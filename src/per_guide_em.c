
/* per_guide_em.c — refined per‑guide EM implementation */
#include "per_guide_em.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double log_nb_pmf(int k, double mu, double r){
    if (mu <= 0) mu = 1e-12;
    if (r  <= 0) r  = 1e-6;
    double logp = 0.0;
    logp  = lgamma(k + r) - lgamma(r) - lgamma(k + 1.0);
    double t = r + mu;
    logp += r * log(r / t) + k * log(mu / t);
    return logp;
}

static double logsumexp2(double a, double b){
    double m = (a > b) ? a : b;
    return m + log(exp(a - m) + exp(b - m));
}

static double safe_clip(double x, double lo, double hi){
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static void set_defaults(PGEMParams *P){
    P->max_iters = 50;
    P->tol = 1e-6;
    P->update_disp = 1;
    P->r_bg_init = 50.0;
    P->r_pos_init = 50.0;
    P->r_lo = 5.0;    P->r_hi = 500.0;
    P->a_bg_lo = 1.0; P->a_bg_hi = 50.0;
    P->a_pos_lo = 1e-6; P->a_pos_hi = 50.0;
    P->pi_lo = 1e-4; P->pi_hi = 0.999;
    P->cap_z = 0.99;
    P->winsor_mult = 5.0;
    P->trim_hi_q = 0.995;
    P->ambient_floor = 1e-9;
}

static double quantile_double(const double *x, int n, double q){
    if (n <= 0) return 0.0;
    int idx = (int)floor((n - 1) * (q<0?0:(q>1?1:q)));
    /* copy, partial sort naive for simplicity */
    double *tmp = (double*)malloc(n*sizeof(double));
    for (int i=0;i<n;++i) tmp[i]=x[i];
    /* full sort */
    for (int i=0;i<n-1;++i){
        for (int j=i+1;j<n;++j){
            if (tmp[j] < tmp[i]){ double t=tmp[i]; tmp[i]=tmp[j]; tmp[j]=t; }
        }
    }
    double v = tmp[idx];
    free(tmp);
    return v;
}

int pgem_fit(const int *c, const int *m, int N, double r_g_in,
             const PGEMParams *params_in, PGEMFit *fit, double *P_pos){
    if (N <= 0 || !c || !m || !fit || !P_pos) return 1;
    PGEMParams P = * (params_in ? params_in : &(PGEMParams){0});
    if (!params_in) set_defaults(&P);

    /* exposures with ridge and quantile for trimming */
    double *E = (double*)malloc(N * sizeof(double));
    if (!E) return 2;
    for (int i=0;i<N;++i) E[i] = (m[i] > 0 ? (double)m[i] : 1.0);
    double E_thr = quantile_double(E, N, P.trim_hi_q);

    double r_g = r_g_in;
    if (r_g < P.ambient_floor) r_g = P.ambient_floor;

    /* init scales: a_bg starts at 1; a_pos from upper quantile of c/E */
    double a_bg = 1.0, a_pos = 0.0;
    double max_rat = 0.0;
    for (int i=0;i<N;++i){
        if (E[i] > 0 && c[i] > 0){
            double r = (double)c[i] / E[i];
            if (r > max_rat) max_rat = r;
        }
    }
    a_pos = safe_clip(fmax(max_rat * 0.5, 1e-3), P.a_pos_lo, P.a_pos_hi);

    double r_bg = safe_clip(P.r_bg_init, P.r_lo, P.r_hi);
    double r_pos = safe_clip(P.r_pos_init, P.r_lo, P.r_hi);
    double pi_pos = 0.05;

    double ll_prev = -INFINITY;
    int it;
    for (it=0; it < P.max_iters; ++it){
        /* E-step */
        double ll = 0.0;
        for (int i=0;i<N;++i){
            double mu_bg = a_bg * E[i] * r_g;
            double mu_pos = a_pos * E[i];
            double l0 = log(1.0 - pi_pos) + log_nb_pmf(c[i], mu_bg, r_bg);
            double l1 = log(pi_pos) + log_nb_pmf(c[i], mu_pos, r_pos);
            double lse = logsumexp2(l0, l1);
            double z = exp(l1 - lse);
            if (z > P.cap_z) z = P.cap_z;
            P_pos[i] = z;
            ll += lse;
        }

        /* M-step with winsorization and exposure trimming */
        double sum_z = 0.0, sum1 = 1e-12, sumE1 = 1e-12;
        double sum0 = 1e-12, sumE0 = 1e-12;

        for (int i=0;i<N;++i){
            int use = (E[i] <= E_thr) ? 1 : 0; /* ignore very large exposures */
            double mu_bg = a_bg * E[i] * r_g;
            double mu_pos = a_pos * E[i];
            double z = P_pos[i];

            double c1_eff = (double)c[i];
            double c0_eff = (double)c[i];
            if (P.winsor_mult > 0){
                double cap1 = P.winsor_mult * fmax(mu_pos, 1e-6);
                double cap0 = P.winsor_mult * fmax(mu_bg, 1e-6);
                if (c1_eff > cap1) c1_eff = cap1;
                if (c0_eff > cap0) c0_eff = cap0;
            }

            sum_z += z * use;
            sum1  += z * c1_eff * use;
            sumE1 += z * E[i] * use;

            sum0  += (1.0 - z) * c0_eff * use;
            sumE0 += (1.0 - z) * (E[i] * r_g) * use;
        }

        pi_pos = safe_clip(sum_z / (double)N, P.pi_lo, P.pi_hi);
        if (sumE1 > 0) a_pos = safe_clip(sum1 / sumE1, P.a_pos_lo, P.a_pos_hi);
        if (sumE0 > 0) a_bg  = safe_clip(sum0 / sumE0,  P.a_bg_lo, P.a_bg_hi);

        if (P.update_disp){
            /* method-of-moments update with weights z and (1-z) */
            double w1=1e-9, mu_bar1=0.0, vnum1=0.0;
            double w0=1e-9, mu_bar0=0.0, vnum0=0.0;
            for (int i=0;i<N;++i){
                int use = (E[i] <= E_thr) ? 1 : 0;
                double z = P_pos[i];
                double mu1 = a_pos * E[i];
                double mu0 = a_bg  * E[i] * r_g;
                w1 += z * use; mu_bar1 += z * mu1 * use;
                w0 += (1.0 - z) * use; mu_bar0 += (1.0 - z) * mu0 * use;
            }
            mu_bar1 /= w1; mu_bar0 /= w0;
            for (int i=0;i<N;++i){
                int use = (E[i] <= E_thr) ? 1 : 0;
                double z = P_pos[i];
                double d1 = (double)c[i] - a_pos * E[i];
                double d0 = (double)c[i] - a_bg  * E[i] * r_g;
                vnum1 += z * d1 * d1 * use;
                vnum0 += (1.0 - z) * d0 * d0 * use;
            }
            double var1 = vnum1 / w1;
            double var0 = vnum0 / w0;
            double r_new1 = r_pos, r_new0 = r_bg;
            if (var1 > mu_bar1 + 1e-9){
                r_new1 = (mu_bar1 * mu_bar1) / (var1 - mu_bar1);
                r_new1 = safe_clip(r_new1, P.r_lo, P.r_hi);
            }
            if (var0 > mu_bar0 + 1e-9){
                r_new0 = (mu_bar0 * mu_bar0) / (var0 - mu_bar0);
                r_new0 = safe_clip(r_new0, P.r_lo, P.r_hi);
            }
            r_pos = r_new1; r_bg = r_new0;
        }

        double drel = (ll_prev == -INFINITY) ? 1.0 : fabs(ll - ll_prev) / (fabs(ll_prev) + 1e-9);
        ll_prev = ll;
        if (drel < P.tol) break;
    }

    fit->pi_pos = pi_pos;
    fit->a_bg   = a_bg;
    fit->a_pos  = a_pos;
    fit->r_bg   = r_bg;
    fit->r_pos  = r_pos;
    fit->ll     = ll_prev;
    fit->iters  = it;

    free(E);
    return 0;
}
