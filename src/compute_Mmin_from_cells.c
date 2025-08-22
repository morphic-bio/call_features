/*
 * compute_Mmin_from_cells.c
 *
 * Implementation of cell-derived ambient threshold (M_min) computation.
 * Provides Otsu, Quantile, and 3-component NB mixture methods.
 */

#define _USE_MATH_DEFINES
#include "compute_Mmin_from_cells.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Safe memory allocation */
static void* safe_malloc(size_t size) {
    void *ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        exit(1);
    }
    return ptr;
}

/* Safe clamp function */
static double safe_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/* Integer comparison for qsort */
static int cmp_int(const void *a, const void *b) {
    int x = *(const int*)a;
    int y = *(const int*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Double comparison for qsort */
static int cmp_double(const void *a, const void *b) {
    double x = *(const double*)a;
    double y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Helper function to parse method string */
MminMethod parse_mmin_method(const char *method_str) {
    if (!method_str) return MMIN_METHOD_OTSU;
    if (strcmp(method_str, "otsu") == 0) return MMIN_METHOD_OTSU;
    if (strcmp(method_str, "quantile") == 0) return MMIN_METHOD_QUANTILE;
    if (strcmp(method_str, "model3") == 0) return MMIN_METHOD_MODEL3;
    return MMIN_METHOD_OTSU;  /* default */
}

/* Helper function to parse init method string */
Mmin3InitMethod parse_mmin3_init(const char *init_str) {
    if (!init_str) return MMIN3_INIT_QUANTILES;
    if (strcmp(init_str, "quantiles") == 0) return MMIN3_INIT_QUANTILES;
    if (strcmp(init_str, "kmeans") == 0) return MMIN3_INIT_KMEANS;
    return MMIN3_INIT_QUANTILES;  /* default */
}

/* Set default Model3 parameters */
void set_default_mmin3_params(Mmin3Params *params, int floor_val) {
    if (!params) return;
    params->max_iters = 50;
    params->tol = 1e-6;
    params->update_disp = 0;
    params->init_method = MMIN3_INIT_QUANTILES;
    params->floor_val = floor_val;
    params->r_lo = 5.0;
    params->r_hi = 500.0;
}

/* ============================================================================
 * OTSU METHOD
 * ============================================================================ */

static int compute_mmin_otsu(const int *cell_tot, int C, int floor_val, int mmin_cap, int *M_min) {
    if (C <= 0) {
        *M_min = floor_val;
        return 0;
    }

    /* Build log1p(m) array */
    double *x = (double*)safe_malloc(C * sizeof(double));
    double min_x = DBL_MAX, max_x = -DBL_MAX;
    
    for (int i = 0; i < C; i++) {
        x[i] = log1p((double)cell_tot[i]);
        if (x[i] < min_x) min_x = x[i];
        if (x[i] > max_x) max_x = x[i];
    }
    
    /* Check for degenerate case */
    if (max_x == min_x) {
        int min_m = cell_tot[0];
        for (int i = 1; i < C; i++) {
            if (cell_tot[i] < min_m) min_m = cell_tot[i];
        }
        *M_min = (min_m > floor_val) ? min_m : floor_val;
        if (*M_min > mmin_cap) *M_min = mmin_cap;
        free(x);
        return 0;
    }
    
    /* Build histogram with 256 bins */
    const int B = 256;
    double bin_width = (max_x - min_x) / B;
    int *hist = (int*)calloc(B, sizeof(int));
    
    for (int i = 0; i < C; i++) {
        int bin = (int)((x[i] - min_x) / bin_width);
        if (bin >= B) bin = B - 1;
        if (bin < 0) bin = 0;
        hist[bin]++;
    }
    
    /* Compute cumulative statistics */
    int *w = (int*)calloc(B, sizeof(int));
    double *mu = (double*)calloc(B, sizeof(double));
    
    w[0] = hist[0];
    mu[0] = hist[0] * (min_x + 0.5 * bin_width);
    
    for (int b = 1; b < B; b++) {
        w[b] = w[b-1] + hist[b];
        double bin_center = min_x + (b + 0.5) * bin_width;
        mu[b] = mu[b-1] + hist[b] * bin_center;
    }
    
    int W = w[B-1];
    
    /* Find optimal threshold by maximizing between-class variance */
    double max_sigma_b = -1.0;
    int best_t = 0;
    
    for (int t = 0; t < B-1; t++) {
        if (w[t] == 0 || w[t] == W) continue;
        
        double p1 = (double)w[t] / W;
        double p2 = 1.0 - p1;
        double mu1 = mu[t] / w[t];
        double mu2 = (mu[B-1] - mu[t]) / (W - w[t]);
        
        double sigma_b = p1 * p2 * (mu1 - mu2) * (mu1 - mu2);
        
        if (sigma_b > max_sigma_b) {
            max_sigma_b = sigma_b;
            best_t = t;
        }
    }
    
    /* Back-transform threshold */
    double thr_log = min_x + (best_t + 0.5) * bin_width;
    double thr = exp(thr_log) - 1.0;
    int M_cells = (int)round(thr);
    
    if (M_cells < floor_val) M_cells = floor_val;
    if (M_cells > mmin_cap) M_cells = mmin_cap;
    
    *M_min = M_cells;
    
    /* Cleanup */
    free(x);
    free(hist);
    free(w);
    free(mu);
    
    return 0;
}

/* ============================================================================
 * QUANTILE METHOD
 * ============================================================================ */

static int compute_mmin_quantile(const int *cell_tot, int C, int floor_val, int mmin_cap, 
                                 double qcells, int *M_min) {
    if (C <= 0) {
        *M_min = floor_val;
        return 0;
    }
    
    /* Copy and sort */
    int *m_sorted = (int*)safe_malloc(C * sizeof(int));
    memcpy(m_sorted, cell_tot, C * sizeof(int));
    qsort(m_sorted, C, sizeof(int), cmp_int);
    
    /* Compute quantile */
    int idx = (int)floor(qcells * (C - 1));
    if (idx < 0) idx = 0;
    if (idx >= C) idx = C - 1;
    
    int M_cells = m_sorted[idx];
    if (M_cells < floor_val) M_cells = floor_val;
    if (M_cells > mmin_cap) M_cells = mmin_cap;
    
    *M_min = M_cells;
    
    free(m_sorted);
    return 0;
}

/* ============================================================================
 * MODEL3 METHOD - 3-component NB mixture
 * ============================================================================ */

/* Gamma function helpers */
static double log_gamma_approx(double x) {
    /* Stirling's approximation for log(Gamma(x)) */
    if (x < 1.0) return log_gamma_approx(x + 1.0) - log(x);
    return 0.5 * log(2.0 * M_PI / x) + x * log(x) - x;
}

/* Log NB PMF using NB2 parameterization */
static double log_nb_pmf(int k, double mu, double r) {
    if (mu <= 0) mu = 1e-12;
    if (r <= 0) r = 1e-6;
    
    double log_pmf = 0.0;
    log_pmf += log_gamma_approx(k + r) - log_gamma_approx(r) - log_gamma_approx(k + 1.0);
    log_pmf += r * log(r / (r + mu)) + k * log(mu / (r + mu));
    
    return log_pmf;
}

/* NB PMF (not log) */
static double nb_pmf(int k, double mu, double r) {
    return exp(log_nb_pmf(k, mu, r));
}

/* 1D k-means for initialization */
static void kmeans_1d(const double *y, int n, int k, double *centers) {
    if (n <= 0 || k <= 0) return;
    
    /* Initialize centers with quantiles */
    double *y_sorted = (double*)safe_malloc(n * sizeof(double));
    memcpy(y_sorted, y, n * sizeof(double));
    qsort(y_sorted, n, sizeof(double), cmp_double);
    
    for (int i = 0; i < k; i++) {
        int idx = (int)((double)i / (k - 1) * (n - 1));
        if (idx >= n) idx = n - 1;
        centers[i] = y_sorted[idx];
    }
    free(y_sorted);
    
    /* Run k-means for 25 iterations */
    int *assign = (int*)safe_malloc(n * sizeof(int));
    
    for (int iter = 0; iter < 25; iter++) {
        /* Assignment step */
        for (int i = 0; i < n; i++) {
            double min_dist = DBL_MAX;
            int best_k = 0;
            for (int j = 0; j < k; j++) {
                double dist = fabs(y[i] - centers[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_k = j;
                }
            }
            assign[i] = best_k;
        }
        
        /* Update centers */
        for (int j = 0; j < k; j++) {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < n; i++) {
                if (assign[i] == j) {
                    sum += y[i];
                    count++;
                }
            }
            if (count > 0) {
                centers[j] = sum / count;
            }
        }
    }
    
    free(assign);
}

static int compute_mmin_model3(const int *cell_tot, int C, int floor_val, int mmin_cap,
                               const Mmin3Params *params, Mmin3Fit *fit, int *M_min) {
    if (C <= 0) {
        *M_min = floor_val;
        if (fit) {
            memset(fit, 0, sizeof(Mmin3Fit));
            fit->converged = 0;
        }
        return 0;
    }
    
    /* Transform to log-space for initialization */
    double *y = (double*)safe_malloc(C * sizeof(double));
    for (int i = 0; i < C; i++) {
        y[i] = log1p((double)cell_tot[i]);
    }
    
    /* Initialize parameters */
    double pi_A = 1.0/3.0, pi_S = 1.0/3.0, pi_D = 1.0/3.0;
    double mu_A, mu_S, mu_D;
    double r_A = 10.0, r_S = 10.0, r_D = 10.0;
    
    /* Initialize means based on method */
    if (params->init_method == MMIN3_INIT_QUANTILES) {
        /* Use quantiles: q20, q50, q80 for A, S, D */
        double *y_sorted = (double*)safe_malloc(C * sizeof(double));
        memcpy(y_sorted, y, C * sizeof(double));
        qsort(y_sorted, C, sizeof(double), cmp_double);
        
        int idx20 = (int)(0.2 * (C - 1));
        int idx50 = (int)(0.5 * (C - 1));
        int idx80 = (int)(0.8 * (C - 1));
        
        /* Back-transform to original space */
        mu_A = exp(y_sorted[idx20]) - 1.0;
        mu_S = exp(y_sorted[idx50]) - 1.0;
        mu_D = exp(y_sorted[idx80]) - 1.0;
        
        free(y_sorted);
    } else {
        /* Use k-means */
        double centers[3];
        kmeans_1d(y, C, 3, centers);
        
        /* Sort centers and back-transform */
        qsort(centers, 3, sizeof(double), cmp_double);
        mu_A = exp(centers[0]) - 1.0;
        mu_S = exp(centers[1]) - 1.0;
        mu_D = exp(centers[2]) - 1.0;
    }
    
    /* Ensure monotonicity */
    if (mu_D < mu_S) mu_D = 2.0 * mu_S;
    if (mu_A <= 0) mu_A = 0.1;
    if (mu_S <= 0) mu_S = 1.0;
    if (mu_D <= 0) mu_D = 2.0;
    
    /* Initialize dispersions using method of moments */
    for (int comp = 0; comp < 3; comp++) {
        double mu_comp = (comp == 0) ? mu_A : (comp == 1) ? mu_S : mu_D;
        double sum_sq = 0.0;
        int count = 0;
        
        /* Simple assignment for initial variance estimate */
        for (int i = 0; i < C; i++) {
            double closest_mu = mu_A;
            double min_dist = fabs(cell_tot[i] - mu_A);
            
            if (fabs(cell_tot[i] - mu_S) < min_dist) {
                closest_mu = mu_S;
                min_dist = fabs(cell_tot[i] - mu_S);
            }
            if (fabs(cell_tot[i] - mu_D) < min_dist) {
                closest_mu = mu_D;
            }
            
            if (fabs(closest_mu - mu_comp) < 1e-6) {
                double diff = cell_tot[i] - mu_comp;
                sum_sq += diff * diff;
                count++;
            }
        }
        
        if (count > 1) {
            double var_comp = sum_sq / (count - 1);
            double r_comp = (mu_comp * mu_comp) / fmax(1e-9, var_comp - mu_comp);
            r_comp = safe_clamp(r_comp, params->r_lo, params->r_hi);
            
            if (comp == 0) r_A = r_comp;
            else if (comp == 1) r_S = r_comp;
            else r_D = r_comp;
        }
    }
    
    /* EM Algorithm */
    double *z_A = (double*)safe_malloc(C * sizeof(double));
    double *z_S = (double*)safe_malloc(C * sizeof(double));
    double *z_D = (double*)safe_malloc(C * sizeof(double));
    
    double prev_ll = -DBL_MAX;
    int converged = 0;
    int iters;
    
    for (iters = 0; iters < params->max_iters; iters++) {
        /* E-step: compute responsibilities */
        double ll = 0.0;
        
        for (int i = 0; i < C; i++) {
            double l_A = log(pi_A) + log_nb_pmf(cell_tot[i], mu_A, r_A);
            double l_S = log(pi_S) + log_nb_pmf(cell_tot[i], mu_S, r_S);
            double l_D = log(pi_D) + log_nb_pmf(cell_tot[i], mu_D, r_D);
            
            /* Numerical stability: subtract max */
            double max_l = l_A;
            if (l_S > max_l) max_l = l_S;
            if (l_D > max_l) max_l = l_D;
            
            double exp_A = exp(l_A - max_l);
            double exp_S = exp(l_S - max_l);
            double exp_D = exp(l_D - max_l);
            double sum_exp = exp_A + exp_S + exp_D;
            
            if (sum_exp > 0) {
                z_A[i] = exp_A / sum_exp;
                z_S[i] = exp_S / sum_exp;
                z_D[i] = exp_D / sum_exp;
            } else {
                z_A[i] = z_S[i] = z_D[i] = 1.0/3.0;
            }
            
            /* Clamp responsibilities */
            z_A[i] = safe_clamp(z_A[i], 1e-6, 1.0 - 1e-6);
            z_S[i] = safe_clamp(z_S[i], 1e-6, 1.0 - 1e-6);
            z_D[i] = safe_clamp(z_D[i], 1e-6, 1.0 - 1e-6);
            
            ll += max_l + log(sum_exp);
        }
        
        /* M-step: update parameters */
        double sum_z_A = 0.0, sum_z_S = 0.0, sum_z_D = 0.0;
        double sum_m_A = 0.0, sum_m_S = 0.0, sum_m_D = 0.0;
        
        for (int i = 0; i < C; i++) {
            sum_z_A += z_A[i];
            sum_z_S += z_S[i];
            sum_z_D += z_D[i];
            
            sum_m_A += z_A[i] * cell_tot[i];
            sum_m_S += z_S[i] * cell_tot[i];
            sum_m_D += z_D[i] * cell_tot[i];
        }
        
        /* Update mixing proportions */
        double total_z = sum_z_A + sum_z_S + sum_z_D;
        if (total_z > 0) {
            pi_A = sum_z_A / total_z;
            pi_S = sum_z_S / total_z;
            pi_D = sum_z_D / total_z;
        }
        
        /* Update means */
        if (sum_z_A > 0) mu_A = sum_m_A / sum_z_A;
        if (sum_z_S > 0) mu_S = sum_m_S / sum_z_S;
        if (sum_z_D > 0) mu_D = sum_m_D / sum_z_D;
        
        /* Ensure means are positive and maintain ordering */
        if (mu_A <= 0) mu_A = 0.1;
        if (mu_S <= 0) mu_S = 1.0;
        if (mu_D <= 0) mu_D = 2.0;
        if (mu_D < mu_S) mu_D = 2.0 * mu_S;
        
        /* Update dispersions if requested */
        if (params->update_disp) {
            /* Method of moments for each component */
            for (int comp = 0; comp < 3; comp++) {
                double mu_comp = (comp == 0) ? mu_A : (comp == 1) ? mu_S : mu_D;
                double *z_comp = (comp == 0) ? z_A : (comp == 1) ? z_S : z_D;
                double sum_z_comp = (comp == 0) ? sum_z_A : (comp == 1) ? sum_z_S : sum_z_D;
                
                double weighted_var = 0.0;
                for (int i = 0; i < C; i++) {
                    double diff = cell_tot[i] - mu_comp;
                    weighted_var += z_comp[i] * diff * diff;
                }
                
                if (sum_z_comp > 0) {
                    weighted_var /= sum_z_comp;
                    if (weighted_var > mu_comp + 1e-9) {
                        double r_new = (mu_comp * mu_comp) / (weighted_var - mu_comp);
                        r_new = safe_clamp(r_new, params->r_lo, params->r_hi);
                        
                        if (comp == 0) r_A = r_new;
                        else if (comp == 1) r_S = r_new;
                        else r_D = r_new;
                    }
                }
            }
        }
        
        /* Check convergence */
        double rel_change = (prev_ll == -DBL_MAX) ? 1.0 : 
                           fabs(ll - prev_ll) / (fabs(prev_ll) + 1e-9);
        
        if (rel_change < params->tol) {
            converged = 1;
            break;
        }
        
        prev_ll = ll;
    }
    
    /* Compute Bayes boundary between A and S */
    int m_star = params->floor_val;
    
    for (int m = 0; m <= 1000; m++) {
        double p_A_given_m = pi_A * nb_pmf(m, mu_A, r_A);
        double p_S_given_m = pi_S * nb_pmf(m, mu_S, r_S);
        
        if (p_S_given_m >= p_A_given_m) {
            m_star = m;
            break;
        }
    }
    
    int M_cells = (m_star > params->floor_val) ? m_star : params->floor_val;
    if (M_cells > mmin_cap) M_cells = mmin_cap;
    
    *M_min = M_cells;
    
    /* Store fit results */
    if (fit) {
        fit->pi_A = pi_A;
        fit->pi_S = pi_S;
        fit->pi_D = pi_D;
        fit->mu_A = mu_A;
        fit->mu_S = mu_S;
        fit->mu_D = mu_D;
        fit->r_A = r_A;
        fit->r_S = r_S;
        fit->r_D = r_D;
        fit->loglik = prev_ll;
        fit->iters = iters;
        fit->converged = converged;
    }
    
    /* Cleanup */
    free(y);
    free(z_A);
    free(z_S);
    free(z_D);
    
    /* Check if EM failed to separate - fallback to quantile */
    if (!converged || fabs(mu_A - mu_S) < 1e-3) {
        fprintf(stderr, "WARNING: Model3 EM failed to converge or separate classes, falling back to quantile method\n");
        return compute_mmin_quantile(cell_tot, C, floor_val, mmin_cap, 0.60, M_min);
    }
    
    return 0;
}

/* ============================================================================
 * MAIN FUNCTION
 * ============================================================================ */

int compute_Mmin_from_cells(
    const int *cell_tot, int C,
    int floor_val, int mmin_cap,
    MminMethod method, double qcells,
    const Mmin3Params *m3_params,
    int *M_min,
    Mmin3Fit *m3_fit
) {
    if (!cell_tot || !M_min || C <= 0) {
        if (M_min) *M_min = floor_val;
        return 1;
    }
    
    switch (method) {
        case MMIN_METHOD_OTSU:
            return compute_mmin_otsu(cell_tot, C, floor_val, mmin_cap, M_min);
            
        case MMIN_METHOD_QUANTILE:
            return compute_mmin_quantile(cell_tot, C, floor_val, mmin_cap, qcells, M_min);
            
        case MMIN_METHOD_MODEL3: {
            Mmin3Params default_params;
            if (!m3_params) {
                set_default_mmin3_params(&default_params, floor_val);
                m3_params = &default_params;
            }
            return compute_mmin_model3(cell_tot, C, floor_val, mmin_cap, m3_params, m3_fit, M_min);
        }
        
        default:
            return compute_mmin_otsu(cell_tot, C, floor_val, mmin_cap, M_min);
    }
}
