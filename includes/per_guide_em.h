
/* per_guide_em.h — refined per‑guide EM (NB mixture with exposure + ambient)
 *
 * Fit for one guide g across N cells:
 *   c_i ~ (1 - pi_pos) * NB(mu_bg_i, r_bg)  +  pi_pos * NB(mu_pos_i, r_pos)
 *   mu_bg_i = a_bg * m_i * r_g
 *   mu_pos_i= a_pos * m_i
 * Returns posteriors P_pos[i] = P(component=positive | c_i).
 *
 * Robust options:
 *   - cap_z:        cap responsibility z_i ≤ cap_z (default 0.99)
 *   - winsor_mult:  winsorize counts in M‑step at winsor_mult × mu_comp (default 5.0)
 *   - trim_hi_q:    ignore top quantile of exposures in M‑step (default 0.995)
 *
 * MIT‑style license.
 */

#ifndef PER_GUIDE_EM_H
#define PER_GUIDE_EM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int    max_iters;      /* default 50 */
    double tol;            /* relative LL tol per cell, default 1e-6 */
    int    update_disp;    /* 1 = update r_bg/r_pos by MoM; 0 = keep fixed */

    /* initial values and clamps */
    double r_bg_init;      /* default 50 */
    double r_pos_init;     /* default 50 */
    double r_lo, r_hi;     /* clamp for r: [5, 500] */
    double a_bg_lo, a_bg_hi;   /* [1, 50] */
    double a_pos_lo, a_pos_hi; /* [1e-6, 50] */
    double pi_lo, pi_hi;   /* [1e-4, 0.999] */

    /* robustness */
    double cap_z;          /* e.g., 0.99 */
    double winsor_mult;    /* e.g., 5.0 */
    double trim_hi_q;      /* e.g., 0.995 */

    /* ambient floor (avoid degenerate r_g=0) */
    double ambient_floor;  /* e.g., 1e-9 */
} PGEMParams;

typedef struct {
    double pi_pos;
    double a_bg;
    double a_pos;
    double r_bg;
    double r_pos;
    double ll;
    int    iters;
} PGEMFit;

/* Fit EM for a single guide.
 * Inputs:
 *   c[N]  : counts for this guide per cell
 *   m[N]  : total guide counts per cell (exposure), N > 0
 *   N     : number of cells
 *   r_g   : ambient composition for this guide in (0,1); floored internally
 *   params: hyperparameters (pass NULL for defaults)
 * Outputs:
 *   fit   : fitted parameters
 *   P_pos : length-N array of posteriors; caller allocates (double P_pos[N])
 * Returns 0 on success; non-zero on failure.
 */
int pgem_fit(const int *c, const int *m, int N, double r_g,
             const PGEMParams *params, PGEMFit *fit, double *P_pos);

#ifdef __cplusplus
}
#endif

#endif /* PER_GUIDE_EM_H */
