/*
 * compute_Mmin_from_cells.h
 *
 * Cell-derived ambient threshold (M_min) computation module.
 * Provides three strategies: Otsu, Quantile, and 3-component NB mixture (Model3).
 * 
 * Replaces the ambient-driven minimum UMI threshold with a floor estimated 
 * from the distribution of per-cell total feature counts.
 */

#ifndef COMPUTE_MMIN_FROM_CELLS_H
#define COMPUTE_MMIN_FROM_CELLS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Method identifiers */
typedef enum {
    MMIN_METHOD_OTSU = 0,
    MMIN_METHOD_QUANTILE = 1,
    MMIN_METHOD_MODEL3 = 2
} MminMethod;

/* Model3 initialization strategies */
typedef enum {
    MMIN3_INIT_QUANTILES = 0,
    MMIN3_INIT_KMEANS = 1
} Mmin3InitMethod;

/* Parameters for Model3 (3-component NB mixture) */
typedef struct {
    int max_iters;          /* default: 50 */
    double tol;             /* default: 1e-6 */
    int update_disp;        /* 0 or 1, default: 0 */
    Mmin3InitMethod init_method; /* default: MMIN3_INIT_QUANTILES */
    int floor_val;          /* minimal allowed outcome */
    /* NB dispersion clamps */
    double r_lo;            /* default: 5.0 */
    double r_hi;            /* default: 500.0 */
} Mmin3Params;

/* Model3 fit results */
typedef struct {
    double pi_A, pi_S, pi_D;    /* mixing proportions for Ambient/Singlet/Doublet */
    double mu_A, mu_S, mu_D;    /* means */
    double r_A, r_S, r_D;       /* NB dispersions */
    double loglik;              /* final log-likelihood */
    int iters;                  /* iterations taken */
    int converged;              /* 1 if converged, 0 if failed */
} Mmin3Fit;

/*
 * Main function to compute M_min from cell totals.
 * 
 * Parameters:
 *   cell_tot: array of per-cell total feature counts (length C)
 *   C: number of cells
 *   floor_val: minimum allowed M_min value
 *   mmin_cap: maximum allowed M_min value
 *   method: MMIN_METHOD_OTSU, MMIN_METHOD_QUANTILE, or MMIN_METHOD_MODEL3
 *   qcells: quantile for QUANTILE method (e.g., 0.60)
 *   m3_params: parameters for Model3 method (can be NULL for defaults)
 *   M_min: output - computed M_min value
 *   m3_fit: output for Model3 results (can be NULL if not needed)
 * 
 * Returns: 0 on success, non-zero on error
 */
int compute_Mmin_from_cells(
    const int *cell_tot, int C,
    int floor_val, int mmin_cap,
    MminMethod method, double qcells,
    const Mmin3Params *m3_params,
    int *M_min,
    Mmin3Fit *m3_fit
);

/* Helper function to parse method string to enum */
MminMethod parse_mmin_method(const char *method_str);

/* Helper function to parse init method string to enum */
Mmin3InitMethod parse_mmin3_init(const char *init_str);

/* Helper function to set default Model3 parameters */
void set_default_mmin3_params(Mmin3Params *params, int floor_val);

#ifdef __cplusplus
}
#endif

#endif /* COMPUTE_MMIN_FROM_CELLS_H */
