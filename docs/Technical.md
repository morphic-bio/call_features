# Technical Notes on Demultiplexing Methods

This document describes the statistical models underpinning the major
assignment modes exposed by the `call_features` toolkit and orchestrated in
`scripts/compareAllMethods.sh`.  The emphasis is on the FLEX (binomial) and EM
(negative-binomial mixture) engines, with supporting notes on the heuristic and
cell-derived *M*<sub>min</sub> variants.  Each section introduces the relevant
probability model, parameter interpretation, and practical implications for
statistical validation.

> **Scope and references.**  The descriptions below draw on public literature
> about feature-barcode demultiplexing (e.g. 10x Genomics whitepapers, Cell
> Ranger documentation, and CRISPR feature selection papers).  Suggested
> reading:
> - 10x Genomics Feature Barcode technology: https://www.10xgenomics.com/resources/white-papers/feature-barcoding-analysis
> - Cell Ranger "Feature Barcode" algorithms overview: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/feature-barcoding
> - Gasperini, M. *et al.* A Genome-wide Framework for mapping gene regulation. *Cell* **176**, 377–390 (2019). https://doi.org/10.1016/j.cell.2018.11.029

## 1. FLEX Binomial Model (`flex`)

### 1.1  Model structure
The FLEX mode assumes that counts for each feature (e.g. hashing oligos or
CRISPR guides) in a given cell largely follow a binomial sampling process once
ambient contamination has been estimated.  Let

- `m_i` be the total feature counts observed in cell `i`.
- `x_{ij}` be the counts for feature `j` in cell `i`.
- `r_j` be the estimated ambient proportion for feature `j` derived from
  background droplets (raw barcodes outside the filtered set).

Under the null hypothesis *H₀* that cell `i` contains only ambient signal, the
counts follow a multinomial distribution with probabilities `r_j`.  FLEX uses a
binomial approximation for the top feature:

```
X_{ij} | H₀ ~ Binomial(m_i, r_j)
```

The tail probability `Pr(X_{ij} ≥ x_{ij})` defines the ambient *p*-value.  After
Benjamini–Hochberg correction (unless `--no-fdr` is used), the smallest
*q*-value must fall below `alpha` to be considered significant.

### 1.2  Parameters and decision rules

- `tau` (`--tau`, default 0.8): minimum fraction `x_{(1)}/m_i` required for the
  top feature `x_{(1)}`.  Intuitively, it demands that at least 80% of counts
  belong to the dominant feature, preventing diffuse ambient reads from driving
  a call.
- `delta` (`--delta`, default 0.4): minimum absolute gap `x_{(1)} - x_{(2)}`
  between the top two features.  This guards against borderline doublets where
  two barcodes are comparable.
- `gamma` (`--gamma`, default 0.9): total dominance threshold
  `(x_{(1)} + x_{(2)}) / m_i`.  A doublet is reported when the two best features
  jointly dominate the cell and both pass the ambient test.
- `alpha` (`--alpha`, default 1e-4): target false-discovery rate for the ambient
  binomial tests.
- `floor` (`--floor`, default 12): minimum total counts to attempt
  classification.  Cells below the floor are marked "low support".
- `ambient-q` (`--ambient-q`, default 0.999): quantile used to set the ambient
  threshold *M*<sub>ambient</sub> from negative droplets.

### 1.3  Statistical workflow
1. Ambient proportions `r_j` are estimated from the background set.
2. A count floor `M_min` is derived as the maximum of:
   - Empirical ambient quantile (via `ambient-q`).
   - Statistical power cutoff obtained by inverting the binomial tail at the
     chosen `alpha`.
   - User-specified `--m-min-fixed`, where applicable.
3. For each cell satisfying `m_i ≥ M_min`:
   - Compute per-feature binomial *p*-values and apply BH correction.
   - Apply thresholds (`tau`, `delta`, `gamma`) to classify as singlet, doublet,
     ambiguous, or low-support.

## 2. EM Negative-Binomial Mixture (`em`, `em-fixed`, `em-noncells`)

### 2.1  Model structure
The EM mode fits a per-feature two-component negative-binomial (NB2) mixture to
model positive cells versus ambient-only cells.  For guide `g` and cell `i`, the
latent assignment `Z_{ig}` ∈ {ambient (A), positive (S)} drives the distribution:

```
C_{ig} | Z_{ig} = A ~ NB(r_bg, p_bg)
C_{ig} | Z_{ig} = S ~ NB(r_pos, p_pos)
```

The EM (Expectation-Maximisation) algorithm iteratively updates the component
parameters `θ = {r_bg, p_bg, r_pos, p_pos, π}`, posterior probabilities
`P(Z_{ig} = S | C_{ig})`, and cell-specific expectations.  Exposure adjustments
use the total counts `m_i` as offsets so that highly sequenced cells are not
penalised.

### 2.2  Parameterisation and defaults
In `call_features`, NB2 is parameterised via `(r, p)` where the mean is
`μ = r (1 - p) / p` and variance `μ + μ² / r`.  The defaults (see
`default_PGEM_params()`) set moderately informative priors on dispersion.

Key user-facing parameters:

- `--use-em`: activates the EM path.
- `--em-fixed-disp`: freeze NB dispersion updates (`r` terms) to stabilise sparse
  datasets.  When omitted, dispersion is optimised per guide.
- `--min-em-counts`: minimum total counts per guide before attempting EM fit;
  low-support guides fall back to ambient-only reporting.

Posterior decision thresholds:

- `tau-pos` (`--tau-pos`): posterior probability `P(Z_{ig} = S | data)` required
  for a primary singlet call.  Default 0.95.
- `tau-pos2` (`--tau-pos2`): looser posterior for a secondary guide when testing
  doublets.  Default `tau_pos - 0.05`.
- `k-min`, `k-min2`: minimum observed counts for primary and secondary guides.
- `gamma-min`, `gamma-min-cand`: dominance ratios for total counts versus
  candidate-only counts; they guard against ambient-heavy cells.
- `doublet-balance`: optional (default on) 20–80% balance constraint to avoid
  calling strongly skewed doublets.
- `k-small`: if the number of allowed features `K` falls below this threshold
  (default 4), EM is disabled to prevent ill-posed fits.

### 2.3  Ambient estimation modes
- `em`: ambient component learned jointly with positives using the filtered
  negatives (STARsolo raw minus filtered cells or MTX-derived background).
- `em-fixed`: same as above but with `--em-fixed-disp` to keep dispersion
  parameters constant.
- `em-noncells`: identical to `em-fixed` but emphasises the use of non-cell
  droplets as the ambient reference; useful when ambient expression differs from
  background estimates derived from candidate cells.

Each iteration proceeds as:
1. E-step: compute posterior weights `w_{ig}` for each cell-guide pair.
2. M-step: update NB parameters via weighted sufficient statistics, subject to
   clipping (caps on `π`, `r` ranges, winsorisation) to avoid degenerate fits.
3. Evaluate convergence using `P.max_iters` and `P.tol` (default 50 iterations,
   tolerance 1e-6).

### 2.4  Decision logic
A cell is called positive for guide `g` when:
- `m_i ≥ M_min` (same cutoff pipeline as FLEX); and
- Posterior thresholds (`tau-pos`, `tau-pos2`) and count floors (`k-min`,
  `k-min2`) are met;
- Dominance rules (gamma thresholds, balance) pass.

Doublets correspond to cells satisfying the secondary thresholds for two guides
simultaneously.  Ambiguous cells fail posterior or dominance checks despite
having sufficient counts, while low-support cells fail `M_min`.

## 3. Simple Ratio Heuristic (`simple`)

The `--simple-assign` mode eschews probabilistic modelling and uses deterministic
criteria based on raw counts:

- `min-count`: minimum count for the top feature; prevents low-abundance noise.
- `min-ratio`: required ratio `x_{(1)} / max(1, x_{(2)})` to declare a singlet.

This is appropriate when features are designed to be mutually exclusive and
strongly expressed (e.g. high-complexity lineage barcodes).  No FDR, ambient, or
posterior calculations are performed; cells failing the conditions are ambiguous
or low-support (`m_i < M_min`).

## 4. Cell-Derived *M*<sub>min</sub> Methods (`otsu`, `quantile`, `model3`)

These methods, invoked via `compareAllMethods.sh`, determine the low-support
threshold directly from filtered cell totals `m_i` rather than relying solely on
ambient droplets.

- `otsu`: applies Otsu’s method on `log1p(m_i)` to split the distribution into
  background and foreground.  This targets bimodal count distributions typical
  of ambient versus true signal.
- `quantile`: sets `M_min` to the specified quantile of the cell-total
  distribution (`--mmin-qcells`, default 0.60).
- `model3`: fits a three-component NB mixture (ambient, singlet, doublet) using
  EM (`compute_Mmin_from_cells`).  Parameters:
  - `--mmin3-max-iters`, `--mmin3-tol`: EM iteration control.
  - `--mmin3-update-disp`: whether to update dispersions.
  - `--mmin3-init`: initialisation strategy (`quantiles` or `kmeans`).
  The resulting `M_min` is chosen at the intersection between ambient and
  singlet components, reflecting the most likely cutoff separating empty and
  occupied droplets.

When `--mmin-from-cells` is active, these data-driven thresholds replace the
standard ambient/power combination.  The toolkit prints summaries such as
mixture weights (`π_A`, `π_S`, `π_D`) and component means for verification.

## 5. Practical Interpretation of Key Parameters

| Parameter | Applies to | Statistical meaning |
|-----------|------------|---------------------|
| `tau` | FLEX | Dominance fraction for top feature: `x_{(1)}/m_i ≥ tau` ensures a strong winner under multinomial sampling. |
| `delta` | FLEX | Absolute gap between top two counts to guard doublets. |
| `gamma` | FLEX | Combined dominance of top two features; doublet rule. |
| `alpha` | FLEX | Target FDR for binomial tail tests (BH-adjusted *q*-values). |
| `ambient-q` | FLEX & EM | High quantile of negative totals to set ambient floor. |
| `tau-pos` | EM | Posterior probability threshold `P(Z = S | data)` for calling a guide, under NB mixture. |
| `tau-pos2` | EM | Secondary posterior for potential doublets. |
| `k-min`, `k-min2` | EM | Minimum observed counts required for trustable EM inference. |
| `gamma-min`, `gamma-min-cand` | EM | Dominance constraints on total counts vs candidate-only counts. |
| `doublet-balance` | EM | Balance requirement (20–80%) for doublet evidence across two guides. |
| `min-count`, `min-ratio` | Simple | Raw-count thresholds for deterministic singlet calls. |
| `mmin-from-cells` (with `otsu`, `quantile`, `model3`) | All modes | Replace default `M_min` with data-driven values derived from filtered cell totals. |

## 6. Possible QC strategies (todos)

1. **Ambient estimation sanity checks.**  Inspect `ambient_colsum` vs EM
   posteriors; a uniform `r_j` implies homogeneous background, while high
   variance suggests dominant features or index hopping.
2. **Binomial tests.**  Recompute `Pr(X ≥ x)` for selected cells using R or
   Python (`scipy.stats.binom.sf`).  Compare BH-adjusted *q*-values to confirm
   implementation.
3. **NB mixture fits.**  For a representative guide:
   - Extract counts `C_{ig}` and total exposures `m_i`.
   - Fit NB mixtures in R (`fitdistrplus::fitdist`) or Python (`statsmodels`) to
     validate `r_bg`, `r_pos`, `π`.  Compare to logged EM summaries.
   - Evaluate posterior convergence; check that log-likelihood increases
     monotonically.
4. **Doublet logic.**  Validate doublet calls by confirming both guides exceed
   `k-min2`/`tau-pos2` and satisfy balance criteria.
5. **Data-driven M_min.**  Plot histograms of `m_i` against reported `M_min`.  In
   the `model3` case, overlay the three NB components and confirm intersection
   points match the printed threshold.

## 7. Apply-All Mode (`--apply-all`)

The `--apply-all` flag modifies the workflow to apply learned thresholds to every barcode in the matrix, not just the filtered list:

### 7.1 Training vs. Application

**Default behavior**: Train thresholds on the filtered barcode list and classify only those barcodes.

**With `--apply-all`**:
1. **Training phase** (unchanged): Ambient model, thresholds (`M_min`, `tau`, `delta`, etc.), and EM parameters are learned using only the filtered barcodes.
2. **Application phase** (new): The learned model is applied to every barcode in the matrix that passes quality filters.

### 7.2 Implementation Details

- **Memory overhead**: Allocates per-column arrays sized to `n_cols` (all barcodes), which can be ~1.6M for large datasets. For FLEX/simple modes, stores feature vectors; for EM mode, stores per-column guide counts.
- **FDR handling**: For FLEX mode, learned q-values are reused without recomputing Benjamini-Hochberg on the expanded set. This preserves the statistical properties learned from the high-quality filtered cells.
- **Duplicate handling**: Barcodes already processed in the filtered set are skipped with a warning.
- **Quality filters**: The same `M_min` threshold applies to all barcodes, ensuring only barcodes with sufficient counts are classified.

### 7.3 Use Cases

- **Universe rescue workflows**: Conservative filtered lists may exclude valid cells with clear feature assignments.
- **Exploratory analysis**: Understand how many additional barcodes could be called if the filtered list were less stringent.
- **Cross-validation**: Compare assignments between filtered and near-filtered barcodes to tune filtering thresholds.

### 7.4 Performance Considerations

For a typical lane with 1.6M barcodes and 10k filtered cells:
- Memory increase: ~6-10 MB per allowed feature (FLEX mode) or ~200 MB total (EM mode with sparse storage)
- Runtime increase: Proportional to the number of unfiltered barcodes that pass `M_min` (typically 10-30% of matrix)

## 8. Summary

- FLEX offers a transparent multinomial/binomial test with interpretable
  dominance thresholds (`tau`, `delta`, `gamma`).
- EM extends this to a per-guide NB mixture, leveraging posterior probabilities
  for finer control in sparse data and doublet detection.
- Simple ratio and cell-derived `M_min` heuristics provide faster alternatives
  when assumptions hold or when count distributions supply clearer cutoffs.
- `--apply-all` extends any method to classify all barcodes in the matrix, enabling
  universe rescue and exploratory workflows.
- All methods share the same output structure; `compareAllMethods.sh` exists to
  benchmark them side-by-side, revealing trade-offs between sensitivity and
  specificity.
