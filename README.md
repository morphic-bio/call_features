# call_features
Produces an executable call_features that calls features from different experiments. Initially used for 10x flex sample barcodes, it has been extended to support different methods suitable for perturb-seq and lineage feature calls. Works with process_features which assigns and error corrects barcodes and UMIs and maps features to cell barcodes (GEMs). Open source fast alternative to 10x Cellranger software.

**Documentation:**
- [docs/Technical.md](docs/Technical.md) - Statistical deep dive into algorithms and parameterization
- [docs/APPLY_ALL_SUMMARY.md](docs/APPLY_ALL_SUMMARY.md) - Quick reference for `--apply-all` flag
- [docs/TESTING_COMPLETE.md](docs/TESTING_COMPLETE.md) - Testing results and validation

---

## 1  Building

```bash
make clean          # remove all executables and objects
make all            # serial build, single-thread runtime

# optional: enable OpenMP (multi-thread EM)
make all OPENMP=1    
```

The build produces `call_features` (streaming demux with optional EM)
and `flex_demux_mtx` (stand-alone FLEX implementation).

---

## 2  Command-line options

### 2.1 Core flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--mtx-dir DIR` | **required** | – | Directory with `matrix.mtx`, `features.tsv`, `barcodes.tsv`, `allowed_features.tsv` (or `features.txt` in `--process-features` mode). |
| `--starsolo-dir DIR` | optional | – | STARsolo *Gene* directory that contains `raw/` and `filtered/` barcode lists. If not provided, filtered barcodes are searched in MTX directory. |
| `--out-prefix PFX` | **required** | – | Prefix for all output files, or directory path (if ends with `/` or is existing directory). |
| `--threads N` | int | **1** | Number of OpenMP threads. When compiled without OpenMP this flag is ignored. |
| `--help`, `-h` | switch | – | Print brief help with all flags and exit. |
| `--cell-list FILE` | optional | *none* | Replace STARsolo's filtered list with your own whitelist. |
| `--tau X` | double | **0.8** | Fraction of counts in top feature required for singlet. |
| `--delta X` | double | **0.4** | Minimum gap `f₁ − f₂` between the two top features. |
| `--gamma X` | double | **0.9** | Dominance `(f₁+f₂)/total` threshold for doublet. |
| `--alpha X` | double | **1e-4** | Per-test error rate (Benjamini–Hochberg when FDR is on). |
| `--floor N` | int | **12** | Hard floor on total feature counts considered for classification. |
| `--ambient-q X` | double | **0.999** | High quantile used by the empirical ambient model. |
| `--no-fdr` | switch | *off* | Use raw p-values instead of FDR-corrected q-values. |
| `--m-min-fixed N` | int | *unset* | **Bypass** all ambient/statistical and cell-derived rules; use `M_min = N` for the low-support cutoff. |
| `--apply-all` | switch | *off* | Apply learned thresholds to **all** barcodes in the matrix (not just filtered list). See section 2.5 below. |

### 2.2 Feature-processing & simple assignment

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--process-features` | switch | off | Treat `features.txt` as the canonical list and ignore `allowed_features.*` – useful for antibody capture where every feature is valid. |
| `--simple-assign` | switch | off | Heuristic one-pass assignment (`ratio`, `min-count`) that avoids statistical tests. |
| `--min-count N` | int | **2** | Minimum raw counts for the top feature in simple mode. |
| `--min-ratio R` | double | **2.0** | Required ratio `f₁/f₂` to call a singlet in simple mode. |

### 2.3 Cell-derived *M*<sub>min</sub> options

| Flag | Type | Default | Notes |
|------|------|---------|-------|
| `--mmin-from-cells` | switch | off | Enable data-driven determination of low-support cutoff. |
| `--mmin-cells-method {em-cell}` | string |
| Initialisation for em-cells. |
| `--mmin3-floor N` | int | `--floor` | Floor passed to em-cells internal negative binomials. |

### 2.4 EM-specific options

| Flag | Type | Default | Rationale |
|------|------|---------|-----------|
| `--use-em` | switch | off | Activate per-guide EM. |
| `--min-em-counts N` | int | **10** | Skip EM for guides whose total count across all cells is `< N`; a background-only fit is written instead. |
| `--em-fixed-disp` | switch | off | Keep NB dispersion parameters fixed during EM. |
| `--tau-pos` | double | 0.95 | Posterior threshold to call a **positive** cell for a guide (primary threshold). |
| `--tau-pos2` | double | `tau-pos − 0.05` | Looser posterior for the 2nd guide. Helps rescue real doublets that have slightly weaker evidence for the second feature. |
| `--k-min` | int | 4 | Min raw counts for a guide to be considered present (primary). |
| `--k-min2` | int | `max(2,k-min-1)` | Min counts for the 2nd guide. Allows doublets with modest secondary counts. |
| `--gamma-min` | double | 0.8 | Dominance threshold `(c1+c2)/total` using **all** features. Ensures the two top guides dominate the library for a doublet. |
| `--gamma-min-cand` | double | 0.85 | Same dominance threshold but denominator is **candidate-only** counts. Captures true doublets in heavy ambient backgrounds where total dominance may be diluted. |
| `--doublet-balance` | 0/1 | 1 | Require 20–80 % balance between the two guides. Set to `0` to keep highly unbalanced doublets. |
| `--debug-amb FILE` | path | – | Write CSV of ambiguous cells and rejection reasons – invaluable when tuning the above knobs. |
| `--k-small` | int | 4 | Minimum number of features required for EM convergence. |

**Why two dominance thresholds?**  In clean data `(c1+c2)/total` ≥ `gamma-min` is robust. Under heavy ambient noise the candidate-only ratio is more reliable; the algorithm considers a cell a doublet if **either** condition passes.

**Balance check** (`--doublet-balance`): avoids misclassifying singlets with trace counts of a second guide as doublets. Disabling it loosens the call when barcode balance is not informative (e.g. targeted panels with skewed expression).

### 2.5 Apply-all mode (`--apply-all`)

By default, `call_features` trains thresholds on the **filtered barcode list** (from `--cell-list`, MTX directory, or STARsolo) and classifies only those barcodes. The `--apply-all` flag extends this workflow:

1. **Training phase**: Thresholds, ambient model, and EM parameters are learned using only the filtered barcodes (unchanged).
2. **Application phase**: The learned model is then applied to **every barcode** in the matrix that passes quality filters (`M_min`).

**Key behaviors:**
- Output files (`assignments.tsv`, `doublets.txt`, etc.) contain results for both filtered and unfiltered barcodes
- Barcodes already processed in the filtered set are skipped with a warning if duplicates are detected
- Quality filters (`M_min`) still apply—barcodes with insufficient counts are excluded
- For **FLEX mode**: Computes p-values vs ambient for each barcode and applies the same decision rules (`f1 ≥ tau`, `(f1-f2) ≥ delta`, `q < alpha`) used in training. Uses raw p-values (no FDR recomputation across expanded set)
- For **EM mode**: Uses stored EM fit parameters (a_bg, a_pos, r_bg, r_pos, pi_pos) to evaluate each guide, then applies the same candidate-counting and dominance logic as the main path
- For **simple-assign**: Trivially extends to all barcodes (no training phase)

**Memory consideration**: `--apply-all` allocates arrays sized to the full matrix (`n_cols`), which can be ~1.6M for large datasets. This is acceptable for correctness-first workflows.

**When to use:**
- You want to "rescue" barcodes that didn't make the filtered list but have clear feature assignments
- Running universe rescue workflows where the filtered list is conservative
- Exploratory analysis to see how many additional cells can be called

**Example:**
```bash
call_features \
  --mtx-dir HTO/matrix \
  --cell-list filtered_barcodes.tsv \
  --out-prefix results/rescue \
  --apply-all \
  --tau 0.8 --delta 0.4 --gamma 0.9
```

---

## 3  Algorithms and when to use them 🌱

| Mode | Key idea | Strengths | Guardrails / when **not** to use |
|------|----------|-----------|----------------------------------|
| **FLEX** (default) | Binomial test vs ambient + dominance rules | Simple, fast, matches Cell Ranger HTO pipeline. | High ambient or low-MOI guide data will inflate false negatives. |
| **EM** | Per-guide NB2 mixture with learnt ambient & positive means | Handles overdispersion, great for low-count CRISPR screens. | Slow on >50 k guides; needs ≥ 2–3 counts in positives (`k-min`). |
| **EM-fixed disp** | Same but freeze dispersion | Stabilises fits in sparse data. | Slight under-fit when guides are highly expressed. |
| **Simple assign** | Pure ratio (`min-ratio`, `min-count`) | Blazing fast, good for cell hashing with distinct HTOs. | No p-values/FDR; unsafe for noisy perturb-seq. |
| **Cell-derived M_min** <br> *EM-cell* | NB mixture (A/S/D components) | Best at separating ambient / singlet / doublet peaks. | More CPU; fails on <10 k cells. |

## 4  Example recipes 🍳

### 4.1 10x **FLEX** cell-hashing (4–12 HTOs)

Standard cell hashing with hashtag oligonucleotides. Each cell should express 1-2 HTOs clearly above ambient.

```bash
call_features \
  --mtx-dir HTO/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/hto_run1 \
  --tau 0.8 --delta 0.4 --gamma 0.9 --floor 12
```

If there are more HTOs than 4 the EM option is possible and may give more coverage – see `scripts/runCallEm.sh` for a more extensive, tuned example.

```bash
call_features \
  --mtx-dir HTO/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/hto_run1 \
  --tau 0.8 --delta 0.4 --gamma 0.9 --floor 12
  --use-em
  --em-fixed-disp
```

### 4.2 Low-MOI CRISPR Perturb-seq (≈3 k guides)

CRISPR screens where most cells are expected to have a single guide. EM handles sparse counts and ambient contamination better than binomial tests.

```bash
call_features \
  --mtx-dir Guides/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/perturb_EM \
  --use-em --k-min 4 --tau-pos 0.95 --tau-pos2 0.9 \
  --gamma-min 0.8 --gamma-min-cand 0.85 \
  --em-fixed-disp \
  --mmin-from-cells --mmin-cells-method model3 --mmin-cap 150
```

### 4.3 High-complexity lineage "Larry" barcodes (250 k possibilities)

Lineage tracing with hundreds of thousands of possible barcodes. Ambient noise tends to be much lower. As a result most cells have unique barcodes or a clearly dominant barcode, so simple ratio-based assignment is sufficient. The ratio can be adjusted. The default ratio seem to be a good compromise between coverage and noise. Alternatively, EM also works very well but is a bit slower with the larger number of barcodes.

```bash
call_features \
  --mtx-dir LB/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/larry \
  --simple-assign --min-count 4 --min-ratio 2 \
  --m-min-fixed 4           # low counts are expected
```

### 4.4 Cell Ranger equivalents 📋

Cell Ranger uses at least *two* different demultiplexing engines, depending on the assay. We can approximate the behavior with the following paramters.

#### 4.4.1 Feature-Barcode (HTO / antibody hashing)

Cell Ranger uses binomial + dominance rules (**FLEX**) with:
- Poisson ambient tail
- Hard floor ≈ 10
- Singlet `τ = 0.9`, doublet if `f₁+f₂ ≥ 0.8`

```bash
call_features \
  --mtx-dir FLEX --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/cellRanger_FLEX \
  --tau 0.9 --gamma 0.8 --floor 10 --no-fdr
```

#### 4.4.2 CRISPR (sgRNA counts)

Cell Ranger uses EM with *Gaussian* mixture on log-counts, hard read cut-off = **3**, dispersion fixed: can use a slightly higher cutoff e.g. 6 to approximate lower CellRanger coverage. 

```bash
call_features \
  --mtx-dir Guides/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/cellRanger_CRISPR_like \
  --use-em \
  --m-min-fixed 3 \  
  --k-min 3 --k-min2 2 \
  --tau-pos 0.9 --tau-pos2 0.85 \
  --gamma-min 0.8 --gamma-min-cand 0.85 \
  --em-fixed-disp
```

#### Gaussian (CellRanger) vs Negative-Binomial (call_features) mixtures — why we use NB on raw counts

| Aspect | Cell Ranger *(Gaussian, log-space)* | call_features *(NB2, count-space)* |
|--------|-------------------------------------|-------------------------------------|
| Data input | `log10(count + 1)` | integer counts |
| Tail behaviour | Symmetric; under-penalises zeros | Heavy-tailed; captures over-dispersion & excess zeros |
| Interpretation | Means/σ in log space | Direct counts ⇒ easier QC |
| Advantage | Very fast | Better fit for low UMIs; no log artifacts |

In practice we observe **higher recall** for low-abundance guides and cleaner ambient
separation with the NB mixture, while still matching Cell Ranger’s hard floor and posterior
thresholds.

---

## 5  Algorithm internals (bird's-eye view)

1. **Cell barcodes**  
   • Load MTX barcodes (columns) into a hash.  
   • Load STARsolo *raw* and *filtered* barcodes.  
   • Optionally, replace the filtered set with `--cell-list` file.

2. **Feature allow-list**  
   Read `allowed_features.tsv`; map matrix rows → 1…K indices. In `--process-features` mode, all features from `features.txt` are considered allowed.

3. **Filtered barcode resolution**  
   Search for filtered cell barcodes in priority order:  
   • **Priority 1**: `--cell-list FILE` (user-provided)  
   • **Priority 2**: `filtered_barcodes.tsv` or `filtered_barcodes.txt` in MTX directory  
   • **Priority 3**: `--starsolo-dir/filtered/barcodes.tsv` (STARsolo fallback)  
   If none found, exit with guidance.

4. **Single scan of `matrix.mtx`**  
   • Ambient counts: features × negatives (raw - filtered barcodes).  
   • Per-cell feature vectors (only for barcodes present in the filtered set).  

5. **Ambient model & thresholds**  
   • Trimmed Poisson pre-filter on negative totals.  
   • Estimate `M_ambient`, compare with NB model if overdispersed.  
   • Statistical power threshold `M_stat` via binomial tail on highest ambient proportion.  
   • Final `M_min = max(M_ambient, M_stat, floor)` unless `--m-min-fixed` or `--mmin-from-cells` override.

6. **Per-cell tests**  
   • **FLEX mode**: For each cell, compute p-value of top feature (`p1`) and second (`p2`) vs ambient.  
   • **EM mode**: Fit per-guide NB mixture models, compute posterior probabilities.  
   • **Simple mode**: Direct ratio and count thresholds.  
   • FDR-correct (unless `--no-fdr`) → `q1`, `q2`.  

7. **Classification**  
   • **Singlet** if `f1 ≥ τ`, `f1-f2 ≥ δ`, `q1 < α`.  
   • **Doublet** if `f1+f2 ≥ γ`, balance check passes, and both `q1,q2 < α`.  
   • Else **Ambiguous**.  
   • Cells below `M_min` → **Low-support**.

Outputs are written as five plain-text files:
- `.assignments.tsv` - Main results with cell classifications
- `.doublets.txt` - List of doublet barcodes  
- `.ambiguous.txt` - Ambiguous assignments
- `.unassignable.txt` - Low-support cells
- `.missing_cells.txt` - Cells in whitelist but not in MTX

### Output Modes

The `--out-prefix` parameter supports two modes:

**Prefix mode (traditional)**: `--out-prefix myrun` creates:
- `myrun.assignments.tsv`, `myrun.doublets.txt`, etc.

**Directory mode**: `--out-prefix results/` or `--out-prefix results` (if `results` exists) creates:
- `results/assignments.tsv`, `results/doublets.txt`, etc.
- The directory is automatically created if it doesn't exist.

---

## 6  Testing and validation

The repository includes test scripts:

```bash
# Test individual methods
./scripts/testCallPerturb.sh {flex|em|em-fixed|simple|em-cell}

# Compare all methods side-by-side
./scripts/compareAllMethods.sh
```

The comparison script generates a detailed table showing how different algorithms perform on the same dataset, helping you choose the best approach for your experiment type.

For a full flag reference and method details, see docs/Methods.md.

---

## 7  Design notes

* **Single-pass MTX** – avoids large memory; scales to millions of barcodes.  
* **Robust ambient estimation** – mixture of empirical and model-based thresholds copes with overdispersion.  
* **Multiple algorithms** – from simple ratio-based to sophisticated EM mixture models.
* **No external dependencies/fast compile** – only glibc and `<math.h>`.  
---

## 8  License

Licensed under the MIT license (see `LICENSE` file).  
Contributions welcome via pull request.
