# assign_features
Produces an executable call_features that calls features from different experiments. Initially used for 10x flex sample barcodes, it has been extended to support different methods suitable for perturb-seq and lineage feature calls. Works with process_features which assigns and error corrects barcodes and UMIs and maps features to cell barcodes (GEMs). Open source fast alternative to 10x Cellranger software.

---

## 1  Building

```bash
make clean          # remove all executables and objects
make all            # compile all executables
```

The build produces `call_features`.

---

## 2  Command-line options

### 2.1 Core flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--mtx-dir DIR` | **required** | – | Directory with `matrix.mtx`, `features.tsv`, `barcodes.tsv`, `allowed_features.tsv` (or `features.txt` in `--process-features` mode). |
| `--starsolo-dir DIR` | **required** | – | STARsolo *Gene* directory that contains `raw/` and `filtered/` barcode lists. |
| `--out-prefix PFX` | **required** | – | Prefix for all output files. |
| `--cell-list FILE` | optional | *none* | Replace STARsolo's filtered list with your own whitelist. |
| `--tau X` | double | **0.8** | Fraction of counts in top feature required for singlet. |
| `--delta X` | double | **0.4** | Minimum gap `f₁ − f₂` between the two top features. |
| `--gamma X` | double | **0.9** | Dominance `(f₁+f₂)/total` threshold for doublet. |
| `--alpha X` | double | **1e-4** | Per-test error rate (Benjamini–Hochberg when FDR is on). |
| `--floor N` | int | **12** | Hard floor on total feature counts considered for classification. |
| `--ambient-q X` | double | **0.999** | High quantile used by the empirical ambient model. |
| `--no-fdr` | switch | *off* | Use raw p-values instead of FDR-corrected q-values. |
| `--m-min-fixed N` | int | *unset* | **Bypass** all ambient/statistical and cell-derived rules; use `M_min = N` for the low-support cutoff. |

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
| `--mmin-cells-method {otsu\|quantile\|model3}` | string | **otsu** | Histogram (Otsu), simple quantile, or three-component NB mixture. |
| `--mmin-qcells Q` | double | **0.60** | Quantile used by the *quantile* method. |
| `--mmin-cap N` | int | **1 000 000** | Safety ceiling: cutoffs above this are clamped. |
| `--mmin3-max-iters` | int | **50** | Model3 EM iterations. |
| `--mmin3-tol` | double | **1e-6** | Convergence tolerance. |
| `--mmin3-update-disp` | 0/1 | **0** | Update NB dispersion parameters during Model3 EM. |
| `--mmin3-init {quantiles\|kmeans}` | string | **quantiles** | Initialisation for Model3. |
| `--mmin3-floor N` | int | `--floor` | Floor passed to Model3 internal negative binomials. |

### 2.4 EM-specific options

| Flag | Type | Default | Rationale |
|------|------|---------|-----------|
| `--use-em` | switch | off | Activate per-guide EM; needed for low-MOI Perturb-seq where ambient NB mixture out-performs the simple binomial model. |
| `--em-fixed-disp` | switch | off | Keeps NB dispersion parameters fixed during EM to avoid over-fitting when guide counts are sparse. |
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

---

## 3  Algorithms and when to use them 🌱

| Mode | Key idea | Strengths | Guardrails / when **not** to use |
|------|----------|-----------|----------------------------------|
| **FLEX** (default) | Binomial test vs ambient + dominance rules | Simple, fast, matches Cell Ranger HTO pipeline. | High ambient or low-MOI guide data will inflate false negatives. |
| **EM** | Per-guide NB2 mixture with learnt ambient & positive means | Handles overdispersion, great for low-count CRISPR screens. | Slow on >50 k guides; needs ≥ 2–3 counts in positives (`k-min`). |
| **EM-fixed disp** | Same but freeze dispersion | Stabilises fits in sparse data. | Slight under-fit when guides are highly expressed. |
| **Simple assign** | Pure ratio (`min-ratio`, `min-count`) | Blazing fast, good for cell hashing with distinct HTOs. | No p-values/FDR; unsafe for noisy perturb-seq. |
| **Cell-derived M_min** <br> *otsu* | Histogram of `log1p(total)` | Data-adaptive cut-off without param tuning. | Unimodal / heavy-tail distributions push cutoff too high. |
| **Cell-derived M_min** <br> *quantile* | Fixed percentile of total counts | Very predictable; works when abundance is heavy-tailed but unimodal. | Needs a chosen quantile (`--mmin-qcells`). |
| **Cell-derived M_min** <br> *model3* | NB mixture (A/S/D components) | Best at separating ambient / singlet / doublet peaks. | More CPU; fails on <10 k cells. |

## 4  Example recipes 🍳

### 4.1 10x **FLEX** cell-hashing (4–12 HTOs)

Standard cell hashing with hashtag oligonucleotides. Each cell should express 1-2 HTOs clearly above ambient.

```bash
call_features \
  --mtx-dir HTO/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/hto_run1 \
  --tau 0.8 --delta 0.4 --gamma 0.9 --floor 12
```

If there are more HTOs than 4 the EM option is possible and may give more coverage - basic command is given but look at scripts/runAssignEm.sh for more extensive parameters

```bash
call_features \
  --mtx-dir HTO/matrix --starsolo-dir Sample/Solo.out/Gene \
  --out-prefix results/hto_run1 \
  --tau 0.8 --delta 0.4 --gamma 0.9 --floor 12
  --use-em
  --em-fixed-disp
```

### 4.2 Low-MOI CRISPR Perturb-seq (≈3 k guides)

CRISPR screens where most cells have 0-1 guides. EM handles sparse counts and ambient contamination better than binomial tests.

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

Lineage tracing with hundreds of thousands of possible barcodes. Most cells have unique barcodes, so simple ratio-based assignment is sufficient. The ratio can be adjusted. The defaults seem to be a good compromise between coverage and noise.

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

3. **Single scan of `matrix.mtx`**  
   • Ambient counts: features × negatives (raw - filtered barcodes).  
   • Per-cell feature vectors (only for barcodes present in the filtered set).  

4. **Ambient model & thresholds**  
   • Trimmed Poisson pre-filter on negative totals.  
   • Estimate `M_ambient`, compare with NB model if overdispersed.  
   • Statistical power threshold `M_stat` via binomial tail on highest ambient proportion.  
   • Final `M_min = max(M_ambient, M_stat, floor)` unless `--m-min-fixed` or `--mmin-from-cells` override.

5. **Per-cell tests**  
   • **FLEX mode**: For each cell, compute p-value of top feature (`p1`) and second (`p2`) vs ambient.  
   • **EM mode**: Fit per-guide NB mixture models, compute posterior probabilities.  
   • **Simple mode**: Direct ratio and count thresholds.  
   • FDR-correct (unless `--no-fdr`) → `q1`, `q2`.  

6. **Classification**  
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

---

## 6  Testing and validation

The repository includes test scripts:

```bash
# Test individual methods
./scripts/testAssignPerturb.sh {flex|em|em-fixed|simple|otsu|quantile|model3}

# Compare all methods side-by-side
./scripts/compareAllMethods.sh
```

The comparison script generates a detailed table showing how different algorithms perform on the same dataset, helping you choose the best approach for your experiment type.

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
