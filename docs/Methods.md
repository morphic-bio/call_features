# Methods and CLI Reference

Authoritative reference for algorithms, flags, and scripts. This is the home for all docs other than `README.md`.

Overview
- Binaries are built into `bin/`:
  - `bin/flex_demux_mtx` — stand-alone FLEX demultiplexing
  - `bin/call_features` — streaming demux with optional per-guide EM and cell-derived M_min
- Build with: `make all` (add `OPENMP=1` for multi-thread EM).
- Output modes: `--out-prefix file` (traditional) or `--out-prefix dir/` (directory mode).

Barcode Sources and Resolution
- Filtered (cells) priority:
  1. `--cell-list FILE`
  2. `<mtx-dir>/filtered_barcodes.tsv` (or `.txt`)
  3. `--starsolo-dir/filtered/barcodes.tsv`
- Raw/negatives:
  - Prefer STARsolo `raw/barcodes.tsv` when `--starsolo-dir` is provided; otherwise use all `barcodes.tsv` from `--mtx-dir` and subtract the filtered set.

Algorithms
- FLEX (default): Binomial vs ambient + dominance rules.
- EM: Per-guide NB2 mixtures with ambient + positive components; optional fixed dispersion.
- Simple: Ratio-based heuristic (`--min-ratio`, `--min-count`).
- Cell-derived M_min: Data-driven low-support cutoff (`otsu`, `quantile`, or `model3`).

Full Flag Reference
- Core
  - `--mtx-dir DIR` (required): directory with `matrix.mtx`, `features.tsv`, `barcodes.tsv`, `allowed_features.tsv` (or `features.txt` with `--process-features`).
  - `--out-prefix PFX` (required): prefix for all outputs.
  - `--starsolo-dir DIR` (optional): STARsolo Gene dir containing `raw/` and `filtered/` barcodes; used if present.
  - `--cell-list FILE` (optional): explicit whitelist; overrides filtered list.
  - `--threads N` (optional): OpenMP threads (when compiled with OpenMP).
  - `--help`/`-h`: print help.
- FLEX-like thresholds
  - `--tau X` (default 0.8), `--delta X` (0.4), `--gamma X` (0.9)
  - `--alpha X` (1e-4), `--ambient-q X` (0.999), `--floor N` (12), `--no-fdr`
  - `--m-min-fixed N`: hard override of the low-support cutoff.
- Simple mode
  - `--simple-assign`, `--min-count N` (2), `--min-ratio R` (2.0)
- EM mode (call_features only)
  - `--use-em`, `--em-fixed-disp`, `--min-em-counts N` (10)
  - Positivity/dominance: `--tau-pos X` (0.95), `--tau-pos2 X` (tau-pos−0.05), `--k-min N` (4), `--k-min2 N` (`max(2,k-min-1)`)
  - Doublet rules: `--gamma-min X` (0.8), `--gamma-min-cand X` (0.85), `--doublet-balance {0|1}` (1)
  - Safety: `--k-small N` (4) minimum features for EM
  - Debug: `--debug-amb FILE`
- Cell-derived M_min (call_features)
  - `--mmin-from-cells`
  - `--mmin-cells-method {otsu|quantile|model3}`; params: `--mmin-qcells`, `--mmin-cap`
  - Model3: `--mmin3-max-iters`, `--mmin3-tol`, `--mmin3-update-disp`, `--mmin3-init {quantiles|kmeans}`, `--mmin3-floor`
- Feature handling
  - `--process-features`: use `features.txt` and treat all as allowed (suppresses enrichment stats); otherwise use `allowed_features.tsv`.

Outputs
- `<PFX>.assignments.tsv`, `<PFX>.doublets.txt`, `<PFX>.ambiguous.txt`, `<PFX>.unassignable.txt`, `<PFX>.missing_cells.txt`.
- In EM mode, optional `--debug-amb` CSV for ambiguous reasons.

Quick Recipes
- FLEX hashing
  - `./bin/flex_demux_mtx --mtx-dir HTO/matrix --out-prefix results/hto --tau 0.8 --delta 0.4 --gamma 0.9 --floor 12`
- Perturb-seq with EM (fixed dispersion)
  - `./bin/call_features --mtx-dir Guides/matrix --out-prefix results/perturb_EM --use-em --em-fixed-disp --k-min 4 --tau-pos 0.95 --tau-pos2 0.90 --gamma-min 0.8 --gamma-min-cand 0.85`
- Lineage with simple ratio
  - `./bin/call_features --mtx-dir LB/matrix --out-prefix results/larry --simple-assign --min-count 4 --min-ratio 2 --m-min-fixed 4`

Scripts (bin-aware)
- `scripts/runCall.sh`: wrapper for `./bin/flex_demux_mtx`.
- `scripts/runCallEm.sh`: wrapper for `./bin/call_features` with EM defaults.
- `scripts/runCallPerturb.sh`: perturbation-flavoured wrapper for `./bin/call_features`.
- `scripts/testSimpleCall.sh`: exercise simple heuristic mode.
- `scripts/testCallPerturb.sh`: multi-method tester (flex/em variants, simple, em-cells).
- `scripts/testLarry.sh`: lineage-oriented variant of the tester.
- `scripts/compareAllMethods.sh`: run all methods and aggregate results.
- `scripts/call_a_directory.sh`: batch over many subdirectories via the tester.

Developer Map
- `src/flex_assign.c`: FLEX classifier (no EM)
- `src/call_features.c`: streaming demux (FLEX, Simple, EM) + cell-derived M_min
- `src/per_guide_em.c`: EM engine used by `call_features`
- `src/compute_Mmin_from_cells.c`: cell-derived M_min estimators
