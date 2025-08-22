#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# testAssignPerturb.sh – comprehensive test script for flex_demux_mtx methods
#
# Usage:
#   ./testAssignPerturb.sh flex                    # Test FLEX method (original binomial)
#   ./testAssignPerturb.sh em                      # Test EM method 
#   ./testAssignPerturb.sh em-fixed                # Test EM with fixed dispersion
#   ./testAssignPerturb.sh em-noncells             # Test EM using non-cells for ambient
#   ./testAssignPerturb.sh simple                  # Test Simple method (cell-derived M_min)
#   ./testAssignPerturb.sh otsu                    # Test Otsu method (cell-derived M_min)
#   ./testAssignPerturb.sh quantile                # Test Quantile method (cell-derived M_min)
#   ./testAssignPerturb.sh model3                  # Test Model3 method (cell-derived M_min)
#   ./testAssignPerturb.sh                         # Show usage and exit
# ---------------------------------------------------------------------------

# Check for method argument
if [[ $# -eq 0 ]] || [[ ! "$1" =~ ^(flex|em|em-fixed|em-noncells|simple|otsu|quantile|model3)$ ]]; then
    echo "Usage: $0 {flex|em|em-fixed|em-noncells|simple|otsu|quantile|model3}"
    echo ""
    echo "Test different demultiplexing methods:"
    echo ""
    echo "Original methods:"
    echo "  flex        - Original FLEX method (binomial vs ambient, no EM)"
    echo "  em          - EM method with dispersion updates"
    echo "  em-fixed    - EM method with fixed dispersion"
    echo "  em-noncells - EM method using non-cells for ambient rate estimates"
    echo ""
    echo "Cell-derived M_min methods (all use EM-fixed):"
    echo "  simple      - Simple ratio-based method (ratio=2, min_counts=4)"
    echo "  otsu        - Histogram-based threshold using Otsu's method"
    echo "  quantile    - Simple quantile-based approach" 
    echo "  model3      - 3-component NB mixture with EM algorithm"
    echo ""
    exit 1
fi

METHOD="$1"
SC_ID="30_KO_ES"

# REQUIRED paths (same as runAssignPerturb.sh)
MTX_DIR="/storage/gene_features/${SC_ID}"
STARSOLO_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/Solo.out/Gene"
OUT_PREFIX="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/features_test_${METHOD}"

# OPTIONAL Flex parameters (unchanged defaults)
TAU="0.8"          # --tau
DELTA="0.4"        # --delta
GAMMA="0.9"        # --gamma
ALPHA="1e-4"       # --alpha
FLOOR="12"         # --floor
AMBIENT_Q="0.999"  # --ambient-q
USE_FDR=true       # set to false to append --no-fdr
MMIN_FIXED="3"    # --m-min-fixed
CELL_LIST="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/non_empty_barcodes.txt"
#CELL_LIST="/storage/gene_features/${SC_ID}/common_barcodes.txt"

# -------- EM-specific knobs --------
TAU_POS="0.85"     # positive-posterior threshold per guide
TAU_POS2="0.80"    # positive-posterior threshold per guide for top2
GAMMA_MIN="0.70"   # min (c1+c2)/total dominance to consider doublet
GAMMA_MIN_CAND=".85" # min candidate-only dominance to consider doublet
K_MIN="3"          # min counts for a guide to be considered present
K_MIN2="2"         # min counts for a guide to be considered present in 2+ cells
DOUBLET_BALANCE="1" # require 20–80% balance for top2 (default 1)
# ------------------------------------

# Cell-derived M_min parameters
MMIN_QCELLS="0.60"     # quantile for quantile method
MMIN_CAP="1000000"     # safety cap

# Simple method parameters
MMIN_SIMPLE_RATIO="2"      # ratio for simple method
MMIN_SIMPLE_MIN_COUNTS="3" # minimum counts for simple method

# Model3-specific parameters
MMIN3_MAX_ITERS="50"   # EM iterations
MMIN3_TOL="1e-6"       # convergence tolerance
MMIN3_UPDATE_DISP="0"  # 0=fixed dispersion, 1=update dispersion
MMIN3_INIT="quantiles" # quantiles|kmeans

# Build the base command (common parameters)
CMD=(./call_features
     --process-features
     --m-min-fixed  "$MMIN_FIXED"
     --mtx-dir      "$MTX_DIR"
     --starsolo-dir "$STARSOLO_DIR"
     --out-prefix   "$OUT_PREFIX"
     --tau          "$TAU"
     --delta        "$DELTA"
     --gamma        "$GAMMA"
     --alpha        "$ALPHA"
     --floor        "$FLOOR"
     --ambient-q    "$AMBIENT_Q")

# Add method-specific parameters
case "$METHOD" in
    "flex")
        echo "=== Testing FLEX method (original binomial) ==="
        echo "Uses binomial test vs ambient rates, no EM"
        # No --use-em flag, uses traditional FLEX approach
        # Remove --process-features for pure FLEX
        CMD=(./call_features
             --process-features
             --mtx-dir      "$MTX_DIR"
             --starsolo-dir "$STARSOLO_DIR"
             --out-prefix   "$OUT_PREFIX"
             --tau          "$TAU"
             --delta        "$DELTA"
             --gamma        "$GAMMA"
             --alpha        "$ALPHA"
             --floor        "$FLOOR"
             --ambient-q    "$AMBIENT_Q")
        ;;
    "em")
        echo "=== Testing EM method ==="
        echo "Uses per-guide EM with dispersion updates"
        CMD+=(--use-em
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE")
        ;;
    "em-fixed")
        echo "=== Testing EM method with fixed dispersion ==="
        echo "Uses per-guide EM without dispersion updates"
        CMD+=(--use-em
              --em-fixed-disp
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE")
        ;;
    "em-noncells")
        echo "=== Testing EM method using non-cells for ambient ==="
        echo "Uses per-guide EM with ambient rates estimated from non-cell barcodes"
        CMD+=(--use-em
              --em-fixed-disp
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE")
        echo "Note: This method estimates ambient rates from non-cell droplets"
        ;;
    "simple")
        echo "=== Testing SIMPLE method (cell-derived M_min) ==="
        echo "Uses ratio-based threshold: ratio=${MMIN_SIMPLE_RATIO}, min_counts=${MMIN_SIMPLE_MIN_COUNTS}"
        CMD+=(--simple-assign
              --min-ratio "$MMIN_SIMPLE_RATIO"
              --min-count "$MMIN_SIMPLE_MIN_COUNTS")
        ;;
    "otsu")
        echo "=== Testing OTSU method (cell-derived M_min) ==="
        echo "Uses histogram of log1p(cell_totals) with between-class variance maximization"
        CMD+=(--use-em
              --em-fixed-disp
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE"
              --mmin-from-cells
              --mmin-cells-method "$METHOD"
              --mmin-cap     "$MMIN_CAP")
        ;;
    "quantile")
        echo "=== Testing QUANTILE method (cell-derived M_min) ==="
        echo "Uses ${MMIN_QCELLS} quantile of cell totals"
        CMD+=(--use-em
              --em-fixed-disp
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE"
              --mmin-from-cells
              --mmin-cells-method "$METHOD"
              --mmin-cap     "$MMIN_CAP"
              --mmin-qcells  "$MMIN_QCELLS")
        ;;
    "model3")
        echo "=== Testing MODEL3 method (cell-derived M_min) ==="
        echo "3-component NB mixture (A/S/D) with EM algorithm"
        CMD+=(--use-em
              --em-fixed-disp
              --tau-pos      "$TAU_POS"
              --tau-pos2     "$TAU_POS2"
              --gamma-min    "$GAMMA_MIN"
              --gamma-min-cand "$GAMMA_MIN_CAND"
              --k-min        "$K_MIN"
              --k-min2       "$K_MIN2"
              --doublet-balance "$DOUBLET_BALANCE"
              --mmin-from-cells
              --mmin-cells-method "$METHOD"
              --mmin-cap     "$MMIN_CAP"
              --mmin3-max-iters "$MMIN3_MAX_ITERS"
              --mmin3-tol       "$MMIN3_TOL"
              --mmin3-update-disp "$MMIN3_UPDATE_DISP"
              --mmin3-init      "$MMIN3_INIT")
        ;;
esac

# Forward optional switches
if [[ -n "$CELL_LIST" ]]; then
    CMD+=(--cell-list "$CELL_LIST")
fi
if [[ "$USE_FDR" == false ]]; then
    CMD+=(--no-fdr)
fi

# Create output directory if needed
mkdir -p "$(dirname "$OUT_PREFIX")"

echo ""
echo "Method: $METHOD"
echo "Output prefix: $OUT_PREFIX"
echo ""
echo "Running: ${CMD[*]}"
echo ""

# Run the command
"${CMD[@]}"

# Check exit status
EXIT_CODE=$?
echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ SUCCESS: $METHOD method completed successfully"
    echo ""
    echo "Output files:"
    ls -la "${OUT_PREFIX}".*
    echo ""
    case "$METHOD" in
        "flex")
            echo "Check the log output above for FLEX binomial test results!"
            ;;
        "em"|"em-fixed"|"em-noncells")
            echo "Check the log output above for EM convergence and per-guide fit results!"
            echo "Look for .per_guide_fit.tsv and .per_guide_poscounts.tsv files."
            ;;
        "simple"|"otsu"|"quantile"|"model3")
            echo "Check the log output above for cell-derived M_min results!"
            ;;
    esac
else
    echo "❌ ERROR: $METHOD method failed with exit code $EXIT_CODE"
fi

echo ""
echo "=== Run completed for method: $METHOD ==="
