#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# testCallPerturb.sh – comprehensive test script for flex_demux_mtx methods
#
# Usage:
#   ./testCallPerturb.sh flex                    # Test FLEX method (original binomial)
#   ./testCallPerturb.sh em                      # Test EM method 
#   ./testCallPerturb.sh em-fixed                # Test EM with fixed dispersion
#   ./testCallPerturb.sh em-noncells             # Test EM using non-cells for ambient
#   ./testCallPerturb.sh simple                  # Test Simple method (cell-derived M_min)
#   ./testCallPerturb.sh em-cells                # Test EM-cells method (cell-derived M_min)
#   ./testCallPerturb.sh                         # Show usage and exit
# ---------------------------------------------------------------------------

# Check for method argument
if [[ $# -eq 0 ]] || [[ ! "$1" =~ ^(flex|em|em-fixed|em-noncells|simple|em-cells)$ ]]; then
    echo "Usage: $0 {flex|em|em-fixed|em-noncells|simple|em-cells}"
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
    echo "  em-cells    - EM-cells method (cell-derived M_min)"
    echo ""
    exit 1
fi

METHOD="$1"
INPUT_DIR="$2"
SC_ID="$3"

# REQUIRED paths (same as runCallPerturb.sh)
#MTX_DIR="/storage/gene_features/${SC_ID}"
MTX_DIR="${INPUT_DIR}/${SC_ID}"
STARSOLO_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/Solo.out/Gene"
OUT_PREFIX="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/larry_features_${METHOD}/"
# make sure the output directory exists
outdir=$(dirname $OUT_PREFIX.add)
mkdir -p "$outdir"
# OPTIONAL Flex parameters (unchanged defaults)
TAU="0.8"          # --tau
DELTA="0.4"        # --delta
GAMMA="0.9"        # --gamma
ALPHA="1e-4"       # --alpha
FLOOR="12"         # --floor
AMBIENT_Q="0.999"  # --ambient-q
USE_FDR=true       # set to false to append --no-fdr
MMIN_FIXED="3"    # --m-min-fixed
THREADS="30"
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
MMIN_SIMPLE_MIN_COUNTS="4" # minimum counts for simple method

# Model3-specific parameters
MMIN3_MAX_ITERS="50"   # EM iterations
MMIN3_TOL="1e-6"       # convergence tolerance
MMIN3_UPDATE_DISP="0"  # 0=fixed dispersion, 1=update dispersion
MMIN3_INIT="quantiles" # quantiles|kmeans

# Build the base command (common parameters)
CMD=(./bin/call_features
     --process-features
     --threads "$THREADS"
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
        CMD=(./bin/call_features
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
    "em-cells")
        echo "=== Testing EM-cells method (cell-derived M_min) ==="
        echo "Uses EM-cells method to estimate M_min"
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
    outdir=$(dirname $OUT_PREFIX.add)
    ls -la "${outdir}"
    echo ""
    case "$METHOD" in
        "flex")
            echo "Check the log output above for FLEX binomial test results!"
            ;;
        "em"|"em-fixed"|"em-noncells")
            echo "Check the log output above for EM convergence and per-guide fit results!"
            echo "Look for .per_guide_fit.tsv and .per_guide_poscounts.tsv files."
            ;;
        "simple"|"em-cells")
            echo "Check the log output above for cell-derived M_min results!"
            ;;
    esac
else
    echo "❌ ERROR: $METHOD method failed with exit code $EXIT_CODE"
fi

echo ""
echo "=== Run completed for method: $METHOD ==="
