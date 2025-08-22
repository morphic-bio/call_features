#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# testSimpleAssign.sh – test script for simple classification mode
#
# Usage:
#   ./testSimpleAssign.sh [min_count] [min_ratio]   # Test simple mode with custom parameters
#   ./testSimpleAssign.sh                           # Test with defaults (min_count=2, min_ratio=2.0)
# ---------------------------------------------------------------------------

# Parse arguments or use defaults
MIN_COUNT=${1:-2}      # default: 2
MIN_RATIO=${2:-2.0}    # default: 2.0

SC_ID="30_KO_ES"

# REQUIRED paths (same as runAssignPerturb.sh)
MTX_DIR="/storage/gene_features/${SC_ID}"
STARSOLO_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/Solo.out/Gene"
OUT_PREFIX="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/features_simple_${MIN_COUNT}_${MIN_RATIO}"

# OPTIONAL Flex parameters (unchanged defaults)
TAU="0.8"          # --tau (not used in simple mode)
DELTA="0.4"        # --delta (not used in simple mode)
GAMMA="0.9"        # --gamma (not used in simple mode)
ALPHA="1e-4"       # --alpha (not used in simple mode)
FLOOR="12"         # --floor (not used in simple mode)
AMBIENT_Q="0.999"  # --ambient-q
USE_FDR=true       # set to false to append --no-fdr

CELL_LIST="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/non_empty_barcodes.txt"
#CELL_LIST="/storage/gene_features/${SC_ID}/common_barcodes.txt"

echo "=== Testing SIMPLE CLASSIFICATION ==="
echo "Simple rules: min_count=${MIN_COUNT}, min_ratio=${MIN_RATIO}"
echo "Logic: top1_count >= min_count AND (top2_count == 0 OR top1/top2 >= min_ratio)"
echo ""

# Build the command - much simpler than other modes
CMD=(./call_features
     --process-features
     --mtx-dir      "$MTX_DIR"
     --starsolo-dir "$STARSOLO_DIR"
     --out-prefix   "$OUT_PREFIX"
     --ambient-q    "$AMBIENT_Q"
     --simple-assign
     --min-count    "$MIN_COUNT"
     --min-ratio    "$MIN_RATIO")

# Forward optional switches
if [[ -n "$CELL_LIST" ]]; then
    CMD+=(--cell-list "$CELL_LIST")
fi
if [[ "$USE_FDR" == false ]]; then
    CMD+=(--no-fdr)
fi

# Create output directory if needed
mkdir -p "$(dirname "$OUT_PREFIX")"

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
    echo "✅ SUCCESS: Simple classification completed successfully"
    echo ""
    echo "Output files:"
    ls -la "${OUT_PREFIX}".*
    echo ""
    echo "Check the log output above for classification results!"
else
    echo "❌ ERROR: Simple classification failed with exit code $EXIT_CODE"
fi

echo ""
echo "=== Run completed for simple classification ==="
echo "Parameters used: min_count=${MIN_COUNT}, min_ratio=${MIN_RATIO}"
