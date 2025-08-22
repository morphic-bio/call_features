#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# runAssign.sh – thin wrapper around flex_demux_mtx
#
# Edit the variables below, then run:
#   ./runAssign.sh
#
# All long-option flag names match those accepted by the program.
# ---------------------------------------------------------------------------

# REQUIRED paths
MTX_DIR="/storage/JAX_features/SC2300771"                  # directory containing matrix.mtx etc.
STARSOLO_DIR="/storage/scRNAseq_output/Alignments/SC2300771/star/Solo.out/Gene"       # STARsolo Gene directory (raw/, filtered/)
OUT_PREFIX="/storage/JAX_features/results/assignments"   # output file prefix (dir will be created)

# OPTIONAL tuning knobs (comment out to keep program defaults)
TAU="0.8"          # --tau
DELTA="0.4"        # --delta
GAMMA="0.9"        # --gamma
ALPHA="1e-4"       # --alpha
FLOOR="12"         # --floor
AMBIENT_Q="0.999"  # --ambient-q
USE_FDR=true       # set to false to append --no-fdr

CELL_LIST="/storage/scRNAseq_output/Alignments/SC2300771/star/non_empty_barcodes.txt"      # optional path to a file with barcodes to use as good cells

# ---------------------------------------------------------------------------
# Build the command
CMD=(./flex_demux_mtx
      --mtx-dir      "$MTX_DIR"
      --starsolo-dir "$STARSOLO_DIR"
      --out-prefix   "$OUT_PREFIX"
      --tau          "$TAU"
      --delta        "$DELTA"
      --gamma        "$GAMMA"
      --alpha        "$ALPHA"
      --floor        "$FLOOR"
      --ambient-q    "$AMBIENT_Q")

# Forward optional switches
if [[ -n "$CELL_LIST" ]]; then
    CMD+=(--cell-list "$CELL_LIST")
fi
if [[ "$USE_FDR" == false ]]; then
    CMD+=(--no-fdr)
fi

# Create output directory if needed
mkdir -p "$(dirname "$OUT_PREFIX")"

# Echo and run
echo "Running: ${CMD[*]}"
"${CMD[@]}"
