#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# runAssignPerturb.sh – wrapper that runs the perturbation features mode (flex_demux_mtx)
#
# Edit the variables below, then execute:
#   ./runAssignPerturb.sh
# ---------------------------------------------------------------------------
SC_ID="30_KO_ES"


# REQUIRED paths (same as runAssignEM.sh)
MTX_DIR="/storage/gene_features/${SC_ID}"
STARSOLO_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/Solo.out/Gene"
OUT_PREFIX="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/features"

# OPTIONAL Flex parameters (unchanged defaults)
TAU="0.8"          # --tau
DELTA="0.4"        # --delta
GAMMA="0.9"        # --gamma
ALPHA="1e-4"       # --alpha
FLOOR="12"         # --floor
AMBIENT_Q="0.1"  # --ambient-q
USE_FDR=true       # set to false to append --no-fdr

CELL_LIST="/mnt/pikachu/storage/MSK-output-2/Alignments/${SC_ID}/star/non_empty_barcodes.txt"

# -------- EM-specific knobs (streams binary only) --------
TAU_POS="0.85"     # positive-posterior threshold per guide
TAU_POS2="0.80"    # positive-posterior threshold per guide for top2
GAMMA_MIN="0.70"   # min (c1+c2)/total dominance to consider doublet
GAMMA_MIN_CAND=".85" # min candidate-only dominance to consider doublet
K_MIN="2"          # min counts for a guide to be considered present
K_MIN2="2"         # min counts for a guide to be considered present in 2+ cells
DOUBLET_BALANCE="1" # require 20–80% balance for top2 (default 1)
# ---------------------------------------------------------

# Build the command
CMD=(./call_features
     --csc
     --use-em
     --em-fixed-disp
     --mtx-dir      "$MTX_DIR"
     --starsolo-dir "$STARSOLO_DIR"
     --out-prefix   "$OUT_PREFIX"
     --tau          "$TAU"
     --delta        "$DELTA"
     --gamma        "$GAMMA"
     --alpha        "$ALPHA"
     --floor        "$FLOOR"
     --ambient-q    "$AMBIENT_Q"
     --tau-pos      "$TAU_POS"
     --tau-pos2     "$TAU_POS2"
     --gamma-min    "$GAMMA_MIN"
     --gamma-min-cand "$GAMMA_MIN_CAND"
     --k-min        "$K_MIN"
     --k-min2      "$K_MIN2"
     --doublet-balance "$DOUBLET_BALANCE")

# Forward optional switches
if [[ -n "$CELL_LIST" ]]; then
    CMD+=(--cell-list "$CELL_LIST")
fi
if [[ "$USE_FDR" == false ]]; then
    CMD+=(--no-fdr)
fi

# Create output directory if needed
mkdir -p "$(dirname "$OUT_PREFIX")"

echo "Running: ${CMD[*]}"
"${CMD[@]}"
