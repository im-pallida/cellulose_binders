#!/bin/bash
set -euo pipefail
export EXP=/home/karina/1cbh_based_linear_binder_cellulose/02_rfd3_linear_symmetry_generation
export RFD3_ENV=/home/karina/miniconda3/envs/foundry
export RFD3_PYTHON="$RFD3_ENV/bin/python"
export RFD3_BIN="$RFD3_ENV/bin/rfd3"
export RFD3_CKPT=/home/karina/.foundry/checkpoints/rfd3_latest.ckpt
export RFD3_OG_PARENT="$RFD3_ENV/lib/python3.12/site-packages"
export RFD3_OVERLAY_PARENT="$EXP/overlay/rfd3_t3_overlay"
# Overlay first, original installed RFD3 second.
export PYTHONPATH="$RFD3_OVERLAY_PARENT:$RFD3_OG_PARENT:${PYTHONPATH:-}"
echo "Using RFD3 overlay:"
echo "  overlay: $RFD3_OVERLAY_PARENT"
echo "  original: $RFD3_OG_PARENT"
echo "  python: $RFD3_PYTHON"
echo "  rfd3: $RFD3_BIN"
