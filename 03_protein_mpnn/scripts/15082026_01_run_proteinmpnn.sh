#!/usr/bin/env bash
# Stage 03, script 1: machine-aware dispatcher for ProteinMPNN.
#
# This script does NOT run ProteinMPNN itself -- it only decides which
# machine we're on and hands off to the script that actually knows how
# to run it there:
#   --ws       -> helping_scripts/15082026_01_run_ws.sh       (bash, direct)
#   --cluster  -> helping_scripts/15082026_01_run_cluster.sbatch (sbatch)
#   (no arg)   -> auto-detect: sbatch on PATH => cluster, else ws
#
# protein_mpnn_run.py is only ever invoked from inside those two helping
# scripts, which is also where MPNN_ROOT/conda-env/paths are resolved.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$(cd "$SCRIPT_DIR/.." && pwd)"
HELPING_SCRIPTS="$SCRIPT_DIR/helping_scripts"
WS_SCRIPT="$HELPING_SCRIPTS/15082026_01_run_ws.sh"
CLUSTER_SBATCH="$HELPING_SCRIPTS/15082026_01_run_cluster.sbatch"

MODE="${1:-auto}"
case "$MODE" in
    --ws) TARGET=ws ;;
    --cluster) TARGET=cluster ;;
    auto|"")
        if command -v sbatch >/dev/null 2>&1; then
            TARGET=cluster
        else
            TARGET=ws
        fi
        ;;
    *)
        echo "ERROR: unrecognized argument: $MODE (expected --ws or --cluster)"
        exit 1
        ;;
esac

if [ "$TARGET" = "cluster" ]; then
    if [ ! -f "$CLUSTER_SBATCH" ]; then
        echo "ERROR: $CLUSTER_SBATCH not found"
        exit 1
    fi
    mkdir -p "$STAGE/logs"
    echo "[dispatch] sbatch found on PATH -- submitting cluster job: $CLUSTER_SBATCH"
    # SLURM copies the .sbatch file into a per-job spool directory before
    # running it, so ${BASH_SOURCE[0]} inside it does NOT resolve to this
    # project's real path -- STAGE/HELPING_SCRIPTS_DIR are passed explicitly
    # here instead, since THIS script's own path resolution (above) is
    # reliable (only the .sbatch file itself gets spool-copied by SLURM,
    # not this dispatcher).
    sbatch \
        --output="$STAGE/logs/mpnn_launcher_%j.out" \
        --error="$STAGE/logs/mpnn_launcher_%j.err" \
        --export="ALL,STAGE=$STAGE,HELPING_SCRIPTS_DIR=$HELPING_SCRIPTS" \
        "$CLUSTER_SBATCH"
else
    if [ ! -f "$WS_SCRIPT" ]; then
        echo "ERROR: $WS_SCRIPT not found"
        exit 1
    fi
    echo "[dispatch] no sbatch on PATH -- running workstation launcher directly: $WS_SCRIPT"
    bash "$WS_SCRIPT"
fi
