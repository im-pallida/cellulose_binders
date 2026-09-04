#!/usr/bin/env bash
# Stage 04 (04_alphafold), script 2: machine-aware dispatcher for
# AlphaFold3.
#
#   --ws       -> helping_scripts/15082026_02_run_ws.sh       (bash, direct, sequential)
#   --cluster  -> helping_scripts/15082026_02_run_cluster.sh  (bash, direct -- submits sbatch jobs itself)
#   (no arg)   -> auto-detect: sbatch on PATH => cluster, else ws
#
# Unlike the ProteinMPNN launcher (one batch job doing the whole
# pipeline), AF3 here means one independent SLURM job PER input json
# (monomer + dimer, 2 per selected sequence -- see
# 15082026_01_af3_json_builder.py). So on the cluster this dispatcher
# does NOT submit anything as a job itself -- it runs
# 15082026_02_run_cluster.sh directly on the login node, which submits
# individual AF3 jobs (15082026_02_run_one_job.sbatch, one per json) and
# throttles them to your account's confirmed 4-concurrent-job cap,
# starting with a one-job canary before queueing the rest.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPING_SCRIPTS="$SCRIPT_DIR/helping_scripts"
WS_SCRIPT="$HELPING_SCRIPTS/15082026_02_run_ws.sh"
CLUSTER_SCRIPT="$HELPING_SCRIPTS/15082026_02_run_cluster.sh"

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
    if [ ! -f "$CLUSTER_SCRIPT" ]; then
        echo "ERROR: $CLUSTER_SCRIPT not found"
        exit 1
    fi
    echo "[dispatch] sbatch found on PATH -- running the cluster submitter directly: $CLUSTER_SCRIPT"
    echo "           This submits one AF3 job per pending json (canary first, then throttled to"
    echo "           4 concurrent) and polls squeue until everything pending is submitted -- it"
    echo "           does no compute itself, so it's safe to leave running in tmux/screen."
    bash "$CLUSTER_SCRIPT"
else
    if [ ! -f "$WS_SCRIPT" ]; then
        echo "ERROR: $WS_SCRIPT not found"
        exit 1
    fi
    echo "[dispatch] no sbatch on PATH -- running workstation launcher directly: $WS_SCRIPT"
    bash "$WS_SCRIPT"
fi
