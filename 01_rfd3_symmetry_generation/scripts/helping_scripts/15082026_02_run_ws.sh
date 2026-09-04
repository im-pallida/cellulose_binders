#!/usr/bin/env bash
# Serial runner - reads the ephemeral run list built by 02_run_rfd3_sizes.sh
# (columns: order, size, global_seq) and calls 02_run_rfd3_one_job.sh once
# per line.
#
# Usage:
#   ./02_run_rfd3_sizes_workstation.sh <experiment> <run_list_path>

set -euo pipefail

STAGE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXPERIMENT="$1"
RUN_LIST="$2"

if [ ! -f "$RUN_LIST" ]; then
    echo "ERROR: missing run list: $RUN_LIST"
    exit 1
fi

while IFS=$'\t' read -r order size global_seq; do
    [ -n "$order" ] || continue
    "$STAGE/scripts/helping_scripts/15082026_02_run_one_job.sh" "$EXPERIMENT" "$order" "$size" "$global_seq" || echo "[fail] $order $size $global_seq"
    echo
done < "$RUN_LIST"

echo "Done."
