#!/usr/bin/env bash
# Stage-01 preflight. Run from the stage root:
#   ./scripts/preflight.sh
# Catches the things that would otherwise fail every job in a queue.
fail=0
say() { printf "  %-46s %s\n" "$1" "$2"; }

if grep -lc $'\r' scripts/*.py scripts/helping_scripts/*.py \
       scripts/helping_scripts/*.sbatch scripts/env/* 2>/dev/null | grep -q .; then
  say "line endings" "CRLF FOUND -> dos2unix those files"; fail=1
else say "line endings" "ok (LF)"; fi

if grep -q '_invoke_rfd3(stage, env, overlay_root, json_path, work_dir)' \
     scripts/helping_scripts/run_one_job.py 2>/dev/null; then
  say "_invoke_rfd3 call site" "ok"
else say "_invoke_rfd3 call site" "NOT FIXED -> every job will fail"; fail=1; fi

if grep -q 'INDEX=' scripts/helping_scripts/*.sbatch 2>/dev/null; then
  say "sbatch argument interface" "OLD 4-arg version -> replace it"; fail=1
else say "sbatch argument interface" "ok (3 args)"; fi

if python3 -c "
import sys; sys.path[:0]=['scripts','scripts/helping_scripts']
from job_paths import resolve_seed_path
from symmetry_check import format_translation_spec, compare_translation
import run_one_job, run_cluster, run_workstation, stage_01_launcher" 2>/dev/null; then
  say "all six modules in sync" "ok"
else
  say "all six modules in sync" "IMPORT FAILS:"; fail=1
  python3 -c "
import sys; sys.path[:0]=['scripts','scripts/helping_scripts']
import stage_01_launcher" 2>&1 | tail -2 | sed 's/^/      /'
fi

# Filename constants in the code vs files on disk. Renaming a file without
# updating its constant fails only at job time, so check it here instead.
python3 - <<'PY' || fail=1
import sys; sys.path[:0]=['scripts','scripts/helping_scripts']
from pathlib import Path
try:
    import run_one_job, run_cluster
except Exception:
    sys.exit(0)          # already reported by the import check above
stage = Path('.').resolve()
rows = [("env file      (run_one_job.ENV_FILE_REL)",  stage / run_one_job.ENV_FILE_REL),
        ("overlay root  (run_one_job.OVERLAY_REL)",   stage / run_one_job.OVERLAY_REL),
        ("sbatch script (run_cluster.SBATCH_SCRIPT)", run_cluster.SBATCH_SCRIPT)]
bad = 0
for label, p in rows:
    ok = p.exists(); bad += not ok
    shown = p.relative_to(stage) if p.is_relative_to(stage) else p
    print(f"  {label:46s} {'ok' if ok else 'MISSING -> ' + str(shown)}")
sys.exit(1 if bad else 0)
PY

echo
[ $fail -eq 0 ] && echo "  PREFLIGHT PASSED" || echo "  PREFLIGHT FAILED - fix the above first"
exit $fail
