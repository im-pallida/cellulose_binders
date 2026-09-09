#!/usr/bin/env python3
"""
Cluster dispatch: submit a run list via sbatch, canary job first.

Each sbatch job runs run_one_job.py as a script on a compute node. The first
row is submitted blocking and verified before the rest are queued, so a broken
environment costs one job instead of the whole run list.

Imported and called by the launcher; not useful on its own.
"""
from __future__ import annotations

import getpass
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence

# job_paths.py lives one level up, in scripts/. Inserted here explicitly so
# this module does not depend on another import having done it first.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import job_paths as jp  # noqa: E402

from run_one_job import DispatchReport  # noqa: E402

SBATCH_SCRIPT = Path(__file__).resolve().parent / "run_cluster.sbatch"
POLL_INTERVAL_S = 30

# squeue can fail transiently while the controller is busy. Tolerate a few
# failures in a row, then give up loudly rather than silently concluding that
# every job has finished.
MAX_CONSECUTIVE_SQUEUE_ERRORS = 5


def _log(*parts: object) -> None:
    print(*parts, flush=True)


def _submit(experiment: str, json_rel: str, global_seq: int, wait: bool, stage: Path) -> Optional[str]:
    """Submit one job. wait=True blocks until it finishes; wait=False returns
    the job ID."""
    log_dir = stage / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
 
    cmd: List[str] = [
        "sbatch",
        f"--chdir={stage}",
        f"--output={log_dir / 'rfd3_%j.out'}",
        f"--error={log_dir / 'rfd3_%j.err'}",
    ]
    if wait:
        cmd.append("--wait")
    cmd += [str(SBATCH_SCRIPT), experiment, json_rel, str(global_seq)]
 
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        what = "job failed on the compute node" if wait else "sbatch could not submit the job"
        raise SystemExit(
            f"ERROR: {what} for {json_rel} seq={global_seq} (exit {proc.returncode}):\n"
            f"{proc.stderr.strip() or proc.stdout.strip()}"
        )
    if wait:
        return None
 
    match = re.search(r"Submitted batch job (\d+)", proc.stdout)
    if not match:
        raise SystemExit(f"ERROR: could not parse job id from sbatch output: {proc.stdout!r}")
    return match.group(1)
 
 
def _poll_until_done(job_ids: Sequence[str]) -> None:
    """Wait for every submitted job to leave the queue."""
    remaining = set(job_ids)
    user = getpass.getuser()
    consecutive_errors = 0
 
    while remaining:
        proc = subprocess.run(
            ["squeue", "--noheader", "--user", user, "--format=%i"],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            consecutive_errors += 1
            if consecutive_errors >= MAX_CONSECUTIVE_SQUEUE_ERRORS:
                raise SystemExit(
                    f"ERROR: squeue failed {consecutive_errors} times in a row; "
                    f"cannot tell which of {len(remaining)} job(s) are still running.\n"
                    f"{proc.stderr.strip()}"
                )
            _log(f"[cluster] WARNING: squeue failed ({proc.returncode}), retrying: "
                 f"{proc.stderr.strip()}")
            time.sleep(POLL_INTERVAL_S)
            continue
 
        consecutive_errors = 0
        still_queued = {line.strip() for line in proc.stdout.split() if line.strip()}
        for job_id in sorted(remaining - still_queued):
            _log(f"[cluster] job {job_id} finished")
        remaining &= still_queued
        if remaining:
            time.sleep(POLL_INTERVAL_S)
 
 
def _missing_outputs(
    experiment: str, rows: Sequence[jp.RunRow], stage: Path
) -> List[jp.RunRow]:
    return [
        (json_rel, global_seq)
        for json_rel, global_seq in rows
        if not jp.expected_raw_cif(stage, experiment, json_rel, global_seq).exists()
    ]
 
 
def dispatch_cluster(experiment: str, run_list: Path, stage: Path) -> DispatchReport:
    if not SBATCH_SCRIPT.is_file():
        raise SystemExit(f"ERROR: sbatch script not found: {SBATCH_SCRIPT}")
 
    rows = jp.read_run_list(run_list)
    if not rows:
        raise SystemExit(f"ERROR: run list {run_list} is empty, nothing to dispatch")
    report = DispatchReport(rows=list(rows))
    _log(f"[dispatch] mode=cluster -> {len(rows)} job(s) via sbatch, canary first")
 
    # Canary: one blocking job, verified, before committing the rest to the queue.
    canary_json_rel, canary_seq = rows[0]
    _log(f"[canary] submitting {canary_json_rel} seq={canary_seq} (sbatch --wait)")
    _submit(experiment, canary_json_rel, canary_seq, wait=True, stage=stage)
    canary_cif = jp.expected_raw_cif(stage, experiment, canary_json_rel, canary_seq)
    if not canary_cif.exists():
        raise SystemExit(
            f"ERROR: canary finished but its output is missing: {canary_cif}\n"
            f"  check the sbatch logs in {stage / 'logs'} before resubmitting."
        )
    _log(f"[canary] ok -> {canary_cif}")
 
    remaining = rows[1:]
    if remaining:
        job_ids: List[str] = []
        for json_rel, global_seq in remaining:
            job_id = _submit(experiment, json_rel, global_seq, wait=False, stage=stage)
            job_ids.append(str(job_id))
            _log(f"[dispatch] submitted {json_rel} seq={global_seq} -> job {job_id}")
 
        _log(f"[dispatch] waiting for {len(job_ids)} job(s)...")
        _poll_until_done(job_ids)
    else:
        _log("[dispatch] canary was the whole run list")
 
    _log("[dispatch] verifying outputs...")
    report.failed = _missing_outputs(experiment, rows, stage)
    _log(f"[dispatch] {report.summary()}")
    for json_rel, global_seq in report.failed:
        missing = jp.expected_raw_cif(stage, experiment, json_rel, global_seq)
        _log(f"  {json_rel} seq={global_seq}: missing {missing}")
    return report
