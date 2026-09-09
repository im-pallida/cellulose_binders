#!/usr/bin/env python3
"""
Workstation dispatch: run every job in a run list serially, in-process.
 
Imported and called by the launcher; not useful on its own.
"""
from __future__ import annotations
 
import sys
from pathlib import Path
 
# job_paths.py lives one level up, in scripts/. Inserted here explicitly so
# this module does not depend on another import having done it first.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import job_paths as jp  # noqa: E402
 
from run_one_job import DispatchReport, JobStatus, run_one_job  # noqa: E402
 
 
def _log(*parts: object) -> None:
    print(*parts, flush=True)
 
 
def dispatch_workstation(experiment: str, run_list: Path, stage: Path) -> DispatchReport:
    rows = jp.read_run_list(run_list)
    report = DispatchReport(rows=list(rows))
    _log(f"[dispatch] mode=workstation -> {len(rows)} job(s), serial, in-process")
 
    for index, (json_rel, global_seq) in enumerate(rows, start=1):
        _log(f"[dispatch] ({index}/{len(rows)}) {json_rel} seq={global_seq}")
        # run_one_job() converts every failure into a FAILED result, so one bad
        # job cannot abandon the rest of the queue.
        result = run_one_job(experiment, json_rel, global_seq, stage=stage)
        if result.status is JobStatus.FAILED:
            _log(f"[fail] {json_rel} seq={global_seq}: {result.message}")
            report.failed.append((json_rel, global_seq))
        _log()
 
    _log(f"[dispatch] {report.summary()}")
    for json_rel, global_seq in report.failed:
        _log(f"  {json_rel} seq={global_seq}")
    return report
