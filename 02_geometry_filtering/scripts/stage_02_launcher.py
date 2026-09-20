#!/usr/bin/env python3
"""
Runs stage 02 end to end:

1. Geometry-filters everything stage 01 handed over that has not been filtered yet;
2. Hands the passed structures to stage 03.

Usage:
    ./stage_02_launcher.py                      # everything outstanding
    ./stage_02_launcher.py --experiment NAME    # just this one
    ./stage_02_launcher.py --no-transfer        # filter only
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent
HELPING_SCRIPTS_DIR = SCRIPT_DIR / "helping_scripts"
# job_paths.py and the shared archive/transfer machinery live at the repo root.
COMMON_DIR = STAGE_ROOT.parent / "common"

sys.path.insert(0, str(COMMON_DIR))
sys.path.insert(0, str(HELPING_SCRIPTS_DIR))

from archives import ArchiveError  # noqa: E402
from geometry_filter import GeometryError, run_geometry_filter  # noqa: E402
from transfer_to_stage03 import TransferError, run_transfer  # noqa: E402


def _log(*parts: object) -> None:
    print(*parts, flush=True)


def filter_structures(stage: Path, experiment: Optional[str]) -> bool:
    """Align and sort everything not already in a results table."""
    _log("[geometry] sorting structures into passed/rejected...")
    try:
        report = run_geometry_filter(stage, experiment)
    except (OSError, ArchiveError, GeometryError) as exc:
        _log(f"[geometry] FAILED: {exc}")
        _log("[geometry] nothing has been lost -- fix the cause and re-run "
             "scripts/helping_scripts/geometry_filter.py")
        return False

    _log(f"[geometry] {report.summary()}")
    for table in report.tables:
        _log(f"[geometry] results -> {table}")
    skipped = report.counts.get("SKIPPED", 0)
    if skipped:
        _log(f"[geometry] {skipped} structure(s) skipped: their stage-01 seed has "
             f"no ligand chain, so there is nothing to check the design against")
    if not report.ok:
        _log(f"[geometry] {report.counts.get('ERROR', 0)} structure(s) could not be "
             f"evaluated -- see the skip_reason column")
    return report.ok


def transfer_structures(stage: Path, experiment: Optional[str]) -> bool:
    """Hand every passed structure not already in stage 03 across to it.

    Runs whatever the filter reported: a SKIPPED or ERROR structure never
    reaches passed/, and the archive merge is atomic, so a filter problem
    cannot leave a half-written source archive for this to read.
    """
    _log("[transfer] handing passed structures to stage 03...")
    try:
        report = run_transfer(stage, experiment)
    except (OSError, ArchiveError, TransferError) as exc:
        _log(f"[transfer] FAILED: {exc}")
        _log("[transfer] filtering is unaffected -- fix the cause and re-run "
             "scripts/helping_scripts/transfer_to_stage03.py")
        return False

    _log(f"[transfer] {report.summary()}")
    if report.incomplete:
        preview = ", ".join(report.incomplete[:5])
        more = f" (+{len(report.incomplete) - 5} more)" if len(report.incomplete) > 5 else ""
        _log(f"[transfer] {len(report.incomplete)} passed structure(s) could not be "
             f"transferred, missing half their file pair: {preview}{more}")
    return report.ok


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--stage", type=Path, default=STAGE_ROOT,
        help=f"stage root directory (default: {STAGE_ROOT})",
    )
    parser.add_argument(
        "--experiment", default=None,
        help="only process this one experiment (default: every outstanding one)",
    )
    parser.add_argument(
        "--no-filter", action="store_true",
        help="skip the geometry filter; transfer whatever already passed",
    )
    parser.add_argument(
        "--no-transfer", action="store_true",
        help="filter only; do not hand anything to stage 03",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    stage = args.stage.resolve()

    filtered_ok = True
    if not args.no_filter:
        filtered_ok = filter_structures(stage, args.experiment)

    transferred_ok = True
    if not args.no_transfer:
        transferred_ok = transfer_structures(stage, args.experiment)

    return 0 if filtered_ok and transferred_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
