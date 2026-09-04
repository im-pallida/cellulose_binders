#!/usr/bin/env python3
"""
Sync: transfer passed proteins from stage-01's sorted_clean/ into
stage-02's inputs/, incrementally -- only proteins not already
transferred are copied on each run.

Source of truth for WHICH structures: whatever's currently sitting under
    01_rfd3_symmetry_generation/sorted_clean/<experiment>/passed/<group>.tar.gz
(rejected/ is never read -- only passed/ proteins become stage-02 inputs.)

Destination:
    02_geometry_filtering/inputs/<experiment>/<group>.tar.gz
Same <group>.tar.gz filename as the source (e.g. "htmedium.tar.gz"),
just relocated under 02_geometry_filtering/inputs/<experiment>/.

Why "incremental" instead of the exact-rebuild-with-removal pattern from
22072026_10_sync_passed_to_stage04.py:
sorted_clean/<experiment>/passed/<group>.tar.gz is NOT cumulative across
stage-01 runs -- 15082026_03_filtering.py's regroup_and_archive writes a
fresh <group>.tar.gz from only that run's newly-routed proteins each
time (see that script's own "known limitation" note). So this script
can't treat sorted_clean as a complete current snapshot the way the
stage04 sync script treats passed/ -- it has to remember, run over run,
which protein_ids it has already pulled in, and only add genuinely new
ones. That memory is job_runs/proteins_transferred_from_01_logs.tsv
(stage-01 side, keyed by protein_id).

Because sorted_clean's archives aren't cumulative, and destination
archives DO need to be cumulative (stage-02 needs to see every passed
protein ever produced, not just the latest stage-01 run's batch), this
script merges into an existing destination archive rather than
overwriting it: for each <group>.tar.gz that has new proteins, it reads
every member already in the destination archive (if any), adds the new
members, and atomically replaces the archive with the combined result.
This is the "decompress-append-recompress" idea from
pipeline_tracking.append_to_growing_archive, generalized to add several
new members at once instead of one file.

Only rows whose members are validated (must have BOTH <protein_id>.json
and <protein_id>.cif/.pdb in the source archive) are transferred; an
incomplete pair is skipped with a printed warning rather than silently
dropped or half-copied.

Usage:
    python 15082026_04_sync_passed_to_stage02.py [--experiment NAME]
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# CONSTANTS -- all relative to this script's own location, so the script
# works unmodified if the whole project folder is copied to another
# machine/user (no hardcoded "/home/karina/..." anywhere).
# Assumes this script lives in 01_rfd3_symmetry_generation/scripts/,
# same convention as 22072026_10_sync_passed_to_stage04.py (the
# producing stage owns the script that pushes its results forward).
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE01_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = STAGE01_ROOT.parent
STAGE02_ROOT = PROJECT_ROOT / "02_geometry_filtering"

SORTED_CLEAN_DIR = STAGE01_ROOT / "sorted_clean"
STAGE02_INPUTS_DIR = STAGE02_ROOT / "inputs"
TRANSFER_LOG_PATH = STAGE01_ROOT / "job_runs" / "proteins_transferred_from_01_logs.tsv"
TRANSFER_LOG_DELIMITER = "\t"

TRANSFER_LOG_FIELDS = [
    "protein_id",
    "experiment",
    "group",
    "source_archive",
    "transferred_at",
]


# ---------------------------------------------------------------------------
# Transfer log -- same load-whole-dict/write-whole-dict pattern as
# stage_01_filter.py's load_table/save_table, keyed on protein_id, but
# written as TSV (tab-delimited) under job_runs/ instead of a CSV under
# tables/, per your naming convention for run/job-tracking logs.
# ---------------------------------------------------------------------------

def load_transfer_log(log_path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    if log_path.exists() and log_path.stat().st_size > 0:
        with log_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=TRANSFER_LOG_DELIMITER)
            for row in reader:
                by_id[row["protein_id"]] = row
    return by_id


def save_transfer_log(by_id: dict[str, dict[str, str]], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    rows_sorted = sorted(by_id.values(), key=lambda row: row.get("protein_id", ""))
    with log_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=TRANSFER_LOG_FIELDS, delimiter=TRANSFER_LOG_DELIMITER, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow({field: row.get(field, "") for field in TRANSFER_LOG_FIELDS})


# ---------------------------------------------------------------------------
# Discovery + grouping (same shape as stage_01_filter.py's helpers)
# ---------------------------------------------------------------------------

def collect_dirs(root_dir: Path) -> list[Path]:
    return sorted(p for p in root_dir.iterdir() if p.is_dir())


def discover_passed_archives(sorted_clean_root: Path, specific_experiment: str | None = None) -> list[Path]:
    """Only passed/ archives -- rejected/ is never a transfer source."""
    if specific_experiment is not None:
        experiments = [sorted_clean_root / specific_experiment]
    else:
        experiments = collect_dirs(sorted_clean_root)

    archives: list[Path] = []
    for experiment in experiments:
        passed_dir = experiment / "passed"
        if passed_dir.is_dir():
            archives.extend(sorted(passed_dir.glob("*.tar.gz")))
    return archives


def group_members_by_protein(member_names: list[str]) -> dict[str, dict[str, str]]:
    """Same logic as stage_01_filter.py's group_members_by_protein --
    directory entries and any non-.json/.cif/.pdb member are ignored so
    they can't create a bogus empty-protein_id group."""
    groups: dict[str, dict[str, str]] = {}
    for member_name in member_names:
        suffix = Path(member_name).suffix.lower()
        if suffix not in (".json", ".cif", ".pdb"):
            continue
        pid = Path(member_name).stem
        entry = groups.setdefault(pid, {})
        if suffix == ".json":
            entry["json"] = member_name
        else:
            entry["structure"] = member_name
    return groups


def group_name_from_archive(archive_path: Path) -> str:
    """'htmedium.tar.gz' -> 'htmedium'. Path.stem only strips ONE suffix
    (would give 'htmedium.tar'), so this strips the full '.tar.gz' explicitly."""
    name = archive_path.name
    if name.endswith(".tar.gz"):
        return name[: -len(".tar.gz")]
    return archive_path.stem


# ---------------------------------------------------------------------------
# Merge: read existing destination members + new source members, write
# combined archive, atomically replace.
# ---------------------------------------------------------------------------

def read_archive_members(archive_path: Path) -> list[tuple[tarfile.TarInfo, bytes]]:
    members: list[tuple[tarfile.TarInfo, bytes]] = []
    with tarfile.open(archive_path, "r:gz") as archive:
        for member_info in archive.getmembers():
            if not member_info.isfile():
                continue
            extracted = archive.extractfile(member_info)
            if extracted is None:
                raise RuntimeError(f"could not read member {member_info.name!r} from {archive_path}")
            members.append((member_info, extracted.read()))
    return members


def merge_group_archive(
    source_path: Path,
    destination_path: Path,
    transferred_ids: set[str],
) -> list[str]:
    """
    Returns the list of protein_ids newly added to destination_path.
    transferred_ids is the GLOBAL set already recorded in the transfer
    log (checked against, not mutated here -- the caller updates it).
    """
    with tarfile.open(source_path, "r:gz") as source:
        groups = group_members_by_protein(source.getnames())
        candidate_ids = sorted(pid for pid in groups if pid not in transferred_ids)
        if not candidate_ids:
            return []

        valid_ids: list[str] = []
        for pid in candidate_ids:
            members = groups[pid]
            if "json" not in members or "structure" not in members:
                print(f"    [SKIP - incomplete pair] {pid}: has {sorted(members)} in {source_path.name}")
                continue
            valid_ids.append(pid)

        if not valid_ids:
            return []

        new_members: list[tuple[tarfile.TarInfo, bytes]] = []
        for pid in valid_ids:
            for kind in ("json", "structure"):
                member_name = groups[pid][kind]
                member_info = source.getmember(member_name)
                extracted = source.extractfile(member_info)
                if extracted is None:
                    raise RuntimeError(f"could not read member {member_name!r} from {source_path}")
                new_members.append((member_info, extracted.read()))

    existing_members: list[tuple[tarfile.TarInfo, bytes]] = []
    if destination_path.exists():
        existing_members = read_archive_members(destination_path)

    existing_names = {member_info.name for member_info, _ in existing_members}
    new_names = {member_info.name for member_info, _ in new_members}
    collision = existing_names & new_names
    if collision:
        raise RuntimeError(
            f"{destination_path}: transfer log inconsistency -- these members are already "
            f"in the destination archive but were about to be re-added: {sorted(collision)}"
        )

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination_path.parent / f".{destination_path.name}.tmp{os.getpid()}"
    with tarfile.open(temp_path, "w:gz") as out:
        for member_info, payload in existing_members + new_members:
            out.addfile(member_info, io.BytesIO(payload))

    with tarfile.open(temp_path, "r:gz") as check:
        names_in_temp = set(check.getnames())
    expected_names = existing_names | new_names
    missing = expected_names - names_in_temp
    if missing:
        temp_path.unlink()
        raise RuntimeError(f"Verification failed building {destination_path}: missing {sorted(missing)}")

    os.replace(temp_path, destination_path)
    return valid_ids


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_transfer(specific_experiment: str | None = None) -> int:
    if not SORTED_CLEAN_DIR.is_dir():
        print(f"RESULT: ERROR (sorted_clean not found: {SORTED_CLEAN_DIR.resolve()})")
        return 1

    log = load_transfer_log(TRANSFER_LOG_PATH)
    transferred_ids = set(log.keys())

    archives = discover_passed_archives(SORTED_CLEAN_DIR, specific_experiment)

    print("===== SYNC: sorted_clean/passed -> stage-02 inputs (incremental) =====")
    print(f"sorted_clean    = {SORTED_CLEAN_DIR.resolve()}")
    print(f"stage02_inputs  = {STAGE02_INPUTS_DIR.resolve()}")
    print(f"transfer_log    = {TRANSFER_LOG_PATH.resolve()}")
    print(f"already transferred (before this run) = {len(transferred_ids)}")
    print()

    newly_transferred_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    new_rows: list[dict[str, str]] = []

    for source_path in archives:
        experiment = source_path.parent.parent.name  # .../<experiment>/passed/<group>.tar.gz
        group = group_name_from_archive(source_path)
        destination_path = STAGE02_INPUTS_DIR / experiment / source_path.name

        new_ids = merge_group_archive(source_path, destination_path, transferred_ids)
        if not new_ids:
            print(f"  {experiment}/{source_path.name}: nothing new")
            continue

        print(f"  {experiment}/{source_path.name}: transferred {len(new_ids)} new protein(s)")
        for protein_id in new_ids:
            new_rows.append({
                "protein_id": protein_id,
                "experiment": experiment,
                "group": group,
                "source_archive": str(source_path.resolve()),
                "transferred_at": newly_transferred_at,
            })
            transferred_ids.add(protein_id)  # keep in sync within this same run

    if new_rows:
        for row in new_rows:
            log[row["protein_id"]] = row
        save_transfer_log(log, TRANSFER_LOG_PATH)

    print()
    print(f"Run complete. Newly transferred: {len(new_rows)}. Total transferred overall: {len(transferred_ids)}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiment", default=None,
                         help="Only scan this one experiment subfolder of sorted_clean/")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return run_transfer(specific_experiment=args.experiment)


if __name__ == "__main__":
    sys.exit(main())
