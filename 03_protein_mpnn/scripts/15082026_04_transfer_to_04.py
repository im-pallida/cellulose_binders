#!/usr/bin/env python3
"""
Stage 03, script 4 (was 15082026_02_transfer_to_04.py -- renumbered now
that 15082026_02_summary.py and 15082026_03_best_seqs.py exist): transfer
the best-3-per-protein sequences from stage-03's outputs_clean/ into
stage-04's (04_alphafold) inputs/, incrementally -- only proteins not
already transferred are copied on each run.

Source of truth for WHICH structures:
    03_protein_mpnn/outputs_clean/<experiment>/<group>.tar.gz
(15082026_03_best_seqs.py writes outputs_clean nested by experiment, same
<experiment>/<group> nesting used everywhere else in this pipeline -- no
"passed/" middle level though, since there's no pass/reject concept for
best-seqs selection.) Destination mirrors the same nesting:
    04_alphafold/inputs/<experiment>/<group>.tar.gz
Same <group>.tar.gz filename as the source, just relocated under its own
<experiment>/ folder.

Members here are also shaped differently: each protein_id is a FOLDER
inside the tar (<protein_id>/<sequence_id>.fa, 3 files -- the 3 lowest-
score designs 15082026_03_best_seqs.py selected), not a flat
<protein_id>.ext file. So the transfer unit here is "every file under one
protein_id/ folder", not a json+structure pair -- there's no separate
kind-of-file to validate for completeness the way the earlier transfer
scripts checked for BOTH json and pdb/cif. If a protein folder has fewer
than 3 files (15082026_03_best_seqs.py already warns about this at
selection time when fewer than 3 sequences were scored), whatever is
there still transfers -- this script doesn't re-block on it, just notes
the count.

Why "incremental" instead of an exact-rebuild pattern: outputs_clean/
<experiment>/<group>.tar.gz IS already cumulative across stage-03 runs on
its own (15082026_03_best_seqs.py's write_merged_archive merges rather
than overwrites) -- but this script still tracks its own transferred set
rather than re-deriving "new" from a diff, for the same reason every
other transfer script in this pipeline does: it's the simplest way to
guarantee idempotency regardless of how the upstream archive evolves,
and it keeps this script's logic identical in shape to
15082026_02_transfer_to_03.py / 15082026_04_transfer_to_02.py. That
memory is job_runs/proteins_transferred_from_03_logs.tsv (stage-03 side,
keyed by protein_id).

Destination archives DO need to be cumulative (stage-04 needs to see
every transferred protein ever produced, not just the latest batch), so
this script merges into an existing destination archive rather than
overwriting it -- same decompress-append-recompress + atomic-replace
pattern as every other transfer script here.

Usage:
    python3 15082026_04_transfer_to_04.py
"""
from __future__ import annotations

import csv
import io
import os
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# CONSTANTS -- all relative to this script's own location, so it works
# unmodified if the whole project folder is copied to another machine/user.
# Assumes this script lives in 03_protein_mpnn/scripts/, same convention
# as every other stage's own transfer-forward script (the producing stage
# owns the script that pushes its results to the next stage).
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE03_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = STAGE03_ROOT.parent
STAGE04_ROOT = PROJECT_ROOT / "04_alphafold"

MPNN_OUTPUTS_CLEAN_DIR = STAGE03_ROOT / "outputs_clean"
STAGE04_INPUTS_DIR = STAGE04_ROOT / "inputs"
TRANSFER_LOG_PATH = STAGE03_ROOT / "job_runs" / "proteins_transferred_from_03_logs.tsv"
TRANSFER_LOG_DELIMITER = "\t"

TRANSFER_LOG_FIELDS = [
    "protein_id",
    "experiment",
    "group",
    "source_archive",
    "file_count",
    "transferred_at",
]


# ---------------------------------------------------------------------------
# Transfer log -- same load-whole-dict/write-whole-dict pattern as every
# other transfer script here, keyed on protein_id, TSV under job_runs/.
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
# Discovery + grouping
# ---------------------------------------------------------------------------

def discover_group_archives(outputs_clean_root: Path) -> list[Path]:
    """outputs_clean/<experiment>/<group>.tar.gz -- no 'passed/' middle level at this stage."""
    if not outputs_clean_root.is_dir():
        return []
    archives: list[Path] = []
    for experiment_dir in sorted(outputs_clean_root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        archives.extend(sorted(experiment_dir.glob("*.tar.gz")))
    return archives


def group_members_by_protein(member_names: list[str]) -> dict[str, list[str]]:
    """
    {'htl000001': ['htl000001/htl00000102.fa', 'htl000001/htl00000109.fa', 'htl000001/htl00000113.fa']}
    Each protein_id is a FOLDER here (unlike earlier stages' flat
    <protein_id>.ext members), so the group is every .fa member under that
    folder, whatever the count.
    """
    groups: dict[str, list[str]] = {}
    for member_name in member_names:
        if "/" not in member_name or not member_name.endswith(".fa"):
            continue
        pid = member_name.split("/", 1)[0]
        groups.setdefault(pid, []).append(member_name)
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
) -> dict[str, int]:
    """
    Returns {protein_id: file_count} for every protein_id newly added to
    destination_path. transferred_ids is the GLOBAL set already recorded
    in the transfer log (checked against, not mutated here -- the caller
    updates it).
    """
    with tarfile.open(source_path, "r:gz") as source:
        groups = group_members_by_protein(source.getnames())
        candidate_ids = sorted(pid for pid in groups if pid not in transferred_ids)
        if not candidate_ids:
            return {}

        new_members: list[tuple[tarfile.TarInfo, bytes]] = []
        file_counts: dict[str, int] = {}
        for pid in candidate_ids:
            member_names = groups[pid]
            if len(member_names) != 3:
                print(f"    [NOTE] {pid}: transferring {len(member_names)} file(s), not the usual 3")
            file_counts[pid] = len(member_names)
            for member_name in member_names:
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
    return file_counts


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_transfer() -> int:
    if not MPNN_OUTPUTS_CLEAN_DIR.is_dir():
        print(f"RESULT: ERROR (outputs_clean not found: {MPNN_OUTPUTS_CLEAN_DIR.resolve()})")
        return 1

    log = load_transfer_log(TRANSFER_LOG_PATH)
    transferred_ids = set(log.keys())

    archives = discover_group_archives(MPNN_OUTPUTS_CLEAN_DIR)

    print("===== SYNC: stage-03 outputs_clean -> stage-04 (04_alphafold) inputs (incremental) =====")
    print(f"outputs_clean   = {MPNN_OUTPUTS_CLEAN_DIR.resolve()}")
    print(f"stage04_inputs  = {STAGE04_INPUTS_DIR.resolve()}")
    print(f"transfer_log    = {TRANSFER_LOG_PATH.resolve()}")
    print(f"already transferred (before this run) = {len(transferred_ids)}")
    print()

    newly_transferred_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    new_rows: list[dict[str, str]] = []

    for source_path in archives:
        experiment = source_path.parent.name  # .../<experiment>/<group>.tar.gz
        group = group_name_from_archive(source_path)
        destination_path = STAGE04_INPUTS_DIR / experiment / source_path.name

        file_counts = merge_group_archive(source_path, destination_path, transferred_ids)
        if not file_counts:
            print(f"  {experiment}/{source_path.name}: nothing new")
            continue

        print(f"  {experiment}/{source_path.name}: transferred {len(file_counts)} new protein(s)")
        for protein_id, file_count in file_counts.items():
            new_rows.append({
                "protein_id": protein_id,
                "experiment": experiment,
                "group": group,
                "source_archive": str(source_path.resolve()),
                "file_count": str(file_count),
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


def main() -> int:
    return run_transfer()


if __name__ == "__main__":
    sys.exit(main())
