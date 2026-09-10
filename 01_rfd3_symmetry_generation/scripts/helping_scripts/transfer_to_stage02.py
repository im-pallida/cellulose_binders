#!/usr/bin/env python3
"""
Hand stage-01's passed structures to stage 02.

    sorted_clean/<experiment>/passed/<group>.tar.gz              <- source
    ../02_geometry_filtering/inputs/<experiment>/<group>.tar.gz  <- destination

Same <group>.tar.gz filename, just relocated. rejected/ is never read: only
passed structures become stage-02 inputs."""
from __future__ import annotations

import argparse
import io
import os
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import job_paths as jp  # noqa: E402

# Shared with the filter on purpose. If the two disagreed about what counts as
# a protein's file pair, structures would pass the filter and then silently
# fail to transfer.
from filter_designs import group_members_by_protein, read_archive_members  # noqa: E402

STAGE = Path(__file__).resolve().parents[2]

PASSED = "passed"


class TransferError(Exception):
    """A transfer could not be completed -- not a structure being skipped."""


def _log(*parts: object) -> None:
    print(*parts, flush=True)


@dataclass
class TransferResult:
    """One destination archive this run added structures to."""
    experiment: str
    group_key: str
    source: Path
    destination: Path
    added: int = 0
    total: int = 0


@dataclass
class TransferReport:
    transferred: int = 0
    already_present: int = 0
    incomplete: List[str] = field(default_factory=list)
    results: List[TransferResult] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """An incomplete pair means a passed structure cannot reach stage 02."""
        return not self.incomplete

    def summary(self) -> str:
        if not self.transferred:
            return f"transfer: nothing new ({self.already_present} already in stage 02)"
        return (f"transfer: {self.transferred} structure(s) -> stage 02, "
                f"{self.already_present} already there, "
                f"{len(self.results)} archive(s) updated")


def passed_archives(stage: Path, experiment: Optional[str] = None) -> Dict[str, List[Path]]:
    """{experiment: [passed/<group>.tar.gz, ...]}.

    Defaults to every filtered experiment, matching the filter's own
    sweep-everything behaviour, so a batch an interrupted run left behind is
    picked up by the next launch whichever experiment that launch was for.
    """
    names = [experiment] if experiment else jp.sorted_experiments(stage)
    found: Dict[str, List[Path]] = {}
    for name in names:
        passed_dir = jp.sorted_clean_dir(stage, name, PASSED)
        if not passed_dir.is_dir():
            continue
        archives = sorted(passed_dir.glob("*.tar.gz"))
        if archives:
            found[name] = archives
    return found


def transferred_ids(destination: Path) -> Set[str]:
    """The protein_ids already in a destination archive.

    This replaces the transfer log 15082026_04 kept: both sides are cumulative
    now, so the destination answers the question directly and cannot drift
    from it. Reuses the filter's grouping so a member naming rule is defined
    in exactly one place.
    """
    if not destination.is_file():
        return set()
    with tarfile.open(destination, "r:gz") as archive:
        return set(group_members_by_protein(archive.getnames()))


def _read_new_members(
    source_path: Path, already: Set[str]
) -> Tuple[List[Tuple[tarfile.TarInfo, bytes]], List[str], int]:
    """Members of `source_path` for proteins not yet in the destination.

    Returns (members, incomplete_ids, skipped_count). A protein missing either
    half of its pair goes to incomplete_ids and contributes no members.
    """
    members: List[Tuple[tarfile.TarInfo, bytes]] = []
    incomplete: List[str] = []

    with tarfile.open(source_path, "r:gz") as source:
        groups = group_members_by_protein(source.getnames())
        candidates = sorted(pid for pid in groups if pid not in already)
        skipped = len(groups) - len(candidates)

        for protein_id in candidates:
            entry = groups[protein_id]
            if "json" not in entry or "structure" not in entry:
                incomplete.append(protein_id)
                _log(f"  [incomplete] {protein_id}: has only {sorted(entry)} "
                     f"in {source_path.name}; not transferred")
                continue
            for kind in ("json", "structure"):
                member_info = source.getmember(entry[kind])
                extracted = source.extractfile(member_info)
                if extracted is None:
                    raise TransferError(
                        f"could not read member {entry[kind]!r} from {source_path}"
                    )
                members.append((member_info, extracted.read()))

    return members, incomplete, skipped


def _merge_members(
    destination: Path, new_members: Sequence[Tuple[tarfile.TarInfo, bytes]]
) -> Tuple[int, int]:
    """Write existing + new members to a temp archive, verify, swap it in.

    Returns (added, total) counted in members, not structures -- each
    structure contributes a json and a structure file.
    """
    existing = read_archive_members(destination) if destination.exists() else []
    existing_names = {info.name for info, _ in existing}
    new_names = {info.name for info, _ in new_members}

    # Cannot happen: candidates were filtered against these very members.
    # Kept as an invariant -- if it ever fires, the grouping rule changed.
    collision = existing_names & new_names
    if collision:
        raise TransferError(
            f"{destination}: already contains {sorted(collision)} but they were "
            f"about to be added again"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_name(f".{destination.name}.tmp{os.getpid()}")
    try:
        with tarfile.open(temp_path, "w:gz") as out:
            for member_info, payload in list(existing) + list(new_members):
                out.addfile(member_info, io.BytesIO(payload))

        with tarfile.open(temp_path, "r:gz") as check:
            present = {Path(name).name for name in check.getnames()}
        missing = (existing_names | new_names) - present
        if missing:
            raise TransferError(
                f"{destination}: verification failed, missing {sorted(missing)}"
            )
        os.replace(temp_path, destination)
    finally:
        temp_path.unlink(missing_ok=True)

    return len(new_names), len(existing_names | new_names)


def run_transfer(stage: Path, experiment: Optional[str] = None) -> TransferReport:
    """Move every passed structure not already in stage 02 across."""
    report = TransferReport()
    by_experiment = passed_archives(stage, experiment)
    if not by_experiment:
        # Distinguish "the filter has not run" from "it ran and nothing passed".
        # regroup_and_archive only creates an outcome directory when something
        # lands in it, so a zero-yield experiment has rejected/ but no passed/,
        # and reporting that as "nothing under sorted_clean" sends people
        # looking for a bug that is not there.
        filtered = jp.sorted_experiments(stage)
        if filtered:
            _log(f"[transfer] no passed/ archive in any of {len(filtered)} filtered "
                 f"experiment(s) ({', '.join(filtered)}) -- nothing has passed the "
                 f"filter yet, so stage 02 has nothing to receive")
        else:
            _log(f"[transfer] nothing filtered yet under {jp.sorted_clean_root(stage)}")
        return report

    _log(f"[transfer] stage 02 inputs -> {jp.stage02_root(stage)}")

    for experiment_name, sources in by_experiment.items():
        for source_path in sources:
            group_key = source_path.name[: -len(".tar.gz")]
            destination = jp.stage02_archive_path(stage, experiment_name, group_key)

            already = transferred_ids(destination)
            members, incomplete, skipped = _read_new_members(source_path, already)
            report.already_present += skipped
            report.incomplete.extend(incomplete)

            if not members:
                continue

            added, total = _merge_members(destination, members)
            structures = len(members) // 2
            report.transferred += structures
            report.results.append(
                TransferResult(experiment_name, group_key, source_path,
                               destination, added, total)
            )
            _log(f"[transfer] {experiment_name}/{group_key}: +{structures} structure(s) "
                 f"({added} new member(s), {total} total) -> {destination}")

    return report


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--stage", type=Path, default=STAGE,
                        help=f"stage root directory (default: {STAGE})")
    parser.add_argument("--experiment", default=None,
                        help="only transfer this one experiment (default: all of them)")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    try:
        report = run_transfer(args.stage.resolve(), args.experiment)
    except (OSError, TransferError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    _log(f"[done] {report.summary()}")
    return 0 if report.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
