"""
Hand a stage's passed structures to the next stage.
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import job_paths as jp  # noqa: E402
from archives import (  # noqa: E402
    ArchiveError,
    group_members_by_protein,
    merge_members_into_archive,
)

PASSED = "passed"

# (stage root, experiment, group_key) -> the archive this group belongs in
DestinationFor = Callable[[Path, str, str], Path]


class TransferError(Exception):
    """A transfer could not be completed -- not a structure being skipped."""


def _log(*parts: object) -> None:
    print(*parts, flush=True)


@dataclass
class TransferResult:
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
        """An incomplete pair means a passed structure cannot move on."""
        return not self.incomplete

    def summary(self) -> str:
        if not self.transferred:
            return f"transfer: nothing new ({self.already_present} already there)"
        return (f"transfer: {self.transferred} structure(s) moved on, "
                f"{self.already_present} already there, "
                f"{len(self.results)} archive(s) updated")


def passed_archives(stage: Path, experiment: Optional[str] = None) -> Dict[str, List[Path]]:
    """{experiment: [passed/<group>.tar.gz, ...]}."""
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
    """The protein_ids already in a destination archive -- the record of what
    has moved. Reuses the shared grouping so the rule for what makes a
    protein's file pair is defined in exactly one place."""
    if not destination.is_file():
        return set()
    with tarfile.open(destination, "r:gz") as archive:
        return set(group_members_by_protein(archive.getnames()))


def _read_new_members(
    source_path: Path, already: Set[str]
) -> Tuple[List[Tuple[tarfile.TarInfo, bytes]], List[str], int]:
    """Members for proteins not yet in the destination, plus the ids whose
    file pair is incomplete and the count of those already transferred."""
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


def run_transfer(
    stage: Path, destination_for: DestinationFor, label: str,
    experiment: Optional[str] = None,
) -> TransferReport:
    """Move every passed structure not already at the destination."""
    report = TransferReport()
    by_experiment = passed_archives(stage, experiment)
    if not by_experiment:
        filtered = jp.sorted_experiments(stage)
        if filtered:
            _log(f"[transfer] no passed/ archive in any of {len(filtered)} filtered "
                 f"experiment(s) ({', '.join(filtered)}) -- nothing has passed the "
                 f"filter yet, so {label} has nothing to receive")
        else:
            _log(f"[transfer] nothing filtered yet under {jp.sorted_clean_root(stage)}")
        return report

    for experiment_name, sources in by_experiment.items():
        for source_path in sources:
            group_key = source_path.name[: -len(".tar.gz")]
            destination = destination_for(stage, experiment_name, group_key)
            if not report.results:
                _log(f"[transfer] {label} inputs -> {destination.parent.parent}")

            already = transferred_ids(destination)
            members, incomplete, skipped = _read_new_members(source_path, already)
            report.already_present += skipped
            report.incomplete.extend(incomplete)
            if not members:
                continue

            added, total = merge_members_into_archive(destination, members)
            structures = len(members) // 2
            report.transferred += structures
            report.results.append(
                TransferResult(experiment_name, group_key, source_path,
                               destination, added, total)
            )
            _log(f"[transfer] {experiment_name}/{group_key}: +{structures} structure(s) "
                 f"({added} new member(s), {total} total) -> {destination}")

    return report


def cli(destination_for: DestinationFor, label: str, default_stage: Path,
        doc: str, argv: Optional[List[str]] = None) -> int:
    """The whole command-line side, so a stage's own file is just its
    destination plus a shebang."""
    parser = argparse.ArgumentParser(
        description=doc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--stage", type=Path, default=default_stage,
                        help=f"stage root directory (default: {default_stage})")
    parser.add_argument("--experiment", default=None,
                        help="only transfer this one experiment (default: all of them)")
    args = parser.parse_args(argv)
    try:
        report = run_transfer(args.stage.resolve(), destination_for, label, args.experiment)
    except (OSError, ArchiveError, TransferError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    _log(f"[done] {report.summary()}")
    return 0 if report.ok else 1

