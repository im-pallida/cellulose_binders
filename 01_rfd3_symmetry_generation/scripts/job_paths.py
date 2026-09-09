"""
The funcs which help to keep the consistency for the paths search.
"""
from __future__ import annotations
 
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple
 
RUN_LIST_SUFFIX = "_run.tsv"
COUNTER_FILENAME = "global_seq_counter.txt"
 
# A run-list row: the json filename relative to json/<experiment>/, and the
# global sequence number assigned to that one structure.
RunRow = Tuple[str, int]
 
 
# Step 1. Naming.
 
 
def group_key_from_json_rel(json_rel: str) -> str:
    """'small.json' -> 'small'. json/<experiment>/ has no subfolders, so this
    is just the filename with the extension stripped."""
    return json_rel.removesuffix(".json")
 
 
def job_name(group_key: str, global_seq: int) -> str:
    """The stem shared by a structure's .cif, its .json, and its archive
    members: 'small' + 42 -> 'small_000042'."""
    return f"{group_key}_{global_seq:06d}"
 
 
#Step 2. Inputs.
 
 
def json_root(stage: Path) -> Path:
    return stage / "json"
 
 
def experiment_json_dir(stage: Path, experiment: str) -> Path:
    return json_root(stage) / experiment
 
 
def json_path(stage: Path, experiment: str, json_rel: str) -> Path:
    return experiment_json_dir(stage, experiment) / json_rel


def design_entries(data: object) -> "dict[str, dict]":
    """The design configs inside an RFD3 input json.
    Raises ValueError if the file is not shaped like either."""
    if not isinstance(data, dict):
        raise ValueError(f"top level must be an object, got {type(data).__name__}")
    if "input" in data or "contig" in data:
        return {"": data}
    if not data:
        raise ValueError("file contains no design entries")
    bad = [name for name, cfg in data.items() if not isinstance(cfg, dict)]
    if bad:
        raise ValueError(
            f"expected each top-level key to name a design config; "
            f"these are not objects: {bad}"
        )
    return dict(data)
 
 
def resolve_seed_path(stage: Path, raw_input: str, json_path_: Path) -> Optional[Path]:
    """Locate the seed structure a json's "input" field points at."""
    candidate = Path(raw_input)
    if candidate.is_absolute():
        return candidate if candidate.is_file() else None
    for base in (stage, json_path_.parent):
        resolved = (base / candidate).resolve()
        if resolved.is_file():
            return resolved
    return None
 
 
# Step 3. Outputs.
 
 
def raw_dir(stage: Path, experiment: str, group_key: str) -> Path:
    return stage / "outputs_raw" / experiment / group_key
 
 
def raw_cif_path(stage: Path, experiment: str, group_key: str, name: str) -> Path:
    return raw_dir(stage, experiment, group_key) / f"{name}.cif"
 
 
def raw_json_path(stage: Path, experiment: str, group_key: str, name: str) -> Path:
    return raw_dir(stage, experiment, group_key) / f"{name}.json"
 
 
def tmp_work_dir(stage: Path, experiment: str, group_key: str, name: str) -> Path:
    return stage / "outputs_raw" / ".tmp" / experiment / group_key / name
 
 
def clean_dir(stage: Path, experiment: str) -> Path:
    return stage / "outputs_clean" / experiment
 
 
def clean_archive_path(stage: Path, experiment: str, group_key: str) -> Path:
    return clean_dir(stage, experiment) / f"{group_key}.tar.gz"
 
 
def expected_raw_cif(stage: Path, experiment: str, json_rel: str, global_seq: int) -> Path:
    """Where a given run-list row's structure will land. Used to check whether
    a job produced its output without re-deriving the naming rules."""
    group_key = group_key_from_json_rel(json_rel)
    return raw_cif_path(stage, experiment, group_key, job_name(group_key, global_seq))
 
 
# Step 4. Run list.
 
 
def run_list_path(stage: Path, experiment: str) -> Path:
    return stage / "job_runs" / f"{experiment}{RUN_LIST_SUFFIX}"
 
 
def counter_path(stage: Path) -> Path:
    return stage / "job_runs" / COUNTER_FILENAME
 
 
def write_run_list(path: Path, rows: Sequence[RunRow]) -> None:
    """One '<json_filename>\\t<global_seq>' row per structure to generate."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for json_rel, global_seq in rows:
            handle.write(f"{json_rel}\t{global_seq}\n")
 
 
def read_run_list(path: Path) -> List[RunRow]:
    """Parse a run list, reporting the offending line number on bad input.
 
    Blank lines and '#' comments are skipped so a run list stays hand-editable.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"could not read run list {path}: {exc}") from exc
 
    rows: List[RunRow] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) != 2:
            raise ValueError(
                f"{path}:{lineno}: expected '<json_filename><TAB><global_seq>', "
                f"got {line!r}"
            )
        json_rel, raw_seq = fields[0].strip(), fields[1].strip()
        if not raw_seq.isdecimal():
            raise ValueError(f"{path}:{lineno}: global_seq must be an integer, got {raw_seq!r}")
        rows.append((json_rel, int(raw_seq)))
    return rows
 
 
def group_rows(rows: Iterable[RunRow]) -> "dict[str, List[RunRow]]":
    """Run-list rows bucketed by group_key, preserving order within a group.
    Archiving works one group at a time because each group has its own tarball."""
    grouped: "dict[str, List[RunRow]]" = {}
    for json_rel, global_seq in rows:
        grouped.setdefault(group_key_from_json_rel(json_rel), []).append((json_rel, global_seq))
    return grouped
