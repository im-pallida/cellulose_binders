#!/usr/bin/env python3
"""
Stage 04 (04_alphafold), helping script: diagnostic only, doesn't write
anything -- reads tables/stage_04_results.csv (written by
15082026_04_compute_rmsd.py) and prints the 10 sequences, among those
with dimer_matched_fraction > MIN_MATCHED_FRACTION, with the best
(lowest) dimer_rmsd / dimer_matched_pairs ratio.

Why this ratio: with 0/1500+ sequences passing the RMSD_THRESHOLD_DIMER /
MATCHED_FRACTION_THRESHOLD_DIMER filter outright, a straight "best RMSD"
or "most matched pairs" ranking alone can each be misleading on their own
(a tiny well-matched fragment can have a great RMSD; a huge alignment can
have a middling RMSD) -- dividing RMSD by matched_pairs rewards designs
that are BOTH accurate AND well-matched, and heavily penalizes a great
RMSD that only came from matching a handful of residues. Lower ratio is
better either way (low RMSD, high matched_pairs). Ranks by the DIMER
figures specifically since the dimer's joint superposition is the
stricter test (it captures relative chain arrangement, not just
per-chain fold accuracy).

MIN_MATCHED_FRACTION (0.30 = 30%) filters out rows before ranking, not
just before printing -- a sequence with a tiny but suspiciously low RMSD
over only a handful of matched residues is exactly the kind of row the
ratio alone can't fully protect against, so this cuts it before it can
occupy a top-10 slot.

A row with dimer_matched_pairs == 0 (or an unparseable/empty dimer_rmsd,
e.g. from a Kabsch fit that never got 3+ matched pairs at all) is skipped
outright -- there's no ratio to compute, and it's not informative anyway.

Usage:
    python3 15082026_04_rank_by_rmsd_ratio.py [--top N] [--min-matched-fraction F]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent.parent   # 04_alphafold/ (helping_scripts/ -> scripts/ -> 04_alphafold/)
RESULTS_TABLE_PATH = STAGE_ROOT / "tables" / "stage_04_results.csv"

MIN_MATCHED_FRACTION = 0.30   # 30% -- only rank sequences with dimer_matched_fraction above this


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        print(f"ERROR: {path} not found -- run 15082026_04_compute_rmsd.py first")
        raise SystemExit(1)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def rank_by_ratio(rows: list[dict[str, str]], min_matched_fraction: float) -> list[tuple[float, dict[str, str]]]:
    ranked = []
    skipped_zero_matched = 0
    skipped_unparseable = 0
    skipped_low_fraction = 0
    for row in rows:
        try:
            dimer_matched = int(row["dimer_matched_pairs"])
            dimer_rmsd = float(row["dimer_rmsd"])
            dimer_fraction = float(row["dimer_matched_fraction"])
        except (KeyError, ValueError):
            skipped_unparseable += 1
            continue
        if dimer_matched == 0:
            skipped_zero_matched += 1
            continue
        if dimer_fraction <= min_matched_fraction:
            skipped_low_fraction += 1
            continue
        ratio = dimer_rmsd / dimer_matched
        ranked.append((ratio, row))
    ranked.sort(key=lambda entry: entry[0])
    print(f"  ({skipped_zero_matched} row(s) skipped: dimer_matched_pairs == 0; "
          f"{skipped_unparseable} row(s) skipped: missing/unparseable dimer_rmsd/fraction; "
          f"{skipped_low_fraction} row(s) skipped: dimer_matched_fraction <= {min_matched_fraction:.0%})")
    return ranked


def print_ranked(ranked: list[tuple[float, dict[str, str]]], top_n: int) -> None:
    header = (
        f"{'sequence_id':<16}{'protein_id':<14}{'monomer_rmsd':>13}{'dimer_rmsd':>12}"
        f"{'mono_matched':>13}{'dimer_matched':>14}{'dimer_matched_%':>16}{'ratio':>10}"
    )
    print(header)
    print("-" * len(header))
    for ratio, row in ranked[:top_n]:
        dimer_fraction = float(row["dimer_matched_fraction"])
        print(
            f"{row['sequence_id']:<16}{row['protein_id']:<14}"
            f"{row.get('monomer_rmsd', ''):>13}{row.get('dimer_rmsd', ''):>12}"
            f"{row.get('monomer_matched_pairs', ''):>13}{row.get('dimer_matched_pairs', ''):>14}"
            f"{dimer_fraction:>15.1%} {ratio:>9.4f}"
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=10, help="how many to print (default 10)")
    ap.add_argument("--min-matched-fraction", type=float, default=MIN_MATCHED_FRACTION,
                     help=f"only rank sequences with dimer_matched_fraction above this (default {MIN_MATCHED_FRACTION})")
    args = ap.parse_args()

    rows = load_rows(RESULTS_TABLE_PATH)
    print(f"Loaded {len(rows)} row(s) from {RESULTS_TABLE_PATH}")
    ranked = rank_by_ratio(rows, args.min_matched_fraction)
    print(f"\nTop {min(args.top, len(ranked))} by dimer_rmsd / dimer_matched_pairs (lower = better), "
          f"restricted to dimer_matched_fraction > {args.min_matched_fraction:.0%}:\n")
    print_ranked(ranked, args.top)


if __name__ == "__main__":
    main()
