#!/usr/bin/env python3
"""
Build the 16 RFD3 "inputs" JSON configs for the SHORT-length variant of
the tyr-hairpin / strand-hairpin seed experiment
(02_rfd3_linear_symmetry_generation).

Identical to 15072026_01_new_seed_tyr_hairpin.py EXCEPT: pre/post size
ranges are scaled to 2/3 of the original (connector ranges unchanged -
those represent the interface loop, not overall protein length).

4 orders (th_str, ht_inv, sh_str, hs_inv) x 4 sizes (small, medium,
large, xlarge) = 16 JSON files.

FIX (24072026): seed paths are now computed relative to THIS SCRIPT'S
OWN FILE LOCATION (via _PROJECT_ROOT, same pattern as _DEFAULT_STAGE_ROOT
below) instead of Path("~/...").expanduser(). The old ~-based path
depends entirely on whoever runs this script and from which machine -
it silently baked in the WRONG path twice now (once as the workstation
user's home dir, once as a different user's cluster home dir), because
neither happened to match where the seed files actually live
(00_reference/ at the project root). Computing it from __file__ instead
means the JSON configs are correct regardless of who generates them or
on which machine, as long as the project's own directory layout is
intact - which is the one thing that's actually guaranteed.

Output layout:
    json/<date>_<experiment>/th_str/small.json
    json/<date>_<experiment>/th_str/medium.json
    json/<date>_<experiment>/th_str/large.json
    json/<date>_<experiment>/th_str/xlarge.json
    json/<date>_<experiment>/ht_inv/small.json
    ... etc for sh_str and hs_inv

Interactive usage:
    python 22072026_01_new_seed_short_seed.py --stage-root ..

    Enter experiment date (DDMMYYYY, e.g. 31012026): 31012026
    Enter experiment name: cellulose_v2

Non-interactive usage:
    python 22072026_01_new_seed_short_seed.py --stage-root .. \
        --date 31012026 --experiment-name cellulose_v2
"""

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------
# Fixed-atom selections (per seed)
# ---------------------------------------------------------------------

_COMMON_FIXED_ATOMS = {
    "A24": "BKBN",
    "A25": "BKBN",
    "A26": "BKBN",
    "A27": "BKBN",
    "A28": "BKBN",
    "A29": "ALL",
    "A30": "BKBN",
    "A31": "ALL",
    "A32": "ALL",
    "A33": "BKBN",
    "A34": "ALL",
    "A35": "BKBN",
    "A36": "BKBN",
}

_TH_FIXED_ATOMS = {"A5": "ALL"}
_SH_FIXED_ATOMS = {"A5": "ALL", "A6": "ALL", "A7": "ALL", "A8": "ALL", "A9": "ALL"}

# ---------------------------------------------------------------------
# Sizes: pre/post ranges (vary by size, not by order)
# Scaled to 2/3 of the original 15072026 values, e.g. small "12-22" ->
# round(12*2/3)=8, round(22*2/3)=15 -> "8-15". Connector ranges (in
# ORDERS below) are left unchanged - those represent the interface loop,
# not overall protein length.
# ---------------------------------------------------------------------

SIZES = {
    "small":  {"pre": "8-15",  "post": "8-15"},
    "medium": {"pre": "13-23", "post": "13-23"},
    "large":  {"pre": "19-30", "post": "19-30"},
    "xlarge": {"pre": "27-40", "post": "27-40"},
}

# ---------------------------------------------------------------------
# This script lives at 02_rfd3_linear_symmetry_generation/scripts/<this
# file>. Its grandparent directory is always the stage root, and the
# stage root's own parent is always the project root, regardless of
# the directory or user this script happens to be run from.
# ---------------------------------------------------------------------
_DEFAULT_STAGE_ROOT = Path(__file__).resolve().parent.parent
_PROJECT_ROOT = _DEFAULT_STAGE_ROOT.parent
_SEED_DIR = _PROJECT_ROOT / "00_reference_structures" / "seed"

# ---------------------------------------------------------------------
# Seeds: path + fixed-atom selection (vary by seed, not by size/order)
# ---------------------------------------------------------------------

SEEDS = {
    "seed_th": {
        "path": str(_SEED_DIR / "1xseed_nocellulose_TH.pdb"),
        "fixed_atoms": {**_COMMON_FIXED_ATOMS, **_TH_FIXED_ATOMS},
    },
    "seed_sh": {
        "path": str(_SEED_DIR / "1xseed_nocellulose_SH.pdb"),
        "fixed_atoms": {**_COMMON_FIXED_ATOMS, **_SH_FIXED_ATOMS},
    },
}

# ---------------------------------------------------------------------
# Orders: which seed, first/second chunk, connector range (vary by
# order, not by size)
# ---------------------------------------------------------------------

ORDERS = {
    "th_str": {"seed": "seed_th", "first_chunk": "A5",     "connector": "34-50", "second_chunk": "A24-36"},
    "ht_inv": {"seed": "seed_th", "first_chunk": "A24-36", "connector": "4-8",   "second_chunk": "A5"},
    "sh_str": {"seed": "seed_sh", "first_chunk": "A5-9",   "connector": "30-46", "second_chunk": "A24-36"},
    "hs_inv": {"seed": "seed_sh", "first_chunk": "A24-36", "connector": "4-8",   "second_chunk": "A5-9"},
}

SELECT_EXPOSED = "A5,A31,A32"
SYMMETRY_ID = "T2EXACT"
IS_SYMMETRIC_MOTIF = False


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def _bound_range(bound_str):
    lo, hi = bound_str.split("-")
    return int(lo), int(hi)


def _chunk_length(chunk_str):
    if "-" in chunk_str:
        lo_str, hi_str = chunk_str.split("-")
        lo = int(lo_str.split("A")[1])  # e.g. "A24" -> 24
        hi = int(hi_str)                # e.g. "36" -> 36 (no "A" prefix on this side)
        length = hi - lo + 1
        return length
    else:
        return 1


def build_length(size_params, order_params):
    pre = _bound_range(size_params["pre"])
    post = _bound_range(size_params["post"])
    conn = _bound_range(order_params["connector"])
    first_chunk = _chunk_length(order_params["first_chunk"])
    second_chunk = _chunk_length(order_params["second_chunk"])
    lo = pre[0] + post[0] + conn[0] + first_chunk + second_chunk
    hi = pre[1] + post[1] + conn[1] + first_chunk + second_chunk
    return f"{lo}-{hi}"


def build_contig(size_params, order_params) -> str:
    return (
        f"{size_params['pre']},"
        f"{order_params['first_chunk']},"
        f"{order_params['connector']},"
        f"{order_params['second_chunk']},"
        f"{size_params['post']}"
    )


def build_config(seed_params, order_params, size_params):
    cfg = {
        "input": seed_params["path"],
        "contig": build_contig(size_params, order_params),
        "length": build_length(size_params, order_params),
        "select_fixed_atoms": seed_params["fixed_atoms"],
        "select_exposed": SELECT_EXPOSED,
        "is_non_loopy": True,
        "plddt_enhanced": True,
        "symmetry": {
            "id": SYMMETRY_ID,
            "is_symmetric_motif": IS_SYMMETRIC_MOTIF,
        },
    }
    return cfg


# ---------------------------------------------------------------------
# Date / experiment-name prompting (same validation rules as the
# original 01_prepare_rfd3_sizes_batch.py script)
# ---------------------------------------------------------------------

def validate_date(date_str: str) -> str:
    if not re.fullmatch(r"\d{8}", date_str):
        raise ValueError(f"Date must be exactly 8 digits (DDMMYYYY), got: {date_str!r}")
    try:
        datetime.strptime(date_str, "%d%m%Y")
    except ValueError as e:
        raise ValueError(f"'{date_str}' is not a real calendar date (DDMMYYYY): {e}")
    return date_str


def validate_experiment_name(name: str) -> str:
    name = name.strip()
    if not name:
        raise ValueError("Experiment name cannot be empty")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise ValueError(
            f"Experiment name can only contain letters, digits, '_' and '-', got: {name!r}"
        )
    return name


# ---------------------------------------------------------------------
# Main: loop over orders x sizes, write 16 JSON files into
# json/<date>_<experiment>/<order_key>/<size_key>.json
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-root", default=str(_DEFAULT_STAGE_ROOT),
                     help="Path to 02_rfd3_linear_symmetry_generation/ "
                          "(default: inferred from this script's own location)")
    ap.add_argument("--date", default=None,
                     help="Experiment date, DDMMYYYY (e.g. 31012026). If omitted, you'll be prompted.")
    ap.add_argument("--experiment-name", default=None,
                     help="Experiment name (e.g. cellulose_v2). If omitted, you'll be prompted.")
    args = ap.parse_args()

    # Sanity check the seed files actually exist BEFORE writing any JSON -
    # fail loudly here rather than baking in a bad path that only shows
    # up as a cryptic FileNotFoundError deep inside an RFD3 job later.
    missing_seeds = [s["path"] for s in SEEDS.values() if not Path(s["path"]).is_file()]
    if missing_seeds:
        print("ERROR: seed file(s) not found on this machine:")
        for path in missing_seeds:
            print(f"  {path}")
        print("Check that 00_reference/seed/ exists at the project root relative to this script.")
        sys.exit(1)

    date_str = args.date
    if date_str is None:
        date_str = input("Enter experiment date (DDMMYYYY, e.g. 31012026): ").strip()
    try:
        date_str = validate_date(date_str)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    experiment_name = args.experiment_name
    if experiment_name is None:
        experiment_name = input("Enter experiment name: ").strip()
    try:
        experiment_name = validate_experiment_name(experiment_name)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    stage_root = Path(args.stage_root).resolve()
    experiment_dir = stage_root / "json" / f"{date_str}_{experiment_name}"

    for order_key, order_params in ORDERS.items():
        order_dir = experiment_dir / order_key
        order_dir.mkdir(parents=True, exist_ok=True)

        seed_params = SEEDS[order_params["seed"]]

        for size_key, size_params in SIZES.items():
            cfg = build_config(seed_params, order_params, size_params)
            # RFD3 expects the JSON to be {example_name: {args...}}, not a
            # flat args dict at the top level - wrap it under the size
            # name (matching the original single-seed script's convention).
            wrapped = {size_key: cfg}
            json_path = order_dir / f"{size_key}.json"
            with open(json_path, "w") as jf:
                json.dump(wrapped, jf, indent=2)
            print(f"Wrote {json_path}")

    print(f"\nExperiment folder: {experiment_dir}")


if __name__ == "__main__":
    main()
