"""
merge_default_identity.py - is the flag-absent DEFAULT PATH byte-identical across two runs?

The validity condition of docs/prereg-combined-head-to-head.md section 5.1. Two branches were
merged onto one tree; both modify benchmark_panoptes_isolated.py and one modifies
workers/math_worker.py. If the merge perturbed the default path, every figure already
published in docs/HANDOFF.md was measured on a different arm than the merged tree runs, the
campaign's `panoptes_default` column is not comparable to anything, and the campaign is VOID.

So the default is PROVEN untouched rather than asserted, exactly as
`worktree-agent-a14c69a80e0afa44a` proved its own change: run the same arms, the same seed
block, the same budget, at the pre-change commit and at the merged commit, and compare EVERY
behavioural cell of round_metrics.csv.

The two wall-clock columns are exempt and named, not silently dropped - `round_seconds` and
`cumulative_seconds` cannot be identical across two executions and a check that demanded it
would always fail. Every other column is compared, including the keys, so the row sets
themselves must agree.

`configuration.seeding.seed_block_sha256` is compared too: if two runs disagree there they
did not see the same starting samples, and identical behavioural columns would be a
coincidence rather than evidence.

Usage:
  python merge_default_identity.py <run_dir_before> <run_dir_after> [--arms a,b,c]
Exit code 0 only if the two runs are identical on every compared cell.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple

#: Cannot be identical across two executions; exempt BY NAME so the exemption is visible.
WALL_CLOCK = ("round_seconds", "cumulative_seconds")


def load(run_dir: Path) -> Tuple[Dict[Tuple[str, int, int], dict], List[str]]:
    path = run_dir / "round_metrics.csv"
    if not path.is_file():
        raise SystemExit(f"no round_metrics.csv at {path}")
    rows: Dict[Tuple[str, int, int], dict] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or [])
        for line in reader:
            key = (line["method"], int(float(line["seed"])), int(float(line["round"])))
            if key in rows:
                raise SystemExit(f"duplicate row {key} in {path}")
            rows[key] = line
    return rows, columns


def seed_digest(run_dir: Path) -> object:
    path = run_dir / "manifest.json"
    if not path.is_file():
        return "NO MANIFEST (a run folder without one is a partial launch, not data)"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    return (((manifest.get("configuration") or {}).get("seeding") or {})
            .get("seed_block_sha256"))


def commit_of(run_dir: Path) -> object:
    path = run_dir / "manifest.json"
    if not path.is_file():
        return None
    manifest = json.loads(path.read_text(encoding="utf-8"))
    code = ((manifest.get("configuration") or {}).get("code") or {})
    return {"commit": code.get("commit"), "dirty_tree": code.get("dirty_tree"),
            "source_tree_sha256": code.get("source_tree_sha256")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("before")
    parser.add_argument("after")
    parser.add_argument("--arms", default="",
                        help="comma-separated arms to compare; default is every arm present "
                             "in BOTH runs")
    args = parser.parse_args()

    before_dir, after_dir = Path(args.before), Path(args.after)
    before, cols_b = load(before_dir)
    after, cols_a = load(after_dir)

    print("=" * 100)
    print("DEFAULT-PATH IDENTITY CHECK")
    print(f"  BEFORE {before_dir}")
    print(f"         code {commit_of(before_dir)}")
    print(f"  AFTER  {after_dir}")
    print(f"         code {commit_of(after_dir)}")
    print("=" * 100)

    if cols_b != cols_a:
        print(f"!! COLUMN SETS DIFFER\n   before: {cols_b}\n   after:  {cols_a}")
        return 2
    compared = [c for c in cols_b if c not in WALL_CLOCK]
    print(f"\ncolumns in round_metrics.csv: {len(cols_b)}")
    print(f"compared: {len(compared)}   exempt (wall clock): {list(WALL_CLOCK)}")

    arms_b = {k[0] for k in before}
    arms_a = {k[0] for k in after}
    if args.arms:
        wanted = {a.strip() for a in args.arms.split(",") if a.strip()}
    else:
        wanted = arms_b & arms_a
    missing = wanted - (arms_b & arms_a)
    if missing:
        print(f"!! arms requested but not present in both runs: {sorted(missing)}")
        return 2
    arms = sorted(wanted)
    print(f"arms compared: {arms}")
    print(f"arms only in BEFORE: {sorted(arms_b - arms_a) or 'none'}")
    print(f"arms only in AFTER:  {sorted(arms_a - arms_b) or 'none'}")

    digest_b, digest_a = seed_digest(before_dir), seed_digest(after_dir)
    same_seeds = digest_b == digest_a and digest_b is not None
    seed_call = ("IDENTICAL" if same_seeds else
                 "DIFFER - the runs did NOT see the same starting samples")
    print(f"\nseed_block_sha256  before {digest_b}  after {digest_a}  -> {seed_call}")

    keys_b = {k for k in before if k[0] in wanted}
    keys_a = {k for k in after if k[0] in wanted}
    only_b, only_a = sorted(keys_b - keys_a), sorted(keys_a - keys_b)
    if only_b or only_a:
        print(f"\n!! ROW SETS DIFFER. only in before: {len(only_b)} e.g. {only_b[:5]}; "
              f"only in after: {len(only_a)} e.g. {only_a[:5]}")

    shared = sorted(keys_b & keys_a)
    seeds = sorted({k[1] for k in shared})
    rounds = sorted({k[2] for k in shared})
    mismatches: List[str] = []
    per_column: Dict[str, int] = {}
    for key in shared:
        for column in compared:
            if before[key][column] != after[key][column]:
                per_column[column] = per_column.get(column, 0) + 1
                if len(mismatches) < 40:
                    mismatches.append(
                        f"{key}  {column}: {before[key][column]!r} -> "
                        f"{after[key][column]!r}")
    cells = len(shared) * len(compared)

    print(f"\n{len(arms)} arm(s) x {len(seeds)} seed(s) x {len(rounds)} round(s) = "
          f"{len(shared)} row(s); {len(compared)} compared column(s) = {cells:,} cells")
    if mismatches or only_b or only_a or not same_seeds:
        print(f"\n!! NOT IDENTICAL: {sum(per_column.values()):,} differing cell(s) in "
              f"{len(per_column)} column(s)")
        for column, count in sorted(per_column.items(), key=lambda kv: -kv[1]):
            print(f"   {column}: {count}")
        print("\nfirst differences:")
        for line in mismatches:
            print(f"   {line}")
        print("\nTHE DEFAULT PATH MOVED. Under the pre-registration's validity condition 5.1 "
              "this is a FINDING and the campaign is VOID.")
        return 1

    print(f"\nIDENTICAL: all {cells:,} cells across {len(compared)} behavioural columns of "
          f"round_metrics.csv agree exactly,")
    print(f"for every one of {len(arms)} arm(s), every one of {len(seeds)} seed(s) and every "
          f"one of {len(rounds)} round(s).")
    print(f"The only exempt columns are {list(WALL_CLOCK)}, which are wall clock.")
    print("The flag-absent default path is UNTOUCHED by the merge.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
