"""Is the flag-absent default path byte-identical across the change?

Compares two run folders cell by cell over every BEHAVIOURAL column of
`round_metrics.csv`, keyed on (method, seed, round). Wall-clock columns are excluded by
name and the exclusion is printed, so a column cannot be quietly dropped to make the
check pass.

The standard is branch worktree-agent-a14c69a80e0afa44a's: 15,840 cells across all 22
behavioural columns agreeing exactly over 3 arms x 48 seeds x 5 rounds. Here it is 4 arms.

Usage:
  identity_check_abf.py <dirA> <dirB>
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys

#: Excluded because they measure the machine, not the algorithm. Named explicitly.
WALL_CLOCK = ("round_seconds", "cumulative_seconds")
KEY = ("method", "seed", "round")


def load(root: str):
    pats = sorted(glob.glob(os.path.join(root, "run_*", "round_metrics.csv")))
    if not pats:
        raise SystemExit(f"no round_metrics.csv under {os.path.join(root, 'run_*')}")
    path = pats[-1]
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    man = os.path.join(os.path.dirname(path), "manifest.json")
    meta = {}
    if os.path.isfile(man):
        m = json.load(open(man, encoding="utf-8"))
        meta = {
            "run_id": m.get("run_id"),
            "fingerprint": m.get("configuration_fingerprint"),
            "commit": (m.get("code") or {}).get("git_commit")
            or (m.get("environment") or {}).get("git_commit"),
            "seed_block": ((m.get("configuration") or {}).get("seeding") or {})
            .get("seed_block_sha256"),
            "wt_anchor": ((m.get("configuration") or {}).get("components") or {})
            .get("wt_anchor"),
        }
    return path, rows, meta


def main() -> int:
    a_root, b_root = sys.argv[1], sys.argv[2]
    a_path, a_rows, a_meta = load(a_root)
    b_path, b_rows, b_meta = load(b_root)

    print("A:", a_path)
    print("   ", json.dumps(a_meta, default=str)[:400])
    print("B:", b_path)
    print("   ", json.dumps(b_meta, default=str)[:400])
    print()

    cols_a = list(a_rows[0].keys())
    cols_b = list(b_rows[0].keys())
    if cols_a != cols_b:
        print("COLUMN SETS DIFFER")
        print("  only in A:", [c for c in cols_a if c not in cols_b])
        print("  only in B:", [c for c in cols_b if c not in cols_a])
    shared = [c for c in cols_a if c in cols_b]
    behavioural = [c for c in shared if c not in WALL_CLOCK and c not in KEY]
    print(f"columns in A: {len(cols_a)}   key: {len(KEY)}   "
          f"excluded as wall clock: {list(WALL_CLOCK)}")
    print(f"BEHAVIOURAL COLUMNS COMPARED ({len(behavioural)}): {behavioural}")
    print()

    def index(rows):
        out = {}
        for r in rows:
            out[tuple(r[k] for k in KEY)] = r
        return out

    ia, ib = index(a_rows), index(b_rows)
    only_a = sorted(set(ia) - set(ib))
    only_b = sorted(set(ib) - set(ia))
    both = sorted(set(ia) & set(ib))
    print(f"rows A={len(a_rows)} B={len(b_rows)} shared keys={len(both)} "
          f"only-A={len(only_a)} only-B={len(only_b)}")
    if only_a[:5] or only_b[:5]:
        print("  only-A sample:", only_a[:5])
        print("  only-B sample:", only_b[:5])

    diffs = []
    cells = 0
    for k in both:
        for c in behavioural:
            cells += 1
            if ia[k][c] != ib[k][c]:
                diffs.append((k, c, ia[k][c], ib[k][c]))

    print()
    print(f"CELLS COMPARED: {cells:,}   DIFFERING: {len(diffs):,}")
    if diffs:
        print("VERDICT: THE DEFAULT PATH MOVED. First 40 differences:")
        for (k, c, va, vb) in diffs[:40]:
            print(f"  {k} {c}: A={va!r} B={vb!r}")
        return 1
    if only_a or only_b:
        print("VERDICT: every shared cell agrees, but the ROW SETS differ - not identity.")
        return 1
    print("VERDICT: IDENTICAL on every behavioural cell. The default path is untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
