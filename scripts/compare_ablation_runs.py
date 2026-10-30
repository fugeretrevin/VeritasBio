"""
scripts/compare_ablation_runs.py - paired comparison of two benchmark run folders.

The benchmark reports each arm's mean independently, which is not enough to decide
whether two configurations differ: the seeds are paired by construction (every arm
sees the identical starting plate per seed), so the paired difference has far less
variance than the two means suggest. Comparing means alone both overstates
differences that are seed noise and hides real ones.

Reads round_metrics.csv from two run folders, pairs them on (method, seed), and
reports the paired t statistic on final best fitness and on top-1% recall, plus how
many seeds moved in each direction.

Usage:
    python scripts/compare_ablation_runs.py <baseline_run_dir> <variant_run_dir>
    python scripts/compare_ablation_runs.py --list
"""
from __future__ import annotations

import argparse
import csv
import math
import pathlib
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmark_comparison_results"


def final_per_seed(run_dir: pathlib.Path) -> dict:
    """Final-round metrics per (method, seed), keyed by the last round present."""
    path = run_dir / "round_metrics.csv"
    if not path.is_file():
        raise FileNotFoundError(f"no round_metrics.csv in {run_dir}")
    latest: dict = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("status") and row["status"].lower() not in ("ok", "completed", ""):
                continue
            try:
                key = (row["method"], int(row["seed"]))
                rnd = int(row["round"])
            except (KeyError, ValueError):
                continue
            prev = latest.get(key)
            if prev is None or rnd > prev[0]:
                try:
                    latest[key] = (
                        rnd,
                        float(row["best_measured_fitness"]),
                        float(row["cum_top1_hits_including_initial"]),
                    )
                except (KeyError, ValueError):
                    continue
    return {k: (v[1], v[2]) for k, v in latest.items()}


def paired_t(diffs: list) -> tuple:
    """Return (mean difference, t statistic). t is nan when it is undefined."""
    n = len(diffs)
    if n < 2:
        return (statistics.fmean(diffs) if diffs else 0.0), float("nan")
    mean = statistics.fmean(diffs)
    sd = statistics.stdev(diffs)
    if sd == 0.0:
        return mean, float("inf") if mean != 0.0 else 0.0
    return mean, mean / (sd / math.sqrt(n))


def describe(run_dir: pathlib.Path) -> str:
    """Pull the arm's ablation note out of report.txt, if it recorded one."""
    report = run_dir / "report.txt"
    if not report.is_file():
        return ""
    for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
        if "ABLATION=" in line:
            return line.strip()
    return "ABLATION=none (no ablation note recorded)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Paired comparison of two benchmark runs.")
    ap.add_argument("baseline", nargs="?", type=pathlib.Path)
    ap.add_argument("variant", nargs="?", type=pathlib.Path)
    ap.add_argument("--list", action="store_true", help="list available run folders")
    args = ap.parse_args(argv)

    if args.list or not (args.baseline and args.variant):
        runs = sorted(RESULTS.glob("run_*"), key=lambda p: p.name, reverse=True)
        if not runs:
            print(f"no run folders under {RESULTS}", file=sys.stderr)
            return 2
        print("available runs (most recent first):")
        for r in runs[:20]:
            note = describe(r)
            print(f"  {r.name}   {note}")
        return 0 if args.list else 2

    base_rows = final_per_seed(args.baseline)
    var_rows = final_per_seed(args.variant)
    shared = sorted(set(base_rows) & set(var_rows))
    if not shared:
        print("no (method, seed) pairs in common - were both runs given the same seeds?",
              file=sys.stderr)
        return 1

    print(f"baseline : {args.baseline.name}")
    print(f"           {describe(args.baseline)}")
    print(f"variant  : {args.variant.name}")
    print(f"           {describe(args.variant)}")
    print(f"paired on {len(shared)} (method, seed) pairs")
    only_base = len(set(base_rows) - set(var_rows))
    only_var = len(set(var_rows) - set(base_rows))
    if only_base or only_var:
        print(f"  note: {only_base} pair(s) only in baseline, {only_var} only in variant - excluded")
    print()

    crit = 2.01  # two-sided 0.05 at n=48; close enough for n near that
    print(f"{'metric':<16}{'baseline':>10}{'variant':>10}{'mean diff':>12}{'t':>9}"
          f"{'verdict':>18}{'won/lost/tied':>16}")
    print("-" * 91)
    for idx, label in ((0, "best fitness"), (1, "top-1% hits")):
        b = [base_rows[k][idx] for k in shared]
        v = [var_rows[k][idx] for k in shared]
        diffs = [vi - bi for bi, vi in zip(b, v)]
        mean, t = paired_t(diffs)
        won = sum(1 for d in diffs if d > 0)
        lost = sum(1 for d in diffs if d < 0)
        tied = sum(1 for d in diffs if d == 0)
        if math.isnan(t):
            verdict = "undetermined"
        elif abs(t) >= crit:
            verdict = "SIGNIFICANT" if mean > 0 else "SIGNIFICANT (worse)"
        else:
            verdict = "not separable"
        print(f"{label:<16}{statistics.fmean(b):>10.4f}{statistics.fmean(v):>10.4f}"
              f"{mean:>+12.4f}{t:>9.2f}{verdict:>18}{f'{won}/{lost}/{tied}':>16}")

    print()
    print(f"Paired on identical seeds; |t| >= {crit} is significant at the 0.05 level.")
    print("'variant' minus 'baseline', so a negative mean difference means the variant is worse.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
