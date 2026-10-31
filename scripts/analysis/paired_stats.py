"""
paired_stats.py - paired per-seed statistics over one benchmark run folder.

The unit of independence is the SEED: every method in a run sees the same shared starting
sample for a given seed, so arms are paired seed by seed and the comparison is within-seed.
Nothing here pools across seeds before differencing.

Reports, for each method against a reference method:
  - mean and sd of the endpoint, per arm
  - the paired difference (mean, sd, 95% bootstrap CI over seeds)
  - Wilcoxon signed-rank (two-sided, exact where scipy chooses it) on the paired differences
  - a sign test (binomial) on wins/losses/ties, which makes no distributional assumption
  - the per-arm mean of n_unique_measurements, which is a FAITHFULNESS check rather than a
    performance figure: an arm that spends its budget re-measuring the same variant will show
    a unique-measurement count far below the nominal budget.

Usage:
  python paired_stats.py <run_dir> [--reference panoptes_default] [--metric cum_top1_hits_including_initial]
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

METRICS = (
    "cum_top1_hits_including_initial",
    "cum_top1_hits_excluding_initial",
    "best_measured_fitness",
    "n_unique_measurements",
)


def load_final_rounds(run_dir: Path):
    """method -> seed -> {metric: value} at the LAST round each (method, seed) reached."""
    rows = defaultdict(dict)
    with (run_dir / "round_metrics.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            # Every round is kept, whatever its status. Filtering on status == "completed"
            # silently discards rounds an arm ran SHORT on - which is exactly the case that
            # matters here, because an arm whose proposals were all rejected as duplicates
            # reports `short_batch` for every round after the first. Dropping those rounds
            # falls back to round 0 and reports the shared starting sample as the endpoint.
            method, seed, rnd = row["method"], int(row["seed"]), int(row["round"])
            previous = rows[method].get(seed)
            if previous is None or rnd > previous["_round"]:
                record = {"_round": rnd, "_status": row.get("status", "")}
                for metric in METRICS:
                    record[metric] = float(row[metric])
                rows[method][seed] = record
    return rows


def paired_bootstrap_ci(diff: np.ndarray, iterations: int = 20000, seed: int = 0):
    rng = np.random.default_rng(seed)
    draws = rng.choice(diff, size=(iterations, diff.size), replace=True).mean(axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument("--reference", default="panoptes_default")
    parser.add_argument("--metric", default="cum_top1_hits_including_initial",
                        choices=METRICS)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    table = load_final_rounds(run_dir)
    if args.reference not in table:
        raise SystemExit(f"reference method {args.reference!r} not in run; have {sorted(table)}")

    print(f"run: {run_dir.name}")
    print(f"metric: {args.metric}   reference: {args.reference}")
    print(f"unit of independence: seed (paired within seed)\n")

    print(f"{'method':<34}{'n':>4}{'mean':>9}{'sd':>8}{'uniq.meas':>11}{'best fit':>10}"
          f"  {'final-round status'}")
    print("-" * 100)
    for method in sorted(table):
        seeds = sorted(table[method])
        values = np.array([table[method][s][args.metric] for s in seeds])
        uniq = np.array([table[method][s]["n_unique_measurements"] for s in seeds])
        best = np.array([table[method][s]["best_measured_fitness"] for s in seeds])
        statuses = sorted({table[method][s]["_status"] for s in seeds})
        print(f"{method:<34}{len(seeds):>4}{values.mean():>9.4f}{values.std(ddof=1):>8.4f}"
              f"{uniq.mean():>11.2f}{best.mean():>10.4f}  {','.join(statuses)}")

    print(f"\nPaired comparisons against {args.reference} on {args.metric}")
    print("-" * 100)
    print(f"{'method':<34}{'paired n':>9}{'mean diff':>11}{'95% CI':>20}"
          f"{'wilcoxon p':>12}{'sign p':>9}{'W/L/T':>12}")
    ref_seeds = set(table[args.reference])
    for method in sorted(table):
        if method == args.reference:
            continue
        shared = sorted(ref_seeds & set(table[method]))
        if not shared:
            print(f"{method:<34}{'0':>9}  no shared seeds")
            continue
        a = np.array([table[method][s][args.metric] for s in shared])
        b = np.array([table[args.reference][s][args.metric] for s in shared])
        diff = a - b
        lo, hi = paired_bootstrap_ci(diff)
        wins = int((diff > 0).sum())
        losses = int((diff < 0).sum())
        ties = int((diff == 0).sum())
        if np.allclose(diff, 0):
            wilcoxon_p = float("nan")
        else:
            wilcoxon_p = float(stats.wilcoxon(a, b, zero_method="wilcox").pvalue)
        sign_p = (float(stats.binomtest(wins, wins + losses, 0.5).pvalue)
                  if wins + losses else float("nan"))
        print(f"{method:<34}{len(shared):>9}{diff.mean():>11.4f}"
              f"{f'[{lo:+.3f}, {hi:+.3f}]':>20}{wilcoxon_p:>12.3g}{sign_p:>9.3g}"
              f"{f'{wins}/{losses}/{ties}':>12}")

    print("\nNote: 'uniq.meas' is the mean number of DISTINCT variants an arm actually measured "
          "out of the\nnominal budget. A value far below the budget means the arm re-selected "
          "variants it had already\nmeasured, and its endpoint is not a search result.")


if __name__ == "__main__":
    main()
