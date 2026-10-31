"""Did every arm actually get the same measurement budget?

The integrity agent noted in passing that a default GB1 run finishes PARTIAL because
GB1 covers 149,360 of the 160,000 possible 4-site combinations, so some generated
proposals are not in the eligible universe, are rejected, and are never backfilled.

If that makes an arm's batches come up SHORT, it has measured fewer than 125 variants
and the comparison is not budget-matched. The review is explicit that every method must
use exactly 125 unique measurements, so this has to be checked rather than assumed.
"""
from __future__ import annotations

import csv
import statistics as st
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")

ARMS: List[Tuple[str, str, str]] = [
    ("gb1 random", "run_20261031_024951_7360aa94", "random"),
    ("gb1 ridge", "run_20261031_024951_7360aa94", "ridge"),
    ("gb1 random_forest", "run_20261031_024951_7360aa94", "random_forest"),
    ("gb1 alde", "run_20261031_025645_7033ab05", "alde"),
    ("gb1 reliability", "run_20261031_024308_86773bff", "panoptes_default"),
    ("gb1 rel-dkl,gfn", "run_20261031_024452_429a9ee4", "panoptes_default"),
    ("gb1 frontier", "run_20261031_025342_33e5cac0", "panoptes_qualify_then_optimize"),
    ("trpb ridge", "run_20261031_025252_7dce6099", "ridge"),
    ("trpb reliability", "run_20261031_024533_c41a0125", "panoptes_default"),
    ("trpb frontier", "run_20261031_025827_8dfabb93", "panoptes_qualify_then_optimize"),
]

print("=" * 86)
print("MEASUREMENT-BUDGET FAIRNESS CHECK")
print("=" * 86)
print(f"  {'arm':<20}{'mean final n':>14}{'min':>7}{'max':>7}{'seeds':>7}{'short of 125':>14}")

for label, run, method in ARMS:
    p = ROOT / run / "round_metrics.csv"
    if not p.exists():
        print(f"  {label:<20}  (missing {run})")
        continue
    last: Dict[int, Tuple[int, int]] = {}
    with open(p, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            try:
                seed = int(float(r["seed"]))
                rnd = int(float(r["round"]))
                n = int(float(r["n_unique_measurements"]))
            except (TypeError, ValueError, KeyError):
                continue
            if seed not in last or rnd >= last[seed][0]:
                last[seed] = (rnd, n)
    vals = [n for _, n in last.values()]
    if not vals:
        continue
    mean = st.fmean(vals)
    print(f"  {label:<20}{mean:>14.2f}{min(vals):>7}{max(vals):>7}{len(vals):>7}"
          f"{125 - mean:>13.2f}")

print()
print("  A budget-matched comparison needs every arm at 125. Any arm materially below it")
print("  measured fewer variants than its competitors and its hit count is not comparable.")
