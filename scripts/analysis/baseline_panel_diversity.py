"""Diversity-adjusted yield for the baseline panel, by diversity.py's own rule.

A raw top-1% count treats N variants in one mutational neighbourhood as N wins. That
matters here specifically because `greedy_measured` is a LOCAL EXPLOIT rule - it is
designed to pick the neighbourhood of the best measured variant - so if it beats the
engine on raw count while returning one dense basin, the raw comparison overstates it.
Checking rather than assuming, because the same argument was asserted without a test in
docs/545.md section 4 and was wrong there.

Rule: sort a seed's hits by fitness descending and greedily keep one only when it is at
least d away from everything already kept. Applied AFTER measurement, so it cannot
advantage a method that knew about it.

Rule 2.1: variant_id in selection_trace.csv is a MUTATION STRING that OMITS wild-type
positions, so every genotype is expanded to full library length with
veritas.calibration.genotypes.canonical_combo and compared with its `hamming`, which
raises on unequal lengths. Nothing is zipped or truncated here.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# The source tree is THIS checkout - the code under test - while run artifacts live in the
# gitignored results root of the REAL repository, which a worktree cannot see (rule 2.8).
# Resolved through veritas.paths rather than hard-coded, and every failure names the paths
# it checked.
CHECKOUT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CHECKOUT))

from veritas import paths as vpaths  # noqa: E402
from veritas.calibration.genotypes import canonical_combo, hamming  # noqa: E402
import benchmark_panoptes_isolated as bench  # noqa: E402

RESULTS_ROOT = vpaths.results_dir("benchmark_comparison_results")
DISTANCES = (2, 3)


def threshold_of(run_dir: Path) -> float:
    """Rule 2.3: the bar comes from the run's OWN manifest, never recomputed."""
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    # results.top1.threshold_value, NOT configuration.components.qualify_threshold. Those
    # are two different bars: qualify_* is the 99.96th-percentile ONLINE qualification bar
    # (9.6878 here, 27 qualifying variants), while the scoring endpoint is the top-1% tier
    # (9.6087, 651 variants). Reading the qualification bar here would silently rescore the
    # whole panel against a 27-variant tier.
    try:
        value = manifest["results"]["top1"]["threshold_value"]
    except (KeyError, TypeError):
        raise SystemExit(
            f"no results.top1.threshold_value in {run_dir / 'manifest.json'}"
        )
    return float(value)


def paired_t(diff: np.ndarray) -> float:
    if diff.size < 2:
        return float("nan")
    sd = float(diff.std(ddof=1))
    if sd == 0.0:
        return float("nan")
    return float(diff.mean() / (sd / math.sqrt(diff.size)))


def main(run_name: str, library_key: str) -> int:
    run_dir = RESULTS_ROOT / run_name
    if not run_dir.is_dir():
        raise SystemExit(f"not a directory: {run_dir}")
    library = bench.resolve_library(library_key)
    positions, parent = library.positions, library.parent_residues
    bar = threshold_of(run_dir)

    # method -> seed -> list of (fitness, combo), hits only, initial sample EXCLUDED so the
    # metric matches cum_top1_hits_excluding_initial.
    hits = defaultdict(lambda: defaultdict(list))
    unparsed = 0
    with (run_dir / "selection_trace.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") == "initial":
                continue
            raw = row.get("revealed_fitness")
            if raw in (None, ""):
                continue
            fitness = float(raw)
            if fitness < bar:
                continue
            combo = canonical_combo(row["variant_id"], positions, parent)
            if combo is None:
                unparsed += 1
                continue
            hits[row["method"]][int(row["seed"])].append((fitness, combo))

    print(f"{library_key}  {run_name}   top-1% bar {bar!r} (from the run's own manifest)")
    print(f"library: {len(positions)} sites; genotypes expanded to full length, WT filled")
    if unparsed:
        print(f"!! {unparsed} variant_id(s) could not be canonicalised and were EXCLUDED")

    def adjusted(entries, d: int) -> int:
        kept = []
        for _, combo in sorted(entries, key=lambda e: -e[0]):
            if all(hamming(combo, other) >= d for other in kept):
                kept.append(combo)
        return len(kept)

    # The seed block comes from round_metrics.csv, NOT from the seeds that happen to appear
    # in `hits`. A seed on which an arm found zero hits has no row in `hits[arm]`, so
    # intersecting the arms' hit-seed sets silently drops every seed where the weakest arm
    # scored nothing - which is most of them for `random`. That conditions the whole table
    # on success and inflates every mean: `random` read 1.60 instead of 0.875. A zero-hit
    # seed is a measurement of zero, not a missing measurement.
    arms = [a for a in bench.BASELINE_PANEL_METHODS if a in hits]
    seeds = set()
    with (run_dir / "round_metrics.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            seeds.add(int(float(row["seed"])))
    seeds = sorted(seeds)
    for arm in arms:
        for seed in seeds:
            hits[arm].setdefault(seed, [])
    print(f"arms {len(arms)}, seeds paired {len(seeds)} (zero-hit seeds counted as 0)\n")

    header = f"{'arm':<20}{'raw':>8}"
    for d in DISTANCES:
        header += f"{'d>=' + str(d):>9}{'redund':>9}"
    print(header)
    print("-" * len(header))
    table = {}
    for arm in arms:
        raw = np.array([len(hits[arm][s]) for s in seeds], dtype=float)
        line = f"{arm:<20}{raw.mean():>8.2f}"
        table[arm] = {"raw": raw}
        for d in DISTANCES:
            adj = np.array([adjusted(hits[arm][s], d) for s in seeds], dtype=float)
            table[arm][d] = adj
            redund = raw.mean() / adj.mean() if adj.mean() else float("nan")
            line += f"{adj.mean():>9.2f}{redund:>8.2f}x"
        print(line)

    for reference in ("ridge", "panoptes_default"):
        if reference not in table:
            continue
        print(f"\npaired against {reference}: mean diff, t, better/worse/tied")
        head = f"{'arm':<20}{'raw':>22}"
        for d in DISTANCES:
            head += f"{'d>=' + str(d):>22}"
        print(head)
        print("-" * len(head))
        for arm in arms:
            if arm == reference:
                continue
            line = f"{arm:<20}"
            for key in ("raw",) + DISTANCES:
                diff = table[arm][key] - table[reference][key]
                b, w, t_ = int((diff > 0).sum()), int((diff < 0).sum()), int((diff == 0).sum())
                line += f"{f'{diff.mean():+.2f} t{paired_t(diff):+.1f} {b}/{w}/{t_}':>22}"
            print(line)
    print("\nredund = raw hits per diversity-adjusted hit; 1.00x means every hit is distinct.")
    print("Hamming separation is a sequence-level proxy, not functional independence.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
