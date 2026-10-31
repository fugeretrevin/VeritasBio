"""
scripts/analysis/shared_substitution_enrichment.py — tail sharing in a variable-depth DMS.

Why this is separate from `tail_motif_narrowness.py`
---------------------------------------------------
That script measures, for each library position, how often the top variants agree on a
residue. It assumes a **combinatorial library**: every variant carries a residue at
every one of a small fixed set of positions, so "agreement at position i" is defined for
all of them.

A random-mutagenesis DMS does not have that shape. avGFP (Sarkisyan et al. 2016) is
51,714 variants spanning **233 distinct mutated positions** at depths 1 through 8+, so
two variants usually touch disjoint position sets and per-position agreement is
undefined for most pairs. Running the combinatorial instrument on it would produce a
number, and the number would not mean what its name says.

The analogous question, asked correctly: **do the best variants share specific
substitutions more than chance?** For the top-k variants, count how often each
(position, residue) substitution recurs, and compare against a permutation null that
samples k variants from the same pool — which holds the depth distribution and the
marginal substitution frequencies fixed.

Not comparable to the combinatorial z
-------------------------------------
`tail_motif_narrowness.py` z values and the `enrichment_z` here measure different
quantities on different supports. A landscape screened by this script cannot be slotted
into the band table in `docs/preregistration-motif-narrowness.md` without a separate
registration, and this script does not assign bands.

Usage
-----
    python scripts/analysis/shared_substitution_enrichment.py --dataset <csv> [--top-k N]
"""
from __future__ import annotations

import argparse
import csv
import random
import re
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_MUT = re.compile(r"^([A-Z])(\d+)([A-Z*])$")

DEFAULT_TOP_K = 50
DEFAULT_PERMUTATIONS = 2000


def load(csv_path: Path) -> Tuple[List[Tuple[Tuple[str, ...], float]], int]:
    """[(substitution tuple, score)], plus a count of unparseable rows."""
    out: List[Tuple[Tuple[str, ...], float]] = []
    skipped = 0
    with csv_path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            toks = [t for t in str(row.get("mutant", "")).split(":") if t]
            if not toks or not all(_MUT.match(t) for t in toks):
                skipped += 1
                continue
            try:
                score = float(row["DMS_score"])
            except (KeyError, ValueError):
                skipped += 1
                continue
            out.append((tuple(sorted(toks)), score))
    return out, skipped


def sharing(variants: Sequence[Tuple[str, ...]]) -> Tuple[float, float, Counter]:
    """(max share, mean share of the top 5 substitutions, the counter itself).

    `share` is the fraction of the sampled variants carrying a given substitution.
    """
    if not variants:
        return 0.0, 0.0, Counter()
    counts: Counter = Counter()
    for subs in variants:
        counts.update(subs)
    n = len(variants)
    top = [c / n for _, c in counts.most_common(5)]
    return (top[0] if top else 0.0), (statistics.fmean(top) if top else 0.0), counts


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Do the best variants of a variable-depth DMS share substitutions?")
    ap.add_argument("--dataset", type=Path, required=True)
    ap.add_argument("--top-k", type=int, default=DEFAULT_TOP_K)
    ap.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    if not args.dataset.is_file():
        print(f"dataset not found: {args.dataset}", file=sys.stderr)
        return 2

    data, skipped = load(args.dataset)
    if len(data) < args.top_k * 4:
        print(f"only {len(data)} parseable variants; need at least {args.top_k * 4}",
              file=sys.stderr)
        return 1

    depths = Counter(len(s) for s, _ in data)
    positions = {int(_MUT.match(t).group(2)) for s, _ in data for t in s}

    print(f"dataset   : {args.dataset.name}")
    print(f"variants  : {len(data):,} parseable" + (f", {skipped:,} skipped" if skipped else ""))
    print(f"positions : {len(positions)} distinct mutated")
    print(f"depths    : {dict(sorted(depths.items())[:8])}")
    print(f"top-k     : {args.top_k}")
    print()

    data.sort(key=lambda kv: kv[1], reverse=True)
    tail = [s for s, _ in data[: args.top_k]]
    obs_max, obs_mean, counts = sharing(tail)

    rng = random.Random(args.seed)
    pool = [s for s, _ in data]
    null_max: List[float] = []
    null_mean: List[float] = []
    for _ in range(args.permutations):
        sample = rng.sample(pool, args.top_k)
        mx, mn, _ = sharing(sample)
        null_max.append(mx)
        null_mean.append(mn)

    def z(obs: float, null: Sequence[float]) -> Tuple[float, float]:
        mu = statistics.fmean(null)
        sd = statistics.pstdev(null) or 1e-12
        p = (sum(1 for v in null if v >= obs) + 1) / (len(null) + 1)
        return (obs - mu) / sd, p

    zmax, pmax = z(obs_max, null_max)
    zmean, pmean = z(obs_mean, null_mean)

    print(f"{'statistic':<34}{'observed':>10}{'null mean':>11}{'z':>9}{'p':>9}")
    print("-" * 73)
    print(f"{'share of most common substitution':<34}{obs_max:>10.3f}"
          f"{statistics.fmean(null_max):>11.3f}{zmax:>9.2f}{pmax:>9.5f}")
    print(f"{'mean share of top 5 substitutions':<34}{obs_mean:>10.3f}"
          f"{statistics.fmean(null_mean):>11.3f}{zmean:>9.2f}{pmean:>9.5f}")

    print()
    print(f"most shared substitutions in the top {args.top_k}:")
    for sub, c in counts.most_common(8):
        print(f"   {sub:<10} {c:>4}/{args.top_k}  ({c / args.top_k:.0%})")

    print()
    verdict = ("CONCENTRATED - the tail shares substitutions far above chance"
               if zmean >= 3.0 else
               "DIFFUSE - tail sharing is not distinguishable from chance"
               if zmean <= 1.0 else "AMBIGUOUS")
    print(f"verdict: {verdict}")
    print()
    print("NOT comparable to tail_motif_narrowness.py z values: different statistic, "
          "different support.")
    print("A landscape screened here needs its own registration before entering any band "
          "table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
