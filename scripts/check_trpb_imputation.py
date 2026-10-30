"""
scripts/check_trpb_imputation.py - can the TrpB tail variants be KNN-imputed?

Johnston et al. (2024) report that 871 of 160,000 fitness values were imputed with
the scikit-learn KNN imputer, and the prepared CSV in this repository carries no
flag distinguishing them. Every conclusion about the >=0.90 tail - the K227/G228
obligate pair, the two-hop path to the optimum, the qualify-then-optimize result -
rests on those five variants being measurements rather than interpolations.

The upstream flags are not available offline, but the question is answerable from
the landscape itself without them.

Argument
--------
KNN imputation replaces a missing value with a (possibly weighted) MEAN of its k
nearest neighbours. A mean cannot exceed the maximum of the values averaged, so the
largest value a k-neighbour imputation can produce is the mean of the k HIGHEST
neighbours. Any variant whose fitness strictly exceeds that ceiling cannot have been
imputed, and is therefore measured.

In this combinatorially complete four-site library the nearest neighbours in one-hot
space are the 76 Hamming-1 variants, and they are all exactly equidistant. Two
consequences:

  * distance weighting degenerates to uniform weighting, so the ceiling above is the
    correct bound for either weighting mode.
  * sklearn breaks the 76-way tie by index order, which is arbitrary with respect to
    fitness. The EXPECTED imputed value is therefore the plain neighbour mean - for
    these variants about 0.32, not 0.91. Reaching 0.91 would require the tie-break to
    have selected precisely the few highest neighbours out of 76.

Assumption, stated because it is load-bearing: that the imputation was performed over
this variant space with sequence-derived features, so that Hamming-1 variants are the
nearest neighbours. If the authors imputed in some other feature space the neighbour
set could differ, and this test would need redoing against that space.

Usage:
    python scripts/check_trpb_imputation.py [--top N] [--k 1,3,5,10]
"""
from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_CSV = ROOT / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv"

SITES = (183, 184, 227, 228)
# Tm9D8* parent residues at the four saturated sites. Cross-checked against every
# notation's wild-type prefix during parsing rather than trusted as a constant.
PARENT = {183: "V", 184: "F", 227: "V", 228: "S"}
_MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")

IMPUTED_COUNT_REPORTED = 871
HAMMING1_NEIGHBOURS = 76  # 4 positions x 19 alternative residues


def load_landscape(csv_path: pathlib.Path) -> dict:
    """Parse the CSV into {4-residue combo: fitness}, verifying the notation."""
    fit = {}
    with csv_path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            combo = dict(PARENT)
            for tok in row["mutant"].strip().split(":"):
                if not tok:
                    continue
                m = _MUT.match(tok)
                if m is None:
                    raise ValueError(f"unparseable mutation {tok!r}")
                wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
                if pos not in PARENT:
                    raise ValueError(f"unexpected position {pos} in {row['mutant']!r}")
                if wt != PARENT[pos]:
                    # A numbering or parent mismatch would silently corrupt every
                    # neighbour set, so fail loudly instead of scoring nonsense.
                    raise ValueError(
                        f"notation {tok!r} disagrees with parent {PARENT[pos]}{pos}"
                    )
                combo[pos] = mut
            fit["".join(combo[p] for p in SITES)] = float(row["DMS_score"])
    return fit


def neighbour_values(combo: str, fit: dict, alphabet: list) -> list:
    """Fitness of every Hamming-1 neighbour present in the landscape."""
    out = []
    for i in range(len(SITES)):
        for aa in alphabet:
            if aa != combo[i]:
                nb = combo[:i] + aa + combo[i + 1:]
                if nb in fit:
                    out.append(fit[nb])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Test whether the highest-fitness TrpB variants could be KNN-imputed.",
    )
    ap.add_argument("--dataset", type=pathlib.Path, default=DEFAULT_CSV)
    ap.add_argument("--top", type=int, default=5,
                    help="how many of the highest-fitness variants to test (default 5)")
    ap.add_argument("--k", default="1,3,5,10",
                    help="comma-separated neighbour counts to bound (default 1,3,5,10)")
    args = ap.parse_args(argv)

    ks = [int(x) for x in args.k.split(",") if x.strip()]
    if not ks:
        print("--k must name at least one neighbour count", file=sys.stderr)
        return 2
    if not args.dataset.is_file():
        print(f"dataset not found: {args.dataset}", file=sys.stderr)
        return 2

    fit = load_landscape(args.dataset)
    alphabet = sorted({c for combo in fit for c in combo})
    space = len(alphabet) ** len(SITES)
    parent_combo = "".join(PARENT[p] for p in SITES)

    print(f"dataset   : {args.dataset.name}")
    print(f"variants  : {len(fit):,} of {space:,} ({space - len(fit):,} absent)")
    print(f"alphabet  : {len(alphabet)} residues at positions {SITES}")
    print(f"parent    : {parent_combo} (present in file: {parent_combo in fit})")
    print(f"reported imputed: {IMPUTED_COUNT_REPORTED:,} values "
          f"({100 * IMPUTED_COUNT_REPORTED / space:.2f}% of the space), unflagged")
    print()

    ranked = sorted(fit.items(), key=lambda kv: kv[1], reverse=True)[: args.top]
    kcols = "".join(f"{'k=' + str(k):>9}" for k in ks)
    print("Highest-fitness variants vs the ceiling on a KNN-imputed value")
    print(f"{'combo':<7}{'fitness':>9}{kcols}{'E[imp]':>9}  verdict")
    print("-" * (25 + 9 * len(ks) + 32))

    default_k = 5
    ruled_out_at_default = 0
    survivors = []
    for combo, f in ranked:
        vals = sorted(neighbour_values(combo, fit, alphabet), reverse=True)
        ceilings = {k: statistics.fmean(vals[:k]) for k in ks}
        expected = statistics.fmean(vals)
        feasible = [k for k in ks if f <= ceilings[k] + 1e-9]
        if default_k in ks:
            if f > ceilings[default_k] + 1e-9:
                ruled_out_at_default += 1
            else:
                survivors.append(combo)
        verdict = ("RULED OUT for every k tested" if not feasible
                   else "possible only at k=" + ",".join(str(k) for k in feasible))
        cells = "".join(f"{ceilings[k]:>9.4f}" for k in ks)
        print(f"{combo:<7}{f:>9.4f}{cells}{expected:>9.4f}  {verdict}")

    # Landscape-wide context: how many variants this test proves are measured.
    exceed = 0
    for combo, f in fit.items():
        vals = neighbour_values(combo, fit, alphabet)
        if vals and f > max(vals):
            exceed += 1

    print()
    if default_k in ks:
        print(f"At sklearn's default n_neighbors={default_k}: "
              f"{ruled_out_at_default} of {len(ranked)} ruled out.")
        if survivors:
            odds = math.comb(HAMMING1_NEIGHBOURS, default_k)
            print(f"Survivor(s) {', '.join(survivors)} additionally require the "
                  f"{HAMMING1_NEIGHBOURS}-way tie-break to have\nselected exactly the "
                  f"{default_k} highest neighbours - 1 in {odds:,} under arbitrary ordering.")
    print()
    print(f"landscape-wide: {exceed:,}/{len(fit):,} ({100 * exceed / len(fit):.2f}%) variants "
          f"exceed their neighbour\nmaximum, so cannot be KNN averages at any k; those are "
          f"measured values.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
