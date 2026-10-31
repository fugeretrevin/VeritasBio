"""
scripts/analysis/tail_motif_narrowness.py — is a landscape's tail one narrow motif?

Why this exists
---------------
On TrpB, all five variants at or above 0.90 carry the same pair, `K227 + G228`, with
tolerant flanks (183 in {A,C,V}, 184 in {I,L}). Removing k-DPP batch diversity *gained*
2.55 top-1% hits there, and ridge — which concentrates implicitly — reaches the ≥0.90
tier in about 17% of campaigns where the diversity-driven engine reaches it in 0%.

Those three observations have one candidate explanation: **diversity is the wrong
objective when the tail is a single narrow motif.** A diversity term spends budget
spreading away from the one productive region.

That explanation makes a prediction, and this script measures the quantity the
prediction is about, so the prediction can be registered *before* any arm is run:

    k-DPP diversity helps where a landscape's tail is multi-modal, and hurts where it
    is a narrow motif.

Nothing here reads an arm, a baseline or a structure, so it runs on every landscape
with a prepared table — including the ones where ALDE cannot function and the ones
with no resolved structure.

The statistic
-------------
For the top-k variants by fitness, at each library position take the frequency of the
most common residue. `concentration` is the mean of those frequencies over positions:
1.0 means every top variant shares a residue at every position, and a low value means
the tail is spread.

Raw concentration is not comparable across landscapes — a library with few observed
residues at a position concentrates by construction. So the observed value is compared
against a **within-landscape null**: fitness is permuted, holding the genotype set and
therefore every marginal residue frequency fixed, and the statistic is recomputed. The
reported `z` and `p` are against that null, which is what makes a 4-site and a 16-site
library comparable.

`conserved_positions` counts positions where the top-k agree at or above
`--conserved-at` (default 0.80); `motif` names them. For TrpB that should recover
`227K` and `228G`.

Genotypes are canonicalised with `veritas.calibration.genotypes.canonical_combo`
(handoff rule 2.1): mutation strings omit wild-type positions, and comparing them
directly lets a 4-tuple truncate to 3.

Usage
-----
    python scripts/analysis/tail_motif_narrowness.py --all
    python scripts/analysis/tail_motif_narrowness.py --library trpb --top-k 5
    python scripts/analysis/tail_motif_narrowness.py --all --json out.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from veritas.calibration.genotypes import canonical_combo  # noqa: E402
from veritas.calibration.landscapes import REGISTRY  # noqa: E402

# The tail is a QUANTILE, floored and capped, not a fixed count: 5 of 1,886 is the
# top 0.27% while 5 of 165,427 is the top 0.003%, so a fixed count compares different
# objects across landscapes. Mirrors the --qualify-threshold -> --qualify-quantile
# lesson, where an absolute bar sat below the global minimum on four landscapes.
DEFAULT_TAIL_QUANTILE = 0.0001   # top 0.01%
DEFAULT_TAIL_FLOOR = 5
DEFAULT_TAIL_CAP = 50
DEFAULT_PERMUTATIONS = 2000
DEFAULT_CONSERVED_AT = 0.80


def _library_spec(key: str) -> Tuple[Sequence[int], str]:
    """(positions, parent_residues) for a registered landscape."""
    spec = REGISTRY[key]
    positions = getattr(spec, "positions", None) or getattr(spec, "sites", None)
    parent = (getattr(spec, "parent_residues", None)
              or getattr(spec, "parent", None)
              or getattr(spec, "wt_residues", None))
    if positions is None or parent is None:
        raise KeyError(
            f"{key}: registry entry exposes no positions/parent_residues "
            f"(has: {sorted(n for n in dir(spec) if not n.startswith('_'))})"
        )
    return list(positions), str(parent)


def load_combos(key: str) -> Tuple[List[str], List[float], Sequence[int], str]:
    """Canonical combos and fitness for a landscape, with unparseable rows reported."""
    spec = REGISTRY[key]
    table = Path(spec.path)
    positions, parent = _library_spec(key)

    combos: List[str] = []
    scores: List[float] = []
    skipped = 0
    with table.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            combo = canonical_combo(row["mutant"], positions, parent)
            if combo is None:
                skipped += 1
                continue
            combo = "".join(combo)
            # Rule 2.1: never compare genotypes of unequal length.
            if len(combo) != len(positions):
                raise AssertionError(
                    f"{key}: canonical combo {combo!r} is {len(combo)} long for "
                    f"{len(positions)} positions"
                )
            combos.append(combo)
            try:
                scores.append(float(row["DMS_score"]))
            except (KeyError, ValueError):
                combos.pop()
                skipped += 1
    if skipped:
        print(f"  note: {key}: {skipped:,} row(s) skipped as unparseable", file=sys.stderr)
    return combos, scores, positions, parent


def concentration(combos: Sequence[str]) -> Tuple[float, List[float]]:
    """Mean over positions of the most-common-residue frequency among `combos`."""
    if not combos:
        return 0.0, []
    n_pos = len(combos[0])
    per_position: List[float] = []
    for i in range(n_pos):
        column = [c[i] for c in combos]
        top = max(column.count(aa) for aa in set(column))
        per_position.append(top / len(column))
    return statistics.fmean(per_position), per_position


def resolve_tail_size(n: int, quantile: float, floor: int, cap: int) -> Tuple[int, str]:
    """(tail size, which bound decided it) - recorded so a bound is never silent."""
    raw_k = int(math.floor(n * quantile))
    if raw_k < floor:
        return min(floor, n), "floor"
    if raw_k > cap:
        return cap, "cap"
    return raw_k, "quantile"


def analyse(key: str, quantile: float, floor: int, cap: int, permutations: int,
            conserved_at: float, seed: int = 0) -> Dict:
    combos, scores, positions, parent = load_combos(key)
    top_k, bound = resolve_tail_size(len(combos), quantile, floor, cap)
    if len(combos) < top_k * 4:
        return {"library": key, "status": "too_small", "n": len(combos)}

    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    tail_idx = order[:top_k]
    tail = [combos[i] for i in tail_idx]
    observed, per_pos = concentration(tail)

    # Within-landscape null: permute fitness, holding the genotype set - and therefore
    # every marginal residue frequency - fixed. A library that concentrates by
    # construction concentrates in the null too, which is the point.
    rng = random.Random(seed)
    null: List[float] = []
    idx = list(range(len(combos)))
    for _ in range(permutations):
        picked = rng.sample(idx, top_k)
        null.append(concentration([combos[i] for i in picked])[0])
    mu = statistics.fmean(null)
    sd = statistics.pstdev(null) or 1e-12
    z = (observed - mu) / sd
    # One-sided: we are asking whether the tail is MORE concentrated than chance.
    p = (sum(1 for v in null if v >= observed) + 1) / (permutations + 1)

    conserved = [
        (positions[i], max(set(col), key=col.count))
        for i, col in ((i, [c[i] for c in tail]) for i in range(len(tail[0])))
        if per_pos[i] >= conserved_at
    ]
    motif = ":".join(f"{pos}{aa}" for pos, aa in conserved) or "(none)"

    return {
        "library": key,
        "status": "ok",
        "n_variants": len(combos),
        "n_positions": len(positions),
        "tail_size": top_k,
        "tail_bound": bound,
        "tail_quantile": quantile,
        "tail_fitness_min": round(min(scores[i] for i in tail_idx), 4),
        "concentration": round(observed, 4),
        "null_mean": round(mu, 4),
        "z": round(z, 2),
        "p": round(p, 5),
        "conserved_positions": len(conserved),
        "motif": motif,
        "per_position": [round(v, 3) for v in per_pos],
        # The prediction this feeds: narrow tail -> diversity should hurt.
        "predicted_dpp_effect": "hurts" if z >= 3.0 else ("helps" if z <= 1.0 else "ambiguous"),
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", action="append", default=None,
                    help="landscape key; repeatable. Default: --all")
    ap.add_argument("--all", action="store_true", help="every registered landscape")
    ap.add_argument("--tail-quantile", type=float, default=DEFAULT_TAIL_QUANTILE,
                    help="tail = this top fraction of the eligible set (default 0.01%%)")
    ap.add_argument("--tail-floor", type=int, default=DEFAULT_TAIL_FLOOR)
    ap.add_argument("--tail-cap", type=int, default=DEFAULT_TAIL_CAP)
    ap.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    ap.add_argument("--conserved-at", type=float, default=DEFAULT_CONSERVED_AT)
    ap.add_argument("--json", type=Path, default=None, help="also write results here")
    args = ap.parse_args(argv)

    keys = sorted(REGISTRY) if (args.all or not args.library) else list(args.library)
    unknown = [k for k in keys if k not in REGISTRY]
    if unknown:
        print(f"unknown landscape(s) {unknown}; registered: {sorted(REGISTRY)}",
              file=sys.stderr)
        return 2

    rows = []
    for key in keys:
        try:
            rows.append(analyse(key, args.tail_quantile, args.tail_floor,
                                args.tail_cap, args.permutations, args.conserved_at))
        except Exception as exc:  # report and continue: one bad table must not hide the rest
            rows.append({"library": key, "status": f"error: {type(exc).__name__}: {exc}"})

    hdr = (f"{'landscape':<15}{'n':>9}{'pos':>5}{'tail':>6}{'bound':>10}"
           f"{'conc':>7}{'null':>7}{'z':>8}{'p':>9}{'cons':>6}  motif / status")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r.get("status") != "ok":
            print(f"{r['library']:<15}{'':>53}  {r.get('status')}")
            continue
        print(f"{r['library']:<15}{r['n_variants']:>9,}{r['n_positions']:>5}"
              f"{r['tail_size']:>6}{r['tail_bound']:>10}{r['concentration']:>7.3f}"
              f"{r['null_mean']:>7.3f}{r['z']:>8.2f}{r['p']:>9.5f}"
              f"{r['conserved_positions']:>6}  {r['motif']}")

    print()
    print(f"concentration = mean over positions of the most-common-residue frequency in")
    print(f"the tail (top {args.tail_quantile:.2%}, floor {args.tail_floor}, cap "
          f"{args.tail_cap}). z and p are against a within-landscape permutation null")
    print(f"that holds the genotype set fixed, so libraries of different width compare.")
    print(f"PREDICTION (register before running arms): k-DPP diversity HURTS where z is")
    print(f"high (narrow tail) and HELPS where z is low (multi-modal tail).")

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
