"""Pre-registered landscape predictors, computed BEFORE the mode comparison exists.

Two questions, two predictors, fixed in advance so neither can be chosen after the fact:

  Q1  Does LOCAL candidate proposal beat ranking the enumerated POOL?
      Predictor: rho_1, the correlation between a variant's fitness and its Hamming-1
      neighbours'. Local expansion carries information exactly when this is high. On
      eqFP611-red the generator proposed candidates averaging 0.1940 fitness against a
      library mean of 0.2316 and found 0 of 82 hits where a random subset of the same
      size would have contained 3.01 - the signature of local expansion on a landscape
      whose neighbourhoods do not predict.

  Q2  Do our prior + Pareto + GP beat ridge?
      Predictor: additivity, the R^2 of a purely additive (one-hot linear) model fit on
      the WHOLE library. Ridge is a linear model: where additivity is high it is already
      near-optimal and nothing can beat it; room exists only where epistasis breaks it.

Both are computed from the full landscape, which a real campaign would not have. The
prospective analogue - the same quantity estimated from the first round's 25 measurements
- is reported alongside, because that is what a customer could actually compute before
deciding. If the two disagree badly, the rule is not deployable even if it is true.
"""
import csv
import random
import statistics
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import benchmark_panoptes_isolated as bench  # noqa: E402

LANDSCAPES = ["trpb", "gb1", "her2", "creilov", "eqfp611_blue", "eqfp611_red",
              "cr9114_h1", "cr6261_h1", "cr6261_h9"]
RNG = random.Random(20261031)
N_PAIRS = 60_000
INITIAL_N = 25


def load(spec):
    """Returns {combo_string: fitness}, combos in declared-site order."""
    wt = dict(zip(spec.positions, spec.parent_residues))
    out = {}
    for row in csv.DictReader(open(spec.dataset)):
        try:
            score = float(row["DMS_score"])
        except (TypeError, ValueError):
            continue
        # HER2 rule 2.3: rows at exactly 0.0 are absent labels, not measurements of zero.
        if spec.key == "her2" and score == 0.0:
            continue
        combo = dict(wt)
        ok = True
        mutant = row["mutant"].strip()
        # GB1 carries a literal "WT" row for the wild type; its combo is the parent.
        if mutant.upper() == "WT":
            out["".join(combo[p] for p in spec.positions)] = score
            continue
        for sub in mutant.split(":"):
            sub = sub.strip()
            if not sub:
                continue
            pos = int(sub[1:-1])
            if pos not in combo:
                ok = False
                break
            combo[pos] = sub[-1]
        if ok:
            out["".join(combo[p] for p in spec.positions)] = score
    return out


def rho_one(table, spec, n_pairs=N_PAIRS):
    """Correlation of fitness across Hamming-1 pairs present in the table."""
    keys = list(table)
    alphabets = [spec.alphabet_for(s) for s in range(spec.n_sites)]
    a, b = [], []
    tries = 0
    while len(a) < n_pairs and tries < n_pairs * 40:
        tries += 1
        base = keys[RNG.randrange(len(keys))]
        slot = RNG.randrange(spec.n_sites)
        options = [c for c in alphabets[slot] if c != base[slot]]
        if not options:
            continue
        neighbour = base[:slot] + RNG.choice(options) + base[slot + 1:]
        other = table.get(neighbour)
        if other is None:
            continue
        a.append(table[base])
        b.append(other)
    if len(a) < 100:
        return float("nan"), len(a)
    return float(np.corrcoef(a, b)[0, 1]), len(a)


def additivity(table, spec, subset=None):
    """R^2 of a one-hot additive model.

    Full-library: ordinary least squares, n >> p, so R^2 is the natural "how much of this
    landscape is additive" measure.

    Prospective (`subset`): 25 measurements against 32-200 one-hot columns is n < p, where
    OLS is not identifiable - which is precisely why the benchmark's own baseline uses
    Ridge(alpha=1.0). The same estimator is used here, fit on the subset and scored
    OUT-OF-SAMPLE on the whole library, since in-sample R^2 on 25 points with a wide
    basis is near 1 by construction and would measure nothing.
    """
    from sklearn.linear_model import LinearRegression, Ridge

    combos = list(subset if subset is not None else table)
    X = bench._one_hot_matrix(combos, spec)
    y = np.array([table[c] for c in combos], dtype=np.float64)

    if subset is None:
        if len(combos) <= X.shape[1] + 1:
            return float("nan")
        return float(LinearRegression().fit(X, y).score(X, y))

    model = Ridge(alpha=1.0).fit(X, y)
    allc = list(table)
    return float(model.score(bench._one_hot_matrix(allc, spec),
                             np.array([table[c] for c in allc])))


print(f"{'landscape':<14}{'n':>8}{'rho_1':>9}{'pairs':>8}"
      f"{'additivity':>12}{'add@25':>9}{'cliff':>8}")
print("-" * 70)
rows = []
for key in LANDSCAPES:
    spec = bench.LIBRARIES[key]
    table = load(spec)
    r1, npairs = rho_one(table, spec)
    add_full = additivity(table, spec)

    # Prospective estimate: fit on 25 random variants, score on the whole library.
    # Median over 9 draws, because a single 25-point draw is noisy.
    keys = list(table)
    estimates = []
    for trial in range(9):
        RNG.shuffle(keys)
        estimates.append(additivity(table, spec, subset=keys[:INITIAL_N]))
    usable = [e for e in estimates if e == e]
    add_25 = statistics.median(usable) if usable else float("nan")

    values = list(table.values())
    k = max(1, round(len(values) * 0.01))
    thresh = sorted(values, reverse=True)[k - 1]
    hits = [v for v in values if v >= thresh]
    cliff = statistics.fmean(hits) / statistics.fmean(values)

    print(f"{key:<14}{len(table):>8}{r1:>9.3f}{npairs:>8}"
          f"{add_full:>12.3f}{add_25:>9.3f}{cliff:>8.2f}")
    rows.append((key, r1, add_full, add_25, cliff))

print()
print("Pre-registered readings:")
print("  Q1  local proposal should beat pool ranking where rho_1 is HIGH and lose where")
print("      it is LOW. eqFP611-red already lost 0 to 15.35, so its rho_1 should be")
print("      among the lowest here; if it is not, the predictor is wrong.")
print("  Q2  our margin over ridge should be LARGEST where additivity is LOWEST.")
print("      CR9114 beat ridge 5.00 to 1.83, so its additivity should be low.")
print()
print("  n = 9 landscapes, and CR6261-H1/H9 are one protein against two antigens, so the")
print("  effective n is 8. A correlation needs |r| > 0.71 for p < 0.05 at n = 8. This is")
print("  hypothesis-generating unless an effect is very large, and must be reported so.")
