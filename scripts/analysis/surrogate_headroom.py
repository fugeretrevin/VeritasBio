"""Will enumerate mode actually help? Decided offline, before implementing anything.

Enumerate mode hands the engine the whole pool instead of 301 self-generated candidates.
That can only pay off if the engine's SURROGATE ranks a pool better than ridge does. So
this measures exactly that, with no campaign:

    train on a random sample the size of a real budget
    predict on 3,000 random held-out variants from the whole library
    report Spearman against truth

Three arms on IDENTICAL one-hot features, so the only thing varying is the model class:

    ridge   sklearn Ridge(alpha=1.0)          - the benchmark's own baseline
    gp      Matern GP, the engine's surrogate family
    rf      RandomForest                      - a nonlinear reference

Reading:
    gp > ridge   enumerate mode should beat ridge; the surrogate is the asset
    gp ~ ridge   enumerate mode gives parity; the selection math is not a differentiator
    gp < ridge   we have a surrogate problem that enumerate mode cannot fix

Also reports the fitness spread of the 301 candidates the engine actually proposed on
eqFP611-red, to show why the in-run correlation of -0.0069 was uninformative.
"""
import csv
import random
import statistics
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import benchmark_panoptes_isolated as bench  # noqa: E402

LANDSCAPES = ["trpb", "gb1", "her2", "creilov", "eqfp611_blue", "eqfp611_red",
              "cr9114_h1", "cr6261_h1", "cr6261_h9"]
HELD_OUT = 3000
N_REPEATS = 5


def load(spec):
    wt = dict(zip(spec.positions, spec.parent_residues))
    out = {}
    for row in csv.DictReader(open(spec.dataset)):
        mutant = row["mutant"].strip()
        try:
            score = float(row["DMS_score"])
        except (TypeError, ValueError):
            continue
        if spec.key == "her2" and score == 0.0:
            continue                      # rule 2.3
        combo = dict(wt)
        if mutant.upper() != "WT":
            ok = True
            for sub in mutant.split(":"):
                sub = sub.strip()
                if not sub:
                    continue
                pos = int(sub[1:-1])
                if pos not in combo:
                    ok = False
                    break
                combo[pos] = sub[-1]
            if not ok:
                continue
        out["".join(combo[p] for p in spec.positions)] = score
    return out


def budget_for(key):
    return {"eqfp611_red": 81, "eqfp611_blue": 81,
            "cr6261_h1": 19, "cr6261_h9": 19, "her2": 113}.get(key, 125)


print(f"{'landscape':<14}{'budget':>7}{'ridge':>9}{'gp':>9}{'rf':>9}"
      f"{'  gp-ridge':>11}  verdict")
print("-" * 76)

summary = []
for key in LANDSCAPES:
    spec = bench.LIBRARIES[key]
    table = load(spec)
    combos = list(table)
    budget = budget_for(key)

    scores = {"ridge": [], "gp": [], "rf": []}
    for repeat in range(N_REPEATS):
        rng = random.Random(1000 + repeat)
        shuffled = combos[:]
        rng.shuffle(shuffled)
        train = shuffled[:budget]
        test = shuffled[budget:budget + HELD_OUT]
        if len(test) < 200:
            test = shuffled[budget:]
        Xtr = bench._one_hot_matrix(train, spec)
        ytr = np.array([table[c] for c in train])
        Xte = bench._one_hot_matrix(test, spec)
        yte = np.array([table[c] for c in test])

        models = {
            "ridge": Ridge(alpha=1.0),
            # Engine's surrogate family: Matern with a learned noise term. normalize_y
            # matters at this sample size, since the assay scales differ by 10x.
            "gp": GaussianProcessRegressor(
                kernel=ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5)
                + WhiteKernel(noise_level=0.1),
                normalize_y=True, random_state=repeat),
            "rf": RandomForestRegressor(n_estimators=200, random_state=repeat, n_jobs=1),
        }
        for name, model in models.items():
            try:
                model.fit(Xtr, ytr)
                pred = model.predict(Xte)
                rho = spearmanr(pred, yte)[0]
                if rho == rho:
                    scores[name].append(float(rho))
            except Exception as exc:
                print(f"  {key} {name} failed: {type(exc).__name__}: {exc}")

    means = {n: (statistics.fmean(v) if v else float("nan")) for n, v in scores.items()}
    delta = means["gp"] - means["ridge"]
    if delta > 0.05:
        verdict = "GP better - surrogate is the asset"
    elif delta < -0.05:
        verdict = "** GP WORSE than ridge **"
    else:
        verdict = "parity - not a differentiator"
    print(f"{key:<14}{budget:>7}{means['ridge']:>9.3f}{means['gp']:>9.3f}"
          f"{means['rf']:>9.3f}{delta:>11.3f}  {verdict}")
    summary.append((key, means["ridge"], means["gp"], means["rf"], delta))

print()
wins = sum(1 for _, _, _, _, d in summary if d > 0.05)
losses = sum(1 for _, _, _, _, d in summary if d < -0.05)
print(f"GP beats ridge on {wins} of {len(summary)} landscapes, loses on {losses}, "
      f"ties on {len(summary) - wins - losses}.")
best_rf = sum(1 for _, r, g, f, _ in summary if f > max(r, g) + 0.05)
print(f"RandomForest beats both on {best_rf} of {len(summary)} - relevant because the "
      f"benchmark already ships random_forest as a baseline.")

# Why the in-run correlation was uninformative: restricted range.
print()
R = "benchmark_comparison_results/run_20261031_225544_8d19b08e"
truth = load(bench.LIBRARIES["eqfp611_red"])
spec = bench.LIBRARIES["eqfp611_red"]
wt = dict(zip(spec.positions, spec.parent_residues))
proposed = set()
for row in csv.DictReader(open(R + "/candidate_trace.csv")):
    proposed.add(row["variant_id"])
vals = []
for notation in proposed:
    combo = dict(wt)
    ok = True
    for sub in notation.split(":"):
        sub = sub.strip()
        if not sub:
            continue
        pos = int(sub[1:-1])
        if pos not in combo:
            ok = False
            break
        combo[pos] = sub[-1]
    if ok:
        key2 = "".join(combo[p] for p in spec.positions)
        if key2 in truth:
            vals.append(truth[key2])
allv = list(truth.values())
print(f"eqFP611-red, fitness spread:")
print(f"  the {len(vals)} candidates the engine proposed: "
      f"mean {statistics.fmean(vals):.4f}  sd {statistics.stdev(vals):.4f}")
print(f"  the whole library:                           "
      f"mean {statistics.fmean(allv):.4f}  sd {statistics.stdev(allv):.4f}")
print("  A correlation measured inside that narrower set is attenuated by restricted")
print("  range, which is why the in-run rho of -0.0069 could not be read as evidence")
print("  about the surrogate.")
