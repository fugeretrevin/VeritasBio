"""Analyse the 2^5 joint component ablation on GB1 and TrpB.

The one-at-a-time sweep could not attribute the engine's behaviour. On TrpB the five
single removals cost more between them than the engine's whole advantage over random,
which means the components are redundant, not additive. On GB1 three of them each
IMPROVED the arm, yet the best recovered only a quarter of the deficit. Both of those
are statements about interactions, which only a full factorial can measure.

What this computes, per landscape:
  - every subset's absolute level, against ridge / random forest / ALDE measured on the
    same seeds in the same run;
  - factorial main effects: for each component, the average effect of removing it taken
    over all 16 subset pairs that differ only in that component, pooled over seeds and
    paired seed-for-seed throughout;
  - two-way interactions, to say whether the components substitute for each other;
  - the extreme tail, P(top 0.01%), because that is the only metric on which the arm's
    deficit replicated across both landscapes.

Every arm ran 48 seeds from base seed 0 with an identical shared starting sample per
seed, so all comparisons are paired.
"""
from __future__ import annotations

import csv
import json
import statistics as st
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio")
import benchmark_panoptes_isolated as iso

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")

# The two reference runs. GB1's was re-measured in the current code state and reproduced
# the earlier baseline's 14.48 exactly, which is what confirms the pool-mode change is
# inert when no pool is supplied; TrpB's is the 48-seed run that also carries the
# baselines, so every comparison below is within-run and within-seed.
SHIPPED = {"gb1": "run_20261031_023520_b5ae03df",
           "trpb": "run_20261030_233220_8e467def"}
BASELINE_RUN = {"gb1": "run_20261030_235620_fd18d20e",
                "trpb": "run_20261030_233220_8e467def"}
COMPONENTS = ("prior", "dkl", "dpp", "gflownet", "pareto")
T_CRIT = 2.01
TOP1 = "cum_top1_hits_including_initial"
BEST = "best_measured_fitness"
QUANTILES = (99.0, 99.9, 99.99)


def final_rows(run: str, method: str) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    with open(ROOT / run / "round_metrics.csv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            seed = int(float(r["seed"]))
            prev = out.get(seed)
            if prev is None or int(float(r["round"])) >= int(float(prev["round"])):
                out[seed] = r
    return out


def series(run: str, method: str, key: str) -> Dict[int, float]:
    return {s: float(r[key]) for s, r in final_rows(run, method).items()
            if r.get(key) not in (None, "")}


def paired(a: Dict[int, float], b: Dict[int, float]) -> Optional[dict]:
    seeds = sorted(set(a) & set(b))
    d = [a[s] - b[s] for s in seeds]
    if len(d) < 3:
        return None
    mean, sd = st.fmean(d), st.stdev(d)
    se = sd / np.sqrt(len(d))
    return {"n": len(d), "mean": mean, "t": (mean / se) if se > 0 else float("nan"),
            "wins": sum(1 for x in d if x > 0), "losses": sum(1 for x in d if x < 0)}


def pooled_paired(pairs: Sequence[Tuple[Dict[int, float], Dict[int, float]]]) -> dict:
    """One t over the concatenated per-seed differences of many subset pairs.

    Each pair contributes 48 paired differences. Treating the pooled differences as one
    sample is the factorial main effect: it averages the component's effect over every
    background configuration instead of only over the shipped one.
    """
    d: List[float] = []
    for a, b in pairs:
        for s in sorted(set(a) & set(b)):
            d.append(a[s] - b[s])
    mean, sd = st.fmean(d), st.stdev(d)
    se = sd / np.sqrt(len(d))
    return {"n": len(d), "mean": mean, "t": (mean / se) if se > 0 else float("nan"),
            "wins": sum(1 for x in d if x > 0), "losses": sum(1 for x in d if x < 0)}


arms_raw = json.loads((SCRATCH / "joint_arms.json").read_text(encoding="utf-8"))

report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


for lib in ("gb1", "trpb"):
    # subset (as a frozenset of removed components) -> run id
    subsets: Dict[FrozenSet[str], str] = {frozenset(): SHIPPED[lib]}
    for label, run in arms_raw[lib].items():
        subsets[frozenset(t for t in label.split(",") if t)] = run
    missing = [s for s in range(32) if frozenset(
        c for b, c in enumerate(COMPONENTS) if s & (1 << b)) not in subsets]
    emit()
    emit("=" * 108)
    emit(f"{lib.upper()}  —  2^5 joint ablation, {len(subsets)}/32 subsets, 48 paired seeds each")
    emit("=" * 108)
    if missing:
        emit(f"  MISSING {len(missing)} subsets")

    # Landscape quantiles, for the extreme-tail columns.
    iso.configure_library(lib)
    universe = iso.load_universe(
        dataset_path=iso.DEFAULT_DATASET, provenance_mode="acknowledge-unflagged",
        imputed_ids_path=None, max_variants=None, subsample_seed=0)
    values = np.array(list(universe.fitness.values()), dtype=float)
    thresh = {q: float(np.percentile(values, q)) for q in QUANTILES}

    top1 = {k: series(v, "panoptes_default", TOP1) for k, v in subsets.items()}
    best = {k: series(v, "panoptes_default", BEST) for k, v in subsets.items()}
    shipped_top1 = top1[frozenset()]

    base_levels = {}
    for m in ("random", "ridge", "random_forest", "alde"):
        vals = series(BASELINE_RUN[lib], m, TOP1)
        if vals:
            base_levels[m] = st.fmean(vals.values())
    emit("  baselines on the same seeds: "
         + "   ".join(f"{m}={v:.2f}" for m, v in base_levels.items()))
    emit()

    # ---------------------------------------------------------------- ranking
    emit(f"  {'removed components':<40}{'top-1%':>8}{'vs ship':>9}{'t':>7}"
         f"{'best fit':>10}{'P1%':>6}{'P.1%':>7}{'P.01%':>7}")
    rows = []
    for ss, run in subsets.items():
        lvl = st.fmean(top1[ss].values())
        cmp_ = paired(top1[ss], shipped_top1) if ss else None
        bvals = list(best[ss].values())
        tail = tuple(100.0 * sum(1 for v in bvals if v >= thresh[q]) / len(bvals)
                     for q in QUANTILES)
        rows.append((lvl, ss, cmp_, st.fmean(bvals), tail))
    for lvl, ss, cmp_, bf, tail in sorted(rows, key=lambda r: -r[0]):
        label = ",".join(sorted(ss)) if ss else "<none: shipped>"
        mark = "" if cmp_ is None else ("*" if cmp_["t"] >= T_CRIT
                                        else ("!" if cmp_["t"] <= -T_CRIT else " "))
        emit(f"  {label:<40}{lvl:>8.2f}"
             + (f"{cmp_['mean']:>+9.2f}{cmp_['t']:>7.2f}" if cmp_ else f"{'-':>9}{'-':>7}")
             + f"{mark:<1}{bf:>9.4f}{tail[0]:>5.0f}%{tail[1]:>6.0f}%{tail[2]:>6.0f}%")
    emit("    * significantly better than shipped    ! significantly worse")

    # ------------------------------------------------------- main effects
    emit()
    emit(f"  FACTORIAL MAIN EFFECTS — average effect of REMOVING each component,")
    emit(f"  over all 16 background configurations (positive = removing it helps):")
    emit(f"    {'component':<12}{'top-1% delta':>14}{'t':>8}{'won/lost':>12}"
         f"{'one-at-a-time':>16}")
    for c in COMPONENTS:
        pairs = []
        for mask in range(32):
            ss = frozenset(x for b, x in enumerate(COMPONENTS) if mask & (1 << b))
            if c in ss:
                continue
            with_c = ss | {c}
            if with_c in top1 and ss in top1:
                pairs.append((top1[with_c], top1[ss]))
        eff = pooled_paired(pairs)
        solo = paired(top1[frozenset({c})], shipped_top1)
        solo_txt = f"{solo['mean']:>+16.2f}" if solo else f"{'-':>16}"
        emit(f"    {c:<12}{eff['mean']:>+14.2f}{eff['t']:>8.2f}"
             f"{str(eff['wins']) + '/' + str(eff['losses']):>12}" + solo_txt)

    # ------------------------------------------------------- interactions
    emit()
    emit("  TWO-WAY INTERACTIONS (positive = the two substitute for each other, so")
    emit("  removing both costs less than the sum of removing each alone):")
    inter = []
    for c, d in combinations(COMPONENTS, 2):
        pairs = []
        for mask in range(32):
            ss = frozenset(x for b, x in enumerate(COMPONENTS) if mask & (1 << b))
            if c in ss or d in ss:
                continue
            quad = [ss | {c, d}, ss | {c}, ss | {d}, ss]
            if all(q in top1 for q in quad):
                seeds = sorted(set.intersection(*(set(top1[q]) for q in quad)))
                both, jc, jd, none_ = quad
                pairs.append({s: (top1[both][s] - top1[jc][s])
                              - (top1[jd][s] - top1[none_][s]) for s in seeds})
        d_all = [v for p in pairs for v in p.values()]
        mean, sd = st.fmean(d_all), st.stdev(d_all)
        t = mean / (sd / np.sqrt(len(d_all)))
        inter.append((abs(t), c, d, mean, t, len(d_all)))
    for _, c, d, mean, t, n in sorted(inter, reverse=True):
        flag = "  <-- significant" if abs(t) >= T_CRIT else ""
        emit(f"    {c + ' x ' + d:<24}{mean:>+9.2f}{t:>8.2f}  (n={n}){flag}")

    # ------------------------------------------------------- the prediction
    all_off = frozenset(COMPONENTS)
    lvl = st.fmean(top1[all_off].values())
    emit()
    emit(f"  ALL FIVE REMOVED: top-1% {lvl:.2f}   "
         + "   ".join(f"vs {m} {lvl - v:+.2f}" for m, v in base_levels.items()))
    for m, v in base_levels.items():
        cmp_ = paired(top1[all_off], series(BASELINE_RUN[lib], m, TOP1))
        if cmp_:
            verdict = ("BETTER" if cmp_["mean"] > 0 else "WORSE") \
                if abs(cmp_["t"]) >= T_CRIT else "not significant"
            emit(f"    vs {m:<16}{cmp_['mean']:>+8.2f}  t={cmp_['t']:>+6.2f}  {verdict}")

(SCRATCH / "joint_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'joint_report.txt'}")
