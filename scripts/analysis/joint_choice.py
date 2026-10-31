"""Pick a configuration from the 2^5 factorial, and test the prior/tail tension.

Two questions the per-landscape rankings do not answer.

1. Is there ONE configuration that is significantly better than the shipped engine on
   both landscapes? Per-landscape optima disagree - GB1's best removes the prior, which
   is catastrophic on TrpB - so the shippable answer is the best configuration that wins
   on both, not the best on either.

2. The extreme tail. Across both landscapes the shipped arm reaches the top 0.01% far
   less often than ALDE, and that is the one deficit that replicated. In the TrpB table
   the subsets with the best tail are the ones that REMOVE the prior, while those same
   subsets lose the most top-1% yield. If that holds up it is a mechanism, not a
   coincidence: the prior concentrates the search, which lifts mean yield and suppresses
   the variance needed to reach the extreme tail. Tested here directly.
"""
from __future__ import annotations

import csv
import json
import statistics as st
import sys
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Tuple

import numpy as np

sys.path.insert(0, r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio")
import benchmark_panoptes_isolated as iso

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
SHIPPED = {"gb1": "run_20261031_023520_b5ae03df",
           "trpb": "run_20261030_233220_8e467def"}
BASELINE_RUN = {"gb1": "run_20261030_235620_fd18d20e",
                "trpb": "run_20261030_233220_8e467def"}
COMPONENTS = ("prior", "dkl", "dpp", "gflownet", "pareto")
T_CRIT = 2.01
TOP1 = "cum_top1_hits_including_initial"
BEST = "best_measured_fitness"


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
    return {"n": len(d), "mean": mean, "sd": sd,
            "t": (mean / se) if se > 0 else float("nan"),
            "wins": sum(1 for x in d if x > 0), "losses": sum(1 for x in d if x < 0)}


arms_raw = json.loads((SCRATCH / "joint_arms.json").read_text(encoding="utf-8"))
report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


land: Dict[str, dict] = {}
for lib in ("gb1", "trpb"):
    subsets: Dict[FrozenSet[str], str] = {frozenset(): SHIPPED[lib]}
    for label, run in arms_raw[lib].items():
        subsets[frozenset(t for t in label.split(",") if t)] = run
    iso.configure_library(lib)
    universe = iso.load_universe(
        dataset_path=iso.DEFAULT_DATASET, provenance_mode="acknowledge-unflagged",
        imputed_ids_path=None, max_variants=None, subsample_seed=0)
    values = np.array(list(universe.fitness.values()), dtype=float)
    land[lib] = {
        "subsets": subsets,
        "top1": {k: series(v, "panoptes_default", TOP1) for k, v in subsets.items()},
        "best": {k: series(v, "panoptes_default", BEST) for k, v in subsets.items()},
        "t9999": float(np.percentile(values, 99.99)),
        "t999": float(np.percentile(values, 99.9)),
        "base": {m: series(BASELINE_RUN[lib], m, TOP1)
                 for m in ("random", "ridge", "random_forest", "alde")},
    }

ALL = [frozenset(c for b, c in enumerate(COMPONENTS) if mask & (1 << b))
       for mask in range(32)]

# ------------------------------------------------- 1. one configuration for both
emit("=" * 104)
emit("A SINGLE CONFIGURATION FOR BOTH LANDSCAPES")
emit("=" * 104)
emit("  Significance is against the shipped engine on the same 48 seeds, per landscape.")
emit("  'worst t' is the weaker of the two, so a configuration only qualifies if it wins")
emit("  on BOTH. Ranked by that worst case.")
emit()
emit(f"  {'removed':<34}{'GB1':>8}{'d':>8}{'t':>7}{'TrpB':>9}{'d':>8}{'t':>7}{'worst t':>10}")
cands = []
for ss in ALL:
    if not ss:
        continue
    g = paired(land["gb1"]["top1"][ss], land["gb1"]["top1"][frozenset()])
    t = paired(land["trpb"]["top1"][ss], land["trpb"]["top1"][frozenset()])
    if not (g and t):
        continue
    cands.append((min(g["t"], t["t"]), ss, g, t))
for worst, ss, g, t in sorted(cands, reverse=True)[:8]:
    gl = st.fmean(land["gb1"]["top1"][ss].values())
    tl = st.fmean(land["trpb"]["top1"][ss].values())
    emit(f"  {','.join(sorted(ss)):<34}{gl:>8.2f}{g['mean']:>+8.2f}{g['t']:>7.2f}"
         f"{tl:>9.2f}{t['mean']:>+8.2f}{t['t']:>7.2f}{worst:>10.2f}")

winner = max(cands)[1]
emit()
emit(f"  CHOICE: remove {','.join(sorted(winner))}  "
     f"(keep {','.join(c for c in COMPONENTS if c not in winner)})")
for lib in ("gb1", "trpb"):
    lvl = st.fmean(land[lib]["top1"][winner].values())
    ship = st.fmean(land[lib]["top1"][frozenset()].values())
    emit(f"    {lib.upper():<6} {ship:.2f} -> {lvl:.2f}")
    for m, base in land[lib]["base"].items():
        c_new = paired(land[lib]["top1"][winner], base)
        c_old = paired(land[lib]["top1"][frozenset()], base)
        def verdict(c):
            return ("BETTER" if c["mean"] > 0 else "WORSE") \
                if abs(c["t"]) >= T_CRIT else "ns"
        emit(f"      vs {m:<14} shipped {c_old['mean']:>+7.2f} ({verdict(c_old):<6})"
             f"   ->  {c_new['mean']:>+7.2f} (t={c_new['t']:>+6.2f}, {verdict(c_new)})")

# ---------------------------------------------- 2. the prior / tail tension
emit()
emit("=" * 104)
emit("DOES THE PRIOR BUY TOP-1% YIELD BY GIVING UP THE EXTREME TAIL?")
emit("=" * 104)
emit("  Across the 16 subsets that keep the prior vs the 16 that remove it, holding the")
emit("  other four components balanced by construction.")
emit()
for lib in ("gb1", "trpb"):
    L = land[lib]
    keep, drop = [], []
    for ss in ALL:
        bvals = list(L["best"][ss].values())
        tail = 100.0 * sum(1 for v in bvals if v >= L["t9999"]) / len(bvals)
        lvl = st.fmean(L["top1"][ss].values())
        (drop if "prior" in ss else keep).append((lvl, tail))
    alde_tail_vals = list(series(BASELINE_RUN[lib], "alde", BEST).values())
    alde_tail = 100.0 * sum(1 for v in alde_tail_vals if v >= L["t9999"]) / len(alde_tail_vals)
    emit(f"  {lib.upper()}   (ALDE reaches the top 0.01% in {alde_tail:.0f}% of campaigns)")
    emit(f"    prior KEPT   ({len(keep):>2} subsets): top-1% {st.fmean(k[0] for k in keep):>6.2f}"
         f"   P(top 0.01%) {st.fmean(k[1] for k in keep):>5.1f}%")
    emit(f"    prior REMOVED({len(drop):>2} subsets): top-1% {st.fmean(d[0] for d in drop):>6.2f}"
         f"   P(top 0.01%) {st.fmean(d[1] for d in drop):>5.1f}%")
    xs = np.array([s[0] for s in keep + drop])
    ys = np.array([s[1] for s in keep + drop])
    r = float(np.corrcoef(xs, ys)[0, 1])
    emit(f"    across all 32 subsets, corr(top-1% yield, P(top 0.01%)) = {r:+.3f}")
    best_tail = max(keep + drop, key=lambda s: s[1])
    emit(f"    best tail any configuration reaches: {best_tail[1]:.0f}% "
         f"(at top-1% {best_tail[0]:.2f}) vs ALDE's {alde_tail:.0f}%")
    emit()

(SCRATCH / "joint_choice.txt").write_text("\n".join(report), encoding="utf-8")
print(f"wrote {SCRATCH / 'joint_choice.txt'}")
