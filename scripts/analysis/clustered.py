"""Seed-clustered re-analysis of the 2^5 factorial main effects and interactions.

The error being corrected: the first analysis pooled every background-pair's per-seed
difference into one sample and ran a t-test over all of them. For a component there are
16 background configurations x 48 seeds = 768 differences, but each SEED contributes 16
of them, and those 16 share a starting sample and a random stream. Treating them as 768
independent observations understates the standard error by up to sqrt(16) = 4x when the
within-seed correlation is high, so the reported t-statistics were inflated.

The correct unit of independence is the seed. For each component:
  1. per seed, average the removal contrast over all 16 backgrounds -> one value per seed
  2. a paired t over those 48 seed-level values
  3. a bootstrap CI resampling whole seeds, which preserves the within-seed correlation

The same correction applies to the two-way interactions (8 backgrounds each).

This also resolves the observation count: 768 differences are computed, of which the
earlier report quoted 734 - that was wins+losses, excluding exact ties, and was
mislabelled as the number of observations.
"""
from __future__ import annotations

import csv
import json
import statistics as st
from itertools import combinations
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
SHIPPED = {"gb1": "run_20261031_023520_b5ae03df",
           "trpb": "run_20261030_233220_8e467def"}
COMPONENTS = ("prior", "dkl", "dpp", "gflownet", "pareto")
TOP1 = "cum_top1_hits_including_initial"
T_CRIT = 2.01
RNG = np.random.default_rng(20261030)


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


def series(run: str, method: str) -> Dict[int, float]:
    return {s: float(r[TOP1]) for s, r in final_rows(run, method).items()
            if r.get(TOP1) not in (None, "")}


def seed_level_t(per_seed: Dict[int, float]) -> dict:
    """Paired t and bootstrap CI over SEEDS, the unit of independence."""
    seeds = sorted(per_seed)
    arr = np.array([per_seed[s] for s in seeds], dtype=float)
    mean, sd = float(arr.mean()), float(arr.std(ddof=1))
    se = sd / np.sqrt(len(arr))
    boots = arr[RNG.integers(0, len(arr), size=(20000, len(arr)))].mean(axis=1)
    return {"n_seeds": len(arr), "mean": mean, "sd": sd,
            "t": mean / se if se > 0 else float("nan"),
            "lo": float(np.percentile(boots, 2.5)),
            "hi": float(np.percentile(boots, 97.5)),
            "wins": int((arr > 0).sum()), "losses": int((arr < 0).sum())}


def naive_t(diffs: Sequence[float]) -> float:
    """The original, incorrect statistic: every difference treated as independent."""
    arr = np.asarray(diffs, dtype=float)
    se = arr.std(ddof=1) / np.sqrt(len(arr))
    return float(arr.mean() / se) if se > 0 else float("nan")


arms_raw = json.loads((SCRATCH / "joint_arms.json").read_text(encoding="utf-8"))
report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


emit("=" * 112)
emit("SEED-CLUSTERED RE-ANALYSIS OF THE 2^5 FACTORIAL")
emit("=" * 112)
emit("  The unit of independence is the seed, not the background-pair difference.")
emit("  'naive t' is the original statistic, kept so the inflation factor is visible.")

for lib in ("gb1", "trpb"):
    subsets: Dict[FrozenSet[str], str] = {frozenset(): SHIPPED[lib]}
    for label, run in arms_raw[lib].items():
        subsets[frozenset(t for t in label.split(",") if t)] = run
    top1 = {k: series(v, "panoptes_default") for k, v in subsets.items()}

    emit()
    emit("=" * 112)
    emit(f"{lib.upper()}  —  main effects of REMOVING each component")
    emit("=" * 112)
    emit(f"  {'component':<10}{'effect':>9}{'clustered t':>13}{'95% CI (seeds)':>22}"
         f"{'seeds +/-':>11}{'naive t':>10}{'inflation':>11}{'n diffs':>9}")
    for c in COMPONENTS:
        backgrounds = []
        for mask in range(32):
            ss = frozenset(x for b, x in enumerate(COMPONENTS) if mask & (1 << b))
            if c in ss:
                continue
            if (ss | {c}) in top1 and ss in top1:
                backgrounds.append((top1[ss | {c}], top1[ss]))
        # per-seed average contrast across all backgrounds
        seeds = sorted(set.intersection(*[set(a) & set(b) for a, b in backgrounds]))
        per_seed = {s: st.fmean([a[s] - b[s] for a, b in backgrounds]) for s in seeds}
        flat = [a[s] - b[s] for a, b in backgrounds for s in seeds]
        r = seed_level_t(per_seed)
        nt = naive_t(flat)
        emit(f"  {c:<10}{r['mean']:>+9.2f}{r['t']:>13.2f}"
             f"{('[' + format(r['lo'], '+.2f') + ', ' + format(r['hi'], '+.2f') + ']'):>22}"
             f"{(str(r['wins']) + '/' + str(r['losses'])):>11}{nt:>10.2f}"
             f"{(nt / r['t'] if r['t'] else float('nan')):>10.2f}x{len(flat):>9}")

    emit()
    emit(f"  {lib.upper()}  —  two-way interactions, seed-clustered")
    emit(f"  {'pair':<22}{'effect':>9}{'clustered t':>13}{'95% CI (seeds)':>22}{'naive t':>10}")
    rows = []
    for c, d in combinations(COMPONENTS, 2):
        quads = []
        for mask in range(32):
            ss = frozenset(x for b, x in enumerate(COMPONENTS) if mask & (1 << b))
            if c in ss or d in ss:
                continue
            quad = (ss | {c, d}, ss | {c}, ss | {d}, ss)
            if all(q in top1 for q in quad):
                quads.append(quad)
        seeds = sorted(set.intersection(
            *[set(top1[q]) for quad in quads for q in quad]))
        per_seed = {}
        flat = []
        for s in seeds:
            vals = [(top1[b][s] - top1[jc][s]) - (top1[jd][s] - top1[n][s])
                    for b, jc, jd, n in quads]
            per_seed[s] = st.fmean(vals)
            flat.extend(vals)
        r = seed_level_t(per_seed)
        rows.append((abs(r["t"]), c, d, r, naive_t(flat)))
    for _, c, d, r, nt in sorted(rows, reverse=True):
        flag = "  significant" if abs(r["t"]) >= T_CRIT else "  ns"
        emit(f"  {c + ' x ' + d:<22}{r['mean']:>+9.2f}{r['t']:>13.2f}"
             f"{('[' + format(r['lo'], '+.2f') + ', ' + format(r['hi'], '+.2f') + ']'):>22}"
             f"{nt:>10.2f}{flag}")

emit()
emit("=" * 112)
emit("WHAT SURVIVES THE CORRECTION")
emit("=" * 112)
emit("  The held-out confirmation of 'remove dkl,gflownet' is UNAFFECTED: it is a single")
emit("  pair of arms compared over 48 independent seeds, correctly paired, one difference")
emit("  per seed. GB1 +11.83 (t=5.75) and TrpB +5.23 (t=5.57) stand as reported.")
emit("  The same is true of every arm-vs-baseline comparison and of the McNemar tail tests.")
emit("  Only the factorial main effects and interactions were affected.")

(SCRATCH / "clustered_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'clustered_report.txt'}")
