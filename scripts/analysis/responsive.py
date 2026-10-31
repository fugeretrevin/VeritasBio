"""Is the degenerate TrpB frontier configuration actually ignoring its data?

Workstream 2.4. All 96 campaigns of `frontier - dkl,gflownet` on TrpB ended at exactly
0.7036. The review is right that zero terminal variance alone does not prove a policy
has stopped learning - a policy that deterministically converges to the same answer
would look identical at the endpoint.

Two explanations had to be separated, and the first is already ruled out:

  (a) every campaign found the global optimum. The TrpB manifest records global_best
      = 1.0, and 0.7036 is well short of it, so this is not what happened.

  (b) the policy's selections do not depend on the measurements. This is what the test
      below examines.

The discriminator does not need new compute. Each seed starts from a DIFFERENT random
25-variant initial sample. A policy that conditions on its data must therefore choose
different candidates in round 1. If the selected sets are identical across seeds whose
training data differ, the selections cannot be a function of that data.

Reported per arm:
  jaccard    mean pairwise Jaccard of the selected (non-initial) sets across seeds.
             ~1.0 means every campaign measured the same variants.
  universal  share of a campaign's selections that appear in EVERY seed's selections.
  union      size of the union of selections over all seeds, against the per-campaign
             count. A responsive policy explores a far larger union than any one run.
  r1_jac     the same Jaccard restricted to round 1, the first decision taken after
             seeing data that differs by construction.

Healthy and baseline arms are included as controls, so "low" and "high" are calibrated
against policies on the same landscape rather than against intuition.
"""
from __future__ import annotations

import csv
import json
import statistics as st
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
DATASET = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data"
               r"\TRPB_THEMA_Johnston_2024.csv")

# (label, [(run, method), ...]) - the suspect arm spans both seed blocks (96 campaigns)
ARMS: List[Tuple[str, List[Tuple[str, str]]]] = [
    ("frontier - dkl,gflownet  [SUSPECT]", [
        ("run_20261031_025657_f6f5b5fd", "panoptes_qualify_then_optimize"),
        ("run_20261031_025926_eada00c6", "panoptes_qualify_then_optimize")]),
    ("frontier - dkl,dpp,gflownet", [
        ("run_20261031_025739_fe0dc20e", "panoptes_qualify_then_optimize"),
        ("run_20261031_030003_74f8caa3", "panoptes_qualify_then_optimize")]),
    ("frontier  (shipped)", [
        ("run_20261031_025608_f51a714f", "panoptes_qualify_then_optimize"),
        ("run_20261031_025827_8dfabb93", "panoptes_qualify_then_optimize")]),
    ("reliability - dkl,gflownet", [
        ("run_20261031_020245_489d3e08", "panoptes_default"),
        ("run_20261031_024735_a7e2cd7c", "panoptes_default")]),
    ("reliability  (shipped)", [
        ("run_20261031_024533_c41a0125", "panoptes_default")]),
    ("ridge      [control]", [("run_20261031_025252_7dce6099", "ridge")]),
    ("random     [floor]", [("run_20261031_025252_7dce6099", "random")]),
]


def load(run: str, method: str
         ) -> Tuple[Dict[Tuple[str, int], Set[str]],
                    Dict[Tuple[str, int], Set[str]],
                    Dict[Tuple[str, int], Set[str]],
                    Dict[Tuple[str, int], float]]:
    """-> (selected, round1, initial, best_fitness) keyed by (run, seed)."""
    sel: Dict[Tuple[str, int], Set[str]] = defaultdict(set)
    r1: Dict[Tuple[str, int], Set[str]] = defaultdict(set)
    init: Dict[Tuple[str, int], Set[str]] = defaultdict(set)
    best: Dict[Tuple[str, int], float] = {}
    with open(ROOT / run / "selection_trace.csv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            try:
                seed = int(float(r["seed"]))
                rnd = int(float(r["round"]))
            except (TypeError, ValueError):
                continue
            vid = str(r.get("variant_id", "")).strip()
            if not vid:
                continue
            key = (run, seed)
            status = (r.get("status") or "").strip()
            if status in ("rejected", "considered", ""):
                continue
            try:
                f = float(r["revealed_fitness"])
                best[key] = max(best.get(key, float("-inf")), f)
            except (TypeError, ValueError):
                pass
            if status == "initial":
                init[key].add(vid)
            else:
                sel[key].add(vid)
                if rnd == 1:
                    r1[key].add(vid)
    return sel, r1, init, best


def mean_jaccard(sets: Sequence[Set[str]], cap: int = 40) -> float:
    sets = [s for s in sets if s][:cap]
    if len(sets) < 2:
        return float("nan")
    vals = [len(a & b) / len(a | b) for a, b in combinations(sets, 2)]
    return st.fmean(vals)


report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


emit("=" * 112)
emit("DATA-RESPONSIVENESS TEST - TrpB")
emit("=" * 112)
emit("  Every seed starts from a different random 25-variant sample. A policy that")
emit("  conditions on its measurements must therefore pick different candidates.")
emit("  TrpB global_best = 1.0, so terminating at 0.7036 is NOT convergence to the optimum.")
emit()
emit(f"  {'arm':<36}{'camps':>7}{'bestF':>9}{'sd':>8}{'sel':>6}{'union':>8}"
     f"{'jaccard':>9}{'universal':>11}{'r1_jac':>9}")

results = {}
for label, sources in ARMS:
    sel: Dict[Tuple[str, int], Set[str]] = {}
    r1: Dict[Tuple[str, int], Set[str]] = {}
    init: Dict[Tuple[str, int], Set[str]] = {}
    best: Dict[Tuple[str, int], float] = {}
    for run, method in sources:
        s, a, i, b = load(run, method)
        sel.update(s); r1.update(a); init.update(i); best.update(b)
    if not sel:
        emit(f"  {label:<36}{'(no trace rows)':>60}")
        continue
    keys = sorted(sel)
    sets = [sel[k] for k in keys]
    bf = [best[k] for k in keys if k in best]
    union: Set[str] = set().union(*sets)
    inter: Set[str] = set.intersection(*sets) if len(sets) > 1 else set(sets[0])
    percamp = st.fmean(len(s) for s in sets)
    jac = mean_jaccard(sets)
    rj = mean_jaccard([r1[k] for k in keys if k in r1])
    results[label] = dict(jac=jac, r1=rj, union=len(union), percamp=percamp,
                          universal=len(inter) / percamp if percamp else float("nan"),
                          n=len(keys), bf=bf, init=[init[k] for k in keys if k in init])
    emit(f"  {label:<36}{len(keys):>7}"
         f"{(st.fmean(bf) if bf else float('nan')):>9.4f}"
         f"{(st.stdev(bf) if len(bf) > 1 else 0.0):>8.4f}"
         f"{percamp:>6.0f}{len(union):>8}{jac:>9.3f}"
         f"{(len(inter) / percamp * 100 if percamp else float('nan')):>10.1f}%{rj:>9.3f}")

emit()
emit("  sel       = variants selected per campaign (excludes the 25 initial)")
emit("  union     = distinct variants any campaign ever selected, across all seeds")
emit("  jaccard   = mean pairwise overlap of selected sets (1.000 = all identical)")
emit("  universal = share of a campaign's picks that EVERY campaign also made")
emit("  r1_jac    = overlap of the round-1 batch, the first decision after seeing data")

# Sanity: the initial samples really do differ across seeds.
emit()
emit("-" * 112)
emit("CONTROL: do the initial samples actually differ?")
emit("-" * 112)
for label in ("frontier - dkl,gflownet  [SUSPECT]", "reliability - dkl,gflownet"):
    r = results.get(label)
    if not r or len(r["init"]) < 2:
        continue
    emit(f"  {label:<36} initial-sample jaccard = {mean_jaccard(r['init']):.4f}"
         f"   (near 0 means each seed starts from genuinely different data)")

# What is the 0.7036 variant, and is it always picked?
emit()
emit("-" * 112)
emit("WHAT IS 0.7036?")
emit("-" * 112)
target: List[Tuple[str, float]] = []
with open(DATASET, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        try:
            f = float(r["DMS_score"])
        except (TypeError, ValueError):
            continue
        if 0.7035 <= f <= 0.7037:
            target.append((str(r["mutant"]), f))
emit(f"  {len(target)} variant(s) in the landscape with fitness in [0.7035, 0.7037]:")
for vid, f in target[:10]:
    emit(f"    {vid}   {f:.6f}")

sus_sel: Dict[Tuple[str, int], Set[str]] = {}
for run, method in ARMS[0][1]:
    s, _a, _i, _b = load(run, method)
    sus_sel.update(s)
if target and sus_sel:
    emit()
    for vid, f in target[:5]:
        hits = sum(1 for s in sus_sel.values() if vid in s)
        emit(f"  '{vid}' selected in {hits} of {len(sus_sel)} suspect campaigns")

# How much of the landscape does the arm ever touch?
emit()
emit("-" * 112)
emit("VERDICT INPUTS")
emit("-" * 112)
sus = results.get("frontier - dkl,gflownet  [SUSPECT]")
hea = results.get("reliability - dkl,gflownet")
if sus and hea:
    emit(f"  suspect : {sus['n']} campaigns, {sus['percamp']:.0f} picks each, "
         f"union {sus['union']} distinct  -> {sus['union'] / sus['percamp']:.2f}x one campaign")
    emit(f"  healthy : {hea['n']} campaigns, {hea['percamp']:.0f} picks each, "
         f"union {hea['union']} distinct  -> {hea['union'] / hea['percamp']:.2f}x one campaign")
    emit(f"  A policy whose union barely exceeds a single campaign is re-deriving the same")
    emit(f"  batch regardless of what it measured.")

(SCRATCH / "responsive_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'responsive_report.txt'}")
