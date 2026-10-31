"""What does the degenerate configuration's reachable set actually contain?

The responsiveness test showed `frontier - dkl,gflownet` is not data-blind - its
round-1 batches do vary with the starting sample - but across 96 campaigns it only
ever selected 294 distinct variants, and it picked V183S:F184L (0.7036) in all 96.

That points at the candidate GENERATOR rather than the acquisition score. If removing
DKL and GFlowNet leaves a generator that proposes only a small neighbourhood of the
wild type, then no acquisition function could do better, because the high-fitness
variants are never on the ballot.

Two things decide it:
  1. the best fitness available ANYWHERE in the reachable set. If that maximum is
     0.7036, the campaigns are not failing to choose well - they are choosing the best
     of what they were offered, and the ceiling is the generator's.
  2. the mutation-depth profile of the reachable set against a healthy arm's. Depth is
     how many of the 4 sites differ from wild type (VFVS).
"""
from __future__ import annotations

import csv
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
DATASET = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data"
               r"\TRPB_THEMA_Johnston_2024.csv")
POSITIONS = (183, 184, 227, 228)
WT = {183: "V", 184: "F", 227: "V", 228: "S"}

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
    ("ridge      [control]", [("run_20261031_025252_7dce6099", "ridge")]),
]


def residues(variant_id: str) -> Optional[Tuple[str, ...]]:
    s = variant_id.strip()
    if not s:
        return None
    if s.upper() == "WT":
        return tuple(WT[p] for p in POSITIONS)
    seen: Dict[int, str] = {}
    for t in s.split(":"):
        if not t:
            continue
        digits = "".join(c for c in t if c.isdigit())
        if not digits or len(t) < 3:
            return None
        pos = int(digits)
        if pos not in WT:
            return None
        seen[pos] = t[-1]
    return tuple(seen.get(p, WT[p]) for p in POSITIONS)


def depth(g: Tuple[str, ...]) -> int:
    return sum(1 for p, a in zip(POSITIONS, g) if a != WT[p])


# landscape fitness by canonical genotype
fitness: Dict[Tuple[str, ...], float] = {}
with open(DATASET, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        g = residues(str(r["mutant"]))
        if g is None:
            continue
        try:
            fitness[g] = float(r["DMS_score"])
        except (TypeError, ValueError):
            pass

by_depth_all: Counter = Counter(depth(g) for g in fitness)
best_by_depth: Dict[int, float] = defaultdict(lambda: float("-inf"))
for g, f in fitness.items():
    best_by_depth[depth(g)] = max(best_by_depth[depth(g)], f)

report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


emit("=" * 110)
emit("WHAT IS IN THE REACHABLE SET?")
emit("=" * 110)
emit(f"  Landscape: {len(fitness):,} genotypes over positions {POSITIONS}, WT "
     f"{''.join(WT[p] for p in POSITIONS)}.")
emit(f"  Depth = how many of the 4 sites differ from WT.")
emit()
emit(f"  {'depth':<8}{'variants':>11}{'best fitness':>15}")
for d in sorted(by_depth_all):
    emit(f"  {d:<8}{by_depth_all[d]:>11,}{best_by_depth[d]:>15.4f}")
emit()

emit("=" * 110)
emit(f"  {'arm':<36}{'reach':>8}{'%land':>8}{'bestAvail':>11}{'meanDep':>9}"
     + "".join(f"{('d' + str(d)):>7}" for d in (0, 1, 2, 3, 4)))
rows = {}
for label, sources in ARMS:
    union: Set[Tuple[str, ...]] = set()
    for run, method in sources:
        with open(ROOT / run / "selection_trace.csv", newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("method") != method:
                    continue
                status = (r.get("status") or "").strip()
                if status in ("", "initial", "rejected", "considered"):
                    continue
                g = residues(str(r.get("variant_id", "")))
                if g is not None:
                    union.add(g)
    if not union:
        continue
    cnt = Counter(depth(g) for g in union)
    avail = [fitness[g] for g in union if g in fitness]
    rows[label] = dict(union=union, best=max(avail) if avail else float("nan"))
    emit(f"  {label:<36}{len(union):>8}{len(union) / len(fitness) * 100:>7.2f}%"
         f"{(max(avail) if avail else float('nan')):>11.4f}"
         f"{st.fmean(depth(g) for g in union):>9.2f}"
         + "".join(f"{cnt.get(d, 0):>7}" for d in (0, 1, 2, 3, 4)))
emit()
emit("  reach     = distinct variants the arm EVER selected across all 96 campaigns")
emit("  bestAvail = the single best fitness present anywhere in that reachable set")
emit("  d0..d4    = how many reachable variants sit at each mutation depth")

emit()
emit("=" * 110)
emit("VERDICT")
emit("=" * 110)
sus = rows.get("frontier - dkl,gflownet  [SUSPECT]")
if sus:
    emit(f"  Best fitness anywhere in the suspect arm's reachable set: {sus['best']:.4f}")
    emit(f"  Observed terminal fitness in all 96 campaigns:            0.7036")
    if abs(sus["best"] - 0.7036) < 1e-3:
        emit()
        emit("  These match. The campaigns are NOT failing to choose well - they are")
        emit("  choosing the best variant that was ever put in front of them. The failure")
        emit("  is in candidate GENERATION, not in acquisition, prediction or batching.")
        emit(f"  Global optimum 1.0000 lies at depth "
             f"{min((depth(g) for g, f in fitness.items() if f >= 0.9999), default=-1)}, "
             f"which the generator reaches "
             f"{Counter(depth(g) for g in sus['union']).get(min((depth(g) for g, f in fitness.items() if f >= 0.9999), default=-1), 0)} "
             f"times in 96 campaigns.")

(SCRATCH / "reachable_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'reachable_report.txt'}")
