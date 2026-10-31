"""Why does reliability win on TrpB and lose on GB1?

What the existing traces CAN and CANNOT answer, stated up front:

  CANNOT. The trace records the 25 accepted picks per round plus a few rejections
  (all 'not_in_eligible_universe'). The ranked candidate pool is not recorded. So
  "did high-fitness variants enter the pool, and how were they ranked" is NOT
  answerable from existing runs - it needs the enriched trace being built separately.

  CAN. Whether the failure is generation-limited, and when in a campaign it happens.
  Both are decided by what was selected, which is recorded.

Three hypotheses, and what separates them:

  H1 generation-limited. Reliability on GB1 cannot propose the right variants, as in
     the TrpB frontier collapse. Signature: small reachable set, truncated depth,
     best-available fitness below the observed best.

  H2 landscape mismatch. GB1's top-1% hits sit at a different mutation depth than
     TrpB's, and a policy that pushes deep overshoots them. Signature: the two
     landscapes' hit-depth profiles differ, and reliability's selection depth sits off
     GB1's hit mass while sitting on TrpB's.

  H3 unguided exploration. Reliability leans on the structural prior, which the
     factorial shows is decisive on TrpB (-16.85) and inert on GB1 (+0.04). Without a
     useful prior its exploration has nothing to steer by. Signature: per-round
     selection quality starts near baseline and decays, rather than starting badly.

The deficit is known to COMPOUND with budget (-2.48 at 50 measurements, -21.92 at 125),
which already argues against a bad start and for something progressive.
"""
from __future__ import annotations

import csv
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
DATA = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data")

LIB = {
    "gb1": dict(csv=DATA / "GB1_STRSG_Wu_2016.csv", pos=(39, 40, 41, 54), wt="VDGV",
                thresh=2.15238),
    "trpb": dict(csv=DATA / "TRPB_THEMA_Johnston_2024.csv", pos=(183, 184, 227, 228),
                 wt="VFVS", thresh=0.349399),
}

ARMS = {
    "gb1": [
        ("random", "run_20261031_024951_7360aa94", "random"),
        ("ridge", "run_20261031_024951_7360aa94", "ridge"),
        ("random_forest", "run_20261031_024951_7360aa94", "random_forest"),
        ("alde", "run_20261031_025645_7033ab05", "alde"),
        ("reliability", "run_20261031_024308_86773bff", "panoptes_default"),
        ("reliability - dkl,gflownet", "run_20261031_024452_429a9ee4", "panoptes_default"),
        ("frontier", "run_20261031_025342_33e5cac0", "panoptes_qualify_then_optimize"),
    ],
    "trpb": [
        ("random", "run_20261031_025252_7dce6099", "random"),
        ("ridge", "run_20261031_025252_7dce6099", "ridge"),
        ("alde", "run_20261031_030547_ae4ac4bb", "alde"),
        ("reliability", "run_20261031_024533_c41a0125", "panoptes_default"),
        ("frontier", "run_20261031_025827_8dfabb93", "panoptes_qualify_then_optimize"),
    ],
}


def make_parser(positions: Tuple[int, ...], wt: str):
    wtmap = dict(zip(positions, wt))

    def residues(variant_id: str) -> Optional[Tuple[str, ...]]:
        s = variant_id.strip()
        if not s:
            return None
        if s.upper() == "WT":
            return tuple(wtmap[p] for p in positions)
        seen: Dict[int, str] = {}
        for t in s.split(":"):
            if not t:
                continue
            digits = "".join(c for c in t if c.isdigit())
            if not digits or len(t) < 3:
                return None
            p = int(digits)
            if p not in wtmap:
                return None
            seen[p] = t[-1]
        return tuple(seen.get(p, wtmap[p]) for p in positions)

    def depth(g: Tuple[str, ...]) -> int:
        return sum(1 for p, a in zip(positions, g) if a != wtmap[p])

    return residues, depth


report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


emit("=" * 104)
emit("WHY RELIABILITY WINS ON TRPB AND LOSES ON GB1")
emit("=" * 104)
emit("  Held-out seeds 48-95. The ranked candidate pool is NOT in the existing traces,")
emit("  so pool-entry and ranking questions are deferred to the enriched trace.")

# ---------------------------------------------------------------- landscape structure
emit()
emit("=" * 104)
emit("H2: do the two landscapes put their top-1% hits at different depths?")
emit("=" * 104)
land: Dict[str, Dict[str, object]] = {}
for lib, spec in LIB.items():
    residues, depth = make_parser(spec["pos"], spec["wt"])
    hits_by_depth: Counter = Counter()
    all_by_depth: Counter = Counter()
    for r in csv.DictReader(open(spec["csv"], newline="", encoding="utf-8")):
        g = residues(str(r["mutant"]))
        if g is None:
            continue
        try:
            f = float(r["DMS_score"])
        except (TypeError, ValueError):
            continue
        d = depth(g)
        all_by_depth[d] += 1
        if f >= spec["thresh"]:
            hits_by_depth[d] += 1
    land[lib] = dict(hits=hits_by_depth, all=all_by_depth)
    tot = sum(hits_by_depth.values())
    emit()
    emit(f"  {lib.upper()}   {tot:,} top-1% hits")
    emit(f"    {'depth':<7}{'variants':>10}{'hits':>8}{'hit rate':>10}{'% of all hits':>15}")
    for d in sorted(all_by_depth):
        h, n = hits_by_depth.get(d, 0), all_by_depth[d]
        emit(f"    {d:<7}{n:>10,}{h:>8,}{h / n * 100:>9.2f}%{(h / tot * 100 if tot else 0):>14.1f}%")
    mean_hit_depth = (sum(d * c for d, c in hits_by_depth.items()) / tot) if tot else float("nan")
    land[lib]["mean_hit_depth"] = mean_hit_depth
    emit(f"    mean depth of a top-1% hit: {mean_hit_depth:.3f}")

# ---------------------------------------------------------------- what arms selected
emit()
emit("=" * 104)
emit("H1/H2: reachable set and selection depth, against each landscape's hit mass")
emit("=" * 104)
for lib in ("trpb", "gb1"):
    spec = LIB[lib]
    residues, depth = make_parser(spec["pos"], spec["wt"])
    emit()
    emit(f"  {lib.upper()}  (mean depth of a top-1% hit = {land[lib]['mean_hit_depth']:.2f})")
    emit(f"    {'arm':<30}{'reach':>7}{'%land':>8}{'selDepth':>10}{'hitDepth':>10}"
         f"{'bestAvail':>11}{'bestSeen':>10}")
    for label, run, method in ARMS[lib]:
        p = ROOT / run / "selection_trace.csv"
        if not p.exists():
            continue
        union: Set[Tuple[str, ...]] = set()
        sel_depths: List[int] = []
        hit_depths: List[int] = []
        best_seen = float("-inf")
        for r in csv.DictReader(open(p, newline="", encoding="utf-8")):
            if r.get("method") != method:
                continue
            if (r.get("status") or "").strip() != "accepted":
                continue
            g = residues(str(r.get("variant_id", "")))
            if g is None:
                continue
            union.add(g)
            d = depth(g)
            sel_depths.append(d)
            try:
                f = float(r["revealed_fitness"])
            except (TypeError, ValueError):
                continue
            best_seen = max(best_seen, f)
            if f >= spec["thresh"]:
                hit_depths.append(d)
        if not union:
            continue
        # best fitness available anywhere in the reachable set
        avail = float("-inf")
        for r in csv.DictReader(open(spec["csv"], newline="", encoding="utf-8")):
            g = residues(str(r["mutant"]))
            if g is None or g not in union:
                continue
            try:
                avail = max(avail, float(r["DMS_score"]))
            except (TypeError, ValueError):
                pass
        n_land = sum(land[lib]["all"].values())
        emit(f"    {label:<30}{len(union):>7}{len(union) / n_land * 100:>7.2f}%"
             f"{st.fmean(sel_depths):>10.2f}"
             f"{(st.fmean(hit_depths) if hit_depths else float('nan')):>10.2f}"
             f"{avail:>11.4f}{best_seen:>10.4f}")
    emit("    selDepth  = mean mutation depth of everything the arm selected")
    emit("    hitDepth  = mean depth of the selections that were top-1% hits")
    emit("    bestAvail = best fitness anywhere in the reachable set (retrospective)")

# ---------------------------------------------------------------- when does it fail
emit()
emit("=" * 104)
emit("H3: WHEN does GB1 go wrong? New top-1% hits found in each round")
emit("=" * 104)
for lib in ("gb1", "trpb"):
    spec = LIB[lib]
    residues, _d = make_parser(spec["pos"], spec["wt"])
    emit()
    emit(f"  {lib.upper()}   new top-1% hits per round (mean over 48 held-out seeds)")
    emit(f"    {'arm':<30}" + "".join(f"{('r' + str(i)):>8}" for i in range(1, 5)))
    for label, run, method in ARMS[lib]:
        p = ROOT / run / "selection_trace.csv"
        if not p.exists():
            continue
        per: Dict[int, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
        seeds: Set[int] = set()
        for r in csv.DictReader(open(p, newline="", encoding="utf-8")):
            if r.get("method") != method:
                continue
            if (r.get("status") or "").strip() != "accepted":
                continue
            try:
                rnd = int(float(r["round"]))
                seed = int(float(r["seed"]))
                f = float(r["revealed_fitness"])
            except (TypeError, ValueError):
                continue
            if rnd == 0:
                continue
            seeds.add(seed)
            if f >= spec["thresh"]:
                per[rnd][seed] += 1
        if not seeds:
            continue
        cells = []
        for rnd in range(1, 5):
            vals = [per[rnd].get(s, 0) for s in sorted(seeds)]
            cells.append(f"{st.fmean(vals):>8.2f}" if vals else f"{'-':>8}")
        emit(f"    {label:<30}" + "".join(cells))

(SCRATCH / "gb1_mechanism_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'gb1_mechanism_report.txt'}")
