"""Does the TrpB advantage represent many distinct solutions, or one dense basin?

Workstream 5 of the review, computed from selection_trace.csv files that already exist.
No new campaigns.

The question matters because a raw top-1% count treats 50 variants in one mutational
neighbourhood as 50 wins. A lab cannot use them that way: if they share a basin they
share whatever liability sits in that basin, and the portfolio is really one candidate.

Method:
  - a hit is a measured variant at or above the run's own top-1% threshold, read from
    the manifest so it is byte-identical to the threshold the benchmark scored with
  - every variant is a 4-site genotype, so a sequence is the 4 mutant residues and
    Hamming distance between two hits is how many of the 4 sites differ (0-4)
  - diversity-adjusted yield: sort hits by fitness descending and greedily keep a hit
    only when it is at least d away from everything already kept. This is the review's
    rule, applied after measurement, so it cannot advantage a method that knew about it
  - the SAME greedy rule is applied to the full top-1% set of the landscape, which gives
    the ceiling: the most distinct hits anybody could report at that threshold. Without
    it, raw counts are being compared against an unknown maximum

Thresholds d=2,3,4 are all reported because on a 4-site library the choice of d is not
obvious and the conclusion should not depend on one pick. Hamming separation is a
sequence-level proxy; it is not evidence of functional or biophysical independence.
"""
from __future__ import annotations

import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
DISTANCES = (2, 3, 4)
T_CRIT = 2.01

HELDOUT = {
    "gb1": {
        "random": ("run_20261031_024951_7360aa94", "random"),
        "ridge": ("run_20261031_024951_7360aa94", "ridge"),
        "random_forest": ("run_20261031_024951_7360aa94", "random_forest"),
        "alde": ("run_20261031_025645_7033ab05", "alde"),
        "reliability": ("run_20261031_024308_86773bff", "panoptes_default"),
        "reliability - dkl,gflownet": ("run_20261031_024452_429a9ee4", "panoptes_default"),
        "frontier": ("run_20261031_025342_33e5cac0", "panoptes_qualify_then_optimize"),
        "frontier - dkl,gflownet": ("run_20261031_025435_b451ee1a", "panoptes_qualify_then_optimize"),
        "frontier - dkl,dpp,gflownet": ("run_20261031_025516_ce3387ff", "panoptes_qualify_then_optimize"),
    },
    "trpb": {
        "random": ("run_20261031_025252_7dce6099", "random"),
        "ridge": ("run_20261031_025252_7dce6099", "ridge"),
        "random_forest": ("run_20261031_025252_7dce6099", "random_forest"),
        "alde": ("run_20261031_030547_ae4ac4bb", "alde"),
        "reliability": ("run_20261031_024533_c41a0125", "panoptes_default"),
        "reliability - dkl,gflownet": ("run_20261031_024735_a7e2cd7c", "panoptes_default"),
        "frontier": ("run_20261031_025827_8dfabb93", "panoptes_qualify_then_optimize"),
        "frontier - dkl,gflownet": ("run_20261031_025926_eada00c6", "panoptes_qualify_then_optimize"),
        "frontier - dkl,dpp,gflownet": ("run_20261031_030003_74f8caa3", "panoptes_qualify_then_optimize"),
    },
}
DATASETS = {
    "gb1": Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data\GB1_STRSG_Wu_2016.csv"),
    "trpb": Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data\TRPB_THEMA_Johnston_2024.csv"),
}


def tokens(variant_id: str) -> Optional[List[Tuple[int, str, str]]]:
    """'V183N:F184R' -> [(183,'V','N'), (184,'F','R')]. 'WT' -> [].

    Mutation strings omit positions that carry the wild-type residue, so
    'V183A:F184A:V227A' is a 4-site genotype whose 4th site is WT. Treating the
    string length as the genotype length would compare a 3-tuple against a
    4-tuple; zip() would silently truncate and understate real distances.
    """
    s = variant_id.strip()
    if not s:
        return None
    if s.upper() == "WT":
        return []
    out: List[Tuple[int, str, str]] = []
    for t in s.split(":"):
        if not t:
            continue
        digits = "".join(c for c in t if c.isdigit())
        if not digits or len(t) < 3:
            return None
        out.append((int(digits), t[0], t[-1]))
    return out


def canonical(variant_id: str, positions: Tuple[int, ...], wt: Dict[int, str]
              ) -> Optional[Tuple[str, ...]]:
    """Residues at every campaign position, wild-type where the string omits one."""
    toks = tokens(variant_id)
    if toks is None:
        return None
    seen = {}
    for pos, _w, mut in toks:
        if pos not in wt:
            return None          # not part of this library
        seen[pos] = mut
    if any(p not in wt for p in seen):
        return None
    return tuple(seen.get(p, wt[p]) for p in positions)


def hamming(a: Sequence[str], b: Sequence[str]) -> int:
    if len(a) != len(b):
        raise ValueError(f"hamming on unequal lengths: {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def dataset_library(lib: str) -> Tuple[Tuple[int, ...], Dict[int, str]]:
    """The library's varied positions and their wild-type residues, from the data.

    Positions are whichever ones the combinatorial library varies; a record that
    touches any other position is not a library member.
    """
    counts: Dict[int, int] = defaultdict(int)
    wt: Dict[int, str] = {}
    with open(DATASETS[lib], newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        vid = next((c for c in (rdr.fieldnames or []) if c.lower() == "mutant"), "mutant")
        for r in rdr:
            toks = tokens(str(r.get(vid, "")))
            if not toks:
                continue
            for pos, w, _m in toks:
                counts[pos] += 1
                wt.setdefault(pos, w)
    # the combinatorial sites are the heavily-varied ones; stray DMS positions are rare
    if not counts:
        return (), {}
    cut = max(counts.values()) * 0.25
    pos = tuple(sorted(p for p, n in counts.items() if n >= cut))
    return pos, {p: wt[p] for p in pos}


def greedy_diverse(hits: List[Tuple[float, Tuple[str, ...]]], d: int) -> int:
    """Hits sorted by fitness desc; keep one only if >= d from all kept."""
    kept: List[Tuple[str, ...]] = []
    for _, g in hits:
        if all(hamming(g, k) >= d for k in kept):
            kept.append(g)
    return len(kept)


def manifest_top1(run: str) -> Tuple[float, int, int]:
    m = json.loads((ROOT / run / "manifest.json").read_text(encoding="utf-8"))
    t = m["top1"]
    return float(t["threshold"]), int(t["n_hits"]), int(m["dataset"]["unique_variants"])


def load_universe(lib: str, threshold: float, positions: Tuple[int, ...],
                  wt: Dict[int, str]) -> Tuple[List[Tuple[float, Tuple[str, ...]]], int]:
    """Every top-1% library member in the landscape, for the ceiling.

    Returns (hits, n_excluded) where excluded records touch positions outside
    the combinatorial library and so are not selectable by any arm.
    """
    out: List[Tuple[float, Tuple[str, ...]]] = []
    excluded = 0
    with open(DATASETS[lib], newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        cols = rdr.fieldnames or []
        vid = next((c for c in cols if c.lower() == "mutant"), None)
        fit = next((c for c in cols if c.lower() == "dms_score"), None)
        if vid is None or fit is None:
            print(f"  ! {lib}: no mutant/DMS_score column in {cols}", file=sys.stderr)
            return out, 0
        for r in rdr:
            try:
                f = float(r[fit])
            except (TypeError, ValueError):
                continue
            if f < threshold:
                continue
            g = canonical(str(r[vid]), positions, wt)
            if g is None:
                excluded += 1
                continue
            out.append((f, g))
    out.sort(key=lambda p: -p[0])
    return out, excluded


def trace_hits(run: str, method: str, threshold: float,
               positions: Tuple[int, ...], wt: Dict[int, str]
               ) -> Tuple[Dict[int, List[Tuple[float, Tuple[str, ...]]]],
                          Dict[int, List[Tuple[float, Tuple[str, ...]]]]]:
    """-> (all measured hits per seed, policy-selected-only hits per seed)"""
    allh: Dict[int, Dict[Tuple[str, ...], float]] = defaultdict(dict)
    selh: Dict[int, Dict[Tuple[str, ...], float]] = defaultdict(dict)
    p = ROOT / run / "selection_trace.csv"
    with open(p, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            status = (r.get("status") or "").strip()
            if status in ("rejected", "considered", ""):
                continue
            try:
                f = float(r["revealed_fitness"])
                seed = int(float(r["seed"]))
            except (TypeError, ValueError):
                continue
            if f < threshold:
                continue
            g = canonical(str(r.get("variant_id", "")), positions, wt)
            if g is None:
                continue
            allh[seed][g] = max(allh[seed].get(g, f), f)
            if status != "initial":
                selh[seed][g] = max(selh[seed].get(g, f), f)
    def fin(d):
        return {s: sorted(((f, g) for g, f in m.items()), key=lambda p: -p[0])
                for s, m in d.items()}
    return fin(allh), fin(selh)


report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


def paired_t(a: Dict[int, float], b: Dict[int, float]) -> Tuple[float, float, int]:
    seeds = sorted(set(a) & set(b))
    d = [a[s] - b[s] for s in seeds]
    if len(d) < 3:
        return float("nan"), float("nan"), len(d)
    m, sd = st.fmean(d), st.stdev(d)
    return m, (m / (sd / np.sqrt(len(d))) if sd > 0 else float("nan")), len(d)


emit("=" * 104)
emit("IS THE ADVANTAGE MANY DISTINCT SOLUTIONS, OR ONE DENSE BASIN?")
emit("=" * 104)
emit("  Held-out seeds 48-95. Hits = measured variants at/above each run's own top-1% threshold.")
emit("  Hamming distance is over the 4 variable sites. Sequence separation is a proxy, not")
emit("  evidence of functional independence.")

ORDER = ["random", "ridge", "random_forest", "alde", "reliability",
         "reliability - dkl,gflownet", "frontier", "frontier - dkl,gflownet",
         "frontier - dkl,dpp,gflownet"]

for lib in ("trpb", "gb1"):
    arms = HELDOUT[lib]
    any_run, any_method = next(iter(arms.values()))
    threshold, n_top1, n_univ = manifest_top1(any_run)
    positions, wt = dataset_library(lib)
    universe, excluded = load_universe(lib, threshold, positions, wt)
    emit()
    emit("=" * 104)
    emit(f"{lib.upper()}   threshold {threshold:.6g}   "
         f"{n_top1:,} of {n_univ:,} variants are top-1%")
    emit(f"  library varies {positions} (WT "
         f"{''.join(wt[p] for p in positions)}); {len(universe):,} top-1% library members"
         + (f", {excluded:,} top-1% records outside the library excluded" if excluded else ""))
    emit("=" * 104)

    ceilings = {d: greedy_diverse(universe, d) for d in DISTANCES}
    emit(f"  CEILING - the same greedy rule on the whole top-1% set:   "
         + "   ".join(f"d>={d}: {ceilings[d]}" for d in DISTANCES))
    emit("  No method can report more diversity-adjusted hits than this, at any budget.")
    emit()

    per_arm_all: Dict[str, Dict[int, List[Tuple[float, Tuple[str, ...]]]]] = {}
    div: Dict[int, Dict[str, Dict[int, float]]] = {d: {} for d in DISTANCES}
    raw: Dict[str, Dict[int, float]] = {}

    emit(f"  {'arm':<30}{'hits':>7}{'meanHam':>9}{'1-clust':>9}{'%largest':>10}"
         + "".join(f"{('d>=' + str(d)):>8}" for d in DISTANCES)
         + f"{'redund':>9}{'%ceil d2':>10}")
    for label in ORDER:
        if label not in arms:
            continue
        run, method = arms[label]
        allh, _selh = trace_hits(run, method, threshold, positions, wt)
        if not allh:
            continue
        per_arm_all[label] = allh
        seeds = sorted(allh)
        raw[label] = {s: float(len(allh[s])) for s in seeds}
        hams, nclust, fraclarge = [], [], []
        for s in seeds:
            gs = [g for _, g in allh[s]]
            if len(gs) >= 2:
                hams.append(st.fmean(hamming(gs[i], gs[j])
                                     for i in range(len(gs)) for j in range(i + 1, len(gs))))
            # single-linkage components joining hits that differ at <= 1 site
            parent = list(range(len(gs)))
            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            for i in range(len(gs)):
                for j in range(i + 1, len(gs)):
                    if hamming(gs[i], gs[j]) <= 1:
                        ri, rj = find(i), find(j)
                        if ri != rj:
                            parent[ri] = rj
            if gs:
                comp = defaultdict(int)
                for i in range(len(gs)):
                    comp[find(i)] += 1
                nclust.append(len(comp))
                fraclarge.append(max(comp.values()) / len(gs))
        for d in DISTANCES:
            div[d][label] = {s: float(greedy_diverse(allh[s], d)) for s in seeds}
        mh = st.fmean(hams) if hams else float("nan")
        mraw = st.fmean(raw[label].values())
        md2 = st.fmean(div[2][label].values())
        emit(f"  {label:<30}{mraw:>7.2f}{mh:>9.2f}"
             f"{(st.fmean(nclust) if nclust else float('nan')):>9.2f}"
             f"{(st.fmean(fraclarge) * 100 if fraclarge else float('nan')):>9.1f}%"
             + "".join(f"{st.fmean(div[d][label].values()):>8.2f}" for d in DISTANCES)
             + f"{(mraw / md2 if md2 else float('nan')):>8.2f}x"
             + f"{(md2 / ceilings[2] * 100 if ceilings[2] else float('nan')):>9.1f}%")

    emit()
    emit("    1-clust = single-linkage groups joining hits that differ at <=1 of 4 sites")
    emit("    redund  = raw hits per diversity-adjusted hit at d>=2; 1.00 means every hit is distinct")

    # Does the lead survive the diversity adjustment?
    best_base = max(("ridge", "random_forest", "alde"),
                    key=lambda c: st.fmean(raw[c].values()) if c in raw else -1)
    emit()
    emit(f"  Lead over the strongest held-out baseline ({best_base}), raw vs diversity-adjusted:")
    emit(f"  {'arm':<30}{'raw':>16}" + "".join(f"{('d>=' + str(d)):>16}" for d in DISTANCES))
    for label in ORDER:
        if label not in per_arm_all or label.startswith(("random", "ridge", "alde")) \
                or label == best_base:
            continue
        cells = []
        m, t, _ = paired_t(raw[label], raw[best_base])
        cells.append(f"{m:>+10.2f}{('t' + format(t, '+.1f')):>6}")
        for d in DISTANCES:
            m2, t2, _ = paired_t(div[d][label], div[d][best_base])
            mark = "*" if abs(t2) >= T_CRIT else " "
            cells.append(f"{m2:>+9.2f}{mark}{('t' + format(t2, '+.1f')):>6}")
        emit(f"  {label:<30}" + "".join(cells))
    emit("    * |t| >= 2.01")

(SCRATCH / "diversity_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'diversity_report.txt'}")
