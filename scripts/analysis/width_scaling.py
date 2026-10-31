"""O-12: is the policy's depth bias an absolute constant or a fraction of library width?

This is the highest-value open question in the project and it costs no compute. Everything
in the mechanism story rests on "reliability selects at mean depth ~2.6", which was
measured on TWO libraries of IDENTICAL width (TrpB and GB1, both 4 sites). Two libraries
of the same width cannot distinguish a constant from a proportion.

  absolute      -> on a 10-site library the policy still selects near depth 2.6, i.e. 26%
                   of sites. It gets effectively shallower as libraries widen, and should
                   degrade monotonically with site count.
  proportional  -> it selects near 65% of sites (2.6/4), i.e. depth ~6.5 of 10. On wide
                   libraries it is NOT a shallow exploiter, the mechanism describes only
                   the two 4-site development landscapes, and it applies to ZERO
                   independent targets. That would void the mechanism's relevance rather
                   than merely weaken the evidence for it.

The HER2 run supplies the missing width: 10 sites, reliability, 48 seeds.

ONE CORRECTION THE NAIVE TEST NEEDS. HER2 runs pool-restricted, so selections come from
an enumerated pool of measured variants rather than free generation. The pool's own depth
distribution constrains what CAN be selected. Comparing raw selected depth across
libraries would therefore confound policy preference with what each universe offers. So
the quantity of interest is selected depth RELATIVE to available depth:

    ratio = mean selected depth / mean universe depth

Under the absolute reading that ratio must FALL as width grows (a constant numerator over
a growing denominator). Under the proportional reading it stays roughly flat.
"""
from __future__ import annotations

import csv
import statistics as st
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
DATA = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data")

LIBS = {
    "TrpB": dict(csv=DATA / "TRPB_THEMA_Johnston_2024.csv", pos=(183, 184, 227, 228),
                 wt="VFVS", stage=None),
    "GB1": dict(csv=DATA / "GB1_STRSG_Wu_2016.csv", pos=(39, 40, 41, 54),
                wt="VDGV", stage=None),
    "HER2": dict(csv=DATA / "HER2_TRASTUZUMAB_Mason_2021.csv",
                 pos=(99, 100, 101, 102, 103, 104, 105, 106, 107, 108),
                 wt="WGGDGFYAMD", stage="train"),
}

# (label, library, run, method)
ARMS: List[Tuple[str, str, str, str]] = [
    ("reliability", "TrpB", "run_20261031_024533_c41a0125", "panoptes_default"),
    ("reliability", "GB1", "run_20261031_024308_86773bff", "panoptes_default"),
    ("reliability", "HER2", "run_20261031_051736_45ae3dd8", "panoptes_default"),
    ("ridge", "TrpB", "run_20261031_025252_7dce6099", "ridge"),
    ("ridge", "GB1", "run_20261031_024951_7360aa94", "ridge"),
    ("ridge", "HER2", "run_20261031_051736_45ae3dd8", "ridge"),
    ("random", "HER2", "run_20261031_051736_45ae3dd8", "random"),
    ("frontier", "HER2", "run_20261031_051736_45ae3dd8",
     "panoptes_qualify_then_optimize"),
]


def make(pos: Tuple[int, ...], wt: str):
    wtmap = dict(zip(pos, wt))

    def residues(vid: str) -> Optional[Tuple[str, ...]]:
        s = vid.strip()
        if not s:
            return None
        if s.upper() == "WT":
            return tuple(wtmap[p] for p in pos)
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
        return tuple(seen.get(p, wtmap[p]) for p in pos)

    def depth(g: Tuple[str, ...]) -> int:
        return sum(1 for p, a in zip(pos, g) if a != wtmap[p])

    return residues, depth


print("=" * 100)
print("O-12: DOES THE DEPTH BIAS SCALE WITH LIBRARY WIDTH?")
print("=" * 100)

# ---- universe depth distribution per library
uni: Dict[str, Dict[str, float]] = {}
for name, spec in LIBS.items():
    residues, depth = make(spec["pos"], spec["wt"])
    depths: List[int] = []
    with open(spec["csv"], newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if spec["stage"] is not None and (r.get("stage") or "") != spec["stage"]:
                continue
            vid = str(r["mutant"])
            if vid.strip().upper() == "WT":
                continue                      # parent is not a library member
            g = residues(vid)
            if g is None:
                continue
            depths.append(depth(g))
    n_sites = len(spec["pos"])
    uni[name] = dict(sites=n_sites, n=len(depths), mean=st.fmean(depths),
                     frac=st.fmean(depths) / n_sites)
    print(f"  {name:<6} {n_sites:>3} sites  universe n={len(depths):>7,}  "
          f"mean depth {st.fmean(depths):>6.3f}  = {st.fmean(depths)/n_sites:>5.1%} of sites")

print()
print(f"  {'arm':<12}{'lib':<6}{'sites':>6}{'selDepth':>10}{'as % sites':>12}"
      f"{'uniDepth':>10}{'sel/uni':>9}{'shortfall':>11}")

rows: Dict[str, List[Tuple[int, float, float]]] = {}
for label, lib, run, method in ARMS:
    spec = LIBS[lib]
    residues, depth = make(spec["pos"], spec["wt"])
    p = ROOT / run / "selection_trace.csv"
    if not p.exists():
        print(f"  {label:<12}{lib:<6}  (missing {run})")
        continue
    sel: List[int] = []
    with open(p, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            if (r.get("status") or "").strip() != "accepted":
                continue
            try:
                if int(float(r["round"])) == 0:
                    continue          # shared initial sample is not a policy choice
            except (TypeError, ValueError):
                continue
            g = residues(str(r.get("variant_id", "")))
            if g is None:
                continue
            sel.append(depth(g))
    if not sel:
        continue
    u = uni[lib]
    m = st.fmean(sel)
    ratio = m / u["mean"]
    print(f"  {label:<12}{lib:<6}{u['sites']:>6}{m:>10.3f}{m/u['sites']:>11.1%}"
          f"{u['mean']:>10.3f}{ratio:>9.3f}{u['mean']-m:>11.3f}")
    rows.setdefault(label, []).append((u["sites"], m, ratio))

print()
print("  selDepth  = mean mutation depth of policy-selected variants (round 0 excluded)")
print("  uniDepth  = mean depth of the eligible universe, i.e. what was on offer")
print("  sel/uni   = the width-invariant quantity. Falls with width => ABSOLUTE.")
print("              Stays flat => PROPORTIONAL.")

print()
print("=" * 100)
print("VERDICT")
print("=" * 100)
for label, obs in rows.items():
    if len(obs) < 2:
        continue
    obs.sort()
    narrow = [o for o in obs if o[0] <= 4]
    wide = [o for o in obs if o[0] > 4]
    if not narrow or not wide:
        continue
    rn, rw = st.fmean(o[2] for o in narrow), st.fmean(o[2] for o in wide)
    dn, dw = st.fmean(o[1] for o in narrow), st.fmean(o[1] for o in wide)
    print(f"  {label}:")
    print(f"    4-site  : selected depth {dn:.3f}, sel/uni ratio {rn:.3f}")
    print(f"    10-site : selected depth {dw:.3f}, sel/uni ratio {rw:.3f}")
    if abs(dw - dn) < 1.0:
        print(f"    -> selected depth barely moved ({dn:.2f} -> {dw:.2f}) while width "
              f"went 4 -> 10: consistent with ABSOLUTE")
    elif abs(rw - rn) < 0.12:
        print(f"    -> the ratio held ({rn:.2f} -> {rw:.2f}) while width went 4 -> 10: "
              f"consistent with PROPORTIONAL")
    else:
        print(f"    -> neither clean: depth {dn:.2f}->{dw:.2f}, ratio {rn:.2f}->{rw:.2f}")
