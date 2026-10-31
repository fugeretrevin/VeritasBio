"""Budget sensitivity, computed from runs that already exist.

The reviewer asks for budget sensitivity at 25/50/75/100/125 measurements as new work.
It needs no new compute: round_metrics.csv already records cumulative top-1% hits at
every round, so the whole curve is in the runs we have.

This matters commercially more than its placement in the plan suggests. 125 measurements
is a large campaign. If the engine's advantage only appears at 125 and is absent at 50,
the result is much less useful than it reads, because a lab deciding whether to buy is
deciding at a budget it chooses, not at ALDE's.

Held-out seeds only (48-95), so nothing here is tuned on.
"""
from __future__ import annotations

import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
TOP1 = "cum_top1_hits_including_initial"
BUDGETS = (25, 50, 75, 100, 125)
T_CRIT = 2.01


def library_of(run: Path) -> Optional[str]:
    try:
        man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    except Exception:
        return None
    p = str(man.get("dataset", {}).get("path", "")).upper()
    return "trpb" if "TRPB" in p else ("gb1" if "GB1" in p else None)


def base_seed(run: Path) -> Optional[int]:
    try:
        man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    except Exception:
        return None
    sv = man.get("config", {}).get("seed_values", [])
    return sv[0] if sv else None


def ablated(run: Path) -> str:
    rt = run / "report.txt"
    if not rt.exists():
        return ""
    txt = rt.read_text(encoding="utf-8", errors="replace")
    for line in txt.splitlines():
        if "ABLATED_COMPONENTS=" in line:
            return line.split("ABLATED_COMPONENTS=", 1)[1].strip()
    return ""


def curve(run: Path, method: str) -> Dict[int, Dict[int, float]]:
    """budget -> {seed: cumulative top-1% hits at or before that budget}"""
    out: Dict[int, Dict[int, float]] = {b: {} for b in BUDGETS}
    p = run / "round_metrics.csv"
    if not p.exists():
        return out
    with open(p, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            try:
                seed = int(float(r["seed"]))
                n = int(float(r["n_unique_measurements"]))
                v = float(r[TOP1])
            except (TypeError, ValueError, KeyError):
                continue
            for b in BUDGETS:
                if n <= b:
                    prev = out[b].get(seed)
                    if prev is None or v >= prev:
                        out[b][seed] = v
    return out


# Collect held-out arms of interest.
arms: Dict[Tuple[str, str], Dict[int, Dict[int, float]]] = {}
for run in sorted(ROOT.glob("run_2026103*")):
    lib = library_of(run)
    if lib is None or base_seed(run) != 48:
        continue
    man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    abl = ablated(run)
    for method, info in man.get("methods", {}).items():
        if info.get("status") != "available":
            continue
        if method == "panoptes_default":
            label = "reliability" if not abl else f"reliability \u2212 {abl}"
        elif method == "panoptes_qualify_then_optimize":
            label = "frontier" if not abl else f"frontier \u2212 {abl}"
        else:
            label = method
        c = curve(run, method)
        if any(c[b] for b in BUDGETS):
            arms.setdefault((lib, label), c)

report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


ORDER = ["random", "ridge", "random_forest", "alde", "reliability",
         "reliability \u2212 dkl,gflownet", "frontier", "frontier \u2212 dkl,gflownet"]

for lib in ("gb1", "trpb"):
    present = [l for l in ORDER if (lib, l) in arms]
    present += sorted(l for (lb, l) in arms if lb == lib and l not in present)
    emit()
    emit("=" * 96)
    emit(f"{lib.upper()}  —  mean cumulative top-1% hits by budget, held-out seeds 48-95")
    emit("=" * 96)
    emit(f"  {'arm':<30}" + "".join(f"{b:>10}" for b in BUDGETS))
    for label in present:
        c = arms[(lib, label)]
        cells = []
        for b in BUDGETS:
            vals = list(c[b].values())
            cells.append(f"{st.fmean(vals):>10.2f}" if vals else f"{'-':>10}")
        emit(f"  {label:<30}" + "".join(cells))

    # The commercially important question: at which budget does our best arm's lead
    # over the best baseline become significant?
    best_base = None
    for cand in ("alde", "ridge", "random_forest"):
        if (lib, cand) in arms:
            v = st.fmean(arms[(lib, cand)][125].values())
            if best_base is None or v > best_base[1]:
                best_base = (cand, v)
    ours = [l for l in present if l.startswith(("reliability", "frontier"))]
    if best_base and ours:
        emit()
        emit(f"  Lead over the strongest held-out baseline ({best_base[0]}, "
             f"{best_base[1]:.2f} at 125), paired per budget:")
        emit(f"  {'arm':<30}" + "".join(f"{b:>14}" for b in BUDGETS))
        for label in ours:
            cells = []
            for b in BUDGETS:
                a, bb = arms[(lib, label)][b], arms[(lib, best_base[0])][b]
                seeds = sorted(set(a) & set(bb))
                d = [a[s] - bb[s] for s in seeds]
                if len(d) < 3:
                    cells.append(f"{'-':>14}")
                    continue
                m, sd = st.fmean(d), st.stdev(d)
                t = m / (sd / np.sqrt(len(d))) if sd > 0 else float("nan")
                mark = "*" if abs(t) >= T_CRIT else " "
                cells.append(f"{m:>+9.2f}{mark}{('t' + format(t, '+.1f')):>4}")
            emit(f"  {label:<30}" + "".join(cells))
        emit("    * |t| >= 2.01")

(SCRATCH / "budget_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"\nwrote {SCRATCH / 'budget_report.txt'}")

