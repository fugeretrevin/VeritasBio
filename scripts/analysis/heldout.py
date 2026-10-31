"""Held-out confirmation of the configuration chosen from the 2^5 factorial, plus
honest interval estimates for the extreme-tail rates.

Two separate jobs.

1. The configuration "remove dkl and gflownet" was picked by scanning 32 subsets on seeds
   0-47, so its advantage on those seeds is optimistically biased by selection. Seeds
   48-95 were never looked at during the scan. The gain has to survive there.

2. The extreme-tail rates have been quoted as point estimates (4%, 0%, 25%, 54%), but at
   48 seeds a 4% rate is TWO campaigns. Clopper-Pearson intervals are computed here so
   the rates are reported with the precision they actually have. This was a fair external
   criticism and it is not cosmetic: it decides whether the tail gap is one finding or
   two coin flips.
"""
from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats as sps

sys.path.insert(0, r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio")
import benchmark_panoptes_isolated as iso

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
SCRATCH = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
               r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
               r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad")
TOP1 = "cum_top1_hits_including_initial"
BEST = "best_measured_fitness"
T_CRIT = 2.01

# Dev-seed (0-47) runs, for the side-by-side.
DEV = {("gb1", "shipped"): "run_20261031_023520_b5ae03df",
       ("trpb", "shipped"): "run_20261030_233220_8e467def"}
DEV_BASE = {"gb1": "run_20261030_235620_fd18d20e",
            "trpb": "run_20261030_233220_8e467def"}


def final_rows(run: str, method: str) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    p = ROOT / run / "round_metrics.csv"
    if not p.exists():
        return out
    with open(p, newline="", encoding="utf-8") as fh:
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
    # Paired bootstrap CI alongside the t, since campaign outcomes are skewed.
    rng = np.random.default_rng(12345)
    arr = np.asarray(d)
    boots = arr[rng.integers(0, len(arr), size=(20000, len(arr)))].mean(axis=1)
    return {"n": len(d), "mean": mean, "t": (mean / se) if se > 0 else float("nan"),
            "lo": float(np.percentile(boots, 2.5)),
            "hi": float(np.percentile(boots, 97.5)),
            "median": float(np.median(arr)),
            "wins": sum(1 for x in d if x > 0), "losses": sum(1 for x in d if x < 0)}


def library_of(run: Path) -> Optional[str]:
    import json
    try:
        man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    except Exception:
        return None
    path = str(man.get("dataset", {}).get("path", "")).upper()
    return "trpb" if "TRPB" in path else ("gb1" if "GB1" in path else None)


def seeds_of(run: Path) -> Tuple[int, ...]:
    import json
    man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    return tuple(man.get("config", {}).get("seed_values", ()))


# --- locate the four held-out runs by their seed block and arm note ---
heldout: Dict[Tuple[str, str], str] = {}
for run in sorted(ROOT.glob("run_2026103*")):
    lib = library_of(run)
    if lib is None:
        continue
    try:
        sv = seeds_of(run)
    except Exception:
        continue
    if not sv or min(sv) != 48:
        continue
    rt = run / "report.txt"
    txt = rt.read_text(encoding="utf-8", errors="replace") if rt.exists() else ""
    if "panoptes_default" not in txt:
        continue
    arm = "minus_dkl_gflownet" if "ABLATED_COMPONENTS=dkl,gflownet" in txt else "shipped"
    heldout[(lib, arm)] = run.name

print("held-out runs found:")
for k, v in sorted(heldout.items()):
    print(f"  {k}  {v}")

report: List[str] = []


def emit(line: str = "") -> None:
    print(line)
    report.append(line)


emit()
emit("=" * 104)
emit("HELD-OUT CONFIRMATION — 'remove dkl,gflownet', seeds 48-95 (never seen during the scan)")
emit("=" * 104)
emit(f"  {'landscape':<12}{'shipped':>9}{'minus':>9}{'delta':>8}{'t':>7}"
     f"{'95% bootstrap CI':>22}{'median d':>10}{'won/lost':>10}")
for lib in ("gb1", "trpb"):
    if (lib, "shipped") not in heldout or (lib, "minus_dkl_gflownet") not in heldout:
        emit(f"  {lib}: runs not found")
        continue
    a = series(heldout[(lib, "minus_dkl_gflownet")], "panoptes_default", TOP1)
    b = series(heldout[(lib, "shipped")], "panoptes_default", TOP1)
    c = paired(a, b)
    emit(f"  {lib.upper():<12}{st.fmean(b.values()):>9.2f}{st.fmean(a.values()):>9.2f}"
         f"{c['mean']:>+8.2f}{c['t']:>7.2f}"
         f"{('[' + format(c['lo'], '+.2f') + ', ' + format(c['hi'], '+.2f') + ']'):>22}"
         f"{c['median']:>+10.1f}{str(c['wins']) + '/' + str(c['losses']):>10}")
emit()
emit("  Dev seeds 0-47 gave GB1 +10.08 (t=5.93) and TrpB +5.35 (t=5.69). If the held-out")
emit("  deltas are close to those, the choice was not an artefact of scanning 32 subsets.")

# ---------------------------------------------------- tail rates with intervals
emit()
emit("=" * 104)
emit("EXTREME TAIL, WITH INTERVALS — P(reach top 0.01% of the eligible universe)")
emit("=" * 104)
emit("  Clopper-Pearson exact 95% intervals on 48 campaigns. k is the raw number of")
emit("  campaigns that got there; the point estimates alone are misleadingly precise.")
emit()
for lib in ("gb1", "trpb"):
    iso.configure_library(lib)
    universe = iso.load_universe(
        dataset_path=iso.DEFAULT_DATASET, provenance_mode="acknowledge-unflagged",
        imputed_ids_path=None, max_variants=None, subsample_seed=0)
    values = np.array(list(universe.fitness.values()), dtype=float)
    t9999 = float(np.percentile(values, 99.99))
    n_above = int((values >= t9999).sum())
    emit(f"  {lib.upper()}  threshold {t9999:.4f}  ({n_above} of {len(values):,} variants "
         f"are at or above it)")
    emit(f"    {'arm':<34}{'k/n':>8}{'rate':>8}{'exact 95% CI':>20}")
    arms: List[Tuple[str, Dict[int, float]]] = []
    for m in ("random", "ridge", "random_forest", "alde",
              "panoptes_default", "panoptes_qualify_then_optimize"):
        s = series(DEV_BASE[lib], m, BEST)
        if s:
            arms.append((m, s))
    if (lib, "minus_dkl_gflownet") in heldout:
        arms.append(("panoptes minus dkl,gflownet (held-out)",
                     series(heldout[(lib, "minus_dkl_gflownet")], "panoptes_default", BEST)))
    for label, s in arms:
        vals = list(s.values())
        k, n = sum(1 for v in vals if v >= t9999), len(vals)
        lo, hi = sps.beta.ppf(0.025, k, n - k + 1) if k else 0.0, \
                 sps.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
        emit(f"    {label:<34}{str(k) + '/' + str(n):>8}{100 * k / n:>7.0f}%"
             f"{('[' + format(100 * lo, '.0f') + '%, ' + format(100 * hi, '.0f') + '%]'):>20}")
    # Does the ALDE-vs-ours tail gap survive a test that respects the small counts?
    ours = [v for v in series(DEV_BASE[lib], "panoptes_default", BEST).values()]
    theirs = [v for v in series(DEV_BASE[lib], "alde", BEST).values()]
    if ours and theirs:
        ko = sum(1 for v in ours if v >= t9999)
        kt = sum(1 for v in theirs if v >= t9999)
        # Paired: same seeds, so McNemar on the discordant pairs is the right test.
        so = series(DEV_BASE[lib], "panoptes_default", BEST)
        stt = series(DEV_BASE[lib], "alde", BEST)
        both = sorted(set(so) & set(stt))
        b_only = sum(1 for s in both if so[s] >= t9999 and stt[s] < t9999)
        c_only = sum(1 for s in both if so[s] < t9999 and stt[s] >= t9999)
        p = sps.binomtest(b_only, b_only + c_only, 0.5).pvalue if (b_only + c_only) else 1.0
        emit(f"    reliability {ko}/48 vs ALDE {kt}/48 — paired McNemar on "
             f"{b_only + c_only} discordant seeds, p={p:.4f}")
    emit()

(SCRATCH / "heldout_report.txt").write_text("\n".join(report), encoding="utf-8")
print(f"wrote {SCRATCH / 'heldout_report.txt'}")
