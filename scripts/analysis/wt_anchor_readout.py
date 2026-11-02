"""The three wild-type anchor modes, paired per landscape, on BOTH endpoints.

Validity FIRST, by construction rather than by judgement. A run folder is admitted only
when every one of these holds; otherwise it is discarded and the reason is printed:

  V1  `manifest.json` exists. An artifact without one is a partial launch, not data -
      a task-stop is not a process kill, and that distinction has already voided ten run
      directories in this project.
  V2  `run.state == "SUCCESS"`.
  V3  48 seeds completed for every arm quoted.
  V4  every seed of every quoted arm spent the FULL budget (125 unique measurements).
  V5  `run.wt_anchor_verified.agrees_with_configuration` is true and `vacuous` is false -
      the mode is read from the engine's own return value, not from the CLI flag
      (rule 2.12).

Then the pairing premise, also checked rather than assumed:

  P1  `configuration.seeding.seed_block_sha256` identical across the three folders.
  P2  `initial_sets_sha256` identical per seed across the three folders.
  P3  the top-1% threshold, READ from each manifest (rule 2.3), identical across folders.
  P4  `ridge` and `random` numerically identical across folders, cell for cell. Those arms
      never read the anchor, so a difference means the arms are not paired.

TWO ENDPOINTS, always both:
  raw   top-1% hits excluding the initial sample, taken from `round_metrics.csv` and
        INDEPENDENTLY recomputed from `selection_trace.csv` as a cross-check.
  d>=2  the pre-registered diversity adjustment: hits fitness-descending, keep one only
        when it is Hamming >= 2 from everything already kept. A raw-only contrast cannot
        tell "found more good variants" from "found the same variants again" - the
        one-hot surrogate arm measured +5.00 raw at t=+5.33 and +0.10 at d>=2 at t=+0.29.

Rule 2.1 is honoured by construction: genotypes go through
`veritas.calibration.genotypes.canonical_combo`, which fills wild type at every omitted
position. Nothing is zipped, no string is compared positionally, and no sixth genotype
implementation is written.

Usage: scripts/analysis/wt_anchor_readout.py <out_root> [landscape ...]

Not to be confused with `wt_anchor_confirmation_readout.py`, which collided with this
file add/add at the integration merge on 2026-11-01 and was renamed. That one reports the
two-arm ablation confirmation and resolves its runs by name from a campaign log; this one
reports all THREE anchor modes on both endpoints from a results root. Separate tools for
separate campaigns, each the provenance for its own results document.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

REPO = str(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from veritas.calibration.endpoint import greedy_diverse  # noqa: E402
from veritas.calibration.genotypes import canonical_combo  # noqa: E402
from veritas.calibration.landscapes import REGISTRY  # noqa: E402

MODES = ("zero", "ablated", "measured")
ARM = "panoptes_default"
CONTROLS = ("ridge", "random")
RAW_COL = "cum_top1_hits_excluding_initial"
BOOTSTRAP = 10_000
RNG = 20261101
EXPECTED_SEEDS = 48
EXPECTED_BUDGET = 125


# ---------------------------------------------------------------- loading + validity ----
def _candidates(root: str, landscape: str, mode: str) -> List[str]:
    return sorted(glob.glob(os.path.join(root, f"{landscape}__{mode}", "run_*")))


def admit(run_dir: str) -> Tuple[Optional[dict], Optional[list], Optional[list], str]:
    man_path = os.path.join(run_dir, "manifest.json")
    if not os.path.isfile(man_path):
        return None, None, None, (
            f"V1 FAIL: no manifest.json - PARTIAL LAUNCH, not data. checked {man_path}")
    man = json.load(open(man_path, encoding="utf-8"))
    state = (man.get("run") or {}).get("state")
    if state != "SUCCESS":
        return None, None, None, f"V2 FAIL: run.state={state!r}, not SUCCESS"
    ver = (man.get("run") or {}).get("wt_anchor_verified") or {}
    if not ver.get("agrees_with_configuration") or ver.get("vacuous"):
        return None, None, None, (
            f"V5 FAIL (rule 2.12): agrees={ver.get('agrees_with_configuration')} "
            f"vacuous={ver.get('vacuous')} observed={ver.get('modes_observed')}")
    rm = os.path.join(run_dir, "round_metrics.csv")
    st = os.path.join(run_dir, "selection_trace.csv")
    for p in (rm, st):
        if not os.path.isfile(p):
            return None, None, None, f"V1 FAIL: missing {p}"
    rows = list(csv.DictReader(open(rm, newline="", encoding="utf-8")))
    trace = list(csv.DictReader(open(st, newline="", encoding="utf-8")))
    return man, rows, trace, "ADMITTED"


def finals(rows) -> Dict[Tuple[str, int], dict]:
    best: Dict[Tuple[str, int], dict] = {}
    for r in rows:
        if r.get("status") != "completed":
            continue
        k = (r["method"], int(r["seed"]))
        if k not in best or int(r["round"]) > int(best[k]["round"]):
            best[k] = r
    return best


def budget_ok(fin, method) -> Tuple[bool, str]:
    seeds = sorted(s for (m, s) in fin if m == method)
    if len(seeds) != EXPECTED_SEEDS:
        return False, f"V3 FAIL: {method} completed {len(seeds)} seeds, not {EXPECTED_SEEDS}"
    bad = [(s, int(fin[(method, s)]["n_unique_measurements"])) for s in seeds
           if int(fin[(method, s)]["n_unique_measurements"]) != EXPECTED_BUDGET]
    if bad:
        return False, (f"V4 FAIL: {method} short of budget on {len(bad)} seed(s), "
                       f"e.g. {bad[:4]} (expected {EXPECTED_BUDGET})")
    return True, f"{method}: {len(seeds)} seeds, all at {EXPECTED_BUDGET} measurements"


def threshold_of(man) -> Optional[float]:
    t = (man.get("results") or {}).get("top1") or {}
    for k in ("threshold_value", "threshold"):
        if t.get(k) is not None:
            return float(t[k])
    return None


# ------------------------------------------------------------------------- endpoints ----
def hits_per_seed(trace, method: str, threshold: float, landscape: str,
                  include_initial: bool) -> Dict[int, List[Tuple[float, tuple]]]:
    """Measured hits per seed, as (fitness, canonical combo). Rule 2.1 via canonical_combo."""
    spec = REGISTRY[landscape]
    out: Dict[int, Dict[tuple, float]] = defaultdict(dict)
    for r in trace:
        if r.get("method") != method:
            continue
        status = (r.get("status") or "").strip()
        if status not in ("accepted", "initial"):
            continue
        if status == "initial" and not include_initial:
            continue
        try:
            f = float(r["revealed_fitness"])
            seed = int(float(r["seed"]))
        except (TypeError, ValueError):
            continue
        if f < threshold:
            continue
        combo = canonical_combo(str(r.get("variant_id", "")), spec.positions,
                                spec.parent_residues)
        if combo is None:
            continue
        prev = out[seed].get(combo)
        out[seed][combo] = f if prev is None else max(prev, f)
    return {s: sorted(((f, g) for g, f in m.items()), key=lambda p: -p[0])
            for s, m in out.items()}


# ------------------------------------------------------------------------ statistics ----
def paired(a: np.ndarray, b: np.ndarray) -> Dict[str, float]:
    """b - a, paired over seeds; seed is the unit of independence."""
    rng = np.random.default_rng(RNG)
    d = b - a
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1)) if n > 1 else 0.0
    t = mean / (sd / math.sqrt(n)) if sd > 0 else (math.inf if mean else 0.0)
    better, worse = int((d > 0).sum()), int((d < 0).sum())
    tied = int((d == 0).sum())
    try:
        t_p = float(stats.ttest_rel(b, a).pvalue)
    except Exception:
        t_p = float("nan")
    try:
        w_p = float(stats.wilcoxon(d).pvalue) if (d != 0).any() else 1.0
    except Exception:
        w_p = float("nan")
    nz = better + worse
    sign_p = float(stats.binomtest(better, nz, 0.5).pvalue) if nz else 1.0
    idx = rng.integers(0, n, size=(BOOTSTRAP, n))
    boots = d[idx].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    # Minimum detectable effect at 80% power, two-sided alpha 0.05, paired.
    mde = 2.8 * sd / math.sqrt(n) if sd > 0 else 0.0
    return {"n": n, "delta": mean, "sd": sd, "t": t, "t_p": t_p, "w_p": w_p,
            "sign_p": sign_p, "better": better, "worse": worse, "tied": tied,
            "lo": float(lo), "hi": float(hi), "mde": mde,
            "mean_a": float(a.mean()), "mean_b": float(b.mean())}


def line(label: str, r: Dict[str, float]) -> str:
    sep = "SEPARABLE" if not (r["lo"] <= 0 <= r["hi"]) else "not separable"
    return (f"    {label:30s} D={r['delta']:+8.3f}  t={r['t']:+7.2f}  p={r['t_p']:9.3g}  "
            f"B/W/T={r['better']:2d}/{r['worse']:2d}/{r['tied']:2d}  "
            f"CI95=[{r['lo']:+7.3f},{r['hi']:+7.3f}]  MDE={r['mde']:5.3f}  {sep}")


# ------------------------------------------------------------------------------ main ----
def analyse(root: str, landscape: str) -> None:
    print("=" * 108)
    print(f"LANDSCAPE  {landscape}")
    print("=" * 108)

    loaded: Dict[str, dict] = {}
    for mode in MODES:
        cands = _candidates(root, landscape, mode)
        if not cands:
            print(f"  {mode:9s} MISSING: no run_* under "
                  f"{os.path.join(root, f'{landscape}__{mode}')}")
            continue
        admitted = None
        for rd in cands:
            man, rows, trace, why = admit(rd)
            print(f"  {mode:9s} {os.path.basename(rd)}  {why}")
            if man is None:
                continue
            fin = finals(rows)
            all_ok = True
            for m in (ARM,) + CONTROLS:
                if not any(k[0] == m for k in fin):
                    continue
                ok, msg = budget_ok(fin, m)
                print(f"            {msg}")
                all_ok = all_ok and ok
            if not all_ok:
                print(f"            DISCARDED by construction")
                continue
            admitted = {"dir": rd, "man": man, "rows": rows, "trace": trace, "fin": fin}
        if admitted:
            loaded[mode] = admitted

    if len(loaded) < 2:
        print("  fewer than two admissible modes; no contrast possible")
        return

    # ---- pairing premises -----------------------------------------------------------
    print("\n  PAIRING PREMISES")
    ok = True
    sb = {m: ((d["man"].get("configuration") or {}).get("seeding") or {})
          .get("seed_block_sha256") for m, d in loaded.items()}
    if len({str(v) for v in sb.values()}) != 1:
        print(f"    P1 FAIL: seed blocks differ {sb}"); ok = False
    else:
        print(f"    P1 OK  seed_block_sha256 = {str(list(sb.values())[0])[:32]} "
              f"in all {len(sb)} folders")
    isets = {m: dict(((d["man"].get("configuration") or {}).get("seeding") or {})
                     .get("initial_sets_sha256") or {}) for m, d in loaded.items()}
    ref = sorted(isets)[0]
    mism = [(m, s) for m in isets if m != ref for s, h in isets[ref].items()
            if isets[m].get(s) != h]
    if mism:
        print(f"    P2 FAIL: starting sample differs in {len(mism)} (mode,seed), "
              f"e.g. {mism[:3]}"); ok = False
    else:
        print(f"    P2 OK  starting sample identical for all {len(isets[ref])} seeds")
    thr = {m: threshold_of(d["man"]) for m, d in loaded.items()}
    if len({round(v, 12) for v in thr.values() if v is not None}) != 1:
        print(f"    P3 FAIL: thresholds differ {thr}"); ok = False
    else:
        print(f"    P3 OK  top-1% threshold from manifests (rule 2.3) = "
              f"{list(thr.values())[0]!r}")
    THR = list(thr.values())[0]

    for c in CONTROLS:
        series = {}
        for m, d in loaded.items():
            series[m] = {s: float(d["fin"][(c, s)][RAW_COL])
                         for (mm, s) in d["fin"] if mm == c}
        present = [m for m in series if series[m]]
        if len(present) < 2:
            print(f"    P4 --  {c} present in {len(present)} folder(s); not checkable")
            continue
        r0 = present[0]
        bad = [(m, s) for m in present[1:] for s in sorted(set(series[r0]) & set(series[m]))
               if series[r0][s] != series[m][s]]
        if bad:
            print(f"    P4 FAIL: {c} differs across folders in {len(bad)} seed(s), "
                  f"e.g. {bad[:3]}"); ok = False
        else:
            mu = float(np.mean(list(series[r0].values())))
            print(f"    P4 OK  {c:7s} identical across all {len(present)} folders "
                  f"(mean {mu:.4f} raw, {len(series[r0])} seeds)")
    print(f"    => PREMISES {'HOLD' if ok else 'DO NOT HOLD - CONTRASTS ARE VOID'}")

    # ---- endpoints -------------------------------------------------------------------
    raw: Dict[str, Dict[int, float]] = {}
    div: Dict[str, Dict[int, float]] = {}
    rawchk: Dict[str, Dict[int, float]] = {}
    for m, d in loaded.items():
        raw[m] = {s: float(d["fin"][(ARM, s)][RAW_COL])
                  for (mm, s) in d["fin"] if mm == ARM}
        hps = hits_per_seed(d["trace"], ARM, THR, landscape, include_initial=False)
        # A seed that found NO hit is absent from `hps` and its endpoint is 0, not
        # missing. Dropping it would silently exclude the worst seeds from the mean and
        # unpair the arms - so the seed set comes from the completed runs, not the hits.
        div[m] = {s: float(greedy_diverse(hps.get(s, []), 2)) for s in raw[m]}
        rawchk[m] = {s: float(len(hps.get(s, []))) for s in raw[m]}
    # cross-check the trace-derived raw count against the harness's own column
    print("\n  ENDPOINT CROSS-CHECK (round_metrics column vs recomputed from the trace)")
    for m in MODES:
        if m not in raw:
            continue
        common = sorted(set(raw[m]) & set(rawchk[m]))
        dif = [s for s in common if abs(raw[m][s] - rawchk[m][s]) > 1e-9]
        print(f"    {m:9s} agree on {len(common) - len(dif)}/{len(common)} seeds"
              + (f"  MISMATCHES e.g. "
                 f"{[(s, raw[m][s], rawchk[m][s]) for s in dif[:3]]}" if dif else ""))

    seeds = sorted(set.intersection(*[set(raw[m]) for m in raw]))
    print(f"\n  {len(seeds)} paired seeds, arm {ARM}")
    hdr = f"    {'mode':9s} {'raw mean':>10s} {'raw sd':>8s} {'d>=2 mean':>10s} {'d>=2 sd':>8s}"
    print(hdr)
    for m in MODES:
        if m not in raw:
            continue
        a = np.array([raw[m][s] for s in seeds])
        b = np.array([div[m][s] for s in seeds])
        print(f"    {m:9s} {a.mean():10.3f} {a.std(ddof=1):8.3f} "
              f"{b.mean():10.3f} {b.std(ddof=1):8.3f}")
    for c in CONTROLS:
        for m, d in loaded.items():
            ss = sorted(s for (mm, s) in d["fin"] if mm == c)
            if not ss:
                continue
            a = np.array([float(d["fin"][(c, s)][RAW_COL]) for s in ss])
            hpsc = hits_per_seed(d["trace"], c, THR, landscape, include_initial=False)
            b = np.array([float(greedy_diverse(hpsc.get(s, []), 2)) for s in ss])
            print(f"    {c:9s} {a.mean():10.3f} {a.std(ddof=1):8.3f} "
                  f"{b.mean():10.3f} {b.std(ddof=1):8.3f}   (reference)")
            break

    pairs = [("zero", "measured"), ("ablated", "measured"), ("zero", "ablated")]
    for name, data in (("RAW top-1% hits excluding initial", raw),
                       ("DIVERSITY-ADJUSTED d>=2", div)):
        print(f"\n  PAIRED CONTRASTS - {name}")
        for a_m, b_m in pairs:
            if a_m not in data or b_m not in data:
                print(f"    {b_m} - {a_m}: UNAVAILABLE")
                continue
            a = np.array([data[a_m][s] for s in seeds], dtype=float)
            b = np.array([data[b_m][s] for s in seeds], dtype=float)
            print(line(f"{b_m} - {a_m}", paired(a, b)))
        # engine vs ridge per mode
        for m in MODES:
            if m not in loaded:
                continue
            d = loaded[m]
            sc = sorted({s for (mm, s) in d["fin"] if mm == ARM} &
                        {s for (mm, s) in d["fin"] if mm == "ridge"})
            if not sc:
                continue
            if data is raw:
                e = np.array([float(d["fin"][(ARM, s)][RAW_COL]) for s in sc])
                rg = np.array([float(d["fin"][("ridge", s)][RAW_COL]) for s in sc])
            else:
                he = hits_per_seed(d["trace"], ARM, THR, landscape, False)
                hr = hits_per_seed(d["trace"], "ridge", THR, landscape, False)
                e = np.array([float(greedy_diverse(he.get(s, []), 2)) for s in sc])
                rg = np.array([float(greedy_diverse(hr.get(s, []), 2)) for s in sc])
            print(line(f"{m}: engine - ridge", paired(rg, e)))
    print()


def main() -> int:
    root = sys.argv[1]
    for ls in (sys.argv[2:] or ["creilov", "trpb", "gb1", "cr9114_h1"]):
        analyse(root, ls)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
