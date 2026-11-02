"""Why does the engine lose when its candidate ballot already contains the hits?

`scripts/analysis/pool_hit_reachability.py` answers the prior question: on eight of the
nine landscapes the pre-veto ballot reaches a depth at which 100% of the top-1% hit set
lives (CreiLOV 99.9%), so the candidate generator's DEPTH is not the constraint
anywhere. On CreiLOV it also shows the post-veto ballot holds enough hits for a perfect
ranker to match one-hot ridge. That localises the loss to RANKING, and this script
measures the ranking, in three independent ways, from artifacts a completed run already
contains plus the label column.

WHAT IT REPORTS
---------------
1. STAGE-WISE RANK CORRELATION, per campaign, from `candidate_trace.csv`:

     resid  the GP residual alone - predicted_fitness minus prior_contribution,
            because mu(x) = m(x) + mu_res(x) (math_worker.py, "Reconstruct full
            predicted fitness")
     pred   mu(x), the composite the optimizer carries forward
     acq    EI, or EI + the Thompson term under hybrid acquisition
     qual   acquisition x Pareto x auxiliary weights, pre-diversity
     prior  the biophysical / structural prior m(x) on its own

   Each against TRUE fitness. The point of separating them is that they compose
   additively on one scale, so a prior that is anti-correlated can cancel and then
   reverse a surrogate that is weakly correct - and the composite is what gets ranked.

2. THE DEPTH-MATCHED BASELINE. An arm's selected-depth histogram (from
   `generation_coverage.json`) is scored against the landscape's own hit rate at each
   depth. That is what a picker with the same depth preference and no other information
   would have achieved. Dividing the arm's real precision by it gives the part of the
   result that is NOT explained by depth preference - which is the number that settles
   whether a depth argument can explain a loss at all. The `random` arm must come out
   near 1.0x; if it does not, the construction is wrong and nothing else here holds.

3. THE SYNTHETIC WILD-TYPE LABEL. `select_next_micro_batch` injects one training row per
   round for the wild type with the literal label `y = 0.0`, whatever the assay's scale
   is, at variance 1e-4 against a real measurement's 1e-2 default. This prints, per
   landscape, how far that label is from the centre of the real labels in label standard
   deviations. It is a property of the label column alone, available before any campaign,
   so it is not a post-hoc fit.

HOW IT STAYS HONEST
-------------------
* Depths, the eligible universe, the hit set and the threshold all come from
  `veritas.calibration.landscapes` - the landscape's own arithmetic. Candidate notations
  are mapped with `veritas.calibration.genotypes.canonical_combo`, which fills wild type
  at omitted positions and asserts full length (HANDOFF rule 2.1). No second parser.
* HANDOFF rule 2.2's anchors are reproduced and printed before any result, and the
  report is stamped unverified if either misses.
* `candidate_trace.csv` is BOUNDED - the top rows by quality plus a sample of rejected
  ones - so every rho here is a range-restricted estimate over a quality-biased
  subsample, and the per-stage base rate it implies is inflated. Both biases flatter the
  engine, so the stage-wise picture is a lower bound on the problem, not an upper one.
  This is printed with the table rather than kept in a docstring.
* Nothing is recomputed that an artifact records: the threshold is the landscape's, the
  yields are read from `round_metrics.csv`, and the depth histograms are the run's own.

USAGE
-----
    .venv\\Scripts\\python.exe -X utf8 scripts/analysis/ranking_decomposition.py \\
        --run <run_dir> --library creilov

    ... --wt-scale-only     # just table 3, no run directory needed

Read-only. It runs no campaign and writes nothing unless --json is given.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from veritas.calibration.genotypes import canonical_combo  # noqa: E402
from veritas.calibration.landscapes import REGISTRY, load_landscape  # noqa: E402

SHALLOW_SHELL_ANCHORS = {"trpb": (168, 2_242), "gb1": (72, 2_167)}
STAGES = ("resid", "pred", "acq", "qual", "prior")
RULE = "=" * 100


def check_anchors() -> bool:
    print()
    print("GENOTYPE VALIDATION (HANDOFF rule 2.2 - shallow shell 1 <= depth <= 2)")
    ok = True
    for key, (want_hits, want_n) in SHALLOW_SHELL_ANCHORS.items():
        try:
            shell = load_landscape(key).shell(2)
        except Exception as exc:
            print(f"  {key:<6} COULD NOT CHECK: {type(exc).__name__}: {exc}")
            ok = False
            continue
        got = (shell.n_hits, shell.n_variants)
        good = got == (want_hits, want_n)
        ok = ok and good
        print(f"  {key:<6} {shell.n_hits:>5,} / {shell.n_variants:>6,} = "
              f"{100.0 * shell.rate:7.4f}%   "
              f"{'PASS' if good else f'FAIL (expected {want_hits}/{want_n})'}")
    if not ok:
        print("  *** anchors did not reproduce; nothing below is trustworthy ***")
    return ok


# ==============================================================================
# 1. Stage-wise rank correlation
# ==============================================================================
def stage_correlations(run: Path, landscape) -> Dict[str, Any]:
    spec = landscape.spec
    index = {c: i for i, c in enumerate(landscape.combos)}
    per_campaign: Dict[str, List[tuple]] = defaultdict(list)
    n_rows = n_off = 0
    with (run / "candidate_trace.csv").open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            n_rows += 1
            combo = canonical_combo(row["variant_id"], spec.positions, spec.parent_residues)
            i = index.get(combo) if combo is not None else None
            if i is None:
                n_off += 1
                continue
            per_campaign[row["campaign"]].append((
                float(landscape.fitness[i]), bool(landscape.hit_mask[i]),
                float(row["predicted_fitness"]), float(row["acquisition_score"]),
                float(row["quality_score"]), int(row["selection_order"]),
                float(row["prior_contribution"]),
            ))

    def rho(a, f):
        return stats.spearmanr(a, f).correlation if np.std(a) > 0 else np.nan

    acc: Dict[str, List[float]] = {s: [] for s in STAGES}
    sel_hits = sel_n = pool_hits = pool_n = 0
    for v in per_campaign.values():
        f = np.array([x[0] for x in v]); h = np.array([x[1] for x in v], dtype=float)
        pred = np.array([x[2] for x in v]); acq = np.array([x[3] for x in v])
        qual = np.array([x[4] for x in v]); sel = np.array([x[5] for x in v])
        prior = np.array([x[6] for x in v])
        acc["resid"].append(rho(pred - prior, f))
        acc["pred"].append(rho(pred, f))
        acc["acq"].append(rho(acq, f))
        acc["qual"].append(rho(qual, f))
        acc["prior"].append(rho(prior, f))
        m = sel > 0
        sel_hits += int(h[m].sum()); sel_n += int(m.sum())
        pool_hits += int(h.sum()); pool_n += len(v)

    out: Dict[str, Any] = {"n_trace_rows": n_rows, "n_off_universe": n_off,
                           "n_campaigns": len(per_campaign), "stages": {}}
    for s in STAGES:
        a = np.array(acc[s], dtype=float)
        a = a[~np.isnan(a)]
        if not a.size:
            continue
        out["stages"][s] = {
            "mean": round(float(a.mean()), 4), "sd": round(float(a.std(ddof=1)), 4)
            if a.size > 1 else None,
            "min": round(float(a.min()), 4), "median": round(float(np.median(a)), 4),
            "max": round(float(a.max()), 4), "n_negative": int((a < 0).sum()),
            "n": int(a.size),
        }
    out["selected_precision"] = (sel_hits / sel_n) if sel_n else None
    out["traced_ballot_base_rate"] = (pool_hits / pool_n) if pool_n else None
    out["lift_over_own_ballot"] = (
        (sel_hits / sel_n) / (pool_hits / pool_n) if sel_n and pool_hits else None
    )
    return out


def print_stage_table(r: Dict[str, Any]) -> None:
    print()
    print("1. STAGE-WISE RANK CORRELATION AGAINST TRUE FITNESS")
    print(f"   {r['n_trace_rows']:,} traced rows ({r['n_off_universe']:,} not eligible-"
          f"universe members), {r['n_campaigns']} campaign(s)")
    print()
    print(f"   {'stage':<8}{'mean rho':>10}{'sd':>8}{'min':>8}{'median':>8}{'max':>8}"
          f"{'rho<0':>8}")
    print("   " + "-" * 58)
    for s in STAGES:
        d = r["stages"].get(s)
        if d is None:
            continue
        sd = "-" if d["sd"] is None else f"{d['sd']:.3f}"
        neg = f"{d['n_negative']}/{d['n']}"
        print(f"   {s:<8}{d['mean']:>10.3f}{sd:>8}{d['min']:>8.3f}"
              f"{d['median']:>8.3f}{d['max']:>8.3f}{neg:>8}")
    if r["lift_over_own_ballot"] is not None:
        print()
        print(f"   committed-batch precision : {100 * r['selected_precision']:.2f}%")
        print(f"   traced-ballot base rate   : {100 * r['traced_ballot_base_rate']:.2f}%")
        print(f"   LIFT OVER ITS OWN BALLOT  : {r['lift_over_own_ballot']:.2f}x"
              + ("   <-- below 1.0 means the full ranking pipeline is WORSE than picking "
                 "uniformly from the ballot it built"
                 if r["lift_over_own_ballot"] < 1.0 else ""))
    print()
    print("   candidate_trace.csv is BOUNDED (top rows by quality plus a sample of")
    print("   rejected ones), so each rho is range-restricted and the base rate is")
    print("   inflated. Both biases flatter the engine.")


# ==============================================================================
# 2. Depth-matched baseline
# ==============================================================================
def depth_matched(run: Path, landscape) -> Dict[str, Any]:
    cov = json.loads((run / "generation_coverage.json").read_text(encoding="utf-8"))
    hist = landscape.depth_histogram()  # {depth: (hits, variants)}
    yields = final_yields(run)
    out: Dict[str, Any] = {"arms": {}, "hit_rate_by_depth": {
        int(d): round(100.0 * h / n, 4) for d, (h, n) in hist.items()}}
    for arm in cov.get("arms", []):
        method = arm["method"]
        sel = arm["coverage"]["union_selected_depth_histogram"]
        total = sum(int(v) for v in sel.values())
        if not total:
            continue
        mean_depth = sum(int(d) * int(c) for d, c in sel.items()) / total
        expected = sum(
            (hist[int(d)][0] / hist[int(d)][1]) * int(c)
            for d, c in sel.items() if int(d) in hist
        ) / total
        got = yields.get(method)
        # Precision over the POST-INITIAL budget, which is what the arm's policy chose.
        precision = None if got is None else got["hits"] / got["post_initial"]
        out["arms"][method] = {
            "mean_selected_depth": round(mean_depth, 3),
            "depth_matched_hit_rate": round(100.0 * expected, 3),
            "observed_precision": None if precision is None else round(100.0 * precision, 3),
            "lift_over_depth_matched": (
                None if precision is None or expected <= 0 else round(precision / expected, 2)
            ),
            "selected_depth_histogram": {int(d): int(c) for d, c in sorted(
                sel.items(), key=lambda kv: int(kv[0]))},
        }
    return out


def final_yields(run: Path) -> Dict[str, Dict[str, float]]:
    """method -> mean final hits (excluding the initial sample) and the post-initial budget."""
    best: Dict[str, Dict[int, tuple]] = defaultdict(dict)
    with (run / "round_metrics.csv").open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            method, seed, rnd = row["method"], int(row["seed"]), int(row["round"])
            prev = best[method].get(seed)
            if prev is None or rnd > prev[0]:
                best[method][seed] = (rnd, float(row["cum_top1_hits_excluding_initial"]),
                                      float(row["n_unique_measurements"]))
    out = {}
    for method, per_seed in best.items():
        hits = float(np.mean([v[1] for v in per_seed.values()]))
        # The initial sample is shared and is not this policy's choice, so precision is
        # taken over the measurements the policy actually selected.
        total_unique = float(np.mean([v[2] for v in per_seed.values()]))
        out[method] = {"hits": hits, "post_initial": _post_initial(run, total_unique)}
    return out


def _post_initial(run: Path, total_unique: float) -> float:
    """Budget minus the shared initial sample, read from the manifest."""
    m = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    seeding = (m.get("configuration") or {}).get("seeding") or {}
    sizes = seeding.get("initial_set_sizes") or {}
    if isinstance(sizes, dict) and sizes:
        initial = float(np.mean([float(v) for v in sizes.values()]))
    else:
        initial = 0.0
    return max(1.0, total_unique - initial)


def print_depth_table(r: Dict[str, Any]) -> None:
    print()
    print("2. DEPTH-MATCHED BASELINE - can a depth argument explain the result?")
    print()
    print(f"   {'arm':<20}{'mean sel. depth':>16}{'depth-matched':>15}"
          f"{'observed':>10}{'lift':>8}")
    print("   " + "-" * 69)
    for method, d in sorted(r["arms"].items()):
        obs = "-" if d["observed_precision"] is None else f"{d['observed_precision']:.2f}%"
        lift = "-" if d["lift_over_depth_matched"] is None else f"{d['lift_over_depth_matched']:.2f}x"
        print(f"   {method:<20}{d['mean_selected_depth']:>16.3f}"
              f"{d['depth_matched_hit_rate']:>14.2f}%{obs:>10}{lift:>8}")
    rnd = r["arms"].get("random", {}).get("lift_over_depth_matched")
    print()
    print("   'depth-matched' is the landscape's own hit rate at the depths the arm chose,")
    print("   weighted by how often it chose them: what a picker with that depth")
    print("   preference and no other information would have scored. 'lift' is the")
    print("   observed precision over it - the part NOT explained by depth preference.")
    if rnd is not None:
        print(f"   CONSTRUCTION CHECK: the random arm's lift is {rnd:.2f}x and must be near")
        print("   1.0, since random has no information beyond its depth profile."
              + ("" if 0.8 <= rnd <= 1.25 else "  *** IT IS NOT - the construction is "
                 "wrong and this table must not be used. ***"))
    else:
        print("   CONSTRUCTION CHECK: no random arm in this run, so the baseline is")
        print("   UNVALIDATED here. Re-run with random included before quoting it.")


# ==============================================================================
# 3. The synthetic wild-type label
# ==============================================================================
def wt_scale_table() -> List[Dict[str, Any]]:
    rows = []
    for key in REGISTRY:
        try:
            ls = load_landscape(key)
        except Exception as exc:
            rows.append({"library": key, "error": f"{type(exc).__name__}: {exc}"})
            continue
        f = ls.fitness
        sd = float(f.std(ddof=1))
        rows.append({
            "library": key, "min": float(f.min()), "mean": float(f.mean()),
            "max": float(f.max()), "sd": sd,
            "zero_inside_range": bool(f.min() <= 0.0 <= f.max()),
            "sds_from_zero": (abs(float(f.mean())) / sd) if sd > 0 else None,
        })
    return rows


def print_wt_table(rows: List[Dict[str, Any]]) -> None:
    print()
    print("3. THE SYNTHETIC WILD-TYPE TRAINING LABEL (y = 0.0, variance 1e-4)")
    print("   Injected into every surrogate fit by select_next_micro_batch, with the")
    print("   assay's scale never consulted.")
    print()
    print(f"   {'library':<14}{'label min':>11}{'mean':>9}{'max':>9}{'sd':>8}"
          f"{'SDs from 0':>12}{'0 in range':>12}")
    print("   " + "-" * 76)
    for r in sorted(rows, key=lambda d: -(d.get("sds_from_zero") or 0.0)):
        if "error" in r:
            print(f"   {r['library']:<14}{r['error']}")
            continue
        z = "-" if r["sds_from_zero"] is None else f"{r['sds_from_zero']:.2f}"
        print(f"   {r['library']:<14}{r['min']:>11.3f}{r['mean']:>9.3f}{r['max']:>9.3f}"
              f"{r['sd']:>8.3f}{z:>12}{('yes' if r['zero_inside_range'] else 'NO'):>12}")
    print()
    print("   Sorted worst first. This is read off the label column alone and is fixed")
    print("   before any campaign runs, so it is not a post-hoc fit. It is NOT a")
    print("   predictor of the win/loss pattern either - CR9114 is second worst and the")
    print("   engine wins there - so treat it as a mechanism to ablate, not a")
    print("   correlation to report. The arm is --ablate-components wt_anchor.")


# ==============================================================================
def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("USAGE")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--run", type=Path, default=None,
                        help="a completed run directory (needs candidate_trace.csv, "
                             "generation_coverage.json, round_metrics.csv, manifest.json)")
    parser.add_argument("--library", default=None, metavar="KEY",
                        help="the landscape that run used.")
    parser.add_argument("--wt-scale-only", action="store_true",
                        help="print only table 3, which needs no run directory.")
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args(argv)

    print(RULE)
    print("RANKING DECOMPOSITION")
    print(RULE)
    anchors_ok = check_anchors()

    payload: Dict[str, Any] = {"anchors_ok": anchors_ok}
    if not args.wt_scale_only:
        if args.run is None or args.library is None:
            parser.error("--run and --library are both required unless --wt-scale-only")
        if args.library not in REGISTRY:
            parser.error(f"unknown landscape {args.library!r}; known: {sorted(REGISTRY)}")
        missing = [n for n in ("candidate_trace.csv", "generation_coverage.json",
                               "round_metrics.csv", "manifest.json")
                   if not (args.run / n).is_file()]
        if missing:
            print()
            print(f"CANNOT PROCEED: {args.run} is missing {missing}.")
            print(f"  checked: {[str(args.run / n) for n in missing]}")
            return 1
        landscape = load_landscape(args.library)
        print()
        print(f"run       : {args.run.name}")
        print(f"landscape : {args.library}  universe {landscape.n_universe:,}  "
              f"hits {landscape.n_hits:,}  threshold {landscape.threshold:.6g} "
              f"[{landscape.threshold_source}]")
        payload["stage_correlations"] = stage_correlations(args.run, landscape)
        print_stage_table(payload["stage_correlations"])
        payload["depth_matched"] = depth_matched(args.run, landscape)
        print_depth_table(payload["depth_matched"])

    payload["wt_scale"] = wt_scale_table()
    print_wt_table(payload["wt_scale"])

    if not anchors_ok:
        print()
        print("  *** rule-2.2 anchors did NOT reproduce; treat every table above as "
              "unverified. ***")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
