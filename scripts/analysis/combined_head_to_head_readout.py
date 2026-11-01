"""
combined_head_to_head_readout.py - one run folder, nine arms, one seed block, BOTH endpoints.

Reads out the campaign pre-registered in docs/prereg-combined-head-to-head.md: the engine
with the one-hot GP surrogate (`panoptes_onehot_surrogate`) against a seven-arm baseline
panel and against `panoptes_default`, measured in ONE run folder on ONE shared starting
sample per seed.

WHY A SIBLING RATHER THAN A FLAG. `baseline_panel_readout.py` reports one endpoint and
`baseline_panel_diversity.py` reports the diversity-adjusted one over
`BASELINE_PANEL_METHODS`, which does not contain the ninth arm. The pre-registration makes
BOTH endpoints co-primary, because the panel already proved the answer depends on which one
is used: on raw TrpB `gp_onehot` beats `ridge` by +9.48 and at d>=2 it LOSES to ridge by
-0.94 (t=-2.0), redundancy 7.06x against 4.54x. So "the strongest baseline" is a different
arm under each endpoint, and a read-out that fixes one comparator across both endpoints is
choosing the comparator while appearing to measure it. That is the structure of the
withdrawn "6.3x the best baseline" claim, whose comparator was fixed at one budget and then
reused at every other.

Rule 2.11 - A VERDICT NEEDS A PAIRED TEST, NOT A MEAN DIFFERENCE. Nothing here prints a mean
difference without the paired t, seeds better/worse/tied, Wilcoxon p, sign-test p and a
whole-seed bootstrap interval beside it, and |t| < 2 with no seed majority is NOT SEPARABLE
with no sign attached.

Rule 2.1 - `variant_id` in selection_trace.csv is a MUTATION STRING that OMITS wild-type
positions. Every genotype is expanded to full library length with
`veritas.calibration.genotypes.canonical_combo` and compared with its `hamming`, which
raises on unequal lengths. Nothing is zipped and nothing is truncated.

Rule 2.3 - the top-1% bar comes from the run's OWN manifest.json, never recomputed.

`load_final_rounds` and `paired_bootstrap_ci` are IMPORTED from paired_stats.py, and
`check_pairing` from baseline_panel_readout.py, rather than reimplemented: a second recipe
for the same quantity is how two published verdicts here came to disagree with their own
artifacts.

Usage:
  python combined_head_to_head_readout.py <run_dir|run_name> <library_key>
      [--engine panoptes_onehot_surrogate] [--also-engine panoptes_default]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import stats

CHECKOUT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CHECKOUT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from veritas import paths as vpaths  # noqa: E402
from veritas.calibration.genotypes import canonical_combo, hamming  # noqa: E402
import benchmark_panoptes_isolated as bench  # noqa: E402
from paired_stats import load_final_rounds, paired_bootstrap_ci  # noqa: E402
from baseline_panel_readout import (  # noqa: E402
    ENGINE_ARMS,
    PANEL_ORDER,
    T_SEPARABLE,
    check_pairing,
    describe_run,
    load_with_initial,
)

#: Diversity separations reported beside the raw endpoint. d>=2 is co-primary with raw under
#: the pre-registration; d>=3 is context.
DISTANCES: Tuple[int, ...] = (2, 3)

RAW = "raw"

#: Development data under rule 2.4: a result on one of these is MECHANISM evidence and must
#: never be presented as generalization evidence, whatever its sign or size.
DEVELOPMENT_LIBRARIES = {"trpb", "gb1", "her2", "cr9114_h1"}


def paired_t(diff: np.ndarray) -> float:
    if diff.size < 2:
        return float("nan")
    sd = float(diff.std(ddof=1))
    if sd == 0.0:
        return float("nan")
    return float(diff.mean() / (sd / math.sqrt(diff.size)))


def verdict(diff: np.ndarray) -> str:
    if diff.size == 0:
        return "no shared seeds"
    if np.allclose(diff, 0.0):
        return "IDENTICAL"
    t = paired_t(diff)
    better = int((diff > 0).sum())
    worse = int((diff < 0).sum())
    if not math.isfinite(t) or abs(t) < T_SEPARABLE:
        return "not separable"
    return "WIN" if better > worse else "LOSS"


HEADER = (f"{'contrast':<40}{'n':>4}{'mean diff':>11}{'t':>8}{'b/w/t':>12}"
          f"{'95% bootstrap CI':>20}{'wilcox p':>11}{'sign p':>10}  verdict")


def row(label: str, diff: np.ndarray) -> str:
    if diff.size == 0:
        return f"{label:<40}{'-':>4}  no shared seeds"
    lo, hi = paired_bootstrap_ci(diff, iterations=10_000)
    better = int((diff > 0).sum())
    worse = int((diff < 0).sum())
    tied = int((diff == 0).sum())
    t = paired_t(diff)
    if np.allclose(diff, 0.0):
        wil = float("nan")
        sign = float("nan")
    else:
        wil = float(stats.wilcoxon(diff, zero_method="wilcox").pvalue)
        sign = (float(stats.binomtest(better, better + worse, 0.5).pvalue)
                if better + worse else float("nan"))
    return (f"{label:<40}{diff.size:>4}{diff.mean():>11.3f}{t:>8.2f}"
            f"{f'{better}/{worse}/{tied}':>12}"
            f"{f'[{lo:+.2f}, {hi:+.2f}]':>20}{wil:>11.3g}{sign:>10.3g}  {verdict(diff)}")


def resolve_run(raw: str) -> Path:
    """A run directory, by absolute path or by name against the REAL results root.

    Rule 2.8: `benchmark_comparison_results/` is gitignored and does not exist in a
    worktree, so a bare name resolves through `veritas.paths` against the real repository
    root. Every failure names every path it checked - a bare "not found" is what let our
    strongest baseline get silently dropped from a 48-seed run.
    """
    checked: List[Path] = []
    direct = Path(raw)
    checked.append(direct.resolve())
    if direct.is_dir():
        return direct.resolve()
    root = vpaths.results_dir("benchmark_comparison_results")
    candidate = root / raw
    checked.append(candidate)
    if candidate.is_dir():
        return candidate
    raise SystemExit(
        "no such run directory. Checked: " + "; ".join(str(p) for p in checked)
    )


def threshold_of(run_dir: Path) -> float:
    """Rule 2.3: the scoring bar comes from the run's OWN manifest.

    `results.top1.threshold_value`, NOT `configuration.components.qualify_threshold`:
    those are two different bars (the online 99.96th-percentile qualification bar against
    the top-1% scoring tier), and reading the wrong one silently rescores the whole panel
    against a tier of a few dozen variants.
    """
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    try:
        return float(manifest["results"]["top1"]["threshold_value"])
    except (KeyError, TypeError):
        raise SystemExit(f"no results.top1.threshold_value in {run_dir / 'manifest.json'}")


def seed_block(run_dir: Path) -> List[int]:
    """Every seed that appears in round_metrics.csv.

    NOT the seeds that appear in the hit table: an arm that found zero hits on a seed has
    no row there, so intersecting hit-seed sets drops every seed where the weakest arm
    scored nothing - which conditions the table on success and inflated `random` from 0.875
    to 1.60 once. A zero-hit seed is a measurement of zero, not a missing measurement.
    """
    seeds = set()
    with (run_dir / "round_metrics.csv").open(newline="", encoding="utf-8") as handle:
        for line in csv.DictReader(handle):
            seeds.add(int(float(line["seed"])))
    return sorted(seeds)


def adjusted_count(entries: Sequence[Tuple[float, str]], d: int) -> int:
    """Hits kept when each must be >= d Hamming from every hit already kept.

    Fitness-descending greedy, applied AFTER measurement so it cannot advantage a method
    that knew about it. `hamming` raises on unequal lengths, which is the rule-2.1 guard:
    both arguments are full-length genotypes with wild type filled.
    """
    kept: List[str] = []
    for _, combo in sorted(entries, key=lambda e: -e[0]):
        if all(hamming(combo, other) >= d for other in kept):
            kept.append(combo)
    return len(kept)


def hit_table(run_dir: Path, library_key: str, arms: Sequence[str], seeds: Sequence[int]):
    """arm -> seed -> [(fitness, full-length genotype)] for hits, initial sample EXCLUDED.

    The initial sample is excluded so the diversity endpoint matches
    `cum_top1_hits_excluding_initial` - the shared starting sample is identical across arms
    by construction and contributes no signal to a paired contrast.
    """
    library = bench.resolve_library(library_key)
    positions, parent = library.positions, library.parent_residues
    bar = threshold_of(run_dir)
    hits: Dict[str, Dict[int, List[Tuple[float, str]]]] = {
        arm: {seed: [] for seed in seeds} for arm in arms
    }
    unparsed = 0
    unknown_arm = defaultdict(int)
    with (run_dir / "selection_trace.csv").open(newline="", encoding="utf-8") as handle:
        for line in csv.DictReader(handle):
            if line.get("status") == "initial":
                continue
            raw_fitness = line.get("revealed_fitness")
            if raw_fitness in (None, ""):
                continue
            if float(raw_fitness) < bar:
                continue
            arm = line["method"]
            if arm not in hits:
                unknown_arm[arm] += 1
                continue
            combo = canonical_combo(line["variant_id"], positions, parent)
            if combo is None:
                unparsed += 1
                continue
            seed = int(line["seed"])
            if seed not in hits[arm]:
                continue
            hits[arm][seed].append((float(raw_fitness), combo))
    return hits, bar, unparsed, dict(unknown_arm), len(positions)


def ordered(arms: Sequence[str]) -> List[str]:
    known = [a for a in PANEL_ORDER if a in arms]
    return known + sorted(a for a in arms if a not in PANEL_ORDER)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument("library")
    parser.add_argument("--engine", default="panoptes_onehot_surrogate")
    parser.add_argument("--also-engine", default="panoptes_default")
    args = parser.parse_args()

    run_dir = resolve_run(args.run_dir)
    if not (run_dir / "manifest.json").is_file():
        raise SystemExit(
            f"no manifest.json in {run_dir}. A run directory without a manifest is a "
            f"PARTIAL LAUNCH, not data (rule 2.3)."
        )
    info = describe_run(run_dir)
    metrics = load_with_initial(run_dir)
    arms = ordered(list(metrics))
    seeds = seed_block(run_dir)

    print("=" * 118)
    print(f"COMBINED HEAD-TO-HEAD   {info['library']}   run {run_dir.name}")
    print(f"  budget {info['budget']} = {info['initial']} initial + {info['rounds']} x "
          f"{info['batch']}   seeds {info['seeds']} from base {info['base_seed']} "
          f"(block {min(seeds)}-{max(seeds)})   state {info['state']}"
          + (f" ({info['state_reason']})" if info["state_reason"] else ""))
    print(f"  arms {len(arms)}: {', '.join(arms)}")
    print(f"  configuration fingerprint {info['fingerprint']}  (quote the COMMIT, not the "
          f"fingerprint: rule 2.12)")
    if args.library in DEVELOPMENT_LIBRARIES:
        print(f"  !! {args.library} is DEVELOPMENT DATA under rule 2.4. Everything below is "
              f"MECHANISM evidence, never generalization evidence.")
    print("=" * 118)

    problems = check_pairing(metrics)
    if problems:
        print("\n  PAIRING / BUDGET PROBLEMS - read these before any number below:")
        for p in problems:
            print(f"    - {p}")
    else:
        print("\n  pairing premises hold: identical round-0 top-1% count per seed across all "
              "arms, budgets matched per arm.")

    hits, bar, unparsed, unknown, n_sites = hit_table(run_dir, args.library, arms, seeds)
    print(f"\n  top-1% bar {bar!r} from the run's own manifest; {n_sites}-site library, "
          f"genotypes expanded to full length with wild type filled")
    if unparsed:
        print(f"  !! {unparsed} variant_id(s) could not be canonicalised and were EXCLUDED")
    if unknown:
        print(f"  !! selection_trace carries arm(s) absent from round_metrics: {unknown}")

    # ---- the endpoints, per arm ---------------------------------------------------------
    values: Dict[str, Dict[str, np.ndarray]] = {}
    for arm in arms:
        per_arm = {RAW: np.array([len(hits[arm][s]) for s in seeds], dtype=float)}
        for d in DISTANCES:
            per_arm[f"d>={d}"] = np.array(
                [adjusted_count(hits[arm][s], d) for s in seeds], dtype=float)
        values[arm] = per_arm

    # Integrity: the raw count rebuilt from selection_trace must equal round_metrics'
    # cum_top1_hits_excluding_initial. They are two independent paths to the same quantity
    # and a disagreement means one of them is not measuring what its name says.
    print()
    print(f"  {'arm':<28}{'raw':>8}{'raw(rm)':>9}", end="")
    for d in DISTANCES:
        print(f"{'d>=' + str(d):>9}{'redund':>9}", end="")
    print(f"{'uniq':>7}{'bestfit':>10}  role")
    print("  " + "-" * 108)
    mismatches: List[str] = []
    for arm in arms:
        rm = np.array([metrics[arm][s]["cum_top1_hits_excluding_initial"]
                       for s in seeds if s in metrics[arm]], dtype=float)
        uniq = np.array([metrics[arm][s]["n_unique_measurements"]
                         for s in seeds if s in metrics[arm]], dtype=float)
        best = np.array([metrics[arm][s]["best_measured_fitness"]
                         for s in seeds if s in metrics[arm]], dtype=float)
        raw = values[arm][RAW]
        if rm.size == raw.size and not np.allclose(rm, raw):
            mismatches.append(
                f"{arm}: selection_trace raw mean {raw.mean():.3f} != round_metrics "
                f"{rm.mean():.3f} on {int((rm != raw).sum())} seed(s)")
        role = ("ENGINE (under test)" if arm in ENGINE_ARMS else "baseline")
        line = f"  {arm:<28}{raw.mean():>8.2f}{rm.mean():>9.2f}"
        for d in DISTANCES:
            adj = values[arm][f'd>={d}']
            redund = raw.mean() / adj.mean() if adj.mean() else float("nan")
            line += f"{adj.mean():>9.2f}{redund:>8.2f}x"
        line += f"{uniq.mean():>7.1f}{best.mean():>10.4f}  {role}"
        print(line)
    if mismatches:
        print("\n  !! ENDPOINT DISAGREEMENT between selection_trace.csv and round_metrics.csv:")
        for m in mismatches:
            print(f"     - {m}")
    else:
        print("\n  raw endpoint agrees exactly between selection_trace.csv and "
              "round_metrics.csv for every arm.")

    baselines = [a for a in arms if a not in ENGINE_ARMS]
    engines = [a for a in (args.engine, args.also_engine) if a in values]
    summary: List[dict] = []

    for endpoint in (RAW,) + tuple(f"d>={d}" for d in DISTANCES):
        print()
        print("=" * 118)
        label = ("raw cum_top1_hits_excluding_initial" if endpoint == RAW
                 else f"diversity-adjusted, {endpoint} Hamming separation")
        print(f"ENDPOINT: {label}")
        print("=" * 118)
        if not baselines:
            print("  no baseline arms in this run")
            continue
        means = {a: float(values[a][endpoint].mean()) for a in arms}
        strongest = max(baselines, key=lambda a: means[a])
        print(f"  strongest baseline under THIS endpoint: `{strongest}` at "
              f"{means[strongest]:.3f}")
        ranked = ", ".join(f"{a} {means[a]:.2f}" for a in
                           sorted(baselines, key=lambda a: -means[a]))
        print(f"  baselines ranked: {ranked}")

        for engine in engines:
            print(f"\n  A. `{engine}` against every baseline   ({engine} - arm)")
            print("  " + HEADER)
            print("  " + "-" * 114)
            for arm in baselines:
                diff = values[engine][endpoint] - values[arm][endpoint]
                print("  " + row(arm, diff))
            diff = values[engine][endpoint] - values[strongest][endpoint]
            print(f"\n  B. HEADLINE: `{engine}` against the STRONGEST baseline "
                  f"`{strongest}` under this endpoint")
            print("  " + HEADER)
            print("  " + row(f"{engine} - {strongest}", diff))
            summary.append({
                "endpoint": endpoint, "engine": engine, "strongest": strongest,
                "baseline_mean": means[strongest], "engine_mean": means[engine],
                "delta": float(diff.mean()), "t": paired_t(diff),
                "better": int((diff > 0).sum()), "worse": int((diff < 0).sum()),
                "tied": int((diff == 0).sum()),
                "ci": paired_bootstrap_ci(diff, iterations=10_000),
                "verdict": verdict(diff),
            })
            print(f"     NOTE: a maximum over {len(baselines)} baseline arms is upward-biased "
                  f"for the baseline by the winner's curse - the bias favours the BASELINE, "
                  f"not the engine. Pre-declared in docs/prereg-combined-head-to-head.md "
                  f"section 2; section A above carries every pairwise row.")

        if len(engines) == 2:
            diff = values[engines[0]][endpoint] - values[engines[1]][endpoint]
            print(f"\n  C. The two engine arms against each other "
                  f"({engines[0]} - {engines[1]})")
            print("  " + HEADER)
            print("  " + row(f"{engines[0]} - {engines[1]}", diff))

        if "random" in values and "ridge" in values:
            diff = values["ridge"][endpoint] - values["random"][endpoint]
            call = verdict(diff)
            print(f"\n  D. Is the registered reference separable from `random`?  "
                  f"ridge - random = {diff.mean():+.3f}, t={paired_t(diff):+.2f}, "
                  f"{int((diff > 0).sum())}/{int((diff < 0).sum())}/"
                  f"{int((diff == 0).sum())} -> {call}")
            model_arms = [a for a in baselines if a != "random"]
            separable = [a for a in model_arms
                         if verdict(values[a][endpoint] - values["random"][endpoint])
                         != "not separable"]
            if not separable:
                print("     => NO model-based baseline is separable from `random` under this "
                      "endpoint. The endpoint cannot discriminate methods here, so EVERY "
                      "verdict on it is weak regardless of sign. That is a property of the "
                      "library and no budget change fixes it.")

    print()
    print("=" * 118)
    print("SUMMARY - engine against the strongest baseline, PER ENDPOINT")
    print("=" * 118)
    print(f"{'endpoint':<10}{'engine arm':<30}{'strongest baseline':<20}{'base':>8}"
          f"{'engine':>8}{'delta':>9}{'t':>8}{'b/w/t':>11}{'95% CI':>20}  verdict")
    print("-" * 118)
    for s in summary:
        bwt = f"{s['better']}/{s['worse']}/{s['tied']}"
        ci = f"[{s['ci'][0]:+.2f}, {s['ci'][1]:+.2f}]"
        print(f"{s['endpoint']:<10}{s['engine']:<30}{s['strongest']:<20}"
              f"{s['baseline_mean']:>8.2f}{s['engine_mean']:>8.2f}{s['delta']:>9.2f}"
              f"{s['t']:>8.2f}{bwt:>11}{ci:>20}  {s['verdict']}")
    print("-" * 118)
    print("delta = engine - strongest baseline UNDER THAT ENDPOINT. The strongest baseline is")
    print("re-selected per endpoint on purpose: the panel measured `gp_onehot` beating `ridge`")
    print("on raw TrpB and LOSING to it at d>=2, so fixing one comparator across both would be")
    print("choosing the comparator while appearing to measure it.")
    print("Hamming separation is a sequence-level proxy for distinctness, NOT functional")
    print("independence.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
