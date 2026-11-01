"""
baseline_panel_readout.py - the baseline panel, read out against the pre-registration.

Answers the two questions docs/prereg-baseline-panel.md asks:

  1. Is `ridge` a fair reference baseline, and where is it weak? Every panel arm is paired
     against `ridge` over seeds, so the question is answered per landscape rather than in
     general - ridge found TrpB's global optimum and was not separable from random on
     CR9114, so "is ridge fair" has no landscape-free answer.

  2. Against the STRONGEST baseline in the panel, how many of the engine's results survive?
     Every arm is also paired against `panoptes_default`.

It does NOT promote anything. `ridge` is the registered reference baseline and is printed
as such; docs/HANDOFF.md 5.4 makes promoting another model to reference status an amendment
to a frozen hypothesis, which needs a dated amendment record and not a read-out script.

Rule 2.11 - A VERDICT NEEDS A PAIRED TEST, NOT A MEAN DIFFERENCE. This is the most repeated
error in the project: eqFP611-blue was published as a WIN at +0.73 with t=0.77, better on
25 of 48 seeds, which is a coin flip. So nothing here prints a mean difference without the
paired t, the seeds better/worse/tied, the Wilcoxon p and a whole-seed bootstrap interval
beside it, and |t| < 2 with roughly half the seeds better is labelled NOT SEPARABLE rather
than given a sign.

`load_final_rounds` and `paired_bootstrap_ci` are IMPORTED from paired_stats.py rather than
reimplemented. A second recipe for the same quantity is how two published verdicts here
came to disagree with their own artifacts.

Usage:
  python baseline_panel_readout.py <run_dir> [<run_dir> ...]
      [--metric cum_top1_hits_excluding_initial] [--reference ridge]
      [--engine panoptes_default]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

# scripts/analysis is deliberately kept off sys.path by the isolated-mode convention this
# directory's README describes, so the sibling import is made explicit.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from paired_stats import load_final_rounds, paired_bootstrap_ci  # noqa: E402

#: Panel order, matching BASELINE_PANEL_METHODS in benchmark_panoptes_isolated.py. Arms
#: present in a run but absent here are appended, so an unexpected arm is never hidden.
PANEL_ORDER = (
    "random",
    "ridge",
    "random_forest",
    "gbm_onehot",
    "knn_hamming",
    "gp_onehot",
    "greedy_measured",
    "panoptes_default",
)

#: Arms that are baselines. `panoptes_default` is the thing under test and is excluded from
#: "the strongest baseline", which is the whole point of the comparison.
ENGINE_ARMS = ("panoptes_default", "panoptes_muse", "panoptes_qualify_then_optimize",
               "panoptes_conditional_probe")

#: |t| below this, with no clear seed majority, is reported as NOT SEPARABLE. Two-sided,
#: 47 degrees of freedom at 48 paired seeds; 2.0 is the convention already used in
#: scripts/analysis/heldout.py (T_CRIT = 2.01) and is kept identical here.
T_SEPARABLE = 2.0

#: Development data under rule 2.4. A result on one of these is mechanism evidence and must
#: never be presented as generalization evidence, whatever its sign or size.
DEVELOPMENT_LIBRARIES = {"trpb", "gb1", "her2", "cr9114_h1"}


def paired_t(diff: np.ndarray) -> float:
    """Paired t over seeds. NaN when the differences are identically zero."""
    if diff.size < 2:
        return float("nan")
    sd = float(diff.std(ddof=1))
    if sd == 0.0:
        return float("nan")
    return float(diff.mean() / (sd / math.sqrt(diff.size)))


def verdict(diff: np.ndarray) -> str:
    """The call, with the paired test built into it rather than inferred from the mean."""
    if diff.size == 0:
        return "no shared seeds"
    if np.allclose(diff, 0.0):
        return "IDENTICAL"
    t = paired_t(diff)
    better = int((diff > 0).sum())
    worse = int((diff < 0).sum())
    if not math.isfinite(t) or abs(t) < T_SEPARABLE:
        return "not separable"
    if better > worse:
        return "WIN"
    return "LOSS"


def manifest_of(run_dir: Path) -> dict:
    path = run_dir / "manifest.json"
    if not path.is_file():
        raise SystemExit(
            f"no manifest at {path}. A run directory without a manifest is not a result "
            f"(rule 2.3: thresholds come from the run's OWN manifest)."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def describe_run(run_dir: Path) -> dict:
    """Library, budget, seed block, state and threshold, read from the run's own manifest."""
    manifest = manifest_of(run_dir)
    resolved = (manifest.get("invocation") or {}).get("resolved_args") or {}
    run = manifest.get("run") or {}
    results = manifest.get("results") or {}
    threshold = None
    for key in ("top1", "top_1_percent", "endpoint"):
        block = results.get(key)
        if isinstance(block, dict) and "threshold" in block:
            threshold = block["threshold"]
            break
    return {
        "run_dir": run_dir,
        "library": resolved.get("library"),
        "initial": resolved.get("initial_size"),
        "batch": resolved.get("batch_size"),
        "rounds": resolved.get("rounds"),
        "seeds": resolved.get("seeds"),
        "base_seed": resolved.get("base_seed"),
        "budget": (
            (resolved.get("initial_size") or 0) + (resolved.get("rounds") or 0) * (resolved.get("batch_size") or 0)
        ),
        "state": run.get("state"),
        "state_reason": run.get("state_reason"),
        "fingerprint": manifest.get("configuration_fingerprint"),
        "threshold": threshold,
    }


def check_pairing(table: Dict[str, Dict[int, dict]]) -> List[str]:
    """
    The paired test's premises, CHECKED rather than assumed.

    Two things have to hold or every difference below is meaningless:

      1. The initial sample must be identical per seed across arms. Round 0's
         `cum_top1_hits_including_initial` is a function of the starting variants alone, so
         any disagreement across arms at the same seed means the arms did not start from the
         same sample and the pairing is fictional.
      2. Budgets must match per arm. An arm that could not fill its budget is
         VALID_BUDGET_LIMITED under HANDOFF 5.8 and must not have its mean compared directly
         against fully budgeted arms - which is exactly the error that produced the void
         "engine arm means at 117-124 measurements against baselines' 125".

    Returned as a list of problem strings, not raised: a reader needs to see which arm broke
    which premise, and on a 48-seed run a single bad arm must not suppress the other seven.
    """
    problems: List[str] = []
    arms = sorted(table)
    seeds = sorted(set.intersection(*(set(table[a]) for a in arms))) if arms else []
    reference = arms[0] if arms else None
    for seed in seeds:
        values = {a: table[a][seed].get("_initial_hits") for a in arms}
        distinct = {v for v in values.values() if v is not None}
        if len(distinct) > 1:
            problems.append(
                f"seed {seed}: round-0 top-1% hits differ across arms ({values}). The arms "
                f"did NOT share an initial sample and the paired test is invalid."
            )
    uniq = {a: float(np.mean([table[a][s]["n_unique_measurements"] for s in sorted(table[a])]))
            for a in arms}
    if uniq:
        top = max(uniq.values())
        for arm, value in sorted(uniq.items()):
            if top - value > 0.5:
                problems.append(
                    f"{arm}: mean unique measurements {value:.2f} against the run's best "
                    f"{top:.2f} - budget NOT matched. Under HANDOFF 5.8 this arm is "
                    f"VALID_BUDGET_LIMITED and its mean must not be compared directly "
                    f"against fully budgeted arms."
                )
    if reference is not None and not seeds:
        problems.append("no seed is shared by every arm; nothing can be paired.")
    return problems


def load_with_initial(run_dir: Path) -> Dict[str, Dict[int, dict]]:
    """`paired_stats.load_final_rounds`, plus round 0's hit count for the pairing check."""
    table = load_final_rounds(run_dir)
    initial: Dict[Tuple[str, int], float] = {}
    with (run_dir / "round_metrics.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(float(row["round"])) == 0:
                initial[(row["method"], int(float(row["seed"])))] = float(
                    row["cum_top1_hits_including_initial"]
                )
    for method, seeds in table.items():
        for seed, record in seeds.items():
            record["_initial_hits"] = initial.get((method, seed))
    return table


def ordered_arms(table: Dict[str, Dict[int, dict]]) -> List[str]:
    known = [a for a in PANEL_ORDER if a in table]
    return known + sorted(a for a in table if a not in PANEL_ORDER)


def pair(table, a: str, b: str, metric: str):
    """`a - b` over the seeds both reached. Returns (diff array, shared seeds)."""
    shared = sorted(set(table[a]) & set(table[b]))
    diff = np.array([table[a][s][metric] - table[b][s][metric] for s in shared], dtype=float)
    return diff, shared


def row(label: str, diff: np.ndarray) -> str:
    if diff.size == 0:
        return f"{label:<22}{'-':>9}  no shared seeds"
    lo, hi = paired_bootstrap_ci(diff, iterations=10_000)
    better = int((diff > 0).sum())
    worse = int((diff < 0).sum())
    tied = int((diff == 0).sum())
    t = paired_t(diff)
    if np.allclose(diff, 0.0):
        wil = float("nan")
    else:
        wil = float(stats.wilcoxon(diff, zero_method="wilcox").pvalue)
    return (f"{label:<22}{diff.size:>5}{diff.mean():>10.3f}"
            f"{(t if math.isfinite(t) else float('nan')):>9.2f}"
            f"{f'{better}/{worse}/{tied}':>12}"
            f"{f'[{lo:+.2f}, {hi:+.2f}]':>18}{wil:>11.3g}  {verdict(diff)}")


HEADER = (f"{'arm':<22}{'n':>5}{'mean diff':>10}{'t':>9}{'bet/wor/tie':>12}"
          f"{'95% CI':>18}{'wilcoxon p':>11}  verdict")


def report_one(run_dir: Path, metric: str, reference: str, engine: str) -> Optional[dict]:
    info = describe_run(run_dir)
    table = load_with_initial(run_dir)
    arms = ordered_arms(table)

    print("=" * 108)
    print(f"{info['library']}   run {run_dir.name}")
    print(f"  budget {info['budget']} = {info['initial']} initial + {info['rounds']} x {info['batch']}"
          f"   seeds {info['seeds']} from base {info['base_seed']}"
          f"   state {info['state']}" + (f" ({info['state_reason']})" if info["state_reason"] else ""))
    print(f"  configuration fingerprint {info['fingerprint']}  (quote the COMMIT, not the "
          f"fingerprint: HANDOFF 2.12)")
    if info["library"] in DEVELOPMENT_LIBRARIES:
        print(f"  !! {info['library']} is DEVELOPMENT DATA under rule 2.4. Anything below is "
              f"MECHANISM evidence, never generalization evidence.")
    print(f"  metric {metric}   reference baseline {reference} (registered; NOT promoted away from)")
    print("=" * 108)

    problems = check_pairing(table)
    if problems:
        print("\n  PAIRING / BUDGET PROBLEMS - read these before any number below:")
        for p in problems:
            print(f"    - {p}")
    else:
        print("\n  pairing premises hold: identical round-0 hit count per seed across arms, "
              "budgets matched per arm.")

    print(f"\n  Per-arm endpoint ({metric}, final round reached)")
    print(f"  {'arm':<22}{'n':>5}{'mean':>10}{'sd':>9}{'uniq.meas':>11}{'best fit':>11}  role")
    print("  " + "-" * 94)
    means: Dict[str, float] = {}
    for arm in arms:
        seeds = sorted(table[arm])
        values = np.array([table[arm][s][metric] for s in seeds], dtype=float)
        uniq = np.array([table[arm][s]["n_unique_measurements"] for s in seeds], dtype=float)
        best = np.array([table[arm][s]["best_measured_fitness"] for s in seeds], dtype=float)
        means[arm] = float(values.mean())
        role = ("ENGINE (under test)" if arm in ENGINE_ARMS
                else "REGISTERED REFERENCE" if arm == reference
                else "additional arm")
        print(f"  {arm:<22}{len(seeds):>5}{values.mean():>10.3f}{values.std(ddof=1):>9.3f}"
              f"{uniq.mean():>11.2f}{best.mean():>11.4f}  {role}")

    baselines = [a for a in arms if a not in ENGINE_ARMS]
    best_baseline = max(baselines, key=lambda a: means[a]) if baselines else None

    print(f"\n  A. Every arm against the registered reference `{reference}`   (arm - {reference})")
    print("  " + HEADER)
    print("  " + "-" * 104)
    for arm in arms:
        if arm == reference:
            continue
        diff, _ = pair(table, arm, reference, metric)
        print("  " + row(arm, diff))

    if engine in table:
        print(f"\n  B. The engine against every baseline   ({engine} - arm)")
        print("  " + HEADER)
        print("  " + "-" * 104)
        for arm in baselines:
            diff, _ = pair(table, engine, arm, metric)
            print("  " + row(arm, diff))

    # Separability of the reference against `random`: the question that makes a win over
    # the reference worth more or less than it reads.
    if "random" in table and reference in table and reference != "random":
        diff, _ = pair(table, reference, "random", metric)
        call = verdict(diff)
        print(f"\n  C. Is the registered reference separable from `random`?  "
              f"{reference} - random = {diff.mean():+.3f}, t={paired_t(diff):+.2f}, "
              f"{int((diff>0).sum())}/{int((diff<0).sum())}/{int((diff==0).sum())} -> {call}")
        if call == "not separable":
            print(f"     => on this landscape a win over `{reference}` is a win over a "
                  f"comparator that is NOT separable from chance. Any claim resting on it "
                  f"must say so.")
        # Landscape informativeness, the HER2 lesson applied before the fact.
        model_arms = [a for a in baselines if a != "random"]
        separable = [a for a in model_arms
                     if verdict(pair(table, a, "random", metric)[0]) != "not separable"]
        if not separable:
            print(f"     => NO model-based arm is separable from `random` here. This "
                  f"landscape's top-1% count cannot discriminate methods, so EVERY verdict "
                  f"on it is weak regardless of sign. That is a property of the library and "
                  f"no budget change fixes it.")

    survival = None
    if engine in table and best_baseline is not None:
        diff, shared = pair(table, engine, best_baseline, metric)
        survival = {
            "library": info["library"],
            "run": run_dir.name,
            "best_baseline": best_baseline,
            "best_baseline_mean": means[best_baseline],
            "engine_mean": means.get(engine),
            "n": diff.size,
            "delta": float(diff.mean()) if diff.size else float("nan"),
            "t": paired_t(diff),
            "better": int((diff > 0).sum()),
            "worse": int((diff < 0).sum()),
            "tied": int((diff == 0).sum()),
            "verdict": verdict(diff),
            "gate": None,
            "development": info["library"] in DEVELOPMENT_LIBRARIES,
            "ridge_vs_random": verdict(pair(table, reference, "random", metric)[0])
            if ("random" in table and reference in table and reference != "random") else None,
        }
        print(f"\n  D. SURVIVAL against the strongest baseline in the panel")
        print(f"     strongest baseline by its OWN mean: `{best_baseline}` at "
              f"{means[best_baseline]:.3f}   (engine {means.get(engine, float('nan')):.3f})")
        print("  " + HEADER)
        print("  " + row(f"{engine} - {best_baseline}", diff))
        print(f"     NOTE: a maximum over {len(baselines)} arms is upward-biased by the "
              f"winner's curse. Pre-declared in docs/prereg-baseline-panel.md section 4; "
              f"section A above carries every pairwise row.")
    return survival


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dirs", nargs="+")
    parser.add_argument("--metric", default="cum_top1_hits_excluding_initial")
    parser.add_argument("--reference", default="ridge")
    parser.add_argument("--engine", default="panoptes_default")
    args = parser.parse_args()

    survivals = []
    for raw in args.run_dirs:
        run_dir = Path(raw)
        if not run_dir.is_dir():
            raise SystemExit(f"not a directory: {run_dir.resolve()}")
        result = report_one(run_dir, args.metric, args.reference, args.engine)
        if result:
            survivals.append(result)
        print()

    if survivals:
        print("=" * 108)
        print("SURVIVAL TABLE - the engine against the strongest baseline in the panel")
        print("=" * 108)
        print(f"{'landscape':<14}{'best baseline':<18}{'base mean':>10}{'engine':>9}"
              f"{'delta':>9}{'t':>8}{'b/w/t':>12}  {'verdict':<14}notes")
        print("-" * 108)
        for s in survivals:
            notes = []
            if s["development"]:
                notes.append("DEVELOPMENT DATA (rule 2.4)")
            if s["ridge_vs_random"] == "not separable":
                notes.append("ridge NOT separable from random")
            bwt = "{}/{}/{}".format(s["better"], s["worse"], s["tied"])
            print(f"{str(s['library']):<14}{s['best_baseline']:<18}"
                  f"{s['best_baseline_mean']:>10.2f}{s['engine_mean']:>9.2f}"
                  f"{s['delta']:>9.2f}{s['t']:>8.2f}{bwt:>12}"
                  f"  {s['verdict']:<14}{'; '.join(notes)}")
        print("-" * 108)
        print("delta = engine - best baseline. A maximum over the panel is upward-biased "
              "for the baseline;\nsee section A of each landscape for every pairwise row.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
