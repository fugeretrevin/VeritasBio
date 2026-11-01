"""Paired read-out for the spatial-veto ablation campaign, ACROSS two run directories.

Why this exists rather than `paired_stats.py` alone. `--ablate-components` is a single
global flag applied to every panoptes arm in a run, so the unablated and ablated engine
arms cannot share a run directory. `paired_stats.py` pairs arms *within* one run folder.
This script pairs the same seed across *two* folders, and it does so by importing
`paired_stats.load_final_rounds` and `paired_stats.paired_bootstrap_ci` rather than
re-deriving them: HANDOFF §4 names one read-out recipe and a second would drift from it
(rule 2.11 is the most repeated error in this project).

Three guards, each of which refuses rather than warns:

  1. **Pairing.** `ridge` is carried in BOTH runs. If ridge's per-seed vector is not
     bit-identical across the pair, the two runs did not see the same initial samples and
     nothing is paired; the comparison is refused. This is the whole reason ridge is paid
     for twice.
  2. **The ablation actually happened.** The ablated run's `manifest.json` must record
     `components.states.spatial_filter == false` and the unablated run's must record
     `true`. HANDOFF rule 2.12: a manifest that reads only the CLI flag is worse than no
     record. `workers.math_worker` is imported and checked for any module-level global
     that could gate `spatial_filter` behind the manifest's back; if one appears, this
     refuses, because the manifest would then be reporting a state it cannot see.
  3. **Same library, same budget.** Compared arms must be the same landscape at the same
     budget (HANDOFF §5.8: only equal-budget arms enter a matched comparison).

Usage:
  python veto_ablation_readout.py --landscape creilov \
      --unablated <run_dir> --ablated <run_dir> [--metric cum_top1_hits_excluding_initial]

  python veto_ablation_readout.py --spec creilov=<unablated>,<ablated> \
      --spec gb1=<unablated>,<ablated>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paired_stats import load_final_rounds, paired_bootstrap_ci  # noqa: E402

DEFAULT_METRIC = "cum_top1_hits_excluding_initial"

#: The published baselines the pre-registration's falsification thresholds are stated
#: against (docs/prereg-spatial-veto-ablation.md §2 and §5). Carried here so the verdict
#: is computed against a committed number, not one retyped at read-out time.
PREREG: Dict[str, Dict[str, float]] = {
    "eqfp611_red": {"panoptes": 0.0000, "ridge": 15.3542, "confirm_at": 7.6771},
    "creilov": {"panoptes": 6.9375, "ridge": 57.2292, "confirm_at": 32.0834},
    "gb1": {"panoptes": 20.3958, "ridge": 33.2917, "confirm_at": 26.8438},
}


class ReadoutError(RuntimeError):
    pass


def _manifest(run_dir: Path) -> dict:
    path = run_dir / "manifest.json"
    if not path.is_file():
        raise ReadoutError(f"no manifest at {path} (checked: {path})")
    return json.loads(path.read_text(encoding="utf-8"))


def _check_no_hidden_global() -> None:
    """Rule 2.12: the manifest is only as honest as what it reads.

    `component_states` consults module globals for `dpp` and `dkl` only. If some future
    global gates the spatial veto as well, the manifest's claim becomes unverifiable from
    the manifest, so this refuses instead of trusting it.
    """
    from workers import math_worker as mw

    suspects = [
        name for name in dir(mw)
        if name.isupper() and "SPATIAL" in name and isinstance(getattr(mw, name), bool)
    ]
    if suspects:
        raise ReadoutError(
            "workers.math_worker now has module-level boolean global(s) "
            f"{suspects} mentioning SPATIAL. `component_states()` reads only "
            "DPP_DIVERSITY_ENABLED and TWO_STAGE_GP_ENABLED, so the manifest's "
            "spatial_filter state can no longer be trusted. Teach component_states "
            "about it before quoting any number from these runs."
        )
    # Positive check: the gate the engine actually uses must still be the ablation set.
    source = Path(mw.__file__).read_text(encoding="utf-8")
    if 'if "spatial_filter" in ablated:' not in source:
        raise ReadoutError(
            "the `if \"spatial_filter\" in ablated:` gate is gone from "
            f"{mw.__file__}; the ablation arm's meaning cannot be verified."
        )


def _verify_pair(landscape: str, unablated: Path, ablated: Path, metric: str):
    mu, ma = _manifest(unablated), _manifest(ablated)

    for tag, man, want in (("unablated", mu, True), ("ablated", ma, False)):
        states = man.get("configuration", {}).get("components", {}).get("states", {})
        if "spatial_filter" not in states:
            raise ReadoutError(
                f"{tag} manifest has no components.states.spatial_filter; the run "
                "predates the ablatable veto and cannot be used as either arm.")
        got = bool(states["spatial_filter"])
        if got is not want:
            raise ReadoutError(
                f"{tag} run {Path(man.get('run_dir', '?')).name or '?'} records "
                f"spatial_filter enabled={got}, expected {want}. The arms are mislabelled.")

    # Same landscape, same budget (rule 5.8).
    for key in ("library", "initial_size", "batch_size", "rounds"):
        a = _find(mu, key)
        b = _find(ma, key)
        if a != b:
            raise ReadoutError(
                f"{landscape}: the two runs disagree on {key!r} ({a!r} vs {b!r}); "
                "they are not a matched-budget pair.")

    tu = load_final_rounds(unablated)
    ta = load_final_rounds(ablated)
    for tag, table in (("unablated", tu), ("ablated", ta)):
        if "ridge" not in table:
            raise ReadoutError(
                f"{tag} run has no `ridge` arm, so the pairing check cannot run. "
                "Both runs must carry ridge.")
    seeds = sorted(set(tu["ridge"]) & set(ta["ridge"]))
    if not seeds:
        raise ReadoutError(f"{landscape}: the two runs share no seeds.")
    ru = np.array([tu["ridge"][s][metric] for s in seeds])
    ra = np.array([ta["ridge"][s][metric] for s in seeds])
    if not np.array_equal(ru, ra):
        bad = [s for s, x, y in zip(seeds, ru, ra) if x != y]
        raise ReadoutError(
            f"{landscape}: ridge is NOT bit-identical across the two runs "
            f"({len(bad)} of {len(seeds)} seeds differ, e.g. {bad[:5]}). The runs did "
            "not see the same initial samples, so the engine arms are not paired. "
            "REFUSING to report. Do not quote this landscape.")
    return tu, ta, seeds


def _find(man: dict, key: str):
    """Locates a scalar config key wherever the manifest schema keeps it."""
    stack = [man]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            if key in node and not isinstance(node[key], (dict, list)):
                return node[key]
            stack.extend(v for v in node.values() if isinstance(v, (dict, list)))
        elif isinstance(node, list):
            stack.extend(v for v in node if isinstance(v, (dict, list)))
    return None


def _paired(a: np.ndarray, b: np.ndarray) -> dict:
    """Paired statistics on a - b, with the seed as the unit of independence."""
    diff = a - b
    n = diff.size
    lo, hi = paired_bootstrap_ci(diff)
    sd = float(diff.std(ddof=1)) if n > 1 else 0.0
    t = float(diff.mean() / (sd / np.sqrt(n))) if sd > 0 else float("nan")
    wins, losses = int((diff > 0).sum()), int((diff < 0).sum())
    wil = (float("nan") if np.allclose(diff, 0)
           else float(stats.wilcoxon(a, b, zero_method="wilcox").pvalue))
    sign = (float(stats.binomtest(wins, wins + losses, 0.5).pvalue)
            if wins + losses else float("nan"))
    return {"n": n, "mean": float(diff.mean()), "sd": sd, "t": t,
            "ci": (lo, hi), "w": wins, "l": losses,
            "tie": int((diff == 0).sum()), "wilcoxon": wil, "sign": sign}


def _row(label: str, s: dict) -> str:
    return (f"  {label:<34}{s['n']:>4}{s['mean']:>10.4f}{s['t']:>9.2f}"
            f"{f'[{s['ci'][0]:+.3f}, {s['ci'][1]:+.3f}]':>22}"
            f"{s['wilcoxon']:>11.3g}{s['sign']:>10.3g}"
            f"{f'{s['w']}/{s['l']}/{s['tie']}':>12}")


def report(landscape: str, unablated: Path, ablated: Path, metric: str) -> None:
    tu, ta, seeds = _verify_pair(landscape, unablated, ablated, metric)

    eng_u = np.array([tu["panoptes_default"][s][metric] for s in seeds])
    eng_a = np.array([ta["panoptes_default"][s][metric] for s in seeds])
    ridge = np.array([tu["ridge"][s][metric] for s in seeds])

    print()
    print("=" * 118)
    print(f"LANDSCAPE {landscape}   metric {metric}   {len(seeds)} paired seeds "
          "(unit of independence: seed)")
    print(f"  unablated run: {unablated.name}")
    print(f"  ablated   run: {ablated.name}")
    print("  manifests verified: spatial_filter enabled=True / False; ridge "
          "bit-identical across both runs")
    print("=" * 118)
    print(f"  {'arm':<34}{'n':>4}{'mean':>10}{'sd':>9}")
    for name, v in (("ridge", ridge), ("panoptes_default (veto ON)", eng_u),
                    ("panoptes_default (veto ABLATED)", eng_a)):
        print(f"  {name:<34}{v.size:>4}{v.mean():>10.4f}{v.std(ddof=1):>9.4f}")

    print()
    print(f"  {'paired comparison':<34}{'n':>4}{'mean diff':>10}{'t':>9}"
          f"{'95% boot CI':>22}{'wilcoxon p':>11}{'sign p':>10}{'W/L/T':>12}")
    print("  " + "-" * 114)
    s_abl = _paired(eng_a, eng_u)
    s_u_ridge = _paired(eng_u, ridge)
    s_a_ridge = _paired(eng_a, ridge)
    print(_row("ablated - unablated", s_abl))
    print(_row("unablated - ridge", s_u_ridge))
    print(_row("ablated - ridge", s_a_ridge))

    pre = PREREG.get(landscape)
    print()
    print("  VERDICT (pre-registered thresholds, docs/prereg-spatial-veto-ablation.md "
          "section 5)")
    gap0 = s_u_ridge["mean"]
    gap1 = s_a_ridge["mean"]
    closed = gap1 - gap0
    frac = closed / abs(gap0) if gap0 != 0 else float("nan")
    print(f"    gap to ridge, veto ON      : {gap0:+.4f}")
    print(f"    gap to ridge, veto ABLATED : {gap1:+.4f}")
    print(f"    closed by ablation         : {closed:+.4f}  "
          f"= {100.0 * frac:+.1f}% of the gap")
    if pre:
        print(f"    pre-registered confirm bar : ablated mean >= {pre['confirm_at']:.4f} "
              f"(half the published gap of {pre['panoptes'] - pre['ridge']:+.4f})")
        print(f"    observed ablated mean      : {eng_a.mean():.4f}")

    if not np.isfinite(frac):
        verdict = "UNDEFINED (no gap to close)"
    elif frac >= 0.50 and s_abl["w"] >= 36 and abs(s_abl["t"]) > 3:
        verdict = "H-VETO CONFIRMED: the veto is THE CAUSE on this landscape"
    elif frac < 0.0:
        verdict = ("NOT A CAUSE - and the veto was NET BENEFICIAL here "
                   "(ablation made it worse)")
    elif frac < 0.10 or (np.isfinite(s_abl["t"]) and abs(s_abl["t"]) < 2):
        verdict = "NOT A CAUSE"
    elif 0.10 <= frac < 0.50 and s_abl["w"] >= 30:
        verdict = "CONTRIBUTING CAUSE, not the main one"
    else:
        verdict = "FALSIFIED as a cause (gap not halved, or no paired separation)"
    print(f"    ==> {verdict}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--landscape")
    ap.add_argument("--unablated")
    ap.add_argument("--ablated")
    ap.add_argument("--spec", action="append", default=[],
                    help="landscape=<unablated_dir>,<ablated_dir>; repeatable")
    ap.add_argument("--metric", default=DEFAULT_METRIC)
    args = ap.parse_args()

    _check_no_hidden_global()

    jobs: List[Tuple[str, Path, Path]] = []
    if args.landscape:
        if not (args.unablated and args.ablated):
            raise SystemExit("--landscape needs both --unablated and --ablated")
        jobs.append((args.landscape, Path(args.unablated), Path(args.ablated)))
    for spec in args.spec:
        key, _, dirs = spec.partition("=")
        u, _, a = dirs.partition(",")
        jobs.append((key, Path(u), Path(a)))
    if not jobs:
        raise SystemExit("nothing to do: pass --landscape/--unablated/--ablated or --spec")

    failures = 0
    for key, u, a in jobs:
        try:
            report(key, u, a, args.metric)
        except ReadoutError as exc:
            failures += 1
            print()
            print(f"!! {key}: REFUSED - {exc}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
