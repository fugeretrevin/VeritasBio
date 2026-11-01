"""Does the engine's REAL candidate ballot contain the hits at all? Per landscape, per round.

WHY THIS EXISTS
---------------
Two different ceilings can produce an identical "the engine lost" number, and no
existing artifact separates them:

  1. **Depth truncation.** `generate_candidate_variants` emits depth 1 and 2 only
     (stages 1, 3, 4 are all `itertools.combinations(..., 2)`); depth >= 3 arrives
     solely from `recombine_lead_variants` and from
     `sample_policy_gradient_candidates`, whose ceiling is
     `max_gfn_depth = max(strategy.max_mutation_depth, 3)` (math_worker.py:2167).
     `docs/NEXT_STEPS.md` section 3.2 measured `max_mutation_depth` to be INERT for
     the combinatorial constructs, so the sampler is what sets the reach.
  2. **Ballot composition.** The pool can reach the right depths and still be a
     sample of the library too small, or too poorly chosen, to contain many hits.

These are different stages and they compound. A depth histogram alone cannot tell
them apart, because "the pool reaches depth 9" says nothing about whether the 930
candidates at those depths include any of the 1,655 hits that also live there.

So this script measures, on a REAL campaign (not a synthetic call to the generator):

  A. `generate_candidate_variants` output            - the combinatorial generator alone
  B. `sample_policy_gradient_candidates` output      - the measured-dominant source
  C. the spatial veto's INPUT                        - the assembled, de-duplicated ballot
  D. the spatial veto's OUTPUT                       - what acquisition actually ranked

and for each stage reports the depth histogram, how many top-1% hits the stage
contains, and how many a uniform sample of the same size would have contained.
C vs D is the spatial veto's contribution, isolated from the generator's; A/B vs C
is the generator's, isolated from the veto's.

WHY THE STAGES HAD TO BE HOOKED SEPARATELY
------------------------------------------
`generation_coverage.json`, which the harness already writes, records the pool
`select_next_micro_batch` hands its trace sink - and that call happens at
math_worker.py:2443, AFTER the veto has replaced `novel_candidates` with
`filtered_candidates` (:2250). Every depth histogram in that file is therefore a
POST-veto histogram, so it cannot attribute a narrow ballot to generation rather
than to filtering. That is the exact confound this script removes, and it is why
stage C is taken from the filter's input rather than from the existing artifact.

HOW IT STAYS HONEST
-------------------
* Depths, the eligible universe and the hit set come from
  `veritas.calibration.landscapes.load_landscape` - the landscape's OWN arithmetic.
  No second genotype parser (HANDOFF rule 2.1); candidate notations are mapped with
  `veritas.calibration.genotypes.canonical_combo`, which fills wild-type at omitted
  positions and asserts full length, so a `zip()` truncation is impossible.
* HANDOFF rule 2.2's two anchors are reproduced and printed BEFORE any result, and
  the whole report is stamped unverified if either misses - the discipline
  `spatial_veto_profile.py` uses.
* The campaign is the production CLI path: this script calls
  `benchmark_panoptes_isolated.main` with the sweep's own budgets, having installed
  three recording wrappers that pass their arguments through unchanged and return the
  real functions' real return values. Nothing about selection is altered, and the
  observed arm result is printed so a perturbed run would be visible.
* `her2` is `pool_restricted`, so stages A and B are legitimately EMPTY there -
  `select_next_micro_batch` is handed the untested universe instead of generating.
  That is reported as `pool_restricted`, never as a generator that produced nothing.
* A candidate that is not a member of the eligible universe is counted as
  `off_universe` rather than silently dropped: it occupies a ballot slot and can
  never be a hit, so hiding it would flatter the ballot.

USAGE
-----
    .venv\\Scripts\\python.exe -X utf8 scripts/analysis/pool_hit_reachability.py
    ... --library creilov            # one landscape (repeatable)
    ... --json out.json              # machine-readable copy
    ... --ablate spatial_filter      # run the arm with the veto disabled

Each landscape runs in its own child process (`--child`), so one campaign's global
state cannot reach the next. Writes nothing except a run directory under
`benchmark_comparison_results/` and, with --json, the report.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from veritas.calibration.genotypes import canonical_combo, depth as combo_depth  # noqa: E402
from veritas.calibration.landscapes import REGISTRY, load_landscape  # noqa: E402

#: HANDOFF rule 2.2. {landscape: (hits, variants) in the shallow shell 1 <= depth <= 2}.
SHALLOW_SHELL_ANCHORS: Dict[str, Tuple[int, int]] = {"trpb": (168, 2_242), "gb1": (72, 2_167)}

#: The sweep's registered budgets (scripts/analysis/generality_sweep.ps1), so a pool
#: measured here is the pool of the campaign shape the sweep actually ran. Landscapes
#: the sweep dropped take 25/25 and are flagged: their ENDPOINT is degenerate at that
#: budget (CR6261's top-1% tier is 19 variants), but the ballot measurement is not an
#: endpoint measurement and stands.
BUDGETS: Dict[str, Tuple[int, int]] = {
    "trpb": (25, 25),
    "gb1": (25, 25),
    "her2": (25, 22),
    "eqfp611_blue": (17, 16),
    "eqfp611_red": (17, 16),
    "creilov": (25, 25),
    "cr9114_h1": (25, 25),
    "cr6261_h1": (25, 25),
    "cr6261_h9": (25, 25),
}
SWEEP_REGISTERED = {"her2", "eqfp611_blue", "eqfp611_red", "creilov", "cr9114_h1", "gb1"}

#: The derived-table acknowledgement the sweep passes, with its reason.
ACK_KEYS = (
    "cr6261_h1", "cr6261_h9", "creilov", "eqfp611_blue",
    "eqfp611_combined", "eqfp611_red", "gb1_combingym",
)
ACK_WHY = (
    "reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth "
    "at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean "
    "workbook which fell to zero at depth 16"
)

STAGES = ("generator", "sampler", "pre_veto", "post_veto")
STAGE_LABEL = {
    "generator": "A generate_candidate_variants",
    "sampler": "B sample_policy_gradient_candidates",
    "pre_veto": "C ballot, pre-veto (filter input)",
    "post_veto": "D ballot, post-veto (what acquisition ranked)",
}


# ==============================================================================
# Anchors
# ==============================================================================
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
        print("  *** The depth/hit pipeline does not reproduce its anchors. Nothing below "
              "is trustworthy. ***")
    return ok


# ==============================================================================
# Child: run one real campaign with the three stages recorded
# ==============================================================================
def run_child(key: str, out_path: Path, ablate: str, seed: int) -> int:
    """Install recording wrappers, run the production CLI, dump raw notations."""
    import benchmark_panoptes_isolated as bench
    from workers import math_worker as mw
    from workers import spatial_filter as sf

    # round_index is not available at these call sites, so rounds are numbered by the
    # order the filter is reached - one call per select_next_micro_batch invocation.
    captured: Dict[str, List[List[str]]] = {s: [] for s in STAGES}
    pending: Dict[str, List[str]] = {"generator": [], "sampler": []}

    real_generate = mw.generate_candidate_variants
    real_sample = mw.sample_policy_gradient_candidates
    real_filter = sf.filter_spatially_viable_variants

    def wrapped_generate(*a: Any, **kw: Any) -> Any:
        out = real_generate(*a, **kw)
        pending["generator"] = [c.notation for c in out]
        return out

    def wrapped_sample(*a: Any, **kw: Any) -> Any:
        out = real_sample(*a, **kw)
        pending["sampler"] = [c.notation for c in out]
        return out

    def wrapped_filter(mutation_strings: Sequence[str], *a: Any, **kw: Any) -> Any:
        approved, pruned = real_filter(list(mutation_strings), *a, **kw)
        # Only the selection call site closes a round. `generate_hamming1_lead_neighbors`
        # also calls this function (math_worker.py:1029) and is never reached in
        # production, but a wrapper that assumed so would silently mis-number rounds.
        caller = sys._getframe(1).f_code.co_name
        if caller == "select_next_micro_batch":
            captured["generator"].append(pending["generator"])
            captured["sampler"].append(pending["sampler"])
            captured["pre_veto"].append(list(mutation_strings))
            captured["post_veto"].append(list(approved))
            pending["generator"] = []
            pending["sampler"] = []
        return approved, pruned

    mw.generate_candidate_variants = wrapped_generate
    mw.sample_policy_gradient_candidates = wrapped_sample
    sf.filter_spatially_viable_variants = wrapped_filter

    initial, batch = BUDGETS[key]
    argv = [
        "--library", key,
        "--methods", "ridge,panoptes_default",
        "--initial-size", str(initial),
        "--batch-size", str(batch),
        "--seeds", "1",
        "--base-seed", str(seed),
        "--rounds", "4",
        "--provenance-mode", "acknowledge-unflagged",
    ]
    if ablate:
        argv += ["--ablate-components", ablate]

    rc = bench.main(argv)

    out_path.write_text(json.dumps({
        "library": key,
        "exit_code": int(rc),
        "initial_size": initial,
        "batch_size": batch,
        "ablate": ablate,
        "seed": seed,
        "pool_restricted": bool(bench.LIBRARIES[key].pool_restricted),
        "rounds": {s: captured[s] for s in STAGES},
    }), encoding="utf-8")
    return 0


# ==============================================================================
# Parent: score the captured notations against the landscape's own arrays
# ==============================================================================
def score_stage(
    notations: Sequence[str],
    *,
    positions: Sequence[int],
    parent: str,
    index: Dict[Tuple[str, ...], int],
    hit_mask: Any,
    n_universe: int,
    n_hits_total: int,
) -> Dict[str, Any]:
    """Depth histogram and hit content of one stage's pool, de-duplicated."""
    seen: set = set()
    depth_hist: Counter = Counter()
    hit_depth_hist: Counter = Counter()
    n_hits = n_off = n_unparseable = 0
    for notation in notations:
        combo = canonical_combo(notation, positions, parent)
        if combo is None:
            n_unparseable += 1
            continue
        if combo in seen:
            continue
        seen.add(combo)
        # The landscape's own arithmetic: WT filled at omitted positions, length
        # asserted inside canonical_combo, Hamming refuses unequal lengths.
        depth_hist[combo_depth(combo, parent)] += 1
        row = index.get(combo)
        if row is None:
            n_off += 1
            continue
        if bool(hit_mask[row]):
            n_hits += 1
            hit_depth_hist[combo_depth(combo, parent)] += 1
    n_distinct = len(seen)
    expected = n_distinct * n_hits_total / n_universe if n_universe else 0.0
    return {
        "n_notations": len(notations),
        "n_distinct": n_distinct,
        "n_off_universe": n_off,
        "n_unparseable": n_unparseable,
        "depth_histogram": {int(k): int(v) for k, v in sorted(depth_hist.items())},
        "depth_min": min(depth_hist) if depth_hist else None,
        "depth_max": max(depth_hist) if depth_hist else None,
        "n_hits_in_pool": n_hits,
        "hit_depth_histogram": {int(k): int(v) for k, v in sorted(hit_depth_hist.items())},
        "expected_hits_if_uniform": round(expected, 3),
        "hit_enrichment": round(n_hits / expected, 3) if expected > 0 else None,
    }


def parse_arm_yields(stdout: str) -> Dict[str, float]:
    """`{method: top1 hits}` from the harness's own summary table.

    Read from the run's own output rather than recomputed, so the yield quoted beside a
    ballot ceiling is the arm's real number and a campaign these wrappers had perturbed
    would show up as a changed yield instead of being assumed away.
    """
    out: Dict[str, float] = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 8 and parts[0] in ("ridge", "panoptes_default", "random") \
                and parts[1] == "available":
            try:
                out[parts[0]] = float(parts[5])
            except ValueError:
                continue
    return out


def analyse(key: str, raw: Dict[str, Any], yields: Dict[str, float]) -> Dict[str, Any]:
    landscape = load_landscape(key)
    spec = landscape.spec
    positions, parent = spec.positions, spec.parent_residues
    index = {combo: i for i, combo in enumerate(landscape.combos)}
    hit_mask = landscape.hit_mask

    hit_depth = Counter(int(d) for d in landscape.depths[hit_mask])
    universe_depth = Counter(int(d) for d in landscape.depths)

    per_stage: Dict[str, Any] = {}
    per_round: Dict[str, List[Dict[str, Any]]] = {}
    for stage in STAGES:
        rounds = raw["rounds"][stage]
        per_round[stage] = [
            score_stage(r, positions=positions, parent=parent, index=index,
                        hit_mask=hit_mask, n_universe=landscape.n_universe,
                        n_hits_total=landscape.n_hits)
            for r in rounds
        ]
        union: List[str] = [n for r in rounds for n in r]
        per_stage[stage] = score_stage(
            union, positions=positions, parent=parent, index=index, hit_mask=hit_mask,
            n_universe=landscape.n_universe, n_hits_total=landscape.n_hits,
        )
        support = set(per_stage[stage]["depth_histogram"])
        reachable = sum(c for d, c in hit_depth.items() if d in support)
        per_stage[stage]["hits_at_a_depth_the_pool_reaches"] = int(reachable)
        per_stage[stage]["pct_hits_at_a_reachable_depth"] = (
            round(100.0 * reachable / landscape.n_hits, 2) if landscape.n_hits else None
        )
        per_stage[stage]["pct_hits_actually_in_pool"] = (
            round(100.0 * per_stage[stage]["n_hits_in_pool"] / landscape.n_hits, 3)
            if landscape.n_hits else None
        )
        # The decisive composition number: the most top-1% hits a PERFECT ranker could
        # have taken from this stage's ballot, given the selection budget. It is
        # min(selection budget, distinct hits that were ever on the ballot) - an upper
        # bound, since it ignores that the same hit may recur across rounds. Compare it
        # against the arm's observed yield: a ceiling far above the yield means ranking
        # is the binding stage, and a ceiling at or below the baseline's yield means the
        # ballot is.
        selection_budget = 4 * raw["batch_size"]
        per_stage[stage]["selection_budget"] = selection_budget
        per_stage[stage]["perfect_ranker_ceiling"] = int(
            min(selection_budget, per_stage[stage]["n_hits_in_pool"])
        )

    return {
        "library": key,
        "label": spec.label,
        "n_sites": spec.n_sites,
        "n_universe": landscape.n_universe,
        "threshold": landscape.threshold,
        "threshold_source": landscape.threshold_source,
        "n_hits": landscape.n_hits,
        "hit_depth_histogram": {int(k): int(v) for k, v in sorted(hit_depth.items())},
        "universe_depth_histogram": {int(k): int(v) for k, v in sorted(universe_depth.items())},
        "hit_depth_median": int(sorted(int(d) for d in landscape.depths[hit_mask])
                                [landscape.n_hits // 2]) if landscape.n_hits else None,
        "pool_restricted": raw["pool_restricted"],
        "initial_size": raw["initial_size"],
        "batch_size": raw["batch_size"],
        "budget": raw["initial_size"] + 4 * raw["batch_size"],
        "budget_registered_by_sweep": key in SWEEP_REGISTERED,
        "ablate": raw["ablate"],
        "exit_code": raw["exit_code"],
        # One entry per select_next_micro_batch CALL, not per protocol round: a round
        # whose batch came up short is BACKFILLED from the same policy
        # (benchmark_panoptes_isolated.py:3690), and each backfill attempt is a second
        # call with the same evidence, hence an identical ballot. Calling these "rounds"
        # would double-count those ballots.
        "n_selector_calls": len(raw["rounds"]["pre_veto"]),
        "arm_yields_top1": yields,
        "stages": per_stage,
        "calls": per_round,
    }


# ==============================================================================
# Reporting
# ==============================================================================
RULE = "=" * 104
THIN = "-" * 104


def _hist(h: Dict[int, int]) -> str:
    if not h:
        return "(empty)"
    return " ".join(f"d{d}:{h[d]}" for d in sorted(h))


def print_library(r: Dict[str, Any]) -> None:
    print()
    print(RULE)
    print(f"LIBRARY  {r['library']}   {r['label']}")
    print(RULE)
    print(f"  universe            : {r['n_universe']:,} eligible, {r['n_sites']} sites")
    print(f"  top-1% threshold    : {r['threshold']:.6g}  [{r['threshold_source']}]")
    print(f"  hits                : {r['n_hits']:,}   median hit depth {r['hit_depth_median']}")
    print(f"  hit depths          : {_hist(r['hit_depth_histogram'])}")
    print(f"  campaign            : initial={r['initial_size']} batch={r['batch_size']} "
          f"budget={r['budget']}"
          + ("" if r["budget_registered_by_sweep"]
             else "   (NOT a sweep-registered budget; endpoint degenerate, ballot stands)")
          + (f"   ablate={r['ablate']}" if r["ablate"] else ""))
    if r["pool_restricted"]:
        print("  POOL-RESTRICTED     : select_next_micro_batch is handed the untested "
              "universe; stages A and B are empty BY DESIGN, not by failure.")
    yields = r["arm_yields_top1"]
    print("  arm yield (top-1%)  : "
          + ("  ".join(f"{m}={v:.2f}" for m, v in sorted(yields.items())) or "(not parsed)")
          + "   [includes the shared initial sample, as the harness reports it]")
    print(f"  selector calls      : {r['n_selector_calls']} "
          f"(4 rounds plus one call per backfill attempt; a backfilled call re-ranks the "
          f"same evidence, so its ballot repeats)")

    print()
    print(f"  {'stage':<42} {'distinct':>9} {'depths':>9} {'offuni':>7} "
          f"{'hits':>6} {'%hitset':>8} {'exp':>8} {'enrich':>7} {'%hits@reachable d':>18}")
    print("  " + THIN)
    for stage in STAGES:
        s = r["stages"][stage]
        span = ("-" if s["depth_min"] is None else f"{s['depth_min']}-{s['depth_max']}")
        enrich = "-" if s["hit_enrichment"] is None else f"{s['hit_enrichment']:.2f}"
        print(f"  {STAGE_LABEL[stage]:<42} {s['n_distinct']:>9,} {span:>9} "
              f"{s['n_off_universe']:>7,} {s['n_hits_in_pool']:>6,} "
              f"{(s['pct_hits_actually_in_pool'] or 0.0):>7.3f}% "
              f"{s['expected_hits_if_uniform']:>8.2f} {enrich:>7} "
              f"{(s['pct_hits_at_a_reachable_depth'] or 0.0):>17.2f}%")
    print()
    for stage in STAGES:
        s = r["stages"][stage]
        print(f"  {stage:<11} depth histogram: {_hist(s['depth_histogram'])}")
    print()
    pre, post = r["stages"]["pre_veto"], r["stages"]["post_veto"]
    if pre["n_distinct"]:
        kept = 100.0 * post["n_distinct"] / pre["n_distinct"]
        print(f"  VETO EFFECT ON THE BALLOT : {post['n_distinct']:,} of {pre['n_distinct']:,} "
              f"distinct candidates survive = {kept:.2f}%; hits in ballot "
              f"{pre['n_hits_in_pool']} -> {post['n_hits_in_pool']}")
    engine = r["arm_yields_top1"].get("panoptes_default")
    base = r["arm_yields_top1"].get("ridge")
    g = stage_gaps(r)
    print("  STAGE GAPS (hits, one unit): "
          + "  ".join(f"{k}={v:.0f}" for k, v in g["gaps"].items())
          + f"   ->  binds first: {g['binds_first']}")
    print(f"  CEILINGS at a {post['selection_budget']}-selection budget: a perfect ranker "
          f"takes <= {pre['perfect_ranker_ceiling']} hits from the PRE-veto ballot, "
          f"<= {post['perfect_ranker_ceiling']} from the POST-veto ballot"
          + (f"; the arm took {engine:.0f}" if engine is not None else "")
          + (f", ridge took {base:.0f}" if base is not None else ""))
    print("  PER SELECTOR CALL (pre-veto distinct / depth span / hits  ->  post-veto)")
    for i, (a, b) in enumerate(zip(r["calls"]["pre_veto"], r["calls"]["post_veto"]), start=1):
        sa = "-" if a["depth_min"] is None else f"{a['depth_min']}-{a['depth_max']}"
        sb = "-" if b["depth_min"] is None else f"{b['depth_min']}-{b['depth_max']}"
        print(f"    call {i}: {a['n_distinct']:>6,} / d{sa:<6} / {a['n_hits_in_pool']:>4} hits"
              f"   ->   {b['n_distinct']:>6,} / d{sb:<6} / {b['n_hits_in_pool']:>4} hits")


def print_summary(results: List[Dict[str, Any]]) -> None:
    print()
    print(RULE)
    print("SUMMARY - pool depth vs hit depth, and which stage binds first")
    print(RULE)
    print(f"  {'library':<14} {'hits':>6} {'hitd':>7} {'preD':>6} {'postD':>6} "
          f"{'%h@pre':>7} {'ceil_pre':>9} {'ceil_post':>10} {'eng':>5} {'ridge':>6} | "
          f"{'gap_dep':>7} {'gap_gen':>7} {'gap_veto':>8} {'gap_acq':>7} {'binds':>22}")
    print("  " + THIN)
    for r in results:
        pre, post = r["stages"]["pre_veto"], r["stages"]["post_veto"]
        hd = r["hit_depth_histogram"]
        hitd = f"{min(hd)}-{max(hd)}" if hd else "-"
        pd = "-" if pre["depth_min"] is None else f"{pre['depth_min']}-{pre['depth_max']}"
        qd = "-" if post["depth_min"] is None else f"{post['depth_min']}-{post['depth_max']}"
        eng = r["arm_yields_top1"].get("panoptes_default")
        rid = r["arm_yields_top1"].get("ridge")
        g = stage_gaps(r)
        gp = g["gaps"]
        s_eng = "-" if eng is None else f"{eng:.0f}"
        s_rid = "-" if rid is None else f"{rid:.0f}"
        s_acq = "-" if "acquisition" not in gp else f"{gp['acquisition']:.0f}"
        print(f"  {r['library']:<14} {r['n_hits']:>6,} {hitd:>7} {pd:>6} {qd:>6} "
              f"{(pre['pct_hits_at_a_reachable_depth'] or 0.0):>6.1f}% "
              f"{pre['perfect_ranker_ceiling']:>9} {post['perfect_ranker_ceiling']:>10} "
              f"{s_eng:>5} {s_rid:>6} | "
              f"{gp['depth']:>7.0f} {gp['generation']:>7.0f} {gp['veto']:>8.0f} "
              f"{s_acq:>7} {g['binds_first']:>22}")
    print()
    print("  hitd/preD/postD are depth SPANS. '%h@pre' is the share of the hit set sitting at")
    print("  a depth the PRE-VETO ballot reaches at least once: the candidate-generation depth")
    print("  ceiling, measured INDEPENDENTLY of the spatial veto so the two findings compose.")
    print("  'ceil_pre'/'ceil_post' are how many hits a PERFECT ranker could have taken from")
    print("  each ballot at this selection budget - the COMPOSITION ceiling.")
    print("  The four gap_* columns are the hits each stage costs, in one unit, and 'binds'")
    print("  names the largest. See stage_gaps() for the definitions.")
    print("  CAVEAT: one seed per landscape. A verdict needs a paired test (HANDOFF rule")
    print("  2.11); these are mechanism measurements, and a gap whose size is comparable to")
    print("  the arm's seed-to-seed spread must not be quoted as an effect.")


def stage_gaps(r: Dict[str, Any]) -> Dict[str, Any]:
    """How many hits each stage costs, in the same unit, so they can be ranked.

    The three are measured on one campaign and are additive by construction, because
    each is the drop in the SAME quantity - the most hits a perfect ranker could still
    take - across one stage boundary:

      depth       hits the pre-veto ballot cannot reach at any depth it ever proposes.
                  Capped at the selection budget, since hits beyond the budget were
                  never takeable anyway.
      generation  hits the pre-veto ballot could have reached by depth but did not
                  contain: budget - ceil_pre, net of the depth term. This is the
                  ballot being a small or badly chosen sample, not a shallow one.
      veto        ceil_pre - ceil_post: headroom the spatial viability filter removes.
      acquisition ceil_post - the arm's observed yield: hits that were on the ballot
                  acquisition ranked, and were not selected.

    Naming the largest is a measurement, not a judgement about which is easiest to fix.
    """
    pre, post = r["stages"]["pre_veto"], r["stages"]["post_veto"]
    budget = post["selection_budget"]
    n_hits = r["n_hits"]
    engine = r["arm_yields_top1"].get("panoptes_default")

    unreachable = n_hits - (pre["hits_at_a_depth_the_pool_reaches"] or 0)
    depth_gap = min(unreachable, budget)
    generation_gap = max(0, min(budget, n_hits) - depth_gap - pre["perfect_ranker_ceiling"])
    veto_gap = pre["perfect_ranker_ceiling"] - post["perfect_ranker_ceiling"]
    acq_gap = (None if engine is None
               else max(0.0, post["perfect_ranker_ceiling"] - engine))

    gaps = {"depth": float(depth_gap), "generation": float(generation_gap),
            "veto": float(veto_gap)}
    if acq_gap is not None:
        gaps["acquisition"] = float(acq_gap)
    largest = max(gaps, key=lambda k: gaps[k]) if any(v > 0 for v in gaps.values()) else "none"
    return {
        "gaps": gaps,
        "binds_first": "pool_restricted+" + largest if r["pool_restricted"] else largest,
    }


def _binding_stage(r: Dict[str, Any]) -> str:
    return stage_gaps(r)["binds_first"]


# ==============================================================================
def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("USAGE")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--library", action="append", default=None, metavar="KEY")
    parser.add_argument("--ablate", default="", metavar="SPEC",
                        help="forwarded as --ablate-components to every campaign.")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--child", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--child-out", type=Path, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.child:
        return run_child(args.child, args.child_out, args.ablate, args.seed)

    requested = args.library or [k for k in REGISTRY]
    unknown = [k for k in requested if k not in REGISTRY or k not in BUDGETS]
    if unknown:
        parser.error(f"unknown landscape(s) {unknown}; known: {sorted(REGISTRY)}")

    print(RULE)
    print("CANDIDATE BALLOT REACHABILITY - depth ceiling vs composition ceiling")
    print(RULE)
    print("  stages recorded  : A generate_candidate_variants, B "
          "sample_policy_gradient_candidates,")
    print("                     C spatial-veto input, D spatial-veto output")
    print("  depths/hit sets  : veritas.calibration.landscapes (HANDOFF rule 2.1 - no new "
          "genotype parser)")
    print(f"  landscapes       : {', '.join(requested)}")
    if args.ablate:
        print(f"  ablation         : --ablate-components {args.ablate}")
    anchors_ok = check_anchors()

    results: List[Dict[str, Any]] = []
    with tempfile.TemporaryDirectory() as tmp:
        for key in requested:
            out = Path(tmp) / f"{key}.json"
            cmd = [sys.executable, "-X", "utf8", str(Path(__file__).resolve()),
                   "--child", key, "--child-out", str(out), "--seed", str(args.seed)]
            if args.ablate:
                cmd += ["--ablate", args.ablate]
            env = dict(os.environ, PYTHONUTF8="1")
            proc = subprocess.run(cmd, cwd=str(REPO_ROOT), env=env,
                                  capture_output=True, text=True, errors="replace")
            if not out.is_file():
                print()
                print(f"FAILED {key}: child exited {proc.returncode} and wrote no capture "
                      f"to {out}")
                tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-12:]
                for line in tail:
                    print(f"    {line}")
                continue
            raw = json.loads(out.read_text(encoding="utf-8"))
            # The arm's own numbers, so a perturbed campaign would be visible rather
            # than assumed away.
            results.append(analyse(key, raw, parse_arm_yields(proc.stdout or "")))
            print_library(results[-1])

    if results:
        print_summary(results)
    if not anchors_ok:
        print()
        print("  *** rule-2.2 anchors did NOT reproduce; treat every table above as "
              "unverified. ***")

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({"anchors_ok": anchors_ok, "libraries": results},
                                        indent=2), encoding="utf-8")
        print(f"\nJSON written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
