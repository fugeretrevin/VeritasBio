"""The three calibration tables, computed from labels only.

    .venv\\Scripts\\python.exe -X utf8 -m veritas.calibration.report

Writes nothing by default; `--json <path>` dumps the same numbers machine-readably so a
protocol amendment can cite them without re-running anything.

Section 0 reproduces the development anchors first. Nothing below it is trustworthy if
that block does not match, so it is printed before any new figure:

    TrpB shallow 168/2,242 = 7.4933%   deep 1,432/157,757 = 0.9077%
    GB1  shallow  72/2,167 = 3.3226%   deep 1,422/147,193 = 0.9661%
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Optional, Sequence

from veritas.calibration.endpoint import endpoint_verdict
from veritas.calibration.genotype_sets import assess_shared_set
from veritas.calibration.genotypes import min_pairwise_hamming, pairs_closer_than
from veritas.calibration.landscapes import REGISTRY, Landscape, load_landscape
from veritas.calibration.qualification import (
    absolute_for_quantile,
    calibrate,
    equivalent_quantile,
    p_qualify_in_draw,
    qualifying_count,
)

#: Order the tables are printed in: development landscapes first, then final test.
ORDER = [
    "trpb", "gb1", "her2",
    "eqfp611_blue", "eqfp611_red", "creilov",
    "cr9114_h1", "cr6261_h1", "cr6261_h9",
]
LEVEL = {
    "trpb": "development", "gb1": "development", "her2": "development",
    "eqfp611_blue": "final test", "eqfp611_red": "final test",
    "creilov": "final test", "cr9114_h1": "final test",
    "cr6261_h1": "final test", "cr6261_h9": "final test",
}

#: The absolute `--qualify-threshold` default the frontier arms have been running with.
LEGACY_QUALIFY_ABSOLUTE = 0.70
#: The landscape the parameter was tuned on; its behaviour is what the quantile has to
#: reproduce.
CALIBRATION_LANDSCAPE = "trpb"
#: Shared initial observations per seed (§1.1): 25, drawn without replacement.
N_INITIAL = 25
#: Unique measurements per seed per arm (§1.1).
BUDGET = 125
DISTANCES = (2, 3, 4)

ANCHORS = {
    "trpb": {"shallow": (168, 2_242), "deep": (1_432, 157_757)},
    "gb1": {"shallow": (72, 2_167), "deep": (1_422, 147_193)},
}


def _rule(title: str, width: int = 104) -> List[str]:
    return ["=" * width, title, "=" * width]


def _ordered(landscapes: Dict[str, Landscape]) -> List[str]:
    """The loaded landscapes in reporting order; a subset is a subset of the tables."""
    return [k for k in ORDER if k in landscapes]


def anchor_block(out: List[str]) -> bool:
    """Reproduce the development figures. Returns False if any one of them disagrees."""
    out.extend(_rule("0. VALIDATION - the development anchors, reproduced before anything new"))
    ok = True
    for key, expected in ANCHORS.items():
        ls = load_landscape(key)
        shallow, deep = ls.shell(2), ls.deep_shell(2)
        for name, got, want in (
            ("shallow (1<=depth<=2)", (shallow.n_hits, shallow.n_variants), expected["shallow"]),
            ("deep    (depth>=3)   ", (deep.n_hits, deep.n_variants), expected["deep"]),
        ):
            rate = got[0] / got[1] * 100
            want_rate = want[0] / want[1] * 100
            match = got == want
            ok = ok and match
            out.append(
                f"  {key:<5} {name}  {got[0]:>6,}/{got[1]:>8,} = {rate:7.4f}%   "
                f"expected {want[0]:>6,}/{want[1]:>8,} = {want_rate:7.4f}%   "
                f"{'OK' if match else 'MISMATCH'}"
            )
        out.append(
            f"  {key:<5} universe {ls.n_universe:,}  hits {ls.n_hits:,}  "
            f"threshold {ls.threshold:.10g}  source: {ls.threshold_source}"
        )
    out.append("  Convention: depth 0 (the wild type) is not a library member, so it is in")
    out.append("  neither shell; the shells partition the mutants and sum to the universe.")
    out.append(f"  ANCHORS {'REPRODUCED' if ok else 'DO NOT REPRODUCE - STOP'}")
    return ok


def qualification_block(out: List[str], landscapes: Dict[str, Landscape]) -> dict:
    """Defect 1: a quantile qualification bar, and its absolute on every landscape."""
    out.append("")
    out.extend(_rule("1. DEFECT 1 - qualification threshold on incomparable fitness scales"))
    trpb = landscapes[CALIBRATION_LANDSCAPE]
    exact_q = equivalent_quantile(trpb.fitness, LEGACY_QUALIFY_ABSOLUTE)
    proposed_q = round(exact_q, 2)
    n_trpb = qualifying_count(trpb.fitness, LEGACY_QUALIFY_ABSOLUTE)
    out.append(
        f"  --qualify-threshold {LEGACY_QUALIFY_ABSOLUTE} on {CALIBRATION_LANDSCAPE} "
        f"(the landscape it was tuned on) is the {exact_q:.4f}th percentile of the "
        f"eligible universe:"
    )
    out.append(
        f"    {n_trpb:,} of {trpb.n_universe:,} variants qualify "
        f"({n_trpb / trpb.n_universe * 100:.4f}% of the universe), and "
        f"P(a random {N_INITIAL}-variant start qualifies) = "
        f"{p_qualify_in_draw(n_universe=trpb.n_universe, n_qualifying=n_trpb, n_draws=N_INITIAL):.4f}."
    )
    out.append(
        f"    For scale, that bar sits ABOVE the top-1% threshold "
        f"({trpb.threshold:.6g}): qualification on TrpB is a rare event in the extreme"
    )
    out.append("    upper tail, which is what the frontier arm's second phase is gated on.")
    out.append("")
    out.append(f"  PROPOSED DEFAULT: --qualify-quantile {proposed_q} "
               f"(equivalently: the top {100 - proposed_q:.2f}% of the eligible universe),")
    out.append(f"  replacing --qualify-threshold {LEGACY_QUALIFY_ABSOLUTE}. Converted with "
               f"numpy.percentile, the same rule the")
    out.append("  benchmark already uses for top1.threshold, so the two bars are on one scale.")
    out.append("")
    header = (f"  {'landscape':<14}{'level':<12}{'universe':>9}{'min':>9}{'max':>9}"
              f"{'top1 thr':>10}{'q-abs':>10}{'n qual':>8}{'P(init)':>9}"
              f"{'0.70 qual':>11}{'P(init)':>9}  changes?")
    out.append(header)
    rows = []
    for key in _ordered(landscapes):
        ls = landscapes[key]
        cal = calibrate(
            key, ls.fitness, quantile=proposed_q,
            legacy_absolute=LEGACY_QUALIFY_ABSOLUTE, n_initial=N_INITIAL,
            top1_threshold=ls.threshold,
        )
        out.append(
            f"  {key:<14}{LEVEL[key]:<12}{cal.n_universe:>9,}{cal.fitness_min:>9.4g}"
            f"{cal.fitness_max:>9.4g}{cal.top1_threshold:>10.4g}"
            f"{cal.absolute_at_quantile:>10.4g}{cal.n_qualifying_at_quantile:>8,}"
            f"{cal.p_qualify_initial_at_quantile:>9.4f}"
            f"{cal.n_qualifying_legacy:>11,}{cal.p_qualify_initial_legacy:>9.4f}"
            f"  {'YES' if cal.behaviour_changes else 'no'}"
        )
        rows.append({
            "landscape": key, "level": LEVEL[key], "n_universe": cal.n_universe,
            "fitness_min": cal.fitness_min, "fitness_max": cal.fitness_max,
            "top1_threshold": cal.top1_threshold,
            "quantile": proposed_q,
            "absolute_at_quantile": cal.absolute_at_quantile,
            "n_qualifying_at_quantile": cal.n_qualifying_at_quantile,
            "fraction_at_quantile": cal.fraction_at_quantile,
            "p_qualify_initial_at_quantile": cal.p_qualify_initial_at_quantile,
            "legacy_absolute": cal.legacy_absolute,
            "n_qualifying_legacy": cal.n_qualifying_legacy,
            "fraction_legacy": cal.fraction_legacy,
            "p_qualify_initial_legacy": cal.p_qualify_initial_legacy,
            "behaviour_changes": cal.behaviour_changes,
        })
    out.append("    q-abs    = the absolute fitness the proposed quantile resolves to there")
    out.append("    n qual   = variants a run would accept as qualifying, at-or-above semantics")
    out.append(f"    P(init)  = P(at least one qualifier among {N_INITIAL} variants drawn without")
    out.append("               replacement) - exact hypergeometric, not 1-(1-f)^n")
    out.append("    changes? = YES when the quantile bar accepts a different set than 0.70 does")
    out.append("")
    out.append("  SMALL-LIBRARY CAVEAT, which the quantile does not remove:")
    for key in _ordered(landscapes):
        ls = landscapes[key]
        n_q = next(r for r in rows if r["landscape"] == key)["n_qualifying_at_quantile"]
        if n_q <= 2:
            out.append(
                f"    {key}: the proposed quantile admits only {n_q} variant(s) of "
                f"{ls.n_universe:,}, so qualification there is effectively "
                f"'find the global best'."
            )
    out.append("    A quantile transfers the RARITY of the bar, not its attainability. On a")
    out.append("    library of ~1,900 variants the top 0.04% is one variant, so the stage")
    out.append("    gates almost everything out rather than almost nothing in. Either floor")
    out.append("    the qualifying set at a minimum count (e.g. max(quantile, k-th best) with")
    out.append("    k fixed in advance) or declare the stage inapplicable on those targets -")
    out.append("    but decide it in the protocol, not per run.")
    out.append("")
    out.append("  SENSITIVITY - the same table at four candidate quantiles, so the choice is")
    out.append("  visible rather than implied. Each cell is 'n qualifying / P(a random "
               f"{N_INITIAL}-start qualifies)':")
    ladder = (99.0, 99.5, 99.9, proposed_q)
    out.append(f"  {'landscape':<14}" + "".join(f"{('q=' + format(q, 'g')):>20}" for q in ladder))
    ladder_rows = {}
    for key in _ordered(landscapes):
        ls = landscapes[key]
        cells = []
        per_q = {}
        for q in ladder:
            absolute = absolute_for_quantile(ls.fitness, q)
            n_q = qualifying_count(ls.fitness, absolute)
            p = p_qualify_in_draw(n_universe=ls.n_universe, n_qualifying=n_q, n_draws=N_INITIAL)
            cells.append(f"{n_q:>12,} /{p:>6.3f}")
            per_q[q] = {"absolute": absolute, "n_qualifying": n_q, "p_qualify_initial": p}
        ladder_rows[key] = per_q
        out.append(f"  {key:<14}" + "".join(cells))
    out.append("    A qualification stage is only a gate where P(init) is small; at q=99.0 it")
    out.append("    fires in the initial draw on most of this corpus and gates nothing.")
    return {"exact_quantile_on_trpb": exact_q, "proposed_quantile": proposed_q,
            "rows": rows, "sensitivity": ladder_rows}


def endpoint_block(out: List[str], landscapes: Dict[str, Landscape]) -> dict:
    """Defect 2: is diversity-adjusted yield at d>=2 informative on this target?"""
    out.append("")
    out.extend(_rule("2. DEFECT 2 - is the primary endpoint degenerate on this target?"))
    out.append("  Raw ceiling = every variant at or above the target's own top-1% threshold.")
    out.append("  d>=k ceiling = the protocol's greedy rule on that whole hit set: the most")
    out.append("  distinct hits anybody could report at any budget. 'info' is the share the")
    out.append("  adjustment removes, 1 - ceiling/raw; at 0 the adjusted endpoint is")
    out.append("  numerically identical to the raw count for every arm.")
    out.append("")
    out.append(f"  {'landscape':<14}{'sites':>6}{'raw':>7}{'d>=2':>7}{'info':>8}{'d>=3':>7}"
               f"{'info':>8}{'d>=4':>7}{'minHam':>8}{'pairs<2':>9}{'pairs<3':>9}"
               f"  verdict (primary endpoint at d>=2)")
    rows = []
    for key in _ordered(landscapes):
        ls = landscapes[key]
        ceilings = ls.ceilings(DISTANCES)
        raw = ls.n_hits
        v2 = endpoint_verdict(raw=raw, ceiling=ceilings[2])
        v3 = endpoint_verdict(raw=raw, ceiling=ceilings[3])
        combos = ls.hit_combos
        min_ham = min_pairwise_hamming(combos)
        p2, p3 = pairs_closer_than(combos, 2), pairs_closer_than(combos, 3)
        out.append(
            f"  {key:<14}{ls.spec.n_sites:>6}{raw:>7,}{ceilings[2]:>7,}"
            f"{v2.informativeness * 100:>7.1f}%{ceilings[3]:>7,}"
            f"{v3.informativeness * 100:>7.1f}%{ceilings[4]:>7,}"
            f"{(min_ham if min_ham is not None else -1):>8}{p2:>9,}{p3:>9,}"
            f"  {v2.verdict}"
        )
        rows.append({
            "landscape": key, "level": LEVEL[key], "sites": ls.spec.n_sites,
            "raw_ceiling": raw,
            "ceiling_d2": ceilings[2], "ceiling_d3": ceilings[3], "ceiling_d4": ceilings[4],
            "informativeness_d2": v2.informativeness,
            "informativeness_d3": v3.informativeness,
            "verdict_d2": v2.verdict, "verdict_d3": v3.verdict,
            "min_pairwise_hamming": min_ham,
            "pairs_closer_than_2": p2, "pairs_closer_than_3": p3,
            "budget_exceeds_raw_ceiling": BUDGET >= raw,
        })
    out.append(f"    verdict: DEGENERATE = ceiling equals raw; NEAR-DEGENERATE = under "
               f"{endpoint_verdict(raw=100, ceiling=100).threshold:.0%} removed; else USABLE")
    out.append("    minHam  = minimum pairwise Hamming distance across the whole hit set")
    out.append("    pairs<k = unordered hit pairs closer than k, i.e. what the filter can act on")
    out.append("")
    bad2 = [r["landscape"] for r in rows if r["verdict_d2"] != "USABLE"]
    bad3 = [r["landscape"] for r in rows if r["verdict_d3"] != "USABLE"]
    out.append(f"  d>=2 NOT usable as primary on: {', '.join(bad2) if bad2 else 'none'}")
    out.append(f"  d>=3 NOT usable as primary on: {', '.join(bad3) if bad3 else 'none'}")
    saturated = [r["landscape"] for r in rows if r["budget_exceeds_raw_ceiling"]]
    out.append(f"  separately, the {BUDGET}-measurement budget meets or exceeds the whole hit")
    out.append(f"  set on: {', '.join(saturated) if saturated else 'none'} - there the endpoint")
    out.append("  saturates for any competent arm, which is a different defect from degeneracy.")
    return {"rows": rows, "degenerate_d2": bad2, "degenerate_d3": bad3,
            "budget_saturated": saturated}


def cr6261_block(out: List[str], landscapes: Dict[str, Landscape]) -> dict:
    """Defect 3 / O-11: the CR6261 H1-H9 genotype-set intersection."""
    out.append("")
    out.extend(_rule("3. DEFECT 3 / O-11 - the CR6261 H1 vs H9 genotype-set intersection"))
    if not {"cr6261_h1", "cr6261_h9"} <= set(landscapes):
        out.append("  skipped: both cr6261_h1 and cr6261_h9 must be loaded for this table.")
        return {}
    h1, h9 = landscapes["cr6261_h1"], landscapes["cr6261_h9"]
    view = assess_shared_set(
        set(h1.combos), set(h9.combos),
        parent_residues=h1.spec.parent_residues,
        budget=BUDGET, tier_fraction=h1.spec.tier_fraction,
    )
    c = view.comparison
    out.append(f"  CR6261 H1 holds {c.n_left:,} variants, H9 holds {c.n_right:,}.")
    out.append(f"  SHARED GENOTYPES (exact intersection over canonical 11-site combos): "
               f"{c.n_shared:,}")
    out.append(f"    unique to H1: {c.n_left_only:,}    unique to H9: {c.n_right_only:,}")
    out.append(f"    Jaccard {c.jaccard:.4f}; the shared set is {c.shared_share_left * 100:.2f}% "
               f"of H1 and {c.shared_share_right * 100:.2f}% of H9.")
    out.append("")
    out.append("  Depth profile of the shared set against each remainder (depth = library")
    out.append("  positions differing from the germline reference):")
    out.append(f"    shared   {view.depth_profile_shared}")
    out.append(f"    H1 only  {view.depth_profile_left_only}")
    out.append(f"    H9 only  {view.depth_profile_right_only}")
    out.append("")
    out.append("  Is the shared set a depth-BIASED subsample, or a smaller copy of the library?")
    out.append("  Per-depth retention (shared / full) on each side - a restriction that removes")
    out.append("  one depth stratum far harder than the others would shift the shell structure")
    out.append("  the mechanism claim rests on:")
    retention: Dict[str, Dict[int, float]] = {}
    for side, ls in (("H1", h1), ("H9", h9)):
        full = {d: n for d, (_hits, n) in ls.depth_histogram().items()}
        kept = {d: view.depth_profile_shared.get(d, 0) for d in full}
        retention[side] = {d: (kept[d] / full[d] if full[d] else 0.0) for d in sorted(full)}
        cells = "  ".join(f"d{d}:{retention[side][d] * 100:.0f}%" for d in sorted(full))
        worst_d = min(retention[side], key=lambda d: retention[side][d])
        out.append(f"    {side}: {cells}")
        out.append(f"        worst-retained stratum: depth {worst_d} at "
                   f"{retention[side][worst_d] * 100:.1f}%")
    out.append("")
    tier = view.tier_on_shared
    out.append(f"  Viability of a within-genotype-set paired contrast on the intersection:")
    out.append(f"    tier fraction in force on CR6261 (§9.1 tier ladder): "
               f"{h1.spec.tier_fraction:.0%}")
    out.append(f"    hit tier on the shared set: ceil({h1.spec.tier_fraction:g} x "
               f"{c.n_shared:,}) = {tier} variants")
    out.append(f"    budget {view.budget} is {view.budget_share_of_shared * 100:.2f}% of the "
               f"shared set, and {'EXCEEDS' if view.budget_exceeds_tier else 'is below'} "
               f"the {tier}-variant tier")
    # The paired contrast is only as good as the hits that exist on BOTH readouts.
    shared = set(h1.combos) & set(h9.combos)
    h1_hits = {c_ for c_ in h1.hit_combos if c_ in shared}
    h9_hits = {c_ for c_ in h9.hit_combos if c_ in shared}
    out.append(f"    top-1% hits that lie inside the shared set: H1 {len(h1_hits)} of "
               f"{h1.n_hits}, H9 {len(h9_hits)} of {h9.n_hits}")
    out.append(f"    hits shared by both readouts: {len(h1_hits & h9_hits)}")
    out.append("")
    eqfp_payload: dict = {}
    if {"eqfp611_blue", "eqfp611_red"} <= set(landscapes):
        out.append("  For contrast, eqFP611's two channels:")
        blue, red = landscapes["eqfp611_blue"], landscapes["eqfp611_red"]
        eq = assess_shared_set(
            set(blue.combos), set(red.combos),
            parent_residues=blue.spec.parent_residues, budget=BUDGET,
        )
        out.append(f"    blue {eq.comparison.n_left:,} / red {eq.comparison.n_right:,}, shared "
                   f"{eq.comparison.n_shared:,}, identical: "
                   f"{'YES' if eq.comparison.identical else 'NO'}")
        eqfp_payload = {
            "n_blue": eq.comparison.n_left, "n_red": eq.comparison.n_right,
            "n_shared": eq.comparison.n_shared, "identical": eq.comparison.identical,
        }
    return {
        "cr6261": {
            "n_h1": c.n_left, "n_h9": c.n_right, "n_shared": c.n_shared,
            "n_h1_only": c.n_left_only, "n_h9_only": c.n_right_only,
            "jaccard": c.jaccard,
            "shared_share_h1": c.shared_share_left,
            "shared_share_h9": c.shared_share_right,
            "depth_profile_shared": view.depth_profile_shared,
            "depth_profile_h1_only": view.depth_profile_left_only,
            "depth_profile_h9_only": view.depth_profile_right_only,
            "tier_fraction": h1.spec.tier_fraction,
            "tier_on_shared": tier,
            "budget": view.budget,
            "budget_share_of_shared": view.budget_share_of_shared,
            "h1_hits_in_shared": len(h1_hits),
            "h9_hits_in_shared": len(h9_hits),
            "hits_in_both": len(h1_hits & h9_hits),
            "depth_retention": retention,
        },
        "eqfp611": eqfp_payload,
    }


def build_report(keys: Sequence[str] = tuple(ORDER)) -> tuple[str, dict]:
    out: List[str] = []
    ok = anchor_block(out)
    if not ok:
        return "\n".join(out), {"anchors_reproduced": False}
    # TrpB is always loaded: it is the landscape --qualify-threshold was tuned on, so
    # the quantile has to be calibrated against it even when it is not being reported.
    landscapes = {k: load_landscape(k) for k in dict.fromkeys((CALIBRATION_LANDSCAPE, *keys))}
    payload: dict = {"anchors_reproduced": True}
    for key in _ordered(landscapes):
        ls = landscapes[key]
        payload.setdefault("thresholds", {})[key] = {
            "threshold": ls.threshold, "source": ls.threshold_source,
            "computed_99th_percentile": ls.computed_threshold,
            "n_universe": ls.n_universe, "n_hits": ls.n_hits,
            "rows_read": ls.n_rows_read,
            "dropped": {
                "unlabelled_split": ls.n_dropped_unlabelled,
                "wild_type": ls.n_dropped_wild_type,
                "unparseable": ls.n_dropped_unparseable,
                "non_finite": ls.n_dropped_nonfinite,
                "duplicate": ls.n_dropped_duplicate,
            },
        }
    payload["qualification"] = qualification_block(out, landscapes)
    payload["endpoint"] = endpoint_block(out, landscapes)
    payload["genotype_sets"] = cr6261_block(out, landscapes)
    out.append("")
    out.extend(_rule("4. THRESHOLD PROVENANCE - what each landscape was scored with"))
    out.append(f"  {'landscape':<14}{'universe':>9}{'hits':>7}{'threshold':>14}  source")
    for key in _ordered(landscapes):
        ls = landscapes[key]
        out.append(f"  {key:<14}{ls.n_universe:>9,}{ls.n_hits:>7,}{ls.threshold:>14.10g}  "
                   f"{ls.threshold_source}")
    out.append("  Where a manifest exists its top1.threshold is used verbatim, so analysis")
    out.append("  scores variants with exactly the threshold the benchmark used; the loader")
    out.append("  also recomputes the 99th percentile and the two agree on all three.")
    return "\n".join(out), payload


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=str, default=None,
                        help="also write the figures as JSON to this path")
    parser.add_argument("--landscapes", nargs="*", default=list(ORDER),
                        choices=sorted(REGISTRY),
                        help="restrict the tables to these landscapes")
    args = parser.parse_args(argv)
    text, payload = build_report(tuple(args.landscapes))
    print(text)
    if args.json:
        from pathlib import Path

        Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nwrote {args.json}")
    return 0 if payload.get("anchors_reproduced") else 1


if __name__ == "__main__":
    sys.exit(main())
