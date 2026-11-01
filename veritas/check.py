"""`veritas check` - decide whether a proposed experiment can succeed, before paying for it.

WHAT THIS IS FOR
----------------
An agent proposes a protein engineering campaign: these positions, that mutation depth.
The proposal is usually sound biology and says nothing about whether the instrument in
hand can execute it. Three properties decide that, none of which is visible in the
sequence, the structure or the annotation the planner reads:

  1. SCOPE     - is every proposed position one the ordered library actually varies?
  2. REACH     - do the good variants live at or below the proposed mutation depth?
  3. VIABILITY - can the engine's own spatial viability veto pass any of them?

Each one is a *ceiling*: it binds before any ranking happens, so a perfect surrogate and
an unlimited budget do not move it. A campaign that fails one of these cannot succeed,
and the failure is invisible in the metrics everyone reports - mean top-1% yield and best
fitness are identical in shape for a search that lost and a search that was never
permitted to begin. This command computes the ceilings and refuses.

WHAT EACH CEILING IS A CEILING *ON* - read this before quoting REACH
--------------------------------------------------------------------
SCOPE and VIABILITY bind on this engine: an unmeasurable position cannot be assayed, and
the veto runs on every candidate inside `select_next_micro_batch`. **REACH does not.**
`max_mutation_depth` is a property of the PLAN, and this engine does not honour it -
`recombine_lead_variants` stacks measured leads and carries the candidate pool past the
declared depth from round 1. Measured 2026-11-01 over the real selection path on all
nine landscapes: the pre-veto ballot spans depths 1-12 on CreiLOV and reaches 99.9% of
its hit set, and the depth ceiling costs **1 hit of 1,655** there and **0** on the other
eight (`scripts/analysis/pool_hit_reachability.py`, on branch
`worktree-agent-a83565bcde59a81f3` @ `bb925c9`; not on `main` at the time of writing).

So a REACH refusal says "a proposer that honours this depth cannot succeed", which is
the right question to ask of a plan and of any system that enforces a depth bound - a
schema validator capping depth at 4 was binding in this repository historically - and is
NOT a prediction about this engine's yield. An earlier version of this module asserted
the stronger reading, and that was wrong. The check prints the distinction itself rather
than relying on a reader having seen this paragraph.

Measured cost of skipping it, on CR9114 before these checks existed: 757,829 candidate
proposals built, 25,712 rejections of which 100% were `not_in_eligible_universe`, 26 of
125 measurements spent, 0 of 651 top-1% hits found.

WHY ALL THREE, AND IN THIS ORDER
--------------------------------
They compose multiplicatively and they are not interchangeable. SCOPE is a property of
the ordered plate. REACH is a property of the measured landscape - it is not in the
protein. VIABILITY is a property of *this engine's own filter*, which no flag could
disable before commit `7415846`, and which on eqFP611-red removes all 82 of the
landscape's hits. An engine can therefore score zero while working perfectly, because
nothing it was permitted to select could have been a hit. Reporting only the composition
would hide which stage binds, so each is reported separately and then composed.

HOW IT STAYS HONEST
-------------------
* Depths and hit sets come from `veritas.calibration.landscapes`, never from a new
  genotype parser. There are already five canonicalisations in this repository and one of
  them silently inverted a published result (HANDOFF rule 2.1).
* The VIABILITY verdict is produced by calling the real
  `filter_spatially_viable_variants` through `scripts/analysis/spatial_veto_profile.py`,
  not a re-implementation of the distance rule.
* **The viability check audits its own input.** A 100% pass rate has two unrelated
  meanings - "the veto found nothing to reject" and "the structure never loaded, so the
  filter returned its input unchanged". Worse, a structure can load and be the *wrong
  structure*: HER2's `1N8Z` resolves to chain A, the trastuzumab LIGHT chain, while that
  library's positions are CDR-H3 on the HEAVY chain, and 9 of 10 mapped positions
  disagree with the library's own parent residue; CreiLOV's `8QI9` is phototropin /
  CrPhotLOV1, a non-aligning homolog, with 11 of 11 disagreeing. Production validates
  neither, because `reference_sequence` is not passed at the call site, so
  `SpatialContactMap`'s mismatch guard never fires. This command reports the mismatch
  rate and marks the viability verdict UNTRUSTED above `MISMATCH_TOLERANCE`, because a
  gate that cannot see its own bad input is the failure it exists to catch.
* No check is silently skipped. A check that could not run says so and why, and that
  propagates to the exit code - an unknown ceiling is not a passing one.

USAGE
-----
    python -m veritas.check --library creilov --depth 3
    python -m veritas.check --library trpb --positions 183,184,227,228 --depth 2
    python -m veritas.check --library creilov --plan plan.json
    python -m veritas.check --library trpb --from-model          # live Meta planner
    python -m veritas.check --library trpb --from-model --all-landscapes
    python -m veritas.check --library gb1 --depth 3 --json

Exit codes: 0 the plan can succeed - 1 the plan cannot succeed, or a ceiling is unknown
- 2 the check itself could not run. A gate that exits 0 on "I don't know" is not a gate.

ONE COSMETIC WART, DOCUMENTED BECAUSE IT CANNOT BE FIXED FROM HERE
------------------------------------------------------------------
Importing the Modal client prints two `DeprecationWarning` lines to **stderr** about
`asyncio.set_event_loop_policy`, before this tool has said anything. They are not
suppressible from this process: `warnings.filterwarnings`, `catch_warnings` +
`simplefilter("ignore", DeprecationWarning)`, and `python -W ignore::DeprecationWarning`
were each tried and each failed, so the client is forcing them past the filter
machinery. Redirect stderr (`2>$null`) if you need clean output - for a recording, say -
and note that doing so also hides this tool's own error messages, which go to stderr
too.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


import numpy as np  # noqa: E402

from veritas.calibration.landscapes import REGISTRY, load_landscape  # noqa: E402

#: A plan is sound only if it can reach essentially the whole hit set. These are the
#: same bands `scripts/demo/gate_catches_planner.py` reports, kept identical so the
#: demo and the tool cannot disagree about what a verdict means.
PASS_PCT = 99.0
DEGRADED_PCT = 50.0

#: Any disagreement at all between the structure's residue and the library's own parent
#: residue is worth printing: it means the deposited structure is not the molecule the
#: plate varies. 15% is `SpatialContactMap`'s own guard - the guard production never
#: fires, because the call site does not pass `reference_sequence`.
MISMATCH_CAVEAT = 0.15

#: Above this fraction the structure is not this library's molecule at all, and the
#: veto's verdict is not evidence about anything. The two are deliberately separate
#: bands, because collapsing them would grade three very different situations the same:
#: TrpB's 8EH0 maps 4 of 4 positions and disagrees at exactly one (184: P in the
#: structure, F in the library parent) - the right protein and the right numbering,
#: carrying one substitution, so the geometry is approximately right and the caveat is
#: a caveat. HER2's 1N8Z resolves to chain A, the trastuzumab LIGHT chain, while that
#: library's positions are CDR-H3 on the HEAVY chain: 9 of 10 disagree. CreiLOV's 8QI9
#: is phototropin / CrPhotLOV1, a non-aligning homolog: 11 of 11 disagree. Reporting
#: TrpB as badly compromised as those would be its own kind of dishonesty.
MISMATCH_UNTRUSTED = 0.50

PASS, DEGRADED, REFUSE, UNKNOWN = "PASS", "DEGRADED", "REFUSE", "UNKNOWN"

#: Structure grades, worst last. `grade_structure` is pure so the three bands can be
#: tested without parsing a PDB or running a campaign.
STRUCTURE_OK, STRUCTURE_CAVEAT, STRUCTURE_UNTRUSTED = "ok", "caveat", "untrusted"


def grade_structure(n_mapped: int, n_mismatch: int) -> str:
    """How much of the veto's verdict is evidence, given its own input?

    Returns one of `STRUCTURE_OK`, `STRUCTURE_CAVEAT`, `STRUCTURE_UNTRUSTED`. Pure, and
    deliberately coarse: the only question is whether the deposited structure is the
    molecule the plate varies. A majority of assayed positions disagreeing with the
    library's own parent residue answers no. A minority disagreeing is consistent with
    the right protein at the right numbering carrying a substitution.

    `n_mapped == 0` is UNTRUSTED rather than OK: nothing was checkable, and an unchecked
    input is not a clean one.
    """
    if n_mapped <= 0:
        return STRUCTURE_UNTRUSTED
    rate = n_mismatch / n_mapped
    if rate > MISMATCH_UNTRUSTED:
        return STRUCTURE_UNTRUSTED
    if rate > MISMATCH_CAVEAT:
        return STRUCTURE_CAVEAT
    return STRUCTURE_OK

RULE = "=" * 78
THIN = "-" * 78


def _render_detail(lines: Sequence[str]) -> List[str]:
    """Indent and wrap a check's detail.

    A line that already begins with whitespace is a pre-aligned table row and is passed
    through untouched; anything else is prose and gets wrapped. Doing it here rather than
    at the call sites means a check author writes a paragraph as one string and cannot
    accidentally hard-wrap it to a width that breaks when a number grows a digit.
    """
    out: List[str] = []
    for line in lines:
        if not line.strip():
            out.append("")
        elif line[:1].isspace():
            out.append(line)
        else:
            out.extend(textwrap.wrap(line, width=74, initial_indent="  ",
                                     subsequent_indent="  "))
    return out


def _verdict_for(pct: float) -> str:
    if not math.isfinite(pct):
        return UNKNOWN
    if pct >= PASS_PCT:
        return PASS
    if pct >= DEGRADED_PCT:
        return DEGRADED
    return REFUSE


# --------------------------------------------------------------------------------------
# the plan under test
# --------------------------------------------------------------------------------------
@dataclass
class Plan:
    """What the proposer asked for, and where it came from.

    `source` is carried through to the report because a plan's provenance changes how a
    refusal should be read: a hand-typed depth is a question, a live model's depth is a
    finding about the model.
    """

    positions: Tuple[int, ...]
    depth: int
    source: str
    hypothesis: str = ""
    rationales: Dict[str, str] = field(default_factory=dict)

    def to_json(self) -> Dict[str, Any]:
        return {
            "positions": list(self.positions), "max_mutation_depth": self.depth,
            "source": self.source, "hypothesis": self.hypothesis,
            "cluster_rationales": self.rationales,
        }


def plan_from_json(path: Path) -> Plan:
    """Read a plan from a file.

    Accepts either this tool's own `--json` output shape or a raw `SearchStrategy` dump,
    because the planner's own serialisation is the thing a user most naturally has on
    disk and making them reshape it invites transcription error.
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    blob = raw.get("plan", raw)
    depth = blob.get("max_mutation_depth", blob.get("depth"))
    if depth is None:
        raise ValueError(f"{path}: no `max_mutation_depth` (or `depth`) field")

    positions: List[int] = []
    if blob.get("positions"):
        positions = [int(p) for p in blob["positions"]]
    else:
        for cluster in blob.get("target_clusters") or []:
            positions.extend(int(p) for p in cluster.get("positions") or [])
    if not positions:
        raise ValueError(f"{path}: no positions, in `positions` or in `target_clusters`")

    rationales = {
        str(c.get("cluster_id", i)): str(c.get("rationale", ""))
        for i, c in enumerate(blob.get("target_clusters") or [])
    }
    return Plan(
        positions=tuple(sorted(set(positions))), depth=int(depth),
        source=f"file:{path.name}", hypothesis=str(blob.get("hypothesis", "")),
        rationales=rationales,
    )


def plan_from_model(library_key: str, *, withhold_scope: bool) -> Plan:
    """Ask the live Meta planner for a plan.

    `withhold_scope` controls whether the orchestrator is told which positions the
    library varies. Withholding it is the honest demonstration of what the gate is for:
    with the scope supplied, `restrict_strategy_to_assayable` clips the plan before
    anything downstream sees it - correct for a scored run, and it makes the catch
    invisible. There is no mock path here on purpose; a mock would prove nothing about
    what a model proposes.
    """
    from controller.orchestrator import Orchestrator, OrchestratorConfig

    bench = _bench()
    library = bench.LIBRARIES[library_key]
    metadata = bench.resolve_library_metadata(library)
    scope = None if withhold_scope else {int(p) for p in library.positions}

    orch = Orchestrator(
        OrchestratorConfig(
            target_accession=(library.accession or f"library:{library.key}"),
            use_mock_llm=False, enable_live_llm=True,
            require_live_llm=True,   # never silently substitute a mock
            random_seed=0,
        ),
        assayable_positions=scope,
    )
    strategy = orch._generate_strategy(metadata, round_number=1, top_evidence=[], current_hop=0)
    if getattr(orch, "strategy_engine_is_mock", None) is not False:
        raise RuntimeError(
            "the planner did not go live, so there is no model plan to gate; check "
            "LLAMA_API_KEY and VERITAS_USE_OS_TRUSTSTORE")

    positions: List[int] = []
    rationales: Dict[str, str] = {}
    for cluster in strategy.target_clusters:
        positions.extend(int(p) for p in cluster.positions)
        rationales[str(cluster.cluster_id)] = str(cluster.rationale or "")
    usage = dict(getattr(orch, "token_usage", {}) or {})
    engine = getattr(orch, "strategy_engine", "?")
    return Plan(
        positions=tuple(sorted(set(positions))), depth=int(strategy.max_mutation_depth),
        source=(f"live:{engine} ({usage.get('total_tokens', '?')} tokens in "
                f"{usage.get('calls', '?')} call(s), scope "
                f"{'withheld' if withhold_scope else 'supplied'})"),
        hypothesis=str(strategy.hypothesis or ""), rationales=rationales,
    )


# --------------------------------------------------------------------------------------
# lazily loaded heavy dependencies
# --------------------------------------------------------------------------------------
_BENCH = None
_PROFILE = None


def _bench():
    """Import `benchmark_panoptes_isolated` by path, once.

    It must be registered in `sys.modules` BEFORE `exec_module`, or dataclass field
    resolution fails with `'NoneType' object has no attribute '__dict__'`: `dataclasses`
    looks the defining module up by name while the class body executes. Loaded lazily
    because it is a large module and the SCOPE and REACH checks do not need it.
    """
    global _BENCH
    if _BENCH is None:
        if "benchmark_panoptes_isolated" in sys.modules:
            _BENCH = sys.modules["benchmark_panoptes_isolated"]
            return _BENCH
        spec = importlib.util.spec_from_file_location(
            "benchmark_panoptes_isolated", REPO_ROOT / "benchmark_panoptes_isolated.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules["benchmark_panoptes_isolated"] = module
        spec.loader.exec_module(module)
        _BENCH = module
    return _BENCH


def _profiler():
    """Import `scripts/analysis/spatial_veto_profile.py` by path, once.

    By path because `scripts/` is not a package. Reused rather than reimplemented: it
    calls the production `filter_spatially_viable_variants` and already validates itself
    against HANDOFF rule 2.2's shallow-shell anchors.
    """
    global _PROFILE
    if _PROFILE is None:
        _bench()   # its `import benchmark_panoptes_isolated` must resolve to our instance
        path = REPO_ROOT / "scripts" / "analysis" / "spatial_veto_profile.py"
        spec = importlib.util.spec_from_file_location("spatial_veto_profile", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules["spatial_veto_profile"] = module
        spec.loader.exec_module(module)
        _PROFILE = module
    return _PROFILE


# --------------------------------------------------------------------------------------
# the three checks
# --------------------------------------------------------------------------------------
@dataclass
class CheckOutcome:
    name: str
    verdict: str
    headline: str
    detail: List[str] = field(default_factory=list)
    numbers: Dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> Dict[str, Any]:
        return {"check": self.name, "verdict": self.verdict,
                "headline": self.headline, "detail": self.detail, **self.numbers}


def check_scope(plan: Plan, library_positions: Sequence[int]) -> CheckOutcome:
    """Is every proposed position one the ordered library actually varies?"""
    measurable: Set[int] = {int(p) for p in library_positions}
    proposed = sorted(set(plan.positions))
    outside = [p for p in proposed if p not in measurable]
    inside = [p for p in proposed if p in measurable]

    detail = [
        f"{p:>10}  {'MEASURABLE' if p in measurable else 'REFUSED':<12} "
        f"{'varied by the ordered library' if p in measurable else 'fixed here - no variant exists to assay'}"
        for p in proposed
    ]
    if not outside:
        verdict, headline = PASS, (
            f"all {len(proposed)} proposed positions are varied by this library")
        detail.append("")
        detail.append(
            "This check passing is worth a sentence rather than a shrug. The planner is "
            "handed this library's annotated sites in its prompt because `8099060` "
            "derived them from the library spec instead of per-target literals. Before "
            "that it proposed distal allosteric residues and whole rounds were discarded "
            "as not_in_eligible_universe. It passes because an earlier fix removed its "
            "cause, which is what a gate that has done its job looks like. It stays, "
            "because nothing guarantees the next target's annotation matches the library "
            "someone actually ordered.")
    else:
        verdict = REFUSE if not inside else DEGRADED
        headline = (f"{len(outside)} of {len(proposed)} proposed positions cannot be "
                    f"measured by this library: {outside}")
        detail.append("")
        detail.append(
            "Not bad biology - a distal residue may genuinely matter. Unmeasurable with "
            "the library in hand, and the model sees a protein, not a plate.")
        detail.append(f"plan as executable: {inside or 'nothing'}")

    return CheckOutcome(
        "scope", verdict, headline, detail,
        {"n_proposed": len(proposed), "n_measurable": len(inside),
         "unmeasurable_positions": outside, "library_positions": sorted(measurable)},
    )


def check_reach(plan: Plan, landscape) -> CheckOutcome:
    """Do the good variants live at or below the proposed mutation depth?

    Depth and the hit set are the landscape's own arrays, so this is its arithmetic
    rather than a second implementation of substitution counting.
    """
    depths = np.asarray(landscape.depths)
    hits = np.asarray(landscape.hit_mask)
    hit_depths = depths[hits]
    n_hits = int(hit_depths.size)
    if not n_hits:
        return CheckOutcome("reach", UNKNOWN,
                            "this landscape has no hits at the registered threshold", [],
                            {"n_hits": 0})

    reachable = int(np.count_nonzero(hit_depths <= plan.depth))
    pct = 100.0 * reachable / n_hits
    verdict = _verdict_for(pct)

    detail = [f"  depth {'hits':>10} {'cumulative':>12} {'% of hit set':>14}"]
    seen = 0
    for d in range(1, int(depths.max()) + 1):
        at_d = int(np.count_nonzero(hit_depths == d))
        if not at_d and d > plan.depth:
            continue
        seen += at_d
        marker = "  <- plan's ceiling" if d == plan.depth else ""
        detail.append(f"  {d:>5} {at_d:>10} {seen:>12} {100.0 * seen / n_hits:>13.1f}%{marker}")

    headline = (f"depth {plan.depth} reaches {reachable} of {n_hits} top-1% hits "
                f"({pct:.1f}%)")
    if verdict == REFUSE:
        detail.append("")
        detail.append(
            f"A proposer that HONOURS this depth finds at most {reachable}, whatever its "
            f"surrogate and whatever its budget. That is arithmetic about the plan.")
    if verdict != PASS:
        detail.append("")
        detail.append(
            "NOT ENFORCED BY THIS ENGINE - measured, 2026-11-01. `max_mutation_depth` "
            "does not cap the candidate pool: `recombine_lead_variants` stacks measured "
            "leads and carries the pool past the declared depth from round 1. "
            "Instrumented over the real selection path on all nine landscapes, the "
            "pre-veto ballot spans depths 1-12 on CreiLOV and reaches 99.9% of its hit "
            "set, and the depth ceiling costs 1 hit of 1,655 there and 0 on the other "
            "eight. So read this row as a property of the PLAN, not as a prediction "
            "about this engine's yield. It is the right question to ask of any proposer "
            "that enforces a depth bound - a schema validator capping depth at 4 was "
            "binding here historically - and it is the wrong question to ask of this "
            "one. Measured by `scripts/analysis/pool_hit_reachability.py` on branch "
            "`worktree-agent-a83565bcde59a81f3`.")
    return CheckOutcome("reach", verdict, headline, detail,
                        {"depth": plan.depth, "n_hits": n_hits,
                         "n_reachable": reachable, "pct_reachable": round(pct, 2)})


def check_viability(plan: Plan, library_key: str) -> CheckOutcome:
    """Can the engine's own spatial viability veto pass any reachable hit?

    This runs the production filter over the landscape's full universe and then asks the
    only question that matters for a plan: of the hits at or below the plan's depth, how
    many survive the veto. A hit the filter rejects is unreachable no matter how large
    the candidate pool is or how good the ranking is.
    """
    prof = _profiler()
    try:
        result = prof.profile_library(
            library_key,
            max_variants=prof.DEFAULT_SAMPLE_TRIGGER,
            sample_trigger=prof.DEFAULT_SAMPLE_TRIGGER,
            sample_size=prof.DEFAULT_SAMPLE_SIZE,
            seed=0,
        )
    except Exception as exc:
        return CheckOutcome(
            "viability", UNKNOWN,
            f"the veto could not be profiled: {type(exc).__name__}: {exc}",
            ["An unmeasured ceiling is not a passing one. Fix this before quoting a "
             "verdict - a silent skip here is how a filter that removes 100% of a hit "
             "set stayed invisible for the life of the project."],
            {},
        )

    struct = result.structure
    n_mapped = len(struct.mapped_positions)
    n_mismatch = len(struct.identity_mismatches)
    mismatch_rate = (n_mismatch / n_mapped) if n_mapped else float("nan")

    # Of the hits the plan can REACH, how many survive the veto? Per-depth rows carry
    # both counts, so this composes with check_reach instead of recomputing depth.
    reach_hits = sum(r.n_hits for r in result.rows if r.depth <= plan.depth)
    reach_pass = sum(r.n_hits_pass for r in result.rows if r.depth <= plan.depth)
    pct = 100.0 * reach_pass / reach_hits if reach_hits else float("nan")

    mismatch_note = (
        f"{n_mismatch} of {n_mapped} mapped positions DISAGREE "
        f"({', '.join(f'{p}:{s}!={l}' for p, s, l in struct.identity_mismatches)})"
        if n_mismatch else f"all {n_mapped} mapped positions agree with the library")
    detail = [
        f"  structure        {struct.verdict}  {struct.detail}",
        f"  positions mapped {n_mapped} of {len(struct.declared_positions)} declared",
        (f"  parent-residue   {mismatch_note}" if n_mapped
         else "  parent-residue   not checkable - no position mapped"),
        "",
        f"  {'depth':>5} {'hits':>8} {'pass veto':>11} {'%':>7}",
    ]
    for row in result.rows:
        if not row.n_hits:
            continue
        mark = "" if row.depth <= plan.depth else "   (beyond the plan)"
        detail.append(f"  {row.depth:>5} {row.n_hits:>8} {row.n_hits_pass:>11} "
                      f"{row.hit_pass_pct:>6.1f}%{mark}")

    verdict = _verdict_for(pct)
    headline = (f"of {reach_hits} hits at depth <= {plan.depth}, {reach_pass} can pass "
                f"the spatial veto ({pct:.1f}%)" if reach_hits
                else f"no hits at depth <= {plan.depth}, so the veto has nothing to pass")

    # The audit of the check's own input. It can only ever DOWNGRADE a verdict, never
    # raise one, so a structure problem can never manufacture a pass.
    grade = grade_structure(n_mapped, n_mismatch)
    untrusted = grade == STRUCTURE_UNTRUSTED
    caveat = grade == STRUCTURE_CAVEAT
    if untrusted and not n_mapped:
        verdict = UNKNOWN
        detail.append("")
        detail.append(
            f"UNTRUSTED: not one of this library's {len(struct.declared_positions)} "
            f"assayed positions could be located in the structure, so the cohesion veto "
            f"had no geometry to apply and `filter_spatially_viable_variants` returned "
            f"its input unchanged. A 100% pass rate here means 'the filter never ran', "
            f"not 'the filter found nothing to reject' - the two are indistinguishable "
            f"in the pass rate alone, which is exactly why this line exists. Structure "
            f"status: {struct.verdict} - {struct.detail}")
    elif untrusted:
        verdict = UNKNOWN
        detail.append("")
        detail.append(
            f"UNTRUSTED: {n_mismatch} of {n_mapped} mapped positions disagree with this "
            f"library's own parent residue ({100.0 * mismatch_rate:.0f}%). A majority "
            f"disagreeing means this is not the molecule the plate varies - the usual "
            f"causes are the wrong chain of a multi-chain entry and a homolog that does "
            f"not align. The veto is therefore measuring geometry that is not this "
            f"library's, and neither a pass nor a refusal here is evidence about "
            f"anything. Production does not catch it: `reference_sequence` is not passed "
            f"at the call site, so SpatialContactMap's mismatch guard never fires. "
            f"Resolve the structure before reading this row.")
    elif caveat:
        if verdict == PASS:
            verdict = DEGRADED
        detail.append("")
        detail.append(
            f"CAVEAT: the structure maps {n_mapped} of "
            f"{len(struct.declared_positions)} assayed positions and disagrees with the "
            f"library's parent residue at {n_mismatch} of them "
            f"({100.0 * mismatch_rate:.0f}%). A minority disagreement is consistent with "
            f"the right protein at the right numbering carrying a substitution - the "
            f"geometry is approximately right and the veto's verdict is usable with this "
            f"stated. Note the arithmetic is coarse on a small library: on four mapped "
            f"positions a single disagreement is 25% by division, not by severity. What "
            f"it rules out is quoting this row as if the structure were the library's "
            f"own parent.")
    if verdict == REFUSE:
        detail.append("")
        detail.append(
            "The veto is the engine's own filter, not a property of the library. It is "
            "ablatable since `7415846` (`--ablate-components spatial_filter`). Disabling "
            "it by default would be a scientific claim about whether distal co-mutations "
            "are jointly viable, not a tuning choice.")

    return CheckOutcome(
        "viability", verdict, headline, detail,
        {"n_hits_at_or_below_depth": reach_hits, "n_passing_veto": reach_pass,
         "pct_passing_veto": round(pct, 2) if math.isfinite(pct) else None,
         "structure_verdict": struct.verdict,
         "n_positions_mapped": n_mapped, "n_parent_residue_mismatches": n_mismatch,
         "parent_residue_mismatches": [
             {"position": p, "structure": s, "library_parent": l}
             for p, s, l in struct.identity_mismatches],
         "structure_untrusted": bool(untrusted),
         "structure_caveat": bool(caveat),
         "whole_universe_hit_pass_pct": round(result.hit_pass_pct, 2)
         if math.isfinite(result.hit_pass_pct) else None,
         "first_cohesion_biting_depth": result.first_cohesion_depth()},
    )


# --------------------------------------------------------------------------------------
# the composed report
# --------------------------------------------------------------------------------------
@dataclass
class Report:
    library: str
    label: str
    plan: Plan
    checks: List[CheckOutcome]
    attainable: Optional[int] = None
    n_hits: Optional[int] = None

    @property
    def verdict(self) -> str:
        """The worst verdict of any check. UNKNOWN never upgrades to PASS."""
        order = [PASS, DEGRADED, UNKNOWN, REFUSE]
        return max((c.verdict for c in self.checks), key=order.index)

    @property
    def exit_code(self) -> int:
        return 0 if self.verdict == PASS else 1

    def to_json(self) -> Dict[str, Any]:
        return {
            "library": self.library, "label": self.label, "plan": self.plan.to_json(),
            "verdict": self.verdict, "exit_code": self.exit_code,
            "attainable_hits": self.attainable, "n_hits": self.n_hits,
            "checks": [c.to_json() for c in self.checks],
        }

    def render(self) -> str:
        out: List[str] = []
        out.append(RULE)
        out.append(f"  veritas check - {self.library}  ({self.label})")
        out.append(RULE)
        out.append(f"  plan source        {self.plan.source}")
        out.append(f"  positions proposed {list(self.plan.positions)}")
        out.append(f"  max mutation depth {self.plan.depth}")
        # Wrapped rather than truncated: a plan's stated reasoning is the thing a reader
        # weighs the refusal against, and a sentence cut mid-word reads as the tool's
        # sloppiness rather than as the model's.
        if self.plan.hypothesis:
            out.extend(textwrap.wrap(
                self.plan.hypothesis, width=74, initial_indent="  hypothesis         ",
                subsequent_indent=" " * 21))
        for cid, why in self.plan.rationales.items():
            if why:
                out.extend(textwrap.wrap(
                    why, width=74, initial_indent=f"    [{cid}] ",
                    subsequent_indent=" " * 6))

        for check in self.checks:
            out.append("")
            out.append(f"{THIN}")
            out.append(f"  {check.name.upper():<10} {check.verdict:<9} {check.headline}")
            out.append(f"{THIN}")
            out.extend(_render_detail(check.detail))

        out.append("")
        out.append(RULE)
        out.append(f"  VERDICT  {self.verdict}")
        closing: List[str] = []
        if self.attainable is not None and self.n_hits:
            pct = 100.0 * self.attainable / self.n_hits
            closing.append(f"Composed ceiling: {self.attainable} of {self.n_hits} top-1% "
                           f"hits are attainable under this plan ({pct:.1f}%) - before "
                           f"any ranking happens.")
        if self.verdict == PASS:
            closing.append("Every ceiling checked is open. This plan is worth paying for.")
        elif self.verdict == UNKNOWN:
            closing.append("A ceiling could not be measured, so this plan is NOT cleared. "
                           "An unknown ceiling is not a passing one.")
        elif self.verdict == DEGRADED:
            closing.append("This plan can partly succeed. Part of the hit set is out of "
                           "its reach before any ranking happens, so a disappointing "
                           "result here is expected rather than informative - read the "
                           "ceiling above as the best this plan could have done.")
        else:
            closing.append("This plan cannot succeed as written. The agent plans; the "
                           "gate decides whether the plan is worth paying for.")
        out.extend(_render_detail(closing))
        out.append(RULE)
        return "\n".join(out)


def run_check(library_key: str, plan: Plan, *, skip_viability: bool = False) -> Report:
    if library_key not in REGISTRY:
        raise KeyError(f"unknown library {library_key!r}; known: {sorted(REGISTRY)}")
    spec = REGISTRY[library_key]
    landscape = load_landscape(library_key)

    checks = [check_scope(plan, spec.positions), check_reach(plan, landscape)]
    attainable: Optional[int] = checks[-1].numbers.get("n_reachable")
    if skip_viability:
        checks.append(CheckOutcome(
            "viability", UNKNOWN, "skipped by --skip-viability",
            ["The engine's spatial veto removes between 37% and 100% of the hit set "
             "depending on the landscape, so a plan cleared without it is not cleared."],
            {"skipped": True}))
        attainable = None
    else:
        viability = check_viability(plan, library_key)
        checks.append(viability)
        passing = viability.numbers.get("n_passing_veto")
        attainable = passing if isinstance(passing, int) else None

    return Report(library=library_key, label=spec.label, plan=plan, checks=checks,
                  attainable=attainable, n_hits=landscape.n_hits)


def render_all_landscapes(plan: Plan) -> str:
    """One plan, every landscape: does its depth reach each one's hit set?

    This is the check that generalises worst and surprises most. A depth a planner picks
    for `a protein` is sound on some measured landscapes and cannot succeed on others,
    and nothing in the sequence, the structure or the annotation says which.
    """
    out = ["", RULE, f"  ONE PLAN, EVERY LANDSCAPE - depth {plan.depth}", RULE,
           f"  {'library':<16}{'sites':>6}{'hits':>8}{'reachable':>11}{'%':>8}  verdict",
           f"  {THIN}"]
    rows: List[Tuple[str, int, int, int, float, str]] = []
    for key in sorted(REGISTRY):
        try:
            ls = load_landscape(key)
        except Exception as exc:
            out.append(f"  {key:<16}{'':>6}{'':>8}{'':>11}{'':>8}  UNAVAILABLE "
                       f"({type(exc).__name__})")
            continue
        hit_depths = np.asarray(ls.depths)[np.asarray(ls.hit_mask)]
        n_hits = int(hit_depths.size)
        reach = int(np.count_nonzero(hit_depths <= plan.depth))
        pct = 100.0 * reach / n_hits if n_hits else float("nan")
        verdict = _verdict_for(pct)
        rows.append((key, ls.spec.n_sites, n_hits, reach, pct, verdict))
        out.append(f"  {key:<16}{ls.spec.n_sites:>6}{n_hits:>8}{reach:>11}{pct:>7.1f}%  "
                   f"{verdict}")

    if rows:
        ok = [r for r in rows if r[5] == PASS]
        doomed = [r for r in rows if r[5] == REFUSE]
        best = max(rows, key=lambda r: r[4])
        tail = [
            f"One plan, {len(rows)} measured landscapes: sound on {len(ok)}, unable to "
            f"succeed on {len(doomed)}. Best case - {best[0]} - reaches {best[4]:.0f}% "
            f"of its own hits at depth {plan.depth}."]
        if doomed:
            worst = min(doomed, key=lambda r: r[4])
            tail.append(f"Worst - {worst[0]} - reaches {worst[3]} of {worst[2]}. A "
                        f"perfect surrogate and unlimited rounds still find {worst[3]}.")
        tail.append("")
        tail.append("The planner is not being careless. A depth like this is a sensible "
                    "exploration budget for a protein, and nobody told it how deep this "
                    "library's answers sit - that is not in the sequence, the structure "
                    "or the annotation. It is a property of the measured landscape, and "
                    "this gate is the only thing holding it.")
        out.append("")
        out.extend(_render_detail(tail))
    out.append(RULE)
    return "\n".join(out)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m veritas.check",
        description=("Decide whether a proposed campaign can succeed before paying for "
                     "it: scope, reach and viability, each a ceiling that binds before "
                     "any ranking happens."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Exit 0 = the plan can succeed; 1 = it cannot, or a ceiling is unknown; "
               "2 = the check could not run.",
    )
    parser.add_argument("--library", required=True, choices=sorted(REGISTRY),
                        help="the library in hand")
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--plan", type=Path, help="a plan as JSON (this tool's output or a "
                                               "raw SearchStrategy dump)")
    src.add_argument("--from-model", action="store_true",
                     help="ask the live Meta planner for a plan (one real API call)")
    parser.add_argument("--positions", default="",
                        help="comma-separated positions; defaults to the library's own")
    parser.add_argument("--depth", type=int, default=None,
                        help="proposed max mutation depth (required without --plan/--from-model)")
    parser.add_argument("--supply-scope", action="store_true",
                        help="with --from-model, tell the planner which positions are "
                             "assayable (as a scored run does). Withheld by default, "
                             "which is what makes the scope check observable.")
    parser.add_argument("--all-landscapes", action="store_true",
                        help="also apply this plan's depth to every registered landscape")
    parser.add_argument("--skip-viability", action="store_true",
                        help="skip the spatial veto profile (faster; the verdict is then "
                             "UNKNOWN, never PASS)")
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    return parser


def _plan_for(args: argparse.Namespace) -> Plan:
    """The plan under test, from whichever of the three sources was asked for."""
    if args.plan:
        return plan_from_json(args.plan)
    if args.from_model:
        if not (os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")):
            raise RuntimeError(
                "--from-model makes a real Meta Model API call and has no mock path, "
                "on purpose: a mock would prove nothing about what a model proposes. "
                "Set LLAMA_API_KEY or META_API_KEY.")
        return plan_from_model(args.library, withhold_scope=not args.supply_scope)
    if args.depth is None:
        raise RuntimeError("without --plan or --from-model, --depth is required.")
    positions = ([int(p) for p in args.positions.split(",") if p.strip()]
                 or list(REGISTRY[args.library].positions))
    return Plan(positions=tuple(sorted(set(positions))), depth=int(args.depth),
                source="hand-specified")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return _run(args)


def _run(args: argparse.Namespace) -> int:
    try:
        from dotenv import load_dotenv

        load_dotenv(REPO_ROOT / ".env", override=False)
    except ImportError:
        pass
    # Before any client opens a TLS connection, and before anything imports a worker
    # that memoises the outcome: `017e707` fixed a cache that kept a "not opted in"
    # answer forever, and `esm_worker` calls this at import time.
    from controller.tls_trust import enable_os_trust_store

    enable_os_trust_store()

    try:
        plan = _plan_for(args)
        report = run_check(args.library, plan, skip_viability=args.skip_viability)
    except (KeyError, ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f"veritas check could not run: {exc}", file=sys.stderr)
        return 2

    if args.json:
        payload = report.to_json()
        if args.all_landscapes:
            payload["all_landscapes_note"] = (
                "rendered in the text report; re-run without --json to read it")
        print(json.dumps(payload, indent=2))
    else:
        print(report.render())
        if args.all_landscapes:
            print(render_all_landscapes(plan))
    return report.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
