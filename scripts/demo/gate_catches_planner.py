"""The agent proposes. The gate decides whether the proposal can be paid for.

Runs the live Meta planner UNCONSTRAINED on a real library, prints what it proposed, then
applies the assayable-scope gate and prints the verdict position by position.

Why unconstrained. In the benchmark the orchestrator is built with
`assayable_positions=library.positions`, so `restrict_strategy_to_assayable` clips the
plan before anything downstream sees it - correct for a scored run, and it means the catch
is invisible. Here the scope is withheld on the first call so the planner proposes freely,
exactly as it would against a library whose geometry nobody declared. That is not a
strawman: it is the documented failure mode. From the benchmark's own limitations block:

    "panoptes_muse is given this library's saturated positions as its assayable scope,
     which is enforced deterministically on the planner's output. Without it the planner
     proposed distal allosteric residues the four-site library cannot measure and whole
     rounds were discarded as not_in_eligible_universe."

The model is not wrong about the biology. A distal allosteric residue may well matter.
It is wrong about the instrument: the library in hand varies four sites, and a proposal
anywhere else cannot be measured, ordered or scored. Only the gate knows that.

Usage:
    .venv\\Scripts\\python.exe -X utf8 scripts/demo/gate_catches_planner.py
    .venv\\Scripts\\python.exe -X utf8 scripts/demo/gate_catches_planner.py --library gb1

Requires LLAMA_API_KEY (or META_API_KEY) in .env, and VERITAS_USE_OS_TRUSTSTORE=1 on a
network that terminates TLS. Costs one live planner call, ~5k tokens.

Antibody libraries are deliberately not offered: a prompt carrying affinity evidence
against influenza hemagglutinin trips the provider's content policy, which is correct
behaviour and not something to design around. See docs/HANDOFF.md.
"""
from __future__ import annotations

import argparse
import importlib.util
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Set

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
os.chdir(REPO)

#: Enzyme and fluorescent-protein libraries only - see the module docstring.
OFFERED = ("trpb", "gb1", "creilov", "eqfp611_blue", "eqfp611_red")

RULE = "=" * 78
THIN = "-" * 78


def banner(text: str) -> None:
    print(f"\n{RULE}\n  {text}\n{RULE}")


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", default="trpb", choices=OFFERED)
    ap.add_argument("--verbose", action="store_true",
                    help="show the orchestrator's own log lines")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")

    try:
        from dotenv import load_dotenv

        load_dotenv(REPO / ".env", override=False)
    except ImportError:
        pass
    from controller.tls_trust import enable_os_trust_store

    enable_os_trust_store()

    if not (os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")):
        print("No LLAMA_API_KEY or META_API_KEY found. This script makes a real call to "
              "the Meta Model API and has no mock path, on purpose: a mock would prove "
              "nothing about what the model proposes.")
        return 2

    spec = importlib.util.spec_from_file_location(
        "bench", REPO / "benchmark_panoptes_isolated.py")
    bench = importlib.util.module_from_spec(spec)
    sys.modules["bench"] = bench
    spec.loader.exec_module(bench)

    from controller.orchestrator import Orchestrator, OrchestratorConfig
    from controller.prompts import restrict_strategy_to_assayable

    library = bench.LIBRARIES[args.library]
    metadata = bench.resolve_library_metadata(library)
    measurable: Set[int] = {int(p) for p in library.positions}

    banner(f"THE LIBRARY IN HAND - {args.library}")
    print(f"  protein            {getattr(metadata, 'name', args.library)}")
    print(f"  sequence length    {len(getattr(metadata, 'sequence', '') or '')} residues")
    print(f"  positions VARIED   {sorted(measurable)}")
    print(f"  everything else    fixed - cannot be ordered, measured or scored")
    print(f"\n  This is the physical constraint. The assay reads "
          f"{len(measurable)} positions and no others.")

    # --- the agent proposes, with no idea what the instrument can read ----------------
    banner("THE AGENT PROPOSES  (live Meta Model API, scope withheld)")
    orch = Orchestrator(
        OrchestratorConfig(
            target_accession=(library.accession or f"library:{library.key}"),
            use_mock_llm=False,
            enable_live_llm=True,
            require_live_llm=True,   # never silently substitute a mock
            random_seed=0,
        ),
        assayable_positions=None,   # withheld: this is the whole point
    )
    print("  calling muse-spark-1.3 ...", flush=True)
    strategy = orch._generate_strategy(
        metadata, round_number=1, top_evidence=[], current_hop=0)

    live = getattr(orch, "strategy_engine_is_mock", None) is False
    print(f"  engine             {getattr(orch, 'strategy_engine', '?')}"
          f"   (mock: {not live})")
    usage = dict(getattr(orch, "token_usage", {}) or {})
    print(f"  tokens             {usage.get('total_tokens', '?')} "
          f"in {usage.get('calls', '?')} call(s)")
    if not live:
        print("\n  The planner did not go live, so there is nothing to gate. Check "
              "VERITAS_USE_OS_TRUSTSTORE and the API key.")
        return 3

    print(f"\n  hypothesis         {(strategy.hypothesis or '')[:200]}")
    print(f"  max_mutation_depth {strategy.max_mutation_depth}")
    print(f"\n  proposed clusters:")
    proposed: List[int] = []
    for cluster in strategy.target_clusters:
        print(f"    [{cluster.cluster_id}]  positions {cluster.positions}")
        print(f"        {(cluster.rationale or '')[:150]}")
        proposed.extend(int(p) for p in cluster.positions)

    # --- check 1: is every proposed position measurable? ------------------------------
    banner("GATE CHECK 1 - CAN THE ASSAY READ WHAT IT PROPOSED?")
    unique = sorted(set(proposed))
    outside = [p for p in unique if p not in measurable]

    for p in unique:
        verdict = "MEASURABLE" if p in measurable else "REFUSED"
        why = ("varied by the ordered library" if p in measurable
               else "fixed in this library - no variant exists to assay")
        print(f"  {p:>10}  {verdict:<14} {why}")

    if outside:
        print(f"\n  {len(outside)} of {len(unique)} refused: {outside}")
        print("  Not bad biology - a distal residue may genuinely matter. Unmeasurable")
        print("  with the library in hand, and the model sees a protein, not a plate.")
        clipped = restrict_strategy_to_assayable(strategy, measurable)
        kept = sorted({int(p) for c in clipped.target_clusters for p in c.positions})
        print(f"  plan as executable: {kept}")
    else:
        print(f"\n  All {len(unique)} measurable - this check PASSES, and it is worth")
        print("  saying why rather than quietly moving on. The planner is handed this")
        print("  library's annotated sites in its prompt, because `8099060` derived the")
        print("  planner's positions from the library spec instead of per-target")
        print("  literals. Before that it proposed distal allosteric residues and whole")
        print("  rounds were discarded as not_in_eligible_universe. The check now passes")
        print("  because an earlier fix removed its cause, which is what a gate that has")
        print("  done its job looks like. It stays, because nothing guarantees the next")
        print("  target's annotation matches the library someone actually ordered.")

    # --- check 2: can the plan's DEPTH reach where the answers live? -------------------
    banner("GATE CHECK 2 - CAN THAT PLAN REACH THE ANSWER?")
    depth = int(strategy.max_mutation_depth)
    print(f"  The model asked for max_mutation_depth = {depth}.")
    print(f"  A proposer with reach {depth} can only build variants carrying at most")
    print(f"  {depth} substitutions. So: where do the good variants actually live?\n")

    from veritas.calibration.landscapes import REGISTRY, load_landscape

    print(f"  {'library':<16}{'sites':>6}{'hits':>8}{'reachable':>11}{'% ':>7}  verdict")
    print(f"  {THIN}")
    rows = []
    for key in ("trpb", "gb1", "eqfp611_blue", "creilov", "cr6261_h1", "cr9114_h1",
                "cr6261_h9"):
        if key not in REGISTRY:
            continue
        ls = load_landscape(key)
        # `depths` is the per-variant substitution count and `hit_mask` marks the
        # top-tier variants, both aligned to the eligible universe - so this is the
        # landscape's own arithmetic rather than a second implementation of depth.
        hit_depths = [int(d) for d, is_hit in zip(ls.depths, ls.hit_mask) if is_hit]
        n_hits = len(hit_depths)
        reach = sum(1 for d in hit_depths if d <= depth)
        pct = 100.0 * reach / n_hits if n_hits else 0.0
        verdict = ("OK" if pct >= 99.0 else
                   "DEGRADED" if pct >= 50.0 else "REFUSE - cannot succeed")
        n_sites = len(ls.spec.positions)
        rows.append((key, n_sites, n_hits, reach, pct, verdict))
        print(f"  {key:<16}{n_sites:>6}{n_hits:>8}{reach:>11}{pct:>6.1f}%  {verdict}")

    ok = [r for r in rows if r[4] >= 99.0]
    doomed = [r for r in rows if r[4] < 50.0]
    best = max(rows, key=lambda r: r[4]) if rows else None
    print(f"\n  One model call, one plan, {len(rows)} libraries: sound on {len(ok)},")
    print(f"  unable to succeed on {len(doomed)}.")
    if best:
        print(f"  Even the best case - {best[0]} - reaches only {best[4]:.0f}% of its own")
        print(f"  good variants at depth {depth}.")
    if doomed:
        worst = min(doomed, key=lambda r: r[4])
        print(f"  Worst - {worst[0]} - reaches {worst[3]} of {worst[2]}. A perfect")
        print(f"  surrogate and unlimited rounds still find {worst[3]}.")
    print(f"\n  The model is not being careless. Depth {depth} is a sensible exploration")
    print(f"  budget for a protein, and nobody told it how deep this library's answers")
    print(f"  sit - that is not in the sequence, the structure or the annotation. It is a")
    print(f"  property of the measured landscape, and the gate is the only thing holding")
    print(f"  it. Which is the point: the plan is not wrong, it is unverifiable from the")
    print(f"  plan alone, and the fix is to check it before paying rather than to ask the")
    print(f"  model to somehow know.")

    # --- what it costs to skip the gate ----------------------------------------------
    banner("WHAT SKIPPING BOTH CHECKS COSTS")
    print("  Measured on CR9114 before these checks existed:")
    print("    757,829 candidate proposals, none above mutation depth 4")
    print("     25,712 rejections, 100% for one reason: not_in_eligible_universe")
    print("         26 of 125 measurements spent; 192 of 192 rounds committed short")
    print("          0 of 651 top-1% hits found")
    print("\n  This harness DID catch it: run state PARTIAL, arm DEGENERATE. Both checks")
    print("  exist because of this episode - a run demotes to PARTIAL on any degenerate")
    print("  arm, and the per-arm classification is veritas/preflight.py check 3.")
    print("\n  What did not catch it is what everyone reports: mean top-1% yield and best")
    print("  fitness. Those are identical in shape for an arm that searched and lost and")
    print("  an arm that never searched. The rejection reason that separates them was in")
    print("  the run's own artifacts and fed no reported metric.")
    print("\n  The agent plans. The gate decides whether the plan is worth paying for.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
