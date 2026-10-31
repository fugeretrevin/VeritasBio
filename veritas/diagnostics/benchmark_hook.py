"""
veritas/diagnostics/benchmark_hook.py - the single call a benchmark makes.

Footprint is the point. `benchmark_panoptes_isolated.py` gains one call:

    from veritas.diagnostics.benchmark_hook import emit_coverage_artifacts
    emit_coverage_artifacts(run_dir, outcomes, eligible_ids=universe.variant_ids,
                            stamp=stamp, top1_predicate=evaluator.is_top1,
                            run_id=run_id)

and everything else lives here, so a concurrent refactor of that file's library
globals merges cleanly. The library definition is INFERRED from the eligible variant
ids (every token carries its own wild-type residue), so this hook references no
module-level `ACTIVE_LIBRARY` / `ACTIVE_SITE_POSITIONS` / `PARENT_RESIDUES` constant.

It also never raises. A diagnostic that can fail a benchmark run is a worse
diagnostic than none: every failure here is logged and returns an empty result.

Outputs written next to the existing CSVs, never over them:
    generation_coverage.csv    per method per round coverage + depth histograms
    generation_coverage.json   the same, plus the collapse diagnosis per arm
    depth_returns.csv          tidy long per-depth returns, one row per
                               run/method/seed/round/depth
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence

from veritas.diagnostics.coverage import (
    CollapseThresholds,
    CoverageAccumulator,
    CoverageDiagnosis,
    assert_generation_coverage,
)
from veritas.diagnostics.depth_returns import (
    DepthReturnLedger,
    MeasuredSelection,
    write_depth_returns_csv,
)
from veritas.diagnostics.genotype import LibrarySpec, infer_library_from_variant_ids
from veritas.diagnostics.trace import write_candidate_trace_csv

logger = logging.getLogger("veritas.diagnostics.benchmark_hook")

# Statuses used by `selection_trace.csv`. `initial` is the shared starting sample,
# which the harness hands to every arm: it is not something candidate generation
# produced, so counting it would flatter every arm by the same 25 variants.
STATUS_SELECTED = "accepted"
STATUS_REJECTED = "rejected"
STATUS_GIVEN = "initial"

# Every row carries the CONTINUOUS quantities alongside the label, so a threshold can
# be revised later from the written artifacts without re-running a single campaign.
# The label is a convenience; `union_share`, `depth_reach` and `deep_enrichment` are
# the data.
COVERAGE_CSV_COLUMNS = (
    "provenance",
    "method",
    "round",
    "campaigns",
    "proposed_distinct",
    "proposed_total",
    "cumulative_distinct",
    "cumulative_share",
    "selected_distinct",
    "cumulative_selected_distinct",
    "cumulative_mean_proposed_depth",
    "proposed_depth_histogram",
    "selected_depth_histogram",
    "pool_source",
    "arm_classification",
    "arm_union_distinct",
    "arm_union_share",
    "arm_mean_proposed_depth",
    "arm_depth_reach",
    "arm_deep_enrichment",
    "arm_deep_proposed",
    "arm_deep_universe_share",
    "threshold_min_union_share",
    "threshold_min_depth_reach",
    "threshold_min_deep_enrichment",
    "threshold_min_proposals_to_judge",
)


def _trace_rows(outcomes: Iterable[Any]) -> Iterable[Any]:
    for outcome in outcomes or ():
        for row in getattr(outcome, "trace", ()) or ():
            yield row


class CandidatePoolRecorder:
    """
    Captures the TRUE generated candidate pool, per (method, seed, round).

    Why this is not optional
    ------------------------
    `selection_trace.csv` holds accepted picks plus harness-level rejections, i.e. the
    selector's committed batch - not its ballot. Coverage computed from it is therefore
    a LOWER BOUND on what was proposed, and a lower bound cannot distinguish the two
    hypotheses that matter:

        (a) a narrow GENERATOR that only ever builds a few hundred genotypes;
        (b) a broad generator whose ACQUISITION function keeps ranking the same shallow
            candidates first.

    Both produce a small, highly repeatable measured set. Only the real pool size and
    its depth profile separate them, which is what this records.

    It is a write-only sink from the selector's point of view: `attach` sets
    `selector.candidate_trace`, the selector hands pools in, and nothing is read back
    until the run is over. It carries no outcome of any kind - identifiers and model
    scores only.
    """

    __slots__ = ("_budget", "_traces", "_pools")

    def __init__(self, budget: Optional[Any] = None) -> None:
        from veritas.diagnostics.trace import TraceBudget

        self._budget = budget or TraceBudget()
        self._traces: Dict[tuple, Any] = {}
        self._pools: Dict[tuple, List[str]] = {}

    def attach(self, selector: Any, *, method: str, seed: Any) -> Optional[Any]:
        """
        Gives one selector a trace sink for one campaign.

        Returns the sink, or None for a selector that has no scoring stage to trace
        (the baselines rank with a one-hot model and generate nothing). Setting the
        attribute is harmless either way: a selector that does not read it is
        unaffected, which is what keeps this additive.
        """
        from veritas.diagnostics.trace import BoundedCandidateTrace

        key = (str(method), str(seed))
        sink = self._traces.get(key)
        if sink is None:

            def _record_pool(round_index: int, ids: Sequence[str], k=key) -> None:
                self._pools.setdefault((k[0], k[1], int(round_index)), []).extend(ids)

            sink = BoundedCandidateTrace(
                campaign=f"{method}|seed{seed}",
                budget=self._budget,
                pool_sink=_record_pool,
            )
            self._traces[key] = sink
        try:
            selector.candidate_trace = sink
        except Exception:  # a selector with __slots__ and no such attribute
            return None
        return sink

    # -- read side, for the post-run hook only ---------------------------------------
    def pools(self) -> Dict[tuple, List[str]]:
        """{(method, seed, round): [every candidate notation proposed that round]}."""
        return {key: list(ids) for key, ids in self._pools.items()}

    def records(self) -> List[Any]:
        rows: List[Any] = []
        for key in sorted(self._traces, key=lambda k: (k[0], k[1])):
            rows.extend(self._traces[key].records)
        return rows

    def rounds_recorded(self) -> int:
        return sum(trace.rounds_recorded for trace in self._traces.values())


def emit_coverage_artifacts(
    run_dir: Path,
    outcomes: Optional[Sequence[Any]],
    *,
    eligible_ids: Optional[Sequence[str]],
    stamp: str,
    top1_predicate: Optional[Callable[[str], bool]] = None,
    run_id: str = "",
    thresholds: Optional[CollapseThresholds] = None,
    library: Optional[LibrarySpec] = None,
    proposed_pools: Optional[CandidatePoolRecorder] = None,
) -> List[CoverageDiagnosis]:
    """
    Computes generation coverage and per-depth returns for every arm of a finished
    run, writes them into the run folder, and warns on any arm whose candidate
    generation collapsed.

    `proposed_pools` is a `CandidatePoolRecorder` the run attached to its selectors.
    When present, coverage is computed from the TRUE generated pool. Without it,
    coverage falls back to the selection trace, which holds committed picks plus
    harness rejections and is therefore a LOWER BOUND on what was proposed - enough to
    show that an arm measures a tiny repeatable slice, but NOT enough to attribute
    that to generation rather than to acquisition. The `pool_source` column records
    which basis each row used, per arm, so no reader has to guess.

    `top1_predicate` is the evaluator's own `is_top1`. It is consulted here, after
    the run, over variants the campaigns already paid to measure - never handed to a
    selector. Without it the per-depth rows still carry proposed, selected and the
    mean/max revealed fitness of the selections; only the hit counts are zero.

    Returns one diagnosis per arm, or [] if anything went wrong.
    """
    try:
        return _emit(
            run_dir=Path(run_dir),
            outcomes=list(outcomes or ()),
            eligible_ids=list(eligible_ids or ()),
            stamp=stamp,
            top1_predicate=top1_predicate,
            run_id=run_id,
            thresholds=thresholds,
            library=library,
            proposed_pools=proposed_pools,
        )
    except Exception as exc:  # a diagnostic must never fail a benchmark run
        logger.error(
            "generation-coverage diagnostics skipped: %s: %s", type(exc).__name__, exc
        )
        return []


def _emit(
    *,
    run_dir: Path,
    outcomes: List[Any],
    eligible_ids: List[str],
    stamp: str,
    top1_predicate: Optional[Callable[[str], bool]],
    run_id: str,
    thresholds: Optional[CollapseThresholds],
    library: Optional[LibrarySpec],
    proposed_pools: Optional[CandidatePoolRecorder],
) -> List[CoverageDiagnosis]:
    rows = list(_trace_rows(outcomes))
    if not rows:
        logger.info("generation-coverage diagnostics skipped: the run produced no trace rows")
        return []
    if library is None:
        library = infer_library_from_variant_ids(
            eligible_ids or [getattr(row, "variant_id", "") for row in rows]
        )
    pools = proposed_pools.pools() if proposed_pools is not None else {}
    # Which arms had their real ballot captured. Reported per arm, never assumed.
    pool_methods = {key[0] for key in pools}

    # One accumulator per arm: the classification is a statement about an arm's
    # candidate generation, and merging arms would hide a collapsed one behind a
    # healthy one.
    accumulators: Dict[str, CoverageAccumulator] = {}
    ledger = DepthReturnLedger(library, run_id=run_id)
    per_round_proposed: Dict[tuple, List[str]] = {}
    per_round_measured: Dict[tuple, List[MeasuredSelection]] = {}

    for row in rows:
        status = getattr(row, "status", "")
        if status == STATUS_GIVEN:
            continue
        method = str(getattr(row, "method", ""))
        seed = getattr(row, "seed", "")
        round_index = int(getattr(row, "round_index", 0) or 0)
        variant_id = getattr(row, "variant_id", "")
        if not isinstance(variant_id, str) or not variant_id:
            continue
        key = (method, str(seed), round_index)
        per_round_proposed.setdefault(key, []).append(variant_id)
        if status == STATUS_SELECTED:
            per_round_measured.setdefault(key, [])
            revealed = getattr(row, "revealed_fitness", None)
            if revealed is not None:
                per_round_measured[key].append(
                    MeasuredSelection(
                        variant_id=variant_id,
                        revealed_fitness=float(revealed),
                        is_top1=bool(top1_predicate(variant_id)) if top1_predicate else False,
                    )
                )

    # Rounds that produced a pool but no surviving trace row still count: an arm whose
    # whole batch was rejected still generated a ballot.
    for key in pools:
        per_round_proposed.setdefault(key, [])

    for key in sorted(per_round_proposed):
        method, seed, round_index = key
        measured = per_round_measured.get(key, [])
        # The TRUE ballot when the sink captured it, otherwise the trace-derived lower
        # bound. The two are unioned rather than one replacing the other: the trace
        # also holds harness rejections, which are proposals the engine committed to.
        proposed = list(per_round_proposed[key])
        if key in pools:
            proposed.extend(pools[key])
        accumulator = accumulators.get(method)
        if accumulator is None:
            accumulator = CoverageAccumulator(library=library, eligible_ids=eligible_ids)
            accumulators[method] = accumulator
        accumulator.record_round(
            campaign=seed,
            round_index=round_index,
            proposed=proposed,
            selected=[selection.variant_id for selection in measured],
        )
        ledger.record_round(
            method=method,
            seed=seed,
            round_index=round_index,
            proposed=proposed,
            measured=measured,
        )

    diagnoses: List[CoverageDiagnosis] = []
    csv_rows: List[Sequence[Any]] = []
    arms: List[Dict[str, Any]] = []
    for method in sorted(accumulators):
        report = accumulators[method].seal()
        diagnosis = assert_generation_coverage(report, thresholds=thresholds, label=method)
        diagnoses.append(diagnosis)
        pool_source = "generated_pool" if method in pool_methods else "selection_trace_lower_bound"
        arms.append(
            {
                "method": method,
                "pool_source": pool_source,
                "diagnosis": diagnosis.to_dict(),
                "coverage": report.to_dict(),
            }
        )
        limits = diagnosis.thresholds
        for round_row in report.rounds:
            csv_rows.append(
                (
                    stamp,
                    method,
                    round_row.round_index,
                    report.campaigns,
                    round_row.proposed_distinct,
                    round_row.proposed_total,
                    round_row.cumulative_distinct,
                    f"{round_row.cumulative_share:.9g}",
                    round_row.selected_distinct,
                    round_row.cumulative_selected_distinct,
                    f"{round_row.cumulative_mean_proposed_depth:.9g}",
                    json.dumps({str(k): v for k, v in round_row.proposed_depth_histogram.items()}),
                    json.dumps({str(k): v for k, v in round_row.selected_depth_histogram.items()}),
                    pool_source,
                    # Label last among the arm columns, continuous values first: the
                    # label is derived, these are the data a revised threshold needs.
                    diagnosis.classification,
                    diagnosis.union_distinct,
                    f"{diagnosis.union_share:.9g}",
                    f"{diagnosis.mean_proposed_depth:.9g}",
                    f"{diagnosis.depth_reach:.9g}",
                    f"{diagnosis.deep_enrichment:.9g}",
                    diagnosis.deep_proposed,
                    f"{diagnosis.deep_universe_share:.9g}",
                    limits.min_union_share,
                    limits.min_depth_reach,
                    limits.min_deep_enrichment,
                    limits.min_proposals_to_judge,
                )
            )

    run_dir.mkdir(parents=True, exist_ok=True)
    with open(run_dir / "generation_coverage.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(COVERAGE_CSV_COLUMNS)
        writer.writerows(csv_rows)

    table = ledger.seal()
    write_depth_returns_csv(run_dir / "depth_returns.csv", table.rows, stamp)

    trace_records = proposed_pools.records() if proposed_pools is not None else []
    if trace_records:
        write_candidate_trace_csv(run_dir / "candidate_trace.csv", trace_records, stamp)

    (run_dir / "generation_coverage.json").write_text(
        json.dumps(
            {
                "provenance": stamp,
                "run_id": run_id,
                "library": dict(library.to_dict()),
                "top1_membership_available": top1_predicate is not None,
                "generated_pool_captured_for": sorted(pool_methods),
                "candidate_trace_rows": len(trace_records),
                "arms": arms,
                "depth_returns": {
                    "methods": list(table.methods()),
                    "mean_selected_depth": {
                        m: table.mean_selected_depth(m) for m in table.methods()
                    },
                    "mean_proposed_depth": {
                        m: table.mean_proposed_depth(m) for m in table.methods()
                    },
                    "new_hits_per_round": {
                        m: {str(k): v for k, v in table.new_hits_per_round(m).items()}
                        for m in table.methods()
                    },
                    "shell_returns": {m: table.shell_returns(m) for m in table.methods()},
                    "rows_written": len(table.rows),
                },
                "note": (
                    "Coverage counts DISTINCT CANONICAL GENOTYPES the arm proposed, "
                    "excluding the shared starting sample. For an arm listed in "
                    "generated_pool_captured_for the basis is the engine's real "
                    "candidate pool; for any other arm it is the selection trace "
                    "(committed picks plus harness rejections), which is a LOWER BOUND "
                    "on the ballot and therefore cannot by itself distinguish narrow "
                    "generation from acquisition repeatedly ranking the same shallow "
                    "candidates. Each row states its basis in pool_source. Per-depth "
                    "returns are derived after the run from outcomes the campaigns "
                    "already paid to measure; top-1% membership is the evaluator's and "
                    "is never exposed to a selector."
                ),
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    return diagnoses
