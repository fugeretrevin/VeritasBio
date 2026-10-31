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

logger = logging.getLogger("veritas.diagnostics.benchmark_hook")

# Statuses used by `selection_trace.csv`. `initial` is the shared starting sample,
# which the harness hands to every arm: it is not something candidate generation
# produced, so counting it would flatter every arm by the same 25 variants.
STATUS_SELECTED = "accepted"
STATUS_REJECTED = "rejected"
STATUS_GIVEN = "initial"

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
    "proposed_depth_histogram",
    "selected_depth_histogram",
)


def _trace_rows(outcomes: Iterable[Any]) -> Iterable[Any]:
    for outcome in outcomes or ():
        for row in getattr(outcome, "trace", ()) or ():
            yield row


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
) -> List[CoverageDiagnosis]:
    """
    Computes generation coverage and per-depth returns for every arm of a finished
    run, writes them into the run folder, and warns on any arm whose candidate
    generation collapsed.

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
) -> List[CoverageDiagnosis]:
    rows = list(_trace_rows(outcomes))
    if not rows:
        logger.info("generation-coverage diagnostics skipped: the run produced no trace rows")
        return []
    if library is None:
        library = infer_library_from_variant_ids(
            eligible_ids or [getattr(row, "variant_id", "") for row in rows]
        )

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

    for (method, seed, round_index), proposed in sorted(per_round_proposed.items()):
        measured = per_round_measured.get((method, seed, round_index), [])
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
        arms.append(
            {
                "method": method,
                "diagnosis": diagnosis.to_dict(),
                "coverage": report.to_dict(),
            }
        )
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
                    json.dumps({str(k): v for k, v in round_row.proposed_depth_histogram.items()}),
                    json.dumps({str(k): v for k, v in round_row.selected_depth_histogram.items()}),
                )
            )

    run_dir.mkdir(parents=True, exist_ok=True)
    with open(run_dir / "generation_coverage.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(COVERAGE_CSV_COLUMNS)
        writer.writerows(csv_rows)

    table = ledger.seal()
    write_depth_returns_csv(run_dir / "depth_returns.csv", table.rows, stamp)

    (run_dir / "generation_coverage.json").write_text(
        json.dumps(
            {
                "provenance": stamp,
                "run_id": run_id,
                "library": dict(library.to_dict()),
                "top1_membership_available": top1_predicate is not None,
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
                    "Coverage is computed from candidates the arms PROPOSED (accepted "
                    "plus harness-rejected), excluding the shared starting sample. "
                    "Per-depth returns are derived after the run from outcomes the "
                    "campaigns already paid to measure; top-1% membership is the "
                    "evaluator's and is never exposed to a selector."
                ),
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    return diagnoses
