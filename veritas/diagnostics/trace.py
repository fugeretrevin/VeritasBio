"""
veritas/diagnostics/trace.py - the enriched candidate-selection trace, bounded from
the start.

Why bounded and not complete
----------------------------
A full per-candidate trace is not viable. GB1 generates on the order of 149,000
candidates per round, over 4 rounds, 48 seeds and several arms: a complete trace is
tens of millions of rows per arm and the write dominates the run. So the bound is
part of the design rather than something bolted on when a disk fills up:

  * EVERY selected candidate is recorded. The batch that was actually committed is
    always fully auditable.
  * The top-N ranked candidates are recorded, so "what did it nearly pick" is
    answerable.
  * A seeded random sample of the rejected candidates is recorded, which is what
    makes the rejected population characterisable (its depth profile, its generation
    source mix) without storing it.

Per round that is `len(selected) + top_n + rejected_sample` rows - a few dozen -
regardless of pool size.

Coverage accounting is separate on purpose. `veritas.diagnostics.coverage` gets the
FULL pool of identifiers (a set insert per candidate, no rows), so the 0.18%-coverage
diagnosis is exact even though the trace is sampled. Pass the accumulator as
`coverage=` and one `record_round` call feeds both.

Outcome discipline
------------------
`predicted_fitness` is a model output. There is no field for a revealed outcome, and
`veritas.diagnostics.coverage` / `.retrospective` enforce the same separation, so a
trace row can never become a channel for the answer key.

Backward compatibility
----------------------
`LEGACY_SELECTION_TRACE_COLUMNS` pins the nine columns of the existing
`selection_trace.csv` written by `benchmark_panoptes_isolated.py`. Today's analysis
scripts read that file, so the enriched trace is written ALONGSIDE it, to
`candidate_trace.csv`, and the legacy schema is asserted in the test suite rather
than trusted.
"""
from __future__ import annotations

import csv
import logging
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Sequence, Tuple

# The existing trace schema. Pinned so a change is a deliberate, reviewed act.
LEGACY_SELECTION_TRACE_COLUMNS: Tuple[str, ...] = (
    "provenance",
    "method",
    "seed",
    "round",
    "rank",
    "variant_id",
    "status",
    "reason",
    "revealed_fitness",
)

logger = logging.getLogger("veritas.diagnostics.trace")

STATUS_SELECTED = "selected"
STATUS_RANKED = "ranked"
STATUS_REJECTED = "rejected"

SOURCE_COMBINATORIAL = "combinatorial"
SOURCE_GFLOWNET = "gflownet"
SOURCE_UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class TraceBudget:
    """
    The bound. `top_n` ranked rows and `rejected_sample` sampled rejections per round,
    on top of every selected candidate. `seed` makes the sample reproducible, so two
    analyses of the same run agree about which rejections they are looking at.
    """

    top_n: int = 50
    rejected_sample: int = 25
    seed: int = 0

    def __post_init__(self) -> None:
        if self.top_n < 0 or self.rejected_sample < 0:
            raise ValueError("trace budget sizes must be non-negative")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CandidateTraceRecord:
    """
    One candidate's scoring detail, as the selector saw it.

    Fields map onto the stages that can each independently ruin a batch:
      generation_source       which generator put it on the ballot (the stage that
                              failed on TrpB)
      predicted_fitness       surrogate posterior mean mu(x)
      posterior_uncertainty   surrogate posterior sigma(x)
      prior_contribution      the biophysical / structural prior m(x) for THIS
                              candidate: the calibrated total the optimizer consumed,
                              including its mechanistic bonus
      prior_structural        the structural-surrogate term alone, on the RAW
                              pre-calibration scale. Reported separately because the
                              calibration is a tanh, so the calibrated total cannot be
                              decomposed after the fact
      prior_raw_total         raw electrostatics + catalytic penalty + structural,
                              pre-calibration, so `prior_structural` has a denominator
      prior_structural_available
                              False when no structure was supplied or it could not be
                              scored - the difference between "the structural term was
                              zero" and "there was no structural term"
      acquisition_score       EI, or EI + the Thompson term under hybrid acquisition
      pareto_rank             multi-objective non-dominated rank (0 = non-dominated)
      pareto_contribution     the multiplier that rank applied to acquisition value
      quality_score           acquisition x pareto x auxiliary weights, pre-diversity
      final_rank              1-based rank by quality_score over the whole pool
      selection_order         1-based position in the committed batch, 0 if unselected
      diversity_adjustment    final_rank - selection_order for a selected candidate:
                              positive means the diversity stage promoted it past
                              better-scoring peers, which is the k-DPP effect made
                              visible instead of inferred
    """

    campaign: str
    round_index: int
    variant_id: str
    generation_source: str
    predicted_fitness: float
    posterior_uncertainty: float
    prior_contribution: float
    prior_structural: float
    prior_raw_total: float
    prior_structural_available: bool
    acquisition_score: float
    pareto_rank: int
    pareto_contribution: float
    quality_score: float
    diversity_adjustment: float
    final_rank: int
    selection_order: int
    status: str
    rejection_reason: str

    @classmethod
    def field_names(cls) -> Tuple[str, ...]:
        return tuple(cls.__dataclass_fields__)

    def to_row(self) -> Tuple[Any, ...]:
        return tuple(getattr(self, name) for name in self.field_names())


class BoundedCandidateTrace:
    """
    The sink `workers.math_worker.select_next_micro_batch` writes a scored pool into.

    Duck-typed on purpose: math_worker calls `record_round(**arrays)` and never
    imports this package, so the production selector keeps no dependency on the
    diagnostics and the two can be merged independently.
    """

    __slots__ = (
        "_campaign",
        "_budget",
        "_coverage",
        "_pool_sink",
        "_records",
        "_rounds_recorded",
    )

    def __init__(
        self,
        campaign: str = "",
        budget: Optional[TraceBudget] = None,
        coverage: Optional[Any] = None,
        pool_sink: Optional[Callable[[int, Sequence[str]], None]] = None,
    ) -> None:
        self._campaign = str(campaign)
        self._budget = budget or TraceBudget()
        self._coverage = coverage
        # The UNBOUNDED pool of identifiers, for coverage accounting. A set insert per
        # candidate, no rows: this is what separates "the generator is narrow" from
        # "acquisition keeps picking the same shallow candidates from a broad pool",
        # and the bounded rows below cannot answer that.
        self._pool_sink = pool_sink
        self._records: list[CandidateTraceRecord] = []
        self._rounds_recorded = 0

    @property
    def records(self) -> Tuple[CandidateTraceRecord, ...]:
        return tuple(self._records)

    @property
    def budget(self) -> TraceBudget:
        return self._budget

    @property
    def rounds_recorded(self) -> int:
        return self._rounds_recorded

    def note_pool(
        self,
        candidates: Sequence[str],
        *,
        round_index: Optional[int] = None,
        source: str = "",
    ) -> None:
        """
        Records a ballot for COVERAGE ONLY, with no per-candidate rows.

        For a generator that does not go through the scored pipeline and so has no
        mu(x)/sigma(x)/acquisition value to report - a greedy Hamming-1 neighbour
        enumeration, a conditional scan. Without this, such a stage is invisible to
        coverage and its arm's ballot is silently under-counted, which is exactly the
        error that made a committed-picks-only coverage number look like evidence
        about generation.
        """
        if not candidates:
            return
        index = self._rounds_recorded + 1 if round_index is None else int(round_index)
        ids = [str(c) for c in candidates]
        if source:
            logger.debug("ballot from %s: %d candidate(s)", source, len(ids))
        if self._coverage is not None:
            self._coverage.record_round(
                campaign=self._campaign, round_index=index, proposed=ids, selected=[]
            )
        if self._pool_sink is not None:
            self._pool_sink(index, ids)
        self._rounds_recorded += 1

    def record_round(
        self,
        *,
        round_index: Optional[int] = None,
        candidates: Sequence[str],
        sources: Sequence[str],
        predicted_fitness: Sequence[float],
        posterior_uncertainty: Sequence[float],
        prior_contribution: Sequence[float],
        acquisition_score: Sequence[float],
        pareto_rank: Sequence[int],
        pareto_contribution: Sequence[float],
        quality_score: Sequence[float],
        selected_indices: Sequence[int],
        rejection_reasons: Optional[Mapping[int, str]] = None,
        prior_structural: Optional[Sequence[Optional[float]]] = None,
        prior_raw_total: Optional[Sequence[Optional[float]]] = None,
        prior_structural_available: Optional[Sequence[bool]] = None,
    ) -> None:
        """
        Records one scored candidate pool, down to the budget.

        Every per-candidate sequence must have the same length as `candidates`:
        a ragged input is refused rather than zipped, because a silently truncated
        zip is how a trace starts describing the wrong candidate.

        `round_index=None` lets the sink number its own rounds, so a caller that does
        not track a round counter (the production selector does not) still produces an
        ordered trace.
        """
        if round_index is None:
            round_index = self._rounds_recorded + 1
        n = len(candidates)
        columns = {
            "sources": sources,
            "predicted_fitness": predicted_fitness,
            "posterior_uncertainty": posterior_uncertainty,
            "prior_contribution": prior_contribution,
            "acquisition_score": acquisition_score,
            "pareto_rank": pareto_rank,
            "pareto_contribution": pareto_contribution,
            "quality_score": quality_score,
        }
        for name, optional in (
            ("prior_structural", prior_structural),
            ("prior_raw_total", prior_raw_total),
            ("prior_structural_available", prior_structural_available),
        ):
            if optional is not None:
                columns[name] = optional
        ragged = {name: len(values) for name, values in columns.items() if len(values) != n}
        if ragged:
            raise ValueError(
                f"candidate trace columns must all be {n} long (one per candidate); "
                f"got {ragged}. Refusing to zip unequal-length columns."
            )

        selected = [int(i) for i in selected_indices]
        bad = [i for i in selected if not 0 <= i < n]
        if bad:
            raise ValueError(f"selected_indices {bad} are outside the candidate pool of {n}")

        # Coverage sees the WHOLE pool; only the rows below are bounded.
        if self._coverage is not None:
            self._coverage.record_round(
                campaign=self._campaign,
                round_index=round_index,
                proposed=list(candidates),
                selected=[candidates[i] for i in selected],
            )
        if self._pool_sink is not None:
            self._pool_sink(int(round_index), [str(c) for c in candidates])

        order = sorted(range(n), key=lambda i: (-float(quality_score[i]), i))
        rank_of = {index: position + 1 for position, index in enumerate(order)}
        selection_order_of = {index: position + 1 for position, index in enumerate(selected)}

        selected_set = set(selected)
        ranked = [i for i in order[: self._budget.top_n] if i not in selected_set]
        ranked_set = set(ranked)
        remainder = [i for i in range(n) if i not in selected_set and i not in ranked_set]
        sample_size = min(self._budget.rejected_sample, len(remainder))
        # Seeded from a STRING, not from hash(): str hashing is salted per process, so a
        # hash-derived seed would give a different sample on every run.
        rng = random.Random(f"{self._budget.seed}|{self._campaign}|{round_index}")
        sampled = sorted(rng.sample(remainder, sample_size)) if sample_size else []

        reasons = dict(rejection_reasons or {})
        for index in sorted(selected_set) + ranked + sampled:
            if index in selected_set:
                status = STATUS_SELECTED
                reason = ""
                selection_order = selection_order_of[index]
                adjustment = float(rank_of[index] - selection_order)
            else:
                status = STATUS_RANKED if index in set(ranked) else STATUS_REJECTED
                selection_order = 0
                adjustment = 0.0
                reason = reasons.get(index) or (
                    f"outranked_by_batch_selection (rank {rank_of[index]} of {n})"
                    if status == STATUS_RANKED
                    else f"not_selected_outside_top_{self._budget.top_n}"
                )
            self._records.append(
                CandidateTraceRecord(
                    campaign=self._campaign,
                    round_index=int(round_index),
                    variant_id=str(candidates[index]),
                    generation_source=str(sources[index] or SOURCE_UNKNOWN),
                    predicted_fitness=float(predicted_fitness[index]),
                    posterior_uncertainty=float(posterior_uncertainty[index]),
                    prior_contribution=float(prior_contribution[index]),
                    prior_structural=(
                        0.0
                        if prior_structural is None or prior_structural[index] is None
                        else float(prior_structural[index])
                    ),
                    prior_raw_total=(
                        0.0
                        if prior_raw_total is None or prior_raw_total[index] is None
                        else float(prior_raw_total[index])
                    ),
                    prior_structural_available=(
                        False
                        if prior_structural_available is None
                        else bool(prior_structural_available[index])
                    ),
                    acquisition_score=float(acquisition_score[index]),
                    pareto_rank=int(pareto_rank[index]),
                    pareto_contribution=float(pareto_contribution[index]),
                    quality_score=float(quality_score[index]),
                    diversity_adjustment=adjustment,
                    final_rank=int(rank_of[index]),
                    selection_order=selection_order,
                    status=status,
                    rejection_reason=reason,
                )
            )
        self._rounds_recorded += 1


def write_candidate_trace_csv(
    path: Path,
    records: Iterable[CandidateTraceRecord],
    stamp: str,
) -> int:
    """
    Writes the enriched trace to its own file, alongside the legacy
    `selection_trace.csv` rather than instead of it. Returns the row count.
    """
    header = ("provenance",) + CandidateTraceRecord.field_names()
    rows = 0
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for record in records:
            writer.writerow((stamp,) + record.to_row())
            rows += 1
    return rows
