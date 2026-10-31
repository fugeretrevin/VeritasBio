"""
veritas/diagnostics/retrospective.py - the ONLY place a coverage metric may reference
a true experimental outcome.

Why this is a separate module with its own gate
-----------------------------------------------
The single most useful number in the TrpB post-mortem was "the best fitness available
anywhere in that 294-variant reachable set is exactly 0.7036". It proved acquisition
was choosing correctly and that candidate generation was at fault.

It is also the most dangerous number in the codebase. A selector that could read it
would be reading the answer key. So it is not merely documented as retrospective, it
is made unreachable during a run:

  1. `CoverageAccumulator` has `__slots__` and accepts bare identifier strings only,
     so a live run cannot hold or record a fitness at all.
  2. `analyse_reachable_set` refuses an UNSEALED accumulator. Sealing is irreversible
     and blocks every further `record_round`, so by the time a fitness has been read,
     the object it was read against can no longer affect a selection that has yet to
     happen.
  3. It also refuses to run without a `RetrospectiveFitnessAccess` token carrying a
     written reason. A truthy flag or a string is not accepted - the caller has to
     construct the type, which is a line of code nobody writes by accident inside a
     selection loop.
  4. The result is a `RetrospectiveCoverageReport`, a type DISJOINT from
     `CoverageReport`. `classify_coverage` rejects it with a TypeError, so an
     outcome-bearing report cannot be mistaken for a live one anywhere downstream.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Mapping, Optional, Tuple

from veritas.diagnostics.genotype import GenotypeParseError, LibrarySpec


class RetrospectiveAccessError(RuntimeError):
    """Raised when hidden experimental outcomes were requested outside a sealed, post-run context."""


@dataclass(frozen=True, slots=True)
class RetrospectiveFitnessAccess:
    """
    An explicit, reasoned permit to read the evaluator's truth table for ANALYSIS.

    Constructing one is a deliberate act and the reason is recorded in the result, so
    a run folder states why the answer key was consulted.
    """

    reason: str

    def __post_init__(self) -> None:
        if not isinstance(self.reason, str) or not self.reason.strip():
            raise ValueError(
                "RetrospectiveFitnessAccess requires a written reason, e.g. "
                "'attributing the 0.7036 ceiling to candidate generation'"
            )


@dataclass(frozen=True, slots=True)
class RetrospectiveCoverageReport:
    """
    Post-run coverage metrics that reference measured fitness.

    Deliberately NOT a subclass of, and not substitutable for, `CoverageReport`. Its
    presence in a data structure means an outcome has been read, which is why the type
    is separate and the live classifier refuses it.
    """

    reason: str
    reachable_distinct: int
    scored_distinct: int
    unscored_distinct: int
    best_reachable_fitness: float
    best_reachable_variant: str
    mean_reachable_fitness: float
    best_eligible_fitness: float
    best_eligible_variant: str
    reachable_fitness_ceiling_gap: float
    best_selected_fitness: Optional[float]
    selected_distinct: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def analyse_reachable_set(
    accumulator: Any,
    *,
    fitness_lookup: Mapping[str, float],
    access: Optional[RetrospectiveFitnessAccess],
) -> RetrospectiveCoverageReport:
    """
    Scores an arm's reachable set against the complete truth table, AFTER the run.

    `accumulator` must be a sealed `CoverageAccumulator`. The headline output is
    `best_reachable_fitness`: the best outcome that was available anywhere on the
    ballot the generator ever produced. When that number sits far below the global
    optimum - 0.7036 against 1.0 on TrpB - the campaign's ceiling was set by candidate
    generation, not by acquisition, and no amount of better ranking could have fixed it.
    """
    if not isinstance(access, RetrospectiveFitnessAccess):
        raise RetrospectiveAccessError(
            "reading measured outcomes requires a RetrospectiveFitnessAccess token with "
            f"a written reason; got {type(access).__name__}. Coverage metrics that "
            "reference true fitness are retrospective only and must never be readable "
            "during selection."
        )
    if not getattr(accumulator, "sealed", False):
        raise RetrospectiveAccessError(
            "the CoverageAccumulator is not sealed. Seal it first: a sealed accumulator "
            "can record no further round, which is what guarantees a number derived from "
            "the truth table cannot flow back into a selection."
        )

    reachable: Tuple[str, ...] = tuple(sorted(accumulator.reachable_ids))
    selected: Tuple[str, ...] = tuple(sorted(accumulator.selected_ids))

    scored = [(v, float(fitness_lookup[v])) for v in reachable if v in fitness_lookup]
    if not scored:
        raise RetrospectiveAccessError(
            f"none of the {len(reachable)} reachable variant(s) appear in the supplied "
            "truth table; check that the identifiers are canonical."
        )
    best_variant, best_value = max(scored, key=lambda pair: pair[1])
    mean_value = sum(value for _, value in scored) / len(scored)

    eligible_best_variant, eligible_best_value = max(
        fitness_lookup.items(), key=lambda pair: float(pair[1])
    )
    selected_scored = [float(fitness_lookup[v]) for v in selected if v in fitness_lookup]

    return RetrospectiveCoverageReport(
        reason=access.reason,
        reachable_distinct=len(reachable),
        scored_distinct=len(scored),
        unscored_distinct=len(reachable) - len(scored),
        best_reachable_fitness=best_value,
        best_reachable_variant=best_variant,
        mean_reachable_fitness=mean_value,
        best_eligible_fitness=float(eligible_best_value),
        best_eligible_variant=eligible_best_variant,
        reachable_fitness_ceiling_gap=float(eligible_best_value) - best_value,
        best_selected_fitness=max(selected_scored) if selected_scored else None,
        selected_distinct=len(selected),
    )


@dataclass(frozen=True, slots=True)
class DepthShellProfile:
    """
    How the top-1% of a landscape is distributed over its mutation-depth shells.

    This is a property of the LANDSCAPE, not of any run: it says where the good
    variants are, which is precisely the information a selection policy must not be
    allowed to read. It is the reference against which a run's own per-depth returns
    (`veritas.diagnostics.depth_returns`) are judged - e.g. TrpB's shallow shell
    (depth <= 2) being several times richer in top-1% variants than its deep shell,
    while GB1's shells are much closer together, is the shape of the GB1-vs-TrpB
    reversal.
    """

    reason: str
    threshold: float
    n_scored: int
    n_hits: int
    per_depth: Mapping[int, Mapping[str, float]]
    shallow_max_depth: int
    shallow: Mapping[str, float]
    deep: Mapping[str, float]
    unparseable: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reason": self.reason,
            "threshold": self.threshold,
            "n_scored": self.n_scored,
            "n_hits": self.n_hits,
            "shallow_max_depth": self.shallow_max_depth,
            "per_depth": {str(k): dict(v) for k, v in sorted(self.per_depth.items())},
            "shallow": dict(self.shallow),
            "deep": dict(self.deep),
            "unparseable": self.unparseable,
        }


def depth_shell_profile(
    library: LibrarySpec,
    *,
    fitness_lookup: Mapping[str, float],
    access: Optional[RetrospectiveFitnessAccess],
    top_percentile: float = 99.0,
    shallow_max_depth: int = 2,
) -> DepthShellProfile:
    """
    Top-1% hit rate per mutation-depth shell over a complete landscape.

    RETROSPECTIVE ONLY, and gated on the same access token as everything else here:
    the threshold is defined over variants nobody measured, so this is answer-key
    information. It exists to explain a result after the fact, never to steer one.
    """
    if not isinstance(access, RetrospectiveFitnessAccess):
        raise RetrospectiveAccessError(
            "depth_shell_profile reads the complete landscape and requires a "
            "RetrospectiveFitnessAccess token with a written reason; got "
            f"{type(access).__name__}."
        )
    if not 0.0 < top_percentile < 100.0:
        raise ValueError(f"top_percentile must be in (0, 100), got {top_percentile}")

    scored: Dict[int, list] = {}
    unparseable = 0
    values: list = []
    for variant_id, value in fitness_lookup.items():
        try:
            depth = library.mutation_depth(variant_id)
        except GenotypeParseError:
            unparseable += 1
            continue
        numeric = float(value)
        scored.setdefault(depth, []).append(numeric)
        values.append(numeric)

    if not values:
        raise RetrospectiveAccessError("the supplied truth table scored no variant")

    values.sort()
    # Linear-interpolated percentile, matching numpy's default so the threshold agrees
    # with the benchmark's evaluator.
    position = (len(values) - 1) * (top_percentile / 100.0)
    lower, upper = int(position), min(int(position) + 1, len(values) - 1)
    threshold = values[lower] + (position - lower) * (values[upper] - values[lower])

    per_depth: Dict[int, Dict[str, float]] = {}
    total_hits = 0
    for depth in sorted(scored):
        shell = scored[depth]
        hits = sum(1 for value in shell if value >= threshold)
        total_hits += hits
        per_depth[depth] = {
            "n": float(len(shell)),
            "hits": float(hits),
            "hit_rate": hits / len(shell),
            "max": max(shell),
        }

    def _shell(predicate) -> Dict[str, float]:
        n = sum(cell["n"] for depth, cell in per_depth.items() if predicate(depth))
        hits = sum(cell["hits"] for depth, cell in per_depth.items() if predicate(depth))
        return {"n": n, "hits": hits, "hit_rate": (hits / n) if n else 0.0}

    return DepthShellProfile(
        reason=access.reason,
        threshold=threshold,
        n_scored=len(values),
        n_hits=total_hits,
        per_depth=per_depth,
        shallow_max_depth=shallow_max_depth,
        shallow=_shell(lambda d: d <= shallow_max_depth),
        deep=_shell(lambda d: d > shallow_max_depth),
        unparseable=unparseable,
    )
