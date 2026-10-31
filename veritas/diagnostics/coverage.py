"""
veritas/diagnostics/coverage.py - candidate-GENERATION coverage, and the collapse
classifier that would have caught the TrpB `frontier` failure.

The failure this was built for
------------------------------
On the published four-site TrpB landscape, the `frontier` configuration with the
`dkl` and `gflownet` components ablated terminated all 96 campaigns at exactly
fitness 0.7036, against a global optimum of 1.0.

Two guards passed it:
  * a terminal-variance check - all campaigns ending at one value looks like a
    deterministic policy, which is legitimate;
  * a data-invariance check - the policy is NOT data-blind. Its round-1 batches
    genuinely vary with the starting sample (Jaccard 0.354, not 1.0).

What was actually wrong was upstream of acquisition. Across those 96 campaigns the
arm only ever reached 294 DISTINCT variants, 0.18% of the 159,999-variant eligible
set, at a mean mutation depth of 2.56, with just 45 distinct depth-4 variants -
while depth 4 holds 130,321 variants, 81% of the landscape, and contains the global
optimum. Plain ridge regression reached 4,028 distinct variants (2.52%) at mean
depth 3.51, 14x more of the landscape. The best fitness anywhere inside that
294-variant reachable set is exactly 0.7036: acquisition was choosing correctly
from an impoverished ballot. The fault is candidate generation alone.

So the metric is coverage of the eligible universe by PROPOSED candidates (not just
selected ones), plus the mutation-depth profile of what was proposed and of what was
selected, per round, cumulatively, and as a union over a whole run.

What this module deliberately does NOT do
-----------------------------------------
It never rejects a run. The review this comes from is explicit that legitimate
deterministic convergence must not be failed, so `assert_generation_coverage` warns
and classifies - healthy / narrow-generation / depth-truncated / collapsed - and
hands back the numbers that justify the label.

It also never touches a measured outcome. Nothing recorded here is anything but a
sequence identifier, and the live accumulator cannot be given a fitness attribute
(it has `__slots__`). Any metric that needs true fitness - "the best fitness
available in the reachable set is 0.7036" - lives in
`veritas.diagnostics.retrospective`, which can only be reached once the accumulator
has been sealed and can no longer influence a selection.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import asdict, dataclass
from typing import Any, Dict, Hashable, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from veritas.diagnostics.genotype import GenotypeParseError, LibrarySpec

logger = logging.getLogger("veritas.diagnostics.coverage")


class GenerationCollapseWarning(UserWarning):
    """Emitted when candidate generation covered too little of the eligible universe."""


# ==============================================================================================
# Reports
# ==============================================================================================
@dataclass(frozen=True, slots=True)
class RoundCoverage:
    """
    One round's generation coverage, aggregated over every campaign of the arm.

    `proposed_distinct` is the number of distinct candidates this round put on the
    ballot anywhere in the run; `cumulative_distinct` is the union over rounds up to
    and including this one. `proposed_total` keeps the un-deduplicated proposal count,
    so re-proposing the same handful of variants is visible as a ratio.
    """

    round_index: int
    proposed_distinct: int
    proposed_total: int
    cumulative_distinct: int
    cumulative_share: float
    selected_distinct: int
    selected_total: int
    cumulative_selected_distinct: int
    cumulative_selected_share: float
    proposed_depth_histogram: Mapping[int, int]
    selected_depth_histogram: Mapping[int, int]

    def to_dict(self) -> Dict[str, Any]:
        payload = asdict(self)
        payload["proposed_depth_histogram"] = {
            str(k): v for k, v in sorted(self.proposed_depth_histogram.items())
        }
        payload["selected_depth_histogram"] = {
            str(k): v for k, v in sorted(self.selected_depth_histogram.items())
        }
        return payload


@dataclass(frozen=True, slots=True)
class CoverageReport:
    """
    The sealed, outcome-free summary of one arm's candidate generation.

    Every field is derived from sequence identifiers and the library definition. There
    is deliberately no field that could hold a measured fitness, and the dataclass is
    frozen with slots, so one cannot be attached later.
    """

    library: LibrarySpec
    library_size: int
    max_library_depth: int
    universe_depth_histogram: Mapping[int, int]
    campaigns: int
    rounds: Tuple[RoundCoverage, ...]
    union_proposed_distinct: int
    union_proposed_total: int
    union_proposed_share: float
    union_proposed_depth_histogram: Mapping[int, int]
    union_selected_distinct: int
    union_selected_total: int
    union_selected_share: float
    union_selected_depth_histogram: Mapping[int, int]
    mean_proposed_depth: float
    mean_selected_depth: float
    unparseable_proposals: int
    off_universe_proposals: int

    @classmethod
    def field_names(cls) -> Tuple[str, ...]:
        return tuple(cls.__dataclass_fields__)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "library": dict(self.library.to_dict()),
            "library_size": self.library_size,
            "max_library_depth": self.max_library_depth,
            "universe_depth_histogram": {
                str(k): v for k, v in sorted(self.universe_depth_histogram.items())
            },
            "campaigns": self.campaigns,
            "union_proposed_distinct": self.union_proposed_distinct,
            "union_proposed_total": self.union_proposed_total,
            "union_proposed_share": self.union_proposed_share,
            "union_proposed_depth_histogram": {
                str(k): v for k, v in sorted(self.union_proposed_depth_histogram.items())
            },
            "union_selected_distinct": self.union_selected_distinct,
            "union_selected_total": self.union_selected_total,
            "union_selected_share": self.union_selected_share,
            "union_selected_depth_histogram": {
                str(k): v for k, v in sorted(self.union_selected_depth_histogram.items())
            },
            "mean_proposed_depth": self.mean_proposed_depth,
            "mean_selected_depth": self.mean_selected_depth,
            "unparseable_proposals": self.unparseable_proposals,
            "off_universe_proposals": self.off_universe_proposals,
            "rounds": [row.to_dict() for row in self.rounds],
        }


# ==============================================================================================
# The live accumulator
# ==============================================================================================
class CoverageAccumulator:
    """
    Append-only, outcome-free tally of what candidate generation put on the ballot.

    Usage during a run:
        accumulator.record_round(campaign=seed, round_index=r, proposed=[...], selected=[...])
    and once, at the end:
        report = accumulator.seal()

    Three properties are structural, not conventional:

      * `__slots__` - no measured outcome can be attached to a live accumulator;
      * `record_round` accepts bare identifier STRINGS only, so an (id, fitness) pair
        or a dict carrying a revealed score is a TypeError rather than something that
        quietly works;
      * sealing is final. Any metric that needs true fitness is reachable only from a
        sealed accumulator (see `veritas.diagnostics.retrospective`), which therefore
        cannot feed anything back into a round that has yet to happen.

    Scale: only identifier SETS are kept, never per-candidate rows, so a 149,000
    candidate pool costs a set insert per candidate and nothing per round. The
    per-candidate detail is the job of `veritas.diagnostics.trace`, which is bounded.
    """

    __slots__ = (
        "_library",
        "_eligible",
        "_universe_depth_histogram",
        "_depth_cache",
        "_campaigns",
        "_round_proposed",
        "_round_selected",
        "_round_proposed_total",
        "_round_selected_total",
        "_unparseable",
        "_off_universe",
        "_sealed",
    )

    def __init__(
        self,
        library: LibrarySpec,
        eligible_ids: Iterable[str],
        *,
        universe_depth_histogram: Optional[Mapping[int, int]] = None,
    ) -> None:
        if not isinstance(library, LibrarySpec):
            raise TypeError(f"library must be a LibrarySpec, got {type(library).__name__}")
        self._library = library
        self._depth_cache: Dict[str, int] = {}
        self._eligible: Set[str] = set()
        histogram: Dict[int, int] = {}
        for variant_id in eligible_ids:
            self._require_identifier(variant_id, "eligible_ids")
            self._eligible.add(variant_id)
            try:
                depth = self._depth(variant_id)
            except GenotypeParseError:
                continue
            histogram[depth] = histogram.get(depth, 0) + 1
        self._universe_depth_histogram: Dict[int, int] = (
            dict(universe_depth_histogram) if universe_depth_histogram is not None else histogram
        )
        self._campaigns: Set[Hashable] = set()
        self._round_proposed: Dict[int, Set[str]] = {}
        self._round_selected: Dict[int, Set[str]] = {}
        self._round_proposed_total: Dict[int, int] = {}
        self._round_selected_total: Dict[int, int] = {}
        self._unparseable = 0
        self._off_universe = 0
        self._sealed = False

    # -- guards ----------------------------------------------------------------------
    @staticmethod
    def _require_identifier(value: Any, where: str) -> None:
        if not isinstance(value, str):
            raise TypeError(
                f"{where} must contain bare variant-identifier strings; got "
                f"{type(value).__name__}. Coverage accounting is outcome-free by "
                "construction: a measured fitness cannot be smuggled in alongside an "
                "identifier. Post-run analysis that needs true fitness belongs in "
                "veritas.diagnostics.retrospective."
            )

    def _depth(self, variant_id: str) -> int:
        cached = self._depth_cache.get(variant_id)
        if cached is None:
            cached = self._library.mutation_depth(variant_id)
            self._depth_cache[variant_id] = cached
        return cached

    @property
    def sealed(self) -> bool:
        return self._sealed

    @property
    def library(self) -> LibrarySpec:
        return self._library

    @property
    def library_size(self) -> int:
        return len(self._eligible)

    # -- recording -------------------------------------------------------------------
    def record_round(
        self,
        *,
        campaign: Hashable,
        round_index: int,
        proposed: Sequence[str],
        selected: Sequence[str],
    ) -> None:
        """
        Records one round of one campaign.

        `proposed` is the full candidate pool the generator put on the ballot;
        `selected` is the batch acquisition committed to. A selected candidate is
        necessarily a proposed one, so `selected` is unioned into the proposed tally -
        which is what makes this usable when all that survives of a run is its
        selection trace.
        """
        if self._sealed:
            raise RuntimeError(
                "this CoverageAccumulator is sealed; it can no longer record a round. "
                "Sealing is what makes retrospective fitness analysis safe: a sealed "
                "accumulator cannot influence a future selection."
            )
        proposed_set = self._round_proposed.setdefault(round_index, set())
        selected_set = self._round_selected.setdefault(round_index, set())
        proposed_count = 0
        selected_count = 0
        offered: Set[str] = set()

        for variant_id in proposed or ():
            self._require_identifier(variant_id, "proposed")
            proposed_count += 1
            offered.add(variant_id)
            if self._accept(variant_id):
                proposed_set.add(variant_id)

        for variant_id in selected or ():
            self._require_identifier(variant_id, "selected")
            selected_count += 1
            # A selected candidate was necessarily on the ballot. Counted as a proposal
            # slot only when the caller did not already list it, so a selection-trace-only
            # reconstruction still reports a truthful proposal total.
            if variant_id not in offered:
                proposed_count += 1
            if self._accept(variant_id):
                selected_set.add(variant_id)
                proposed_set.add(variant_id)

        self._round_proposed_total[round_index] = (
            self._round_proposed_total.get(round_index, 0) + proposed_count
        )
        self._round_selected_total[round_index] = (
            self._round_selected_total.get(round_index, 0) + selected_count
        )
        self._campaigns.add(campaign)

    def _accept(self, variant_id: str) -> bool:
        """Validates a proposal, counting (not raising on) the two ways it can be unusable."""
        try:
            self._depth(variant_id)
        except GenotypeParseError:
            self._unparseable += 1
            return False
        if self._eligible and variant_id not in self._eligible:
            self._off_universe += 1
            return False
        return True

    # -- sealing ---------------------------------------------------------------------
    @property
    def reachable_ids(self) -> frozenset:
        """
        Every distinct candidate generation ever proposed - the arm's reachable set.

        Identifiers only. This is what `veritas.diagnostics.retrospective` scores
        against the truth table, after sealing.
        """
        union: Set[str] = set()
        for ids in self._round_proposed.values():
            union |= ids
        return frozenset(union)

    @property
    def selected_ids(self) -> frozenset:
        union: Set[str] = set()
        for ids in self._round_selected.values():
            union |= ids
        return frozenset(union)

    def seal(self) -> CoverageReport:
        """Freezes the accumulator and returns its outcome-free report."""
        self._sealed = True
        size = max(len(self._eligible), 1)

        rows: List[RoundCoverage] = []
        cumulative_proposed: Set[str] = set()
        cumulative_selected: Set[str] = set()
        for round_index in sorted(self._round_proposed):
            proposed = self._round_proposed[round_index]
            selected = self._round_selected.get(round_index, set())
            cumulative_proposed |= proposed
            cumulative_selected |= selected
            rows.append(
                RoundCoverage(
                    round_index=round_index,
                    proposed_distinct=len(proposed),
                    proposed_total=self._round_proposed_total.get(round_index, 0),
                    cumulative_distinct=len(cumulative_proposed),
                    cumulative_share=len(cumulative_proposed) / size,
                    selected_distinct=len(selected),
                    selected_total=self._round_selected_total.get(round_index, 0),
                    cumulative_selected_distinct=len(cumulative_selected),
                    cumulative_selected_share=len(cumulative_selected) / size,
                    proposed_depth_histogram=self._histogram(proposed),
                    selected_depth_histogram=self._histogram(selected),
                )
            )

        union_proposed = cumulative_proposed
        union_selected = cumulative_selected
        return CoverageReport(
            library=self._library,
            library_size=len(self._eligible),
            max_library_depth=self._library.max_depth,
            universe_depth_histogram=dict(sorted(self._universe_depth_histogram.items())),
            campaigns=len(self._campaigns),
            rounds=tuple(rows),
            union_proposed_distinct=len(union_proposed),
            union_proposed_total=sum(self._round_proposed_total.values()),
            union_proposed_share=len(union_proposed) / size,
            union_proposed_depth_histogram=self._histogram(union_proposed),
            union_selected_distinct=len(union_selected),
            union_selected_total=sum(self._round_selected_total.values()),
            union_selected_share=len(union_selected) / size,
            union_selected_depth_histogram=self._histogram(union_selected),
            mean_proposed_depth=self._mean_depth(union_proposed),
            mean_selected_depth=self._mean_depth(union_selected),
            unparseable_proposals=self._unparseable,
            off_universe_proposals=self._off_universe,
        )

    def _histogram(self, ids: Iterable[str]) -> Dict[int, int]:
        histogram: Dict[int, int] = {}
        for variant_id in ids:
            depth = self._depth(variant_id)
            histogram[depth] = histogram.get(depth, 0) + 1
        return dict(sorted(histogram.items()))

    def _mean_depth(self, ids: Iterable[str]) -> float:
        depths = [self._depth(v) for v in ids]
        return (sum(depths) / len(depths)) if depths else 0.0


# ==============================================================================================
# Classification
# ==============================================================================================
@dataclass(frozen=True, slots=True)
class CollapseThresholds:
    """
    Calibrated against the measured TrpB figures, with the arithmetic written out so a
    future change to a number is an argument about evidence rather than about taste.

      min_union_share = 0.005
          The collapsed arm reached 294 / 159,999 = 0.18%, 2.7x below this line.
          Ridge reached 4,028 / 159,999 = 2.52%, 5x above it. Anything that covers
          less than half a percent of an eligible universe over a whole 96-campaign
          run was not searching it.

      min_depth_reach = 0.75   (mean proposed depth / deepest depth the library allows)
          The collapsed arm: 2.56 / 4 = 0.64. Ridge: 3.51 / 4 = 0.88. The optimum of
          this landscape is a depth-4 genotype, and 81% of the landscape is depth 4,
          so a generator averaging well below full depth is not reaching the part of
          the space where the answer lives.

      min_deep_enrichment = 0.5
          Share of proposals at maximum depth, divided by that depth's share of the
          eligible universe. 1.0 means the generator samples full-depth genotypes in
          proportion to how many there are. The collapsed arm scored 0.19 (45 of 294
          = 15% against 81% of the landscape); ridge scored 0.66. This is the metric
          that names the MECHANISM - depth truncation - rather than just the symptom.

      min_proposals_to_judge = 50
          Below this there is not enough evidence to call anything; the diagnosis says
          so instead of guessing.
    """

    min_union_share: float = 0.005
    min_depth_reach: float = 0.75
    min_deep_enrichment: float = 0.5
    min_proposals_to_judge: int = 50

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


HEALTHY = "healthy"
NARROW_GENERATION = "narrow-generation"
DEPTH_TRUNCATED = "depth-truncated"
COLLAPSED = "collapsed"
INSUFFICIENT_DATA = "insufficient-data"


@dataclass(frozen=True, slots=True)
class CoverageDiagnosis:
    """A classification plus every number that justifies it."""

    classification: str
    narrow_generation: bool
    depth_truncated: bool
    union_distinct: int
    library_size: int
    union_share: float
    mean_proposed_depth: float
    max_library_depth: int
    depth_reach: float
    deep_proposed: int
    deep_universe_share: float
    deep_enrichment: float
    selected_distinct: int
    campaigns: int
    thresholds: CollapseThresholds
    justification: str
    label: str = ""

    @property
    def is_collapse(self) -> bool:
        return self.classification in (NARROW_GENERATION, DEPTH_TRUNCATED, COLLAPSED)

    def to_dict(self) -> Dict[str, Any]:
        payload = asdict(self)
        payload["thresholds"] = self.thresholds.to_dict()
        payload["is_collapse"] = self.is_collapse
        return payload


def classify_coverage(
    report: CoverageReport,
    thresholds: Optional[CollapseThresholds] = None,
    *,
    label: str = "",
) -> CoverageDiagnosis:
    """
    Classifies one arm's candidate generation. Pure: it warns about nothing and raises
    only when handed the wrong type of report.

    `TypeError` on anything that is not a `CoverageReport` is deliberate: a
    `RetrospectiveCoverageReport` carries measured outcomes, and the classifier is a
    live-run component. The two must never be interchangeable.
    """
    if not isinstance(report, CoverageReport):
        raise TypeError(
            f"classify_coverage takes an outcome-free CoverageReport, got "
            f"{type(report).__name__}. Retrospective reports reference true fitness "
            "and must never be fed to a live-run classifier."
        )
    limits = thresholds or CollapseThresholds()

    size = max(report.library_size, 1)
    union = report.union_proposed_distinct
    share = report.union_proposed_share
    max_depth = report.max_library_depth or 1
    mean_depth = report.mean_proposed_depth
    depth_reach = mean_depth / max_depth if max_depth else 0.0

    deep_proposed = report.union_proposed_depth_histogram.get(max_depth, 0)
    deep_universe = report.universe_depth_histogram.get(max_depth, 0)
    deep_universe_share = deep_universe / size
    proposed_deep_share = deep_proposed / union if union else 0.0
    deep_enrichment = (
        proposed_deep_share / deep_universe_share if deep_universe_share > 0 else 0.0
    )

    if union < limits.min_proposals_to_judge:
        justification = (
            f"only {union} distinct candidate(s) recorded over {report.campaigns} "
            f"campaign(s); fewer than {limits.min_proposals_to_judge} is not enough "
            "evidence to classify generation coverage."
        )
        return CoverageDiagnosis(
            classification=INSUFFICIENT_DATA,
            narrow_generation=False,
            depth_truncated=False,
            union_distinct=union,
            library_size=report.library_size,
            union_share=share,
            mean_proposed_depth=mean_depth,
            max_library_depth=max_depth,
            depth_reach=depth_reach,
            deep_proposed=deep_proposed,
            deep_universe_share=deep_universe_share,
            deep_enrichment=deep_enrichment,
            selected_distinct=report.union_selected_distinct,
            campaigns=report.campaigns,
            thresholds=limits,
            justification=justification,
            label=label,
        )

    narrow = share < limits.min_union_share
    truncated = (
        depth_reach < limits.min_depth_reach or deep_enrichment < limits.min_deep_enrichment
    )
    if narrow and truncated:
        classification = COLLAPSED
    elif narrow:
        classification = NARROW_GENERATION
    elif truncated:
        classification = DEPTH_TRUNCATED
    else:
        classification = HEALTHY

    reasons: List[str] = [
        f"reached {union} distinct candidates ({share * 100:.2f}% of "
        f"{report.library_size} eligible) over {report.campaigns} campaign(s)",
        f"mean proposed depth {mean_depth:.2f} of {max_depth} "
        f"(depth reach {depth_reach:.2f})",
        f"{deep_proposed} proposal(s) at depth {max_depth} = "
        f"{proposed_deep_share * 100:.1f}% of what was proposed, against "
        f"{deep_universe_share * 100:.1f}% of the eligible set "
        f"(enrichment {deep_enrichment:.2f}x)",
    ]
    if narrow:
        reasons.append(
            f"NARROW GENERATION: coverage {share * 100:.2f}% is below the "
            f"{limits.min_union_share * 100:.2f}% floor, so acquisition was ranking an "
            "impoverished ballot"
        )
    if truncated:
        reasons.append(
            f"DEPTH TRUNCATED: depth reach {depth_reach:.2f} / enrichment "
            f"{deep_enrichment:.2f}x against floors {limits.min_depth_reach} / "
            f"{limits.min_deep_enrichment}, so the deepest region of the library was "
            "effectively unreachable"
        )
    if classification == HEALTHY:
        reasons.append("coverage and depth profile are both within calibrated bounds")

    return CoverageDiagnosis(
        classification=classification,
        narrow_generation=narrow,
        depth_truncated=truncated,
        union_distinct=union,
        library_size=report.library_size,
        union_share=share,
        mean_proposed_depth=mean_depth,
        max_library_depth=max_depth,
        depth_reach=depth_reach,
        deep_proposed=deep_proposed,
        deep_universe_share=deep_universe_share,
        deep_enrichment=deep_enrichment,
        selected_distinct=report.union_selected_distinct,
        campaigns=report.campaigns,
        thresholds=limits,
        justification="; ".join(reasons),
        label=label,
    )


def assert_generation_coverage(
    report: CoverageReport,
    *,
    thresholds: Optional[CollapseThresholds] = None,
    label: str = "",
    log: Optional[logging.Logger] = None,
) -> CoverageDiagnosis:
    """
    The run assertion. It WARNS AND CLASSIFIES; it never raises on a collapse.

    Deliberately not a hard gate. An arm that converges to the same answer from a
    broad, deep ballot every time is a legitimate deterministic result, and a guard
    that failed it would be wrong. What goes wrong silently is the opposite case -
    a varying, data-responsive policy searching 0.18% of its own library - so this
    names that case, attaches the arithmetic, and leaves the verdict to a human.
    """
    diagnosis = classify_coverage(report, thresholds, label=label)
    emitter = log or logger
    prefix = f"{label}: " if label else ""

    if diagnosis.classification == INSUFFICIENT_DATA:
        emitter.info("%sgeneration coverage not classified - %s", prefix, diagnosis.justification)
        return diagnosis
    if diagnosis.classification == HEALTHY:
        emitter.info(
            "%sgeneration coverage healthy - %s", prefix, diagnosis.justification
        )
        return diagnosis

    message = (
        f"{prefix}CANDIDATE GENERATION {diagnosis.classification.upper()}. "
        f"{diagnosis.justification}. This is a diagnosis, not a rejection: "
        "deterministic convergence from a broad ballot is legitimate, but these "
        "numbers say the ballot itself was restricted, so attribute the result to "
        "candidate generation before attributing it to acquisition."
    )
    emitter.warning(message)
    warnings.warn(message, GenerationCollapseWarning, stacklevel=2)
    return diagnosis
