"""
tests/test_generation_coverage_diagnostics.py - the diagnostic that would have caught
the TrpB `frontier` collapse, and the bounded candidate-selection trace.

The failure this locks in (measured, not hypothetical):

    On the published four-site TrpB landscape, the `frontier` configuration with the
    `dkl` and `gflownet` components ablated terminated ALL 96 campaigns at exactly
    fitness 0.7036 against a global optimum of 1.0.

    A terminal-variance check PASSES that configuration, and so does a data-invariance
    check: the round-1 batches genuinely vary with the starting sample (Jaccard 0.354,
    not 1.0). The policy is not data-blind, so neither guard fires. What was actually
    broken is candidate GENERATION: across 96 campaigns the arm only ever reached 294
    DISTINCT variants, 0.18% of the 159,999-variant eligible set, with a mean mutation
    depth of 2.56 and just 45 distinct depth-4 variants - while depth 4 holds 130,321
    variants (81% of the landscape) and contains the global optimum. Plain ridge
    regression reached 4,028 distinct variants (2.52%) at mean depth 3.51.

    The best fitness anywhere in that 294-variant reachable set is exactly 0.7036, so
    acquisition was choosing correctly from an impoverished ballot.

Acceptance: the classifier must label the collapse arm `collapsed` and the ridge
control `healthy`, from those figures, and must NEVER reject a run - legitimate
deterministic convergence has to stay passable.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import tempfile
import unittest
import warnings
from pathlib import Path

from veritas.diagnostics.coverage import (
    CollapseThresholds,
    CoverageAccumulator,
    CoverageReport,
    GenerationCollapseWarning,
    assert_generation_coverage,
    classify_coverage,
)
from veritas.diagnostics.depth_returns import (
    DepthReturnLedger,
    DepthReturnRow,
    MeasuredSelection,
    OutcomeDisclosureError,
    write_depth_returns_csv,
)
from veritas.diagnostics.genotype import (
    GenotypeParseError,
    LibrarySpec,
    infer_library_from_variant_ids,
)
from veritas.diagnostics.retrospective import (
    RetrospectiveAccessError,
    RetrospectiveCoverageReport,
    RetrospectiveFitnessAccess,
    analyse_reachable_set,
    depth_shell_profile,
)
from veritas.diagnostics.trace import (
    LEGACY_SELECTION_TRACE_COLUMNS,
    BoundedCandidateTrace,
    CandidateTraceRecord,
    TraceBudget,
    write_candidate_trace_csv,
)

# --- The real TrpB four-site library -----------------------------------------------------------
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
TRPB = LibrarySpec(positions=(183, 184, 227, 228), parent_residues="VFVS")

# The saturated 20^4 library: depth d holds C(4, d) * 19^d variants.
#   d=0: 1   d=1: 76   d=2: 2,166   d=3: 27,436   d=4: 130,321   -> 160,000 total
# Dropping the parent itself leaves the 159,999 eligible variants of the analysis.
UNIVERSE_SIZE = 159_999
UNIVERSE_DEPTH_4 = 130_321


def _notation(genotype: str) -> str:
    """Position-ordered notation for a genotype, omitting wild-type sites (as the data does)."""
    return ":".join(
        f"{TRPB.parent_residues[i]}{pos}{genotype[i]}"
        for i, pos in enumerate(TRPB.positions)
        if genotype[i] != TRPB.parent_residues[i]
    )


def _eligible_ids() -> list[str]:
    ids = []
    for combo in itertools.product(AA_ALPHABET, repeat=4):
        genotype = "".join(combo)
        if genotype != TRPB.parent_residues:
            ids.append(_notation(genotype))
    return ids


def _ids_at_depth(depth: int, count: int, offset: int = 0) -> list[str]:
    """`count` distinct eligible variant ids whose mutation depth is exactly `depth`."""
    out: list[str] = []
    for slots in itertools.combinations(range(4), depth):
        alternatives = [
            [aa for aa in AA_ALPHABET if aa != TRPB.parent_residues[s]] for s in slots
        ]
        for picks in itertools.product(*alternatives):
            genotype = list(TRPB.parent_residues)
            for slot, aa in zip(slots, picks):
                genotype[slot] = aa
            out.append(_notation("".join(genotype)))
            if len(out) >= count + offset:
                return out[offset:]
    raise AssertionError(f"library has fewer than {count + offset} variants at depth {depth}")


def _reachable_set(depth_counts: dict[int, int]) -> list[str]:
    ids: list[str] = []
    for depth, count in sorted(depth_counts.items()):
        ids.extend(_ids_at_depth(depth, count))
    assert len(ids) == len(set(ids)) == sum(depth_counts.values())
    return ids


# The measured collapse: 294 distinct variants, mean depth 2.5578 (reported 2.56),
# 45 of them at depth 4.
COLLAPSE_DEPTHS = {1: 26, 2: 123, 3: 100, 4: 45}
# The ridge control: 4,028 distinct variants at mean depth 3.51, 2,154 at depth 4.
RIDGE_DEPTHS = {2: 100, 3: 1_774, 4: 2_154}


class TestGenotypeCanonicalisation(unittest.TestCase):
    """
    Mutation strings OMIT wild-type positions. Canonicalising to full library length,
    filling the parent residue where a position is absent, is the step we got wrong
    once: `V183A:F184A:V227A` is a FOUR-site genotype whose fourth site is wild type.
    """

    def test_omitted_wild_type_position_is_filled_not_truncated(self) -> None:
        genotype = TRPB.canonical_genotype("V183A:F184A:V227A")
        self.assertEqual(genotype, "AAAS")
        self.assertEqual(len(genotype), len(TRPB.positions))
        self.assertEqual(TRPB.mutation_depth("V183A:F184A:V227A"), 3)

    def test_wild_type_notation_has_depth_zero(self) -> None:
        self.assertEqual(TRPB.canonical_genotype("WT"), "VFVS")
        self.assertEqual(TRPB.mutation_depth("WT"), 0)
        self.assertEqual(TRPB.mutation_depth(""), 0)

    def test_full_depth_genotype(self) -> None:
        self.assertEqual(TRPB.canonical_genotype("V183A:F184I:V227K:S228G"), "AIKG")
        self.assertEqual(TRPB.mutation_depth("V183A:F184I:V227K:S228G"), 4)

    def test_token_order_does_not_matter(self) -> None:
        self.assertEqual(
            TRPB.canonical_genotype("V227K:V183A"), TRPB.canonical_genotype("V183A:V227K")
        )

    def test_position_outside_the_library_is_rejected(self) -> None:
        with self.assertRaises(GenotypeParseError):
            TRPB.canonical_genotype("S65T")

    def test_parent_residue_mismatch_is_rejected(self) -> None:
        with self.assertRaises(GenotypeParseError):
            TRPB.canonical_genotype("A183V")

    def test_duplicate_position_is_rejected(self) -> None:
        with self.assertRaises(GenotypeParseError):
            TRPB.canonical_genotype("V183A:V183K")

    def test_unequal_length_genotype_is_never_zipped_or_truncated(self) -> None:
        with self.assertRaises(GenotypeParseError):
            TRPB.depth_of_genotype("AAA")
        with self.assertRaises(GenotypeParseError):
            TRPB.depth_of_genotype("AAAAA")

    def test_library_spec_rejects_inconsistent_definition(self) -> None:
        with self.assertRaises(ValueError):
            LibrarySpec(positions=(183, 184), parent_residues="VFV")
        with self.assertRaises(ValueError):
            LibrarySpec(positions=(184, 183), parent_residues="FV")

    def test_library_is_inferable_from_variant_ids_alone(self) -> None:
        spec = infer_library_from_variant_ids(
            ["V183A:F184I", "V227K:S228G", "V183C:S228A"]
        )
        self.assertEqual(spec, TRPB)

    def test_inference_rejects_conflicting_parent_residues(self) -> None:
        with self.assertRaises(GenotypeParseError):
            infer_library_from_variant_ids(["V183A", "A183V"])


class TestCoverageAccounting(unittest.TestCase):
    """Task 1: distinct PROPOSED candidates per round and cumulatively, with depth histograms."""

    def setUp(self) -> None:
        self.accumulator = CoverageAccumulator(
            library=TRPB,
            eligible_ids=[
                _notation("AFVS"),   # depth 1
                _notation("AIVS"),   # depth 2
                _notation("AIKS"),   # depth 3
                _notation("AIKG"),   # depth 4
                _notation("CIKG"),   # depth 4
            ],
        )

    def test_counts_proposed_not_just_selected(self) -> None:
        self.accumulator.record_round(
            campaign="seed0",
            round_index=1,
            proposed=[_notation("AFVS"), _notation("AIVS"), _notation("AIKS")],
            selected=[_notation("AIKS")],
        )
        report = self.accumulator.seal()
        row = report.rounds[0]
        self.assertEqual(row.proposed_distinct, 3)
        self.assertEqual(row.selected_distinct, 1)
        self.assertEqual(report.union_proposed_distinct, 3)
        self.assertEqual(report.union_selected_distinct, 1)
        self.assertAlmostEqual(report.union_proposed_share, 3 / 5)
        self.assertAlmostEqual(report.union_selected_share, 1 / 5)

    def test_selected_candidates_count_as_proposed(self) -> None:
        self.accumulator.record_round(
            campaign="seed0", round_index=1, proposed=[], selected=[_notation("AIKG")]
        )
        report = self.accumulator.seal()
        self.assertEqual(report.union_proposed_distinct, 1)

    def test_cumulative_union_is_deduplicated_across_rounds_and_campaigns(self) -> None:
        a, b, c = _notation("AFVS"), _notation("AIVS"), _notation("AIKS")
        self.accumulator.record_round(campaign="s0", round_index=1, proposed=[a, b], selected=[a])
        self.accumulator.record_round(campaign="s0", round_index=2, proposed=[b, c], selected=[c])
        self.accumulator.record_round(campaign="s1", round_index=1, proposed=[a], selected=[a])
        report = self.accumulator.seal()
        self.assertEqual(report.campaigns, 2)
        by_round = {row.round_index: row for row in report.rounds}
        self.assertEqual(by_round[1].proposed_distinct, 2)      # a, b over both campaigns
        self.assertEqual(by_round[1].cumulative_distinct, 2)
        self.assertEqual(by_round[2].proposed_distinct, 2)      # b, c
        self.assertEqual(by_round[2].cumulative_distinct, 3)    # a, b, c
        self.assertEqual(report.union_proposed_distinct, 3)
        self.assertEqual(by_round[1].proposed_total, 3)         # 2 + 1 proposals, with repeats

    def test_depth_histograms_for_proposed_and_selected(self) -> None:
        self.accumulator.record_round(
            campaign="s0",
            round_index=1,
            proposed=[_notation("AFVS"), _notation("AIKG"), _notation("CIKG")],
            selected=[_notation("AIKG")],
        )
        report = self.accumulator.seal()
        self.assertEqual(dict(report.union_proposed_depth_histogram), {1: 1, 4: 2})
        self.assertEqual(dict(report.union_selected_depth_histogram), {4: 1})
        self.assertAlmostEqual(report.mean_proposed_depth, (1 + 4 + 4) / 3)
        self.assertAlmostEqual(report.mean_selected_depth, 4.0)

    def test_universe_depth_histogram_is_recorded(self) -> None:
        report = self.accumulator.seal()
        self.assertEqual(
            dict(report.universe_depth_histogram), {1: 1, 2: 1, 3: 1, 4: 2}
        )
        self.assertEqual(report.library_size, 5)
        self.assertEqual(report.max_library_depth, 4)

    def test_recording_rejects_anything_that_is_not_a_bare_identifier(self) -> None:
        with self.assertRaises(TypeError):
            self.accumulator.record_round(
                campaign="s0",
                round_index=1,
                proposed=[(_notation("AFVS"), 0.7036)],   # an (id, fitness) pair
                selected=[],
            )
        with self.assertRaises(TypeError):
            self.accumulator.record_round(
                campaign="s0", round_index=1, proposed=[], selected=[{"id": "x", "y": 1.0}]
            )

    def test_unknown_variant_is_counted_but_flagged(self) -> None:
        self.accumulator.record_round(
            campaign="s0", round_index=1, proposed=[_notation("AFVS"), "V183A:F184Z"], selected=[]
        )
        report = self.accumulator.seal()
        self.assertEqual(report.unparseable_proposals, 1)
        self.assertEqual(report.union_proposed_distinct, 1)

    def test_sealing_is_final(self) -> None:
        self.accumulator.record_round(
            campaign="s0", round_index=1, proposed=[_notation("AFVS")], selected=[]
        )
        self.assertFalse(self.accumulator.sealed)
        self.accumulator.seal()
        self.assertTrue(self.accumulator.sealed)
        with self.assertRaises(RuntimeError):
            self.accumulator.record_round(
                campaign="s0", round_index=2, proposed=[_notation("AIVS")], selected=[]
            )


class TestClassifierOnTheRealFailure(unittest.TestCase):
    """
    Task 2 / exit gate. Calibrated against the measured figures, and it WARNS - it
    never rejects, because legitimate deterministic convergence must stay passable.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.eligible = _eligible_ids()
        assert len(cls.eligible) == UNIVERSE_SIZE, len(cls.eligible)

    def _report(self, depth_counts: dict[int, int]) -> CoverageReport:
        accumulator = CoverageAccumulator(library=TRPB, eligible_ids=self.eligible)
        ids = _reachable_set(depth_counts)
        # 96 campaigns x 4 rounds, all drawing from the same impoverished reachable set.
        for campaign in range(96):
            for round_index in range(1, 5):
                start = ((campaign * 4) + round_index) * 25 % len(ids)
                window = [ids[(start + i) % len(ids)] for i in range(25)]
                accumulator.record_round(
                    campaign=f"seed{campaign}",
                    round_index=round_index,
                    proposed=window,
                    selected=window[:25],
                )
        # Guarantee every member of the reachable set is seen, as in the real run.
        accumulator.record_round(campaign="seed0", round_index=4, proposed=ids, selected=[])
        return accumulator.seal()

    def test_universe_depth_structure_matches_the_analysis(self) -> None:
        report = CoverageAccumulator(library=TRPB, eligible_ids=self.eligible).seal()
        self.assertEqual(report.library_size, UNIVERSE_SIZE)
        self.assertEqual(report.universe_depth_histogram[4], UNIVERSE_DEPTH_4)
        self.assertEqual(report.max_library_depth, 4)
        self.assertAlmostEqual(UNIVERSE_DEPTH_4 / UNIVERSE_SIZE, 0.8145, places=4)

    def test_the_real_trpb_collapse_is_classified_collapsed(self) -> None:
        report = self._report(COLLAPSE_DEPTHS)
        # The figures the diagnosis has to reproduce.
        self.assertEqual(report.union_proposed_distinct, 294)
        self.assertAlmostEqual(report.union_proposed_share, 0.0018, places=4)
        self.assertAlmostEqual(report.mean_proposed_depth, 2.56, places=2)
        self.assertEqual(report.union_proposed_depth_histogram[4], 45)

        diagnosis = classify_coverage(report)
        self.assertEqual(diagnosis.classification, "collapsed")
        self.assertTrue(diagnosis.narrow_generation)
        self.assertTrue(diagnosis.depth_truncated)
        # The numbers that justify the label travel with it.
        self.assertAlmostEqual(diagnosis.union_share, 0.0018, places=4)
        self.assertAlmostEqual(diagnosis.mean_proposed_depth, 2.56, places=2)
        self.assertLess(diagnosis.depth_reach, 0.75)
        self.assertLess(diagnosis.deep_enrichment, 0.5)
        self.assertIn("294", diagnosis.justification)
        self.assertIn("0.18", diagnosis.justification)

    def test_the_ridge_control_is_classified_healthy(self) -> None:
        report = self._report(RIDGE_DEPTHS)
        self.assertEqual(report.union_proposed_distinct, 4_028)
        self.assertAlmostEqual(report.union_proposed_share, 0.0252, places=4)
        self.assertAlmostEqual(report.mean_proposed_depth, 3.51, places=2)

        diagnosis = classify_coverage(report)
        self.assertEqual(diagnosis.classification, "healthy")
        self.assertFalse(diagnosis.narrow_generation)
        self.assertFalse(diagnosis.depth_truncated)
        self.assertGreater(diagnosis.depth_reach, 0.75)
        self.assertGreater(diagnosis.deep_enrichment, 0.5)

    def test_narrow_but_deep_generation_is_distinguished_from_depth_truncation(self) -> None:
        # Broad coverage, but the generator never leaves depth <= 2: depth-truncated only.
        report = self._report({1: 76, 2: 2_166})
        diagnosis = classify_coverage(report)
        self.assertTrue(diagnosis.depth_truncated)
        self.assertEqual(diagnosis.classification, "depth-truncated")

        # Deep enough, but almost nothing is ever reached: narrow-generation only.
        report = self._report({4: 300})
        diagnosis = classify_coverage(report)
        self.assertTrue(diagnosis.narrow_generation)
        self.assertFalse(diagnosis.depth_truncated)
        self.assertEqual(diagnosis.classification, "narrow-generation")

    def test_the_run_assertion_warns_and_never_rejects(self) -> None:
        report = self._report(COLLAPSE_DEPTHS)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            diagnosis = assert_generation_coverage(report, label="frontier/ablate=dkl,gflownet")
        self.assertEqual(diagnosis.classification, "collapsed")
        self.assertTrue(
            any(issubclass(w.category, GenerationCollapseWarning) for w in caught),
            [str(w.message) for w in caught],
        )
        self.assertTrue(any("frontier" in str(w.message) for w in caught))

    def test_a_healthy_run_emits_no_warning(self) -> None:
        report = self._report(RIDGE_DEPTHS)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            diagnosis = assert_generation_coverage(report, label="ridge")
        self.assertEqual(diagnosis.classification, "healthy")
        self.assertEqual(
            [w for w in caught if issubclass(w.category, GenerationCollapseWarning)], []
        )

    def test_deterministic_convergence_is_not_failed(self) -> None:
        """
        An arm that converges to one answer every time, from a broad and deep ballot,
        is a legitimate result. The diagnostic classifies; it must not reject.
        """
        accumulator = CoverageAccumulator(library=TRPB, eligible_ids=self.eligible)
        pool = _reachable_set(RIDGE_DEPTHS)
        winner = _ids_at_depth(4, 1, offset=9_000)[0]
        for campaign in range(96):
            for round_index in range(1, 5):
                start = (campaign * 97 + round_index * 29) % len(pool)
                window = [pool[(start + i) % len(pool)] for i in range(25)]
                accumulator.record_round(
                    campaign=f"seed{campaign}",
                    round_index=round_index,
                    proposed=window + [winner],
                    selected=[winner],      # identical terminal pick in every campaign
                )
        accumulator.record_round(campaign="seed0", round_index=4, proposed=pool, selected=[])
        report = accumulator.seal()
        self.assertEqual(report.union_selected_distinct, 1)     # zero terminal variance
        diagnosis = assert_generation_coverage(report, label="deterministic")
        self.assertEqual(diagnosis.classification, "healthy")

    def test_too_little_evidence_is_reported_as_such_not_as_collapse(self) -> None:
        accumulator = CoverageAccumulator(library=TRPB, eligible_ids=self.eligible)
        accumulator.record_round(
            campaign="s0", round_index=1, proposed=_ids_at_depth(4, 3), selected=[]
        )
        diagnosis = classify_coverage(accumulator.seal())
        self.assertEqual(diagnosis.classification, "insufficient-data")

    def test_thresholds_are_explicit_and_reported(self) -> None:
        report = self._report(COLLAPSE_DEPTHS)
        thresholds = CollapseThresholds(min_union_share=1e-9, min_depth_reach=0.0,
                                        min_deep_enrichment=0.0)
        diagnosis = classify_coverage(report, thresholds=thresholds)
        self.assertEqual(diagnosis.classification, "healthy")
        self.assertEqual(diagnosis.thresholds, thresholds)

    def test_diagnosis_serialises_for_the_run_folder(self) -> None:
        diagnosis = classify_coverage(self._report(COLLAPSE_DEPTHS))
        payload = json.loads(json.dumps(diagnosis.to_dict()))
        self.assertEqual(payload["classification"], "collapsed")
        self.assertEqual(payload["union_distinct"], 294)
        self.assertNotIn("fitness", json.dumps(payload).lower())


class TestHiddenOutcomesCannotReachSelection(unittest.TestCase):
    """
    Task 4. Coverage metrics that reference true fitness - "the best fitness available
    in the reachable set is 0.7036" - are RETROSPECTIVE ONLY. This is enforced by the
    types, not by convention.
    """

    def setUp(self) -> None:
        self.accumulator = CoverageAccumulator(
            library=TRPB, eligible_ids=[_notation("AFVS"), _notation("AIKG"), _notation("CIKG")]
        )
        self.accumulator.record_round(
            campaign="s0", round_index=1, proposed=[_notation("AFVS")], selected=[_notation("AFVS")]
        )
        self.truth = {
            _notation("AFVS"): 0.7036,
            _notation("AIKG"): 1.0,
            _notation("CIKG"): 0.42,
        }

    def test_live_accumulator_cannot_be_given_a_fitness_attribute(self) -> None:
        with self.assertRaises(AttributeError):
            self.accumulator.fitness = self.truth          # type: ignore[attr-defined]
        with self.assertRaises(AttributeError):
            self.accumulator.revealed_fitness = 0.7036     # type: ignore[attr-defined]

    def test_live_report_holds_no_outcome_field_and_is_immutable(self) -> None:
        report = self.accumulator.seal()
        self.assertEqual(
            [f for f in report.field_names() if "fit" in f or "score" in f or "truth" in f], []
        )
        with self.assertRaises(Exception):
            report.best_available_fitness = 0.7036         # type: ignore[attr-defined]

    def test_retrospective_analysis_refuses_an_unsealed_accumulator(self) -> None:
        with self.assertRaises(RetrospectiveAccessError):
            analyse_reachable_set(
                self.accumulator,
                fitness_lookup=self.truth,
                access=RetrospectiveFitnessAccess(reason="post-run analysis"),
            )

    def test_retrospective_analysis_requires_an_explicit_access_token(self) -> None:
        self.accumulator.seal()
        with self.assertRaises(RetrospectiveAccessError):
            analyse_reachable_set(self.accumulator, fitness_lookup=self.truth, access=None)
        with self.assertRaises(RetrospectiveAccessError):
            analyse_reachable_set(
                self.accumulator, fitness_lookup=self.truth, access="yes please"
            )
        with self.assertRaises(ValueError):
            RetrospectiveFitnessAccess(reason="")

    def test_retrospective_analysis_reports_the_best_reachable_fitness(self) -> None:
        self.accumulator.seal()
        result = analyse_reachable_set(
            self.accumulator,
            fitness_lookup=self.truth,
            access=RetrospectiveFitnessAccess(reason="attributing the 0.7036 ceiling"),
        )
        self.assertIsInstance(result, RetrospectiveCoverageReport)
        self.assertAlmostEqual(result.best_reachable_fitness, 0.7036)
        self.assertAlmostEqual(result.best_eligible_fitness, 1.0)
        self.assertAlmostEqual(result.reachable_fitness_ceiling_gap, 1.0 - 0.7036)
        self.assertEqual(result.reason, "attributing the 0.7036 ceiling")

    def test_a_sealed_accumulator_can_never_influence_another_selection(self) -> None:
        self.accumulator.seal()
        analyse_reachable_set(
            self.accumulator,
            fitness_lookup=self.truth,
            access=RetrospectiveFitnessAccess(reason="post-run"),
        )
        with self.assertRaises(RuntimeError):
            self.accumulator.record_round(
                campaign="s0", round_index=2, proposed=[_notation("AIKG")], selected=[]
            )

    def test_the_live_classifier_rejects_a_retrospective_report(self) -> None:
        self.accumulator.seal()
        result = analyse_reachable_set(
            self.accumulator,
            fitness_lookup=self.truth,
            access=RetrospectiveFitnessAccess(reason="post-run"),
        )
        with self.assertRaises(TypeError):
            classify_coverage(result)                      # type: ignore[arg-type]

    def test_retrospective_report_is_a_disjoint_type(self) -> None:
        self.assertFalse(issubclass(RetrospectiveCoverageReport, CoverageReport))
        self.assertFalse(issubclass(CoverageReport, RetrospectiveCoverageReport))


class TestBoundedCandidateTrace(unittest.TestCase):
    """
    Task 3. A full per-candidate trace is not viable - GB1 is ~149,000 candidates x 4
    rounds x 48 seeds x many arms - so the trace is bounded FROM THE START: every
    selected candidate, the top-N ranked, and a seeded random sample of the rejected.
    """

    def _round(self, trace: BoundedCandidateTrace, n: int = 400) -> None:
        ids = _ids_at_depth(4, n)
        trace.record_round(
            round_index=1,
            candidates=ids,
            sources=["gflownet" if i % 3 else "combinatorial" for i in range(n)],
            predicted_fitness=[0.001 * i for i in range(n)],
            posterior_uncertainty=[0.05] * n,
            prior_contribution=[0.01 * (i % 7) for i in range(n)],
            acquisition_score=[0.002 * i for i in range(n)],
            pareto_rank=[i % 4 for i in range(n)],
            pareto_contribution=[1.0 / (1.0 + 0.3 * (i % 4)) for i in range(n)],
            quality_score=[0.002 * i / (1.0 + 0.3 * (i % 4)) for i in range(n)],
            selected_indices=[n - 1, n - 2, 5],
        )

    def test_trace_is_bounded_by_the_budget(self) -> None:
        trace = BoundedCandidateTrace(
            campaign="seed0", budget=TraceBudget(top_n=10, rejected_sample=5, seed=1)
        )
        self._round(trace, n=400)
        self.assertLessEqual(len(trace.records), 3 + 10 + 5)
        self.assertGreater(len(trace.records), 10)

    def test_every_selected_candidate_is_recorded(self) -> None:
        trace = BoundedCandidateTrace(
            campaign="seed0", budget=TraceBudget(top_n=2, rejected_sample=1, seed=1)
        )
        self._round(trace, n=400)
        selected = [r for r in trace.records if r.status == "selected"]
        self.assertEqual(len(selected), 3)

    def test_records_carry_every_required_field(self) -> None:
        trace = BoundedCandidateTrace(campaign="seed0", budget=TraceBudget(top_n=5, rejected_sample=3))
        self._round(trace)
        required = {
            "variant_id", "generation_source", "predicted_fitness", "posterior_uncertainty",
            "prior_contribution", "acquisition_score", "pareto_contribution",
            "diversity_adjustment", "final_rank", "status", "rejection_reason",
        }
        self.assertTrue(required.issubset(set(CandidateTraceRecord.field_names())))
        for record in trace.records:
            self.assertTrue(record.variant_id)
            self.assertIn(record.generation_source, {"gflownet", "combinatorial"})
            self.assertIn(record.status, {"selected", "ranked", "rejected"})
            self.assertGreaterEqual(record.final_rank, 1)
            if record.status == "selected":
                self.assertEqual(record.rejection_reason, "")
            else:
                self.assertTrue(record.rejection_reason)

    def test_final_rank_follows_the_quality_score(self) -> None:
        trace = BoundedCandidateTrace(campaign="s", budget=TraceBudget(top_n=20, rejected_sample=0))
        self._round(trace, n=100)
        ranked = sorted(trace.records, key=lambda r: r.final_rank)
        self.assertEqual(ranked[0].final_rank, 1)
        for earlier, later in zip(ranked, ranked[1:]):
            self.assertGreaterEqual(earlier.quality_score, later.quality_score - 1e-12)

    def test_diversity_adjustment_records_the_batch_selection_shift(self) -> None:
        trace = BoundedCandidateTrace(campaign="s", budget=TraceBudget(top_n=5, rejected_sample=0))
        self._round(trace, n=100)
        promoted = [r for r in trace.records if r.status == "selected" and r.variant_id]
        # Candidate index 5 has a poor quality score but was selected by the diversity
        # stage, so its adjustment must be positive (promoted past better-ranked peers).
        self.assertTrue(any(r.diversity_adjustment > 0 for r in promoted))

    def test_sampling_is_reproducible_for_a_fixed_seed(self) -> None:
        def sample() -> list[str]:
            trace = BoundedCandidateTrace(
                campaign="s", budget=TraceBudget(top_n=5, rejected_sample=7, seed=99)
            )
            self._round(trace, n=300)
            return [r.variant_id for r in trace.records if r.status == "rejected"]

        self.assertEqual(sample(), sample())
        self.assertEqual(len(sample()), 7)

    def test_trace_feeds_the_coverage_accumulator_with_the_unbounded_pool(self) -> None:
        accumulator = CoverageAccumulator(library=TRPB, eligible_ids=_ids_at_depth(4, 500))
        trace = BoundedCandidateTrace(
            campaign="seed0",
            budget=TraceBudget(top_n=5, rejected_sample=5),
            coverage=accumulator,
        )
        self._round(trace, n=400)
        report = accumulator.seal()
        # Coverage sees all 400 proposals even though only ~13 rows were traced.
        self.assertEqual(report.union_proposed_distinct, 400)
        self.assertEqual(report.union_selected_distinct, 3)
        self.assertLess(len(trace.records), 30)

    def test_record_round_rejects_ragged_inputs(self) -> None:
        trace = BoundedCandidateTrace(campaign="s")
        with self.assertRaises(ValueError):
            trace.record_round(
                round_index=1,
                candidates=["V183A", "F184I"],
                sources=["combinatorial"],
                predicted_fitness=[0.1, 0.2],
                posterior_uncertainty=[0.1, 0.2],
                prior_contribution=[0.0, 0.0],
                acquisition_score=[0.1, 0.2],
                pareto_rank=[0, 1],
                pareto_contribution=[1.0, 0.8],
                quality_score=[0.1, 0.2],
                selected_indices=[0],
            )

    def test_csv_round_trip(self) -> None:
        trace = BoundedCandidateTrace(campaign="seed0", budget=TraceBudget(top_n=5, rejected_sample=3))
        self._round(trace, n=200)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate_trace.csv"
            write_candidate_trace_csv(path, trace.records, stamp="provenance=TEST")
            with open(path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), len(trace.records))
        self.assertEqual(rows[0]["provenance"], "provenance=TEST")
        self.assertIn("generation_source", rows[0])
        self.assertIn("prior_contribution", rows[0])

    def test_trace_never_carries_a_revealed_outcome(self) -> None:
        fields = set(CandidateTraceRecord.field_names())
        self.assertNotIn("revealed_fitness", fields)
        self.assertNotIn("true_fitness", fields)
        # `predicted_fitness` is a model output, not an outcome; it must be the only
        # fitness-shaped column.
        self.assertEqual({f for f in fields if "fitness" in f}, {"predicted_fitness"})


class TestLegacySelectionTraceCompatibility(unittest.TestCase):
    """Task 5: `selection_trace.csv` keeps its exact columns - analysis scripts read it."""

    def test_legacy_columns_are_pinned(self) -> None:
        self.assertEqual(
            LEGACY_SELECTION_TRACE_COLUMNS,
            ("provenance", "method", "seed", "round", "rank", "variant_id", "status",
             "reason", "revealed_fitness"),
        )

    def test_benchmark_writer_still_emits_exactly_those_columns(self) -> None:
        import benchmark_panoptes_isolated as iso

        outcome = iso.RunOutcome(method="ridge", seed=0, status="completed")
        outcome.trace.append(
            iso.TraceRow(method="ridge", seed=0, round_index=1, rank=0,
                         variant_id="V183A:F184I", status="accepted", reason="",
                         revealed_fitness=0.7036)
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "selection_trace.csv"
            iso.write_trace_csv(path, [outcome], stamp="provenance=TEST")
            with open(path, newline="", encoding="utf-8") as handle:
                reader = csv.reader(handle)
                header = tuple(next(reader))
                first = next(reader)
        self.assertEqual(header, LEGACY_SELECTION_TRACE_COLUMNS)
        self.assertEqual(first[5], "V183A:F184I")
        self.assertEqual(first[8], "0.7036")


class TestBenchmarkHook(unittest.TestCase):
    """The single hook the benchmark calls, derived from the legacy trace rows."""

    def test_hook_builds_coverage_from_legacy_trace_rows(self) -> None:
        from veritas.diagnostics.benchmark_hook import emit_coverage_artifacts

        import benchmark_panoptes_isolated as iso

        eligible = _ids_at_depth(1, 76) + _ids_at_depth(2, 2_166) + _ids_at_depth(4, 40_000)
        reachable = _ids_at_depth(1, 60)
        outcome = iso.RunOutcome(method="panoptes_default", seed=0, status="completed")
        outcome.trace.append(
            iso.TraceRow(method="panoptes_default", seed=0, round_index=0, rank=0,
                         variant_id=_ids_at_depth(4, 1)[0], status="initial", reason="",
                         revealed_fitness=0.1)
        )
        for rank, vid in enumerate(reachable):
            outcome.trace.append(
                iso.TraceRow(method="panoptes_default", seed=0, round_index=1, rank=rank,
                             variant_id=vid, status="accepted", reason="", revealed_fitness=0.2)
            )
        outcome.trace.append(
            iso.TraceRow(method="panoptes_default", seed=0, round_index=1, rank=99,
                         variant_id=reachable[0], status="rejected",
                         reason="duplicate_within_batch", revealed_fitness=None)
        )

        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            diagnoses = emit_coverage_artifacts(
                run_dir,
                [outcome],
                eligible_ids=eligible,
                stamp="provenance=TEST",
                top1_predicate=lambda vid: vid == reachable[0],
                run_id="run_test",
            )
            payload = json.loads((run_dir / "generation_coverage.json").read_text(encoding="utf-8"))
            with open(run_dir / "generation_coverage.csv", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            with open(run_dir / "depth_returns.csv", newline="", encoding="utf-8") as handle:
                depth_rows = list(csv.DictReader(handle))

        self.assertEqual(len(diagnoses), 1)
        self.assertEqual(diagnoses[0].classification, "collapsed")
        # The shared starting sample is not selector-generated, so it is excluded.
        self.assertEqual(diagnoses[0].union_distinct, 60)
        self.assertEqual(payload["arms"][0]["method"], "panoptes_default")
        self.assertTrue(rows)
        self.assertEqual(rows[0]["provenance"], "provenance=TEST")

        # Per-depth returns travel with the coverage artifacts, tidy long.
        self.assertEqual(len(depth_rows), 1)              # one depth present
        self.assertEqual(depth_rows[0]["method"], "panoptes_default")
        self.assertEqual(depth_rows[0]["depth"], "1")
        self.assertEqual(depth_rows[0]["selected"], "60")
        self.assertEqual(depth_rows[0]["proposed"], "61")  # 60 accepted + 1 rejected
        self.assertEqual(depth_rows[0]["top1_hits"], "1")
        self.assertEqual(payload["depth_returns"]["mean_selected_depth"]["panoptes_default"], 1.0)
        self.assertEqual(
            payload["depth_returns"]["new_hits_per_round"]["panoptes_default"], {"1": 1.0}
        )

    def test_hook_never_raises_on_a_malformed_run(self) -> None:
        from veritas.diagnostics.benchmark_hook import emit_coverage_artifacts

        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(
                emit_coverage_artifacts(Path(tmp), [], eligible_ids=[], stamp="x"), []
            )
            self.assertEqual(
                emit_coverage_artifacts(Path(tmp) / "nope", None, eligible_ids=None, stamp="x"),
                [],
            )


class TestDepthReturnLedger(unittest.TestCase):
    """
    Per-round, per-depth RETURNS - the measurement the depth-adaptive-allocation
    hypothesis needs. Recorded during a run, readable only once the run is over.
    """

    def _ledger(self) -> DepthReturnLedger:
        ledger = DepthReturnLedger(TRPB, run_id="run_test")
        d1 = _ids_at_depth(1, 2)
        d2 = _ids_at_depth(2, 1)
        d3 = _ids_at_depth(3, 1)
        d4 = _ids_at_depth(4, 3)
        ledger.record_round(
            method="reliability",
            seed=48,
            round_index=1,
            proposed=d1 + d2 + d3 + d4,
            measured=[
                MeasuredSelection(d1[0], 0.10, True),
                MeasuredSelection(d1[1], 0.20, False),
                MeasuredSelection(d3[0], 0.30, False),
            ],
        )
        ledger.record_round(
            method="reliability",
            seed=48,
            round_index=2,
            proposed=d4,
            measured=[
                MeasuredSelection(d4[0], 0.70, True),
                MeasuredSelection(d4[1], 0.90, True),
            ],
        )
        ledger.record_round(
            method="reliability",
            seed=49,
            round_index=1,
            proposed=d2,
            measured=[MeasuredSelection(d2[0], 0.40, False)],
        )
        ledger.record_round(
            method="reliability",
            seed=49,
            round_index=2,
            proposed=d4,
            measured=[MeasuredSelection(d4[2], 0.50, False)],
        )
        return ledger

    def test_ledger_is_write_only_until_sealed(self) -> None:
        ledger = self._ledger()
        with self.assertRaises(OutcomeDisclosureError):
            _ = ledger.rows
        self.assertFalse(ledger.sealed)

    def test_recording_discloses_nothing(self) -> None:
        ledger = DepthReturnLedger(TRPB)
        returned = ledger.record_round(
            method="m", seed=0, round_index=1,
            proposed=_ids_at_depth(4, 1),
            measured=[MeasuredSelection(_ids_at_depth(4, 1)[0], 0.7036, True)],
        )
        self.assertIsNone(returned)

    def test_ledger_cannot_be_given_another_attribute(self) -> None:
        ledger = DepthReturnLedger(TRPB)
        with self.assertRaises(AttributeError):
            ledger.landscape = {"x": 1.0}      # type: ignore[attr-defined]

    def test_outcomes_must_arrive_as_measured_selections(self) -> None:
        ledger = DepthReturnLedger(TRPB)
        variant = _ids_at_depth(4, 1)[0]
        for bad in ([(variant, 0.7036)], [{"variant_id": variant, "fitness": 0.7}], [variant]):
            with self.assertRaises(TypeError):
                ledger.record_round(method="m", seed=0, round_index=1, measured=bad)

    def test_sealing_is_final(self) -> None:
        ledger = self._ledger()
        ledger.seal()
        self.assertTrue(ledger.sealed)
        with self.assertRaises(RuntimeError):
            ledger.record_round(method="reliability", seed=48, round_index=3)

    def test_rows_are_tidy_long_keyed_by_run_method_seed_round_depth(self) -> None:
        table = self._ledger().seal()
        self.assertEqual(
            DepthReturnRow.field_names()[:5],
            ("run_id", "method", "seed", "round_index", "depth"),
        )
        keys = [(r.run_id, r.method, r.seed, r.round_index, r.depth) for r in table.rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertTrue(all(r.run_id == "run_test" for r in table.rows))

        cell = next(r for r in table.rows if r.seed == "48" and r.round_index == 1 and r.depth == 1)
        self.assertEqual((cell.proposed, cell.selected, cell.top1_hits), (2, 2, 1))
        self.assertAlmostEqual(cell.mean_selected_fitness, 0.15)
        self.assertAlmostEqual(cell.max_selected_fitness, 0.20)

        deep = next(r for r in table.rows if r.seed == "48" and r.round_index == 2 and r.depth == 4)
        self.assertEqual((deep.selected, deep.top1_hits), (2, 2))
        self.assertAlmostEqual(deep.mean_selected_fitness, 0.80)
        self.assertAlmostEqual(deep.max_selected_fitness, 0.90)

    def test_cumulative_columns_accumulate_within_a_campaign_and_depth(self) -> None:
        ledger = DepthReturnLedger(TRPB)
        d4 = _ids_at_depth(4, 2)
        ledger.record_round(method="m", seed=1, round_index=1,
                            measured=[MeasuredSelection(d4[0], 0.5, True)])
        ledger.record_round(method="m", seed=1, round_index=2,
                            measured=[MeasuredSelection(d4[1], 0.6, True)])
        rows = {r.round_index: r for r in ledger.seal().rows}
        self.assertEqual((rows[1].cumulative_selected, rows[1].cumulative_top1_hits), (1, 1))
        self.assertEqual((rows[2].cumulative_selected, rows[2].cumulative_top1_hits), (2, 2))

    def test_roll_ups_reproduce_the_depth_and_return_trajectory(self) -> None:
        table = self._ledger().seal()
        self.assertAlmostEqual(table.mean_selected_depth("reliability"), 19 / 7)
        self.assertEqual(table.new_hits_per_round("reliability"), {1: 0.5, 2: 1.0})
        shells = table.shell_returns("reliability", shallow_max_depth=2)
        self.assertEqual(shells["shallow"]["selected"], 3.0)
        self.assertAlmostEqual(shells["shallow"]["hit_rate"], 1 / 3)
        self.assertEqual(shells["deep"]["selected"], 4.0)
        self.assertAlmostEqual(shells["deep"]["hit_rate"], 0.5)

    def test_unparseable_selection_is_counted_not_raised(self) -> None:
        ledger = DepthReturnLedger(TRPB)
        ledger.record_round(
            method="m", seed=0, round_index=1,
            proposed=["S65T"],
            measured=[MeasuredSelection("S65T", 0.5, False)],
        )
        self.assertEqual(ledger.seal().unparseable_selections, 2)

    def test_csv_is_tidy_long(self) -> None:
        table = self._ledger().seal()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "depth_returns.csv"
            written = write_depth_returns_csv(path, table.rows, stamp="provenance=TEST")
            with open(path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(written, len(table.rows))
        self.assertEqual(
            tuple(rows[0]), ("provenance",) + DepthReturnRow.field_names()
        )

    def test_the_coverage_accumulator_still_refuses_an_outcome(self) -> None:
        """The two ledgers are separate: the outcome-free one stays outcome-free."""
        accumulator = CoverageAccumulator(library=TRPB, eligible_ids=_ids_at_depth(4, 2))
        with self.assertRaises(TypeError):
            accumulator.record_round(
                campaign="s0",
                round_index=1,
                proposed=[],
                selected=[MeasuredSelection(_ids_at_depth(4, 1)[0], 0.7036, True)],
            )


_DATA_DIR = Path(__file__).resolve().parent.parent / "scripts" / "data"
_TRPB_CSV = _DATA_DIR / "TRPB_THEMA_Johnston_2024.csv"
_GB1_CSV = _DATA_DIR / "GB1_STRSG_Wu_2016.csv"


def _load_landscape(path: Path) -> dict[str, float]:
    table: dict[str, float] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            mutant = (row.get("mutant") or "").strip()
            raw = (row.get("DMS_score") or "").strip()
            if not mutant or not raw or mutant == "WT":
                continue
            try:
                table[mutant] = float(raw)
            except ValueError:
                continue
    return table


class TestRealLandscapeDepthShells(unittest.TestCase):
    """
    Validates the depth parsing against the two real landscapes, by reproducing the
    measured shell hit rates exactly:

        TrpB shallow (depth<=2) 7.49%, deep (depth>=3) 0.91%
        GB1  shallow (depth<=2) 3.32%, deep (depth>=3) 0.97%

    This is where the GB1-vs-TrpB reversal lives: TrpB's shallow shell is 8.3x richer
    in top-1% variants than its deep shell, GB1's only 3.4x, while the two deep
    shells are nearly identical. A shallow-biased selector therefore wins on TrpB and
    loses on GB1 without anything about the selector changing.

    These numbers reference true fitness over variants nobody measured, so they are
    RETROSPECTIVE and gated behind an access token.
    """

    def test_access_token_is_required(self) -> None:
        with self.assertRaises(RetrospectiveAccessError):
            depth_shell_profile(TRPB, fitness_lookup={"V183A": 0.1}, access=None)
        with self.assertRaises(RetrospectiveAccessError):
            depth_shell_profile(TRPB, fitness_lookup={"V183A": 0.1}, access=True)

    @unittest.skipUnless(_TRPB_CSV.is_file(), "TrpB landscape CSV not present")
    def test_trpb_shell_hit_rates(self) -> None:
        profile = depth_shell_profile(
            TRPB,
            fitness_lookup=_load_landscape(_TRPB_CSV),
            access=RetrospectiveFitnessAccess(reason="validating depth parsing on TrpB"),
        )
        self.assertEqual(profile.n_scored, UNIVERSE_SIZE)
        self.assertEqual(profile.unparseable, 0)
        self.assertEqual(profile.per_depth[4]["n"], float(UNIVERSE_DEPTH_4))
        self.assertAlmostEqual(profile.shallow["hit_rate"] * 100, 7.49, places=2)
        self.assertAlmostEqual(profile.deep["hit_rate"] * 100, 0.91, places=2)
        # The global optimum is a depth-4 genotype: the deep shell is where the answer is.
        self.assertAlmostEqual(profile.per_depth[4]["max"], 1.0, places=6)

    @unittest.skipUnless(_GB1_CSV.is_file(), "GB1 landscape CSV not present")
    def test_gb1_shell_hit_rates(self) -> None:
        table = _load_landscape(_GB1_CSV)
        library = infer_library_from_variant_ids(table)
        self.assertEqual(library.positions, (39, 40, 41, 54))
        self.assertEqual(library.parent_residues, "VDGV")
        profile = depth_shell_profile(
            library,
            fitness_lookup=table,
            access=RetrospectiveFitnessAccess(reason="validating depth parsing on GB1"),
        )
        self.assertEqual(profile.unparseable, 0)
        self.assertAlmostEqual(profile.shallow["hit_rate"] * 100, 3.32, places=2)
        self.assertAlmostEqual(profile.deep["hit_rate"] * 100, 0.97, places=2)


class TestMathWorkerTraceSink(unittest.TestCase):
    """The production selector accepts the sink and fills it per candidate."""

    def test_select_next_micro_batch_populates_the_sink(self) -> None:
        from controller.prompts import ResidueCluster, SearchStrategy
        from workers.math_worker import select_next_micro_batch
        from workers.rag_worker import fetch_protein_metadata

        metadata = fetch_protein_metadata("P42212", use_cache=True)
        strategy = SearchStrategy(
            target_protein="avGFP",
            target_accession="P42212",
            hypothesis="coverage trace sink",
            target_clusters=[
                ResidueCluster(cluster_id="c1", positions=[64, 65, 68],
                               rationale="probe"),
            ],
            forbidden_positions=[],
            max_mutation_depth=2,
        )
        trace = BoundedCandidateTrace(
            campaign="seed0", budget=TraceBudget(top_n=5, rejected_sample=5, seed=0)
        )
        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=metadata.sequence,
            evidence=None,
            batch_size=10,
            random_seed=0,
            candidate_trace=trace,
        )
        self.assertEqual(len(batch.variants), 10)
        self.assertTrue(trace.records)
        selected = {r.variant_id for r in trace.records if r.status == "selected"}
        self.assertEqual(selected, {v.notation for v in batch.variants})
        self.assertTrue(
            {r.generation_source for r in trace.records}
            <= {"combinatorial", "gflownet"}
        )
        for record in trace.records:
            self.assertTrue(math.isfinite(record.predicted_fitness))
            self.assertTrue(math.isfinite(record.posterior_uncertainty))
            self.assertTrue(math.isfinite(record.prior_contribution))
            self.assertTrue(math.isfinite(record.acquisition_score))


if __name__ == "__main__":
    unittest.main()
