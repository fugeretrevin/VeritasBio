"""
tests/test_assayable_scope.py — Keeping proposals inside what the assay can measure.

The problem this covers:
    A planner reasoning about biophysics legitimately proposes distal allosteric
    residues. On a library that saturates only a handful of sites those proposals are
    unmeasurable: they are discarded downstream as not_in_eligible_universe and the
    round's batch silently runs short. Observed on the TrpB four-site landscape, the
    live planner proposed R51/P52/L77 and two of four rounds produced zero usable
    selections.

The fix has two halves, both tested here:
  1. Inference: the measurable scope is recovered by diffing measured variants against
     wild type, so nothing is hardcoded and the scope is discovered from evidence.
  2. Enforcement: prompt guidance does not bind a model, so the scope is also applied
     deterministically to whatever the planner returns.
"""

from __future__ import annotations

import unittest

from controller.prompts import (
    BiophysicalHypothesis,
    ResidueCluster,
    SearchStrategy,
    infer_assayable_positions,
    restrict_strategy_to_assayable,
)

WT = "MKGYFGPYGGQYVPEILMGALEELEAAYEG"


def _strategy(cluster_positions, hyp_residues, chemistries, depth=4, forbidden=None):
    return SearchStrategy(
        target_protein="Test",
        target_accession="TEST",
        hypothesis="test strategy",
        target_clusters=[
            ResidueCluster(cluster_id=f"c{i}", positions=list(pos), rationale="test")
            for i, pos in enumerate(cluster_positions)
        ],
        forbidden_positions=list(forbidden or []),
        max_mutation_depth=depth,
        mechanistic_hypotheses=[
            BiophysicalHypothesis(
                pocket_name="p",
                target_residues=list(hyp_residues),
                mechanism="m",
                structural_rationale="r",
                suggested_chemistries=dict(chemistries),
                epistatic_coupling="e",
            )
        ],
    )


class TestInferAssayablePositions(unittest.TestCase):
    def test_recovers_saturated_sites_from_evidence(self) -> None:
        """Positions that vary across measured variants are the library's scope."""
        # Three variants, each substituting within positions 5, 9 and 14 only.
        def mutate(pairs):
            seq = list(WT)
            for pos, aa in pairs:
                seq[pos - 1] = aa
            return "".join(seq)

        evidence = [
            {"variant_sequence": mutate([(5, "A")]), "computed_score": 0.1},
            {"variant_sequence": mutate([(9, "L"), (14, "W")]), "computed_score": 0.4},
            {"variant_sequence": mutate([(5, "T"), (9, "V")]), "computed_score": 0.2},
        ]
        self.assertEqual(infer_assayable_positions(WT, evidence), [5, 9, 14])

    def test_no_evidence_yields_no_scope(self) -> None:
        """Round 1 has nothing to infer from, so the planner must stay unconstrained."""
        self.assertEqual(infer_assayable_positions(WT, []), [])
        self.assertEqual(infer_assayable_positions(WT, None), [])

    def test_ignores_length_changing_constructs(self) -> None:
        """An indel says nothing about substitution scope and must not widen it."""
        evidence = [
            {"variant_sequence": WT[:-3], "computed_score": 0.1},
            {"variant_sequence": WT[:4] + "A" + WT[5:], "computed_score": 0.2},
        ]
        self.assertEqual(infer_assayable_positions(WT, evidence), [5])


class TestRestrictStrategyToAssayable(unittest.TestCase):
    def test_clips_out_of_scope_positions(self) -> None:
        """Off-library positions are removed from clusters, residues and chemistries."""
        strategy = _strategy(
            cluster_positions=[[51, 52], [183, 184]],
            hyp_residues=[51, 77, 183, 184],
            chemistries={"51": ["N"], "77": ["Y"], "183": ["A"], "184": ["I"]},
        )
        clipped = restrict_strategy_to_assayable(strategy, {183, 184, 227, 228})

        kept = sorted({p for c in clipped.target_clusters for p in c.positions})
        self.assertEqual(kept, [183, 184])
        self.assertEqual(clipped.mechanistic_hypotheses[0].target_residues, [183, 184])
        self.assertEqual(
            sorted(clipped.mechanistic_hypotheses[0].suggested_chemistries),
            ["183", "184"],
        )
        # Depth cannot exceed the number of measurable positions actually retained.
        self.assertLessEqual(clipped.max_mutation_depth, 2)

    def test_empty_cluster_is_dropped_not_kept_empty(self) -> None:
        """A cluster entirely out of scope disappears rather than becoming empty."""
        strategy = _strategy(
            cluster_positions=[[51, 52], [227, 228]],
            hyp_residues=[227, 228],
            chemistries={"227": ["K"], "228": ["G"]},
        )
        clipped = restrict_strategy_to_assayable(strategy, {227, 228})
        self.assertEqual(len(clipped.target_clusters), 1)
        self.assertEqual(clipped.target_clusters[0].positions, [227, 228])

    def test_no_overlap_raises_rather_than_substituting(self) -> None:
        """
        A wholly unmeasurable strategy must fail loudly. Silently swapping in a
        different search space would misreport what the campaign explored.
        """
        strategy = _strategy(
            cluster_positions=[[51, 52]],
            hyp_residues=[51, 52],
            chemistries={"51": ["N"]},
        )
        with self.assertRaises(ValueError) as ctx:
            restrict_strategy_to_assayable(strategy, {183, 184, 227, 228})
        self.assertIn("no assayable positions", str(ctx.exception))

    def test_absent_scope_leaves_strategy_untouched(self) -> None:
        """With no scope supplied the planner stays unconstrained by design."""
        strategy = _strategy(
            cluster_positions=[[51, 52], [183]],
            hyp_residues=[51, 183],
            chemistries={"51": ["N"], "183": ["A"]},
            depth=4,
        )
        for empty in (None, set(), []):
            same = restrict_strategy_to_assayable(strategy, empty)
            self.assertEqual(
                sorted({p for c in same.target_clusters for p in c.positions}), [51, 52, 183]
            )
            self.assertEqual(same.max_mutation_depth, 4)

    def test_forbidden_positions_are_clipped_to_scope(self) -> None:
        """
        A forbidden list referring to unreachable sites is noise; drop it. Catalytic
        position 83 is outside the measurable scope, so withholding it says nothing;
        184 is in scope and must stay withheld.
        """
        strategy = _strategy(
            cluster_positions=[[183]],
            hyp_residues=[183],
            chemistries={"183": ["A"]},
            forbidden=[83, 184],
        )
        clipped = restrict_strategy_to_assayable(strategy, {183, 184})
        self.assertEqual(clipped.forbidden_positions, [184])


if __name__ == "__main__":
    unittest.main()
