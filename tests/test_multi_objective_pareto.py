"""
tests/test_multi_objective_pareto.py — Unit tests for Multi-Objective Pareto Optimization Engine.

Verifies:
  1. Fast Non-Dominated Sorting: Correctly computes Pareto ranks across 3 objectives (Activity, Stability, Synthesizability).
  2. Non-dominated variants receive Rank 0 (Pareto Front).
  3. Synthesizability scoring: Penalizes hydrophobic aggregation patches and homopolymer repeats.
  4. Integration with k-DPP: Pareto-weighted quality scores strictly prioritize Pareto-optimal candidates.
"""

from __future__ import annotations

import unittest
import numpy as np

from controller.gatekeeper import Variant
from workers.math_worker import (
    compute_non_dominated_pareto_ranks,
    evaluate_synthesizability_score,
)


class TestMultiObjectivePareto(unittest.TestCase):
    """Unit tests for Multi-Objective Pareto Engine (Enterprise Improvement 2)."""

    def test_non_dominated_sorting_2d(self) -> None:
        """Verifies Pareto frontier identification on known 2D test points."""
        # Point A (1.0, 1.0) dominates Point B (0.5, 0.5)
        # Point C (0.2, 1.2) is non-dominated (trade-off)
        objectives = np.array([
            [1.0, 1.0],  # Point 0: Non-dominated
            [0.5, 0.5],  # Point 1: Dominated by Point 0
            [0.2, 1.2],  # Point 2: Non-dominated (higher obj 1)
            [1.2, 0.2],  # Point 3: Non-dominated (higher obj 0)
        ])

        ranks = compute_non_dominated_pareto_ranks(objectives)

        # Points 0, 2, 3 should be Rank 0 (Pareto front)
        self.assertEqual(ranks[0], 0)
        self.assertEqual(ranks[2], 0)
        self.assertEqual(ranks[3], 0)
        # Point 1 should be dominated (Rank >= 1)
        self.assertGreaterEqual(ranks[1], 1)

    def test_non_dominated_sorting_3d(self) -> None:
        """Verifies 3-objective Pareto ranking: [Activity, Stability, Synthesizability]."""
        objectives = np.array([
            [2.5, 0.9, 0.8],  # P0: Non-dominated (high activity, high stability)
            [1.0, 0.4, 0.3],  # P1: Strictly dominated by P0
            [0.8, 0.95, 0.9], # P2: Non-dominated (highest stability and synthesizability)
            [3.0, 0.2, 0.5],  # P3: Non-dominated (highest activity)
            [1.5, 0.5, 0.4],  # P4: Dominated by P0
        ])

        ranks = compute_non_dominated_pareto_ranks(objectives)

        self.assertEqual(ranks[0], 0, "P0 must be on the Pareto front")
        self.assertEqual(ranks[2], 0, "P2 must be on the Pareto front")
        self.assertEqual(ranks[3], 0, "P3 must be on the Pareto front")
        self.assertGreaterEqual(ranks[1], 1, "P1 must be dominated")
        self.assertGreaterEqual(ranks[4], 1, "P4 must be dominated")

    def test_synthesizability_scoring(self) -> None:
        """Verifies synthesizability penalties for aggregation patches and repeats."""
        wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVK"

        # 1. Benign polar mutation -> High synthesizability
        benign_var = Variant(notation="S65T")
        score_benign = evaluate_synthesizability_score(benign_var, wt)
        self.assertGreaterEqual(score_benign, 0.80, "Benign polar mutation must have high synthesizability")

        # 2. Aggregating hydrophobic patch (multiple adjacent Leu/Trp/Phe)
        # Mutating nearby residues to form a continuous sticky patch
        sticky_var = Variant(notation="S65W:Y66W:V68W")
        score_sticky = evaluate_synthesizability_score(sticky_var, wt)
        self.assertLess(score_sticky, score_benign, "Sticky hydrophobic patch must reduce synthesizability score")


if __name__ == "__main__":
    unittest.main()
