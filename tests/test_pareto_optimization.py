"""
tests/test_pareto_optimization.py — Test suite for Multi-Objective Pareto Optimization (Track B).

Tests Seam 1:
- compute_multi_objective_pareto_front
- ParetoPoint and ParetoAnalysisResult data contracts
- Fast non-dominated sorting (NSGA-II) trade-off surface
- Dynamic multi-objective weighting (Fitness vs. Stability vs. Synthesizability vs. Developability)
- Integration with select_next_micro_batch
"""

from __future__ import annotations

import unittest
import numpy as np

from controller.gatekeeper import SingleMutation, Variant
from controller.prompts import ResidueCluster, SearchStrategy
from workers.math_worker import (
    ParetoPoint,
    ParetoAnalysisResult,
    compute_multi_objective_pareto_front,
    select_next_micro_batch,
)


class TestParetoOptimization(unittest.TestCase):
    """Test suite for multi-objective Pareto optimization and developability trade-offs."""

    def setUp(self):
        # Model antibody heavy chain variable fragment (VH)
        self.wt_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSVIYSGGSYTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDGYYYGMDVWGQGTTVTVSS"

        # Construct 4 variants with known conflicting tradeoffs
        # Variant 1: High fitness, severe developability liabilities (NG deamidation at 25-26 + free Cys at 24)
        self.v_high_fit_bad_dev = Variant(
            notation="A24C:S25N",
            mutations=[
                SingleMutation(wt_aa="A", position=24, mut_aa="C"),  # Odd cysteine (3 Cys total)
                SingleMutation(wt_aa="S", position=25, mut_aa="N"),  # Creates NG at 25-26
            ],
        )

        # Variant 2: Moderate fitness, clean developability, zero liabilities
        self.v_balanced = Variant(
            notation="A24V",
            mutations=[
                SingleMutation(wt_aa="A", position=24, mut_aa="V"),
            ],
        )

        # Variant 3: Low fitness, poor synthesizability (homopolymer repeat GGG at 25-27)
        self.v_dominated = Variant(
            notation="S25G:F27G",
            mutations=[
                SingleMutation(wt_aa="S", position=25, mut_aa="G"),
                SingleMutation(wt_aa="F", position=27, mut_aa="G"),  # creates GGG homopolymer (25=G, 26=G, 27=G)
            ],
        )

        # Variant 4: High stability & synthesizability, moderate fitness
        self.v_high_stab = Variant(
            notation="V5L",
            mutations=[
                SingleMutation(wt_aa="V", position=5, mut_aa="L"),
            ],
        )

        self.candidates = [
            self.v_high_fit_bad_dev,
            self.v_balanced,
            self.v_dominated,
            self.v_high_stab,
        ]

    def test_compute_multi_objective_pareto_front_contracts_and_tradeoffs(self):
        """Verifies Pareto front extraction, NSGA-II non-dominated sorting, and utility calculation."""
        predicted_fitness = np.array([2.50, 1.80, 0.30, 1.65], dtype=np.float64)
        stability_priors = np.array([0.10, 0.40, -0.50, 0.90], dtype=np.float64)

        # Default balanced weights: 40% fitness, 30% developability, 15% synthesizability, 15% stability
        result: ParetoAnalysisResult = compute_multi_objective_pareto_front(
            candidates=self.candidates,
            wt_sequence=self.wt_seq,
            predicted_fitness=predicted_fitness,
            stability_priors=stability_priors,
        )

        # 1. Structure check
        self.assertIsInstance(result, ParetoAnalysisResult)
        self.assertEqual(len(result.all_points), 4)
        self.assertTrue(result.rank_0_count >= 2)
        self.assertTrue(len(result.frontier_points) == result.rank_0_count)

        # 2. Check points metadata
        p_bad_dev = next(p for p in result.all_points if p.notation == self.v_high_fit_bad_dev.notation)
        p_balanced = next(p for p in result.all_points if p.notation == self.v_balanced.notation)
        p_dominated = next(p for p in result.all_points if p.notation == self.v_dominated.notation)

        # Variant 1 must have liabilities detected
        self.assertTrue(p_bad_dev.developability_score < p_balanced.developability_score)
        self.assertTrue(len(p_bad_dev.liabilities.get("deamidation_sites", [])) >= 1 or p_bad_dev.liabilities.get("unpaired_cysteine", False))

        # Variant 3 should be dominated by Variant 2 and 4 (higher rank index means dominated)
        self.assertTrue(p_dominated.rank > 0)

        # Both Variant 1 (supreme fitness) and Variant 2 (high developability) should be on the Rank 0 front
        self.assertEqual(p_bad_dev.rank, 0)
        self.assertEqual(p_balanced.rank, 0)

        # 3. Dynamic utility weights shifting
        # Developability-heavy regime
        dev_result = compute_multi_objective_pareto_front(
            candidates=self.candidates,
            wt_sequence=self.wt_seq,
            predicted_fitness=predicted_fitness,
            stability_priors=stability_priors,
            weights={"fitness": 0.1, "developability": 0.7, "synthesizability": 0.1, "stability": 0.1},
        )
        dev_bad = next(p for p in dev_result.all_points if p.notation == self.v_high_fit_bad_dev.notation)
        dev_bal = next(p for p in dev_result.all_points if p.notation == self.v_balanced.notation)
        self.assertTrue(dev_bal.weighted_utility > dev_bad.weighted_utility)

        # Fitness-heavy regime
        fit_result = compute_multi_objective_pareto_front(
            candidates=self.candidates,
            wt_sequence=self.wt_seq,
            predicted_fitness=predicted_fitness,
            stability_priors=stability_priors,
            weights={"fitness": 0.8, "developability": 0.1, "synthesizability": 0.05, "stability": 0.05},
        )
        fit_bad = next(p for p in fit_result.all_points if p.notation == self.v_high_fit_bad_dev.notation)
        fit_bal = next(p for p in fit_result.all_points if p.notation == self.v_balanced.notation)
        self.assertTrue(fit_bad.weighted_utility > fit_bal.weighted_utility)

    def test_select_next_micro_batch_with_custom_pareto_weights(self):
        """Verifies that select_next_micro_batch accepts and respects custom pareto_weights_config."""
        strategy = SearchStrategy(
            target_protein="TestVH",
            target_accession="CUSTOM_VH",
            hypothesis="Test custom Pareto weights",
            target_clusters=[
                ResidueCluster(cluster_id="cdr1", positions=[24, 25], rationale="CDR-H1"),
            ],
            forbidden_positions=[],
            max_mutation_depth=1,
        )

        # Pass custom pareto weights strongly penalizing poor developability
        custom_weights = {
            "fitness": 0.2,
            "developability": 0.6,
            "synthesizability": 0.1,
            "stability": 0.1,
        }

        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=self.wt_seq,
            evidence=[],
            batch_size=5,
            pareto_weights_config=custom_weights,
        )

        self.assertEqual(len(batch.variants), 5)
        self.assertEqual(len(batch.pareto_ranks), 5)
        self.assertTrue(all(r >= 0 for r in batch.pareto_ranks))


if __name__ == "__main__":
    unittest.main()
