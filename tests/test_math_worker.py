"""
tests/test_math_worker.py — Unit test suite for workers.math_worker (Quant Engine).

Tests:
    1. Physicochemical property featurization (delta vectors, dimension, canonical AAs).
    2. Combinatorial candidate generation within clusters (forbidden filtering, depth limits).
    3. Gaussian Process fitting and Expected Improvement (EI) calculation.
    4. Batch selection sizing and integration with Gatekeeper schemas.
"""

import math
import unittest
import numpy as np

from controller.gatekeeper import Variant
from controller.prompts import ResidueCluster, SearchStrategy
from workers.math_worker import (
    AA_FEATURE_DIM,
    compute_expected_improvement,
    featurize_single_delta,
    featurize_variant_on_positions,
    generate_candidate_variants,
    select_next_micro_batch,
)
from workers.rag_worker import fetch_protein_metadata


class TestMathWorker(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.metadata = fetch_protein_metadata("P42212", use_cache=True)
        cls.wt_sequence = cls.metadata.sequence
        cls.strategy = SearchStrategy(
            target_protein="avGFP",
            target_accession="P42212",
            hypothesis="Target chromophore core and interacting pocket",
            target_clusters=[
                ResidueCluster(
                    cluster_id="chromophore",
                    positions=[64, 65, 66, 67, 68],
                    rationale="Core fluorophore and immediate environment",
                ),
                ResidueCluster(
                    cluster_id="pocket",
                    positions=[203, 222],
                    rationale="Fluorescence shift and stability contacts",
                ),
            ],
            forbidden_positions=[1, 2, 238],
            max_mutation_depth=2,
        )

    # -------------------------------------------------------------------------
    # 1. Featurization Tests
    # -------------------------------------------------------------------------
    def test_featurize_single_delta_shape(self):
        delta = featurize_single_delta("S", "T")
        self.assertEqual(len(delta), AA_FEATURE_DIM)
        self.assertTrue(np.all(np.isfinite(delta)))

    def test_featurize_variant_vector_length(self):
        # 5 active positions: 64, 65, 66, 67, 68 -> vector length = 5 * 5 = 25
        active_positions = [64, 65, 66, 67, 68]
        v = Variant(notation="S65T")
        vec = featurize_variant_on_positions(v, active_positions, self.wt_sequence)
        self.assertEqual(len(vec), 25)
        self.assertTrue(np.all(np.isfinite(vec)))

        # Position 64 is unmutated -> first 5 entries should be 0.0
        np.testing.assert_array_equal(vec[0:5], np.zeros(5))

        # Position 65 is S->T -> entries 5:10 should be non-zero
        self.assertFalse(np.all(vec[5:10] == 0.0))

    # -------------------------------------------------------------------------
    # 2. Candidate Generation Tests
    # -------------------------------------------------------------------------
    def test_generate_candidates_forbidden_exclusion(self):
        candidates = generate_candidate_variants(self.strategy, self.wt_sequence)
        self.assertGreater(len(candidates), 50)

        # Ensure NO candidate mutates any of the forbidden positions [1, 2, 238]
        forbidden_set = {1, 2, 238}
        for cand in candidates:
            cand_positions = {m.position for m in cand.mutations}
            self.assertTrue(cand_positions.isdisjoint(forbidden_set))

    def test_combinatorial_double_mutants_generated(self):
        candidates = generate_candidate_variants(self.strategy, self.wt_sequence)
        double_mutants = [c for c in candidates if len(c.mutations) == 2]
        self.assertGreater(len(double_mutants), 0)
        # Verify valid notation format
        self.assertIn(":", double_mutants[0].notation)

    # -------------------------------------------------------------------------
    # 3. Expected Improvement & GP Active Learning Selection
    # -------------------------------------------------------------------------
    def test_compute_expected_improvement_vectorized(self):
        means = np.array([1.5, 0.5, 0.0])
        stds = np.array([0.5, 0.8, 0.1])
        f_best = 1.0
        eis = compute_expected_improvement(means, stds, f_best=f_best, xi=0.01)

        self.assertEqual(len(eis), 3)
        self.assertTrue(np.all(eis >= 0.0))
        # Candidate 1 (mean 1.5 > f_best 1.0) should have the highest EI
        self.assertGreater(eis[0], eis[1])
        self.assertGreater(eis[0], eis[2])

    def test_select_next_micro_batch_cold_start(self):
        # Round 1: No evidence yet
        batch = select_next_micro_batch(
            strategy=self.strategy,
            wt_sequence=self.wt_sequence,
            evidence=None,
            batch_size=50,
        )

        self.assertEqual(len(batch.variants), 50)
        self.assertEqual(len(batch.expected_improvements), 50)
        self.assertEqual(len(batch.predicted_means), 50)
        self.assertEqual(len(batch.predicted_uncertainties), 50)

        # All variants should pass biophysical validation against WT
        for v in batch.variants:
            v.validate_against_wt(self.wt_sequence)

    def test_select_next_micro_batch_with_evidence(self):
        # Round 2: Historical evidence present
        sample_evidence = [
            {"variant_sequence": "S65T", "computed_score": 0.85},
            {"variant_sequence": "Y66H", "computed_score": -0.40},
            {"variant_sequence": "T203Y", "computed_score": 1.45},  # High-performer
        ]
        # In evidence, variant_sequence can be full sequence or mutation notation
        # Here we test providing the mutated sequence
        full_seq_high = Variant(notation="T203Y").apply_to(self.wt_sequence)
        full_seq_low = Variant(notation="Y66H").apply_to(self.wt_sequence)
        evidence_full = [
            {"variant_sequence": full_seq_high, "computed_score": 1.45},
            {"variant_sequence": full_seq_low, "computed_score": -0.40},
        ]

        batch = select_next_micro_batch(
            strategy=self.strategy,
            wt_sequence=self.wt_sequence,
            evidence=evidence_full,
            batch_size=25,
        )

        self.assertEqual(len(batch.variants), 25)
        # All selected EIs must be non-negative
        for ei in batch.expected_improvements:
            self.assertGreaterEqual(ei, 0.0)


if __name__ == "__main__":
    unittest.main()
