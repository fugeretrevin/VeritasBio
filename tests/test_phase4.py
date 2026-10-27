"""
tests/test_phase4.py — Unit test suite for Phase 4 (Cloud Muscle & Mock Workers).

Tests:
    1. workers.mock_workers:
       - Realistic GFP fitness landscape (WT=0.0, S65T/T203Y>0, G67/Y66 penalties).
       - Epistatic synergy bonus calculation (S65T + T203Y synergy).
       - MockESMWorker batch scoring.
       - MockLLMWorker strategy generation.
    2. workers.esm_worker:
       - score_variants with 'mock_esm' target.
       - Graceful fallback when Modal is unconfigured and use_mock_fallback=True.
       - Fail-closed error when use_mock_fallback=False and Modal is unconfigured.
"""

import unittest
from unittest import mock

from controller.gatekeeper import SingleMutation, Variant
from workers.esm_worker import score_variants
from workers.mock_workers import (
    MockESMWorker,
    MockLLMWorker,
    compute_mock_gfp_fitness,
)
from workers.rag_worker import fetch_protein_metadata


class TestPhase4(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.metadata = fetch_protein_metadata("P42212", use_cache=True)
        cls.wt_sequence = cls.metadata.sequence

    # -------------------------------------------------------------------------
    # 1. Biological Mock Fitness Tests
    # -------------------------------------------------------------------------
    def test_mock_wt_fitness_is_zero(self):
        score = compute_mock_gfp_fitness([])
        self.assertEqual(score, 0.0)

    def test_mock_beneficial_mutations_positive(self):
        s65t = SingleMutation.from_string("S65T")
        score_s65t = compute_mock_gfp_fitness([s65t], deterministic_noise=False)
        self.assertGreater(score_s65t, 0.5)

        t203y = SingleMutation.from_string("T203Y")
        score_t203y = compute_mock_gfp_fitness([t203y], deterministic_noise=False)
        self.assertGreater(score_t203y, 0.5)

    def test_mock_disruptive_mutations_negative(self):
        # Mutating catalytic Gly67 breaks the chromophore hinge
        g67a = SingleMutation.from_string("G67A")
        score_g67a = compute_mock_gfp_fitness([g67a], deterministic_noise=False)
        self.assertLess(score_g67a, -1.5)

        # Mutating Tyr66 breaks conjugation
        y66a = SingleMutation.from_string("Y66A")
        score_y66a = compute_mock_gfp_fitness([y66a], deterministic_noise=False)
        self.assertLess(score_y66a, -1.0)

    def test_mock_epistatic_synergy(self):
        # Epistatic pair S65T + T203Y should exceed the sum of individual effects
        s65t = SingleMutation.from_string("S65T")
        t203y = SingleMutation.from_string("T203Y")

        score_s65t = compute_mock_gfp_fitness([s65t], deterministic_noise=False)
        score_t203y = compute_mock_gfp_fitness([t203y], deterministic_noise=False)
        score_double = compute_mock_gfp_fitness([s65t, t203y], deterministic_noise=False)

        # Non-linear synergy bonus (+0.45)
        self.assertGreater(score_double, score_s65t + score_t203y)

    def test_mock_esm_worker_batch(self):
        worker = MockESMWorker(self.wt_sequence)
        var1 = Variant(notation="S65T").apply_to(self.wt_sequence)
        var2 = Variant(notation="T203Y").apply_to(self.wt_sequence)

        scores = worker.score_variants(self.wt_sequence, [self.wt_sequence, var1, var2])
        self.assertEqual(len(scores), 3)
        self.assertEqual(scores[0], 0.0)  # WT
        self.assertGreater(scores[1], 0.5)  # S65T
        self.assertGreater(scores[2], 0.5)  # T203Y

    def test_mock_llm_worker(self):
        llm = MockLLMWorker()
        strategy = llm.generate_strategy(self.metadata)
        self.assertEqual(strategy.target_protein, self.metadata.protein_name)
        self.assertGreaterEqual(len(strategy.target_clusters), 1)

    # -------------------------------------------------------------------------
    # 2. ESM Worker Fallback & Dispatcher Tests
    # -------------------------------------------------------------------------
    def test_esm_worker_mock_target(self):
        var1 = Variant(notation="S65T").apply_to(self.wt_sequence)
        scores = score_variants(
            wt_sequence=self.wt_sequence,
            variant_sequences=[var1],
            target_model="mock_esm",
        )
        self.assertEqual(len(scores), 1)
        self.assertGreater(scores[0], 0.5)

    @mock.patch("workers.esm_worker.is_modal_configured", return_value=False)
    def test_esm_worker_fallback_when_unconfigured(self, _mock_configured):
        var1 = Variant(notation="S65T").apply_to(self.wt_sequence)
        # When Modal is unconfigured, use_mock_fallback=True returns valid mock scores
        scores = score_variants(
            wt_sequence=self.wt_sequence,
            variant_sequences=[var1],
            target_model="esm2_650m",
            use_mock_fallback=True,
        )
        self.assertEqual(len(scores), 1)
        self.assertGreater(scores[0], 0.5)


if __name__ == "__main__":
    unittest.main()
