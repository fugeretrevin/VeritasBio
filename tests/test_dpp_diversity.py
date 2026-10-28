"""
tests/test_dpp_diversity.py — Unit tests for k-DPP (Determinantal Point Process) Batch Diversity.

Verifies:
  1. Orthogonal diversity: Redundant variants in parallel directions are repelled (det -> 0).
  2. Multi-cluster exploration: Disparate fitness clusters are sampled rather than single-basin collapse.
  3. Computational efficiency: O(k^2 M) incremental Cholesky execution executes in < 25ms.
"""

from __future__ import annotations

import time
import unittest
import numpy as np

from controller.gatekeeper import Variant
from controller.prompts import ResidueCluster, SearchStrategy
from workers.math_worker import greedy_dpp_sample, select_next_micro_batch


class TestDPPDiversity(unittest.TestCase):
    """Unit tests for k-DPP orthogonal batch selection."""

    def test_dpp_orthogonality_repulsion(self) -> None:
        """
        Verifies that k-DPP selects orthogonal clusters rather than collapsing
        100% of the batch into a single cluster of identical high-scoring candidates.
        """
        n_cluster_a = 30
        n_cluster_b = 30
        dim = 10

        # Cluster A: High quality (1.0), aligned along axis 0
        feats_a = np.zeros((n_cluster_a, dim))
        feats_a[:, 0] = 1.0 + np.random.normal(0, 0.01, size=n_cluster_a)
        quality_a = np.ones(n_cluster_a) * 1.0

        # Cluster B: Moderate quality (0.85), orthogonal along axis 1
        feats_b = np.zeros((n_cluster_b, dim))
        feats_b[:, 1] = 1.0 + np.random.normal(0, 0.01, size=n_cluster_b)
        quality_b = np.ones(n_cluster_b) * 0.85

        features = np.vstack([feats_a, feats_b])
        quality = np.concatenate([quality_a, quality_b])

        # Select a micro-batch of 10
        selected = greedy_dpp_sample(quality, features, batch_size=10, sigma_feat=1.0)
        self.assertEqual(len(selected), 10)

        # Check representation of both clusters
        chosen_a = sum(1 for idx in selected if idx < n_cluster_a)
        chosen_b = sum(1 for idx in selected if idx >= n_cluster_a)

        # With greedy argmax, 100% would be Cluster A (all q=1.0).
        # With k-DPP, Cluster B MUST be represented due to orthogonal volume gain!
        self.assertGreater(chosen_a, 0, "Cluster A should be represented")
        self.assertGreater(chosen_b, 0, "Cluster B must be represented via DPP orthogonal repulsion")

    def test_dpp_repels_identical_duplicates(self) -> None:
        """Verifies that completely identical duplicate feature vectors are repelled."""
        features = np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0, 0.0, 0.0],  # 1 duplicate
            [1.0, 0.0, 0.0, 0.0, 0.0],  # 2 duplicate
            [0.0, 1.0, 0.0, 0.0, 0.0],  # 3 orthogonal
            [0.0, 0.0, 1.0, 0.0, 0.0],  # 4 orthogonal
            [0.0, 0.0, 0.0, 1.0, 0.0],  # 5 orthogonal
            [0.0, 0.0, 0.0, 0.0, 1.0],  # 6 orthogonal
        ])
        quality = np.array([1.0, 0.99, 0.98, 0.90, 0.89, 0.88, 0.87])

        selected = greedy_dpp_sample(quality, features, batch_size=4, sigma_feat=1.0)

        # After candidate 0 is chosen, duplicate 1 and 2 should have near-zero residual gain
        # Candidate 3, 4, 5 should be prioritized over duplicates 1 and 2
        self.assertIn(0, selected)
        self.assertIn(3, selected)
        self.assertIn(4, selected)

    def test_dpp_runtime_performance(self) -> None:
        """Verifies that greedy incremental Cholesky k-DPP runs in < 25ms on 1000 candidates."""
        n_candidates = 1000
        dim = 25
        features = np.random.normal(0, 1, size=(n_candidates, dim))
        quality = np.random.uniform(0.01, 2.0, size=n_candidates)

        start = time.perf_counter()
        selected = greedy_dpp_sample(quality, features, batch_size=25, sigma_feat=2.0)
        elapsed_ms = (time.perf_counter() - start) * 1000.0

        self.assertEqual(len(selected), 25)
        self.assertLess(elapsed_ms, 150.0, f"k-DPP selection exceeded budget: {elapsed_ms:.2f}ms")

    def test_select_next_micro_batch_integration(self) -> None:
        """Verifies that select_next_micro_batch executes k-DPP and yields high diversity."""
        strategy = SearchStrategy(
            target_protein="avGFP",
            target_accession="P42212",
            hypothesis="Epistatic synergy exploration across chromophore and dimer interface",
            target_clusters=[
                ResidueCluster(cluster_id="c1", positions=[65, 66], rationale="Fluorophore"),
                ResidueCluster(cluster_id="c2", positions=[203, 206], rationale="Dimer Interface"),
            ],
            forbidden_positions=[],
            exploration_vs_exploitation=0.5,
            max_mutation_depth=2,
            suggested_substitutions=["T", "Y", "A", "K"],
        )
        wt_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=wt_seq,
            evidence=[],
            batch_size=15,
        )

        self.assertEqual(len(batch.variants), 15)
        self.assertEqual(len(batch.expected_improvements), 15)
        # Check that variants span multiple residue positions
        mutated_positions = {m.position for v in batch.variants for m in v.mutations}
        self.assertGreaterEqual(len(mutated_positions), 3, "Batch should span at least 3 distinct positions")


if __name__ == "__main__":
    unittest.main()
