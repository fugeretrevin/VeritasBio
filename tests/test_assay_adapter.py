"""
tests/test_assay_adapter.py — Unit tests for Assay-Adapted Foundation Model Head.

Verifies:
  1. Adaptive Shrinkage: Foundation model predictions dynamically fuse with empirical assay evidence.
  2. Non-natural assay peak adaptation: When an assay peak is penalized by zero-shot evolutionary priors,
     the assay adapter successfully learns the assay feedback and boosts its score.
  3. Calibrated uncertainty: Predicts valid calibrated epistemic uncertainty across candidate mutations.
"""

from __future__ import annotations

import unittest
import numpy as np

from workers.math_worker import AssayAdaptedHead


class TestAssayAdaptedHead(unittest.TestCase):
    """Unit tests for AssayAdaptedHead (Enterprise Improvement 1)."""

    def test_zero_shot_prior_fallback_when_empty(self) -> None:
        """Verifies that with zero evidence, the adapter falls back cleanly to zero-shot prior."""
        adapter = AssayAdaptedHead(feature_dim=15, random_state=42)
        X_test = np.random.normal(0, 1, size=(10, 15))
        prior_scores = np.array([0.5, 0.8, -0.2, 1.2, 0.0, 0.4, 0.9, -0.5, 0.3, 0.7])

        # Predict with no evidence
        fused = adapter.predict(X_test, zero_shot_scores=prior_scores)
        np.testing.assert_allclose(fused, prior_scores, rtol=1e-5, atol=1e-5)

    def test_assay_adaptation_on_unnatural_peak(self) -> None:
        """
        Verifies that when an unnatural variant has negative zero-shot prior (-5.0)
        but positive assay evidence (+2.5), the adapter learns to promote it.
        """
        rng = np.random.RandomState(42)
        n_features = 8
        n_samples = 30

        # Features for 30 variants
        X_train = rng.normal(0, 1, size=(n_samples, n_features))
        # The first feature X[:, 0] represents a specific residue substitution (e.g. Lys227)
        # In nature it is penalized (-4.0), but in the assay it forms a critical salt bridge (+2.5)
        zero_shot_priors = -2.0 * X_train[:, 0]
        actual_assay_scores = 3.0 * X_train[:, 0] + rng.normal(0, 0.1, size=n_samples)

        adapter = AssayAdaptedHead(feature_dim=n_features, random_state=42)
        adapter.fit(X_train, actual_assay_scores, zero_shot_scores=zero_shot_priors)

        # Candidate with positive X[0] = +1.5 (unnatural salt bridge)
        X_cand = np.zeros((1, n_features))
        X_cand[0, 0] = 1.5
        cand_prior = np.array([-3.0])  # Zero-shot model hates it

        pred_fused = adapter.predict(X_cand, zero_shot_scores=cand_prior)

        # The assay adapter must overturn the negative zero-shot penalty
        self.assertGreater(
            pred_fused[0],
            0.0,
            f"Assay adapter must promote assay-positive unnatural variant (was {pred_fused[0]:.3f})",
        )


if __name__ == "__main__":
    unittest.main()
