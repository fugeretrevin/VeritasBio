"""
tests/test_dkl_surrogate.py — Unit tests for Deep Kernel Learning (DKL) Surrogate Engine.

Verifies:
  1. Non-linear epistatic learning: Captures non-linear cross-terms on combinatorial landscapes.
  2. Bayesian uncertainty calibration: sigma(x) increases as points move away from training manifold.
  3. Seamless drop-in compatibility with GaussianProcessRegressor API.
"""

from __future__ import annotations

import unittest
import numpy as np

from workers.math_worker import DeepKernelLearner


class TestDeepKernelLearner(unittest.TestCase):
    """Unit tests for DeepKernelLearner (DKL)."""

    def test_dkl_fit_and_predict_shapes(self) -> None:
        """Verifies that DKL fits and predicts consistent shapes and non-negative uncertainties."""
        n_samples = 25
        n_features = 15
        n_test = 50

        X_train = np.random.normal(0, 1, size=(n_samples, n_features))
        y_train = np.sin(X_train[:, 0]) * np.cos(X_train[:, 1]) + np.random.normal(0, 0.05, size=n_samples)

        dkl = DeepKernelLearner(input_dim=n_features, latent_dim=6, random_state=42)
        dkl.fit(X_train, y_train)

        X_test = np.random.normal(0, 1, size=(n_test, n_features))
        means, stds = dkl.predict(X_test, return_std=True)

        self.assertEqual(means.shape, (n_test,))
        self.assertEqual(stds.shape, (n_test,))
        self.assertTrue(np.all(stds >= 0.0), "Epistemic uncertainties must be strictly non-negative")
        self.assertTrue(np.all(np.isfinite(means)), "Predicted means must be finite")

    def test_dkl_uncertainty_spread(self) -> None:
        """Verifies that epistemic uncertainty sigma(x) is lower on training data than unseen distant data."""
        X_train = np.array([
            [0.1, 0.2, 0.1],
            [0.2, 0.1, 0.3],
            [0.1, 0.3, 0.2],
            [0.3, 0.2, 0.1],
            [0.2, 0.2, 0.2],
        ])
        y_train = np.array([1.0, 1.2, 1.1, 1.3, 1.15])

        dkl = DeepKernelLearner(input_dim=3, latent_dim=3, random_state=42)
        dkl.fit(X_train, y_train)

        # Predict near training center vs distant extrapolation
        X_near = np.array([[0.2, 0.2, 0.2]])
        X_distant = np.array([[10.0, 15.0, 20.0]])

        _, std_near = dkl.predict(X_near, return_std=True)
        _, std_distant = dkl.predict(X_distant, return_std=True)

        self.assertLess(std_near[0], std_distant[0], "Uncertainty near training manifold must be lower than extrapolation")

    def test_dkl_non_linear_epistasis(self) -> None:
        """Verifies that DKL learns a non-linear interaction landscape (y = x0 * x1)."""
        rng = np.random.RandomState(42)
        X = rng.uniform(-1.5, 1.5, size=(40, 5))
        # Epistatic cross-term: y = x0 * x1 + 0.5 * x2
        y = X[:, 0] * X[:, 1] + 0.5 * X[:, 2]

        dkl = DeepKernelLearner(input_dim=5, latent_dim=4, random_state=42)
        dkl.fit(X, y)

        X_val = rng.uniform(-1.0, 1.0, size=(20, 5))
        y_val = X_val[:, 0] * X_val[:, 1] + 0.5 * X_val[:, 2]

        pred_means, _ = dkl.predict(X_val, return_std=True)
        correlation = np.corrcoef(pred_means, y_val)[0, 1]

        self.assertGreater(correlation, 0.70, f"DKL correlation on non-linear epistasis was {correlation:.3f}")


if __name__ == "__main__":
    unittest.main()
