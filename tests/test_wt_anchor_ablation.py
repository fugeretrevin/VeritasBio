"""The synthetic wild-type training row must be omittable, and omitting it must be the
only thing that changes.

Why this exists. `select_next_micro_batch` injects one training row per round for the
wild type with the literal label `y = 0.0`, whatever the assay's scale is, at variance
1e-4 while a real measurement defaults to 1e-2. The distance from 0.0 to the centre of
each landscape's real labels, in label standard deviations, is:

    creilov 15.68 | cr9114_h1 14.14 | cr6261_h1 11.00 | cr6261_h9 8.35
    eqfp611_blue 1.19 | eqfp611_red 0.80 | her2 0.56 | trpb 0.35 | gb1 0.20

TrpB and GB1 - the two landscapes every quotable engine figure comes from - really are
centred near zero, so the row is nearly true there and the defect is invisible. On
CreiLOV the real labels occupy 2.793..4.196 and the row sits 11.6 units below the
minimum, which is 8.3x the whole label range.

Measured consequence, CreiLOV, 12 paired seeds, ridge held constant at 59.17:
omitting this one row moves the engine from 9.08 to 38.83 top-1% hits
(+29.75, t = 12.17, better on 12 of 12 seeds) - 59% of a 50.08-hit deficit.

These tests assert the MECHANISM: that the switch removes exactly one row, that the row
it removes is the synthetic one, and that the default path still adds it. They do not
assert a yield. A test asserting "the ablated arm scores 38.83" would pin one landscape's
number into the suite, which is the error this repository has made three times.
"""
from __future__ import annotations

import unittest
from unittest import mock

import numpy as np

from controller.gatekeeper import Variant
from workers import math_worker


class TestWtAnchorIsAblatable(unittest.TestCase):
    def test_it_is_in_the_ablation_vocabulary(self) -> None:
        self.assertIn("wt_anchor", math_worker.ABLATABLE_COMPONENTS)

    def test_the_benchmark_accepts_the_name(self) -> None:
        import benchmark_panoptes_isolated as bench

        self.assertIn("wt_anchor", bench.ABLATABLE_COMPONENTS)
        self.assertEqual(bench.parse_ablate_components("wt_anchor"), ("wt_anchor",))

    # ---- the mechanism ---------------------------------------------------------------
    def _fit_calls(self, ablate):
        """Run one real round and capture every (X, y) the surrogate was fitted on."""
        from controller.prompts import ResidueCluster, SearchStrategy

        wt = "M" + "ACDEFGHIKLMNPQRSTVWY" * 12
        positions = [10, 20, 30, 40]
        strategy = SearchStrategy(
            target_protein="synthetic", target_accession="TEST",
            hypothesis="exercise the training-set assembly",
            target_clusters=[ResidueCluster(
                cluster_id="sites", positions=positions,
                rationale="a synthetic four-site library")],
            max_mutation_depth=3,
        )
        pool = [Variant(notation=n) for n in (
            f"{wt[9]}10A", f"{wt[19]}20C", f"{wt[9]}10A:{wt[19]}20C", f"{wt[29]}30D")]
        # Labels deliberately far from zero, which is the case the component is about.
        evidence = [
            {"variant_sequence": Variant(notation=f"{wt[9]}10G").apply_to(wt),
             "computed_score": 3.90},
            {"variant_sequence": Variant(notation=f"{wt[19]}20G").apply_to(wt),
             "computed_score": 4.05},
            {"variant_sequence": Variant(notation=f"{wt[29]}30G").apply_to(wt),
             "computed_score": 3.80},
        ]
        seen = []
        real_build = math_worker.build_gp_model

        def spy_build(*a, **kw):
            model = real_build(*a, **kw)
            real_fit = model.fit

            def fit(X, y, *fa, **fkw):
                seen.append((np.asarray(X), np.asarray(y)))
                return real_fit(X, y, *fa, **fkw)

            model.fit = fit
            return model

        with mock.patch.object(math_worker, "build_gp_model", side_effect=spy_build):
            np.random.seed(5)
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=evidence, batch_size=2,
                candidate_pool=pool, ablate=ablate, random_seed=3,
            )
        self.assertTrue(seen, "the surrogate was never fitted, so this test proves nothing")
        return seen[0]

    def test_the_default_path_trains_on_a_zero_labelled_row(self) -> None:
        """Guards the test below against passing vacuously, and documents the behaviour
        being made switchable: with evidence at 3.80-4.05, a 0.0 label is present."""
        _X, y = self._fit_calls(())
        self.assertIn(0.0, list(y), "the default path no longer injects a y=0.0 row; the "
                                    "ablation test below is no longer measuring anything")
        self.assertEqual(len(y), 4, "expected the WT row plus three evidence rows")

    def test_ablating_it_removes_exactly_that_row(self) -> None:
        X_default, y_default = self._fit_calls(())
        X_ablated, y_ablated = self._fit_calls(("wt_anchor",))
        self.assertEqual(len(y_ablated), len(y_default) - 1,
                         "wt_anchor removed something other than exactly one row")
        self.assertNotIn(0.0, list(y_ablated),
                         "the synthetic zero-labelled row survived its own ablation")
        self.assertEqual(X_ablated.shape[1], X_default.shape[1],
                         "the feature width moved, so more than the row changed")

    def test_the_remaining_rows_are_the_measured_ones_unchanged(self) -> None:
        """Only the synthetic row goes. Every real label must survive untouched, or the
        arm would be measuring a different training set rather than one row fewer."""
        _Xd, y_default = self._fit_calls(())
        _Xa, y_ablated = self._fit_calls(("wt_anchor",))
        kept = sorted(v for v in y_default if v != 0.0)
        self.assertEqual(sorted(y_ablated.tolist()), kept)

    def test_the_sample_variance_vector_stays_aligned(self) -> None:
        """The WT row also contributed a 1e-4 variance entry. If the row goes and that
        entry stays, `sample_variances_arr` is one longer than X_train, which
        math_worker silently turns into None (:2139) - so a misalignment would not
        raise, it would quietly disable heteroscedastic noise for the whole arm.
        """
        from controller.prompts import ResidueCluster, SearchStrategy

        wt = "M" + "ACDEFGHIKLMNPQRSTVWY" * 12
        strategy = SearchStrategy(
            target_protein="synthetic", target_accession="TEST",
            hypothesis="variance alignment",
            target_clusters=[ResidueCluster(cluster_id="s", positions=[10, 20, 30, 40],
                                            rationale="synthetic")],
            max_mutation_depth=3,
        )
        pool = [Variant(notation=f"{wt[9]}10A"), Variant(notation=f"{wt[19]}20C")]
        evidence = [
            {"variant_sequence": Variant(notation=f"{wt[9]}10G").apply_to(wt),
             "computed_score": 3.9, "variance": 0.01},
            {"variant_sequence": Variant(notation=f"{wt[19]}20G").apply_to(wt),
             "computed_score": 4.0, "variance": 0.02},
        ]
        captured = {}
        real_dkl = math_worker.DeepKernelLearner

        class Spy(real_dkl):  # type: ignore[misc,valid-type]
            def fit(self, X, y, sample_variance=None, **kw):
                captured["n_rows"] = len(y)
                captured["n_var"] = None if sample_variance is None else len(sample_variance)
                return super().fit(X, y, sample_variance=sample_variance, **kw)

        with mock.patch.object(math_worker, "DeepKernelLearner", Spy), \
             mock.patch.object(math_worker, "TWO_STAGE_GP_ENABLED", True):
            np.random.seed(5)
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=evidence, batch_size=1,
                candidate_pool=pool, ablate=("wt_anchor",), random_seed=3,
            )
        self.assertIn("n_rows", captured, "the DKL surrogate was never fitted, so the "
                                          "variance alignment was not exercised")
        self.assertEqual(
            captured["n_var"], captured["n_rows"],
            f"sample_variance has {captured['n_var']} entries for "
            f"{captured['n_rows']} training rows; math_worker would silently drop it")


if __name__ == "__main__":
    unittest.main()
