"""The wild-type anchor must be REPAIRABLE, not only deletable - and the default must not move.

Three mutually exclusive behaviours, so one run folder can hold all three and they pair
within seed:

  * `zero`     - today's default: the literal `y = 0.0` at variance 1e-4, whatever the
                 assay's scale is.
  * `ablated`  - `--ablate-components wt_anchor`: no anchor row at all.
  * `measured` - the repair: the wild type's measured label at a real measurement's
                 variance.

Ablation is not the fix. The anchor exists to pin the surrogate at a known point so early
rounds are not fitting an unconstrained GP to 25 observations; deleting it may cost
accuracy wherever `y = 0.0` happens to be nearly true of the parent. The repair is to
anchor at the parent's MEASURED label, which
`scripts/analysis/wt_label_availability.py` establishes exists on all nine landscapes -
though on seven of them it is not in the prepared table and has to come from the pinned
upstream workbook or an in-repo constant.

What these tests assert is MECHANISM, never a yield. A test pinning "the measured arm
scores X on CreiLOV" would freeze one landscape's number into the suite, which is the
mistake this repository has made three times.
"""
from __future__ import annotations

import unittest
from unittest import mock

import numpy as np

from controller.gatekeeper import Variant
from workers import math_worker


def _strategy():
    from controller.prompts import ResidueCluster, SearchStrategy

    return SearchStrategy(
        target_protein="synthetic", target_accession="TEST",
        hypothesis="exercise the training-set assembly",
        target_clusters=[ResidueCluster(
            cluster_id="sites", positions=[10, 20, 30, 40],
            rationale="a synthetic four-site library")],
        max_mutation_depth=3,
    )


WT = "M" + "ACDEFGHIKLMNPQRSTVWY" * 12

#: Labels deliberately far from zero - the case the whole component is about. A fixture
#: centred on zero would let a measured-label bug pass as a rounding difference.
EVIDENCE = [
    {"variant_sequence": Variant(notation=f"{WT[9]}10G").apply_to(WT),
     "computed_score": 3.90},
    {"variant_sequence": Variant(notation=f"{WT[19]}20G").apply_to(WT),
     "computed_score": 4.05},
    {"variant_sequence": Variant(notation=f"{WT[29]}30G").apply_to(WT),
     "computed_score": 3.80},
]

POOL = [Variant(notation=n) for n in (
    f"{WT[9]}10A", f"{WT[19]}20C", f"{WT[9]}10A:{WT[19]}20C", f"{WT[29]}30D")]


def _run(**kwargs):
    """One real round. Returns (the (X, y) the GP was fitted on, the returned batch)."""
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
        batch = math_worker.select_next_micro_batch(
            strategy=_strategy(), wt_sequence=WT, evidence=EVIDENCE, batch_size=2,
            candidate_pool=POOL, random_seed=3, **kwargs,
        )
    assert seen, "the surrogate was never fitted, so this test proves nothing"
    return seen[0], batch


class TestTheVocabulary(unittest.TestCase):
    def test_three_modes_and_only_three(self) -> None:
        self.assertEqual(set(math_worker.WT_ANCHOR_MODES),
                         {"zero", "ablated", "measured"})

    def test_the_default_mode_is_still_zero(self) -> None:
        """A default flip re-bases every number in this repository. It is pinned here so
        the one-line edit that would do it silently fails CI instead."""
        self.assertEqual(math_worker.DEFAULT_WT_ANCHOR_MODE, "zero")

    def test_the_measured_variance_is_a_real_readings_variance(self) -> None:
        """1e-4 asserts the fabricated row is 100x more certain than a real unreplicated
        assay reading, which this module puts at 1e-2. Parity, not a tuning knob."""
        self.assertEqual(math_worker.WT_ANCHOR_MEASURED_VARIANCE, 1e-2)
        self.assertGreater(math_worker.WT_ANCHOR_MEASURED_VARIANCE, 1e-4)


class TestTheDefaultPathIsUntouched(unittest.TestCase):
    """The strongest requirement: with no flag, nothing may change."""

    def test_omitting_every_new_argument_is_identical_to_passing_mode_zero(self) -> None:
        (X_omitted, y_omitted), b1 = _run()
        (X_zero, y_zero), b2 = _run(wt_anchor_mode="zero")
        np.testing.assert_array_equal(X_omitted, X_zero)
        np.testing.assert_array_equal(y_omitted, y_zero)
        self.assertEqual([v.notation for v in b1.variants],
                         [v.notation for v in b2.variants])
        self.assertEqual(b1.expected_improvements, b2.expected_improvements)

    def test_the_default_still_injects_the_fabricated_zero_row(self) -> None:
        """Guards every test below against passing vacuously. If the default stops
        injecting y=0.0, there is no defect left to repair and these tests measure
        nothing."""
        (_X, y), batch = _run()
        self.assertIn(0.0, list(y))
        self.assertEqual(len(y), 4, "expected the anchor row plus three evidence rows")
        self.assertEqual(batch.wt_anchor_state["mode"], "zero")
        self.assertEqual(batch.wt_anchor_state["label"], 0.0)
        self.assertEqual(batch.wt_anchor_state["variance"], 1e-4)


class TestTheMeasuredMode(unittest.TestCase):
    def test_it_replaces_the_zero_with_the_supplied_label(self) -> None:
        (_Xd, y_default), _ = _run()
        (_Xm, y_measured), batch = _run(
            wt_anchor_mode="measured", wt_anchor_label=4.0232)
        self.assertEqual(len(y_measured), len(y_default),
                         "the measured mode changed the ROW COUNT; it must change only "
                         "the label and variance of the one anchor row")
        self.assertNotIn(0.0, list(y_measured),
                         "the fabricated 0.0 survived the repair")
        self.assertEqual(batch.wt_anchor_state["mode"], "measured")
        self.assertAlmostEqual(batch.wt_anchor_state["label"], 4.0232)

    def test_the_measured_rows_are_left_alone(self) -> None:
        """Only the anchor row moves. Every real label must survive untouched, or the arm
        is measuring a different training set rather than a different anchor."""
        (_Xd, y_default), _ = _run()
        (_Xm, y_measured), _ = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232)
        kept = sorted(v for v in y_default if v != 0.0)
        self.assertEqual(sorted(v for v in y_measured.tolist() if v != 4.0232), kept)

    def test_the_anchor_features_are_the_wild_type_in_all_three_modes(self) -> None:
        """The repair must move the LABEL, not the point being anchored."""
        (X_zero, _yz), _ = _run()
        (X_meas, _ym), _ = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232)
        np.testing.assert_array_equal(X_zero, X_meas)

    def test_it_refuses_to_run_without_a_label(self) -> None:
        """No defensible default exists: the defect being repaired IS a fabricated label,
        so a silent fallback to 0.0 would be the same defect wearing the fix's name."""
        for bad in (None, float("nan"), float("inf")):
            with self.subTest(label=bad):
                with self.assertRaises(ValueError) as ctx:
                    _run(wt_anchor_mode="measured", wt_anchor_label=bad)
                self.assertIn("wt_anchor_label", str(ctx.exception))

    def test_the_variance_defaults_to_a_real_measurements_and_is_overridable(self) -> None:
        _xy, batch = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232)
        self.assertEqual(batch.wt_anchor_state["variance"],
                         math_worker.WT_ANCHOR_MEASURED_VARIANCE)
        _xy2, batch2 = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232,
                            wt_anchor_variance=1.76e-3)
        self.assertAlmostEqual(batch2.wt_anchor_state["variance"], 1.76e-3)

    def test_an_override_cannot_reintroduce_the_1e_4_overconfidence(self) -> None:
        _xy, batch = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232,
                          wt_anchor_variance=1e-12)
        self.assertGreaterEqual(batch.wt_anchor_state["variance"], 1e-4)

    def test_the_anchor_is_in_residual_space_like_every_measured_row(self) -> None:
        """The GP is fitted on `y - m(x)`. The `zero` path hardcodes a prior of 0.0 for
        the anchor; the repair EVALUATES the prior at the wild type, so the anchor sits
        in the same space as the data beside it."""
        prior = mock.Mock(return_value=0.25)
        (_X, y), batch = _run(wt_anchor_mode="measured", wt_anchor_label=4.0,
                              prior_mean_fn=prior)
        self.assertAlmostEqual(batch.wt_anchor_state["prior_at_wt"], 0.25)
        self.assertIn(4.0 - 0.25, [round(v, 10) for v in y.tolist()])
        self.assertEqual(batch.wt_anchor_state["label"], 4.0,
                         "the recorded label must be the MEASUREMENT, not the residual")


class TestTheThreeModesAreMutuallyExclusive(unittest.TestCase):
    def test_ablating_and_requesting_measured_is_refused(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            _run(ablate=("wt_anchor",), wt_anchor_mode="measured",
                 wt_anchor_label=4.0232)
        self.assertIn("ablated", str(ctx.exception).lower())

    def test_an_unknown_mode_raises_rather_than_falling_back(self) -> None:
        with self.assertRaises(ValueError):
            _run(wt_anchor_mode="wild_guess")

    def test_the_ablated_arm_reports_itself_as_ablated(self) -> None:
        (_X, y), batch = _run(ablate=("wt_anchor",))
        self.assertEqual(batch.wt_anchor_state["mode"], "ablated")
        self.assertEqual(batch.wt_anchor_state["rows"], 0)
        self.assertEqual(len(y), 3)

    def test_all_three_modes_give_three_distinct_training_sets(self) -> None:
        """If any two coincided, a three-arm campaign would be a two-arm campaign."""
        (_, y_zero), _ = _run()
        (_, y_abl), _ = _run(ablate=("wt_anchor",))
        (_, y_meas), _ = _run(wt_anchor_mode="measured", wt_anchor_label=4.0232)
        sets = [tuple(sorted(v.tolist())) for v in (y_zero, y_abl, y_meas)]
        self.assertEqual(len(set(sets)), 3, f"two modes produced the same fit: {sets}")


class TestTheVarianceVectorStaysAligned(unittest.TestCase):
    """A misalignment used to become `None` rather than raise: the array was built only
    `if len(sample_variances_list) == len(X_train)`, and `DeepKernelLearner.fit` ALSO
    drops a mismatched vector silently - so the round quietly became homoscedastic while
    `y_train`, misaligned by the same event, stayed misaligned."""

    def _captured(self, **kwargs):
        captured = {}
        real_dkl = math_worker.DeepKernelLearner

        class Spy(real_dkl):  # type: ignore[misc,valid-type]
            def fit(self, X, y, sample_variance=None, **kw):
                captured["n_rows"] = len(y)
                captured["n_var"] = (None if sample_variance is None
                                     else len(sample_variance))
                captured["variances"] = (None if sample_variance is None
                                         else list(sample_variance))
                return super().fit(X, y, sample_variance=sample_variance, **kw)

        with mock.patch.object(math_worker, "DeepKernelLearner", Spy), \
             mock.patch.object(math_worker, "TWO_STAGE_GP_ENABLED", True):
            np.random.seed(5)
            math_worker.select_next_micro_batch(
                strategy=_strategy(), wt_sequence=WT, evidence=EVIDENCE, batch_size=1,
                candidate_pool=POOL, random_seed=3, **kwargs,
            )
        self.assertIn("n_rows", captured, "the DKL surrogate was never fitted, so the "
                                         "alignment was not exercised")
        return captured

    def test_aligned_in_all_three_modes(self) -> None:
        for name, kwargs in (
            ("zero", {}),
            ("ablated", {"ablate": ("wt_anchor",)}),
            ("measured", {"wt_anchor_mode": "measured", "wt_anchor_label": 4.0232}),
        ):
            with self.subTest(mode=name):
                c = self._captured(**kwargs)
                self.assertIsNotNone(
                    c["n_var"],
                    "sample_variance arrived as None, so heteroscedastic noise was "
                    "silently disabled for the whole arm")
                self.assertEqual(c["n_var"], c["n_rows"])

    def test_the_measured_anchor_carries_a_real_measurements_variance(self) -> None:
        c = self._captured(wt_anchor_mode="measured", wt_anchor_label=4.0232)
        self.assertIn(math_worker.WT_ANCHOR_MEASURED_VARIANCE, c["variances"])
        self.assertNotIn(1e-4, c["variances"],
                         "the 1e-4 over-confidence survived the repair")

    def test_a_desynchronised_append_raises_instead_of_degrading(self) -> None:
        """The guard itself. Forcing one variance entry to go missing must VOID the
        round, not quietly turn the vector into None - because `y_train`, misaligned by
        the same event, would stay misaligned and nothing in the artifacts would say so.
        """
        real_featurize = math_worker.featurize_variant_on_positions
        calls = {"n": 0}

        def flaky_prior(variant):
            # Raises on the second evidence row, AFTER its feature vector is appended.
            calls["n"] += 1
            if calls["n"] == 2:
                raise RuntimeError("synthetic prior failure between the two appends")
            return 0.0

        with mock.patch.object(math_worker, "featurize_variant_on_positions",
                               side_effect=real_featurize):
            np.random.seed(5)
            # The evidence loop swallows the exception, so X gains a row that y and the
            # variance list never get. The guard must catch the resulting mismatch.
            with self.assertRaises(RuntimeError) as ctx:
                math_worker.select_next_micro_batch(
                    strategy=_strategy(), wt_sequence=WT, evidence=EVIDENCE,
                    batch_size=1, candidate_pool=POOL, random_seed=3,
                    prior_mean_fn=flaky_prior,
                )
        self.assertIn("misalignment", str(ctx.exception))


class TestTheHarnessRecordsWhatTheEngineDid(unittest.TestCase):
    """Rule 2.12: a manifest's component state is only as honest as what it reads. Two
    manifests wrongly claimed dpp and dkl were enabled because the recorder read the CLI
    flag rather than the module-level global."""

    def test_the_state_comes_back_on_the_batch(self) -> None:
        _xy, batch = _run()
        self.assertIsNotNone(batch.wt_anchor_state)
        self.assertIn("mode", batch.wt_anchor_state)
        self.assertIn("reason", batch.wt_anchor_state)

    def test_the_verification_declares_itself_vacuous_when_no_arm_called_the_engine(
            self) -> None:
        """A check that cannot fail must say so. Rule 2.10: the pre-flight threshold
        check's first version passed the exact defect it was written for."""
        import benchmark_panoptes_isolated as bench

        empty = bench.RunOutcome(method="ridge", seed=0, status="completed")
        block = bench._wt_anchor_verification("zero", [empty])
        self.assertTrue(block["vacuous"])
        self.assertFalse(block["agrees_with_configuration"])
        self.assertIn("VACUOUS", block["note"])

    def test_the_verification_is_not_vacuous_when_the_engine_ran(self) -> None:
        import benchmark_panoptes_isolated as bench

        ran = bench.RunOutcome(
            method="panoptes_default", seed=0, status="completed",
            wt_anchor_observed=[{"round": 1, "mode": "measured", "label": 4.0232}],
        )
        block = bench._wt_anchor_verification("measured", [ran])
        self.assertFalse(block["vacuous"])
        self.assertTrue(block["agrees_with_configuration"])
        self.assertEqual(block["modes_observed"], ["measured"])

    def test_a_disagreement_is_reported_rather_than_smoothed_over(self) -> None:
        import benchmark_panoptes_isolated as bench

        lying = bench.RunOutcome(
            method="panoptes_default", seed=0, status="completed",
            wt_anchor_observed=[{"round": 1, "mode": "zero", "label": 0.0}],
        )
        block = bench._wt_anchor_verification("measured", [lying])
        self.assertFalse(block["agrees_with_configuration"])
        self.assertEqual(block["modes_observed"], ["zero"])

    def test_the_selector_refuses_measured_on_a_library_with_no_recorded_label(
            self) -> None:
        """Never invent one."""
        import benchmark_panoptes_isolated as bench
        from dataclasses import replace

        library = replace(bench.LIBRARIES["trpb"], parent_measured_fitness=None,
                          parent_fitness_source=None)
        with self.assertRaises(bench.BenchmarkError) as ctx:
            bench.PanoptesSelector(
                0, "M" + "A" * 300, object(), False, None, library,
                wt_anchor_mode="measured")
        self.assertIn("parent_measured_fitness", str(ctx.exception))

    def test_the_cli_rejects_ablation_plus_measured(self) -> None:
        import benchmark_panoptes_isolated as bench

        with self.assertRaises(bench.BenchmarkError):
            bench.PanoptesSelector(
                0, "M" + "A" * 300, object(), False, None, bench.LIBRARIES["trpb"],
                ablate_components={"wt_anchor"}, wt_anchor_mode="measured")


class TestEveryLibraryRecordsItsParentLabel(unittest.TestCase):
    def test_all_nine_have_one_with_a_source(self) -> None:
        """`scripts/analysis/wt_label_availability.py` resolved all nine from three
        independent sources and cross-checked each against the bracket the prepared
        table's own `DMS_score_bin` column implies. A library carrying a value with no
        recorded provenance is worse than one carrying none, because it will be trusted.
        """
        import benchmark_panoptes_isolated as bench

        for key, library in bench.LIBRARIES.items():
            with self.subTest(library=key):
                self.assertIsNotNone(library.parent_measured_fitness,
                                     f"{key} records no parent label")
                self.assertTrue(library.parent_fitness_source,
                                f"{key} records a parent label with no provenance")

    def test_the_recorded_labels_match_the_landscape_loaders_view(self) -> None:
        """The two landscapes whose prepared table KEEPS a parent row must agree with it,
        so a typo in the registry cannot pass."""
        import benchmark_panoptes_isolated as bench

        self.assertEqual(bench.LIBRARIES["gb1"].parent_measured_fitness, 1.0)
        self.assertEqual(bench.LIBRARIES["her2"].parent_measured_fitness, 0.1)

    def test_creilov_is_the_largest_violation_and_trpb_is_not_a_null(self) -> None:
        """The reframe this work turns on. `|0 - label mean| / sd` says TrpB is a
        negative control at 0.35 sd. `|wt_label - 0| / sd` - how wrong the anchor is
        about the point it claims to pin - says 6.81, and the parent is itself a top-1%
        hit there. The second quantity is the one the repair changes.
        """
        import benchmark_panoptes_isolated as bench
        from veritas.calibration.landscapes import load_landscape

        errors = {}
        for key in ("creilov", "trpb", "gb1", "eqfp611_red"):
            ls = load_landscape(key)
            label = bench.LIBRARIES[key].parent_measured_fitness
            errors[key] = abs(float(label)) / float(ls.fitness.std(ddof=0))
        self.assertGreater(errors["creilov"], 15.0)
        self.assertGreater(errors["trpb"], 5.0,
                           "TrpB is being treated as a negative control for this repair; "
                           "it is not")
        self.assertLess(errors["eqfp611_red"], 1.0,
                        "eqFP611-red is the genuine negative control and has stopped "
                        "being one")


if __name__ == "__main__":
    unittest.main()
