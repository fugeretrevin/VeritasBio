"""
Regression tests for the prior fitness-scale calibration, the de-double-counted
desolvation penalty, and the opt-in auxiliary objective weighting.

These lock in three properties that were defects when measured on the TrpB
landscape:

  * the prior is consumed as a GP prior MEAN over measured fitness, so it has to
    live on the fitness scale. It previously spanned -15.3 to +0.8 against an
    observable range of -0.16 to 1.00.
  * the buried-core and cofactor-shell terms describe one desolvation event and
    both fired on it, so V227K paid -2.0 and -2.5 for a single physical cause.
  * the synthesizability and developability proxies are anti-correlated with
    measured fitness, and multiplied acquisition value by as little as 0.4.
"""
from __future__ import annotations

import math
import unittest

from controller.gatekeeper import Variant


class TestPriorFitnessScaleCalibration(unittest.TestCase):
    """
    Mechanism tests for the calibration. The SHIPPED DEFAULT is off, because a paired
    48-seed benchmark measured it worse (33.31 vs 45.60 top-1% hits, t = 12.50, worse on
    46 of 48 seeds), so each test here enables it explicitly. The default itself is
    asserted separately in TestCalibrationDefaultIsOff.
    """

    def setUp(self) -> None:
        import workers.math_worker as mw
        self._saved = mw.PRIOR_CALIBRATION_ENABLED
        mw.PRIOR_CALIBRATION_ENABLED = True

    def tearDown(self) -> None:
        import workers.math_worker as mw
        mw.PRIOR_CALIBRATION_ENABLED = self._saved
    def test_calibration_is_bounded_by_the_configured_spread(self) -> None:
        from workers.math_worker import PRIOR_FITNESS_SPREAD, calibrate_prior_to_fitness_scale

        # The raw heuristic reaches about -15 on real variants; measured fitness
        # spans roughly -0.16 to 1.0, so the prior must not exceed its spread.
        for raw in (-1e6, -100.0, -15.3, -4.5, 0.0, 0.8, 100.0, 1e6):
            out = calibrate_prior_to_fitness_scale(raw)
            self.assertLessEqual(abs(out), PRIOR_FITNESS_SPREAD + 1e-9,
                                 msg=f"raw={raw} escaped the spread bound")

    def test_calibration_preserves_ranking(self) -> None:
        """Monotone, so it fixes the magnitude defect without altering any ranking."""
        from workers.math_worker import calibrate_prior_to_fitness_scale

        raws = [-15.3, -8.0, -4.5, -2.75, -1.0, -0.3, 0.0, 0.15, 0.3, 0.8]
        outs = [calibrate_prior_to_fitness_scale(r) for r in raws]
        self.assertEqual(outs, sorted(outs),
                         "calibration must be monotone in the raw score")
        # Strictly increasing, so distinct raw scores stay distinguishable.
        for a, b in zip(outs, outs[1:]):
            self.assertLess(a, b)

    def test_zero_spread_ablates_the_prior(self) -> None:
        from workers.math_worker import calibrate_prior_to_fitness_scale

        # The hand-built prior is fitted to nothing; ablating it is a live question,
        # so the switch has to work.
        for raw in (-15.3, -4.5, 0.0, 0.8):
            self.assertEqual(calibrate_prior_to_fitness_scale(raw, spread=0.0), 0.0)

    def test_non_positive_softening_is_rejected(self) -> None:
        from workers.math_worker import calibrate_prior_to_fitness_scale

        with self.assertRaises(ValueError):
            calibrate_prior_to_fitness_scale(-1.0, softening=0.0)

    def test_calibration_can_be_bypassed_for_ablation_arms(self) -> None:
        import workers.math_worker as mw

        # setUp enabled calibration, so it is bounded here.
        self.assertLessEqual(abs(mw.calibrate_prior_to_fitness_scale(-15.3)),
                             mw.PRIOR_FITNESS_SPREAD + 1e-9)
        # Disabling it passes the raw score through. Raising the spread cannot achieve
        # this, because tanh saturates, which is why the bypass is a separate switch.
        mw.PRIOR_CALIBRATION_ENABLED = False
        self.assertAlmostEqual(mw.calibrate_prior_to_fitness_scale(-15.3), -15.3)

    def test_biological_prior_output_stays_on_the_fitness_scale(self) -> None:
        from workers.math_worker import PRIOR_FITNESS_SPREAD, compute_biological_prior

        # A multi-mutant hitting catalytic invariants used to return about -13.5.
        variant = Variant(notation="V183K:F184R:V227D:S228E")
        prior = compute_biological_prior(variant, catalytic_invariants={183, 227})
        self.assertLessEqual(abs(prior), PRIOR_FITNESS_SPREAD + 1e-9)


class TestDesolvationIsBilledOnce(unittest.TestCase):
    def test_buried_and_cofactor_terms_do_not_both_charge_full_price(self) -> None:
        from pathlib import Path

        from workers.spatial_filter import compute_structural_surrogate_score

        pdb = Path(__file__).resolve().parent.parent / "workers" / "data" / "8EH0.pdb"
        if not pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        # V227 is buried AND within 4.5 A of the PLP cofactor, so V227K used to take
        # -2.0 from the buried-core term plus -2.5 from the cofactor-shell term for
        # one desolvation event, totalling the -4.5 that ranked the global optimum as
        # catastrophic.
        score = compute_structural_surrogate_score("V227K", pdb_path=pdb)
        self.assertGreater(score, -4.5,
                           "V227K is still paying the full double penalty")

    def test_formal_charge_costs_more_than_a_neutral_amide(self) -> None:
        from pathlib import Path

        from workers.spatial_filter import compute_structural_surrogate_score

        pdb = Path(__file__).resolve().parent.parent / "workers" / "data" / "8EH0.pdb"
        if not pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        # The old "KRDENQ" class charged K/R/D/E and N/Q identically.
        charged = compute_structural_surrogate_score("V227K", pdb_path=pdb)
        amide = compute_structural_surrogate_score("V227N", pdb_path=pdb)
        self.assertLess(charged, amide,
                        "burying a formal charge should cost more than an amide")


class TestAuxiliaryObjectiveWeightingIsOptIn(unittest.TestCase):
    def test_default_does_not_apply_the_proxy_multipliers(self) -> None:
        import inspect

        from workers.math_worker import select_next_micro_batch

        sig = inspect.signature(select_next_micro_batch)
        param = sig.parameters.get("apply_auxiliary_objective_weights")
        self.assertIsNotNone(param, "the opt-in parameter is missing")
        # Measured on TrpB: synthesizability rho = -0.104 against fitness, and the
        # combined multiplier averages 0.722 on the five best variants against 0.782
        # for random ones. It must not be applied unless a caller asks.
        self.assertFalse(param.default)

    def test_proxy_scores_are_still_computed_for_reporting(self) -> None:
        from workers.math_worker import (
            analyze_biologics_developability,
            evaluate_synthesizability_score,
        )

        # Turning the weighting off must not remove the measurements themselves.
        wt = "M" + "A" * 240
        variant = Variant(notation="A100V")
        synth = evaluate_synthesizability_score(variant, wt, active_site_positions={100})
        dev = analyze_biologics_developability(variant.apply_to(wt))
        self.assertIsInstance(float(synth), float)
        self.assertIn("developability_score", dev)


class TestCalibrationDefaultIsOff(unittest.TestCase):
    def test_shipped_default_leaves_the_prior_uncalibrated(self) -> None:
        import workers.math_worker as mw

        # Measured, paired, 48 identical seeds:
        #   calibrated (spread 0.25)  0.7149 best / 33.31 top-1%
        #   raw magnitudes           0.7462 best / 45.60 top-1%
        #   t = 4.77 and 12.50; the calibrated arm lost on 46 of 48 seeds.
        # The term behaves as a weighted ranking feature whose weight is its
        # magnitude, not as a prior mean, so bounding it removed most of its
        # contribution. Flip this only with a paired benchmark run to back it.
        self.assertFalse(
            mw.PRIOR_CALIBRATION_ENABLED,
            "prior calibration is measured worse; defaulting it on needs new evidence",
        )


if __name__ == "__main__":
    unittest.main()
