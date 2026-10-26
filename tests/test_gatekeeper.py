"""
tests/test_gatekeeper.py — Unit test suite for controller.gatekeeper.

Built using standard library `unittest` (zero external test runner dependencies).
Tests cover:
    1. Biophysical validation (SingleMutation, multi-point Variant, canonical AAs, WT match).
    2. Fail-closed behavior on malformed inputs, unknown models, and out-of-bounds positions.
    3. Two-state cost calculations (cold start vs. warm keep-alive amortization).
    4. ROI ($EI/$) evaluations.
    5. Integration with ledger.database.authorize_spend() (approvals, over-budget rejections, audit logs).
"""

import math
import os
import tempfile
import unittest

from controller.gatekeeper import (
    BatchEvaluationRequest,
    Gatekeeper,
    GatekeeperConfig,
    SingleMutation,
    Variant,
)
from ledger import database as ledger_db

# GFP wild-type sequence snippet (chromophore region around residues 64-68 is LSYGV)
GFP_MINI_WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHM"


class TestGatekeeper(unittest.TestCase):

    # -------------------------------------------------------------------------
    # 1. Biophysical & Pydantic Schema Tests
    # -------------------------------------------------------------------------
    def test_single_mutation_parsing(self):
        m = SingleMutation.from_string("S65T")
        self.assertEqual(m.wt_aa, "S")
        self.assertEqual(m.position, 65)
        self.assertEqual(m.mut_aa, "T")
        self.assertEqual(m.to_string(), "S65T")

    def test_single_mutation_invalid_syntax(self):
        with self.assertRaises(ValueError):
            SingleMutation.from_string("invalid_mutation")
        with self.assertRaises(ValueError):
            SingleMutation.from_string("65T")
        with self.assertRaises(ValueError):
            SingleMutation.from_string("S65")

    def test_single_mutation_non_canonical_aa(self):
        with self.assertRaises(ValueError):
            SingleMutation.from_string("S65Z")  # Z is not canonical
        with self.assertRaises(ValueError):
            SingleMutation.from_string("B65T")  # B is not canonical

    def test_single_mutation_same_wt_and_mut(self):
        with self.assertRaises(ValueError):
            SingleMutation.from_string("S65S")  # No substitution

    def test_variant_multi_point_parsing(self):
        v = Variant(notation="S65T:Y66H:V68A")
        self.assertEqual(len(v.mutations), 3)
        self.assertEqual([m.to_string() for m in v.mutations], ["S65T", "Y66H", "V68A"])

    def test_variant_colliding_positions(self):
        with self.assertRaises(ValueError):
            Variant(notation="S65T:S65A")  # Double mutation at same site

    def test_variant_validation_against_wt(self):
        v = Variant(notation="S65T:Y66H")
        v.validate_against_wt(GFP_MINI_WT)  # Should pass without error

        # Apply mutation and check full mutated string
        mutated = v.apply_to(GFP_MINI_WT)
        self.assertEqual(len(mutated), len(GFP_MINI_WT))
        self.assertEqual(mutated[64], "T")  # 1-indexed 65
        self.assertEqual(mutated[65], "H")  # 1-indexed 66

    def test_variant_mismatch_against_wt(self):
        # If mutation claims position 65 is 'A', but WT has 'S' -> must raise!
        v = Variant(notation="A65T")
        with self.assertRaises(ValueError) as cm:
            v.validate_against_wt(GFP_MINI_WT)
        self.assertIn("Biophysical mismatch", str(cm.exception))

    def test_variant_position_exceeds_wt(self):
        v = Variant(notation="S999T")
        with self.assertRaises(ValueError) as cm:
            v.validate_against_wt(GFP_MINI_WT)
        self.assertIn("exceeds WT sequence length", str(cm.exception))

    # -------------------------------------------------------------------------
    # 2. Cost Engine & Two-State Amortization Tests
    # -------------------------------------------------------------------------
    def test_two_state_cost_model(self):
        gk = Gatekeeper(
            config=GatekeeperConfig(
                a10g_usd_per_sec=0.000306,
                gpu_sec_per_variant=0.02,
                cold_start_sec=20.0,
                idle_timeout_sec=120.0,
            )
        )

        # 1. Cold state: 50 variants
        # Total GPU sec = 50 * 0.02 (1.0s) + 20.0s = 21.0s
        # Cost = 21.0 * 0.000306 = $0.006426
        cost_cold, is_warm_cold = gk.calculate_cost(
            variant_count=50, target_model="esm2_650m", force_warm=False
        )
        self.assertFalse(is_warm_cold)
        self.assertTrue(math.isclose(cost_cold, 0.006426, rel_tol=1e-4))

        # 2. Warm state: 50 variants
        # Total GPU sec = 50 * 0.02 (1.0s) + 0.0s boot = 1.0s
        # Cost = 1.0 * 0.000306 = $0.000306
        cost_warm, is_warm_warm = gk.calculate_cost(
            variant_count=50, target_model="esm2_650m", force_warm=True
        )
        self.assertTrue(is_warm_warm)
        self.assertTrue(math.isclose(cost_warm, 0.000306, rel_tol=1e-4))

        # Warm execution should be ~21x cheaper!
        self.assertGreater(cost_cold, cost_warm * 20)

    def test_mock_cost(self):
        gk = Gatekeeper()
        cost_mock, is_warm = gk.calculate_cost(variant_count=50, target_model="mock_esm")
        self.assertEqual(cost_mock, gk.config.mock_nominal_cost)
        self.assertGreater(cost_mock, 0.0)

    # -------------------------------------------------------------------------
    # 3. ROI Evaluation Tests
    # -------------------------------------------------------------------------
    def test_roi_calculation(self):
        gk = Gatekeeper()
        roi = gk.calculate_roi(expected_improvements=[0.5, 0.3, 0.2], cost_usd=0.005)
        self.assertEqual(roi, 200.0)
        self.assertIsNone(gk.calculate_roi(None, 0.005))

    # -------------------------------------------------------------------------
    # 4. Fail-Closed Integration with Ledger
    # -------------------------------------------------------------------------
    def test_gatekeeper_authorize_success(self):
        tmp_dir = tempfile.mkdtemp()
        db_path = os.path.join(tmp_dir, "test_ledger.db")
        ledger_db.init_db(db_path)
        ledger_db.create_budget(total_allocation=1.0, db_path=db_path)

        gk = Gatekeeper(db_path=db_path)

        request = BatchEvaluationRequest(
            job_id="job-gatekeeper-01",
            target_model="esm2_650m",
            wildtype_sequence=GFP_MINI_WT,
            variants=[
                Variant(notation="S65T"),
                Variant(notation="Y66H"),
            ],
            expected_improvements=[0.45, 0.35],
        )

        decision = gk.evaluate_and_authorize(request)

        self.assertTrue(decision.approved)
        self.assertEqual(decision.job_id, "job-gatekeeper-01")
        self.assertGreater(decision.estimated_cost, 0.0)
        self.assertIsNotNone(decision.roi_ei_per_dollar)
        self.assertEqual(len(decision.mutated_sequences), 2)

        # Verify ledger was updated
        job_in_ledger = ledger_db.get_job("job-gatekeeper-01", db_path=db_path)
        self.assertIsNotNone(job_in_ledger)
        self.assertEqual(job_in_ledger["status"], "APPROVED")

    def test_gatekeeper_overbudget_rejected(self):
        tmp_dir = tempfile.mkdtemp()
        db_path = os.path.join(tmp_dir, "test_ledger.db")
        ledger_db.init_db(db_path)
        # Tiny budget: $0.000005 (less than cold start cost)
        ledger_db.create_budget(total_allocation=0.000005, db_path=db_path)

        gk = Gatekeeper(db_path=db_path)

        request = BatchEvaluationRequest(
            job_id="job-overbudget",
            target_model="esm2_650m",
            wildtype_sequence=GFP_MINI_WT,
            variants=[Variant(notation="S65T")],
        )

        decision = gk.evaluate_and_authorize(request, force_warm=False)

        self.assertFalse(decision.approved)
        self.assertIn("rejected", decision.rejection_reason.lower())

        # Rejection logged in ledger for audit trail
        job_in_ledger = ledger_db.get_job("job-overbudget", db_path=db_path)
        self.assertIsNotNone(job_in_ledger)
        self.assertEqual(job_in_ledger["status"], "REJECTED")

    def test_gatekeeper_fail_closed_on_bad_payload(self):
        gk = Gatekeeper()
        bad_payload = {"random_key": "junk"}
        decision = gk.evaluate_and_authorize(bad_payload)

        self.assertFalse(decision.approved)
        self.assertIn("validation error", decision.rejection_reason.lower())

    def test_gatekeeper_fail_closed_on_unknown_model(self):
        gk = Gatekeeper()
        bad_payload = {
            "job_id": "bad-model-job",
            "target_model": "gpt4_physics_engine",  # disallowed model
            "wildtype_sequence": GFP_MINI_WT,
            "variants": [{"notation": "S65T"}],
        }
        decision = gk.evaluate_and_authorize(bad_payload)
        self.assertFalse(decision.approved)
        self.assertIn("allowed models", decision.rejection_reason.lower())


if __name__ == "__main__":
    unittest.main()
