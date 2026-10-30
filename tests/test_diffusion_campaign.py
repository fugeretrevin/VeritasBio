"""
tests/test_diffusion_campaign.py — Verification of Autonomous De Novo Diffusion CLI Pipeline.

Test coverage:
  1. Full end-to-end campaign execution with mock diffusion.
  2. Verification of physical lab file emissions (Twist CSV, Opentrons OT-2, FDA certificate, PDBs).
  3. Fail-closed budget rejection when funds are insufficient.
  4. Cryptographic FDA 21 CFR Part 11 audit integrity check on campaign ledger.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ledger import database as ledger_db
from scripts.run_diffusion_campaign import run_campaign


class TestDiffusionCampaign(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.out_path = Path(self.temp_dir.name)
        self.db_path = self.out_path / "test_campaign_ledger.db"

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_cli_campaign_execution(self):
        """Verifies full execution of de novo diffusion campaign and file artifact emissions."""
        result = run_campaign(
            target_pdb="2YBR",
            num_designs=3,
            binder_length=65,
            budget_usd=0.05,
            target_model="mock_diffusion",
            output_dir=self.out_path,
            db_path=self.db_path,
        )

        self.assertTrue(result.approved)
        self.assertEqual(result.status, "COMPLETED")
        self.assertEqual(len(result.designs), 3)

        # Check emitted files
        twist_file = self.out_path / "twist_de_novo_orders.csv"
        self.assertTrue(twist_file.is_file(), "Missing twist_de_novo_orders.csv")
        self.assertIn("cDNA_Sequence", twist_file.read_text(encoding="utf-8"))

        robot_file = self.out_path / "opentrons_de_novo_protocol.py"
        self.assertTrue(robot_file.is_file(), "Missing opentrons_de_novo_protocol.py")
        self.assertIn("protocol_api", robot_file.read_text(encoding="utf-8"))

        cert_file = self.out_path / "fda_part11_de_novo_certificate.json"
        self.assertTrue(cert_file.is_file(), "Missing fda_part11_de_novo_certificate.json")
        self.assertIn("ledger_hash_sha256", cert_file.read_text(encoding="utf-8"))

        pdb_file = self.out_path / "pdb" / "VB-DIFF-001.pdb"
        self.assertTrue(pdb_file.is_file(), "Missing VB-DIFF-001.pdb")
        self.assertIn("ATOM", pdb_file.read_text(encoding="utf-8"))

        # Check cryptographic audit chain in campaign ledger
        is_valid, count, msg = ledger_db.verify_audit_integrity(db_path=self.db_path)
        self.assertTrue(is_valid, f"Audit chain verification failed: {msg}")
        self.assertGreater(count, 0)

    def test_insufficient_budget_rejection(self):
        """Verifies fail-closed halt when budget is below required trajectory cost."""
        result = run_campaign(
            target_pdb="2YBR",
            num_designs=5,
            binder_length=65,
            budget_usd=0.0001,
            # Implemented engine, so this covers budget gating rather than
            # engine availability. mock_diffusion costs $0.0005 > $0.0001 budget.
            target_model="mock_diffusion",
            output_dir=self.out_path,
            db_path=self.db_path,
        )

        self.assertFalse(result.approved)
        self.assertEqual(result.status, "REJECTED")
        self.assertEqual(len(result.designs), 0)


if __name__ == "__main__":
    unittest.main()
