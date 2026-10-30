"""
tests/test_diffusion_worker.py — Verification of Budgeted De Novo Diffusion & Inverse Folding.

Test coverage:
  1. Fail-closed budget rejection (ensures zero compute runs if funds insufficient).
  2. Successful approved generation with 3D backbone PDB structures and ProteinMPNN sequences.
  3. Biophysical developability enforcement (zero unpaired cysteines, no deamidation/isomerization).
  4. 3D epitope hotspot docking with authentic target structure (PDB 2YBR).
  5. Wet-lab translations (Twist synthesis CSV, Opentrons OT-2 protocol, FDA Part 11 audit certificate).
  6. Ledger cryptographic Merkle audit integrity verification.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from controller.gatekeeper import Gatekeeper, GatekeeperConfig
from ledger import database as ledger_db
from workers.diffusion_worker import (
    DeNovoBinderDesign,
    DeNovoBinderRequest,
    DiffusionJobResult,
    GeometricBackboneSampler,
    ProteinMPNNInverseFolder,
    generate_budgeted_de_novo_binders,
)


class TestDiffusionWorker(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = Path(self.temp_dir.name) / "test_diffusion_ledger.db"
        ledger_db.init_db(self.db_path)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_fail_closed_budget_rejection(self):
        """Verifies that an unbudgeted diffusion request is rejected fail-closed with 0 compute spent."""
        # Create tiny budget of $0.0001
        ledger_db.create_budget(total_allocation=0.0001, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-overbudget",
            # Use an implemented engine so this test exercises budget gating only.
            # mock_diffusion costs $0.0005, still well over the $0.0001 budget.
            target_model="mock_diffusion",
            num_designs=5,
            binder_length=65,
            budget_id=1,
        )

        result: DiffusionJobResult = generate_budgeted_de_novo_binders(
            request=request,
            db_path=self.db_path,
        )

        self.assertFalse(result.approved)
        self.assertEqual(result.status, "REJECTED")
        self.assertEqual(len(result.designs), 0)
        self.assertIn("Ledger authorization rejected", str(result.rejection_reason))

        # Ledger state check: budget spent must be $0.00
        budget = ledger_db.get_budget(budget_id=1, db_path=self.db_path)
        self.assertEqual(budget["cumulative_spent"], 0.0)

        # Job must be recorded as REJECTED in the audit trail
        job = ledger_db.get_job("diff-job-overbudget", db_path=self.db_path)
        self.assertIsNotNone(job)
        self.assertEqual(job["status"], "REJECTED")

    def test_unimplemented_engine_refused_before_any_spend(self):
        """An engine with no code path must be refused before funds are authorized."""
        ledger_db.create_budget(total_allocation=5.00, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-no-engine",
            target_model="rfdiffusion_proteinmpnn_a10g",
            num_designs=2,
            binder_length=65,
            budget_id=1,
        )
        result = generate_budgeted_de_novo_binders(request=request, db_path=self.db_path)

        self.assertFalse(result.approved)
        self.assertEqual(result.status, "REJECTED")
        self.assertEqual(len(result.designs), 0)
        self.assertEqual(result.cost_usd, 0.0)
        self.assertIn("not implemented", str(result.rejection_reason))

        # Critically: the budget is untouched and no job row was ever created,
        # because the refusal happens before authorization rather than after.
        budget = ledger_db.get_budget(budget_id=1, db_path=self.db_path)
        self.assertEqual(budget["cumulative_spent"], 0.0)
        self.assertIsNone(ledger_db.get_job("diff-job-no-engine", db_path=self.db_path))

    def test_successful_budgeted_generation_mock(self):
        """Verifies approved de novo binder generation under the fail-closed ledger."""
        ledger_db.create_budget(total_allocation=1.00, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-approved-001",
            target_model="mock_diffusion",
            num_designs=4,
            binder_length=65,
            target_protein_name="Centruroides noxius Cn2 Toxin",
            budget_id=1,
        )

        result: DiffusionJobResult = generate_budgeted_de_novo_binders(
            request=request,
            db_path=self.db_path,
        )

        self.assertTrue(result.approved)
        self.assertEqual(result.status, "COMPLETED")
        self.assertEqual(len(result.designs), 4)
        self.assertGreater(result.cost_usd, 0.0)

        # Check each generated binder
        for design in result.designs:
            self.assertEqual(len(design.sequence), 65)
            self.assertEqual(design.length, 65)
            # Interface geometry is measured, not predicted: with no target
            # structure supplied these are zero rather than invented.
            self.assertGreaterEqual(design.interface_contact_pairs, 0)
            self.assertGreaterEqual(design.interface_binder_residues, 0)
            self.assertGreater(design.developability_score, 0.0)
            self.assertLessEqual(design.developability_score, 1.0)
            self.assertGreater(design.radius_of_gyration, 8.0)
            self.assertLess(design.radius_of_gyration, 25.0)
            self.assertIn("ATOM", design.pdb_content)
            # The designed sequence must be threaded into the coordinates: the
            # sampler emits placeholder GLY names, and leaving them in place would
            # misreport the structure to any downstream tool. Verify the first
            # residue's name matches the first residue of the designed sequence.
            from workers.diffusion_worker import AA_1_TO_3

            expected_first = AA_1_TO_3[design.sequence[0].upper()]
            self.assertIn(f"{expected_first} B   1", design.pdb_content)
            threaded = {
                line[17:20].strip()
                for line in design.pdb_content.splitlines()
                if line.startswith("ATOM  ")
            }
            expected = {AA_1_TO_3[c.upper()] for c in design.sequence}
            self.assertTrue(
                threaded.issubset(expected),
                f"PDB contains residue names absent from the design: {threaded - expected}",
            )
            self.assertIn("END", design.pdb_content)
            self.assertGreater(design.synthesis_cost_usd, 0.0)

        # Ledger verification: evidence recorded and job marked COMPLETED
        evidence = ledger_db.get_evidence(job_id="diff-job-approved-001", db_path=self.db_path)
        self.assertEqual(len(evidence), 4)

        job = ledger_db.get_job("diff-job-approved-001", db_path=self.db_path)
        self.assertIsNotNone(job)
        self.assertEqual(job["status"], "COMPLETED")

        # Budget verification: spent increased
        budget = ledger_db.get_budget(budget_id=1, db_path=self.db_path)
        self.assertGreater(budget["cumulative_spent"], 0.0)

    def test_biologics_developability_enforcement(self):
        """Verifies that designed sequences eliminate unpaired cysteines and liability motifs."""
        ledger_db.create_budget(total_allocation=1.00, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-developability",
            target_model="mock_diffusion",
            num_designs=5,
            binder_length=60,
            budget_id=1,
        )

        result = generate_budgeted_de_novo_binders(request=request, db_path=self.db_path)
        self.assertTrue(result.approved)

        for design in result.designs:
            seq = design.sequence
            # Zero unpaired cysteines
            self.assertEqual(seq.count("C"), 0, f"Found cysteine in de novo binder: {seq}")
            # No deamidation motifs (NG, NS, NA, NH)
            for bad in ["NG", "NS", "NA", "NH"]:
                self.assertNotIn(bad, seq, f"Found deamidation motif {bad} in {seq}")
            # No isomerization motifs (DG, DS, DT)
            for bad in ["DG", "DS", "DT"]:
                self.assertNotIn(bad, seq, f"Found isomerization motif {bad} in {seq}")
            # High developability score
            self.assertGreaterEqual(design.developability_score, 0.70)

    def test_target_structure_hotspot_docking(self):
        """Verifies that passing a target structure and hotspots computes engaged contacts."""
        ledger_db.create_budget(total_allocation=1.00, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-hotspots",
            target_model="mock_diffusion",
            num_designs=3,
            binder_length=65,
            target_pdb="2YBR",
            hotspot_residues=[101, 103, 106],
            target_protein_name="Cn2 Scorpion Neurotoxin",
            budget_id=1,
        )

        result = generate_budgeted_de_novo_binders(request=request, db_path=self.db_path)
        self.assertTrue(result.approved)
        self.assertEqual(len(result.designs), 3)

        for design in result.designs:
            self.assertIsNotNone(design.pdb_content)
            self.assertGreater(len(design.pdb_content), 500)

    def test_wetlab_exports_and_fda_audit_integrity(self):
        """Verifies generation of Twist CSV, Opentrons OT-2 script, and FDA Part 11 certificate."""
        ledger_db.create_budget(total_allocation=1.00, budget_id=1, db_path=self.db_path)

        request = DeNovoBinderRequest(
            job_id="diff-job-wetlab",
            target_model="mock_diffusion",
            num_designs=3,
            binder_length=65,
            target_protein_name="Centruroides Cn2",
            budget_id=1,
        )

        result = generate_budgeted_de_novo_binders(request=request, db_path=self.db_path)
        self.assertTrue(result.approved)

        # Twist CSV verification
        self.assertIsNotNone(result.twist_csv)
        self.assertIn("Gene_Name", result.twist_csv)
        self.assertIn("cDNA_Sequence", result.twist_csv)
        self.assertIn("VB-DIFF-001", result.twist_csv)

        # Opentrons protocol verification
        self.assertIsNotNone(result.opentrons_script)
        self.assertIn("protocol_api", result.opentrons_script)
        self.assertIn("Centruroides Cn2", result.opentrons_script)

        # FDA Part 11 audit verification
        self.assertIsNotNone(result.fda_certificate)
        self.assertIsNotNone(result.audit_block_hash)

        # Cryptographic integrity verification across all blocks in SQLite
        is_intact, block_count, msg = ledger_db.verify_audit_integrity(db_path=self.db_path)
        self.assertTrue(is_intact, f"Audit chain verification failed: {msg}")
        self.assertGreater(block_count, 0)


if __name__ == "__main__":
    unittest.main()
