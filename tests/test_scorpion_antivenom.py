"""
tests/test_scorpion_antivenom.py — Test suite for Scorpion Antivenom Neutralizing Antibody Campaign.

Verifies:
1. Metadata & structural PDB integrity for SCFV_9004G_CN2 and PDB 2YBR.
2. Spatial distance mapping between scFv CDR-H3 paratope and Centruroides noxius Cn2 toxin.
3. Biologics developability liability screening on CDR-H3 variants.
4. Autonomous active learning campaign execution with ledger authorization.
5. Automated wet-lab exports (Twist CSV, Opentrons OT-2 script, FDA Part 11 certificate).
"""

import json
import os
import tempfile
import unittest
from pathlib import Path

from controller.gatekeeper import SingleMutation, Variant
from ledger import database as ledger_db
from scripts.run_scorpion_antivenom_campaign import (
    CDR_H3_POSITIONS,
    evaluate_ground_truth_affinity,
    generate_candidate_paratope_library,
    run_scorpion_antivenom_campaign,
)
from workers.math_worker import analyze_biologics_developability
from workers.rag_worker import fetch_protein_metadata
from workers.spatial_filter import SpatialContactMap

ROOT_DIR = Path(__file__).resolve().parent.parent


class TestScorpionAntivenom(unittest.TestCase):
    """Test suite for the scorpion antivenom engineering pipeline."""

    def test_metadata_and_pdb_integrity(self):
        """Verifies metadata ingestion and PDB 2YBR crystal structure loading."""
        meta = fetch_protein_metadata("SCFV_9004G_CN2")
        self.assertEqual(meta.accession, "SCFV_9004G_CN2")
        self.assertEqual(meta.pdb_id, "2YBR")
        self.assertEqual(meta.length, 117)
        self.assertIn(22, meta.catalytic_positions)
        self.assertIn(96, meta.catalytic_positions)
        self.assertIn(101, meta.binding_positions)

        pdb_file = ROOT_DIR / "workers" / "data" / "2YBR.pdb"
        self.assertTrue(pdb_file.is_file(), f"Missing 2YBR.pdb at {pdb_file}")
        self.assertGreater(pdb_file.stat().st_size, 100_000)

    def test_spatial_contact_calculation(self):
        """Verifies 3D spatial distance extraction between scFv VH and Cn2 toxin."""
        pdb_file = ROOT_DIR / "workers" / "data" / "2YBR.pdb"
        contact_map = SpatialContactMap(pdb_path=pdb_file, chain_id="A")
        self.assertGreater(len(contact_map.ca_coordinates), 100)

        # Check CDR-H3 positions exist in C-alpha coordinates
        for pos in [99, 100, 101, 102, 103, 104, 105, 106]:
            self.assertIn(pos, contact_map.ca_coordinates)

    def test_developability_filtering_on_paratope(self):
        """Verifies developability liability detection on scFv CDR-H3 mutations."""
        meta = fetch_protein_metadata("SCFV_9004G_CN2")
        wt_seq = meta.sequence

        # Clean mutation: V101F (hydrophobic pi-stacking)
        clean_var = Variant(notation="V101F", mutations=[SingleMutation(wt_aa="V", position=101, mut_aa="F")])
        clean_seq = clean_var.apply_to(wt_seq)
        clean_eval = analyze_biologics_developability(clean_seq, wt_seq)
        self.assertFalse(clean_eval.get("unpaired_cysteine", False))
        self.assertGreaterEqual(clean_eval.get("developability_score", 0.0), 0.65)

        # Dangerous mutation: V101C (introduces unpaired cysteine, breaking C22-C96)
        bad_var = Variant(notation="V101C", mutations=[SingleMutation(wt_aa="V", position=101, mut_aa="C")])
        bad_seq = bad_var.apply_to(wt_seq)
        bad_eval = analyze_biologics_developability(bad_seq, wt_seq)
        self.assertTrue(bad_eval.get("unpaired_cysteine", False))

        # Check ground truth affinity reflects penalty
        clean_affinity = evaluate_ground_truth_affinity(clean_seq, wt_seq)
        bad_affinity = evaluate_ground_truth_affinity(bad_seq, wt_seq)
        self.assertGreater(clean_affinity, 1.10)
        self.assertLess(bad_affinity, 0.20)

    def test_candidate_library_generation(self):
        """Verifies focused paratope library creation excludes unpaired cysteines."""
        meta = fetch_protein_metadata("SCFV_9004G_CN2")
        cands = generate_candidate_paratope_library(meta.sequence, max_candidates=100, seed=42)
        self.assertEqual(len(cands), 100)
        for c in cands:
            for m in c.mutations:
                self.assertNotEqual(m.mut_aa, "C", "Cysteine should not be introduced into paratope")
                self.assertIn(m.position, CDR_H3_POSITIONS)

    def test_scorpion_antivenom_campaign_fast(self):
        """Executes a 2-round antivenom active learning campaign end-to-end."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_db = os.path.join(tmp_dir, "test_scorp.db")
            results = run_scorpion_antivenom_campaign(
                total_budget_usd=0.25,
                rounds=2,
                batch_size=10,
                seed=42,
                db_path=temp_db,
                export_dir=tmp_dir,
            )

            # Verification of results
            self.assertEqual(results["target"], "SCFV_9004G_CN2")
            self.assertEqual(results["total_evaluated"], 30)  # 10 seed + 20 in 2 rounds
            self.assertGreater(results["best_score"], 0.90)
            self.assertGreaterEqual(results["developability_score"], 0.60)

            # Verification of physical export artifacts
            twist_csv = Path(results["twist_csv"])
            self.assertTrue(twist_csv.is_file(), "Twist CSV order must exist")
            self.assertGreater(twist_csv.stat().st_size, 200)

            ot2_script = Path(results["opentrons_script"])
            self.assertTrue(ot2_script.is_file(), "Opentrons OT-2 script must exist")
            self.assertGreater(ot2_script.stat().st_size, 500)

            fda_cert = Path(results["fda_certificate"])
            self.assertTrue(fda_cert.is_file(), "FDA Part 11 certificate must exist")
            with open(fda_cert, "r", encoding="utf-8") as f:
                cert_data = json.load(f)
            # The verifiable claim is that the chain re-verified, not that the
            # system is Part 11 compliant.
            self.assertTrue(cert_data.get("audit_chain_verified"))
            self.assertIn("verified", cert_data.get("standard_compliance", "").lower())
            self.assertIn("ledger_hash_sha256", cert_data)
            self.assertTrue(cert_data.get("merkle_audit_chain_valid", False))

    def test_ui_and_3d_viewer_integration(self):
        """Verifies UI target resolution and 3Dmol viewer generation for scorpion antivenom."""
        from ui.structure_viewer import resolve_target_pdb_context, generate_3d_viewer_html

        pdb_path, title, label, cat_res, cat_label, color = resolve_target_pdb_context("SCFV_9004G_CN2")
        self.assertEqual(pdb_path.name, "2YBR.pdb")
        self.assertIn("2YBR", title)
        self.assertEqual(color, "#d97706")
        self.assertIn(101, cat_res)

        html = generate_3d_viewer_html(target_accession="SCFV_9004G_CN2")
        self.assertIn("2YBR", html)
        self.assertIn("#d97706", html)

    def test_orchestrator_execution_on_scorpion_target(self):
        """Verifies autonomous 6-step loop orchestrator runs on SCFV_9004G_CN2."""
        from controller.orchestrator import Orchestrator, OrchestratorConfig

        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_db = os.path.join(tmp_dir, "test_orch_scorp.db")
            ledger_db.init_db(temp_db)
            config = OrchestratorConfig(
                target_accession="SCFV_9004G_CN2",
                total_budget_usd=0.50,
                max_iterations=1,
                batch_size=5,
                target_model="mock_esm",
                use_mock_llm=True,
                db_path=temp_db,
            )
            orch = Orchestrator(config)
            summary = orch.run()
            self.assertEqual(summary.iterations_completed, 1)
            evidence = ledger_db.get_evidence(db_path=temp_db)
            self.assertGreaterEqual(len(evidence), 5)


if __name__ == "__main__":
    unittest.main()
