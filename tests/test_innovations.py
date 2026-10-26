"""
tests/test_innovations.py — Unit tests for Hackathon & Startup innovations.

Tests:
    1. Codon optimization for E. coli expression and restriction site removal.
    2. Twist / IDT CSV order sheet generation.
    3. Opentrons (OT-2) 96-well automated liquid-handling protocol generation.
    4. FDA 21 CFR Part 11 cryptographic audit certificate (SHA-256 ledger proof).
    5. 3Dmol WebGL HTML generation with authentic PDB coordinates.
    6. Plotly active-learning latent landscape generation.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from ledger import database as ledger_db
from ui.landscape_plot import generate_landscape_figure
from ui.structure_viewer import generate_3d_viewer_html, load_pdb_content
from workers.wetlab_worker import (
    codon_optimize_sequence,
    generate_cryptographic_audit_certificate,
    generate_opentrons_protocol,
    generate_twist_order_csv,
)


class TestInnovations(unittest.TestCase):
    """Unit tests for enterprise wet-lab translation and visual modules."""

    def test_codon_optimization(self) -> None:
        """Verifies amino acid to DNA codon translation and restriction site suppression."""
        prot = "MSKGEELFT"
        cdna = codon_optimize_sequence(prot)
        # Should be len(prot)*3 + 3 (stop codon) = 30 bp
        self.assertEqual(len(cdna), 30)
        self.assertTrue(cdna.endswith("TAA"))
        # Forbidden cloning restriction sites should not be present
        self.assertNotIn("GAATTC", cdna)
        self.assertNotIn("CTCGAG", cdna)

    def test_twist_order_csv_generation(self) -> None:
        """Verifies CSV generation for Twist Bioscience / IDT portals."""
        variants = [
            {"variant_sequence": "MSKGEELFT", "computed_score": 3.85, "mutations": "S65T"},
            {"variant_sequence": "MSKGEELFT", "computed_score": 1.20, "mutations": "T203Y"},
        ]
        csv_str = generate_twist_order_csv(variants)
        self.assertIn("Gene_Name,Target_Protein,Mutations,Fitness_Score", csv_str)
        self.assertIn("VB-GFP-001_S65T", csv_str)
        self.assertIn("3.8500", csv_str)

    def test_opentrons_protocol_generation(self) -> None:
        """Verifies executable Python protocol generation for Opentrons OT-2 robots."""
        names = ["GFP_Var_01", "GFP_Var_02", "GFP_Var_03"]
        script = generate_opentrons_protocol(names)
        self.assertIn("from opentrons import protocol_api", script)
        self.assertIn("opentrons_96_tiprack_20ul", script)
        self.assertIn("p20_single_gen2", script)
        self.assertIn("p300_multi_gen2", script)
        self.assertIn("GFP_Var_01", script)

    def test_fda_audit_certificate(self) -> None:
        """Verifies SHA-256 cryptographic audit certificate generation."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = Path(tmp_dir) / "test_audit.db"
            ledger_db.init_db(db_path)
            ledger_db.create_budget(1.0, budget_id=1, db_path=db_path)

            cert = generate_cryptographic_audit_certificate(db_path)
            self.assertTrue(cert["certificate_id"].startswith("VB-AUDIT-"))
            self.assertEqual(len(cert["ledger_hash_sha256"]), 64)
            self.assertIn("FDA 21 CFR Part 11", cert["standard_compliance"])
            self.assertEqual(cert["total_budget_usd"], 1.0)

    def test_3dmol_html_generation(self) -> None:
        """Verifies 3Dmol.js HTML block generation with 1EMA crystal structure."""
        pdb_content = load_pdb_content()
        self.assertTrue(len(pdb_content) > 50000, "PDB file should contain ATOM records")
        html_str = generate_3d_viewer_html([65, 163, 203])
        self.assertIn("<iframe srcdoc=", html_str)
        self.assertIn("3Dmol-min.js", html_str)
        self.assertIn("1EMA", html_str)

    def test_plotly_landscape_figure(self) -> None:
        """Verifies Plotly latent landscape figure generation."""
        fig = generate_landscape_figure([])
        self.assertIsNotNone(fig)
        self.assertEqual(fig.layout.paper_bgcolor, "#0b1329")


if __name__ == "__main__":
    unittest.main()
