"""
tests/test_mechanistic_cot.py — Unit tests for Meta Llama 3.1 Mechanistic Biophysical CoT.

Tests:
    1. BiophysicalHypothesis schema validation and chemistry sanitizer.
    2. SearchStrategy JSON deserialization containing mechanistic hypotheses.
    3. MockLLMWorker production of authentic structural hypotheses for TrpB & avGFP.
    4. compute_mechanistic_cot_bonus calculation and prior reward biasing.
    5. Prioritization of suggested chemistries in candidate variant generation.
    6. SHA-256 cryptographic audit certificate inclusion of biophysical rationale.
    7. format_mechanistic_hypotheses_html card rendering for Gradio dashboard.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from controller.gatekeeper import SingleMutation, Variant
from controller.prompts import (
    BiophysicalHypothesis,
    ResidueCluster,
    SearchStrategy,
    parse_strategy_response,
)
from ledger import database as ledger_db
from ui.dashboard import format_mechanistic_hypotheses_html
from workers.math_worker import compute_mechanistic_cot_bonus, generate_candidate_variants
from workers.mock_workers import MockLLMWorker
from workers.rag_worker import fetch_protein_metadata
from workers.wetlab_worker import generate_cryptographic_audit_certificate


class TestMechanisticCoT(unittest.TestCase):
    """Test suite for Meta Llama 3.1 Mechanistic Chain-of-Thought biophysical integration."""

    def test_biophysical_hypothesis_validation(self) -> None:
        """Verifies schema validation and chemistry sanitizer."""
        hyp = BiophysicalHypothesis(
            pocket_name="substrate_channel",
            target_residues=[184, 227],
            mechanism="Steric volume expansion to accommodate bulkier non-canonical indoles.",
            structural_rationale="Phe184 is 6.2A from the Lys87 PLP cofactor. Mutating to Ile relieves steric clash.",
            suggested_chemistries={"184": ["i", "v", "l"], "227": "mq"},
            epistatic_coupling="Positive synergy across loops 1 and 2.",
        )
        self.assertEqual(hyp.pocket_name, "substrate_channel")
        self.assertEqual(hyp.target_residues, [184, 227])
        # Chemistry sanitizer should uppercase and convert string "mq" into ["M", "Q"]
        self.assertEqual(hyp.suggested_chemistries["184"], ["I", "V", "L"])
        self.assertEqual(hyp.suggested_chemistries["227"], ["M", "Q"])

    def test_search_strategy_json_parsing_with_cot(self) -> None:
        """Verifies parse_strategy_response deserializes mechanistic_hypotheses cleanly."""
        raw_json = json.dumps({
            "target_protein": "TrpB",
            "target_accession": "TRPB_TM9D8",
            "hypothesis": "Modulating entrance loop dynamics around PLP cofactor.",
            "target_clusters": [
                {
                    "cluster_id": "loop_1",
                    "positions": [183, 184],
                    "rationale": "Active site entrance loop.",
                }
            ],
            "forbidden_positions": [83],
            "max_mutation_depth": 2,
            "mechanistic_hypotheses": [
                {
                    "pocket_name": "loop_1_entrance",
                    "target_residues": [183, 184],
                    "mechanism": "Active site channel dilation.",
                    "structural_rationale": "Phe184 forms a steric gate at the mouth of the catalytic cleft.",
                    "suggested_chemistries": {"184": ["I", "V"]},
                    "epistatic_coupling": "Gating residue.",
                }
            ],
        })

        strat = parse_strategy_response(raw_json, sequence_length=400)
        self.assertEqual(len(strat.mechanistic_hypotheses), 1)
        hyp = strat.mechanistic_hypotheses[0]
        self.assertEqual(hyp.pocket_name, "loop_1_entrance")
        self.assertEqual(hyp.suggested_chemistries["184"], ["I", "V"])

    def test_mock_llm_worker_produces_cot_hypotheses(self) -> None:
        """Verifies MockLLMWorker generates authentic structural hypotheses for both targets."""
        worker = MockLLMWorker()

        # 1. TrpB
        trpb_meta = fetch_protein_metadata("TRPB_TM9D8")
        strat_trpb = worker.generate_strategy(trpb_meta, round_number=2)
        self.assertTrue(len(strat_trpb.mechanistic_hypotheses) >= 1)
        trpb_hyp = strat_trpb.mechanistic_hypotheses[0]
        self.assertIn("184", trpb_hyp.suggested_chemistries)
        self.assertIn("I", trpb_hyp.suggested_chemistries["184"])
        self.assertIn("PLP", trpb_hyp.structural_rationale)

        # 2. avGFP
        gfp_meta = fetch_protein_metadata("P42212")
        strat_gfp = worker.generate_strategy(gfp_meta, round_number=2)
        self.assertTrue(len(strat_gfp.mechanistic_hypotheses) >= 1)
        gfp_hyp = strat_gfp.mechanistic_hypotheses[0]
        self.assertIn("203", gfp_hyp.suggested_chemistries)
        self.assertIn("Y", gfp_hyp.suggested_chemistries["203"])
        self.assertIn("pi-", gfp_hyp.mechanism.lower())

    def test_compute_mechanistic_cot_bonus(self) -> None:
        """Verifies that candidates matching Meta Llama 3.1 suggested chemistries receive prior reward bonus."""
        hyp = BiophysicalHypothesis(
            pocket_name="pocket",
            target_residues=[184, 227],
            mechanism="Steric relief",
            structural_rationale="Cofactor proximity",
            suggested_chemistries={"184": ["I", "V"], "227": ["M"]},
        )

        # Variant matching 1 suggested chemistry
        v1 = Variant(notation="F184I", mutations=[SingleMutation(wt_aa="F", position=184, mut_aa="I")])
        bonus1 = compute_mechanistic_cot_bonus(v1, [hyp])
        self.assertAlmostEqual(bonus1, 0.12, places=3)

        # Variant matching 2 suggested chemistries
        v2 = Variant(
            notation="F184I:V227M",
            mutations=[
                SingleMutation(wt_aa="F", position=184, mut_aa="I"),
                SingleMutation(wt_aa="V", position=227, mut_aa="M"),
            ],
        )
        bonus2 = compute_mechanistic_cot_bonus(v2, [hyp])
        self.assertAlmostEqual(bonus2, 0.24, places=3)

        # Non-matching variant
        v_none = Variant(notation="F184R", mutations=[SingleMutation(wt_aa="F", position=184, mut_aa="R")])
        bonus_none = compute_mechanistic_cot_bonus(v_none, [hyp])
        self.assertEqual(bonus_none, 0.0)

    def test_candidate_generation_prioritizes_cot_chemistries(self) -> None:
        """Verifies that suggested chemistries appear early in candidate generation."""
        hyp = BiophysicalHypothesis(
            pocket_name="substrate_loop",
            target_residues=[184],
            mechanism="Volume expansion",
            structural_rationale="Proximity to PLP",
            suggested_chemistries={"184": ["I", "V"]},
        )
        strat = SearchStrategy(
            target_protein="TrpB",
            target_accession="TRPB_TM9D8",
            hypothesis="Test CoT candidate prioritization",
            target_clusters=[ResidueCluster(cluster_id="c1", positions=[184], rationale="test")],
            forbidden_positions=[],
            max_mutation_depth=1,
            mechanistic_hypotheses=[hyp],
        )
        wt = "M" * 183 + "F" + "M" * 200

        candidates = generate_candidate_variants(strat, wt, max_candidates=50)
        top_notations = [c.notation for c in candidates[:4]]
        # F184I and F184V must be generated among the very first single mutants
        self.assertIn("F184I", top_notations)
        self.assertIn("F184V", top_notations)

    def test_fda_audit_certificate_incorporates_cot(self) -> None:
        """Verifies that FDA 21 CFR Part 11 audit certificates cryptographically hash and store CoT hypotheses."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = Path(tmp_dir) / "test_cot_audit.db"
            ledger_db.init_db(db_path)
            ledger_db.create_budget(1.0, budget_id=1, db_path=db_path)

            hyp = BiophysicalHypothesis(
                pocket_name="catalytic_gate",
                target_residues=[184],
                mechanism="Indole channel dilation",
                structural_rationale="Steric gating 6.2A from cofactor",
                suggested_chemistries={"184": ["I"]},
            )

            cert = generate_cryptographic_audit_certificate(
                db_path=db_path,
                target_protein="TRPB_TM9D8",
                mechanistic_hypotheses=[hyp],
            )

            self.assertIn("mechanistic_biophysical_rationale", cert)
            rationale_list = cert["mechanistic_biophysical_rationale"]
            self.assertEqual(len(rationale_list), 1)
            self.assertEqual(rationale_list[0]["pocket_name"], "catalytic_gate")
            self.assertEqual(len(cert["ledger_hash_sha256"]), 64)

    def test_dashboard_html_rendering(self) -> None:
        """Verifies format_mechanistic_hypotheses_html produces valid, rich HTML."""
        worker = MockLLMWorker()
        meta = fetch_protein_metadata("TRPB_TM9D8")
        strat = worker.generate_strategy(meta, round_number=2)

        html = format_mechanistic_hypotheses_html(strat.model_dump(), target_protein="TRPB_TM9D8")
        self.assertIn("META LLAMA 3.1 INSTRUCT", html)
        self.assertIn("Biophysical Mechanistic Chain-of-Thought", html)
        self.assertIn("inter_loop_lid_complex", html)
        self.assertIn("Pos 184", html)
        self.assertIn("PLP cofactor", html)
        self.assertIn("Zero Direct LLM Math", html)


if __name__ == "__main__":
    unittest.main()
