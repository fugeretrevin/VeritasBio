"""
tests/test_phase2.py — Unit tests for Phase 2 (Biological Priors & Combinatorial Prompts).

Tests:
    1. workers.rag_worker: UniProt metadata retrieval, sequence length, chromophore parsing, LLM summary generation.
    2. controller.prompts: SearchStrategy validation, forbidden overlap checks, out-of-bounds position detection, prompt rendering.
"""

import json
import unittest

from controller.prompts import (
    FEW_SHOT_EXAMPLE_JSON,
    ResidueCluster,
    SearchStrategy,
    build_strategy_user_prompt,
    get_strategy_system_prompt,
    parse_strategy_response,
)
from workers.rag_worker import (
    FeatureAnnotation,
    ProteinMetadata,
    fetch_protein_metadata,
    parse_uniprot_response,
)


class TestPhase2(unittest.TestCase):

    # -------------------------------------------------------------------------
    # 1. RAG Worker Tests (UniProt avGFP P42212)
    # -------------------------------------------------------------------------
    def test_fetch_metadata_cached(self):
        metadata = fetch_protein_metadata("P42212", use_cache=True)
        self.assertEqual(metadata.accession, "P42212")
        self.assertEqual(metadata.length, 238)
        self.assertEqual(len(metadata.sequence), 238)
        self.assertTrue(metadata.sequence.startswith("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT"))

        # Verify chromophore residue indices (65, 66, 67)
        self.assertIn(65, metadata.chromophore_residues)
        self.assertIn(66, metadata.chromophore_residues)
        self.assertIn(67, metadata.chromophore_residues)

    def test_metadata_summary_rendering(self):
        metadata = fetch_protein_metadata("P42212", use_cache=True)
        summary = metadata.to_llm_summary()
        self.assertIn("P42212", summary)
        self.assertIn("Sequence Length: 238", summary)
        self.assertIn("Chromophore Triad", summary)

    def test_parse_uniprot_missing_sequence(self):
        bad_data = {"primaryAccession": "BAD", "sequence": {}}
        with self.assertRaises(ValueError):
            parse_uniprot_response(bad_data)

    # -------------------------------------------------------------------------
    # 2. Prompts & Strategy Schema Tests
    # -------------------------------------------------------------------------
    def test_parse_valid_strategy_json(self):
        raw_json = json.dumps(FEW_SHOT_EXAMPLE_JSON)
        strategy = parse_strategy_response(raw_json, sequence_length=238)
        self.assertEqual(strategy.target_protein, "avGFP")
        self.assertEqual(len(strategy.target_clusters), 2)
        self.assertEqual(strategy.target_clusters[0].cluster_id, "chromophore_core")
        self.assertIn(65, strategy.target_clusters[0].positions)

    def test_parse_strategy_with_code_fences(self):
        # Even if LLM wraps in ```json ... ```, parser should clean it
        raw_with_fences = f"```json\n{json.dumps(FEW_SHOT_EXAMPLE_JSON)}\n```"
        strategy = parse_strategy_response(raw_with_fences, sequence_length=238)
        self.assertEqual(strategy.target_accession, "P42212")

    def test_parse_strategy_out_of_bounds(self):
        payload = dict(FEW_SHOT_EXAMPLE_JSON)
        # Position 999 exceeds sequence length 238
        payload["target_clusters"] = [
            {"cluster_id": "out_of_bounds", "positions": [999], "rationale": "invalid"}
        ]
        with self.assertRaises(ValueError) as cm:
            parse_strategy_response(json.dumps(payload), sequence_length=238)
        self.assertIn("exceeds sequence length", str(cm.exception))

    def test_forbidden_position_overlap(self):
        payload = dict(FEW_SHOT_EXAMPLE_JSON)
        # Cluster contains position 65, which is also set as forbidden
        payload["forbidden_positions"] = [65]
        with self.assertRaises(ValueError) as cm:
            parse_strategy_response(json.dumps(payload), sequence_length=238)
        self.assertIn("forbidden positions", str(cm.exception))

    def test_user_prompt_construction(self):
        metadata = fetch_protein_metadata("P42212", use_cache=True)
        summary = metadata.to_llm_summary()
        top_evidence = [
            {"variant_sequence": "S65T:Y66H", "computed_score": 1.45},
            {"variant_sequence": "T203Y", "computed_score": 1.22},
        ]
        user_prompt = build_strategy_user_prompt(
            protein_summary=summary,
            round_number=2,
            top_evidence=top_evidence,
        )
        self.assertIn("Round 2", user_prompt)
        self.assertIn("S65T:Y66H", user_prompt)
        self.assertIn("P42212", user_prompt)

    def test_system_prompt_immutability(self):
        sys_prompt = get_strategy_system_prompt()
        self.assertIn("Strategic Biological Planner", sys_prompt)
        self.assertIn("raw JSON", sys_prompt)


if __name__ == "__main__":
    unittest.main()
