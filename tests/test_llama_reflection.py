"""
tests/test_llama_reflection.py — Unit Tests for Gap 1: Closed-Loop Residual Reflection.

Tests:
    1. Validation of ResidualOutlier and ReflectionDiagnosis Pydantic schemas.
    2. Prompt construction and fail-closed JSON parsing for Llama 3.1 reflections.
    3. MockLLMWorker deterministic reflection generation.
    4. Closed-loop Orchestrator integration: verifies that prediction residuals are computed,
       reflection diagnoses are formed, and dynamic hypotheses are merged into the next round.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from controller.orchestrator import Orchestrator, OrchestratorConfig
from controller.prompts import (
    BiophysicalHypothesis,
    ReflectionDiagnosis,
    ResidualOutlier,
    build_reflection_prompt,
    get_reflection_system_prompt,
    parse_reflection_response,
)
from workers.mock_workers import MockLLMWorker
from workers.rag_worker import fetch_protein_metadata


class TestLlamaReflection(unittest.TestCase):
    """Verifies Closed-Loop Residual Reflection schemas, prompts, and orchestration."""

    def test_residual_outlier_and_reflection_diagnosis_schemas(self) -> None:
        """Verifies schema validation, positional chemistries, and cleaning logic."""
        outlier = ResidualOutlier(
            variant_notation="F184I:S228G",
            predicted_score=0.4200,
            observed_score=0.8005,
            residual=0.3805,
            distance_angstroms=18.95,
            mutation_details=["F184I (aromatic->branched)", "S228G (polar->flexible)"],
        )
        self.assertEqual(outlier.variant_notation, "F184I:S228G")
        self.assertAlmostEqual(outlier.residual, 0.3805, places=4)
        self.assertEqual(len(outlier.mutation_details), 2)

        diagnosis = ReflectionDiagnosis(
            target_protein="TrpB TM9D8",
            round_number=2,
            positive_outliers=[outlier],
            negative_outliers=[],
            mechanistic_diagnosis="F184I relieved steric hindrance in the substrate entrance channel.",
            updated_suggested_chemistries={"184": ["I", "V"], "228": ["G", "A"]},
            forbidden_chemistries={"183": ["K", "R"]},
            structural_cooperativity_insight="Coupled loops 183-184 and 227-228 relax to open active-site channel.",
        )
        self.assertEqual(diagnosis.round_number, 2)
        self.assertIn("184", diagnosis.updated_suggested_chemistries)
        self.assertIn("I", diagnosis.updated_suggested_chemistries["184"])
        self.assertIn("K", diagnosis.forbidden_chemistries["183"])

    def test_reflection_prompt_construction_and_parsing(self) -> None:
        """Verifies prompt construction and strict JSON parsing."""
        sys_prompt = get_reflection_system_prompt()
        self.assertIn("Closed-Loop Residual Reflection", sys_prompt)

        pos_outliers = [
            ResidualOutlier(
                variant_notation="S65T:T203Y",
                predicted_score=1.40,
                observed_score=2.25,
                residual=0.85,
                distance_angstroms=7.8,
                mutation_details=["S65T", "T203Y"],
            )
        ]
        neg_outliers = [
            ResidualOutlier(
                variant_notation="Y66K",
                predicted_score=0.20,
                observed_score=-1.80,
                residual=-2.00,
                distance_angstroms=None,
                mutation_details=["Y66K (chromophore aromatic to charged)"],
            )
        ]

        user_prompt = build_reflection_prompt(
            target_protein="avGFP",
            round_number=1,
            positive_outliers=pos_outliers,
            negative_outliers=neg_outliers,
            protein_summary="avGFP 238 AA fluorescent protein",
            pdb_id="1EMA",
        )
        self.assertIn("S65T:T203Y", user_prompt)
        self.assertIn("Y66K", user_prompt)
        self.assertIn("Required ReflectionDiagnosis JSON", user_prompt)

        # Test fail-closed parser
        sample_json = json.dumps({
            "target_protein": "avGFP",
            "round_number": 1,
            "mechanistic_diagnosis": "S65T accelerates chromophore maturation; T203Y provides pi-stacking.",
            "updated_suggested_chemistries": {"65": ["T"], "203": ["Y", "F"]},
            "forbidden_chemistries": {"66": ["K", "R", "D", "E"]},
            "structural_cooperativity_insight": "Chromophore environment stabilizes fluorophore conjugated ring.",
        })
        diag = parse_reflection_response(sample_json)
        self.assertEqual(diag.target_protein, "avGFP")
        self.assertEqual(diag.updated_suggested_chemistries["65"], ["T"])

        # Test malformed JSON fail-closed rejection
        with self.assertRaises(ValueError):
            parse_reflection_response("Not valid JSON here!")

    def test_mock_llm_worker_reflection(self) -> None:
        """Verifies MockLLMWorker deterministic reflection synthesis."""
        mock = MockLLMWorker()
        metadata = fetch_protein_metadata("P42212", use_cache=True)
        pos = [
            ResidualOutlier(
                variant_notation="S65T",
                predicted_score=0.3,
                observed_score=0.85,
                residual=0.55,
                mutation_details=["S65T"],
            )
        ]
        neg = [
            ResidualOutlier(
                variant_notation="G67A",
                predicted_score=0.1,
                observed_score=-2.2,
                residual=-2.3,
                mutation_details=["G67A"],
            )
        ]
        ref = mock.generate_reflection(metadata, round_number=1, positive_outliers=pos, negative_outliers=neg)
        self.assertIsInstance(ref, ReflectionDiagnosis)
        self.assertEqual(ref.round_number, 1)
        self.assertIn("S65T", ref.mechanistic_diagnosis)
        self.assertIn("65", ref.updated_suggested_chemistries)
        self.assertIn("67", ref.forbidden_chemistries)

    def test_orchestrator_closed_loop_reflection_execution(self) -> None:
        """Verifies that Orchestrator executes reflection and yields reflection_complete event."""
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
            db_path = f.name

        config = OrchestratorConfig(
            target_accession="P42212",
            total_budget_usd=1.00,
            batch_size=10,
            max_iterations=2,
            target_model="mock_esm",
            use_mock_llm=True,
            enable_reflection=True,
            reflection_top_k=2,
            db_path=db_path,
        )

        orch = Orchestrator(config)
        steps_seen = []
        for event in orch.run_stream():
            steps_seen.append(event.step)

        # Confirm that reflection completed and was recorded in telemetry
        self.assertIn("reflection_complete", steps_seen)
        self.assertIsNotNone(orch.last_reflection)
        self.assertIsInstance(orch.last_reflection, ReflectionDiagnosis)
        self.assertTrue(len(orch.last_reflection.mechanistic_diagnosis) > 10)

        # Confirm that in subsequent strategy, hypotheses are enriched by reflection
        strategy = orch._generate_strategy(orch.last_metadata, round_number=2, top_evidence=[])
        hyp_names = [h.pocket_name for h in strategy.mechanistic_hypotheses]
        self.assertTrue(any(f"reflection_round_{orch.last_reflection.round_number}" in name for name in hyp_names))


if __name__ == "__main__":
    unittest.main()
