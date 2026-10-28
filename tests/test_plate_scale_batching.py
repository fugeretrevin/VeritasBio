"""
tests/test_plate_scale_batching.py — Unit tests for Gap 2: Dual-Mode Plate-Scale Asynchronous Batching.
"""

from __future__ import annotations

import csv
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np

from controller.orchestrator import CampaignMode, Orchestrator, OrchestratorConfig
from ledger import database as ledger_db
from scripts.ingest_plate_data import ingest_plate_csv, parse_plate_csv
from workers.math_worker import greedy_dpp_sample
from workers.wetlab_worker import generate_plate_layout


class TestPlateScaleBatching(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = Path(self.temp_dir.name) / "test_plate.db"
        ledger_db.init_db(self.db_path)
        ledger_db.create_budget(total_allocation=10.00, db_path=self.db_path)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_campaign_mode_and_plate_capacity(self) -> None:
        """Verifies plate capacity math across 96-well and 384-well formats with anti-edge barriers."""
        # 96-well with barrier: 36 barrier wells + 4 control wells -> 56 candidate wells
        cfg_96_barrier = OrchestratorConfig(
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="96_well",
            anti_edge_barrier=True,
        )
        self.assertEqual(cfg_96_barrier.get_plate_candidate_capacity(), 56)

        # 96-well without barrier: 8 control wells -> 88 candidate wells
        cfg_96_nobarrier = OrchestratorConfig(
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="96_well",
            anti_edge_barrier=False,
        )
        self.assertEqual(cfg_96_nobarrier.get_plate_candidate_capacity(), 88)

        # 384-well with barrier: 76 barrier wells -> 308 candidate wells
        cfg_384_barrier = OrchestratorConfig(
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="384_well",
            anti_edge_barrier=True,
        )
        self.assertEqual(cfg_384_barrier.get_plate_candidate_capacity(), 308)

        # 384-well without barrier: 16 control wells -> 368 candidate wells
        cfg_384_nobarrier = OrchestratorConfig(
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="384_well",
            anti_edge_barrier=False,
        )
        self.assertEqual(cfg_384_nobarrier.get_plate_candidate_capacity(), 368)

    def test_greedy_dpp_sample_plate_scale_diversity(self) -> None:
        """Verifies k-DPP subset selection selects 56 and 88 mutually repelling, high-quality candidates."""
        np.random.seed(42)
        n_candidates = 600
        feature_dim = 15

        # Create clustered features
        features = np.random.randn(n_candidates, feature_dim)
        quality_scores = np.random.uniform(0.1, 1.0, size=n_candidates)

        # Test selecting 56 candidates (96-well plate with barrier)
        selected_56 = greedy_dpp_sample(
            quality_scores=quality_scores,
            features=features,
            batch_size=56,
            sigma_feat=1.5,
        )
        self.assertEqual(len(selected_56), 56)
        self.assertEqual(len(set(selected_56)), 56)

        # Test selecting 88 candidates (96-well full plate)
        selected_88 = greedy_dpp_sample(
            quality_scores=quality_scores,
            features=features,
            batch_size=88,
            sigma_feat=1.5,
        )
        self.assertEqual(len(selected_88), 88)
        self.assertEqual(len(set(selected_88)), 88)

        # Verify quality of selected candidates is significantly higher than random
        mean_selected_q = float(np.mean(quality_scores[selected_56]))
        mean_random_q = float(np.mean(quality_scores))
        self.assertGreater(mean_selected_q, mean_random_q)

    def test_orchestrator_plate_scale_execution_and_artifact_export(self) -> None:
        """Verifies Orchestrator runs in PLATE_SCALE_CRO mode, sizing to 56 variants and exporting wetlab files."""
        config = OrchestratorConfig(
            target_accession="P42212",
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="96_well",
            anti_edge_barrier=True,
            target_model="mock_esm",
            db_path=str(self.db_path),
            max_iterations=1,
            enable_live_llm=False,
        )
        orch = Orchestrator(config)
        events = list(orch.run_stream())

        # Check completed round
        update_events = [e for e in events if e.step == "update"]
        self.assertEqual(len(update_events), 1)

        # Verify wet-lab artifacts were generated
        self.assertIsNotNone(orch.last_wetlab_artifacts)
        artifacts = orch.last_wetlab_artifacts
        self.assertIn("plate_layout", artifacts)
        self.assertIn("twist_order_csv", artifacts)
        self.assertIn("opentrons_script", artifacts)
        self.assertIn("plate_csv", artifacts)

        layout = artifacts["plate_layout"]
        self.assertEqual(layout["barrier_wells_count"], 36)
        self.assertEqual(layout["assigned_variants_count"], 56)
        self.assertIn("PBS_Barrier", artifacts["plate_csv"])

    def test_orchestrator_pause_for_wetlab_ingestion(self) -> None:
        """Verifies that pause_for_wetlab_ingestion halts loop at waiting_for_assay."""
        config = OrchestratorConfig(
            target_accession="P42212",
            campaign_mode=CampaignMode.PLATE_SCALE_CRO,
            plate_format="96_well",
            anti_edge_barrier=True,
            pause_for_wetlab_ingestion=True,
            target_model="mock_esm",
            db_path=str(self.db_path),
            max_iterations=2,
        )
        orch = Orchestrator(config)
        events = list(orch.run_stream())

        # Must halt at waiting_for_assay without executing expensive cloud model
        wait_events = [e for e in events if e.step == "waiting_for_assay"]
        self.assertEqual(len(wait_events), 1)
        self.assertIn("Paused for wet-lab assay", wait_events[0].message)

        # Ensure no execution or update events occurred
        exec_events = [e for e in events if e.step in ("execution", "update")]
        self.assertEqual(len(exec_events), 0)

    def test_ingest_plate_csv_workflow(self) -> None:
        """Verifies end-to-end ingestion of plate reader experimental data into ledger."""
        # Generate a mock 96-well plate layout with experimental readings
        mock_variants = [
            {"mutations": f"S65T:V68L", "variant_sequence": "dummy", "computed_score": 1.20},
            {"mutations": f"T203Y", "variant_sequence": "dummy", "computed_score": 1.15},
            {"mutations": f"F64L:S65T:T203Y", "variant_sequence": "dummy", "computed_score": 2.40},
        ]
        layout = generate_plate_layout(mock_variants, plate_format="96_well", anti_edge_barrier=True)
        raw_csv = layout["csv_layout"]

        # Simulate adding a 'Reading' column to the plate CSV using csv writer
        in_reader = csv.reader(io.StringIO(raw_csv))
        out_buf = io.StringIO()
        out_writer = csv.writer(out_buf)
        header = next(in_reader) + ["Reading"]
        out_writer.writerow(header)

        for row in in_reader:
            if not row:
                continue
            sample_type = row[3]
            sample_name = row[4]
            if sample_type == "evaporation_barrier":
                val = 0.05
            elif sample_type == "control_wt":
                val = 1.00  # WT parent baseline
            elif sample_type == "control_blank":
                val = 0.02
            elif "F64L" in sample_name:
                val = 2.85  # High-activity winner
            elif "S65T" in sample_name:
                val = 1.65
            elif "T203Y" in sample_name:
                val = 1.40
            else:
                val = 0.85
            out_writer.writerow(row + [str(val)])

        plate_file = Path(self.temp_dir.name) / "mock_plate_reader.csv"
        plate_file.write_text(out_buf.getvalue(), encoding="utf-8")

        # Ingest the plate data
        ingest_result = ingest_plate_csv(
            csv_path=plate_file,
            db_path=self.db_path,
            target_accession="P42212",
            normalize_to_wt=True,
        )

        self.assertEqual(ingest_result["variants_ingested"], 3)
        self.assertEqual(ingest_result["barrier_wells_ignored"], 36)
        self.assertAlmostEqual(ingest_result["wt_baseline_score"], 1.00)
        self.assertAlmostEqual(ingest_result["top_fitness_discovered"], 2.85)
        self.assertIn("F64L", ingest_result["top_variant"])

        # Verify evidence is recorded in SQLite ledger
        evidence = ledger_db.get_evidence(db_path=self.db_path)
        self.assertEqual(len(evidence), 3)

        # Verify jobs in ledger
        jobs = ledger_db.list_jobs(db_path=self.db_path)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0]["status"], "COMPLETED")
        self.assertEqual(jobs[0]["target_model"], "wetlab_plate_reader")


if __name__ == "__main__":
    unittest.main()
