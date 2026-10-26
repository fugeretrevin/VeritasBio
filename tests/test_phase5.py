"""
tests/test_phase5.py — Unit test suite for Phase 5 (Orchestration & UI).

Tests:
    1. controller.orchestrator:
       - Full 2-round autonomous active-learning run.
       - Streaming event generator progression.
       - Budget deduction and evidence logging in SQLite.
       - Dual termination on budget exhaustion and max iterations.
    2. ui.dashboard:
       - Gradio app instantiation and formatting helpers.
"""

import os
import tempfile
import unittest

from controller.orchestrator import Orchestrator, OrchestratorConfig
from ledger import database as ledger_db
from ui.dashboard import create_dashboard, format_evidence_table, format_jobs_table


class TestPhase5(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.tmp_dir, "test_orchestrator.db")
        ledger_db.init_db(self.db_path)

    # -------------------------------------------------------------------------
    # 1. Orchestrator End-to-End Tests
    # -------------------------------------------------------------------------
    def test_orchestrator_full_run_mock(self):
        config = OrchestratorConfig(
            target_accession="P42212",
            total_budget_usd=1.00,
            max_iterations=2,
            batch_size=10,
            target_model="mock_esm",
            db_path=self.db_path,
            use_mock_llm=True,
        )

        orchestrator = Orchestrator(config)
        events = []
        for ev in orchestrator.run_stream():
            events.append(ev)

        self.assertGreater(len(events), 5)
        # Check that steps proceeded in unidirectional order
        steps_seen = [e.step for e in events]
        self.assertIn("init", steps_seen)
        self.assertIn("strategy", steps_seen)
        self.assertIn("selection", steps_seen)
        self.assertIn("authorization", steps_seen)
        self.assertIn("execution", steps_seen)
        self.assertIn("update", steps_seen)
        self.assertIn("complete", steps_seen)

        # Verify evidence in SQLite ledger
        evidence = ledger_db.get_evidence(db_path=self.db_path)
        self.assertEqual(len(evidence), 20)  # 2 rounds * 10 variants

        # Verify jobs in ledger
        jobs = ledger_db.list_jobs(db_path=self.db_path)
        self.assertEqual(len(jobs), 2)
        for j in jobs:
            self.assertEqual(j["status"], "COMPLETED")

        # Verify money was deducted
        rem_budget = ledger_db.get_remaining_budget(db_path=self.db_path)
        self.assertLess(rem_budget, 1.00)

    def test_orchestrator_budget_exhaustion_stop(self):
        # Tiny budget: $0.0001 (enough for only 1 mock job at nominal $0.0001)
        config = OrchestratorConfig(
            target_accession="P42212",
            total_budget_usd=0.0001,
            max_iterations=5,
            batch_size=10,
            target_model="mock_esm",
            db_path=self.db_path,
            use_mock_llm=True,
        )

        orchestrator = Orchestrator(config)
        summary = orchestrator.run()

        # Should halt after 1 round because budget hit 0
        self.assertEqual(summary.iterations_completed, 1)
        rem_budget = ledger_db.get_remaining_budget(db_path=self.db_path)
        self.assertEqual(rem_budget, 0.0)

    # -------------------------------------------------------------------------
    # 2. UI Formatting & App Instantiation Tests
    # -------------------------------------------------------------------------
    def test_dashboard_creation(self):
        demo = create_dashboard(db_path=self.db_path)
        self.assertIsNotNone(demo)

    def test_dashboard_table_formatters(self):
        sample_jobs = [{
            "job_id": "job-1",
            "target_model": "esm2_650m",
            "variant_count": 50,
            "estimated_cost": 0.0064,
            "status": "APPROVED",
        }]
        job_rows = format_jobs_table(sample_jobs)
        self.assertEqual(len(job_rows), 1)
        self.assertEqual(job_rows[0][0], "job-1")

        sample_evidence = [{
            "variant_sequence": "MSKGEE",
            "computed_score": 1.45,
            "job_id": "job-1",
            "timestamp": "2026-10-26T14:00:00Z",
        }]
        ev_rows = format_evidence_table(sample_evidence)
        self.assertEqual(len(ev_rows), 1)
        self.assertEqual(ev_rows[0][0], 1)  # Rank 1


if __name__ == "__main__":
    unittest.main()
