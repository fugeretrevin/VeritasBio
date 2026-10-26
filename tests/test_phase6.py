"""
tests/test_phase6.py — Unit test suite for Phase 6 (Benchmark & Validation Engine).

Tests:
    1. scripts.evaluate:
       - Sequence to mutation string parser (singles, multi-point, wildtype).
       - Sarkisyan 2016 GFP DMS benchmark ingestion.
       - Spearman rank correlation (ρ) calculation against ledger evidence.
       - Top-decile hit rate and enrichment factor.
       - Cloud compute savings and efficiency multiple vs. brute-force.
"""

import json
import os
import tempfile
import unittest

from ledger import database as ledger_db
from scripts.evaluate import (
    DEFAULT_BENCHMARK_PATH,
    BenchmarkMetrics,
    evaluate_run,
    load_ground_truth_benchmark,
    sequence_to_mutation_str,
)
from workers.rag_worker import fetch_protein_metadata


class TestPhase6(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.tmp_dir, "test_eval.db")
        ledger_db.init_db(self.db_path)

        meta = fetch_protein_metadata("P42212")
        self.wt = meta.sequence

    # -------------------------------------------------------------------------
    # 1. Sequence-to-Mutation Parser Tests
    # -------------------------------------------------------------------------
    def test_sequence_to_mutation_wt(self):
        self.assertEqual(sequence_to_mutation_str(self.wt, self.wt), "WT")

    def test_sequence_to_mutation_single(self):
        # Mutate residue 65 (index 64): S -> T
        chars = list(self.wt)
        chars[64] = "T"
        mut_seq = "".join(chars)
        self.assertEqual(sequence_to_mutation_str(self.wt, mut_seq), "S65T")

    def test_sequence_to_mutation_double(self):
        # Mutate residue 65: S -> T, residue 203 (index 202): T -> Y
        chars = list(self.wt)
        chars[64] = "T"
        chars[202] = "Y"
        mut_seq = "".join(chars)
        self.assertEqual(sequence_to_mutation_str(self.wt, mut_seq), "S65T:T203Y")

    # -------------------------------------------------------------------------
    # 2. Benchmark Ingestion Tests
    # -------------------------------------------------------------------------
    def test_load_benchmark_data(self):
        bench = load_ground_truth_benchmark(str(DEFAULT_BENCHMARK_PATH))
        self.assertIn("S65T", bench)
        self.assertIn("S65T:T203Y", bench)
        self.assertIn("Y66A", bench)
        self.assertGreater(bench["S65T:T203Y"], bench["S65T"])
        self.assertLess(bench["Y66A"], 0.0)

    # -------------------------------------------------------------------------
    # 3. Evaluation Metrics Engine Tests
    # -------------------------------------------------------------------------
    def test_evaluate_run_with_correlated_evidence(self):
        # Set up a ledger run with budget
        ledger_db.create_budget(10.0, db_path=self.db_path)

        # Mutate sequences
        chars1 = list(self.wt)
        chars1[64] = "T"  # S65T (true score 1.85)
        seq_s65t = "".join(chars1)

        chars2 = list(self.wt)
        chars2[64] = "T"
        chars2[202] = "Y"  # S65T:T203Y (true score 4.45)
        seq_double = "".join(chars2)

        chars3 = list(self.wt)
        chars3[65] = "A"  # Y66A (true score -4.50)
        seq_y66a = "".join(chars3)

        # Authorize and record evidence
        job_id = "job-eval-001"
        ledger_db.authorize_spend(job_id, "mock_esm", 3, 0.0003, db_path=self.db_path)
        ledger_db.record_evidence(
            job_id,
            [
                (seq_double, 4.30),
                (seq_s65t, 1.90),
                (seq_y66a, -4.20),
            ],
            db_path=self.db_path,
        )
        ledger_db.mark_job_completed(job_id, db_path=self.db_path)

        metrics = evaluate_run(
            db_path=self.db_path,
            wt_sequence=self.wt,
            total_candidate_universe_size=1000,
        )

        self.assertIsInstance(metrics, BenchmarkMetrics)
        self.assertEqual(metrics.total_ledger_evidence, 3)
        self.assertEqual(metrics.matched_benchmark_variants, 3)
        # Perfect rank order -> Spearman rho = 1.0
        self.assertAlmostEqual(metrics.spearman_rho, 1.0, places=4)
        self.assertGreater(metrics.top_decile_hit_rate, 0.0)
        self.assertGreater(metrics.cloud_efficiency_multiple, 1.0)
        self.assertGreater(metrics.cloud_cost_savings_pct, 50.0)
        self.assertEqual(metrics.best_discovered_variant, "S65T:T203Y")

    def test_evaluate_run_empty_ledger_handles_gracefully(self):
        metrics = evaluate_run(db_path=self.db_path, wt_sequence=self.wt)
        self.assertEqual(metrics.total_ledger_evidence, 0)
        self.assertEqual(metrics.matched_benchmark_variants, 0)
        self.assertEqual(metrics.spearman_rho, 0.0)


if __name__ == "__main__":
    unittest.main()
