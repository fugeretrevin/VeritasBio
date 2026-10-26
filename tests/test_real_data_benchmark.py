"""
tests/test_real_data_benchmark.py — Unit & Integration tests for authentic ProteinGym real data.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.evaluate import load_ground_truth_benchmark, sequence_to_mutation_str
from scripts.run_real_data_benchmark import DATA_CSV, load_dms_library, run_real_active_learning


class TestRealDataBenchmark(unittest.TestCase):
    def test_csv_file_exists_and_parses(self):
        """Verifies that the authentic Sarkisyan 2016 ProteinGym dataset exists and has >50k rows."""
        self.assertTrue(DATA_CSV.exists(), "ProteinGym GFP CSV not found")
        lib = load_dms_library(DATA_CSV)
        self.assertGreaterEqual(len(lib), 50000, "Should load at least 50,000 real empirical variants")
        # Landmark single substitutions
        self.assertIn("S65T", lib)
        self.assertIn("T203A", lib)
        self.assertIn("V163A", lib)
        self.assertAlmostEqual(lib["S65T"], 3.5504, places=2)

    def test_evaluate_csv_loader(self):
        """Verifies that evaluate.py load_ground_truth_benchmark parses CSV files."""
        gt = load_ground_truth_benchmark(str(DATA_CSV))
        self.assertIn("S65T", gt)
        self.assertIn("K3R", gt)

    def test_sequence_to_mutation_str_sorting(self):
        """Verifies canonical position sorting in mutation strings."""
        wt = "ACDEFGHIKLMNPQRSTVWY"  # 20 aa mock
        # Mutate pos 10 (L->A) and pos 2 (C->G)
        mut = "AGDEFGHIKAMNPQRSTVWY"
        notation = sequence_to_mutation_str(wt, mut)
        # Should be sorted: C2G:L10A
        self.assertEqual(notation, "C2G:L10A")

    def test_real_active_learning_simulation(self):
        """Runs a fast 2-round real-data active learning test on a temporary ledger."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_db = Path(tmp_dir) / "test_real.db"
            metrics = run_real_active_learning(
                total_budget_usd=0.50,
                rounds=2,
                batch_size=10,
                db_path=str(temp_db),
            )
            self.assertEqual(metrics["screened"], 20)
            self.assertGreater(metrics["cost_savings_pct"], 80.0)
            self.assertGreaterEqual(metrics["spearman_rho"], -1.0)
            self.assertLessEqual(metrics["spearman_rho"], 1.0)


if __name__ == "__main__":
    unittest.main()
