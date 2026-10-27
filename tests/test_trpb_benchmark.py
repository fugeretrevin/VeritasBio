"""
tests/test_trpb_benchmark.py — Unit & Integration tests for the authentic TrpB fitness landscape.

Guards the integrity of the prepared dataset as much as the loop that consumes it: if the
upstream file or the conversion ever drifts, the combinatorial-completeness and parent-lineage
assertions below fail loudly rather than quietly producing a plausible-looking benchmark.
"""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from scripts.download_trpb_benchmark import (
    ACTIVE_SITE_POSITIONS,
    EXPECTED_PARENT_SUBSTITUTIONS,
    PARENT_CACHE_JSON,
)
from scripts.run_trpb_benchmark import (
    DATA_CSV,
    load_fitness_library,
    load_parent_fitness,
    run_trpb_active_learning,
)
from workers.rag_worker import fetch_protein_metadata


class TestTrpBDatasetIntegrity(unittest.TestCase):
    def test_csv_file_exists_and_parses(self):
        """The prepared TrpB benchmark exists and holds the full published library."""
        self.assertTrue(
            DATA_CSV.exists(),
            "TrpB benchmark CSV not found — run: python scripts/download_trpb_benchmark.py",
        )
        lib = load_fitness_library(DATA_CSV)
        # 20^4 combinations minus the parent itself (which carries no mutations).
        self.assertEqual(len(lib), 159_999)

    def test_library_is_combinatorially_complete(self):
        """
        Every variant must be 1-4 substitutions confined to the saturated sites, and the counts
        must match 19^k exactly. This is the strongest available check that the conversion
        neither dropped nor invented a single measurement.
        """
        order_counts: Counter[int] = Counter()
        positions_seen = set()
        with open(DATA_CSV, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                muts = row["mutant"].split(":")
                order_counts[len(muts)] += 1
                for m in muts:
                    positions_seen.add(int(m[1:-1]))

        self.assertEqual(sorted(positions_seen), ACTIVE_SITE_POSITIONS)
        self.assertEqual(order_counts[1], 4 * 19)           # 76
        self.assertEqual(order_counts[2], 6 * 19 ** 2)      # 2,166
        self.assertEqual(order_counts[3], 4 * 19 ** 3)      # 27,436
        self.assertEqual(order_counts[4], 19 ** 4)          # 130,321
        self.assertEqual(sum(order_counts.values()), 159_999)

    def test_scores_are_unaltered_and_in_published_range(self):
        """Fitness values stay inside the published range; nothing is rescaled or clipped."""
        lib = load_fitness_library(DATA_CSV)
        vals = list(lib.values())
        self.assertAlmostEqual(min(vals), -0.164070661, places=9)
        self.assertAlmostEqual(max(vals), 1.0, places=9)

    def test_bin_column_follows_documented_rule(self):
        """DMS_score_bin is exactly 'score >= parent fitness', with no other thresholding."""
        parent = load_parent_fitness()
        with open(DATA_CSV, "r", encoding="utf-8") as f:
            for i, row in enumerate(csv.DictReader(f)):
                expected = 1 if float(row["DMS_score"]) >= parent else 0
                self.assertEqual(int(row["DMS_score_bin"]), expected, f"row {i}: {row['mutant']}")
                if i >= 5000:
                    break


class TestTrpBParentMetadata(unittest.TestCase):
    def test_parent_record_is_labelled_as_derived(self):
        """The offline cache must declare its provenance rather than pose as a UniProt entry."""
        record = json.loads(PARENT_CACHE_JSON.read_text(encoding="utf-8"))
        prov = record["veritasbio_provenance"]
        self.assertIn("derived", prov["record_type"])
        self.assertEqual(len(prov["substitutions_vs_uniprot_wild_type"]), EXPECTED_PARENT_SUBSTITUTIONS)
        self.assertEqual(prov["active_site_positions"], ACTIVE_SITE_POSITIONS)

    def test_rag_worker_resolves_parent_offline(self):
        """fetch_protein_metadata must serve the Tm9D8* parent from the bundled cache."""
        meta = fetch_protein_metadata("TRPB_TM9D8", use_cache=True)
        self.assertEqual(len(meta.sequence), 397)
        self.assertTrue(meta.sequence.endswith("LEHHHHHH"))
        self.assertIn("Thermotoga", meta.organism)

    def test_mutation_notation_matches_parent_sequence(self):
        """
        Each mutant's stated wild-type residue must match the parent sequence at that position,
        which is what Variant.validate_against_wt enforces downstream.
        """
        wt = fetch_protein_metadata("TRPB_TM9D8", use_cache=True).sequence
        with open(DATA_CSV, "r", encoding="utf-8") as f:
            for i, row in enumerate(csv.DictReader(f)):
                for m in row["mutant"].split(":"):
                    self.assertEqual(wt[int(m[1:-1]) - 1], m[0], f"row {i}: {row['mutant']}")
                if i >= 5000:
                    break


class TestTrpBActiveLearning(unittest.TestCase):
    def test_trpb_active_learning_simulation(self):
        """Runs a fast 2-round TrpB active learning test on a temporary ledger."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_db = Path(tmp_dir) / "test_trpb.db"
            metrics = run_trpb_active_learning(
                total_budget_usd=0.50,
                rounds=2,
                batch_size=10,
                pool_size=3000,
                db_path=str(temp_db),
            )
            self.assertEqual(metrics["screened"], 20)
            self.assertGreater(metrics["cost_savings_pct"], 80.0)
            self.assertGreaterEqual(metrics["spearman_rho"], -1.0)
            self.assertLessEqual(metrics["spearman_rho"], 1.0)
            # The held-out set must be genuinely disjoint from what was screened.
            self.assertEqual(metrics["holdout_n"], 3000 - 20)
            # A random baseline is required for the hit count to mean anything.
            self.assertGreater(metrics["random_baseline_hits"], 0.0)


if __name__ == "__main__":
    unittest.main()
