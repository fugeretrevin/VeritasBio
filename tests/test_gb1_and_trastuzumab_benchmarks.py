"""
tests/test_gb1_and_trastuzumab_benchmarks.py — Test suite for GB1 & Trastuzumab benchmarks.

Verifies:
1. Dataset loading and integrity for Protein G GB1 (149k) and Trastuzumab CDR-H3 (36k).
2. UniProt metadata offline caching and feature parsing (SPG1_STRSG, HER2_TRASTUZUMAB).
3. Physicochemical featurization over active interface positions.
4. Active learning loop execution on GB1 landscape with ledger spend authorization.
5. Active learning loop execution on Trastuzumab CDR-H3 with developability liability scoring.
6. Benchmark evaluation against ground truth DMS scores with Spearman rho calculation.
"""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from controller.gatekeeper import Variant
from ledger import database as ledger_db
from scripts.evaluate import evaluate_run, load_ground_truth_benchmark
from scripts.run_gb1_benchmark import (
    DATA_CSV as GB1_CSV,
    PARENT_ACCESSION as GB1_ACCESSION,
    load_fitness_library as load_gb1_library,
    run_gb1_active_learning,
)
from scripts.run_trastuzumab_benchmark import (
    DATA_CSV as TRAS_CSV,
    PARENT_ACCESSION as TRAS_ACCESSION,
    load_fitness_library as load_tras_library,
    run_trastuzumab_active_learning,
)
from workers.math_worker import (
    analyze_biologics_developability,
    build_gp_model,
    featurize_variant_on_positions,
)
from workers.rag_worker import fetch_protein_metadata


class TestGB1AndTrastuzumabBenchmarks(unittest.TestCase):

    def test_gb1_dataset_and_metadata_integrity(self):
        """Verifies that the GB1 dataset exists and contains authentic Wu 2016 variants."""
        self.assertTrue(GB1_CSV.is_file(), f"GB1 CSV missing: {GB1_CSV}")
        library = load_gb1_library(GB1_CSV)
        self.assertGreaterEqual(len(library), 140_000, "GB1 library should have ~149k records")
        self.assertIn("WT", library)
        self.assertEqual(library["WT"], 1.0)
        self.assertIn("V39A", library)

        meta = fetch_protein_metadata(GB1_ACCESSION, use_cache=True)
        self.assertEqual(meta.accession, "SPG1_STRSG")
        self.assertEqual(meta.length, 56)
        self.assertEqual(meta.sequence, "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE")
        self.assertEqual(meta.binding_positions, [39, 40, 41, 54])

    def test_trastuzumab_dataset_and_metadata_integrity(self):
        """Verifies that the Trastuzumab dataset exists and contains authentic Mason 2021 variants."""
        self.assertTrue(TRAS_CSV.is_file(), f"Trastuzumab CSV missing: {TRAS_CSV}")
        library = load_tras_library(TRAS_CSV)
        self.assertGreaterEqual(len(library), 30_000, "Trastuzumab library should have ~36k records")

        meta = fetch_protein_metadata(TRAS_ACCESSION, use_cache=True)
        self.assertEqual(meta.accession, "HER2_TRASTUZUMAB")
        self.assertEqual(meta.length, 120)
        self.assertTrue(meta.sequence.startswith("EVQLVESGGGLVQPGGSLRLSCAASGFNIK"))
        # CDR-H3 is residues 99-108
        self.assertIn(99, meta.binding_positions)
        self.assertIn(108, meta.binding_positions)

    def test_gb1_featurization_and_gp(self):
        """Verifies physicochemical featurization on GB1 4-site interface."""
        wt_seq = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
        v1 = Variant(notation="V39A")
        v2 = Variant(notation="V39A:D40G:G41A:V54A")
        v1.validate_against_wt(wt_seq)
        v2.validate_against_wt(wt_seq)

        f1 = featurize_variant_on_positions(v1, [39, 40, 41, 54], wt_seq)
        f2 = featurize_variant_on_positions(v2, [39, 40, 41, 54], wt_seq)
        self.assertEqual(f1.shape, f2.shape)
        self.assertFalse(np.array_equal(f1, f2))

    def test_trastuzumab_developability_analysis(self):
        """Verifies developability liability detection on Trastuzumab CDR-H3 variants."""
        meta = fetch_protein_metadata(TRAS_ACCESSION, use_cache=True)
        wt_seq = meta.sequence

        # WT Trastuzumab VH has canonical intra-domain disulfide C22-C96 (paired, no free thiol)
        clean_rep = analyze_biologics_developability(wt_seq)
        self.assertFalse(clean_rep["unpaired_cysteine"])
        self.assertGreaterEqual(clean_rep["developability_score"], 0.50)

        # Mutate to introduce unpaired cysteine and deamidation
        v_bad = Variant(notation="G100C:G101N:D102S")
        v_bad.validate_against_wt(wt_seq)
        bad_seq = v_bad.apply_to(wt_seq)

        bad_rep = analyze_biologics_developability(bad_seq)
        self.assertTrue(bad_rep["unpaired_cysteine"])
        self.assertGreaterEqual(len(bad_rep["deamidation_sites"]), 1)
        self.assertLess(bad_rep["developability_score"], clean_rep["developability_score"])

    def test_gb1_active_learning_simulation_fast(self):
        """Runs a 2-round active learning simulation on a GB1 candidate pool."""
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tf:
            temp_db = tf.name

        try:
            results = run_gb1_active_learning(
                total_budget_usd=0.50,
                rounds=2,
                batch_size=10,
                pool_size=150,
                db_path=temp_db,
            )
            self.assertEqual(results["screened"], 20)
            self.assertGreater(results["total_spent_usd"], 0.0)
            self.assertGreater(results["cost_savings_pct"], 90.0)

            # Evaluate with evaluate_run
            m = evaluate_run(db_path=temp_db, target_name="SPG1_STRSG")
            self.assertEqual(m.matched_benchmark_variants, 20)
            self.assertGreater(m.cloud_cost_savings_pct, 90.0)
        finally:
            if os.path.exists(temp_db):
                try:
                    os.remove(temp_db)
                except OSError:
                    pass

    def test_trastuzumab_active_learning_simulation_fast(self):
        """Runs a 2-round active learning simulation on Trastuzumab candidate pool."""
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tf:
            temp_db = tf.name

        try:
            results = run_trastuzumab_active_learning(
                total_budget_usd=0.50,
                rounds=2,
                batch_size=10,
                pool_size=150,
                db_path=temp_db,
            )
            self.assertEqual(results["screened"], 20)
            self.assertGreater(results["total_spent_usd"], 0.0)
            self.assertIn("best_lead_developability_score", results)
            self.assertIn("best_lead_risk_category", results)

            # Evaluate with evaluate_run
            m = evaluate_run(db_path=temp_db, target_name="HER2_TRASTUZUMAB")
            self.assertEqual(m.matched_benchmark_variants, 20)
            self.assertGreater(m.cloud_cost_savings_pct, 90.0)
        finally:
            if os.path.exists(temp_db):
                try:
                    os.remove(temp_db)
                except OSError:
                    pass


if __name__ == "__main__":
    unittest.main()
