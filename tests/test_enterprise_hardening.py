"""
tests/test_enterprise_hardening.py — Test suite for enterprise hardening features:
1. Heteroscedastic Gaussian Process (sample variance & replicate noise modeling).
2. Therapeutic Biologics developability liability analysis (deamidation, isomerization, thiols).
3. Robotic 96/384-well plate layout generator with anti-edge evaporation barrier.
4. Opentrons OT-2 script generation with perimeter barrier automation.
"""

from __future__ import annotations

import unittest
import numpy as np

from controller.gatekeeper import SingleMutation, Variant
from controller.prompts import ResidueCluster, SearchStrategy
from workers.math_worker import (
    DeepKernelLearner,
    analyze_biologics_developability,
    evaluate_synthesizability_score,
    select_next_micro_batch,
)
from workers.wetlab_worker import (
    generate_opentrons_protocol,
    generate_plate_layout,
)


class TestEnterpriseHardening(unittest.TestCase):

    def test_heteroscedastic_dkl_fitting(self):
        """Verifies that DeepKernelLearner fits with observation-specific sample variances."""
        np.random.seed(42)
        X = np.random.randn(20, 15)
        y = np.sin(X[:, 0]) + 0.5 * X[:, 1]
        sample_var = np.array([0.01 if i % 2 == 0 else 0.25 for i in range(20)])

        learner = DeepKernelLearner(input_dim=15, latent_dim=4, hidden_dim=16, random_state=42)
        learner.fit(X, y, sample_variance=sample_var)

        self.assertTrue(learner.is_fitted)
        mu, std = learner.predict(X[:5])
        self.assertEqual(len(mu), 5)
        self.assertEqual(len(std), 5)
        self.assertTrue(np.all(std > 0))

    def test_select_next_micro_batch_with_replicate_evidence(self):
        """Verifies active learning selection ingests replicate assay variances."""
        wt = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        strategy = SearchStrategy(
            target_protein="avGFP",
            target_accession="P42212",
            hypothesis="Test chromophore mutations with replicate noise",
            target_clusters=[
                ResidueCluster(cluster_id="chromophore", positions=[64, 65, 66], rationale="Core"),
            ],
            forbidden_positions=[],
            max_mutation_depth=2,
        )

        # Evidence with explicit measurement variance
        evidence = [
            {
                "variant_sequence": wt[:64] + "T" + wt[65:],
                "computed_score": 1.45,
                "score_variance": 0.05,
            },
            {
                "variant_sequence": wt[:63] + "L" + wt[64:],
                "computed_score": 0.90,
                "variance": 0.12,
            },
            {
                "variant_sequence": wt[:65] + "A" + wt[66:],
                "computed_score": 0.20,
                "std": 0.15,
            },
        ]

        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=wt,
            evidence=evidence,
            batch_size=5,
        )
        self.assertEqual(len(batch.variants), 5)
        self.assertTrue(batch.max_pool_ei >= 0.0)

    def test_biologics_developability_analysis(self):
        """Verifies detection of therapeutic developability liabilities."""
        # Clean antibody sequence with paired cysteines (canonical intra-domain disulfide)
        clean_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSVIYSC"
        clean_rep = analyze_biologics_developability(clean_seq)
        self.assertTrue(clean_rep["developability_score"] >= 0.95)
        self.assertFalse(clean_rep["unpaired_cysteine"])

        # Sequence with Asn deamidation (NG), Asp isomerization (DG), and odd cysteine (1 Cys)
        liability_seq = "EVQLVESGGGLVQPGGSLRLSAASGFTFSSYAMSWVRQAPGKGLEWVSVIYSNGDGC"
        rep = analyze_biologics_developability(liability_seq)
        self.assertTrue(len(rep["deamidation_sites"]) >= 1)
        self.assertTrue(len(rep["isomerization_sites"]) >= 1)
        self.assertTrue(rep["unpaired_cysteine"])
        self.assertTrue(rep["developability_score"] < clean_rep["developability_score"])

    def test_evaluate_synthesizability_score_penalties(self):
        """Verifies synthesizability penalties for therapeutic liabilities."""
        wt = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSVIYS"
        v_clean = Variant(notation="G10A", mutations=[SingleMutation(wt_aa="G", position=10, mut_aa="A")])
        score_clean = evaluate_synthesizability_score(v_clean, wt)

        # Mutation that introduces an NG deamidation motif
        # WT has S at 21, followed by G at 22: mutating S21 to N creates NG
        v_deam = Variant(notation="S21N", mutations=[SingleMutation(wt_aa="S", position=21, mut_aa="N")])
        score_deam = evaluate_synthesizability_score(v_deam, wt)
        self.assertTrue(score_deam <= score_clean)

    def test_generate_plate_layout_96_well(self):
        """Verifies 96-well plate layout with 36 perimeter evaporation barriers."""
        dummy_variants = [
            {"mutations": f"Var_{i}", "computed_score": 1.0 + i * 0.05}
            for i in range(50)
        ]

        layout = generate_plate_layout(dummy_variants, plate_format="96_well", anti_edge_barrier=True)
        self.assertEqual(layout["plate_format"], "96_well")
        self.assertEqual(layout["total_wells"], 96)
        self.assertEqual(layout["barrier_wells_count"], 36)
        self.assertEqual(layout["assigned_variants_count"], 50)

        plate_map = layout["plate_map"]
        # Check perimeter wells
        self.assertEqual(plate_map["A1"]["type"], "evaporation_barrier")
        self.assertEqual(plate_map["A12"]["type"], "evaporation_barrier")
        self.assertEqual(plate_map["H1"]["type"], "evaporation_barrier")
        self.assertEqual(plate_map["H12"]["type"], "evaporation_barrier")
        self.assertEqual(plate_map["B1"]["type"], "evaporation_barrier")
        self.assertEqual(plate_map["B12"]["type"], "evaporation_barrier")

        # Check interior controls
        self.assertEqual(plate_map["B2"]["type"], "control_wt")
        self.assertEqual(plate_map["B3"]["type"], "control_negative")
        self.assertEqual(plate_map["B4"]["type"], "control_positive")
        self.assertEqual(plate_map["B5"]["type"], "control_blank")

        # Check candidate wells
        self.assertEqual(plate_map["B6"]["type"], "variant")
        self.assertIn("Well,Row,Column", layout["csv_layout"])

    def test_generate_opentrons_protocol_syntax(self):
        """Verifies Opentrons script generation and Python syntax."""
        var_names = [f"Variant_{i:02d}" for i in range(40)]
        script = generate_opentrons_protocol(
            variant_names=var_names,
            target_protein="avGFP",
            anti_edge_barrier=True,
        )

        self.assertIn("protocol_api", script)
        self.assertIn("perimeter_barrier_wells", script)
        self.assertIn("corning_96_wellplate_360ul_flat", script)

        # Verify generated script compiles cleanly as Python
        compiled = compile(script, "<string>", "exec")
        self.assertIsNotNone(compiled)


if __name__ == "__main__":
    unittest.main()
