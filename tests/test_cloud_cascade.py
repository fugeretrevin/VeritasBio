"""
tests/test_cloud_cascade.py — Unit and integration tests for Multi-Fidelity Cloud Cascade & Lead Refinement.

Verifies:
  1. Fail-closed tiered pricing calculations in Gatekeeper for 'esm_cascade' and 'mock_cascade'.
  2. Mathematical guarantee of Hamming-1 neighborhood generation (e.g. discovering V183A:F184I:V227M:S228G from F184I:V227M:S228G).
  3. Non-greedy acquisition filtering combining Tier-1 scores with DKL uncertainty and CoT bonuses.
  4. Multi-fidelity scoring behavior on MockESMWorker.
  5. End-to-end active learning campaign with cloud cascade and Tier-3 lead polish.
"""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from controller.gatekeeper import Gatekeeper, GatekeeperConfig, SingleMutation, Variant
from controller.orchestrator import Orchestrator, OrchestratorConfig
from ledger import database as ledger_db
from workers.math_worker import (
    filter_cascade_promoted_candidates,
    generate_hamming1_lead_neighbors,
)
from workers.mock_workers import MockESMWorker
from workers.rag_worker import fetch_protein_metadata


class TestCloudCascadeAndLeadRefinement(unittest.TestCase):
    """Test suite for guaranteed-value cloud cascade and lead refinement."""

    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.db_path = Path(self.tmp_dir.name) / "test_cascade.db"
        ledger_db.init_db(self.db_path)
        ledger_db.create_budget(total_allocation=2.00, budget_id=1, db_path=self.db_path)

        self.gk_config = GatekeeperConfig(
            a10g_usd_per_sec=0.000306,
            gpu_sec_per_variant=0.02,
            gpu_sec_per_variant_150m=0.005,
            cascade_top_k_ratio=0.20,
            cold_start_sec=20.0,
            idle_timeout_sec=45.0,
            mock_nominal_cost=0.0001,
        )
        self.gatekeeper = Gatekeeper(config=self.gk_config, db_path=self.db_path)

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

    def test_tiered_cost_calculation(self) -> None:
        """Verify tiered cascade pricing is strictly cheaper than brute-force 650M."""
        n_variants = 100
        # Warm container cost
        cost_cascade, is_warm = self.gatekeeper.calculate_cost(
            variant_count=n_variants,
            target_model="esm_cascade",
            force_warm=True,
        )
        cost_650m, _ = self.gatekeeper.calculate_cost(
            variant_count=n_variants,
            target_model="esm2_650m",
            force_warm=True,
        )

        # Expected:
        # Cascade: (100 * 0.005) + (20 * 0.02) = 0.5s + 0.4s = 0.9s * 0.000306 = $0.000275
        # 650M: 100 * 0.02 = 2.0s * 0.000306 = $0.000612
        self.assertTrue(is_warm)
        self.assertLess(cost_cascade, cost_650m)
        ratio = cost_cascade / cost_650m
        self.assertAlmostEqual(ratio, 0.9 / 2.0, places=2)

        # Mock cascade cost verification
        mock_cost, _ = self.gatekeeper.calculate_cost(
            variant_count=n_variants,
            target_model="mock_cascade",
        )
        self.assertEqual(mock_cost, 0.00015)

    def test_hamming1_lead_neighbors_discovers_trpb_optimum(self) -> None:
        """
        Verify that given the 3-mutation lead F184I:V227M:S228G (fitness 0.8005),
        Hamming-1 Polish systematically generates the 1.0000 global optimum V183A:F184I:V227M:S228G!
        """
        meta = fetch_protein_metadata("TRPB_TM9D8", use_cache=True)
        wt_seq = meta.sequence

        lead = Variant(
            notation="F184I:V227M:S228G",
            mutations=[
                SingleMutation(wt_aa="F", position=184, mut_aa="I"),
                SingleMutation(wt_aa="V", position=227, mut_aa="M"),
                SingleMutation(wt_aa="S", position=228, mut_aa="G"),
            ],
        )
        lead.validate_against_wt(wt_seq)

        active_positions = [183, 184, 227, 228]
        neighbors = generate_hamming1_lead_neighbors(
            lead_variant=lead,
            active_positions=active_positions,
            wt_sequence=wt_seq,
        )

        notations = {n.notation for n in neighbors}

        # 1. Mathematical guarantee: 4 positions * 19 alternative AAs = 76 neighbors
        self.assertEqual(len(neighbors), 76)
        self.assertEqual(len(notations), 76)

        # 2. Mathematical guarantee: V183A:F184I:V227M:S228G is in the generated set!
        self.assertIn("V183A:F184I:V227M:S228G", notations)

        # 3. Every neighbor must be strictly 1 mutation away from lead
        for n in neighbors:
            n.validate_against_wt(wt_seq)
            # Compare mutations
            pos_diffs = 0
            lead_map = {184: "I", 227: "M", 228: "G"}
            for p in active_positions:
                expected_lead_aa = lead_map.get(p, wt_seq[p - 1])
                curr_neighbor_aa = wt_seq[p - 1]
                for m in n.mutations:
                    if m.position == p:
                        curr_neighbor_aa = m.mut_aa
                if curr_neighbor_aa != expected_lead_aa:
                    pos_diffs += 1
            self.assertEqual(pos_diffs, 1)

    def test_filter_cascade_promoted_candidates_non_greedy(self) -> None:
        """Verify non-greedy promotion protects high-uncertainty / high-CoT variants."""
        cands = [
            Variant(notation="S65T", mutations=[SingleMutation(wt_aa="S", position=65, mut_aa="T")]),
            Variant(notation="T203Y", mutations=[SingleMutation(wt_aa="T", position=203, mut_aa="Y")]),
            Variant(notation="V163A", mutations=[SingleMutation(wt_aa="V", position=163, mut_aa="A")]),
            Variant(notation="F64L", mutations=[SingleMutation(wt_aa="F", position=64, mut_aa="L")]),
        ]

        # Tier 1 raw scores: S65T is highest, F64L is lowest
        t1_scores = [1.5, 1.2, 0.9, 0.8]
        # But F64L has very high surrogate uncertainty (e.g. 1.8) and CoT bonus (0.36)
        dkl_stds = [0.1, 0.1, 0.2, 1.8]
        cot_bonuses = [0.0, 0.0, 0.0, 0.36]

        promoted, indices = filter_cascade_promoted_candidates(
            candidates=cands,
            tier1_scores=t1_scores,
            top_k_ratio=0.50,  # Promote top 2
            dkl_stds=dkl_stds,
            cot_bonuses=cot_bonuses,
            kappa=0.6,
        )

        self.assertEqual(len(promoted), 2)
        # F64L priority: 0.8 + 0.6 * 1.8 + 0.36 = 2.24 -> Ranked #1!
        self.assertEqual(promoted[0].notation, "F64L")
        self.assertEqual(indices[0], 3)

    def test_mock_cascade_worker_scoring(self) -> None:
        """Verify MockESMWorker.score_variants_cascade returns valid tiered results."""
        meta = fetch_protein_metadata("P42212", use_cache=True)
        wt_seq = meta.sequence

        worker = MockESMWorker(wt_sequence=wt_seq)
        variants = [
            Variant(notation="S65T", mutations=[SingleMutation(wt_aa="S", position=65, mut_aa="T")]),
            Variant(notation="T203Y", mutations=[SingleMutation(wt_aa="T", position=203, mut_aa="Y")]),
            Variant(notation="V163A", mutations=[SingleMutation(wt_aa="V", position=163, mut_aa="A")]),
            Variant(notation="F64L", mutations=[SingleMutation(wt_aa="F", position=64, mut_aa="L")]),
            Variant(notation="E222G", mutations=[SingleMutation(wt_aa="E", position=222, mut_aa="G")]),
        ]
        seqs = [v.apply_to(wt_seq) for v in variants]

        res = worker.score_variants_cascade(wt_seq, seqs, top_k_ratio=0.40)
        self.assertEqual(len(res["tier1_scores"]), 5)
        self.assertEqual(len(res["promoted_indices"]), 2)
        self.assertEqual(len(res["promoted_sequences"]), 2)
        self.assertEqual(len(res["tier2_scores"]), 2)

        for s in res["tier1_scores"]:
            self.assertTrue(math.isfinite(s))
        for s in res["tier2_scores"]:
            self.assertTrue(math.isfinite(s))

    def test_orchestrator_cascade_and_lead_polish_loop(self) -> None:
        """Integration test: runs active learning loop with cascade and lead polish enabled."""
        config = OrchestratorConfig(
            target_accession="P42212",
            total_budget_usd=1.00,
            max_iterations=2,
            batch_size=10,
            target_model="mock_esm",
            db_path=str(self.db_path),
            use_mock_llm=True,
            enable_live_llm=False,
            enable_cloud_cascade=True,
            cascade_top_k_ratio=0.20,
            enable_lead_refinement=True,
        )

        orchestrator = Orchestrator(config)
        events = list(orchestrator.run_stream())

        event_steps = [e.step for e in events]
        self.assertIn("init", event_steps)
        self.assertIn("authorization", event_steps)
        self.assertIn("execution", event_steps)
        self.assertIn("update", event_steps)
        self.assertIn("refinement", event_steps)
        self.assertIn("complete", event_steps)

        # Check ledger recorded jobs and evidence
        jobs = ledger_db.list_jobs(db_path=self.db_path)
        self.assertGreaterEqual(len(jobs), 2)
        # Check that mock_cascade was authorized for the rounds
        round_jobs = [j for j in jobs if "refine" not in j["job_id"]]
        for rj in round_jobs:
            self.assertEqual(rj["target_model"], "mock_cascade")

        evidence = ledger_db.get_evidence(db_path=self.db_path)
        self.assertGreaterEqual(len(evidence), 20)

        # Final event checks
        last_event = events[-1]
        self.assertEqual(last_event.step, "complete")
        self.assertGreater(last_event.top_variants[0]["computed_score"], 0.0)


if __name__ == "__main__":
    unittest.main()
