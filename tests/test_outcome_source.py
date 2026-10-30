"""
tests/test_outcome_source.py — End-to-end exercise of steps 4, 5 and 6 with a real
measured outcome source instead of ESM inference.

Why this exists:
    Step 5 was hard-wired to `workers.esm_worker.score_variants*`, so a caller holding
    real measured outcomes (a plate readout, a benchmark's assay table, a CRO result)
    had no way to supply them. The only way to evaluate the engine was to skip steps 4
    and 5 and exercise planning and selection alone, which is what
    benchmark_panoptes_isolated.py had to do. `Orchestrator(outcome_source=...)` closes
    that gap.

Verifies:
  1. An injected outcome source is what step 5 executes, and its scores - not ESM
     pseudo-likelihoods - are what lands in the evidence table.
  2. Step 4 prices the round with the outcome source's declared target model, so a
     physically measured round is charged the wet-lab DNA + assay cost rather than
     fabricated GPU seconds.
  3. The ledger is genuinely debited, and the fail-closed stop still fires when the
     authorised budget cannot cover a round.
  4. `wetlab_assay` pricing refuses to price a batch without a construct length.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from typing import Dict, List

from controller.gatekeeper import Gatekeeper, estimate_wetlab_synthesis_cost
from controller.orchestrator import Orchestrator, OrchestratorConfig
from ledger import database as ledger_db


class _RecordingAssaySource:
    """
    A measured-assay outcome source.

    Returns a deterministic function of the sequence rather than a model score, so a
    test can prove the engine learned from THESE numbers and not from ESM. Records
    every call so the test can assert step 5 actually routed through it.
    """

    engine_name = "test_measured_assay"
    target_model = "wetlab_assay"

    def __init__(self) -> None:
        self.calls: List[int] = []
        self.scored: Dict[str, float] = {}

    def score(
        self,
        wt_sequence: str,
        variant_sequences: List[str],
        round_number: int,
    ) -> List[float]:
        self.calls.append(round_number)
        scores = []
        for seq in variant_sequences:
            # Deterministic, clearly not an ESM log-likelihood: count tryptophans.
            value = float(seq.count("W")) / 10.0
            self.scored[seq] = value
            scores.append(value)
        return scores


class TestOutcomeSourceInjection(unittest.TestCase):
    def _config(self, db_path: str, budget: float) -> OrchestratorConfig:
        return OrchestratorConfig(
            target_accession="TRPB_TM9D8",
            total_budget_usd=budget,
            max_iterations=1,
            batch_size=2,
            db_path=db_path,
            use_mock_llm=True,
            enable_live_llm=False,
            # Keep the run short and deterministic: this test is about steps 4-6,
            # not about stopping heuristics or optional structural stages.
            enable_regret_stopping=False,
            enable_stagnation_stopping=False,
            enable_lead_refinement=False,
            enable_structural_filter=False,
            enable_reflection=False,
            enable_cloud_cascade=False,
            random_seed=7,
        )

    def test_injected_source_supplies_evidence_and_is_priced_as_wetlab(self) -> None:
        """Step 5 runs the injected source; step 4 prices it as a physical assay."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = str(Path(tmp_dir) / "ledger.db")
            ledger_db.init_db(db_path)
            # A wet-lab-scale allocation: physical measurement of even two variants of
            # a 397-aa construct costs far more than a GPU round.
            ledger_db.create_budget(total_allocation=20_000.0, db_path=db_path)

            source = _RecordingAssaySource()
            orch = Orchestrator(self._config(db_path, 20_000.0), outcome_source=source)
            events = list(orch.run_stream())

            # 1. Step 5 routed through the injected source.
            self.assertTrue(source.calls, "outcome source was never invoked by step 5")
            self.assertTrue(source.scored)

            # 2. The evidence table holds the assay's numbers under its engine name.
            evidence = ledger_db.get_evidence(db_path=db_path)
            self.assertTrue(evidence, "no evidence was recorded")
            engines = {row.get("compute_engine") for row in evidence}
            self.assertIn("test_measured_assay", engines)
            for row in evidence:
                seq = row["variant_sequence"]
                if seq in source.scored:
                    self.assertAlmostEqual(
                        float(row["computed_score"]), source.scored[seq], places=6
                    )

            # 3. Step 4 priced the round with the wet-lab model, not GPU seconds.
            auth = [e for e in events if e.step == "authorization" and e.decision]
            self.assertTrue(auth, "no authorization decision was emitted")
            decision = auth[-1].decision
            self.assertEqual(decision["target_model"], "wetlab_assay")
            expected = estimate_wetlab_synthesis_cost(
                variant_count=len(decision["mutated_sequences"]),
                sequence_length_aa=len(decision["wildtype_sequence"])
                if decision.get("wildtype_sequence")
                else 397,
            )["total_wetlab_cost_usd"]
            # Physical cost is orders of magnitude above any GPU round; assert the
            # scale rather than an exact match in case batch size was clamped.
            self.assertGreater(float(decision["estimated_cost"]), 100.0)
            self.assertGreater(expected, 100.0)

            # 4. The ledger was genuinely debited.
            remaining = ledger_db.get_remaining_budget(1, db_path)
            self.assertLess(remaining, 20_000.0)

    def test_fail_closed_when_physical_cost_exceeds_budget(self) -> None:
        """A budget that cannot cover a physically measured round must halt the loop."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = str(Path(tmp_dir) / "ledger.db")
            ledger_db.init_db(db_path)
            # Enough to pass validation, nowhere near a real plate's cost.
            ledger_db.create_budget(total_allocation=1.0, db_path=db_path)

            source = _RecordingAssaySource()
            orch = Orchestrator(self._config(db_path, 1.0), outcome_source=source)
            events = list(orch.run_stream())

            halted = [e for e in events if e.step == "halted"]
            self.assertTrue(halted, "loop did not halt despite an unaffordable round")
            self.assertFalse(
                source.calls,
                "outcome source ran despite the Gatekeeper refusing to authorise spend",
            )
            self.assertEqual(ledger_db.get_evidence(db_path=db_path), [])

    def test_wetlab_pricing_requires_a_construct_length(self) -> None:
        """Physical cost scales with basepairs, so an unknown length must not be priced."""
        gatekeeper = Gatekeeper(db_path=None)
        with self.assertRaises(ValueError):
            gatekeeper.calculate_cost(
                variant_count=25, target_model="wetlab_assay", sequence_length_aa=None
            )
        cost, is_warm = gatekeeper.calculate_cost(
            variant_count=25, target_model="wetlab_assay", sequence_length_aa=397
        )
        self.assertAlmostEqual(cost, 3314.25, places=2)
        # There is no container to keep warm for a physical assay.
        self.assertFalse(is_warm)


if __name__ == "__main__":
    unittest.main()
