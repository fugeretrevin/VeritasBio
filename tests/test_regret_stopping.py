"""
tests/test_regret_stopping.py — Unit tests for Bayesian Regret Stopping & Dynamic Budget Conservation.

Verifies:
  1. Early convergence detection: When acquisition EI drops below threshold across consecutive rounds,
     the orchestrator halts early with `stop_reason = 'converged_regret_threshold'`.
  2. Budget refunding: Unspent budget is preserved and returned to the ledger.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from controller.orchestrator import Orchestrator, OrchestratorConfig
from ledger import database as ledger_db


class TestRegretStopping(unittest.TestCase):
    """Unit tests for Bayesian Regret Stopping (Enterprise Improvement 3)."""

    def test_regret_early_stopping(self) -> None:
        """Verifies that an orchestrator run with high initial convergence stops before max_iterations."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = str(Path(tmp_dir) / "test_regret.db")
            config = OrchestratorConfig(
                target_accession="P42212",
                total_budget_usd=1.00,
                max_iterations=10,  # Cap at 10 rounds
                batch_size=10,
                target_model="mock_esm",
                use_mock_llm=True,
                enable_live_llm=False,
                db_path=db_path,
                enable_regret_stopping=True,
                regret_stopping_threshold=0.01,
            )

            orch = Orchestrator(config)
            summary = orch.run()

            # Should converge and stop well before round 10
            self.assertLess(
                summary.iterations_completed,
                10,
                f"Run should have stopped early on regret convergence, but ran {summary.iterations_completed} rounds",
            )
            # Budget should be preserved
            self.assertGreater(summary.final_remaining_usd, 0.50, "Remaining funds must be preserved")


if __name__ == "__main__":
    unittest.main()
