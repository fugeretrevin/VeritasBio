"""
Wet-lab spend must be priced by the server, never by the caller.

Both /api/cloudlab/dispatch and /api/campaign/ingest-curves used to pass
estimated_cost=request.nominal_cost_usd straight into ledger authorize_spend, so a
caller set the price its own spend was charged at. Sending 0.000001 authorised
effectively unlimited physical dispatches against the budget: the cap was denominated
in numbers the client supplied.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import ledger.database as ledger_db


class TestWetlabPricingIsServerSide(unittest.TestCase):
    def test_request_models_do_not_accept_a_cost(self) -> None:
        """A cost field on the request model is the vulnerability; it must be absent."""
        import server

        for model_name in ("CloudLabDispatchRequest", "IngestCurvesRequest"):
            model = getattr(server, model_name, None)
            if model is None:
                continue
            fields = set(getattr(model, "model_fields", {}))
            self.assertNotIn(
                "nominal_cost_usd", fields,
                f"{model_name} still lets the caller price its own spend",
            )

    def test_price_scales_with_variant_count_and_construct_length(self) -> None:
        import server

        # DNA synthesis plus robotic assay: cost scales with both how many variants are
        # built and how long each construct is. A flat per-round constant would miss both.
        small = server._price_wetlab_round(2, 238)
        more_variants = server._price_wetlab_round(25, 238)
        longer_protein = server._price_wetlab_round(2, 397)

        self.assertGreater(more_variants, small)
        self.assertGreater(longer_protein, small)

        # And the figure is wet-lab scale, not GPU scale. Pricing a physically measured
        # round as GPU seconds understated it by about five orders of magnitude.
        self.assertGreater(server._price_wetlab_round(25, 397), 100.0)

    def test_dispatch_refuses_when_the_budget_cannot_cover_the_real_price(self) -> None:
        from fastapi.testclient import TestClient

        import server as server_mod

        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "pricing.db"
            ledger_db.init_db(db_path)
            # Enough for the old fake price, nowhere near a real round.
            ledger_db.create_budget(total_allocation=5.00, db_path=db_path)
            saved = server_mod.DEFAULT_DB_PATH
            server_mod.DEFAULT_DB_PATH = db_path
            self.addCleanup(setattr, server_mod, "DEFAULT_DB_PATH", saved)

            client = TestClient(server_mod.app)
            resp = client.post(
                "/api/cloudlab/dispatch",
                json={
                    "provider": "emerald_cloud_lab",
                    "target_accession": "P42212",
                    "variants": ["S65T", "F64L:S65T"],
                    "fantasy_scores": [1.75, 2.20],
                    # Even if a caller tries the old field, pydantic drops it and the
                    # server prices the round itself.
                    "nominal_cost_usd": 0.000001,
                },
            )
            self.assertEqual(resp.status_code, 402, resp.text)
            # The refusal states the real price, so an operator can see what to fund.
            self.assertIn("prices at", resp.text)

            budget = ledger_db.get_budget(1, db_path=db_path)
            self.assertEqual(budget["cumulative_spent"], 0.0,
                             "a refused dispatch must not move the budget")


if __name__ == "__main__":
    unittest.main()
