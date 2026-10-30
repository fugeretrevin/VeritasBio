"""
tests/test_server_api.py — Unit Tests for PanoptesBio / VeritasBio FastAPI Backend.

Verifies:
    1. /api/health returns online status and budget allocation telemetry.
    2. /api/targets lists validated benchmarks (TrpB, GB1, HER2, avGFP, Scorpion Antivenom).
    3. /api/structure/{id} resolves authentic PDB coordinates.
    4. /api/campaign/stream streams Server-Sent Events (SSE) in real time to web clients.
    5. /api/ledger returns Merkle block audit integrity.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from server import app


class TestServerAPI(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.client = TestClient(app)

    def test_health_endpoint(self) -> None:
        """Verifies engine health and budget status."""
        resp = self.client.get("/api/health")
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertEqual(data.get("status"), "online")
        self.assertIn("budget", data)
        self.assertIn("remaining_usd", data["budget"])

    def test_targets_list_endpoint(self) -> None:
        """Verifies that all 5 target proteins and benchmarks are registered."""
        resp = self.client.get("/api/targets")
        self.assertEqual(resp.status_code, 200)
        targets = resp.json()
        self.assertGreaterEqual(len(targets), 5)
        target_ids = [t["id"] for t in targets]
        self.assertIn("TRPB_TM9D8", target_ids)
        self.assertIn("SCFV_9004G_CN2", target_ids)
        self.assertIn("HER2_TRASTUZUMAB", target_ids)
        self.assertIn("SPG1_STRSG", target_ids)
        self.assertIn("P42212", target_ids)

    def test_structure_pdb_endpoint(self) -> None:
        """Verifies PDB coordinate delivery for 3Dmol viewer."""
        resp = self.client.get("/api/structure/TRPB_TM9D8")
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertEqual(data.get("pdb_filename"), "8EH0.pdb")
        self.assertIn("ATOM", data.get("pdb_content", ""))

    def test_campaign_stream_requires_post_authorization(self) -> None:
        """A GET alone must not be able to start a campaign that spends money."""
        # No ticket: the old endpoint took budget and batch size straight off the
        # query string, so an <img> tag could start a paid run.
        resp = self.client.get("/api/campaign/stream/not-a-real-ticket")
        self.assertEqual(resp.status_code, 404)

        prepared = self.client.post(
            "/api/campaign/prepare-stream",
            json={"target_accession": "P42212", "total_budget_usd": 0.50,
                  "max_iterations": 1, "batch_size": 5},
        )
        self.assertEqual(prepared.status_code, 200, prepared.text)
        ticket = prepared.json()["stream_ticket"]

        # Tickets are single-use: consume it, then confirm a replay is refused.
        with self.client.stream("GET", f"/api/campaign/stream/{ticket}") as resp:
            self.assertEqual(resp.status_code, 200)
            for line in resp.iter_lines():
                if line.startswith("data:"):
                    if json.loads(line[5:].strip()).get("step") == "complete":
                        break
        replay = self.client.get(f"/api/campaign/stream/{ticket}")
        self.assertEqual(replay.status_code, 404)

    def test_campaign_sse_stream_endpoint(self) -> None:
        """Verifies real-time Server-Sent Events (SSE) streaming of active learning events."""
        prepared = self.client.post(
            "/api/campaign/prepare-stream",
            json={"target_accession": "P42212", "total_budget_usd": 0.50,
                  "max_iterations": 1, "batch_size": 5},
        )
        self.assertEqual(prepared.status_code, 200, prepared.text)
        ticket = prepared.json()["stream_ticket"]

        with self.client.stream("GET", f"/api/campaign/stream/{ticket}") as resp:
            self.assertEqual(resp.status_code, 200)
            self.assertEqual(resp.headers.get("content-type"), "text/event-stream; charset=utf-8")
            
            data_events = []
            for line in resp.iter_lines():
                if line.startswith("data:"):
                    payload = json.loads(line[5:].strip())
                    data_events.append(payload)
                    if payload.get("step") == "complete":
                        break
            
            self.assertGreater(len(data_events), 0)
            steps = [e.get("step") for e in data_events]
            self.assertIn("init", steps)

    def test_ledger_state_endpoint(self) -> None:
        """Verifies ledger state and cryptographic Merkle chain verification."""
        resp = self.client.get("/api/ledger")
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertIn("budget", data)
        self.assertIn("audit_integrity", data)
        self.assertTrue(data["audit_integrity"]["is_valid"])

    def test_ingest_plate_json_endpoint(self) -> None:
        """Verifies microplate reader CSV ingestion via JSON body into SQLite ledger."""
        csv_payload = (
            "Well,Sample_Type,Sample_Name,Reading\n"
            "A01,evaporation_barrier,H2O,0.01\n"
            "A02,evaporation_barrier,H2O,0.01\n"
            "B02,control_wt,WT,1.00\n"
            "B03,variant,S65T:V68L,1.85\n"
            "B04,variant,T203Y,1.40\n"
        )
        resp = self.client.post(
            "/api/campaign/ingest-plate",
            json={
                "csv_content": csv_payload,
                "target_accession": "P42212",
                "normalize_to_wt": True,
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data.get("success"))
        self.assertEqual(data.get("variants_ingested"), 3)
        self.assertEqual(data.get("barrier_wells_ignored"), 2)
        self.assertAlmostEqual(data.get("wt_baseline_score"), 1.00)
        self.assertAlmostEqual(data.get("top_fitness_discovered"), 1.85)
        self.assertIn("block_hash", data)

    def test_ingest_plate_multipart_endpoint(self) -> None:
        """Verifies microplate reader CSV ingestion via multipart file upload."""
        csv_payload = (
            "Well,Sample_Type,Sample_Name,Reading\n"
            "A01,evaporation_barrier,H2O,0.02\n"
            "B02,control_wt,WT,1.00\n"
            "B05,variant,F64L:S65T,2.20\n"
        )
        files = {"file": ("tecan_plate_data.csv", csv_payload.encode("utf-8"), "text/csv")}
        data = {"target_accession": "P42212", "normalize_to_wt": "true"}
        resp = self.client.post("/api/campaign/ingest-plate", files=files, data=data)
        self.assertEqual(resp.status_code, 200)
        res_data = resp.json()
        self.assertTrue(res_data.get("success"))
        self.assertEqual(res_data.get("variants_ingested"), 2)
        self.assertEqual(res_data.get("barrier_wells_ignored"), 1)
        self.assertAlmostEqual(res_data.get("top_fitness_discovered"), 2.20)

    def test_wetlab_export_endpoint(self) -> None:
        """Verifies physical biology and LIMS export files generation."""
        resp = self.client.post("/api/wetlab/export?target_name=Test_Leads&plate_format=96_well")
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertIn("twist_csv", data)
        self.assertIn("opentrons_script", data)
        self.assertIn("fda_certificate_json", data)
        self.assertIn("benchling_manifest", data)
        self.assertIn("plate_csv", data)
        self.assertIn("Test_Leads", data["twist_csv"])
        self.assertIn("opentrons", data["opentrons_script"].lower())
        self.assertIn("Plate_Name", data["benchling_manifest"])

    def test_pareto_frontier_endpoint(self) -> None:
        """Verifies multi-objective Pareto trade-off analysis API (Track B)."""
        payload = {
            "target_protein": "SCFV_9004G_CN2",
            "variants": ["G24V", "G24V:S25N"],
            "weights": {
                "fitness": 0.40,
                "developability": 0.30,
                "synthesizability": 0.15,
                "stability": 0.15,
            },
        }
        resp = self.client.post("/api/campaign/pareto-frontier", json=payload)
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertEqual(data.get("target_protein"), "SCFV_9004G_CN2")
        self.assertIn("points", data)
        self.assertIn("frontier_points", data)
        self.assertIn("chart_data", data)
        self.assertGreaterEqual(len(data["points"]), 2)
        self.assertGreaterEqual(data["rank_0_count"], 1)
        self.assertTrue(any(pt["is_frontier"] for pt in data["chart_data"]))




    def test_relax_structure_endpoint(self) -> None:
        """Verifies Track C sidechain rotamer relaxation endpoint."""
        payload = {
            "target_protein": "P42212",
            "mutations": ["S65T"],
            "max_iterations": 100,
            "harmonic_k": 50.0,
        }
        resp = self.client.post("/api/campaign/relax-structure", json=payload)
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data.get("success"))
        self.assertEqual(data.get("pdb_filename"), "1EMA.pdb")
        res = data.get("result", {})
        self.assertEqual(res.get("variant_notation"), "S65T")
        self.assertTrue(res.get("converged"))
        self.assertIn("ATOM", res.get("relaxed_pdb_str", ""))

    def test_auth_register_and_login_endpoints(self) -> None:
        """
        Verifies operator account registration and sign-in via /api/auth/register and
        /api/auth/login.

        The account store is redirected to a temporary file: these routes persist real
        credential material (emails, PBKDF2 hashes and salts), so the suite must not
        write into ledger/users.json. Isolation also keeps the test idempotent, since
        re-registering an existing email is now correctly refused.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            with patch("controller.auth.USERS_STORE_PATH", Path(tmp_dir) / "users.json"):
                reg_resp = self.client.post(
                    "/api/auth/register",
                    json={
                        "email": "test.operator@panoptesbio.org",
                        "password": "securepassword123",
                        "display_name": "Dr. Test Operator",
                        "organization_id": "PanoptesBio Lab",
                    },
                )
                self.assertEqual(reg_resp.status_code, 200)
                reg_data = reg_resp.json()
                self.assertTrue(reg_data.get("authenticated"))
                self.assertIn("id_token", reg_data)
                self.assertEqual(reg_data["user"]["email"], "test.operator@panoptesbio.org")
                # Self-service registration mints the lowest-privilege role and cannot sign.
                self.assertEqual(reg_data["user"]["role"], "lead_biologist")
                self.assertFalse(reg_data["user"]["is_part11_certified"])

                login_resp = self.client.post(
                    "/api/auth/login",
                    json={
                        "email": "test.operator@panoptesbio.org",
                        "password": "securepassword123",
                    },
                )
                self.assertEqual(login_resp.status_code, 200)
                login_data = login_resp.json()
                self.assertTrue(login_data.get("authenticated"))
                self.assertEqual(login_data["user"]["display_name"], "Dr. Test Operator")

    def test_auth_register_ignores_client_supplied_role(self) -> None:
        """
        A privileged role in the request body must not reach the account record.
        `role` is no longer part of RegisterUserRequest, so pydantic drops it and the
        account is minted as lead_biologist rather than 400-ing on the refusal path.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            with patch("controller.auth.USERS_STORE_PATH", Path(tmp_dir) / "users.json"):
                resp = self.client.post(
                    "/api/auth/register",
                    json={
                        "email": "climber@panoptesbio.org",
                        "password": "securepassword123",
                        "display_name": "Role Climber",
                        "role": "lab_director",
                    },
                )
                self.assertEqual(resp.status_code, 200)
                user = resp.json()["user"]
                self.assertEqual(user["role"], "lead_biologist")
                self.assertFalse(user["is_part11_certified"])


if __name__ == "__main__":
    unittest.main()

