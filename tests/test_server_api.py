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
import unittest

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

    def test_campaign_sse_stream_endpoint(self) -> None:
        """Verifies real-time Server-Sent Events (SSE) streaming of active learning events."""
        with self.client.stream(
            "GET",
            "/api/campaign/stream?target_accession=P42212&total_budget_usd=0.50&max_iterations=1&batch_size=5",
        ) as resp:
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


if __name__ == "__main__":
    unittest.main()
