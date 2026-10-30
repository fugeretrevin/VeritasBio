"""
tests/test_firebase_auth.py — Unit tests for PanoptesBio Firebase Authentication & FDA 21 CFR Part 11 Electronic Identity.

Tests:
  1. FirebaseConfig default properties & client dictionary export.
  2. Mock and development token verification across RBAC roles.
  3. Fail-closed rejection of empty, whitespace, or invalid tokens.
  4. FDA 21 CFR §11.200 HMAC-SHA256 electronic signature generation.
  5. Gatekeeper & SQLite Ledger integration: verifying operator identity is stamped into immutable Merkle blocks.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from controller.auth import (
    AuthenticatedUser,
    AuthenticationError,
    create_part11_electronic_signature,
    verify_firebase_id_token,
)
from controller.firebase_config import FirebaseConfig, get_firebase_config
from controller.gatekeeper import BatchEvaluationRequest, Gatekeeper, SingleMutation, Variant
from ledger import database as ledger_db


class TestFirebaseAuth(unittest.TestCase):
    """Test suite for Firebase Authentication & FDA Part 11 Electronic Identity."""

    def test_firebase_config_defaults(self) -> None:
        """Verifies that PanoptesBio project configuration loads default credentials."""
        cfg = get_firebase_config()
        self.assertEqual(cfg.project_id, "panoptesbio")
        self.assertEqual(cfg.auth_domain, "panoptesbio.firebaseapp.com")
        self.assertEqual(cfg.storage_bucket, "panoptesbio.firebasestorage.app")
        self.assertEqual(cfg.measurement_id, "G-3CQM0RGVH3")

        client_dict = cfg.to_client_dict()
        self.assertIn("apiKey", client_dict)
        self.assertIn("projectId", client_dict)
        self.assertEqual(client_dict["projectId"], "panoptesbio")

    def test_verify_mock_tokens_and_roles(self) -> None:
        """Verifies that mock and test tokens correctly resolve AuthenticatedUser with roles."""
        user1 = verify_firebase_id_token("mock_token_director_john")
        self.assertEqual(user1.uid, "usr_john")
        self.assertEqual(user1.email, "john@panoptesbio.com")
        self.assertEqual(user1.role, "lab_director")
        self.assertEqual(user1.organization_id, "panoptesbio")
        # Mock identities are deliberately NOT Part 11 certified: their role is
        # parsed from the caller-supplied token string, so they must never be able
        # to mint an electronic signature.
        self.assertFalse(user1.is_part11_certified)

        user2 = verify_firebase_id_token("mock_token_compliance_alice")
        self.assertEqual(user2.uid, "usr_alice")
        self.assertEqual(user2.email, "alice@panoptesbio.com")
        self.assertEqual(user2.role, "compliance_officer")

        user_default = verify_firebase_id_token("panoptes_demo_token")
        self.assertEqual(user_default.uid, "usr_panoptes_001")
        self.assertEqual(user_default.role, "lead_biologist")

    def test_mock_tokens_rejected_when_not_explicitly_enabled(self) -> None:
        """Mock identities must be opt-in; the shipped default rejects them."""
        import os
        from controller.auth import mock_auth_enabled

        saved = os.environ.pop("VERITAS_ALLOW_MOCK_AUTH", None)
        try:
            self.assertFalse(mock_auth_enabled())
            # Role is parsed out of the token string, so accepting these by default
            # let anyone self-assign lab_director from a browser.
            for forged in ("panoptes_demo_token", "mock_token_director_attacker",
                           "mock_token_compliance_attacker"):
                with self.assertRaises(AuthenticationError):
                    verify_firebase_id_token(forged)
        finally:
            if saved is not None:
                os.environ["VERITAS_ALLOW_MOCK_AUTH"] = saved

    def test_unsigned_jwt_claims_are_not_trusted(self) -> None:
        """A JWT whose signature cannot be verified must be rejected, not decoded."""
        import base64, json

        payload = base64.urlsafe_b64encode(
            json.dumps({"sub": "usr_attacker", "email": "attacker@evil.test",
                        "role": "compliance_officer"}).encode()
        ).decode().rstrip("=")
        forged = f"eyJhbGciOiJub25lIn0.{payload}.signature_not_checked"

        # Previously this token's claims were base64-decoded and trusted verbatim,
        # including `role`, whenever live verification raised.
        with self.assertRaises(AuthenticationError):
            verify_firebase_id_token(forged, allow_mock=True)

    def test_signing_refuses_without_secret_or_certification(self) -> None:
        """Electronic signatures require a configured secret AND a certified identity."""
        import os
        import controller.auth as auth_mod

        user = verify_firebase_id_token("mock_token_director_john", allow_mock=True)
        self.assertFalse(user.is_part11_certified)
        # Mock identity -> refused even though a secret is configured in tests.
        with self.assertRaises(AuthenticationError):
            create_part11_electronic_signature(user, action="APPROVE", meaning="test")

        # Certified identity but no secret -> also refused, because a signature keyed
        # on an empty or committed secret is forgeable offline from the source.
        certified = user.model_copy(update={"is_part11_certified": True})
        saved = auth_mod.PART11_SIGNATURE_SECRET
        auth_mod.PART11_SIGNATURE_SECRET = ""
        try:
            with self.assertRaises(AuthenticationError):
                create_part11_electronic_signature(certified, action="APPROVE", meaning="test")
        finally:
            auth_mod.PART11_SIGNATURE_SECRET = saved

    def test_fail_closed_invalid_tokens(self) -> None:
        """Verifies fail-closed rejection on missing, empty, or forged tokens."""
        with self.assertRaises(AuthenticationError):
            verify_firebase_id_token(None)

        with self.assertRaises(AuthenticationError):
            verify_firebase_id_token("")

        with self.assertRaises(AuthenticationError):
            verify_firebase_id_token("   ")

        # Invalid token without allow_mock should raise AuthenticationError
        with self.assertRaises(AuthenticationError):
            verify_firebase_id_token("mock_token_fake", allow_mock=False)

    def test_part11_electronic_signature(self) -> None:
        """Verifies that electronic signatures satisfy FDA 21 CFR §11.200 dual-component binding."""
        user = AuthenticatedUser(
            uid="usr_director_01",
            email="director@panoptesbio.com",
            display_name="Dr. Director",
            organization_id="panoptesbio",
            role="lab_director",
        )
        sig = create_part11_electronic_signature(
            user=user,
            action="Approval of Modal A10G Compute Spend",
            meaning="Authorized compute budget spend under 21 CFR §11.200",
        )

        self.assertEqual(sig["operator_uid"], "usr_director_01")
        self.assertEqual(sig["operator_email"], "director@panoptesbio.com")
        self.assertEqual(sig["role"], "lab_director")
        self.assertIn("electronic_signature_hash", sig)
        self.assertEqual(len(sig["electronic_signature_hash"]), 64)  # SHA-256 hex length
        self.assertEqual(sig["signature_algorithm"], "HMAC-SHA256 (21 CFR §11.200)")

    def test_gatekeeper_stamping_operator_into_ledger(self) -> None:
        """Verifies that Gatekeeper verifies the auth token and stamps operator identity into SQLite Merkle blocks."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = str(Path(tmp_dir) / "test_auth_ledger.db")
            ledger_db.init_db(db_path)
            ledger_db.create_budget(total_allocation=10.0, db_path=db_path)

            gatekeeper = Gatekeeper(db_path=Path(db_path))
            wt_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
            variant = Variant(
                notation="S65T",
                mutations=[SingleMutation(wt_aa="S", position=65, mut_aa="T")],
            )

            req = BatchEvaluationRequest(
                job_id="job-auth-test-001",
                target_model="mock_esm",
                wildtype_sequence=wt_seq,
                variants=[variant],
                auth_token="mock_token_director_sarah",
            )

            decision = gatekeeper.evaluate_and_authorize(req)
            self.assertTrue(decision.approved, f"Job must be approved: {decision.rejection_reason}")

            # Inspect the cryptographic audit blocks in SQLite
            audit_blocks = ledger_db.get_audit_blocks(db_path=db_path)
            spend_blocks = [b for b in audit_blocks if b["action_type"] == "SPEND_AUTHORIZED"]
            self.assertEqual(len(spend_blocks), 1)

            payload = json.loads(spend_blocks[0]["payload_json"])
            self.assertEqual(payload["operator_uid"], "usr_sarah")
            self.assertEqual(payload["operator_email"], "sarah@panoptesbio.com")

            # Verify entire FDA Part 11 Merkle hash chain remains intact
            is_valid, count, msg = ledger_db.verify_audit_integrity(db_path=db_path)
            self.assertTrue(is_valid, f"Merkle chain must be intact: {msg}")
            self.assertGreaterEqual(count, 3)


if __name__ == "__main__":
    unittest.main()
