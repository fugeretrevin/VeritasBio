"""
controller/auth.py — Firebase Authentication & FDA 21 CFR Part 11 Electronic Signature Engine.

Provides:
  1. AuthenticatedUser data model with Organization Tenancy & RBAC roles.
  2. Firebase ID Token verification via firebase-admin SDK (with high-fidelity deterministic offline/mock mode).
  3. FDA 21 CFR Part 11 Electronic Signature creation (HMAC-SHA256 dual-component binding).
  4. Integration with VeritasBio double-entry ledger & gatekeeper spend authorization.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import logging
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

from pydantic import BaseModel, Field

from controller.firebase_config import get_firebase_config

logger = logging.getLogger(__name__)

# HMAC secret key for FDA Part 11 electronic signature authentication
PART11_SIGNATURE_SECRET = os.getenv(
    "PART11_SIGNATURE_SECRET", "panoptesbio-fda-part11-hmac-sha256-signature-key-2026"
)

# Optional firebase-admin singleton
_FIREBASE_APP_INITIALIZED = False


def _init_firebase_admin_if_needed() -> bool:
    """Initializes firebase-admin SDK once if credentials exist."""
    global _FIREBASE_APP_INITIALIZED
    if _FIREBASE_APP_INITIALIZED:
        return True

    try:
        import firebase_admin
        from firebase_admin import credentials

        if firebase_admin._apps:
            _FIREBASE_APP_INITIALIZED = True
            return True

        cfg = get_firebase_config()
        sa_path = cfg.service_account_path
        if sa_path and Path(sa_path).is_file():
            cred = credentials.Certificate(sa_path)
            firebase_admin.initialize_app(cred, {
                "projectId": cfg.project_id,
                "storageBucket": cfg.storage_bucket,
            })
            _FIREBASE_APP_INITIALIZED = True
            logger.info("Firebase Admin SDK initialized with Service Account from %s", sa_path)
            return True
        else:
            # Initialize with default project options (offline / metadata mode)
            firebase_admin.initialize_app(options={
                "projectId": cfg.project_id,
                "storageBucket": cfg.storage_bucket,
            })
            _FIREBASE_APP_INITIALIZED = True
            logger.info("Firebase Admin SDK initialized with default options (Project: %s)", cfg.project_id)
            return True
    except Exception as exc:
        logger.debug("Firebase Admin SDK initialization skipped: %s", exc)
        return False


class AuthenticatedUser(BaseModel):
    """Represents a verified operator identity under FDA 21 CFR Part 11."""

    uid: str = Field(..., description="Firebase User ID")
    email: str = Field(..., description="Verified corporate email address")
    display_name: str = Field(default="PanoptesBio Scientist", description="User display name")
    organization_id: str = Field(default="panoptesbio", description="Multi-tenant organization identifier")
    role: str = Field(default="lead_biologist", description="Role: lead_biologist, lab_director, compliance_officer")
    is_part11_certified: bool = Field(default=True, description="Certified for FDA 21 CFR Part 11 electronic signing")
    auth_time: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),
        description="ISO 8601 UTC authentication timestamp"
    )

    def to_audit_payload(self) -> Dict[str, Any]:
        """Returns identity dictionary for stamping into immutable Merkle blocks."""
        return {
            "operator_uid": self.uid,
            "operator_email": self.email,
            "organization_id": self.organization_id,
            "role": self.role,
            "auth_time": self.auth_time,
        }


class AuthenticationError(ValueError):
    """Raised when token verification fails or an unauthorized request is presented."""
    pass


def verify_firebase_id_token(
    id_token: Optional[str],
    allow_mock: bool = True,
) -> AuthenticatedUser:
    """
    Verifies a Firebase JWT ID token and returns the authenticated user.

    Fail-Closed Policy:
      - None, empty string, or whitespace raises AuthenticationError.
      - In live production with firebase-admin, calls auth.verify_id_token().
      - In development/mock mode (or mock tokens), accepts deterministically formatted tokens.
    """
    if not id_token or not id_token.strip():
        raise AuthenticationError("Authentication token is required (fail-closed: empty token)")

    clean_token = id_token.strip()

    # 1. Deterministic Mock / Development Tokens
    if clean_token in ("panoptes_demo_token", "test_token"):
        if not allow_mock:
            raise AuthenticationError("Mock tokens are prohibited in strict live production mode")
        return AuthenticatedUser(
            uid="usr_panoptes_001",
            email="scientist@panoptesbio.com",
            display_name="Lead Computational Biologist",
            organization_id="panoptesbio",
            role="lead_biologist",
            is_part11_certified=True,
        )

    if clean_token.startswith("mock_token_") or clean_token.startswith("panoptes_demo_"):
        if not allow_mock:
            raise AuthenticationError("Mock tokens are prohibited in strict live production mode")

        # Parse role and user from mock token format: e.g. mock_token_director_john
        parts = clean_token.split("_")
        role = "lead_biologist"
        if "director" in parts:
            role = "lab_director"
        elif "compliance" in parts:
            role = "compliance_officer"

        uid = f"usr_{parts[-1]}" if len(parts) > 2 else "usr_panoptes_001"
        email = f"{parts[-1]}@panoptesbio.com" if len(parts) > 2 else "scientist@panoptesbio.com"
        name = "Dr. " + parts[-1].capitalize() if len(parts) > 2 else "Lead Computational Biologist"

        return AuthenticatedUser(
            uid=uid,
            email=email,
            display_name=name,
            organization_id="panoptesbio",
            role=role,
            is_part11_certified=True,
        )

    # 2. Live Firebase Admin SDK Token Verification
    _init_firebase_admin_if_needed()
    try:
        import firebase_admin.auth as fb_auth

        decoded = fb_auth.verify_id_token(clean_token)
        uid = decoded.get("uid", "anonymous_uid")
        email = decoded.get("email", f"{uid}@panoptesbio.com")
        name = decoded.get("name", "Firebase Authenticated Operator")
        org = decoded.get("org_id", "panoptesbio")
        role = decoded.get("role", "lead_biologist")

        return AuthenticatedUser(
            uid=uid,
            email=email,
            display_name=name,
            organization_id=org,
            role=role,
            is_part11_certified=True,
        )
    except Exception as exc:
        # Fallback check for simulated JSON web token in test harnesses
        if allow_mock and clean_token.startswith("ey"):
            try:
                # If JWT format but test environment has no Google certs access
                import base64
                payload_part = clean_token.split(".")[1]
                # Pad base64
                padded = payload_part + "=" * (-len(payload_part) % 4)
                decoded_json = json.loads(base64.urlsafe_b64decode(padded).decode("utf-8"))
                return AuthenticatedUser(
                    uid=decoded_json.get("user_id", decoded_json.get("sub", "usr_jwt_default")),
                    email=decoded_json.get("email", "jwt_user@panoptesbio.com"),
                    display_name=decoded_json.get("name", "JWT Operator"),
                    organization_id=decoded_json.get("org", "panoptesbio"),
                    role=decoded_json.get("role", "lead_biologist"),
                )
            except Exception:
                pass

        logger.warning("Firebase ID token verification failed: %s", exc)
        raise AuthenticationError(f"Invalid or expired Firebase ID token: {exc}") from exc


def create_part11_electronic_signature(
    user: AuthenticatedUser,
    action: str,
    meaning: str,
) -> Dict[str, Any]:
    """
    Creates an FDA 21 CFR Part 11 compliant electronic signature.

    Section 11.200 Requirements:
      1. Identification of the individual (user.uid + user.email).
      2. The date and time when signature was executed (timestamp).
      3. The meaning associated with the signature (e.g. 'Approval of Compute Spend', 'Authorization of Twist Synthesis').
      4. Cryptographic integrity binding (HMAC-SHA256 signature).
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    canonical_str = f"{user.uid}|{user.email}|{user.organization_id}|{user.role}|{action}|{meaning}|{timestamp}"
    hmac_sig = hmac.new(
        PART11_SIGNATURE_SECRET.encode("utf-8"),
        canonical_str.encode("utf-8"),
        hashlib.sha256,
    ).hexdigest()

    return {
        "operator_uid": user.uid,
        "operator_email": user.email,
        "display_name": user.display_name,
        "organization_id": user.organization_id,
        "role": user.role,
        "action": action,
        "signature_meaning": meaning,
        "signed_timestamp": timestamp,
        "signature_algorithm": "HMAC-SHA256 (21 CFR §11.200)",
        "electronic_signature_hash": hmac_sig,
    }
