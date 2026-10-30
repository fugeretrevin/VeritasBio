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
from typing import Any, Dict, Optional, Tuple

from pydantic import BaseModel, Field

from controller.firebase_config import get_firebase_config

logger = logging.getLogger(__name__)

# HMAC secret for electronic signatures. There is deliberately NO committed
# default: a hardcoded secret in a public repo means every signature it produces
# is forgeable offline by anyone with the source. When unset, signing refuses.
PART11_SIGNATURE_SECRET = os.getenv("PART11_SIGNATURE_SECRET", "")

# Mock identities are a development convenience and must be opted into explicitly.
# Previously `allow_mock` defaulted to True and no caller ever passed False, so any
# token beginning `mock_token_` was accepted AND its role was parsed out of the
# attacker-supplied string (mock_token_director_x -> lab_director).
MOCK_AUTH_ENV_VAR = "VERITAS_ALLOW_MOCK_AUTH"


def mock_auth_enabled() -> bool:
    """True only when mock identities have been explicitly enabled by env var."""
    return os.getenv(MOCK_AUTH_ENV_VAR, "").strip().lower() in ("1", "true", "yes", "on")

USERS_STORE_PATH = Path(__file__).resolve().parent.parent / "ledger" / "users.json"


def _load_users_store() -> Dict[str, Dict[str, Any]]:
    """Loads persisted user accounts from ledger/users.json."""
    if not USERS_STORE_PATH.is_file():
        return {}
    try:
        return json.loads(USERS_STORE_PATH.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _save_users_store(users: Dict[str, Dict[str, Any]]) -> None:
    """Saves persisted user accounts to ledger/users.json."""
    USERS_STORE_PATH.parent.mkdir(parents=True, exist_ok=True)
    USERS_STORE_PATH.write_text(json.dumps(users, indent=2), encoding="utf-8")


def _hash_password(password: str, salt: str) -> str:
    """Computes PBKDF2-HMAC-SHA256 password digest."""
    dk = hashlib.pbkdf2_hmac(
        "sha256",
        password.encode("utf-8"),
        salt.encode("utf-8"),
        100_000,
    )
    return dk.hex()


def register_user_account(
    email: str,
    password: str,
    display_name: str,
    organization_id: str = "panoptesbio",
    role: str = "lead_biologist",
) -> Tuple["AuthenticatedUser", str]:
    """
    Registers a new operator account with PBKDF2-HMAC-SHA256 password hashing
    and returns (AuthenticatedUser, session_token).
    """
    clean_email = email.strip().lower()
    if not clean_email or "@" not in clean_email:
        raise AuthenticationError("Please provide a valid email address.")
    if not password or len(password) < 6:
        raise AuthenticationError("Password must be at least 6 characters.")
    clean_name = display_name.strip() or clean_email.split("@")[0]

    users = _load_users_store()
    if clean_email in users:
        if hmac.compare_digest(_hash_password(password, users[clean_email].get("password_salt", "")), users[clean_email].get("password_hash", "")): return authenticate_user_account(clean_email, password)
        raise AuthenticationError("An account with this email already exists. Please sign in.")

    salt = hashlib.sha256(f"{clean_email}:{time.time()}".encode("utf-8")).hexdigest()[:24]
    pw_hash = _hash_password(password, salt)
    uid = "usr_" + hashlib.sha256(clean_email.encode("utf-8")).hexdigest()[:12]

    record = {
        "uid": uid,
        "email": clean_email,
        "display_name": clean_name,
        "organization_id": organization_id.strip() or "panoptesbio",
        "role": role if role in ("lead_biologist", "lab_director", "compliance_officer") else "lead_biologist",
        "password_salt": salt,
        "password_hash": pw_hash,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    users[clean_email] = record
    _save_users_store(users)

    user = AuthenticatedUser(
        uid=record["uid"],
        email=record["email"],
        display_name=record["display_name"],
        organization_id=record["organization_id"],
        role=record["role"],
        is_part11_certified=True,
    )
    token = f"panoptes_session_{uid}"
    return user, token


def authenticate_user_account(email: str, password: str) -> Tuple["AuthenticatedUser", str]:
    """
    Authenticates an existing operator account against ledger/users.json
    and returns (AuthenticatedUser, session_token).
    """
    clean_email = email.strip().lower()
    users = _load_users_store()
    record = users.get(clean_email)
    if not record:
        raise AuthenticationError("No account found for this email. Please create an account first.")

    expected_hash = record.get("password_hash", "")
    salt = record.get("password_salt", "")
    if not hmac.compare_digest(_hash_password(password, salt), expected_hash):
        raise AuthenticationError("Incorrect password. Please try again.")

    user = AuthenticatedUser(
        uid=record["uid"],
        email=record["email"],
        display_name=record["display_name"],
        organization_id=record.get("organization_id", "panoptesbio"),
        role=record.get("role", "lead_biologist"),
        is_part11_certified=True,
    )
    token = f"panoptes_session_{user.uid}"
    return user, token

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
    allow_mock: Optional[bool] = None,
) -> AuthenticatedUser:
    """
    Verify a Firebase JWT ID token and return the authenticated user.

    Fail-closed policy:
      - None, empty, or whitespace raises AuthenticationError.
      - Live verification goes through firebase-admin `auth.verify_id_token()`,
        which checks the signature.
      - Mock identities are accepted ONLY when explicitly enabled, via
        `allow_mock=True` or the VERITAS_ALLOW_MOCK_AUTH environment variable.
        `allow_mock=None` (the default) consults the environment, which is off
        unless set. Mock identities are never Part 11 certified, so they cannot
        mint electronic signatures.
      - A token that fails live verification is rejected. An earlier revision
        base64-decoded any token starting with "ey" and trusted its claims,
        including `role`, without verifying the signature; that path is removed.
    """
    if not id_token or not id_token.strip():
        raise AuthenticationError("Authentication token is required (fail-closed: empty token)")

    clean_token = id_token.strip()
    permit_mock = mock_auth_enabled() if allow_mock is None else bool(allow_mock)

    # 0. Registered User Session Tokens
    if clean_token.startswith("panoptes_session_"):
        target_uid = clean_token.replace("panoptes_session_", "", 1)
        users = _load_users_store()
        for rec in users.values():
            if rec.get("uid") == target_uid:
                return AuthenticatedUser(
                    uid=rec["uid"],
                    email=rec["email"],
                    display_name=rec.get("display_name", rec["email"]),
                    organization_id=rec.get("organization_id", "panoptesbio"),
                    role=rec.get("role", "lead_biologist"),
                    is_part11_certified=True,
                )

    # 1. Deterministic Mock / Development Tokens
    if clean_token in ("panoptes_demo_token", "test_token"):
        if not permit_mock:
            raise AuthenticationError(
                "Mock tokens are rejected. Set VERITAS_ALLOW_MOCK_AUTH=1 for local "
                "development, or supply a real Firebase ID token."
            )
        return AuthenticatedUser(
            uid="usr_panoptes_001",
            email="scientist@panoptesbio.com",
            display_name="Lead Computational Biologist (mock)",
            organization_id="panoptesbio",
            role="lead_biologist",
            is_part11_certified=False,
        )

    if clean_token.startswith("mock_token_") or clean_token.startswith("panoptes_demo_"):
        if not permit_mock:
            raise AuthenticationError(
                "Mock tokens are rejected. Set VERITAS_ALLOW_MOCK_AUTH=1 for local "
                "development, or supply a real Firebase ID token."
            )

        # Role is parsed from the token string, which the caller controls. This is
        # acceptable only because mock auth is development-only and grants no
        # Part 11 certification; it must never be enabled in a deployment.
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
            is_part11_certified=False,
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
        # No unverified fallback. An earlier revision caught this exception and, for
        # any token starting with "ey", base64-decoded the JWT payload and trusted
        # its claims - including `role` - without checking the signature. That let
        # anyone mint an arbitrary identity by crafting a token body, so the path is
        # gone: a token that fails signature verification is rejected.
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

    Refuses when PART11_SIGNATURE_SECRET is unset - a signature keyed on an empty
    or committed secret is forgeable by anyone with the source, which is worse than
    no signature because it looks authoritative. Also refuses for identities that
    are not Part 11 certified, which includes every mock identity.
    """
    if not PART11_SIGNATURE_SECRET:
        raise AuthenticationError(
            "PART11_SIGNATURE_SECRET is not configured; refusing to mint an electronic "
            "signature. Set it in the environment to enable signing."
        )
    if not user.is_part11_certified:
        raise AuthenticationError(
            f"Identity {user.uid!r} is not Part 11 certified (mock or unverified); "
            "refusing to mint an electronic signature."
        )

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
