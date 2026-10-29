"""
controller/firebase_config.py — Firebase / PanoptesBio Cloud Configuration.

Provides validated client and server credentials for:
  - Firebase Authentication (Google OAuth, Email/Password, FDA Part 11 electronic signatures)
  - Firebase Cloud Storage (panoptesbio.firebasestorage.app) for wet-lab manifests & PDB structures
  - Firebase Analytics (G-3CQM0RGVH3)
"""

from __future__ import annotations

import os
from typing import Any, Dict
from pydantic import BaseModel, Field


class FirebaseConfig(BaseModel):
    """Configuration for PanoptesBio Firebase instance."""

    api_key: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_API_KEY", "AIzaSyB4nBuxqnKFqFeyS_RZVCNlK6ox0iXEbq8")
    )
    auth_domain: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_AUTH_DOMAIN", "panoptesbio.firebaseapp.com")
    )
    project_id: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_PROJECT_ID", "panoptesbio")
    )
    storage_bucket: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_STORAGE_BUCKET", "panoptesbio.firebasestorage.app")
    )
    messaging_sender_id: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_MESSAGING_SENDER_ID", "561444549088")
    )
    app_id: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_APP_ID", "1:561444549088:web:e1be29926bdaaf9766fb85")
    )
    measurement_id: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_MEASUREMENT_ID", "G-3CQM0RGVH3")
    )
    service_account_path: str = Field(
        default_factory=lambda: os.getenv("FIREBASE_SERVICE_ACCOUNT_PATH", "")
    )

    def to_client_dict(self) -> Dict[str, str]:
        """Returns public client dictionary suitable for Firebase JS SDK."""
        return {
            "apiKey": self.api_key,
            "authDomain": self.auth_domain,
            "projectId": self.project_id,
            "storageBucket": self.storage_bucket,
            "messagingSenderId": self.messaging_sender_id,
            "appId": self.app_id,
            "measurementId": self.measurement_id,
        }

    def to_js_snippet(self) -> str:
        """Returns JavaScript snippet for embedding in web dashboards and Gradio HTML."""
        cfg = self.to_client_dict()
        return f"""
        const firebaseConfig = {{
            apiKey: "{cfg['apiKey']}",
            authDomain: "{cfg['authDomain']}",
            projectId: "{cfg['projectId']}",
            storageBucket: "{cfg['storageBucket']}",
            messagingSenderId: "{cfg['messagingSenderId']}",
            appId: "{cfg['appId']}",
            measurementId: "{cfg['measurementId']}"
        }};
        """


# Global singleton instance
_DEFAULT_FIREBASE_CONFIG: FirebaseConfig | None = None


def get_firebase_config() -> FirebaseConfig:
    """Returns the active FirebaseConfig singleton."""
    global _DEFAULT_FIREBASE_CONFIG
    if _DEFAULT_FIREBASE_CONFIG is None:
        _DEFAULT_FIREBASE_CONFIG = FirebaseConfig()
    return _DEFAULT_FIREBASE_CONFIG
