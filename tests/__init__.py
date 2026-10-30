"""VeritasBio test suite.

Marked as a package so the suite can be discovered from any working directory with:
    python -m unittest discover -s tests -t .

Mock identities and a signing secret are enabled here, for the test process only.
Production code defaults both off: `controller.auth` rejects mock tokens unless
VERITAS_ALLOW_MOCK_AUTH is set, and refuses to mint electronic signatures unless
PART11_SIGNATURE_SECRET is configured. Setting them in this module keeps the suite
exercising the auth paths without weakening the shipped defaults.
"""

import os

os.environ.setdefault("VERITAS_ALLOW_MOCK_AUTH", "1")
os.environ.setdefault("PART11_SIGNATURE_SECRET", "test-only-signing-secret-not-for-deployment")
os.environ.setdefault("VERITAS_CLOUDLAB_WEBHOOK_SECRET", "test-only-webhook-secret-not-for-deployment")
