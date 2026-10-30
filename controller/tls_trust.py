"""
controller/tls_trust.py — opt-in OS certificate store for outbound TLS.

Why this exists
---------------
On networks that terminate TLS (corporate proxies, and consumer security suites such as
Norton or Kaspersky), outbound HTTPS presents a certificate signed by a locally
installed root CA. Windows trusts that CA; `certifi` — which `requests`, `urllib3` and
the Modal client use by default — does not. The result is:

    SSLError(SSLCertVerificationError(1, '[SSL: CERTIFICATE_VERIFY_FAILED]
    certificate verify failed: unable to get local issuer certificate'))

on every call to the Meta/Llama API and to Modal, which silently degrades a campaign to
the mock planner and mock scoring. `docs/545.md` recorded the live end-to-end run as
impossible for this reason.

`truststore` resolves it the correct way: certificate verification stays fully enabled,
and trust anchors are read from the operating system store instead of certifi. This is
NOT a verification bypass — `verify=False` would be, and must never be used here.

Why opt-in
----------
Two reasons it is not applied unconditionally:

  * pip ships a vendored `ssl` shim that recurses infinitely when truststore patches
    SSLContext, so an always-on injection breaks package installation in the same venv.
  * On a network that is not intercepting TLS, certifi is the better default: it is
    independent of whatever happens to be in the machine's certificate store.

So callers set `VERITAS_USE_OS_TRUSTSTORE=1` when they are on an intercepted network.
`enable_os_trust_store()` is then safe to call repeatedly and from several entry points.

Verified
--------
With this enabled, a 48-campaign end-to-end run completed against live Modal A10G
ESM-2 inference (410.71 GPU-seconds, engine `esm2_650m`) and 384 live Muse planner
calls (2,276,858 provider-reported tokens). Without it, both fail closed.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

logger = logging.getLogger(__name__)

TRUSTSTORE_ENV_VAR = "VERITAS_USE_OS_TRUSTSTORE"

_ENABLED: Optional[bool] = None


def os_trust_store_requested() -> bool:
    """True when the caller has opted in via the environment."""
    return os.getenv(TRUSTSTORE_ENV_VAR, "").strip().lower() in ("1", "true", "yes", "on")


def enable_os_trust_store(force: bool = False) -> bool:
    """
    Route TLS verification through the OS certificate store, if opted in.

    Returns True when the OS store is in effect, False otherwise. Idempotent: the
    injection is performed at most once per process. `force=True` skips the env-var
    check, for callers that have established out of band that interception is present.

    Never raises: a missing `truststore` package or a failed injection leaves certifi in
    place and logs the reason, because failing to start is worse than failing to verify
    against the OS store.
    """
    global _ENABLED
    if _ENABLED is not None:
        return _ENABLED
    if not (force or os_trust_store_requested()):
        _ENABLED = False
        return False
    try:
        import truststore

        truststore.inject_into_ssl()
    except ImportError:
        logger.warning(
            "%s is set but the 'truststore' package is not installed, so TLS still "
            "verifies against certifi. Install it with: pip install truststore",
            TRUSTSTORE_ENV_VAR,
        )
        _ENABLED = False
        return False
    except Exception as exc:
        logger.warning(
            "truststore injection failed (%s: %s); continuing with certifi",
            type(exc).__name__, exc,
        )
        _ENABLED = False
        return False

    logger.info(
        "TLS verification now uses the OS certificate store (%s is set). "
        "Verification remains enabled.", TRUSTSTORE_ENV_VAR,
    )
    _ENABLED = True
    return True
