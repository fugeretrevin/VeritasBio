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

How it reaches Modal specifically
---------------------------------
Modal does not use `requests` for its control plane; it speaks gRPC over HTTP/2 through
`grpclib`. `truststore.inject_into_ssl()` still covers it, but by a different route, so
it is worth writing down:

  * `grpclib.client.Channel._get_default_ssl_context()` calls
    `ssl.create_default_context(cafile=certifi.where())` — that `certifi.where()` is the
    hardcoded default that fails on an intercepted network.
  * `ssl.create_default_context` instantiates the *module global* `ssl.SSLContext`,
    which is exactly what `inject_into_ssl()` rebinds. Because the global is resolved at
    call time, the injection takes effect as long as it happens before the channel is
    built — the `cafile` argument is then satisfied by `truststore`'s context while the
    OS store supplies the anchors.

So the ordering requirement is "inject before the first channel is constructed," which
is looser than it sounds: `modal` is safe to import first. What is *not* safe is calling
`enable_os_trust_store()` before the opt-in env var exists; see the note on its return
value not being memoized.

Verified
--------
With this enabled, a 48-campaign end-to-end run completed against live Modal A10G
ESM-2 inference (410.71 GPU-seconds, engine `esm2_650m`) and 384 live Muse planner
calls (2,276,858 provider-reported tokens). Without it, both fail closed.

Re-verified 2026-11-01: with the env var set, `modal.Function.from_name` hydrates
`veritasbio-esm/score_variants_remote` and a 5-variant A10G call returns in ~23s with
`engine=esm2_650m, is_mock=False`. With it unset, the same call raises
`ConnectionError: Could not connect to the Modal server.`
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

TRUSTSTORE_ENV_VAR = "VERITAS_USE_OS_TRUSTSTORE"

# The process-wide injection has happened and is in effect.
_INJECTED: bool = False
# Injection was attempted and cannot succeed in this process (missing package, or it
# raised). Cached so a per-call retry does not re-log the same warning forever.
_UNAVAILABLE: bool = False


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

    A "not opted in" answer is deliberately **not** memoized. Only the successful
    injection is, because that is the thing that must happen at most once. Caching the
    negative broke Modal for days: `workers/esm_worker.py` calls this at import time,
    before any entry point has run `load_dotenv()`, so the first call saw an unset
    VERITAS_USE_OS_TRUSTSTORE, cached False, and every later call -- including the one
    the entry point makes *after* loading .env -- returned that stale False. certifi
    stayed in place and the gRPC handshake kept failing.
    """
    global _INJECTED, _UNAVAILABLE
    if _INJECTED:
        return True
    if _UNAVAILABLE:
        return False
    if not (force or os_trust_store_requested()):
        # Not cached: the opt-in may arrive later, via load_dotenv() or os.environ.
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
        _UNAVAILABLE = True
        return False
    except Exception as exc:
        logger.warning(
            "truststore injection failed (%s: %s); continuing with certifi",
            type(exc).__name__, exc,
        )
        _UNAVAILABLE = True
        return False

    logger.info(
        "TLS verification now uses the OS certificate store (%s is set). "
        "Verification remains enabled.", TRUSTSTORE_ENV_VAR,
    )
    _INJECTED = True
    return True


def os_trust_store_active() -> bool:
    """True when the injection has actually been performed in this process."""
    return _INJECTED
