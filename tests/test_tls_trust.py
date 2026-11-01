"""
Regression tests for `controller.tls_trust`, the OS-certificate-store opt-in.

The defect these lock in kept live Modal inference dead for days on an intercepted
network, with `VERITAS_USE_OS_TRUSTSTORE=1` sitting in `.env` the whole time.

`enable_os_trust_store()` memoized *every* outcome in a single `_ENABLED` tri-state,
including "the caller has not opted in". `workers/esm_worker.py` calls it at import
time, and that import happens before any entry point runs `load_dotenv()`. So the
first call saw an unset environment, cached `False`, and the call the entry point
made *after* loading `.env` returned that stale `False` instead of injecting. certifi
stayed in place and Modal's gRPC handshake kept failing with
CERTIFICATE_VERIFY_FAILED -- while the function reported, truthfully, that the OS
trust store was "not enabled".

The invariant that actually needs memoizing is the *injection*, which must happen at
most once per process. A negative answer must not be sticky, because the opt-in can
arrive later. These tests pin that distinction, and pin that nothing here ever
weakens verification.

`truststore.inject_into_ssl` is mocked throughout: it rebinds `ssl.SSLContext`
process-wide, which would leak into every other test in the suite.
"""
from __future__ import annotations

import unittest
from unittest import mock

import controller.tls_trust as tls_trust
from controller.tls_trust import (
    TRUSTSTORE_ENV_VAR,
    enable_os_trust_store,
    os_trust_store_active,
    os_trust_store_requested,
)


class _TrustStoreTestCase(unittest.TestCase):
    """Restores the module's process-wide injection state around each test."""

    def setUp(self) -> None:
        self._saved = (tls_trust._INJECTED, tls_trust._UNAVAILABLE)
        tls_trust._INJECTED = False
        tls_trust._UNAVAILABLE = False

    def tearDown(self) -> None:
        tls_trust._INJECTED, tls_trust._UNAVAILABLE = self._saved


class TestOptInParsing(_TrustStoreTestCase):

    def test_accepted_truthy_spellings(self) -> None:
        for raw in ("1", "true", "TRUE", "yes", "On", "  1  "):
            with self.subTest(raw=raw):
                with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: raw}):
                    self.assertTrue(os_trust_store_requested())

    def test_absent_or_falsy_is_not_opted_in(self) -> None:
        for raw in ("", "0", "false", "no", "off", "maybe"):
            with self.subTest(raw=raw):
                with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: raw}):
                    self.assertFalse(os_trust_store_requested())

    def test_missing_variable_is_not_opted_in(self) -> None:
        with mock.patch.dict("os.environ", {}, clear=True):
            self.assertFalse(os_trust_store_requested())


class TestNegativeAnswerIsNotMemoized(_TrustStoreTestCase):
    """
    THE REGRESSION. This is the exact sequence `workers/esm_worker.py` produces.
    """

    def test_opting_in_after_a_negative_call_still_injects(self) -> None:
        injected = mock.Mock()
        with mock.patch.dict("os.environ", {}, clear=True):
            # Import-time call, before load_dotenv(): correctly declines...
            self.assertFalse(enable_os_trust_store())
        self.assertFalse(os_trust_store_active())

        # ...and must NOT have poisoned the decision for the call that follows
        # load_dotenv(). Before the fix this second call returned a stale False.
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch("truststore.inject_into_ssl", injected):
                self.assertTrue(enable_os_trust_store())
        injected.assert_called_once()
        self.assertTrue(os_trust_store_active())

    def test_declining_never_injects(self) -> None:
        injected = mock.Mock()
        with mock.patch.dict("os.environ", {}, clear=True):
            with mock.patch("truststore.inject_into_ssl", injected):
                for _ in range(3):
                    self.assertFalse(enable_os_trust_store())
        injected.assert_not_called()
        self.assertFalse(os_trust_store_active())


class TestInjectionHappensAtMostOnce(_TrustStoreTestCase):

    def test_repeated_calls_inject_once(self) -> None:
        injected = mock.Mock()
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch("truststore.inject_into_ssl", injected):
                results = [enable_os_trust_store() for _ in range(5)]
        self.assertEqual(results, [True] * 5)
        injected.assert_called_once()

    def test_force_bypasses_the_env_check(self) -> None:
        injected = mock.Mock()
        with mock.patch.dict("os.environ", {}, clear=True):
            with mock.patch("truststore.inject_into_ssl", injected):
                self.assertTrue(enable_os_trust_store(force=True))
        injected.assert_called_once()

    def test_an_injected_process_stays_injected_if_the_var_is_later_unset(self) -> None:
        """The patch is process-wide and cannot be undone, so the report must not lie."""
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch("truststore.inject_into_ssl", mock.Mock()):
                self.assertTrue(enable_os_trust_store())
        with mock.patch.dict("os.environ", {}, clear=True):
            self.assertTrue(enable_os_trust_store())
            self.assertTrue(os_trust_store_active())


class TestFailureLeavesCertifiInPlace(_TrustStoreTestCase):
    """A failure must degrade to certifi, never to a weakened verification mode."""

    def test_missing_package_returns_false_and_does_not_raise(self) -> None:
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch.dict("sys.modules", {"truststore": None}):
                # `import truststore` against a None entry raises ImportError.
                self.assertFalse(enable_os_trust_store())
        self.assertFalse(os_trust_store_active())

    def test_a_raising_injection_returns_false_and_does_not_raise(self) -> None:
        boom = mock.Mock(side_effect=RuntimeError("recursive ssl shim"))
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch("truststore.inject_into_ssl", boom):
                self.assertFalse(enable_os_trust_store())
        self.assertFalse(os_trust_store_active())

    def test_a_hard_failure_is_not_retried(self) -> None:
        boom = mock.Mock(side_effect=RuntimeError("recursive ssl shim"))
        with mock.patch.dict("os.environ", {TRUSTSTORE_ENV_VAR: "1"}):
            with mock.patch("truststore.inject_into_ssl", boom):
                for _ in range(4):
                    self.assertFalse(enable_os_trust_store())
        boom.assert_called_once()


class TestNoVerificationBypass(unittest.TestCase):
    """
    `tls_trust` exists to move trust anchors, not to stop checking them. A future
    edit that reaches for the easy fix should fail here.
    """

    def test_module_source_never_disables_verification(self) -> None:
        import inspect
        source = inspect.getsource(tls_trust)
        # Strip the docstring, which names these constructs in order to forbid them.
        body = source.split('"""', 2)[-1]
        for forbidden in (
            "verify=False",
            "CERT_NONE",
            "check_hostname = False",
            "check_hostname=False",
            "_create_unverified_context",
            "trusted-host",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, body)


class TestEsmWorkerEnablesTheTrustStoreBeforeModal(unittest.TestCase):
    """
    `workers/esm_worker.py` is imported directly by scripts that never go through
    `controller.orchestrator`, so it has to load `.env` and opt in by itself. If this
    ordering regresses, Modal silently falls back to the mock worker again.
    """

    def test_esm_worker_loads_dotenv_before_enabling_the_trust_store(self) -> None:
        # Parsed, not grepped: the module docstring names all three of these while
        # documenting the defect, so a text search finds the prose, not the code.
        import ast
        import inspect
        import workers.esm_worker as esm_worker

        tree = ast.parse(inspect.getsource(esm_worker))

        def first_call(name: str):
            lines = [n.lineno for n in ast.walk(tree)
                     if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Name) and n.func.id == name]
            return min(lines) if lines else None

        def first_import(module: str):
            lines = [n.lineno for n in ast.walk(tree)
                     if isinstance(n, ast.Import)
                     and any(a.name == module for a in n.names)]
            return min(lines) if lines else None

        load_at = first_call("load_dotenv")
        enable_at = first_call("enable_os_trust_store")
        modal_at = first_import("modal")

        self.assertIsNotNone(load_at, "esm_worker must load .env itself")
        self.assertIsNotNone(enable_at, "esm_worker must opt in to the OS trust store")
        self.assertIsNotNone(modal_at, "esm_worker must import modal")
        self.assertLess(load_at, enable_at,
                        ".env must be loaded before the opt-in var is read")
        self.assertLess(enable_at, modal_at,
                        "the trust store must be injected before modal is imported")


if __name__ == "__main__":
    unittest.main()
