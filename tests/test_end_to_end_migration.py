"""
tests/test_end_to_end_migration.py - benchmark_end_to_end.py's manifest and provenance.

This file had NO test coverage at all while the isolated harness was modernised around
it, which is how it kept two defects the other harness had already fixed:

  * it wrote its manifest with `json.dumps(..., default=str)`, so any value that could
    not be serialised became an opaque string instead of a loud failure, and the
    manifest stopped being machine-readable;
  * it compared `provenance_status == "UNVERIFIED_INCLUDES_IMPUTED"` against a bare
    string. `LabelProvenance` is deliberately not a `str` subclass precisely so one
    opaque stamp cannot keep standing in for six materially different claims, and a bare
    string comparison against a member is always False - the limitation would have gone
    silently missing from the report.

A full campaign needs Modal credentials and a GPU, so these tests cover what can be
verified without them: argument parsing, the library threaded explicitly, the derived
scope checked against the declared design, the provenance mapping, and a manifest that
round-trips without a `default=` fallback.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import benchmark_end_to_end as e2e
import benchmark_panoptes_isolated as iso


_CONTEXT: dict = {}


def _context() -> dict:
    """
    The dataset, evaluator and derived scope, loaded once for the whole module.

    The prepared TrpB table is ~160k rows; loading it per test would make a cheap set of
    assertions about manifest shape cost more than the assertions are worth.
    """
    if not _CONTEXT:
        args = e2e.build_parser().parse_args([
            "--seeds", "1", "--rounds", "1", "--batch-size", "5", "--initial-size", "10",
            "--max-variants", "200", "--provenance-mode", "acknowledge-unflagged",
            "--dry-run",
        ])
        library = iso.resolve_library(args.library)
        universe = iso.load_universe(
            library=library,
            dataset_path=Path(args.dataset),
            provenance_mode=args.provenance_mode,
            imputed_ids_path=None,
            max_variants=args.max_variants,
            subsample_seed=args.subsample_seed,
        )
        scope = e2e.derive_library_scope(
            universe.variant_ids, _synthetic_parent(library),
        )
        _CONTEXT.update(
            args=args,
            library=library,
            universe=universe,
            evaluator=iso.Evaluator(universe),
            scope=scope,
            scope_check=e2e.check_scope_against_library(scope, library),
        )
    return _CONTEXT


def _synthetic_parent(library) -> str:
    """
    A parent sequence carrying the library's declared residues at its declared sites.

    Built here rather than fetched through `_load_panoptes_context`, so these tests do
    not depend on protein metadata, a structure file or the network.
    """
    residues = ["A"] * (max(library.positions) + 1)
    for position, residue in zip(library.positions, library.parent_residues):
        residues[position - 1] = residue
    return "".join(residues)


class TestEndToEndManifest(unittest.TestCase):
    """The manifest: parseable without `default=str`, and explicit about provenance."""

    _cached_manifest = None

    @classmethod
    def setUpClass(cls) -> None:
        context = _context()
        cls.args = context["args"]
        cls.library = context["library"]
        cls.universe = context["universe"]
        cls.evaluator = context["evaluator"]
        cls.scope = context["scope"]
        cls.scope_check = context["scope_check"]

    def _manifest(self) -> dict:
        """
        The manifest, built once for the class.

        Deep-copied on every call, so a test that deletes a key to check the
        completeness refusal cannot affect another test. Built once because
        `repository_state()` hashes the source tree, which is not something to repeat
        six times in a unit test.
        """
        import copy

        if getattr(type(self), "_cached_manifest", None) is None:
            type(self)._cached_manifest = self._build_manifest()
        return copy.deepcopy(type(self)._cached_manifest)

    def _build_manifest(self) -> dict:
        from datetime import datetime, timezone

        return e2e.build_manifest(
            run_id="e2e_test",
            args=self.args,
            library=self.library,
            universe=self.universe,
            scope=self.scope,
            scope_check=self.scope_check,
            evaluator=self.evaluator,
            status={"random": {"status": "available", "detail": "baseline: random"}},
            rows=[],
            production=e2e.summarise_production([]),
            production_facts=[],
            started=datetime.now(timezone.utc),
        )

    def test_the_manifest_serialises_without_a_default_fallback(self) -> None:
        manifest = self._manifest()
        # The defect: `default=str` silently rescued anything unserialisable. Encoding
        # with no `default=` at all is the assertion - it raises if a value needs one.
        encoded = json.dumps(iso._json_safe(manifest), indent=2)
        self.assertEqual(json.dumps(json.loads(encoded), indent=2), encoded)

    def test_the_manifest_passes_its_own_completeness_check(self) -> None:
        e2e.validate_manifest(self._manifest())

    def test_an_incomplete_manifest_is_refused_by_name(self) -> None:
        manifest = self._manifest()
        del manifest["label_provenance"]
        with self.assertRaises(iso.ManifestIncompleteError) as caught:
            e2e.validate_manifest(manifest)
        self.assertIn("label_provenance", str(caught.exception))

    def test_both_vocabularies_are_recorded_side_by_side(self) -> None:
        manifest = self._manifest()
        legacy = self.universe.provenance_status
        resolved = iso.map_legacy_provenance_status(legacy)
        # The new stamp leads...
        self.assertEqual(manifest["label_provenance"]["status"], resolved.value)
        # ...and the legacy string is KEPT, both inside the block and at the top level,
        # so documented greps and existing readers keep matching.
        self.assertEqual(manifest["label_provenance"]["legacy_status"], legacy)
        self.assertEqual(manifest["provenance_status"], legacy)
        self.assertIn(resolved.value, manifest["label_provenance"]["vocabulary"])

    def test_the_manifest_names_the_library_it_scored(self) -> None:
        manifest = self._manifest()
        self.assertEqual(manifest["library"]["key"], self.library.key)
        self.assertEqual(
            manifest["library"]["declared_positions"], list(self.library.positions)
        )
        self.assertEqual(
            manifest["library_scope"]["declared_vs_derived"]["derived_positions"],
            list(self.scope.positions),
        )

    def test_the_manifest_records_which_roots_the_run_resolved(self) -> None:
        manifest = self._manifest()
        self.assertIn("repository_root", manifest["environment"])
        self.assertIn("in_worktree", manifest["environment"])


class TestProvenanceIsAnEnumNotAString(unittest.TestCase):
    """The comparison that used to be made against a bare string."""

    def test_a_member_never_compares_equal_to_its_own_legacy_string(self) -> None:
        # This is why the old `universe.provenance_status == "UNVERIFIED_..."` branch had
        # to be replaced rather than wrapped: the two vocabularies cannot be compared.
        member = e2e.LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE
        self.assertNotEqual(member, "HISTORICAL_UNKNOWN_PROVENANCE")
        self.assertNotEqual(member, "UNVERIFIED_INCLUDES_IMPUTED")
        self.assertEqual(member.value, "HISTORICAL_UNKNOWN_PROVENANCE")

    def test_the_stamp_carries_the_new_claim_and_the_legacy_string(self) -> None:
        universe = _context()["universe"]
        stamp = e2e.provenance_stamp(universe)
        self.assertIn("provenance=HISTORICAL_UNKNOWN_PROVENANCE", stamp)
        self.assertIn("legacy=UNVERIFIED_INCLUDES_IMPUTED", stamp)

        # And the limitation fires off the enum, naming both.
        limitations = e2e._limitations(
            e2e.build_parser().parse_args(["--dry-run"]), universe, {},
        )
        self.assertTrue(
            any("HISTORICAL_UNKNOWN_PROVENANCE" in item
                and "UNVERIFIED_INCLUDES_IMPUTED" in item for item in limitations),
            limitations,
        )


class TestLibraryIsThreadedNotGlobal(unittest.TestCase):
    """`--library` already worked; what is checked here is that nothing reads globals."""

    def test_the_deleted_module_globals_are_not_referenced(self) -> None:
        for name in ("ACTIVE_LIBRARY", "ACTIVE_SITE_POSITIONS", "PARENT_RESIDUES"):
            self.assertFalse(hasattr(e2e, name), name)
            self.assertFalse(hasattr(iso, name), name)

    def test_every_registered_library_can_be_selected_and_has_its_own_dataset(self) -> None:
        for key in sorted(iso.LIBRARIES):
            args = e2e.build_parser().parse_args(["--library", key])
            self.assertEqual(args.library, key)
            spec = iso.resolve_library(key)
            self.assertEqual(spec.key, key)
            # Each library carries its own prepared table; main() substitutes it when
            # --library names one and --dataset was left at the default.
            self.assertTrue(Path(spec.dataset).is_file(), spec.dataset)

    def test_a_derived_scope_outside_the_declared_design_is_fatal(self) -> None:
        library = iso.resolve_library("trpb")
        rogue = e2e.LibraryScope(
            positions=(183, 184, 227, 228, 999),
            parent_residues={183: "V", 184: "F", 227: "V", 228: "S", 999: "A"},
            source="synthetic",
        )
        with self.assertRaises(e2e.E2EError) as caught:
            e2e.check_scope_against_library(rogue, library)
        self.assertIn("999", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
