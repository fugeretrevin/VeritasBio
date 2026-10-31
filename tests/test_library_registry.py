"""
tests/test_library_registry.py - the benchmark is library-aware, with no active-library state.

Guards the refactor that made a new combinatorial library benchmarkable through data and
configuration:

  1. The `ACTIVE_LIBRARY` / `ACTIVE_SITE_POSITIONS` / `PARENT_RESIDUES` module globals are
     gone, and every registry entry round-trips into the canonical
     `veritas.datasets.LibraryDesign` against its real wild-type sequence.
  2. Two libraries are usable in one process. Under the globals this was impossible:
     each `configure_library()` call clobbered the previous one, so a cross-protein study
     could not hold two universes at once.
  3. THE EXIT GATE. trpb, gb1 and her2 still produce exactly the eligible universes and
     top-1% thresholds they produced before the refactor. These numbers are asserted
     literally, so a future change to ingestion cannot move them unnoticed.
  4. HER2's 25,114 unlabelled rows stay out of the universe and its 5,832 floor
     measurements are recorded as censored rather than exact.
"""

from __future__ import annotations

import unittest
from typing import Dict

import benchmark_panoptes_isolated as bench
from veritas.datasets import FitnessDirection, SelectionMode
from veritas.datasets.errors import WildTypeMismatchError

#: Captured from the pre-refactor code (commit 38e98be) and re-verified after it. Moving
#: any of these means the eligible universe or the top-1% definition changed, which
#: invalidates every previously published number for that library.
EXPECTED: Dict[str, Dict[str, object]] = {
    "trpb": {
        "n_eligible": 159_999,
        "top1_threshold": 0.3493993794399995,
        "n_top1": 1_600,
        "global_best": 1.0,
        "rows_parsed": 159_999,
        "one_hot_width": 80,
        "first_id": "F184A",
        "last_id": "V227Y:S228Y",
        "unlabelled_excluded": 0,
        "censored_at_floor": 0,
    },
    "gb1": {
        "n_eligible": 149_360,
        "top1_threshold": 2.152377110000002,
        "n_top1": 1_494,
        "global_best": 8.761966,
        "rows_parsed": 149_361,
        "one_hot_width": 80,
        "first_id": "D40A",
        "last_id": "V54Y",
        "unlabelled_excluded": 0,
        "censored_at_floor": 0,
    },
    "her2": {
        "n_eligible": 11_276,
        "top1_threshold": 2.4882,
        "n_top1": 113,
        "global_best": 15.0842,
        "rows_parsed": 11_277,
        "one_hot_width": 200,
        "first_id": "A106T:M107H:D108V",
        "last_id": "W99Y:G101V:Y105F:A106Q:M107F:D108W",
        # The held-out split: absent labels, never admitted as measurements of zero.
        "unlabelled_excluded": 25_114,
        # Left-censored at the assay floor of 0.1, value preserved.
        "censored_at_floor": 5_832,
    },
}


def _load(key: str):
    spec = bench.LIBRARIES[key]
    return spec, bench.load_universe(
        library=spec,
        dataset_path=spec.dataset,
        provenance_mode="acknowledge-unflagged",
        imputed_ids_path=None,
    )


class TestNoActiveLibraryState(unittest.TestCase):
    def test_the_active_library_globals_are_gone(self):
        for name in ("ACTIVE_LIBRARY", "ACTIVE_SITE_POSITIONS", "PARENT_RESIDUES",
                     "PARENT_ACCESSION", "PARENT_METADATA_JSON", "STRUCTURE_PDB",
                     "PUBLISHED_MEASURED_COUNT", "PUBLISHED_IMPUTED_COUNT",
                     "PUBLISHED_GLOBAL_OPTIMUM", "FITNESS_DIRECTION",
                     "configure_library"):
            with self.subTest(name=name):
                self.assertFalse(
                    hasattr(bench, name),
                    f"{name} is back at module scope; the library must be threaded "
                    f"through explicitly so two libraries can coexist in one process",
                )

    def test_every_registry_entry_is_internally_consistent(self):
        for key, spec in bench.LIBRARIES.items():
            with self.subTest(library=key):
                self.assertEqual(spec.key, key)
                self.assertEqual(len(spec.parent_residues), spec.n_sites)
                self.assertEqual(sorted(spec.positions), list(spec.positions))
                self.assertEqual(len(set(spec.positions)), spec.n_sites)
                self.assertEqual(spec.one_hot_width, spec.n_sites * 20)
                self.assertIn(spec.fitness_direction,
                              [d.value for d in FitnessDirection])
                self.assertIn(spec.selection_mode, (SelectionMode.POOL_RESTRICTED,
                                                    SelectionMode.GENERATIVE))
                self.assertTrue(spec.dataset.is_file(), f"missing dataset {spec.dataset}")

    def test_registry_entries_convert_to_the_canonical_schema(self):
        """to_design() validates the declared design against the real WT sequence."""
        for key in sorted(bench.LIBRARIES):
            with self.subTest(library=key):
                _metadata, sequence, _pdb, design = bench._load_panoptes_context(
                    bench.LIBRARIES[key], pdb=False,
                )
                spec = bench.LIBRARIES[key]
                self.assertEqual(design.assay_id, key)
                self.assertEqual(design.site_positions, spec.positions)
                self.assertEqual(design.wild_type_combo, spec.parent_residues)
                self.assertEqual(design.wild_type_sequence, sequence)
                self.assertEqual(design.selection_mode, spec.selection_mode)
                self.assertEqual(design.n_sites, spec.n_sites)

    def test_a_library_declaring_the_wrong_parent_is_refused(self):
        import dataclasses

        spec = bench.LIBRARIES["trpb"]
        wrong = dataclasses.replace(spec, parent_residues="AAAA")
        _metadata, sequence, _pdb, _design = bench._load_panoptes_context(spec, pdb=False)
        with self.assertRaises(WildTypeMismatchError):
            wrong.to_design(sequence)

    def test_the_benchmark_encoders_refill_omitted_wild_type_sites(self):
        """A mutation string omits WT positions; the genotype is always full length."""
        trpb = bench.LIBRARIES["trpb"]
        self.assertEqual(bench._combo_from_notation("V183A:F184A:V227A", trpb), "AAAS")
        self.assertEqual(bench._combo_from_notation("V183A", trpb), "AFVS")
        self.assertEqual(bench._combo_from_notation("WT", trpb), "VFVS")
        self.assertEqual(
            bench._canonical_variant_id("AAAS", trpb), "V183A:F184A:V227A",
        )
        # And an unequal-length genotype is refused rather than truncated or padded.
        with self.assertRaises(bench.BenchmarkError):
            bench._canonical_variant_id("AAA", trpb)
        with self.assertRaises(bench.BenchmarkError):
            bench._canonical_variant_id("AAAAA", trpb)

    def test_two_libraries_encode_independently_in_one_process(self):
        trpb, gb1 = bench.LIBRARIES["trpb"], bench.LIBRARIES["gb1"]
        self.assertEqual(bench._combo_from_notation("V183A", trpb), "AFVS")
        self.assertEqual(bench._combo_from_notation("V39A", gb1), "ADGV")
        # Interleaved, to catch any cached or global per-library state.
        self.assertEqual(bench._combo_from_notation("V183A", trpb), "AFVS")
        self.assertIsNone(bench._combo_from_notation("V183A", gb1))
        self.assertIsNone(bench._combo_from_notation("V39A", trpb))
        self.assertEqual(bench._one_hot_matrix(["AFVS"], trpb).shape, (1, 80))
        self.assertEqual(
            bench._one_hot_matrix(["WGGDGFYAMD"], bench.LIBRARIES["her2"]).shape, (1, 200),
        )

    def test_selection_modes_are_declared_and_distinct(self):
        self.assertEqual(bench.LIBRARIES["her2"].selection_mode,
                         SelectionMode.POOL_RESTRICTED)
        self.assertTrue(bench.LIBRARIES["her2"].pool_restricted)
        for key in ("trpb", "gb1"):
            self.assertEqual(bench.LIBRARIES[key].selection_mode, SelectionMode.GENERATIVE)
            self.assertFalse(bench.LIBRARIES[key].pool_restricted)


class TestEligibleUniversesAreUnchanged(unittest.TestCase):
    """THE EXIT GATE: the three shipped libraries must not move."""

    @classmethod
    def setUpClass(cls) -> None:
        # Loaded once and held together. Three real universes coexisting in one process
        # is itself the property the globals made impossible.
        cls.loaded = {key: _load(key) for key in sorted(EXPECTED)}

    def test_eligible_universe_and_top1_threshold_are_identical(self):
        for key, expected in EXPECTED.items():
            spec, universe = self.loaded[key]
            evaluator = bench.Evaluator(universe)
            with self.subTest(library=key):
                self.assertEqual(len(universe), expected["n_eligible"])
                self.assertEqual(universe.fingerprint.unique_variants, expected["n_eligible"])
                self.assertEqual(universe.fingerprint.rows_parsed, expected["rows_parsed"])
                self.assertEqual(evaluator.top1_threshold, expected["top1_threshold"])
                self.assertEqual(evaluator.n_top1, expected["n_top1"])
                self.assertEqual(evaluator.global_best, expected["global_best"])
                self.assertEqual(universe.variant_ids[0], expected["first_id"])
                self.assertEqual(universe.variant_ids[-1], expected["last_id"])
                self.assertEqual(universe.one_hot.shape,
                                 (expected["n_eligible"], expected["one_hot_width"]))
                self.assertEqual(universe.duplicates_dropped, 0)
                self.assertEqual(universe.library.key, key)
                self.assertIs(universe.library, spec)

    def test_every_genotype_is_full_length_over_its_library(self):
        for key, (spec, universe) in self.loaded.items():
            with self.subTest(library=key):
                lengths = {len(c) for c in universe.combos.values()}
                self.assertEqual(lengths, {spec.n_sites})

    def test_her2_absent_labels_and_censoring_are_recorded_not_conflated(self):
        spec, universe = self.loaded["her2"]
        self.assertEqual(universe.unlabelled_excluded, EXPECTED["her2"]["unlabelled_excluded"])
        self.assertEqual(universe.censored_at_floor, EXPECTED["her2"]["censored_at_floor"])
        # No absent label leaked into the fitness table as a measured zero.
        self.assertNotIn(0.0, set(universe.fitness.values()))
        self.assertEqual(min(universe.fitness.values()), spec.detection_floor)
        # Censored values keep their measured value; only the flag records the limit.
        at_floor = sum(1 for v in universe.fitness.values() if v == spec.detection_floor)
        self.assertEqual(at_floor, universe.censored_at_floor)

    def test_trpb_and_gb1_declare_no_label_split(self):
        for key in ("trpb", "gb1"):
            _spec, universe = self.loaded[key]
            with self.subTest(library=key):
                self.assertEqual(universe.unlabelled_excluded, 0)
                self.assertEqual(universe.censored_at_floor, 0)


class TestRealDatasetsIngestThroughTheSchema(unittest.TestCase):
    """
    The deliverable claim: a shipped library ingests through `veritas.datasets` from its
    registry entry and the real wild-type sequence alone - no benchmark code, no
    optimiser code, and no per-library branch in either.
    """

    @classmethod
    def setUpClass(cls) -> None:
        from veritas.datasets import ingest_dataset

        cls.data = {}
        for key in sorted(bench.LIBRARIES):
            spec = bench.LIBRARIES[key]
            _metadata, sequence, _pdb, _design = bench._load_panoptes_context(spec, pdb=False)
            cls.data[key] = ingest_dataset(spec.ingest_config(sequence))

    def test_record_counts_agree_with_the_benchmark_universe(self):
        for key, expected in EXPECTED.items():
            with self.subTest(library=key):
                self.assertEqual(len(self.data[key].records), expected["n_eligible"])
                self.assertEqual(self.data[key].one_hot().shape,
                                 (expected["n_eligible"], expected["one_hot_width"]))

    def test_no_duplicate_genotypes_in_any_shipped_dataset(self):
        """The default policy REFUSES duplicates, so ingesting at all proves there are none."""
        for key, data in self.data.items():
            with self.subTest(library=key):
                self.assertEqual(data.report.duplicate_aggregation, "reject")
                self.assertEqual(data.report.duplicate_groups, 0)
                ids = [r.variant_id for r in data.records]
                self.assertEqual(len(set(ids)), len(ids))

    def test_coverage_is_reported_rather_than_filled_in(self):
        self.assertTrue(self.data["trpb"].report.library_is_complete)
        self.assertAlmostEqual(self.data["trpb"].report.coverage_fraction, 1.0)
        # GB1 is missing ~6.6% of its product. Those combinations are ABSENT from the
        # eligible set, not imputed into it.
        gb1 = self.data["gb1"].report
        self.assertFalse(gb1.library_is_complete)
        self.assertEqual(gb1.combinations_missing, 20 ** 4 - 1 - 149_360)
        self.assertAlmostEqual(gb1.coverage_fraction, 0.933506, places=6)

    def test_her2_labels_censoring_and_pool_are_all_distinct(self):
        from veritas.datasets import (
            Censoring,
            LabelStatus,
            SchemaError,
            SelectionMode as Mode,
        )

        data = self.data["her2"]
        self.assertEqual(data.report.rows_absent_label, 25_114)
        self.assertEqual(data.report.left_censored, 5_832)
        self.assertEqual(data.report.exact, 11_276 - 5_832)
        # Absent labels carry no number at all, so they cannot be read as measured zeros.
        for record in data.absent_labels[:50]:
            self.assertIs(record.fitness, None)
            self.assertEqual(record.label_status, LabelStatus.ABSENT)
            self.assertEqual(record.censoring, Censoring.NOT_APPLICABLE)
        table = data.fitness_table()
        self.assertEqual(len(table), 11_276)
        self.assertNotIn(0.0, set(table.values()))
        # A floor measurement keeps its value and is flagged, not clamped or dropped.
        censored = [r for r in data.records if r.censoring is Censoring.LEFT_CENSORED]
        self.assertEqual(len(censored), 5_832)
        self.assertTrue(all(r.fitness == 0.1 for r in censored))
        # Pool mode is explicit, and a pool is only defined in pool mode.
        self.assertEqual(data.selection_mode, Mode.POOL_RESTRICTED)
        self.assertEqual(len(data.candidate_pool()), 11_276)
        # And a generative library refuses to hand out one, so a result on an enumerated
        # pool can never be reported as a result on unconstrained generation.
        with self.assertRaises(SchemaError):
            self.data["trpb"].candidate_pool()

    def test_generative_libraries_expose_constraints_not_a_pool(self):
        from veritas.datasets import SchemaError

        for key in ("trpb", "gb1"):
            with self.subTest(library=key):
                data = self.data[key]
                constraints = data.mutation_constraints()
                self.assertEqual(constraints["positions"],
                                 list(bench.LIBRARIES[key].positions))
                self.assertEqual(constraints["max_mutation_depth"], 4)
                with self.assertRaises(SchemaError):
                    data.candidate_pool()


if __name__ == "__main__":
    unittest.main()
