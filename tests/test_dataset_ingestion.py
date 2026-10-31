"""
tests/test_dataset_ingestion.py - the library-aware dataset schema and ingestion contract.

These tests fix the behaviour that makes a NEW combinatorial library benchmarkable through
data and configuration alone:

  * the canonical schema separates REQUIRED design facts from OPTIONAL provenance;
  * every validation REJECTS rather than coercing - wrong wild type, wrong alphabet, an
    indel, an inconsistent genotype length and a position outside the declared design are
    all errors, never a silently repaired record;
  * a mutation string OMITS wild-type positions, so 'V183A:F184A:V227A' is a FOUR-site
    genotype whose fourth site is wild type. Everything canonicalises to full length over
    the library's positions before anything compares two genotypes;
  * an absent label, a true measured zero and a detection floor are three different things
    and are never conflated;
  * POOL_RESTRICTED and GENERATIVE are declared modes, not an accident of the data.
"""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from veritas.datasets import (
    AssayMeta,
    Censoring,
    ColumnMap,
    DuplicateGenotypeError,
    FitnessDirection,
    IngestConfig,
    LabelStatus,
    LibraryDesign,
    SelectionMode,
    VariablePosition,
    ingest_dataset,
)
from veritas.datasets.errors import (
    AlphabetError,
    ConflictingMeasurementError,
    GenotypeLengthError,
    IndelError,
    LabelPolicyError,
    PositionOutsideLibraryError,
    SchemaError,
    WildTypeMismatchError,
)

AA20 = "ACDEFGHIKLMNPQRSTVWY"


def _wt_sequence(length: int = 240, **residues: str) -> str:
    """A filler sequence with specific residues pinned at specific 1-indexed positions."""
    seq = ["G"] * length
    for key, residue in residues.items():
        seq[int(key.lstrip("p")) - 1] = residue
    return "".join(seq)


def trpb_like_design(selection_mode: SelectionMode = SelectionMode.POOL_RESTRICTED) -> LibraryDesign:
    """A complete 20-letter four-site library, shaped like the TrpB active site."""
    wt = _wt_sequence(240, p183="V", p184="F", p227="V", p228="S")
    return LibraryDesign(
        target_id="TRPB_TM9D8",
        assay_id="trpb_four_site_johnston_2024",
        wild_type_sequence=wt,
        positions=tuple(
            VariablePosition(position=p, wild_type_aa=wt[p - 1], allowed_aa=tuple(AA20))
            for p in (183, 184, 227, 228)
        ),
        fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
        selection_mode=selection_mode,
    )


def binary_design() -> LibraryDesign:
    """A two-site binary library: the sparse-alphabet case."""
    wt = _wt_sequence(20, p5="A", p9="A")
    return LibraryDesign(
        target_id="TOY_BINARY",
        assay_id="toy_binary_v1",
        wild_type_sequence=wt,
        positions=(
            VariablePosition(position=5, wild_type_aa="A", allowed_aa=("A", "V")),
            VariablePosition(position=9, wild_type_aa="A", allowed_aa=("A", "V")),
        ),
        fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
        selection_mode=SelectionMode.POOL_RESTRICTED,
    )


def _write_csv(directory: Path, name: str, header: str, rows: list[str]) -> Path:
    path = directory / name
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


# ==================================================================================
# The design itself
# ==================================================================================
class TestLibraryDesign(unittest.TestCase):
    def test_complete_twenty_letter_four_site_library(self):
        design = trpb_like_design()
        self.assertEqual(design.n_sites, 4)
        self.assertEqual(design.site_positions, (183, 184, 227, 228))
        self.assertEqual(design.wild_type_combo, "VFVS")
        self.assertEqual(design.one_hot_width, 4 * 20)
        self.assertEqual(design.combinatorial_size, 20 ** 4)
        self.assertTrue(design.is_saturated_alphabet)

    def test_sparse_alphabet_library_is_not_saturated(self):
        design = binary_design()
        self.assertEqual(design.one_hot_width, 2 * 2)
        self.assertEqual(design.combinatorial_size, 4)
        self.assertFalse(design.is_saturated_alphabet)

    def test_mismatched_wild_type_reference_is_rejected(self):
        wt = _wt_sequence(240, p183="V")
        with self.assertRaises(WildTypeMismatchError):
            LibraryDesign(
                target_id="T",
                assay_id="a",
                wild_type_sequence=wt,
                positions=(VariablePosition(position=183, wild_type_aa="A", allowed_aa=tuple(AA20)),),
                fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
                selection_mode=SelectionMode.POOL_RESTRICTED,
            )

    def test_position_outside_the_sequence_is_rejected(self):
        with self.assertRaises(SchemaError):
            LibraryDesign(
                target_id="T",
                assay_id="a",
                wild_type_sequence="GGG",
                positions=(VariablePosition(position=9, wild_type_aa="G", allowed_aa=("G", "A")),),
                fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
                selection_mode=SelectionMode.POOL_RESTRICTED,
            )

    def test_duplicate_and_unordered_positions_are_rejected(self):
        wt = _wt_sequence(20, p5="A", p9="A")
        site_a = VariablePosition(position=5, wild_type_aa="A", allowed_aa=("A", "V"))
        site_b = VariablePosition(position=9, wild_type_aa="A", allowed_aa=("A", "V"))
        with self.assertRaises(SchemaError):
            LibraryDesign(
                target_id="T", assay_id="a", wild_type_sequence=wt,
                positions=(site_a, site_a),
                fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
                selection_mode=SelectionMode.POOL_RESTRICTED,
            )
        with self.assertRaises(SchemaError):
            LibraryDesign(
                target_id="T", assay_id="a", wild_type_sequence=wt,
                positions=(site_b, site_a),
                fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
                selection_mode=SelectionMode.POOL_RESTRICTED,
            )

    def test_allowed_alphabet_must_be_declared_and_legal(self):
        wt = _wt_sequence(20, p5="A")
        with self.assertRaises(SchemaError):
            VariablePosition(position=5, wild_type_aa="A", allowed_aa=())
        with self.assertRaises(AlphabetError):
            VariablePosition(position=5, wild_type_aa="A", allowed_aa=("A", "B"))
        # The wild-type residue must be representable, or a mutation string that omits
        # the site cannot be canonicalised back to full length.
        with self.assertRaises(SchemaError):
            LibraryDesign(
                target_id="T", assay_id="a", wild_type_sequence=wt,
                positions=(VariablePosition(position=5, wild_type_aa="A", allowed_aa=("V", "L")),),
                fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
                selection_mode=SelectionMode.POOL_RESTRICTED,
            )

    def test_generative_mode_declares_its_mutation_constraints(self):
        design = trpb_like_design(SelectionMode.GENERATIVE)
        self.assertTrue(design.is_generative)
        self.assertFalse(design.is_pool_restricted)
        constraints = design.mutation_constraints()
        self.assertEqual(constraints["positions"], [183, 184, 227, 228])
        self.assertEqual(constraints["max_mutation_depth"], 4)
        self.assertEqual(constraints["allowed_aa"]["183"], list(AA20))

    def test_pool_restricted_mode_is_distinct_from_generative(self):
        design = trpb_like_design(SelectionMode.POOL_RESTRICTED)
        self.assertTrue(design.is_pool_restricted)
        self.assertFalse(design.is_generative)

    def test_optional_metadata_is_separate_from_required_design(self):
        design = trpb_like_design()
        # Optional provenance defaults to 'unknown', never to a plausible guess.
        self.assertIsNone(design.assay.units)
        self.assertIsNone(design.assay.detection_floor)
        self.assertIsNone(design.assay.source_publication)
        enriched = design.with_assay(AssayMeta(
            units="normalised enrichment",
            transformation="none",
            detection_floor=0.1,
            source_publication="Mason et al. 2021",
            source_doi="10.1038/s41551-021-00699-9",
            dataset_version="v1",
            replicate_count=3,
        ))
        self.assertEqual(enriched.assay.detection_floor, 0.1)
        self.assertEqual(design.assay.detection_floor, None)  # frozen; original untouched


# ==================================================================================
# Notation: wild-type positions are OMITTED
# ==================================================================================
class TestNotationCanonicalisation(unittest.TestCase):
    def setUp(self):
        self.design = trpb_like_design()

    def test_mutation_string_omitting_wild_type_sites_is_full_length(self):
        # THE case that was got wrong: three tokens, four sites, fourth site is WT.
        combo = self.design.combo_from_mutations("V183A:F184A:V227A")
        self.assertEqual(combo, "AAAS")
        self.assertEqual(len(combo), self.design.n_sites)

    def test_variable_mutation_depth_all_canonicalise_to_full_length(self):
        cases = {
            "": "VFVS",
            "WT": "VFVS",
            "V183A": "AFVS",
            "V183A:S228G": "AFVG",
            "V183A:F184I:V227K": "AIKS",
            "V183A:F184I:V227K:S228G": "AIKG",
        }
        for notation, expected in cases.items():
            with self.subTest(notation=notation):
                combo = self.design.combo_from_mutations(notation)
                self.assertEqual(combo, expected)
                self.assertEqual(len(combo), 4)

    def test_canonical_id_round_trip_drops_wild_type_tokens(self):
        self.assertEqual(self.design.canonical_variant_id("AAAS"), "V183A:F184A:V227A")
        self.assertEqual(self.design.canonical_variant_id("AIKG"), "V183A:F184I:V227K:S228G")
        self.assertEqual(self.design.canonical_variant_id("VFVS"), "")
        for combo in ("AAAS", "AIKG", "AFVS"):
            with self.subTest(combo=combo):
                notation = self.design.canonical_variant_id(combo)
                self.assertEqual(self.design.combo_from_mutations(notation), combo)

    def test_token_order_does_not_change_the_genotype(self):
        self.assertEqual(
            self.design.combo_from_mutations("S228G:V183A"),
            self.design.combo_from_mutations("V183A:S228G"),
        )

    def test_explicit_wild_type_token_is_the_same_genotype(self):
        self.assertEqual(self.design.combo_from_mutations("V183A:F184F"), "AFVS")
        self.assertEqual(
            self.design.canonical_variant_id(self.design.combo_from_mutations("V183A:F184F")),
            "V183A",
        )

    def test_position_outside_the_declared_design_is_rejected(self):
        with self.assertRaises(PositionOutsideLibraryError):
            self.design.combo_from_mutations("S65T")

    def test_wrong_wild_type_residue_in_the_token_is_rejected(self):
        with self.assertRaises(WildTypeMismatchError):
            self.design.combo_from_mutations("A183V")

    def test_residue_outside_the_declared_alphabet_is_rejected(self):
        with self.assertRaises(AlphabetError):
            self.design.combo_from_mutations("V183X")
        # A sparse library does NOT infer its alphabet from what happens to be observed.
        binary = binary_design()
        self.assertEqual(binary.combo_from_mutations("A5V"), "VA")
        with self.assertRaises(AlphabetError):
            binary.combo_from_mutations("A5L")

    def test_two_tokens_at_one_position_are_rejected(self):
        with self.assertRaises(SchemaError):
            self.design.combo_from_mutations("V183A:V183L")

    def test_indels_are_detected_not_treated_as_substitutions(self):
        for notation in ("V183del", "183insA", "V183-", "V183_", "V183fs", "V183*", "V183"):
            with self.subTest(notation=notation):
                with self.assertRaises(IndelError):
                    self.design.combo_from_mutations(notation)

    def test_sequence_of_the_wrong_length_is_an_indel_not_a_substitution_set(self):
        wt = self.design.wild_type_sequence
        with self.assertRaises(IndelError):
            self.design.combo_from_sequence(wt[:-1])
        with self.assertRaises(IndelError):
            self.design.combo_from_sequence(wt + "A")

    def test_substitution_positions_are_detected_from_a_sequence(self):
        wt = self.design.wild_type_sequence
        mutated = list(wt)
        mutated[182] = "A"
        mutated[227] = "G"
        self.assertEqual(self.design.combo_from_sequence("".join(mutated)), "AFVG")

    def test_substitution_outside_the_library_positions_is_rejected(self):
        mutated = list(self.design.wild_type_sequence)
        mutated[64] = "T"
        with self.assertRaises(PositionOutsideLibraryError):
            self.design.combo_from_sequence("".join(mutated))

    def test_inconsistent_genotype_length_is_rejected_never_truncated(self):
        with self.assertRaises(GenotypeLengthError):
            self.design.canonical_variant_id("AAA")
        with self.assertRaises(GenotypeLengthError):
            self.design.canonical_variant_id("AAAAA")
        with self.assertRaises(GenotypeLengthError):
            self.design.hamming_distance("AAAS", "AAA")
        self.assertEqual(self.design.hamming_distance("AAAS", "VFVS"), 3)

    def test_one_hot_is_library_shaped_and_rejects_bad_genotypes(self):
        matrix = self.design.one_hot(["AAAS", "VFVS"])
        self.assertEqual(matrix.shape, (2, 80))
        self.assertTrue(np.all(matrix.sum(axis=1) == 4))
        with self.assertRaises(GenotypeLengthError):
            self.design.one_hot(["AAA"])
        with self.assertRaises(AlphabetError):
            self.design.one_hot(["XAAS"])

    def test_two_designs_do_not_share_state(self):
        """No module-level active library: two designs are usable in one process."""
        a, b = trpb_like_design(), binary_design()
        self.assertEqual(a.combo_from_mutations("V183A"), "AFVS")
        self.assertEqual(b.combo_from_mutations("A5V"), "VA")
        self.assertEqual(a.combo_from_mutations("V183A"), "AFVS")


# ==================================================================================
# Ingestion
# ==================================================================================
class TestIngestion(unittest.TestCase):
    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_ingests_a_complete_four_site_library(self):
        design = trpb_like_design()
        path = _write_csv(self.dir, "trpb.csv", "mutant,DMS_score", [
            "V183A:F184A:V227A,0.5",
            "V183A:F184I:V227K:S228G,1.0",
            "S228G,0.25",
        ])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(len(data.records), 3)
        self.assertEqual(
            [r.variant_id for r in data.records],
            ["S228G", "V183A:F184A:V227A", "V183A:F184I:V227K:S228G"],
        )
        for record in data.records:
            self.assertEqual(len(record.combo), 4)
            self.assertEqual(record.label_status, LabelStatus.MEASURED)
            self.assertEqual(record.censoring, Censoring.EXACT)
        self.assertEqual(data.fitness_table()["S228G"], 0.25)
        self.assertEqual(data.one_hot().shape, (3, 80))

    def test_binary_library_round_trips(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", [
            "A5V,0.2", "A9V,0.3", "A5V:A9V,0.9",
        ])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual({r.combo for r in data.records}, {"VA", "AV", "VV"})

    def test_missing_but_permissible_combinations_are_reported_not_invented(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A9V,0.3"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(len(data.records), 2)
        # Coverage excludes the wild type from both sides: it is a permitted genotype
        # but not a library member, so a complete library must not read as short by one.
        self.assertEqual(data.report.combinatorial_size, 4)
        self.assertEqual(data.report.combinations_present, 2)
        self.assertEqual(data.report.combinations_missing, 1)   # only A5V:A9V is absent
        self.assertAlmostEqual(data.report.coverage_fraction, 2 / 3)
        self.assertNotIn("A5V:A9V", data.fitness_table())
        self.assertFalse(data.report.library_is_complete)

    def test_a_complete_small_library_reads_as_complete(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", [
            "A5V,0.2", "A9V,0.3", "A5V:A9V,0.9",
        ])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(data.report.combinations_missing, 0)
        self.assertTrue(data.report.library_is_complete)
        self.assertAlmostEqual(data.report.coverage_fraction, 1.0)

    def test_wild_type_row_is_not_a_library_member(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["WT,1.0", "A5V,0.2"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual([r.variant_id for r in data.records], ["A5V"])
        self.assertEqual(data.report.wild_type_rows, 1)

    # ---- duplicates and conflicting measurements ---------------------------------
    def test_duplicate_genotype_is_rejected_by_default(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A5V,0.2"])
        with self.assertRaises(DuplicateGenotypeError):
            ingest_dataset(IngestConfig(design=design, source=path))

    def test_duplicate_detected_across_different_spellings_of_one_genotype(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A5V:A9A,0.2"])
        with self.assertRaises(DuplicateGenotypeError) as ctx:
            ingest_dataset(IngestConfig(design=design, source=path))
        self.assertIn("A5V", str(ctx.exception))

    def test_conflicting_measurements_are_rejected_under_aggregation(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A5V,0.9"])
        with self.assertRaises(ConflictingMeasurementError):
            ingest_dataset(IngestConfig(
                design=design, source=path,
                duplicate_policy="mean", conflict_tolerance=0.01,
            ))

    def test_replicate_aggregation_is_deterministic_and_documented(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", [
            "A5V,0.20", "A5V,0.30", "A5V,0.40", "A9V,0.5",
        ])
        config = IngestConfig(
            design=design, source=path,
            duplicate_policy="mean", conflict_tolerance=0.5,
        )
        first = ingest_dataset(config)
        second = ingest_dataset(config)
        record = {r.variant_id: r for r in first.records}["A5V"]
        self.assertAlmostEqual(record.fitness, 0.30)
        self.assertEqual(record.n_replicates, 3)
        self.assertEqual(record.replicate_values, (0.20, 0.30, 0.40))
        self.assertEqual(first.report.duplicate_aggregation, "mean")
        self.assertEqual(first.fitness_table(), second.fitness_table())

        median = ingest_dataset(IngestConfig(
            design=design, source=path,
            duplicate_policy="median", conflict_tolerance=0.5,
        ))
        self.assertAlmostEqual({r.variant_id: r for r in median.records}["A5V"].fitness, 0.30)

        first_wins = ingest_dataset(IngestConfig(
            design=design, source=path,
            duplicate_policy="first", conflict_tolerance=0.5,
        ))
        self.assertAlmostEqual({r.variant_id: r for r in first_wins.records}["A5V"].fitness, 0.20)

    def test_unknown_duplicate_policy_is_rejected(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2"])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path, duplicate_policy="whatever"))

    # ---- labels, placeholders and censoring --------------------------------------
    def test_unlabelled_rows_with_numeric_placeholders_are_not_measurements(self):
        """The HER2 case: 'test'-stage rows all carry DMS_score exactly 0.0."""
        design = binary_design().with_assay(AssayMeta(detection_floor=0.1))
        path = _write_csv(self.dir, "her2like.csv", "mutant,DMS_score,stage", [
            "A5V,0.5,train",
            "A9V,0.1,train",
            "A5V:A9V,0.0,test",
        ])
        data = ingest_dataset(IngestConfig(
            design=design, source=path,
            labelled_stage="train", unlabelled_placeholder=0.0,
        ))
        self.assertEqual({r.variant_id for r in data.records}, {"A5V", "A9V"})
        self.assertEqual([r.variant_id for r in data.absent_labels], ["A5V:A9V"])
        absent = data.absent_labels[0]
        self.assertIs(absent.fitness, None)
        self.assertEqual(absent.label_status, LabelStatus.ABSENT)
        self.assertEqual(absent.censoring, Censoring.NOT_APPLICABLE)
        self.assertNotIn("A5V:A9V", data.fitness_table())
        self.assertEqual(data.report.rows_absent_label, 1)

    def test_a_real_measurement_in_the_unlabelled_split_is_rejected(self):
        design = binary_design()
        path = _write_csv(self.dir, "her2like.csv", "mutant,DMS_score,stage", [
            "A5V,0.5,train", "A9V,0.77,test",
        ])
        with self.assertRaises(LabelPolicyError):
            ingest_dataset(IngestConfig(
                design=design, source=path,
                labelled_stage="train", unlabelled_placeholder=0.0,
            ))

    def test_missing_stage_column_is_rejected_when_the_policy_needs_it(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.5"])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path, labelled_stage="train"))

    def test_a_true_measured_zero_is_not_an_absent_label(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.0", "A9V,0.4"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        zero = {r.variant_id: r for r in data.records}["A5V"]
        self.assertEqual(zero.fitness, 0.0)
        self.assertEqual(zero.label_status, LabelStatus.MEASURED)
        self.assertEqual(zero.censoring, Censoring.EXACT)
        self.assertEqual(data.report.rows_absent_label, 0)

    def test_values_at_the_assay_floor_are_left_censored_not_exact(self):
        design = binary_design().with_assay(AssayMeta(detection_floor=0.1, units="enrichment"))
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.1", "A9V,0.4"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        by_id = {r.variant_id: r for r in data.records}
        self.assertEqual(by_id["A5V"].censoring, Censoring.LEFT_CENSORED)
        self.assertEqual(by_id["A5V"].fitness, 0.1)  # value preserved, not clamped away
        self.assertEqual(by_id["A9V"].censoring, Censoring.EXACT)
        self.assertEqual(data.report.left_censored, 1)
        self.assertEqual(data.censored_ids(), ("A5V",))

    def test_a_value_below_the_declared_detection_floor_is_rejected(self):
        design = binary_design().with_assay(AssayMeta(detection_floor=0.1))
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.05"])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path))

    def test_values_at_the_saturation_ceiling_are_right_censored(self):
        design = binary_design().with_assay(AssayMeta(detection_ceiling=1.0))
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,1.0", "A9V,0.4"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        by_id = {r.variant_id: r for r in data.records}
        self.assertEqual(by_id["A5V"].censoring, Censoring.RIGHT_CENSORED)
        self.assertEqual(data.report.right_censored, 1)

    def test_non_numeric_and_non_finite_scores_are_rejected(self):
        design = binary_design()
        for bad in ("abc", "nan", "inf", ""):
            with self.subTest(bad=bad):
                path = _write_csv(self.dir, f"bad_{bad or 'empty'}.csv",
                                  "mutant,DMS_score", [f"A5V,{bad}"])
                with self.assertRaises(SchemaError):
                    ingest_dataset(IngestConfig(design=design, source=path))

    # ---- transformations ---------------------------------------------------------
    def test_fitness_transformation_is_applied_and_recorded(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,100.0", "A9V,10.0"])
        data = ingest_dataset(IngestConfig(design=design, source=path, transform="log10"))
        by_id = {r.variant_id: r for r in data.records}
        self.assertAlmostEqual(by_id["A5V"].fitness, 2.0)
        self.assertAlmostEqual(by_id["A5V"].raw_fitness, 100.0)
        self.assertEqual(data.report.transform, "log10")
        self.assertEqual(data.design.assay.transformation, "log10")

    def test_negate_transform_flips_the_declared_fitness_direction(self):
        wt = _wt_sequence(20, p5="A", p9="A")
        design = LibraryDesign(
            target_id="T", assay_id="a", wild_type_sequence=wt,
            positions=binary_design().positions,
            fitness_direction=FitnessDirection.LOWER_IS_BETTER,
            selection_mode=SelectionMode.POOL_RESTRICTED,
        )
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,2.0", "A9V,5.0"])
        data = ingest_dataset(IngestConfig(design=design, source=path, transform="negate"))
        self.assertEqual(data.design.fitness_direction, FitnessDirection.HIGHER_IS_BETTER)
        self.assertAlmostEqual({r.variant_id: r for r in data.records}["A5V"].fitness, -2.0)

    def test_an_impossible_transform_is_rejected_not_silently_nan(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.0"])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path, transform="log10"))
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path, transform="not_a_transform"))

    # ---- structural rejections ---------------------------------------------------
    def test_indel_rows_are_rejected(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5del,0.2"])
        with self.assertRaises(IndelError):
            ingest_dataset(IngestConfig(design=design, source=path))

    def test_sequence_column_length_mismatch_is_an_indel(self):
        design = binary_design()
        short = design.wild_type_sequence[:-1]
        path = _write_csv(self.dir, "bin.csv", "mutant,mutated_sequence,DMS_score",
                          [f"A5V,{short},0.2"])
        with self.assertRaises(IndelError):
            ingest_dataset(IngestConfig(design=design, source=path))

    def test_sequence_column_must_agree_with_the_mutation_string(self):
        design = binary_design()
        wrong = list(design.wild_type_sequence)
        wrong[8] = "V"  # the sequence says A9V, the notation says A5V
        path = _write_csv(self.dir, "bin.csv", "mutant,mutated_sequence,DMS_score",
                          [f"A5V,{''.join(wrong)},0.2"])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path))

    def test_missing_required_columns_are_rejected(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "variant,score", ["A5V,0.2"])
        with self.assertRaises(SchemaError) as ctx:
            ingest_dataset(IngestConfig(design=design, source=path))
        self.assertIn("mutant", str(ctx.exception))
        renamed = ingest_dataset(IngestConfig(
            design=design, source=path,
            column_map=ColumnMap(variant="variant", fitness="score", sequence=None, stage=None),
        ))
        self.assertEqual(len(renamed.records), 1)

    def test_empty_dataset_is_rejected(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", [])
        with self.assertRaises(SchemaError):
            ingest_dataset(IngestConfig(design=design, source=path))

    # ---- modes -------------------------------------------------------------------
    def test_pool_restricted_dataset_exposes_its_enumerated_pool(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A9V,0.3"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(data.selection_mode, SelectionMode.POOL_RESTRICTED)
        self.assertEqual(data.candidate_pool(), ("A5V", "A9V"))
        self.assertEqual(data.candidate_pool(exclude=("A5V",)), ("A9V",))

    def test_generative_dataset_reports_constraints_and_refuses_a_pool(self):
        design = LibraryDesign(
            target_id="T", assay_id="a",
            wild_type_sequence=binary_design().wild_type_sequence,
            positions=binary_design().positions,
            fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
            selection_mode=SelectionMode.GENERATIVE,
        )
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2", "A9V,0.3"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(data.selection_mode, SelectionMode.GENERATIVE)
        with self.assertRaises(SchemaError):
            data.candidate_pool()
        self.assertEqual(data.mutation_constraints()["max_mutation_depth"], 2)

    def test_source_is_fingerprinted(self):
        design = binary_design()
        path = _write_csv(self.dir, "bin.csv", "mutant,DMS_score", ["A5V,0.2"])
        data = ingest_dataset(IngestConfig(design=design, source=path))
        self.assertEqual(len(data.source_sha256), 64)
        self.assertEqual(data.source_size_bytes, path.stat().st_size)
        self.assertEqual(data.report.rows_read, 1)


if __name__ == "__main__":
    unittest.main()
