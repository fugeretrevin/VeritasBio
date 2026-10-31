"""
benchmark_panoptes_isolated.py - Standalone comparison benchmark on a combinatorial library.

Purpose:
    Compares PanoptesBio's production selection path against credible alternatives on a
    registered combinatorial protein library, under a single shared hidden-results protocol.
    Nothing in the existing project is imported for its side effects, and nothing in the
    existing project is modified: this file writes only into a fresh run folder under
    `benchmark_comparison_results/` and imports production modules read-only.

Libraries:
    `LIBRARIES` registers the four-site TrpB active site (Johnston et al. 2024), the
    four-site GB1 domain (Wu et al. 2016) and the ten-site trastuzumab CDR-H3 library
    (Mason et al. 2021). The library is chosen with --library and then THREADED through
    every call path as an argument; there is no module-level active library, so two
    libraries can be loaded and compared in one process.

    A new library is added by registering a `LibrarySpec` - positions, parent residues,
    permissible alphabet, label policy, assay limits, provenance - and pointing it at a
    prepared CSV. No selection or optimisation code changes. The canonical schema and the
    validating ingestion path live in `veritas/datasets/`, which imports nothing from
    `workers` or `controller`.

Methods compared:
    A. random          - uniform sampling of untested eligible variants, without replacement.
    B. ridge           - one-hot (4 positions x 20 AA = 80 features) ridge regression, refit each round.
    C. random_forest   - same encoding, RandomForestRegressor, refit each round.
    D. alde            - the published ALDE implementation (jsunn-y/ALDE). Availability-probed;
                         never silently replaced by a homemade method. Requesting it when the
                         checkout, interpreter or dependencies are absent ABORTS the run (exit 4)
                         rather than printing an UNAVAILABLE row, which once cost two runs.
    E. panoptes_default- the real production planning + selection path with its supported default
                         (non-Muse) strategy.
    F. panoptes_muse   - the same production path with live Muse planning and reflection feedback.
                         Requires explicit --enable-live-muse opt-in.

Protocol (shared by every method):
    - One evaluator owns the complete measured fitness table. Selectors never see it.
    - Per seed: the same 25 starting variants for every method, then 4 rounds of 25 newly
      selected variants. 125 unique measurements per method per seed, including the start.
    - A method commits its batch before the evaluator reveals those outcomes.
    - Duplicate, already-observed, unknown and invalid selections are rejected and logged.
      A short batch is BACKFILLED from the same policy: the selector is asked for a longer
      ranked list and the next eligible picks are taken, up to a bounded number of attempts.
      It is never padded with random draws, and any secondary rule is named in the manifest
      and in the round metrics. A batch that still cannot fill stays short, loudly.

Selection modes (kept explicit and distinct per library):
    POOL-RESTRICTED - the selector is handed the enumerated untested measured variants and
    ranks them. The default for retrospective benchmarking on a library that is a sparse
    sample of its product (HER2: 11,277 measured genotypes against a 20^10 product), where
    free generation proposes variants nobody assayed.
    GENERATIVE - the selector proposes under the library's mutation constraints, which is
    what the saturated libraries (TrpB, GB1) support. Ranking a supplied pool is strictly
    easier, so a pool-mode result never validates unconstrained generation, and every
    pool-mode run says so in its report.

Labels and censoring (never conflated):
    A missing label, a true measured zero and a detection floor are three different things.
    HER2 ships 36,391 rows across a train/test split; the 25,114 'test' rows all carry
    DMS_score exactly 0.0, which is an ABSENT LABEL, not a measurement of no binding - they
    are excluded, and each is checked against the declared placeholder so a real measurement
    cannot be discarded as an absent one. Of the 11,276 eligible HER2 variants, 5,832 sit
    exactly at the assay floor of 0.1; those keep their measured value and are COUNTED as
    left-censored, so a block of ties reads as a detection limit rather than a real plateau.

Provenance (read this before quoting any number from this benchmark):
    The published TrpB four-site landscape is NOT fully measured. Johnston et al. (2024) state
    that 159,129 variants (99.45%) had sufficient sequencing coverage for quantification and
    that "the missing 871 fitness values were imputed" with the scikit-learn KNN imputer.
    The prepared CSV in this repository carries no flag distinguishing the two, so by default
    this benchmark refuses to run on it and tells you how to supply the imputed-variant list.
    See --provenance-mode. Each library carries its own provenance record; nothing quotes
    one library's published split while loading another's.

    Every run resolves its labels onto the explicit `LabelProvenance` vocabulary
    (HISTORICAL_MEASURED / HISTORICAL_IMPUTED / HISTORICAL_UNKNOWN_PROVENANCE /
    SIMULATED_NOISY_OUTCOME / PROSPECTIVE_PHYSICAL_MEASUREMENT / MODEL_PREDICTION) and records
    it in the manifest together with the legacy stamp. An acknowledged-unflagged TrpB run is
    HISTORICAL_UNKNOWN_PROVENANCE, not 'measured with a caveat', and it carries the
    differential-bias limitation verbatim: a shared contaminated universe does NOT make the
    comparison between methods unbiased.

The run manifest is the record:
    manifest.json is authoritative and machine-readable. Its `configuration` subtree states
    what was run - acquisition policy, components enabled AND disabled, dataset identity and
    eligibility rules, model and acquisition hyperparameters, seed block, commit and dirty
    state, environment, budget and schedule, label provenance, assay transformation and
    threshold definition - and carries a sha256 `configuration_fingerprint`. It contains no
    outcome, deliberately: 62 completed arms once had to be told apart by fingerprinting their
    RESULTS, which cannot be allowed to identify a configuration. Outcomes live under
    `results`. The console transcript, the CLI string and report.txt are convenience.

Usage:
    python benchmark_panoptes_isolated.py --help
    python benchmark_panoptes_isolated.py --smoke-test
    python benchmark_panoptes_isolated.py --library gb1 --imputed-ids <file>
"""

from __future__ import annotations

# --- Isolation preamble: must precede every project import. --------------------------------
import sys

sys.dont_write_bytecode = True  # no __pycache__ anywhere, including under the project tree

import argparse
import csv
import hashlib
import inspect
import json
import logging
import os
import platform
import random
import statistics
import subprocess
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "benchmark_comparison_results"

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}

BENCHMARK_VERSION = "1.0.0"

# Bumped whenever the manifest's shape changes. A consumer reads this first.
#   2.0.0 - the manifest became the authoritative, machine-readable record of a run:
#           a self-contained `configuration` subtree with a fingerprint, the explicit
#           label-provenance vocabulary, component ablation state, seed-block identity,
#           code/environment identity, run state and warnings.
MANIFEST_SCHEMA_VERSION = "2.0.0"

# The selector components `--ablate-components` may disable. Mirrors
# workers.math_worker.ABLATABLE_COMPONENTS; run_self_checks asserts the two agree, so a
# name this CLI accepts can never be one the engine rejects half way through a round.
ABLATABLE_COMPONENTS: Tuple[str, ...] = ("dkl", "dpp", "gflownet", "pareto", "prior")

# Process exit codes. Distinct so CI can tell a configuration refusal from a crash.
EXIT_OK = 0
EXIT_ERROR = 1
EXIT_REPOSITORY_CHANGED = 2
EXIT_PROVENANCE = 3
EXIT_BASELINE_UNAVAILABLE = 4

logger = logging.getLogger("panoptes_isolated_benchmark")

# The canonical library/dataset schema lives outside this file, so a new combinatorial
# library is describable without touching the benchmark or the optimiser. Imported
# read-only; nothing in veritas.datasets imports workers or controller.
from veritas.datasets import (  # noqa: E402  (after the isolation preamble, deliberately)
    AssayMeta,
    ColumnMap,
    FitnessDirection,
    IngestConfig,
    LibraryDesign,
    SelectionMode,
    VariablePosition,
)


# ==============================================================================================
# Library registry — the protocol is not TrpB-specific
# ==============================================================================================
@dataclass(frozen=True)
class LibrarySpec:
    """
    A combinatorial library this protocol can run on, described by data alone.

    These were TrpB constants at module scope (`ACTIVE_SITE_POSITIONS`,
    `PARENT_RESIDUES`, `PARENT_ACCESSION`, `PUBLISHED_*`), read directly by the encoders
    and reassigned by a `configure_library()` that mutated nine globals. That made every
    result single-target - the binding limitation on every claim in docs/NEXT_STEPS.md,
    because a method tuned and measured on one landscape demonstrates nothing about
    generalisation - and it made two libraries in one process impossible, which is
    exactly what a cross-protein study needs.

    There is now no active-library state anywhere. Every function that needs a library
    takes one as an argument, and `to_design()` converts this registry entry into the
    canonical `veritas.datasets.LibraryDesign` so a new library is ingested and validated
    through data and configuration without touching the optimiser.
    """

    # ---- REQUIRED --------------------------------------------------------------------
    key: str
    label: str
    dataset: Path
    accession: str
    positions: Tuple[int, ...]
    parent_residues: str
    structure_pdb: Optional[Path]
    published_optimum: Optional[str]
    published_measured_count: Optional[int]
    published_imputed_count: Optional[int]
    provenance_note: str
    # ---- OPTIONAL --------------------------------------------------------------------
    # When set, only CSV rows whose 'stage' column holds this value enter the universe.
    # HER2 ships 36,391 rows of which 25,114 are an unlabelled held-out split carrying
    # DMS_score exactly 0.0 - not a measurement of zero binding but an absent label.
    # Admitting them would put 69% of the universe at an invented value and make every
    # quantile meaningless.
    labelled_stage: Optional[str] = None
    # Saturated libraries (TrpB, GB1) cover ~93-100% of the combinatorial product, so a
    # proposed Hamming neighbour is almost always a variant somebody measured. A
    # pool-restricted library is a sparse sample of its product, so free generation
    # proposes variants that do not exist and the batch comes up short. The selector is
    # handed the untested measured variants instead and ranks them.
    pool_restricted: bool = False
    # The value the unlabelled split carries where a measurement would be. Declared so
    # ingestion can VERIFY that every discarded row really is unlabelled, rather than
    # trusting the split column and silently throwing a measurement away.
    unlabelled_placeholder: Optional[float] = None
    # Per-position permissible alphabets, one string per declared position. None means
    # all 20 at every position. Never inferred from the observed sequences: an assay
    # covering 11,277 of 20^10 genotypes would otherwise "declare" whatever it sampled.
    allowed_aa_per_position: Optional[Tuple[str, ...]] = None
    # Assay metadata. Each default is "unknown", never a plausible guess.
    fitness_direction: str = FitnessDirection.HIGHER_IS_BETTER.value
    assay_units: Optional[str] = None
    assay_transformation: str = "none"
    detection_floor: Optional[float] = None
    detection_ceiling: Optional[float] = None
    replicate_count: Optional[int] = None
    source_publication: Optional[str] = None
    source_doi: Optional[str] = None
    dataset_version: Optional[str] = None
    # Deterministic duplicate handling for ingestion through veritas.datasets.
    duplicate_policy: str = "reject"
    conflict_tolerance: float = 0.0

    # -- derived -----------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        return len(self.positions)

    @property
    def selection_mode(self) -> SelectionMode:
        """
        POOL_RESTRICTED and GENERATIVE are kept explicit and distinct. Ranking an
        enumerated feasible pool is a strictly easier problem than generating under
        mutation constraints, so success in pool mode does not validate generation.
        """
        return SelectionMode.POOL_RESTRICTED if self.pool_restricted else SelectionMode.GENERATIVE

    @property
    def metadata_json(self) -> Path:
        return REPO_ROOT / "workers" / "data" / f"{self.accession}.json"

    def alphabet_for(self, slot: int) -> str:
        if self.allowed_aa_per_position is None:
            return AA_ALPHABET
        return self.allowed_aa_per_position[slot]

    @property
    def one_hot_width(self) -> int:
        return sum(len(self.alphabet_for(s)) for s in range(self.n_sites))

    def assay_meta(self) -> AssayMeta:
        return AssayMeta(
            units=self.assay_units,
            transformation=(
                self.assay_transformation if self.assay_transformation != "none" else None
            ),
            detection_floor=self.detection_floor,
            detection_ceiling=self.detection_ceiling,
            replicate_count=self.replicate_count,
            measurement_provenance=self.provenance_note,
            source_publication=self.source_publication,
            source_doi=self.source_doi,
            dataset_version=self.dataset_version,
        )

    def to_design(self, wild_type_sequence: str) -> LibraryDesign:
        """
        The canonical schema entry for this library.

        Validates the declared positions and parent residues against the real wild-type
        sequence: a mismatch raises `WildTypeMismatchError` rather than aligning the
        dataset against the wrong reference.
        """
        if len(self.parent_residues) != self.n_sites:
            raise BenchmarkError(
                f"library {self.key!r} declares {self.n_sites} positions but "
                f"{len(self.parent_residues)} parent residues ({self.parent_residues!r})"
            )
        if self.allowed_aa_per_position is not None and (
            len(self.allowed_aa_per_position) != self.n_sites
        ):
            raise BenchmarkError(
                f"library {self.key!r} declares {self.n_sites} positions but "
                f"{len(self.allowed_aa_per_position)} per-position alphabets"
            )
        return LibraryDesign(
            target_id=self.accession,
            assay_id=self.key,
            wild_type_sequence=wild_type_sequence,
            positions=tuple(
                VariablePosition(
                    position=position,
                    wild_type_aa=self.parent_residues[slot],
                    allowed_aa=tuple(self.alphabet_for(slot)),
                )
                for slot, position in enumerate(self.positions)
            ),
            fitness_direction=FitnessDirection(self.fitness_direction),
            selection_mode=self.selection_mode,
            assay=self.assay_meta(),
            label=self.label,
        )

    def ingest_config(
        self,
        wild_type_sequence: str,
        dataset_path: Optional[Path] = None,
    ) -> IngestConfig:
        """The configuration that ingests this library through `veritas.datasets`."""
        return IngestConfig(
            design=self.to_design(wild_type_sequence),
            source=dataset_path or self.dataset,
            column_map=ColumnMap(),
            labelled_stage=self.labelled_stage,
            unlabelled_placeholder=self.unlabelled_placeholder,
            duplicate_policy=self.duplicate_policy,
            conflict_tolerance=self.conflict_tolerance,
            transform=self.assay_transformation,
        )


LIBRARIES: Dict[str, LibrarySpec] = {
    "trpb": LibrarySpec(
        key="trpb",
        label="TrpB Tm9D8* four-site active site (Johnston et al. 2024)",
        dataset=REPO_ROOT / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv",
        accession="TRPB_TM9D8",
        positions=(183, 184, 227, 228),
        parent_residues="VFVS",
        structure_pdb=REPO_ROOT / "workers" / "data" / "8EH0.pdb",
        published_optimum="V183A:F184I:V227K:S228G",
        published_measured_count=159_129,
        published_imputed_count=871,
        provenance_note=(
            "Johnston et al. report 159,129 of 160,000 variants quantified and 871 "
            "fitness values imputed with the scikit-learn KNN imputer; the prepared CSV "
            "carries no flag separating them."
        ),
        assay_units="fitness normalised to the landscape maximum",
        source_publication="Johnston et al. 2024",
        source_doi="10.1101/2024.06.23.600144",
    ),
    "gb1": LibrarySpec(
        key="gb1",
        label="GB1 domain four-site library (Wu et al. 2016)",
        dataset=REPO_ROOT / "scripts" / "data" / "GB1_STRSG_Wu_2016.csv",
        accession="SPG1_STRSG",
        positions=(39, 40, 41, 54),
        parent_residues="VDGV",
        structure_pdb=REPO_ROOT / "workers" / "data" / "2QMT.pdb",
        # Not stated in our prepared CSV. Left None rather than guessed: the dataset
        # holds 149,360 of the 160,000 possible combinations, so roughly 6.6% of the
        # library is simply absent rather than imputed, and absent is not the same
        # problem as silently model-filled.
        published_optimum=None,
        published_measured_count=None,
        published_imputed_count=None,
        provenance_note=(
            "This prepared CSV carries 149,360 of 160,000 combinations and no "
            "provenance column. Missing combinations are absent from the eligible "
            "universe rather than imputed into it, but the measured/imputed split of "
            "what IS present is unverified."
        ),
        assay_units="relative fitness (parent = 1.0)",
        source_publication="Wu et al. 2016",
        source_doi="10.7554/eLife.16965",
    ),
    "her2": LibrarySpec(
        key="her2",
        label="Trastuzumab CDR-H3 ten-site library, HER2 binding (Mason et al. 2021)",
        dataset=REPO_ROOT / "scripts" / "data" / "HER2_TRASTUZUMAB_Mason_2021.csv",
        accession="HER2_TRASTUZUMAB",
        positions=(99, 100, 101, 102, 103, 104, 105, 106, 107, 108),
        parent_residues="WGGDGFYAMD",
        structure_pdb=REPO_ROOT / "workers" / "data" / "1N8Z.pdb",
        published_optimum=None,
        published_measured_count=11_277,
        published_imputed_count=0,
        labelled_stage="train",
        # 11,277 labelled variants against a 20^10 product: 0.00000001% saturation.
        # Free generation would propose unmeasured variants essentially every time.
        pool_restricted=True,
        # Every one of the 25,114 held-out rows carries DMS_score exactly 0.0. Declaring
        # the placeholder lets ingestion VERIFY that - a held-out row carrying anything
        # else would be a real measurement about to be discarded, and is refused.
        unlabelled_placeholder=0.0,
        # 5,833 of the 11,277 labelled rows sit exactly at 0.1. That is the assay's
        # detection floor, so those are left-censored observations, not exact values.
        detection_floor=0.1,
        assay_units="normalised HER2 binding enrichment",
        source_publication="Mason et al. 2021",
        source_doi="10.1038/s41551-021-00699-9",
        provenance_note=(
            "The prepared CSV carries 36,391 rows across a 'train'/'test' split; the "
            "25,114 'test' rows all carry DMS_score exactly 0.0, which is an absent "
            "label rather than a measurement of no binding, so only the 11,277 "
            "labelled rows enter the eligible universe. Of those, 5,833 (51.7%) sit at "
            "the assay floor of 0.1 - this library is censored from below, so a large "
            "mass of ties at the floor is expected and is not a defect."
        ),
    ),
}

#: The default library for the CLI's --dataset default only. It is a default value, not
#: process state: nothing reads it to decide which library is running.
DEFAULT_LIBRARY_KEY = "trpb"
DEFAULT_DATASET = LIBRARIES[DEFAULT_LIBRARY_KEY].dataset


def resolve_library(key: str) -> LibrarySpec:
    """
    Looks one library up by key. Replaces the old `configure_library()`.

    `configure_library()` reassigned nine module globals (`ACTIVE_LIBRARY`,
    `ACTIVE_SITE_POSITIONS`, `PARENT_RESIDUES`, ...) that the encoders read at call time.
    That made the library implicit process state: one library per process, no way to
    ingest two in one run, and every encoder silently coupled to whichever call had run
    last. This function only looks the spec up; the caller threads it through explicitly.
    """
    if key not in LIBRARIES:
        raise BenchmarkError(
            f"unknown library {key!r}; choose from {', '.join(sorted(LIBRARIES))}"
        )
    spec = LIBRARIES[key]
    logger.info("Library: %s (%s), positions %s, parent %s, mode %s",
                spec.key, spec.label, list(spec.positions), spec.parent_residues,
                spec.selection_mode.value)
    return spec


# ==============================================================================================
# Errors
# ==============================================================================================
class BenchmarkError(RuntimeError):
    """Base class for actionable benchmark configuration failures."""


class ProvenanceError(BenchmarkError):
    """Raised when measurement provenance cannot be established for the eligible dataset."""


class ConfigurationError(BenchmarkError):
    """
    Raised for a CLI configuration that cannot mean what it appears to mean.

    Specifically: an `--ablate-components` value that names no component, names an empty
    one, or names something the production engine does not recognise. Each of those used
    to run the full unablated engine - or to fail inside the first round - while the run
    folder was labelled as an ablation arm.
    """


class BaselineUnavailableError(BenchmarkError):
    """
    Raised when a method that was explicitly requested cannot run.

    This is deliberately fatal. A requested ALDE baseline was once reported as a single
    'UNAVAILABLE' row in a report because ALDE_PATH was unset, and two complete benchmark
    runs were published from that report before anyone read the row. A comparison that is
    missing a requested arm is not a comparison.
    """


class ManifestIncompleteError(BenchmarkError):
    """Raised when a manifest is missing a field the schema requires."""


# ==============================================================================================
# Label provenance: an explicit, non-interchangeable vocabulary
# ==============================================================================================
class LabelProvenance(Enum):
    """
    Where the number a run optimises against actually came from.

    Deliberately NOT a `str` subclass: these must never compare equal to a bare string
    or to each other, because the whole failure mode being fixed is one opaque stamp
    ('UNVERIFIED_INCLUDES_IMPUTED') standing in for six materially different claims.
    Serialise with `.value`; parse with `LabelProvenance(value)`.
    """

    HISTORICAL_MEASURED = "HISTORICAL_MEASURED"
    HISTORICAL_IMPUTED = "HISTORICAL_IMPUTED"
    HISTORICAL_UNKNOWN_PROVENANCE = "HISTORICAL_UNKNOWN_PROVENANCE"
    SIMULATED_NOISY_OUTCOME = "SIMULATED_NOISY_OUTCOME"
    PROSPECTIVE_PHYSICAL_MEASUREMENT = "PROSPECTIVE_PHYSICAL_MEASUREMENT"
    MODEL_PREDICTION = "MODEL_PREDICTION"


LABEL_PROVENANCE_DEFINITIONS: Dict[LabelProvenance, str] = {
    LabelProvenance.HISTORICAL_MEASURED: (
        "A value physically measured in a previously published or previously run experiment, "
        "and identified as measured by the dataset itself (a provenance column) or by an "
        "external list of the records that were not measured."
    ),
    LabelProvenance.HISTORICAL_IMPUTED: (
        "A value filled in by a model from neighbouring records in a historical dataset - for "
        "the TrpB landscape, the scikit-learn KNN imputer the authors used. Not a measurement."
    ),
    LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE: (
        "A historical dataset whose records cannot be separated into measured and imputed. "
        "The universe is known to contain both and no per-record flag exists, so no individual "
        "label can be claimed as measured."
    ),
    LabelProvenance.SIMULATED_NOISY_OUTCOME: (
        "A generated value from a simulator or synthetic landscape. Useful for validating "
        "mechanism; says nothing about biological performance."
    ),
    LabelProvenance.PROSPECTIVE_PHYSICAL_MEASUREMENT: (
        "A value measured physically after this system selected the variant - the only kind "
        "that supports a prospective claim. No run in this benchmark produces one."
    ),
    LabelProvenance.MODEL_PREDICTION: (
        "A model's own prediction (for example an ESM pseudo-likelihood) used in place of an "
        "outcome. Never comparable with a measurement."
    ),
}

# The legacy strings this benchmark has written into traces, reports and plots. They are
# KEPT alongside the new vocabulary - existing run folders and the greps in docs/545.md
# still carry them - and mapped forwards rather than reinterpreted.
LEGACY_PROVENANCE_STATUS_TO_ENUM: Dict[str, LabelProvenance] = {
    "VERIFIED_COLUMN": LabelProvenance.HISTORICAL_MEASURED,
    "VERIFIED_EXTERNAL_LIST": LabelProvenance.HISTORICAL_MEASURED,
    "UNVERIFIED_INCLUDES_IMPUTED": LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE,
    "SYNTHETIC_NOT_BIOLOGICAL": LabelProvenance.SIMULATED_NOISY_OUTCOME,
}

# Why HISTORICAL_UNKNOWN_PROVENANCE is not a bookkeeping detail. This text goes into the
# manifest verbatim, so no run can be read without it.
DIFFERENTIAL_BIAS_WARNING = (
    "Do NOT assume that because every method searches the same contaminated universe the "
    "comparison between methods is unbiased. The contamination is not independent of the "
    "methods: the published KNN-imputed values sit wherever sequencing coverage was poorest, "
    "which is not uniform over the landscape. If imputed records cluster in the region one "
    "method prefers to explore, that method is scored against model output while another is "
    "scored against measurements, and the bias is DIFFERENTIAL, not common-mode. Nothing in "
    "this run bounds that effect, and no shared-universe argument removes it."
)

# Kept until both are closed; the manifest names them so a reader of any single run folder
# learns the state of the open work, not just the stamp.
PROVENANCE_OPEN_BLOCKERS: Tuple[str, ...] = (
    "TrpB: the published imputed-variant identifiers (Johnston et al. 2024, SI Appendix) have "
    "not been obtained, so the 871 KNN-imputed records cannot be excluded or even located. "
    "Supply them with --imputed-ids to move this run to HISTORICAL_MEASURED.",
    "GB1: label provenance is unresolved - the held GB1 dataset has not been checked for "
    "imputed or inferred labels, so it must not be run under a measured-provenance claim.",
)


def map_legacy_provenance_status(status: str) -> LabelProvenance:
    """
    Maps a legacy provenance string onto the explicit vocabulary.

    Unknown strings raise: guessing is exactly how one stamp came to stand for six claims.
    """
    try:
        return LEGACY_PROVENANCE_STATUS_TO_ENUM[status]
    except KeyError:
        raise ProvenanceError(
            f"no mapping from legacy provenance status '{status}' onto the explicit "
            f"vocabulary {[m.value for m in LabelProvenance]}. Add one deliberately in "
            "LEGACY_PROVENANCE_STATUS_TO_ENUM; it must not be inferred."
        ) from None


def provenance_stamp(universe: "Universe") -> str:
    """
    The one-line stamp written into every CSV, plot and report.

    Carries BOTH vocabularies: the new enum leads, and the legacy string is retained so
    existing run folders and documented greps keep matching.
    """
    resolved = map_legacy_provenance_status(universe.provenance_status)
    stamp = f"provenance={resolved.value} legacy={universe.provenance_status}"
    if universe.synthetic:
        stamp += " SYNTHETIC_NOT_BIOLOGICAL-mechanism-validation-only"
    return stamp


# ==============================================================================================
# Dataset loading, provenance and the eligible universe
# ==============================================================================================
@dataclass(frozen=True)
class DatasetFingerprint:
    """Identity of the dataset file the run was scored against."""

    path: str
    sha256: str
    size_bytes: int
    rows_parsed: int
    unique_variants: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "rows_parsed": self.rows_parsed,
            "unique_variants": self.unique_variants,
        }


@dataclass
class Universe:
    """
    The eligible candidate universe and its evaluator-only truth table.

    `fitness` is deliberately the ONLY place complete outcomes live. Selectors are handed
    `variant_ids` and `one_hot` (both sequence-derived) plus their own revealed measurements.
    """

    variant_ids: List[str]                      # canonical, position-ordered mutation notation
    combos: Dict[str, str]                      # variant id -> one residue per variable site
    fitness: Dict[str, float]                   # EVALUATOR ONLY
    one_hot: np.ndarray                         # (n_variants, n_sites * 20), rows == variant_ids
    index_of: Dict[str, int]
    # The library this universe was built from. Carried on the object rather than read
    # from module state, so a universe is self-describing and two can coexist.
    library: LibrarySpec
    fingerprint: DatasetFingerprint
    provenance_status: str
    provenance_detail: str
    excluded_imputed: int
    duplicates_dropped: int
    #: Rows that exist in the file but carry no label (HER2's held-out split). Counted so
    #: a reader can see the eligible universe is the labelled split, not the whole file.
    unlabelled_excluded: int = 0
    #: Measurements sitting at the library's declared detection limits. Their values are
    #: unchanged - what is recorded is that they are censored, not exact.
    censored_at_floor: int = 0
    censored_at_ceiling: int = 0
    synthetic: bool = False
    # Ground truth for a planted synthetic interaction, for the TEST HARNESS only.
    # Never read by a selector: a policy that could see this would be scored on a
    # landscape whose answer it had been handed, which is the failure mode the
    # synthetic landscapes exist to avoid.
    planted_interaction: Optional[Dict[str, Any]] = None
    # The two remaining silent drops, counted rather than only logged, because 'any
    # excluded records' is part of what the manifest has to state: a drop changes the
    # eligible universe, and the eligible universe defines the top-1% threshold.
    # `unmapped_dropped` is only reachable under --lenient-rows, where a record that
    # does not match the declared design is counted instead of refused; it must not
    # then vanish. `malformed_rows_dropped` covers blank and non-finite scores, which
    # were skipped with a bare `continue`.
    unmapped_dropped: int = 0
    malformed_rows_dropped: int = 0

    def __len__(self) -> int:
        return len(self.variant_ids)


def _sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _canonical_variant_id(combo: str, library: LibrarySpec) -> str:
    """
    Builds position-ordered ProteinGym notation for a FULL-LENGTH combo over
    `library.positions`, dropping the sites that carry the parent residue.

    The combo has one residue per declared position - 4 characters for TrpB and GB1, 10
    for the HER2 CDR-H3 library - and that is enforced, because a short combo would
    otherwise be zipped against the position list and silently produce a different
    variant id.
    """
    if len(combo) != library.n_sites:
        raise BenchmarkError(
            f"genotype {combo!r} has {len(combo)} site(s) but library {library.key!r} "
            f"declares {library.n_sites} ({list(library.positions)}); genotypes are never "
            f"truncated or padded to match"
        )
    parts = [
        f"{library.parent_residues[i]}{pos}{combo[i]}"
        for i, pos in enumerate(library.positions)
        if combo[i] != library.parent_residues[i]
    ]
    return ":".join(parts)


def _combo_from_notation(notation: str, library: LibrarySpec) -> Optional[str]:
    """
    Maps 'V183A:F184I' to the FULL-LENGTH combo 'AIVS', one residue per position in
    `library.positions`.

    A mutation string OMITS wild-type positions, so 'V183A:F184A:V227A' is a four-site
    genotype whose fourth site is the parent residue; the returned combo is always
    `library.n_sites` characters, with the parent filled in wherever the string is silent.

    Returns None - rather than raising - when the notation touches a position outside the
    library, disagrees with the parent residue, or is not substitution notation. The
    engine's own proposals are passed through here, and an off-library proposal must be
    recorded as a rejected selection in the trace rather than abort the campaign. The
    strict, raising equivalent is `veritas.datasets.LibraryDesign.combo_from_mutations`,
    which is what dataset ingestion uses.
    """
    residues = list(library.parent_residues)
    if notation in ("", "WT"):
        return "".join(residues)
    for token in notation.split(":"):
        token = token.strip()
        if len(token) < 3:
            return None
        wt_aa, mut_aa, pos_text = token[0], token[-1], token[1:-1]
        if not pos_text.isdigit():
            return None
        position = int(pos_text)
        if position not in library.positions:
            return None
        slot = library.positions.index(position)
        if wt_aa != library.parent_residues[slot] or mut_aa not in AA_INDEX:
            return None
        residues[slot] = mut_aa
    return "".join(residues)


def _one_hot_matrix(combos: Sequence[str], library: LibrarySpec) -> np.ndarray:
    """One-hot encodes combos as len(library.positions) x 20 binary features
    (80 for the four-site libraries, 200 for the ten-site HER2 library)."""
    matrix = np.zeros((len(combos), library.n_sites * len(AA_ALPHABET)), dtype=np.float32)
    for row, combo in enumerate(combos):
        if len(combo) != library.n_sites:
            raise BenchmarkError(
                f"genotype {combo!r} at row {row} has {len(combo)} site(s) but library "
                f"{library.key!r} declares {library.n_sites}"
            )
        for slot, residue in enumerate(combo):
            matrix[row, slot * len(AA_ALPHABET) + AA_INDEX[residue]] = 1.0
    return matrix


def _load_imputed_ids(
    dataset_path: Path,
    explicit_path: Optional[Path],
) -> Tuple[Optional[Set[str]], str]:
    """
    Resolves the list of imputed (non-measured) variants, in order of preference:
      1. --imputed-ids FILE               (one variant id or 4-residue combo per line, or a CSV column)
      2. <dataset>.provenance.json        ({"imputed_variants": [...]})
    Returns (ids, source-description), or (None, reason) when no list could be resolved.
    """
    def _parse_lines(path: Path) -> Set[str]:
        ids: Set[str] = set()
        with open(path, "r", encoding="utf-8") as handle:
            for raw in handle:
                token = raw.strip().split(",")[0].strip().strip('"')
                if not token or token.lower() in ("mutant", "variant", "variant_id", "combo"):
                    continue
                ids.add(token)
        return ids

    if explicit_path is not None:
        if not explicit_path.is_file():
            raise ProvenanceError(f"--imputed-ids file not found: {explicit_path}")
        if explicit_path.suffix.lower() == ".json":
            payload = json.loads(explicit_path.read_text(encoding="utf-8"))
            ids = set(payload.get("imputed_variants") or payload.get("imputed") or [])
        else:
            ids = _parse_lines(explicit_path)
        return ids, f"--imputed-ids {explicit_path}"

    sidecar = dataset_path.with_suffix(dataset_path.suffix + ".provenance.json")
    if sidecar.is_file():
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
        ids = set(payload.get("imputed_variants") or payload.get("imputed") or [])
        return ids, f"sidecar {sidecar.name}"

    return None, "no imputed-variant list supplied and the dataset carries no provenance column"


def _provenance_help(library: LibrarySpec) -> str:
    """
    The actionable refusal text, written for whichever library is being loaded.

    Library-specific because the numbers are: TrpB has a published measured/imputed split,
    GB1 has none recorded, and quoting TrpB's 871 imputed values while loading GB1 would
    be a fabricated provenance claim.
    """
    if library.published_measured_count is not None and library.published_imputed_count:
        why = f"""
  The upstream authors report that {library.published_measured_count:,} variants had
  sufficient sequencing coverage to be quantified, and that {library.published_imputed_count}
  fitness values were IMPUTED. The prepared CSV in this repository re-encodes the published
  label column verbatim and carries no flag separating measured from imputed records, so a
  run over it would mix {library.published_imputed_count} model-filled values into a
  benchmark presented as experimental."""
    else:
        why = f"""
  {library.provenance_note}"""
    return f"""
Measurement provenance could not be established for library '{library.key}'
({library.dataset.name}).

Why this matters:{why}

How to fix it (pick one):
  1. Supply the imputed-variant list:
       --imputed-ids <file>
     The file may be one variant per line (either 'V183A:F184I:V227K:S228G' notation or a
     4-residue combo such as 'AIKG'), a single-column CSV, or a JSON object shaped
     {{"imputed_variants": [...]}}. The published list is derivable from the authors'
     supplementary material (SI Appendix, Fig. S16) or the upstream per-replicate read counts.
  2. Drop a sidecar next to the dataset named '<dataset>.csv.provenance.json' with the same
     JSON shape; it is picked up automatically.
  3. Add a provenance column to a copy of the dataset: one of 'provenance', 'is_imputed',
     'imputed' or 'measured', and point --dataset at the copy.
  4. Acknowledge the limitation explicitly and proceed anyway:
       --provenance-mode acknowledge-unflagged
     Every output file, the plot and the report are then stamped
     'provenance=HISTORICAL_UNKNOWN_PROVENANCE legacy=UNVERIFIED_INCLUDES_IMPUTED', and the
     manifest records the differential-bias limitation in full. Results from such a run must
     not be described as a clean experimentally-grounded comparison, and must not be defended
     on the grounds that every method saw the same contaminated universe.
"""


def load_universe(
    library: LibrarySpec,
    dataset_path: Path,
    provenance_mode: str,
    imputed_ids_path: Optional[Path],
    max_variants: Optional[int] = None,
    subsample_seed: int = 0,
    strict_rows: bool = True,
) -> Universe:
    """
    Loads one library's prepared CSV into an eligible universe, enforcing the provenance
    and label policies declared on its `LibrarySpec`.

    `library` is an explicit argument, not module state: that is what makes this loader
    reusable for a library it has never seen, and what lets a cross-protein study hold
    two universes in one process.

    Filtering rules, all reported in the manifest:
      - a record must map onto the declared positions and agree with the parent residues.
        With `strict_rows` (the default) a record that does not is a REFUSAL naming the
        row, because a dataset whose rows do not match the declared design is a
        configuration error, not noise to be dropped;
      - rows outside `library.labelled_stage` carry no label and never enter the universe.
        When the library declares `unlabelled_placeholder`, each such row is checked
        against it, so a real measurement cannot be discarded as an absent label;
      - duplicate variant ids keep the first occurrence and are counted;
      - non-finite scores are dropped;
      - imputed records are dropped when a provenance list is available.

    Censoring is RECORDED, never applied: a value at the declared detection floor stays
    at its measured value and is counted in `censored_at_floor`, so a downstream reader
    can see that a block of ties is a detection limit rather than a real plateau.
    """
    if not dataset_path.is_file():
        raise BenchmarkError(
            f"Dataset not found: {dataset_path}\n"
            f"Prepare it first with:  python scripts/download_trpb_benchmark.py"
        )

    sha, size = _sha256_file(dataset_path)

    rows: List[Tuple[str, float, Optional[str]]] = []
    malformed_rows = 0
    provenance_column: Optional[str] = None
    unlabelled_skipped = 0
    required_stage = library.labelled_stage
    with open(dataset_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]
        if required_stage is not None and "stage" not in fieldnames:
            raise BenchmarkError(
                f"library {library.key!r} requires rows with stage "
                f"{required_stage!r}, but {dataset_path} has no 'stage' column "
                f"(found {fieldnames}). Refusing to admit unlabelled rows as "
                f"measurements."
            )
        if "mutant" not in fieldnames or "DMS_score" not in fieldnames:
            raise BenchmarkError(
                f"{dataset_path} does not look like a prepared benchmark CSV for library "
                f"{library.key!r} (expected 'mutant' and 'DMS_score' columns, found "
                f"{fieldnames})."
            )
        for candidate in ("provenance", "is_imputed", "imputed", "measured"):
            if candidate in fieldnames:
                provenance_column = candidate
                break
        for offset, row in enumerate(reader):
            line = offset + 2  # 1-indexed, header included
            if required_stage is not None and (row.get("stage") or "").strip() != required_stage:
                # An absent label, not a measurement of zero. When the library declares
                # what the placeholder is, verify it: a held-out row carrying anything
                # else is a real measurement about to be thrown away.
                if library.unlabelled_placeholder is not None:
                    placeholder_text = (row.get("DMS_score") or "").strip()
                    try:
                        observed = float(placeholder_text)
                    except ValueError:
                        observed = float("nan")
                    if not (abs(observed - library.unlabelled_placeholder) < 1e-9):
                        raise BenchmarkError(
                            f"{dataset_path.name} line {line}: row is outside the labelled "
                            f"stage {required_stage!r} but carries DMS_score "
                            f"{placeholder_text!r}, not the declared placeholder "
                            f"{library.unlabelled_placeholder!r}. That looks like a real "
                            f"measurement in the held-out split; refusing to discard it as "
                            f"an absent label."
                        )
                unlabelled_skipped += 1
                continue
            notation = (row.get("mutant") or "").strip()
            raw_score = (row.get("DMS_score") or "").strip()
            if not notation or not raw_score:
                malformed_rows += 1
                continue
            try:
                score = float(raw_score)
            except ValueError:
                malformed_rows += 1
                continue
            if not np.isfinite(score):
                malformed_rows += 1
                continue
            flag = (row.get(provenance_column) or "").strip() if provenance_column else None
            rows.append((notation, score, flag))

    if not rows:
        raise BenchmarkError(f"No usable rows parsed from {dataset_path}")

    # ---- Provenance resolution -------------------------------------------------------------
    imputed_lookup: Optional[Set[str]] = None
    provenance_detail: str
    if provenance_column:
        provenance_status = "VERIFIED_COLUMN"
        provenance_detail = f"dataset column '{provenance_column}'"
    else:
        imputed_lookup, source = _load_imputed_ids(dataset_path, imputed_ids_path)
        if imputed_lookup is not None:
            provenance_status = "VERIFIED_EXTERNAL_LIST"
            provenance_detail = f"{len(imputed_lookup)} imputed ids from {source}"
        elif provenance_mode == "acknowledge-unflagged":
            provenance_status = "UNVERIFIED_INCLUDES_IMPUTED"
            # The imputed count is library-specific and may be unknown, so fall back to
            # the library's own provenance note rather than printing "None".
            if library.published_imputed_count is not None:
                caveat = (
                    f"up to {library.published_imputed_count} published KNN-imputed values "
                    f"are included"
                )
            else:
                caveat = f"measured/imputed split unverified ({library.provenance_note})"
            provenance_detail = (
                f"{source}; {caveat} because --provenance-mode acknowledge-unflagged was passed"
            )
        else:
            raise ProvenanceError(_provenance_help(library).strip())

    # ---- Build the universe ----------------------------------------------------------------
    seen: Set[str] = set()
    duplicates = 0
    excluded_imputed = 0
    unmapped = 0
    censored_at_floor = 0
    censored_at_ceiling = 0
    ids: List[str] = []
    combos: List[str] = []
    scores: List[float] = []
    floor = library.detection_floor
    ceiling = library.detection_ceiling

    for notation, score, flag in rows:
        combo = _combo_from_notation(notation, library)
        if combo is None:
            if strict_rows:
                raise BenchmarkError(
                    f"{dataset_path.name}: record {notation!r} does not map onto library "
                    f"{library.key!r}'s declared positions {list(library.positions)} with "
                    f"parent residues {library.parent_residues!r}. Either the dataset or the "
                    f"library design is wrong; pass --lenient-rows to drop such records and "
                    f"count them instead."
                )
            unmapped += 1
            continue
        variant_id = _canonical_variant_id(combo, library)
        if not variant_id:
            continue  # the parent itself is not a library member
        if variant_id in seen:
            duplicates += 1
            continue
        if provenance_column:
            is_imputed = flag is not None and (
                flag.lower() in ("1", "true", "yes", "imputed")
                if provenance_column != "measured"
                else flag.lower() in ("0", "false", "no")
            )
            if is_imputed:
                excluded_imputed += 1
                continue
        elif imputed_lookup is not None and (variant_id in imputed_lookup or combo in imputed_lookup):
            excluded_imputed += 1
            continue
        seen.add(variant_id)
        ids.append(variant_id)
        combos.append(combo)
        scores.append(score)
        # Censoring is COUNTED, never applied. A value at the declared floor keeps its
        # measured value; what changes is that the report can say the ties are a
        # detection limit rather than a real plateau.
        if floor is not None and score <= floor + 1e-12:
            censored_at_floor += 1
        elif ceiling is not None and score >= ceiling - 1e-12:
            censored_at_ceiling += 1

    if unlabelled_skipped:
        logger.warning(
            "%d record(s) outside stage %r were dropped as unlabelled; the eligible "
            "universe is the labelled split only.", unlabelled_skipped, required_stage
        )
    if unmapped:
        logger.warning("%d record(s) did not map onto positions %s and were dropped.",
                       unmapped, list(library.positions))
    if censored_at_floor:
        logger.info(
            "%d of %d eligible variants (%.1f%%) sit at the declared detection floor %s "
            "and are LEFT-CENSORED, not exact measurements.",
            censored_at_floor, len(ids), 100.0 * censored_at_floor / max(len(ids), 1), floor,
        )
    if len(ids) < 1000:
        raise BenchmarkError(
            f"Only {len(ids)} eligible variants survived filtering; refusing to benchmark on so few records."
        )

    if max_variants is not None and max_variants < len(ids):
        rng = random.Random(subsample_seed)
        keep = sorted(rng.sample(range(len(ids)), max_variants))
        ids = [ids[i] for i in keep]
        combos = [combos[i] for i in keep]
        scores = [scores[i] for i in keep]
        provenance_detail += f"; candidate universe subsampled to {max_variants} variants (seed {subsample_seed})"

    order = sorted(range(len(ids)), key=lambda i: ids[i])
    ids = [ids[i] for i in order]
    combos = [combos[i] for i in order]
    scores = [scores[i] for i in order]

    universe = Universe(
        variant_ids=ids,
        combos=dict(zip(ids, combos)),
        fitness=dict(zip(ids, scores)),
        one_hot=_one_hot_matrix(combos, library),
        index_of={v: i for i, v in enumerate(ids)},
        library=library,
        fingerprint=DatasetFingerprint(
            path=str(dataset_path),
            sha256=sha,
            size_bytes=size,
            rows_parsed=len(rows),
            unique_variants=len(ids),
        ),
        provenance_status=provenance_status,
        provenance_detail=provenance_detail,
        excluded_imputed=excluded_imputed,
        duplicates_dropped=duplicates,
        unlabelled_excluded=unlabelled_skipped,
        censored_at_floor=censored_at_floor,
        censored_at_ceiling=censored_at_ceiling,
        unmapped_dropped=unmapped,
        malformed_rows_dropped=malformed_rows,
    )

    best_id = max(universe.fitness, key=lambda v: universe.fitness[v])
    if library.published_optimum is None:
        # No published optimum recorded for this library, so there is nothing to check
        # against. Log the observed best so a fitness-direction error is still visible.
        logger.info("Global optimum in the eligible set: %s (%.4f); no published value "
                    "on record for library %r to cross-check against.",
                    best_id, universe.fitness[best_id], library.key)
    elif best_id != library.published_optimum:
        logger.warning(
            "Global optimum in the eligible set is %s, not the published %s. Check the dataset and "
            "the fitness direction before interpreting results.",
            best_id, library.published_optimum,
        )
    return universe


def build_synthetic_universe(
    library: LibrarySpec,
    n_variants: int = 0,
    seed: int = 7,
    plant_sign_epistasis: bool = False,
    solo_penalty: float = 0.60,
    pair_bonus: float = 1.40,
) -> Universe:
    """
    Builds a clearly-labelled SYNTHETIC landscape over `library`'s design, for mechanism
    validation.

    Values here are generated additive-plus-pairwise numbers, NOT biology. Every output of a
    synthetic run is stamped so these numbers can never be mistaken for assay performance.

    `n_variants <= 0` (the default) covers the complete 20^n_sites - 1 combo space, which is what
    the production adapter needs: it proposes variants from the full saturated library, so a sparse
    synthetic subsample would reject nearly every one of its selections and leave that code path
    untested. A positive value subsamples instead, for a faster but less representative check.
    """
    import itertools

    n_sites = library.n_sites
    rng = np.random.default_rng(seed)
    singles = rng.normal(0.0, 0.25, size=(n_sites, len(AA_ALPHABET)))
    pair = rng.normal(0.0, 0.12, size=(n_sites, len(AA_ALPHABET), n_sites, len(AA_ALPHABET)))

    # Enumerating the complete product is only tractable for a small library: it is
    # 160,000 genotypes for a four-site one and 20^10 for the HER2 CDR-H3 library. A
    # library that large has to be subsampled explicitly rather than silently.
    full_size = len(AA_ALPHABET) ** n_sites
    if not n_variants and full_size > 1_000_000:
        raise BenchmarkError(
            f"the complete product of library {library.key!r} is {len(AA_ALPHABET)}^{n_sites} "
            f"genotypes, which is too large to enumerate. Pass --synthetic-size to subsample "
            f"a synthetic landscape over this design."
        )
    complete = [
        "".join(parts)
        for parts in itertools.product(AA_ALPHABET, repeat=n_sites)
        if "".join(parts) != library.parent_residues
    ]
    if n_variants and n_variants < len(complete):
        combos: Sequence[str] = random.Random(seed).sample(complete, n_variants)
    else:
        combos = complete

    # Optionally plant a SIGN-EPISTATIC pair with a known, hidden structure.
    #
    # Why this exists. The engine's measured weakness is conditional interactions: on
    # TrpB, K227 has mean marginal effect -0.021 and is better than parent in only 39%
    # of backgrounds, yet 100% of variants above 0.90 fitness carry it, because its
    # benefit is conditional on G228. Any model fitted on marginals correctly learns
    # "K227 is bad" and stops sampling it.
    #
    # Developing a fix on TrpB is unfalsifiable: we already know the answer, so a method
    # tuned until it recovers that pair proves nothing about a landscape we have not
    # seen. This plants an equivalent structure at a position/residue pair chosen by the
    # generator's own RNG, so a policy can be tested on landscapes whose interaction it
    # has never been told about. `planted` records the ground truth for the test
    # harness; it is never placed in the Universe handed to a selector.
    planted: Optional[Dict[str, Any]] = None
    if plant_sign_epistasis:
        slot_a, slot_b = rng.choice(n_sites, size=2, replace=False)
        res_a = AA_ALPHABET[int(rng.integers(len(AA_ALPHABET)))]
        res_b = AA_ALPHABET[int(rng.integers(len(AA_ALPHABET)))]
        planted = {
            "slot_a": int(slot_a), "residue_a": res_a,
            "slot_b": int(slot_b), "residue_b": res_b,
            "position_a": library.positions[int(slot_a)],
            "position_b": library.positions[int(slot_b)],
            "solo_penalty": float(solo_penalty),
            "pair_bonus": float(pair_bonus),
        }

    ids: List[str] = []
    combo_list: List[str] = []
    scores: List[float] = []
    for combo in sorted(combos):
        value = 0.0
        for i, residue_i in enumerate(combo):
            value += singles[i, AA_INDEX[residue_i]]
            for j, residue_j in enumerate(combo):
                if j > i:
                    value += pair[i, AA_INDEX[residue_i], j, AA_INDEX[residue_j]]
        if planted is not None:
            has_a = combo[planted["slot_a"]] == planted["residue_a"]
            has_b = combo[planted["slot_b"]] == planted["residue_b"]
            # Each residue alone is deleterious; together they are strongly beneficial.
            # This is the shape of the TrpB K227/G228 effect, which makes the marginal
            # effect of each one negative while the combination carries the optimum.
            if has_a and has_b:
                value += pair_bonus
            elif has_a or has_b:
                value -= solo_penalty
        ids.append(_canonical_variant_id(combo, library))
        combo_list.append(combo)
        scores.append(float(value))

    order = sorted(range(len(ids)), key=lambda i: ids[i])
    ids = [ids[i] for i in order]
    combo_list = [combo_list[i] for i in order]
    scores = [scores[i] for i in order]

    return Universe(
        variant_ids=ids,
        combos=dict(zip(ids, combo_list)),
        fitness=dict(zip(ids, scores)),
        one_hot=_one_hot_matrix(combo_list, library),
        index_of={v: i for i, v in enumerate(ids)},
        library=library,
        fingerprint=DatasetFingerprint(
            path=f"<synthetic:seed={seed}>",
            sha256=hashlib.sha256(f"synthetic-{seed}-{n_variants}".encode()).hexdigest(),
            size_bytes=0,
            rows_parsed=len(ids),
            unique_variants=len(ids),
        ),
        provenance_status="SYNTHETIC_NOT_BIOLOGICAL",
        provenance_detail=(
            "generated additive + pairwise values for mechanism validation only"
            + (
                f"; PLANTED sign epistasis at "
                f"{planted['residue_a']}{planted['position_a']} x "
                f"{planted['residue_b']}{planted['position_b']} "
                f"(each alone -{solo_penalty}, together +{pair_bonus})"
                if planted else ""
            )
        ),
        planted_interaction=planted,
        excluded_imputed=0,
        duplicates_dropped=0,
        synthetic=True,
    )


# ==============================================================================================
# Evaluator and per-method observation state
# ==============================================================================================
#: The qualification bar's default, expressed as a QUANTILE of the eligible universe
#: (0-100) rather than as an absolute fitness.
#:
#: It used to be the absolute fitness 0.70. On TrpB, the landscape it was tuned on, that
#: is the 99.9602nd percentile of the eligible universe: 64 of 159,999 variants, with
#: P(a random 25-variant start already qualifies) = 0.0100. It was an extreme-tail gate.
#:
#: Fitness scales are not comparable across this corpus, so the same constant admits
#: 31.6% of GB1, 9.11% of HER2 (1,027 of 11,276 - which is why every HER2 seed logged
#: "qualified before round 1"), and 100% of CreiLOV, CR9114 H1 and both CR6261
#: landscapes, where 0.70 lies below the global minimum because the bnAb assay floor is
#: 7.0 on a -log10(Kd) scale. On 8 of 9 landscapes the bar sat BELOW that target's own
#: top-1% threshold, so "qualifying" was an easier bar than being a hit, which inverts
#: the stage's purpose.
#:
#: 99.96 is chosen to reproduce 0.70's behaviour on the landscape the parameter was tuned
#: on: it resolves to 0.6998 on TrpB and accepts the IDENTICAL 64 variants, so the
#: development result does not move. Resolved absolutes elsewhere: GB1 5.362, HER2 6.379,
#: eqFP611 blue 1.588, eqFP611 red 1.521, CreiLOV 4.166, CR9114 H1 9.687, CR6261 H1
#: 9.655, CR6261 H9 10.150.
#:
#: KNOWN CAVEAT, deliberately unresolved: at this quantile the qualifying set on CR6261
#: H1/H9 is ONE variant, so qualification there becomes "find the global best". A
#: quantile transfers the RARITY of the bar, not its ATTAINABILITY. No minimum-qualifying-
#: set floor is applied; instead every run records `qualify_qualifying_count` in its
#: manifest so a reader sees the attainability of the bar it actually ran under. A floor
#: would have to be chosen by the protocol in advance, not by this module.
DEFAULT_QUALIFY_QUANTILE = 99.96

#: The pre-fix absolute. Kept so a run can reproduce the old bar deliberately via
#: --qualify-threshold, and so the manifest can name what changed.
LEGACY_QUALIFY_ABSOLUTE = 0.70

#: How many extra candidate requests a round may make to fill a short batch.
#:
#: Bounded on purpose. Each attempt asks the policy for a longer ranked list, and the
#: margin grows geometrically, so four attempts cover a rejection rate far above the
#: ~6.65% that GB1's incomplete coverage produces. A batch that is still short after this
#: is reported as short, loudly, rather than retried forever.
BACKFILL_MAX_ATTEMPTS = 4


@dataclass
class RevealedState:
    """A single method's private view: only what the evaluator has already revealed to it."""

    measured: Dict[str, float] = field(default_factory=dict)

    @property
    def ids(self) -> List[str]:
        return list(self.measured.keys())

    def best(self) -> float:
        return max(self.measured.values()) if self.measured else float("nan")


class Evaluator:
    """
    Owns the complete measured fitness table and the top-1% definition.

    Selection methods receive this object's `reveal()` output only after they have committed a
    batch. The threshold and top-1% membership stay evaluator-side and are never passed to a
    selector, logged into selector-visible structures, or exposed through candidate metadata.
    """

    def __init__(
        self,
        universe: Universe,
        *,
        qualify_quantile: Optional[float] = None,
        qualify_threshold: Optional[float] = None,
    ) -> None:
        if qualify_quantile is not None and qualify_threshold is not None:
            raise BenchmarkError(
                "qualify_quantile and qualify_threshold are mutually exclusive: the first "
                "states the bar as a quantile of the eligible universe, the second as an "
                "absolute fitness. Passing both leaves it ambiguous which one gated the run."
            )
        self._universe = universe
        self.harness_exports: List[str] = []
        values = np.array([universe.fitness[v] for v in universe.variant_ids], dtype=np.float64)
        # Defined once over the eligible measured dataset. Ties: every variant at or above the
        # threshold counts as a hit, so the hit set can exceed 1% of the pool when values repeat.
        self.top1_threshold = float(np.percentile(values, 99.0))
        self._top1: Set[str] = {v for v in universe.variant_ids if universe.fitness[v] >= self.top1_threshold}
        self.global_best = float(values.max())
        self.n_top1 = len(self._top1)

        # The qualification bar is resolved HERE, from the same `values` array that defines
        # top1_threshold just above. Using that exact array is the point: it guarantees ONE
        # definition of the eligible universe - including a library's `labelled_stage`
        # exclusion, which is what removes HER2's 25,114 unlabelled rows - so the
        # qualification bar and the top-1% bar are expressed on one scale and cannot drift
        # apart. np.percentile with linear interpolation is the same machinery used for
        # top1_threshold and for veritas.calibration.qualification.absolute_for_quantile.
        if qualify_threshold is not None:
            self.qualify_quantile: Optional[float] = None
            self.qualify_threshold_absolute = float(qualify_threshold)
            self.qualify_threshold_source = "explicit-override"
        else:
            quantile = (
                DEFAULT_QUALIFY_QUANTILE if qualify_quantile is None else float(qualify_quantile)
            )
            if not 0.0 <= quantile <= 100.0:
                raise BenchmarkError(
                    f"qualify_quantile must be a percentile in [0, 100], got {quantile}"
                )
            self.qualify_quantile = quantile
            self.qualify_threshold_absolute = float(np.percentile(values, quantile))
            self.qualify_threshold_source = "quantile"
        # Recorded because a quantile transfers the RARITY of the bar, not its
        # ATTAINABILITY: at the default quantile this is 1 on the CR6261 landscapes, which
        # turns qualification into "find the global best". A run states this rather than
        # having it inferred.
        self.n_qualifying = int(np.count_nonzero(values >= self.qualify_threshold_absolute))

    @property
    def eligible_ids(self) -> List[str]:
        return self._universe.variant_ids

    def reveal(self, variant_ids: Sequence[str]) -> Dict[str, float]:
        """Reveals measured outcomes for an already-committed batch."""
        return {v: self._universe.fitness[v] for v in variant_ids}

    def is_top1(self, variant_id: str) -> bool:
        return variant_id in self._top1

    def export_truth_table(self, method: str) -> Dict[str, float]:
        """
        Hands the complete fitness table to a method that runs its own simulation loop.

        This is the one documented exception to the hidden-results separation, and it is
        recorded here so it shows up in the manifest rather than happening quietly. It exists
        because ALDE is a simulation harness: its loop owns the lookup table and reads an entry
        only after its acquisition function has committed to that point, which preserves the
        commit-then-reveal order even though the array is in its memory. The top-1% threshold
        and membership are NOT exported, and this benchmark re-scores every returned pick with
        its own evaluator.
        """
        self.harness_exports.append(method)
        logger.info("Evaluator exported the full fitness table to harness method '%s'.", method)
        return dict(self._universe.fitness)


@dataclass
class CandidateView:
    """Everything a selector is allowed to see. Deliberately free of outcome information."""

    variant_ids: List[str]
    one_hot: np.ndarray
    index_of: Dict[str, int]
    combos: Dict[str, str]
    #: Library DESIGN - positions, parent residues, alphabet, selection mode. Design
    #: facts are not outcome information, and a selector that cannot see them cannot be
    #: library-aware.
    library: LibrarySpec
    batch_size: int
    n_rounds: int

    def untested(self, revealed: RevealedState) -> List[str]:
        observed = set(revealed.measured)
        return [v for v in self.variant_ids if v not in observed]


# ==============================================================================================
# Selection methods
# ==============================================================================================
class Selector:
    """Interface every method implements. `propose` must not consult the evaluator."""

    name = "selector"
    version = BENCHMARK_VERSION
    config: Dict[str, Any] = {}

    def __init__(self, seed: int) -> None:
        self.seed = seed
        self.model_calls = 0
        self.token_usage: Optional[int] = None  # None == unknown, never guessed
        self.notes: List[str] = []
        # Optional, write-only diagnostics sink set by the harness
        # (veritas.diagnostics.benchmark_hook.CandidatePoolRecorder.attach). A selector
        # that has no scoring stage simply ignores it; nothing here reads it back, so
        # it cannot influence a selection.
        self.candidate_trace: Optional[Any] = None
        # Which selection rule produced the last proposed batch. Written into the round
        # metrics so a reader can see which measurements came from the arm's own policy
        # and which came from a secondary rule - a fallback that is not visible per round
        # is indistinguishable from the policy it replaced.
        self.last_batch_source: str = "policy"

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        raise NotImplementedError


def _rank_select(
    predictions: np.ndarray,
    candidate_ids: Sequence[str],
    n: int,
) -> List[str]:
    """
    Picks the n highest-predicted candidates.

    Tie-breaking is reproducible: ties are broken by the candidate's position in the universe's
    fixed lexicographic order, which `candidate_ids` preserves. np.lexsort is stable on that key.
    """
    order = np.lexsort((np.arange(len(candidate_ids)), -predictions))
    return [candidate_ids[int(i)] for i in order[:n]]


class RandomSelector(Selector):
    """A. Uniform random sampling of untested eligible variants, without replacement."""

    name = "random"
    config = {"strategy": "uniform without replacement"}

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        pool = view.untested(revealed)
        rng = random.Random((self.seed + 1) * 100_003 + round_index)
        return rng.sample(pool, min(n, len(pool)))


class RidgeSelector(Selector):
    """
    B. One-hot ridge regression over the 80 binary features, refit on revealed measurements
    each round, selecting the highest predicted untested variants.

    alpha=1.0 is scikit-learn's documented default and was NOT tuned against held-out results.
    """

    name = "ridge"
    config = {"model": "sklearn.linear_model.Ridge", "alpha": 1.0, "encoding": "onehot_4x20"}

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from sklearn.linear_model import Ridge

        pool = view.untested(revealed)
        if not revealed.measured:
            return RandomSelector(self.seed).propose(round_index, n, view, revealed)
        train_rows = [view.index_of[v] for v in revealed.measured]
        model = Ridge(alpha=1.0)
        model.fit(view.one_hot[train_rows], np.array([revealed.measured[v] for v in revealed.measured]))
        self.model_calls += 1
        predictions = model.predict(view.one_hot[[view.index_of[v] for v in pool]])
        return _rank_select(predictions, pool, min(n, len(pool)))


class RandomForestSelector(Selector):
    """
    C. One-hot random forest over the same encoding and the same revealed observations,
    refit each round. Nonlinear, so it can represent interactions between the four sites.

    Settings are scikit-learn defaults apart from n_estimators=300 (variance reduction) and a
    per-seed random_state for reproducibility. Nothing here was tuned on hidden results.
    """

    name = "random_forest"
    config = {
        "model": "sklearn.ensemble.RandomForestRegressor",
        "n_estimators": 300,
        "max_depth": None,
        "min_samples_leaf": 1,
        "encoding": "onehot_4x20",
    }

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from sklearn.ensemble import RandomForestRegressor

        pool = view.untested(revealed)
        if not revealed.measured:
            return RandomSelector(self.seed).propose(round_index, n, view, revealed)
        train_rows = [view.index_of[v] for v in revealed.measured]
        model = RandomForestRegressor(
            n_estimators=300,
            random_state=self.seed,
            n_jobs=-1,
        )
        model.fit(view.one_hot[train_rows], np.array([revealed.measured[v] for v in revealed.measured]))
        self.model_calls += 1
        predictions = model.predict(view.one_hot[[view.index_of[v] for v in pool]])
        return _rank_select(predictions, pool, min(n, len(pool)))


# ----------------------------------------------------------------------------------------------
# D. Published ALDE
# ----------------------------------------------------------------------------------------------
ALDE_REQUIRED_MODULES = ("torch", "gpytorch", "botorch", "xgboost", "sklearn", "pandas")
ALDE_REQUIRED_FILES = ("execute_simulation.py", "src/optimize.py", "src/objectives.py",
                       "src/acquisition.py", "src/encoding_utils.py")

# Worker executed by the isolated ALDE interpreter. It is written to the isolated workspace at
# run time so that the ALDE checkout itself stays byte-for-byte unmodified, and so that all
# integration logic lives in this one benchmark file.
#
# Faithfulness to the authors' simulation (execute_simulation.py, commit d0b3593):
#   - encoding   : src.encoding_utils.generate_onehot, their own function, reshaped to [n, 80]
#   - objective  : src.objectives.Combo, unmodified (it max-normalises y exactly as in their runs)
#   - model      : DNN_ENSEMBLE with architecture [d, 30, 30, 1], the branch their script selects
#                  for onehot encodings; lr 1e-3, activation lrelu, train_iter 300, dropout 0
#   - acquisition: TS (Thompson sampling)
#   - loop       : src.optimize.BayesianOptimization.run, unmodified
# The only substitution is the initial sample: their utils.samp_discrete(n, obj, seed) draw is
# replaced by this benchmark's shared starting variants, so every method starts from the same 25.
ALDE_WORKER_SOURCE = '''
"""Generated by benchmark_panoptes_isolated.py. Runs the authors' ALDE loop unmodified."""
import json, os, sys, time, warnings
warnings.filterwarnings("ignore")

payload_path, result_path = sys.argv[1], sys.argv[2]
with open(payload_path, "r", encoding="utf-8") as handle:
    cfg = json.load(handle)

sys.path.insert(0, cfg["alde_checkout"])
os.chdir(cfg["work_dir"])

import numpy as np
import pandas as pd
import torch
import random

from src.encoding_utils import generate_onehot
from src.optimize import BayesianOptimization, BO_ARGS
import src.objectives as objectives

combos = np.load(cfg["data_npz"], allow_pickle=True)
combo_list = [str(c) for c in combos["combos"]]
fitness = combos["fitness"].astype(np.float64)

# Lay the data out the way Combo expects, inside the isolated work dir.
data_dir = os.path.join("data", cfg["protein"])
os.makedirs(data_dir, exist_ok=True)
pd.DataFrame({"Combo": combo_list, "fitness": fitness}).to_csv(
    os.path.join(data_dir, "fitness.csv"), index=False)
if not os.path.exists(os.path.join(data_dir, "onehot_x.pt")):
    encoded = generate_onehot(combo_list)
    encoded = torch.reshape(encoded, (encoded.shape[0], -1))
    torch.save(encoded, os.path.join(data_dir, "onehot_x.pt"))

obj = objectives.Combo(cfg["protein"], "onehot")
domain = obj.get_domain()
disc_X = obj.get_points()[0]
disc_y = obj.get_points()[1]

seed = int(cfg["seed"])
torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

# The shared starting sample replaces utils.samp_discrete, so all methods begin identically.
start_indices = torch.tensor(cfg["start_indices"], dtype=torch.long)
start_x = disc_X[start_indices]
start_y = disc_y[start_indices]

arc = [domain[0].size(-1), 30, 30, 1]   # the onehot branch of execute_simulation.py
savedir = os.path.join(cfg["savedir_subdir"], "run_")
os.makedirs(os.path.dirname(savedir), exist_ok=True)

args = BO_ARGS(
    mtype="DNN_ENSEMBLE",
    kernel="RBF",
    acq_fn="TS",
    xi=4,
    architecture=arc,
    activation="lrelu",
    min_noise=1e-6,
    trainlr=1e-3,
    train_iter=int(cfg["train_iter"]),
    dropout=0,
    mcdropout=0,
    verbose=int(cfg["verbose"]),
    bb_fn=obj.objective,
    domain=domain,
    disc_X=disc_X,
    disc_y=disc_y,
    noise_std=0,
    n_rand_init=0,
    budget=int(cfg["budget"]),
    query_cost=1,
    queries_x=start_x,
    queries_y=start_y,
    indices=start_indices,
    savedir=savedir,
    batch_size=int(cfg["batch_size"]),
)

started = time.time()
BayesianOptimization.run(args, seed)
runtime = time.time() - started

selected = torch.load(savedir + "indices.pt").reshape(-1).tolist()
n_initial = len(cfg["start_indices"])
queried = [int(i) for i in selected[n_initial:]]

with open(result_path, "w", encoding="utf-8") as handle:
    json.dump({
        "queried_indices": queried,
        "n_initial": n_initial,
        "runtime_seconds": runtime,
        "versions": {
            "torch": torch.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }, handle)
print("ALDE_WORKER_OK", len(queried))
'''


def _repository_family_roots() -> List[Path]:
    """
    Directories that could hold sibling checkouts of this repository.

    A git worktree lives at `<repo>/.claude/worktrees/<name>`, so `REPO_ROOT.parent` is
    the worktrees folder, not the directory the real checkout sits in. An agent running
    in a worktree therefore looked for a sibling ALDE clone two levels too deep,
    concluded it was "not on disk anywhere" and skipped the arm. This walks back out to
    the real checkout's parent as well, so the search covers both layouts.
    """
    roots: List[Path] = [REPO_ROOT.parent]
    parts = REPO_ROOT.parts
    if ".claude" in parts and "worktrees" in parts:
        index = parts.index(".claude")
        if index >= 1:
            roots.append(Path(*parts[:index]).parent)
    seen: Set[str] = set()
    ordered: List[Path] = []
    for root in roots:
        key = str(root).lower()
        if key not in seen:
            seen.add(key)
            ordered.append(root)
    return ordered


def resolve_alde_paths(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Resolves the ALDE checkout and interpreter, RECORDING every path it looked at.

    ALDE is not vendored into this repository: it is a separate clone with its own
    interpreter, so "unavailable" can mean an unset variable, a path that does not exist,
    or a missing environment. An overnight agent hit exactly this and reported ALDE "not
    on disk anywhere" when it was present at a path the search never considered - which
    is a silent environment problem dressed up as a result.

    Every candidate and its outcome is returned, so the fatal message can say "looked in
    X, Y, Z; not found" instead of a bare refusal, and so a discovered path is marked as
    discovered rather than passing for one the operator pinned.
    """
    def _record(source: str, value: Optional[Any], explicit: bool, is_dir: bool) -> Dict[str, Any]:
        if value is None or str(value).strip() == "":
            return {"source": source, "path": None, "exists": False, "explicit": explicit,
                    "note": "not set"}
        path = Path(str(value)).expanduser()
        try:
            path = path.resolve()
        except OSError:
            pass
        exists = path.is_dir() if is_dir else path.is_file()
        return {
            "source": source,
            "path": str(path),
            "exists": exists,
            "explicit": explicit,
            "note": "" if exists else ("no such directory" if is_dir else "no such file"),
        }

    def _select(candidates: Sequence[Dict[str, Any]]) -> Tuple[Optional[Path], str, str]:
        """
        First existing candidate wins - but an EXPLICIT candidate that does not exist
        stops the search.

        Falling through from a missing `--alde-path` to a discovered sibling clone would
        silently run a different checkout from the one the operator named, and report it
        as though it were theirs. An operator who names a path gets that path or a
        refusal, never a substitute.
        """
        for candidate in candidates:
            if candidate["exists"]:
                candidate["selected"] = True
                return Path(candidate["path"]), candidate["source"], ""
            if candidate["explicit"] and candidate["path"] is not None:
                candidate["blocking"] = True
                return None, "none", (
                    f"{candidate['source']} was set to {candidate['path']} but "
                    f"{candidate['note']}; refusing to fall back to a discovered path, "
                    f"because that would run a different checkout than the one named"
                )
        return None, "none", ""

    candidates: List[Dict[str, Any]] = [
        _record("--alde-path", getattr(args, "alde_path", None), True, True),
        _record("ALDE_PATH", os.getenv("ALDE_PATH"), True, True),
    ]
    for root in _repository_family_roots():
        candidates.append(
            _record(f"sibling of the repository ({root})", root / "ALDE", False, True)
        )
    checkout, checkout_source, checkout_block = _select(candidates)

    interpreters: List[Dict[str, Any]] = [
        _record("--alde-python", getattr(args, "alde_python", None), True, False),
        _record("ALDE_PYTHON", os.getenv("ALDE_PYTHON"), True, False),
    ]
    if checkout is not None:
        for parent in (checkout.parent, checkout):
            for relative in (Path(".venv-alde") / "Scripts" / "python.exe",
                             Path(".venv-alde") / "bin" / "python",
                             Path(".venv") / "Scripts" / "python.exe",
                             Path(".venv") / "bin" / "python"):
                interpreters.append(_record(
                    f"{parent / relative.parts[0]} convention", parent / relative, False, False,
                ))
    python_exe, interpreter_source, interpreter_block = _select(interpreters)

    return {
        "checkout": checkout,
        "checkout_source": checkout_source,
        "python_exe": python_exe,
        "interpreter_source": interpreter_source,
        "searched_checkouts": candidates,
        "searched_interpreters": interpreters,
        "blocked": "; ".join(b for b in (checkout_block, interpreter_block) if b),
    }


def _searched_summary(entries: Sequence[Dict[str, Any]]) -> str:
    """One line naming every path that was tried and why it did not qualify."""
    parts = []
    for entry in entries:
        if entry.get("path") is None:
            parts.append(f"{entry['source']}: not set")
        else:
            parts.append(
                f"{entry['source']}: {entry['path']}"
                + ("" if entry["exists"] else f" ({entry['note']})")
            )
    return "; ".join(parts) or "nothing to search"


def probe_alde(
    checkout: Optional[Path],
    python_exe: Optional[Path],
    resolution: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Reports precisely what is and is not available for running the authors' ALDE implementation
    (https://github.com/jsunn-y/ALDE, Nat Commun 2025, doi:10.1038/s41467-025-55987-8).

    Dependencies are probed in ALDE's OWN interpreter, not this one: the two environments are
    deliberately separate, so importing torch here would prove nothing. No part of ALDE is
    reimplemented: when something is missing, the method is reported UNAVAILABLE with the
    specific cause (checkout, dependency, data or adapter) and nothing is substituted.

    `resolution` carries the record of every path that was searched (see
    `resolve_alde_paths`). It goes into the reason verbatim, because "not found" without
    the list of places looked is unactionable - and because ALDE lives outside this
    repository, which is how a worktree run came to report it missing when it was there.
    """
    searched_checkouts = list((resolution or {}).get("searched_checkouts") or [])
    searched_interpreters = list((resolution or {}).get("searched_interpreters") or [])
    blocked = str((resolution or {}).get("blocked") or "")
    report: Dict[str, Any] = {
        "checkout_path": str(checkout) if checkout else None,
        "checkout_source": (resolution or {}).get("checkout_source", "unknown"),
        "python_exe": str(python_exe) if python_exe else None,
        "interpreter_source": (resolution or {}).get("interpreter_source", "unknown"),
        "searched_checkouts": searched_checkouts,
        "searched_interpreters": searched_interpreters,
        "blocked": blocked,
        "checkout_present": False,
        "interpreter_present": False,
        "missing_files": [],
        "dependencies": {},
        "missing_dependencies": [],
        "imports_alde": None,
        "import_error": None,
        "commit": None,
        "available": False,
        "reason": "",
        "pinned_versions": "alde.yml pins python=3.11, pytorch=2.1.1, gpytorch==1.11, botorch==0.9.4",
        "note": (
            "ALDE is NOT vendored into this repository: it is a separate clone with its "
            "own interpreter, so every path below is outside this tree."
        ),
    }

    if checkout is None:
        report["reason"] = (
            (blocked + ". " if blocked else "")
            + "No usable ALDE checkout. ALDE is a separate clone of "
            "https://github.com/jsunn-y/ALDE, not part of this repository. Looked in: "
            + _searched_summary(searched_checkouts)
            + ". Pass --alde-path /path/to/ALDE or set ALDE_PATH."
        )
        return report
    if not checkout.is_dir():
        report["reason"] = (
            f"ALDE checkout {checkout} (from {report['checkout_source']}) is not a "
            f"directory. Looked in: {_searched_summary(searched_checkouts)}."
        )
        return report

    report["checkout_present"] = True
    report["missing_files"] = [f for f in ALDE_REQUIRED_FILES if not (checkout / f).exists()]
    try:
        report["commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(checkout),
            capture_output=True, text=True, timeout=15,
        ).stdout.strip() or None
    except Exception:
        report["commit"] = None

    if python_exe is None or not Path(python_exe).is_file():
        report["reason"] = (
            (blocked + ". " if blocked else "")
            + f"ALDE's own interpreter was not found (checkout {checkout} came from "
            f"{report['checkout_source']}). Looked in: "
            + _searched_summary(searched_interpreters)
            + ". Create an isolated environment for the checkout and point --alde-python "
              "at its python executable, or set ALDE_PYTHON."
        )
        return report
    report["interpreter_present"] = True

    probe_source = (
        "import json,sys\n"
        "out={'versions':{},'missing':[],'import_error':None}\n"
        "for m in ('torch','gpytorch','botorch','xgboost','sklearn','pandas','numpy'):\n"
        "    try:\n"
        "        out['versions'][m]=getattr(__import__(m),'__version__','unknown')\n"
        "    except Exception as e:\n"
        "        out['missing'].append(m); out['versions'][m]=None\n"
        "try:\n"
        "    import src.optimize, src.objectives, src.acquisition, src.encoding_utils\n"
        "    out['imports_alde']=True\n"
        "except Exception as e:\n"
        "    out['imports_alde']=False; out['import_error']='%s: %s' % (type(e).__name__, e)\n"
        "print(json.dumps(out))\n"
    )
    try:
        completed = subprocess.run(
            [str(python_exe), "-c", probe_source],
            cwd=str(checkout), capture_output=True, text=True, timeout=300,
        )
        parsed = json.loads(completed.stdout.strip().splitlines()[-1])
        report["dependencies"] = parsed["versions"]
        report["missing_dependencies"] = parsed["missing"]
        report["imports_alde"] = parsed["imports_alde"]
        report["import_error"] = parsed["import_error"]
    except Exception as exc:
        report["reason"] = f"could not probe the ALDE interpreter: {type(exc).__name__}: {exc}"
        return report

    problems = []
    if report["missing_files"]:
        problems.append(f"missing repository file(s): {report['missing_files']}")
    if report["missing_dependencies"]:
        problems.append("missing Python dependenc(ies) in the ALDE environment: "
                        + ", ".join(report["missing_dependencies"]))
    if report["imports_alde"] is False:
        problems.append(f"ALDE's own modules do not import: {report['import_error']}")
    report["reason"] = "; ".join(problems) or "checkout, interpreter and dependencies all present"
    report["available"] = not problems
    return report


class AldeSelector(Selector):
    """
    D. The authors' ALDE implementation, run unmodified in its own interpreter.

    ALDE owns its active-learning loop: BayesianOptimization.optimize() retrains and then picks
    `batch_size` points per round until `budget` queries are spent. So this adapter runs one
    campaign per seed through `run_campaign` rather than proposing round by round, and this
    benchmark's evaluator independently re-scores whatever ALDE picked.

    Hidden-results discipline inside ALDE (verified in src/optimize.py at commit d0b3593):
      - the surrogate trains only on `queries_x` / `norm_y`, the points already queried;
      - `disc_y[idx]` is read in `update_trajectory` only AFTER a point has been selected, which
        is the same commit-then-reveal order this benchmark enforces for every other method;
      - the full `disc_y` is otherwise used only for regret logging (`obj_max`).

    Budget convention: ALDE's `budget` EXCLUDES the initial sample, so 25 shared starting
    measurements plus budget=100 gives the same 125 total measurements as every other method.
    """

    name = "alde"

    def __init__(
        self,
        seed: int,
        checkout: Path,
        python_exe: Path,
        workspace: Path,
        probe: Dict[str, Any],
        timeout_seconds: int,
        train_iter: int = 300,
        verbose: int = 1,
    ) -> None:
        super().__init__(seed)
        self.checkout = checkout
        self.python_exe = python_exe
        self.workspace = workspace
        self.probe = probe
        self.timeout_seconds = timeout_seconds
        self.train_iter = train_iter
        self.verbose = verbose
        self._truth: Optional[Dict[str, float]] = None
        self.config = {
            "source": "https://github.com/jsunn-y/ALDE",
            "commit": probe.get("commit"),
            "paper": "doi:10.1038/s41467-025-55987-8",
            "encoding": "onehot (src.encoding_utils.generate_onehot)",
            "model": "DNN_ENSEMBLE, architecture [d, 30, 30, 1]",
            "acquisition": "TS",
            "train_iter": train_iter,
            "lr": 1e-3,
            "activation": "lrelu",
            "budget_convention": "budget excludes the initial sample",
            "execution": "subprocess, separate interpreter",
            "interpreter": str(python_exe),
            "dependency_versions": probe.get("dependencies"),
            "deviation_from_alde_yml": probe.get("deviation", ""),
        }

    def run_campaign(
        self,
        view: CandidateView,
        initial_ids: Sequence[str],
        n_rounds: int,
        batch_size: int,
    ) -> List[List[str]]:
        """Runs one full ALDE campaign and returns its picks split into per-round batches."""
        # The encoded design space is identical across seeds, so it is built once in a shared
        # work dir; only each seed's outputs and logs are per-run. Without this, every seed
        # would rewrite the same ~100 MB one-hot tensor.
        work_dir = self.workspace / "alde_work"
        run_root = work_dir / "runs" / f"seed_{self.seed}_{uuid.uuid4().hex[:6]}"
        run_root.mkdir(parents=True, exist_ok=True)

        data_npz = work_dir / "universe.npz"
        universe_key = hashlib.sha256(
            "".join(view.variant_ids).encode()
        ).hexdigest()
        key_path = work_dir / "universe.key"
        if not data_npz.is_file() or not key_path.is_file() or key_path.read_text().strip() != universe_key:
            np.savez_compressed(
                data_npz,
                combos=np.array([view.combos[v] for v in view.variant_ids], dtype=object),
                fitness=np.array([self._truth_placeholder(v) for v in view.variant_ids], dtype=np.float64),
            )
            key_path.write_text(universe_key, encoding="utf-8")
            for stale in (work_dir / "data").glob("**/onehot_x.pt"):
                stale.unlink()  # the design space changed, so the cached encoding is invalid
        payload = {
            "alde_checkout": str(self.checkout),
            "work_dir": str(work_dir),
            "savedir_subdir": str(run_root),
            "data_npz": str(data_npz),
            # ALDE uses this only to name its own cache and output folders; it must track
            # the library being run or two libraries share one encoded design space on
            # disk. Taken from the candidate view, which carries its own library.
            "protein": f"{view.library.key}_panoptes",
            "seed": self.seed,
            "budget": n_rounds * batch_size,   # excludes the initial sample, per ALDE's convention
            "batch_size": batch_size,
            "train_iter": self.train_iter,
            "verbose": self.verbose,
            "start_indices": [view.index_of[v] for v in initial_ids],
        }
        payload_path = run_root / "payload.json"
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        worker_path = run_root / "alde_worker_generated.py"
        worker_path.write_text(ALDE_WORKER_SOURCE, encoding="utf-8")
        result_path = run_root / "result.json"

        environment = dict(os.environ)
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["OMP_NUM_THREADS"] = environment.get("OMP_NUM_THREADS", "4")
        started = time.perf_counter()
        completed = subprocess.run(
            [str(self.python_exe), str(worker_path), str(payload_path), str(result_path)],
            cwd=str(run_root), capture_output=True, text=True,
            timeout=self.timeout_seconds, env=environment,
        )
        elapsed = time.perf_counter() - started
        (run_root / "stdout.log").write_text(completed.stdout or "", encoding="utf-8")
        (run_root / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")

        if completed.returncode != 0 or not result_path.is_file():
            tail = (completed.stderr or completed.stdout or "").strip().splitlines()[-12:]
            raise BenchmarkError(
                f"ALDE subprocess failed (exit {completed.returncode}) after {elapsed:.1f}s. "
                f"Logs: {run_root}. Last output:\n  " + "\n  ".join(tail)
            )

        result = json.loads(result_path.read_text(encoding="utf-8"))
        queried = result["queried_indices"]
        self.model_calls += n_rounds  # one surrogate refit per round
        self.notes.append(
            f"ALDE commit {self.probe.get('commit', 'unknown')[:12]} returned {len(queried)} "
            f"queries in {result['runtime_seconds']:.1f}s (torch {result['versions']['torch']})"
        )
        batches: List[List[str]] = []
        for start in range(0, n_rounds * batch_size, batch_size):
            chunk = queried[start:start + batch_size]
            batches.append([
                view.variant_ids[i] for i in chunk if 0 <= i < len(view.variant_ids)
            ])
        return batches

    def _truth_placeholder(self, variant_id: str) -> float:
        """
        Supplies the measured value for ALDE's `disc_y`.

        ALDE is a simulation harness: its loop reads `disc_y[idx]` to reveal the outcome of a
        point it has just chosen, exactly as this benchmark's evaluator reveals a committed
        batch. The surrogate never trains on unqueried entries. The injected table is set by
        `attach_truth` before the campaign runs.
        """
        if self._truth is None:
            raise BenchmarkError("ALDE campaign started before the evaluator exported its table.")
        return self._truth[variant_id]

    def attach_truth(self, truth: Dict[str, float]) -> None:
        """Called by the runner only for harness-style methods that run their own loop."""
        self._truth = truth

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        raise BenchmarkError("ALDE runs as a campaign; run_campaign() should have been used.")


# ----------------------------------------------------------------------------------------------
# E/F. The existing PanoptesBio production path
# ----------------------------------------------------------------------------------------------
class _LiveCallRecorder(logging.Handler):
    """
    Observes the orchestrator's own log stream to record whether Muse calls actually went live.

    Production falls back to the mock planner when the Meta API call fails. That fallback is
    silent to the caller, so a run could otherwise be mislabelled 'live'. This handler only
    reads log records; it changes no production behaviour.
    """

    def __init__(self) -> None:
        super().__init__(level=logging.INFO)
        self.live_attempts = 0
        self.fallbacks = 0
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        try:
            message = record.getMessage()
        except Exception:
            return
        if "Querying Meta Model API" in message:
            self.live_attempts += 1
        if "Meta Model API" in message and "failed" in message:
            self.fallbacks += 1
            self.messages.append(message[:200])


class _AnnotationDrivenPlanner:
    """
    The non-live planner used by `panoptes_default`.

    Why this exists instead of `workers.mock_workers.MockLLMWorker`:
        MockLLMWorker is a scripted demo planner. For any accession matching TrpB it
        returns the literal positions 183/184/227/228 - exactly the four saturated
        sites of this benchmark's library - and `suggested_chemistries` naming three
        of the four residues of the published global optimum V183A:F184I:V227K:S228G
        (183:[A,V,C], 184:[I,V,L], 228:[G,A,S]). That is appropriate for a product
        demo and disqualifying for a baseline: it encodes the answer this benchmark
        is supposed to measure the discovery of.

    What this does instead:
        Derives the search space from the target's own annotations, which is how a
        real campaign knows what its ordered library can assay. Positions come from
        `metadata.binding_positions` (populated from the entry's annotated active
        sites), catalytic machinery from `metadata.catalytic_positions` is forbidden,
        and chemistries are a fixed neutral alphabet-spanning set chosen for physical
        diversity, never for fitness on this landscape. Nothing here mentions TrpB,
        so the same code path serves any target with annotated sites.
    """

    # A neutral, target-independent probe set spanning small/polar/charged/branched/
    # aromatic classes. Deliberately NOT tuned to any landscape: it must not encode
    # which residues win. GFlowNet sampling in step 3 explores beyond these.
    NEUTRAL_CHEMISTRIES: Tuple[str, ...] = ("A", "S", "D", "K", "L", "F")

    def generate_strategy(
        self,
        metadata: Any,
        round_number: int = 1,
        top_evidence: Optional[Sequence[Dict[str, Any]]] = None,
        current_hop: int = 0,
    ) -> Any:
        from controller.prompts import BiophysicalHypothesis, ResidueCluster, SearchStrategy

        annotated = sorted({int(p) for p in (getattr(metadata, "binding_positions", None) or [])})
        catalytic = sorted({int(p) for p in (getattr(metadata, "catalytic_positions", None) or [])})
        targetable = [p for p in annotated if p not in set(catalytic)]
        if not targetable:
            raise BenchmarkError(
                "annotation-driven planner found no targetable annotated positions for "
                f"{getattr(metadata, 'accession', '<unknown>')}: binding_positions={annotated}, "
                f"catalytic_positions={catalytic}. panoptes_default cannot plan without "
                "annotated sites; supply them in the target's metadata."
            )

        # Split the annotated sites into two coupled clusters so the production
        # selector sees interacting groups rather than one undifferentiated block.
        midpoint = (len(targetable) + 1) // 2
        groups = [targetable[:midpoint], targetable[midpoint:]]
        clusters = [
            ResidueCluster(
                cluster_id=f"annotated_site_group_{i + 1}",
                positions=group,
                rationale="Annotated active/binding-site positions reported for this target.",
            )
            for i, group in enumerate(groups)
            if group
        ]

        chemistries = {str(p): list(self.NEUTRAL_CHEMISTRIES) for p in targetable}
        hypothesis = (
            f"Round {round_number}: combinatorial exploration of the {len(targetable)} "
            "annotated active-site positions, preserving annotated catalytic machinery."
        )
        return SearchStrategy(
            target_protein=metadata.protein_name,
            target_accession=metadata.accession,
            hypothesis=hypothesis,
            target_clusters=clusters,
            forbidden_positions=catalytic,
            # Depth tracks the annotated site count so a combinatorially complete
            # library is reachable, rather than being fixed to a literal.
            max_mutation_depth=min(len(targetable), 2 if round_number == 1 else 4),
            mechanistic_hypotheses=[
                BiophysicalHypothesis(
                    pocket_name="annotated_active_site",
                    target_residues=targetable,
                    mechanism="Combinatorial substitution across annotated active-site positions.",
                    structural_rationale=(
                        f"Positions {targetable} are the annotated substrate/cofactor contact "
                        f"sites for this entry; catalytic positions {catalytic} are withheld."
                    ),
                    suggested_chemistries=chemistries,
                    epistatic_coupling="Annotated co-located sites are treated as potentially coupled.",
                )
            ],
        )


def parse_ablate_components(spec: Optional[str]) -> Tuple[str, ...]:
    """
    Parses --ablate-components into a normalised, validated tuple.

    The old implementation was `{tok.strip() for tok in (spec or "").split(",") if tok.strip()}`,
    which silently swallowed every malformed value:
      * '--ablate-components ,'  and  '--ablate-components " "'  produced the empty set, so a
        folder labelled as an ablation arm actually ran the complete engine;
      * 'dpp,,gflownet' quietly dropped the empty token, hiding a typo;
      * an unknown name was passed straight through to select_next_micro_batch, which raises -
        but only once the first round of the first seed was already running, by which time the
        run folder existed and the failure surfaced as a per-seed 'failed' row.

    Rules now:
      * None, or the argument's own default empty string, means ablate nothing;
      * any other value that names no component is a hard error;
      * empty tokens (leading, trailing or doubled separators) are a hard error;
      * unknown names are rejected here, before anything runs, naming the valid set;
      * names are de-duplicated and sorted so the manifest's arm identity is canonical.
    """
    if spec is None:
        return ()
    raw = str(spec)
    if raw == "":
        return ()  # the parser default: explicitly 'disable no components'
    if not raw.strip():
        raise ConfigurationError(
            "--ablate-components was given a whitespace-only value, which names no component. "
            "Omit the flag (or pass an empty string) to ablate nothing, or name components "
            f"from: {', '.join(ABLATABLE_COMPONENTS)}."
        )

    names: List[str] = []
    for token in raw.split(","):
        name = token.strip()
        if not name:
            raise ConfigurationError(
                f"--ablate-components {raw!r} contains an empty component name. Remove the "
                "stray separator; an empty entry used to be dropped silently, so a typo "
                "produced an unablated run reported as an ablation arm."
            )
        names.append(name)

    unknown = sorted({n for n in names if n not in ABLATABLE_COMPONENTS})
    if unknown:
        raise ConfigurationError(
            f"unknown --ablate-components name(s) {unknown}. Recognised components: "
            f"{', '.join(ABLATABLE_COMPONENTS)}. (Names are case-sensitive and must match "
            "workers.math_worker.ABLATABLE_COMPONENTS.)"
        )
    return tuple(sorted(set(names)))


def component_states(disabled: Sequence[str]) -> Dict[str, bool]:
    """Maps every ablatable component to whether it was ENABLED in this run."""
    disabled_set = set(disabled)
    return {name: name not in disabled_set for name in ABLATABLE_COMPONENTS}


def resolve_ablate_components(args: argparse.Namespace) -> Tuple[str, ...]:
    """
    Returns the validated component list for a parsed argument namespace.

    main() parses once, before any folder is created, and caches the result on the
    namespace; everything downstream reads that cache so the value cannot drift.
    """
    cached = getattr(args, "ablate_components_resolved", None)
    if cached is None:
        cached = parse_ablate_components(getattr(args, "ablate_components", ""))
        args.ablate_components_resolved = cached
    return tuple(cached)


def enforce_requested_methods_available(
    methods: Sequence[str],
    method_status: Dict[str, Dict[str, Any]],
    allow_unavailable: bool = False,
) -> None:
    """
    Fails the run when a method that was explicitly requested cannot execute.

    This is the fix for the incident that cost two benchmark runs: `--methods ... ,alde`
    with ALDE_PATH unset produced `status=UNAVAILABLE` as one row of a report, the run
    completed, and the comparison was quoted as though the baseline had been there. Any
    requested method whose status is not 'available' - UNAVAILABLE, UNSUPPORTED or SKIPPED
    - now aborts before a single measurement is spent.

    `allow_unavailable` (--allow-unavailable-baselines) restores the report-only behaviour
    for the one legitimate use: probing an environment to see what it can run.
    """
    broken: List[Tuple[str, str, str]] = []
    for method in methods:
        status = method_status.get(method, {"status": "UNSUPPORTED", "detail": "no status recorded"})
        state = str(status.get("status", "unknown"))
        if state.lower() != "available":
            broken.append((method, state, str(status.get("detail", "")).strip()))
    if not broken:
        return

    lines = [
        f"{len(broken)} requested method(s) cannot run, so this comparison would be missing "
        "an arm it claims to measure:",
        "",
    ]
    for method, state, detail in broken:
        lines.append(f"  {method}: {state}")
        for chunk in (detail.split("; ") if detail else []):
            if chunk.strip():
                lines.append(f"      - {chunk.strip()}")
    lines += [
        "",
        "Fix the environment, or drop the method from --methods. To probe availability",
        "without running a comparison, pass --allow-unavailable-baselines, which restores",
        "the old behaviour of reporting the arm as a row instead of failing.",
    ]
    message = "\n".join(lines)
    if allow_unavailable:
        logger.warning(
            "--allow-unavailable-baselines: continuing without %s. THIS RUN IS NOT A "
            "COMPLETE COMPARISON.", ", ".join(m for m, _, _ in broken),
        )
        return
    raise BaselineUnavailableError(message)


class PanoptesSelector(Selector):
    """
    E/F. The real production planning + selection path.

    Production functions called, in production order:
        workers.rag_worker.fetch_protein_metadata      (step 1, cached, no network)
        _AnnotationDrivenPlanner.generate_strategy     (step 2, annotation-derived search space) [E]
        controller.orchestrator.Orchestrator._generate_strategy (step 2, live Muse planner)  [F]
        workers.math_worker.select_next_micro_batch    (step 3, DKL + GFlowNet + Pareto + k-DPP)
        controller.orchestrator.Orchestrator._generate_reflection (step 6b, feedback)        [F]

    Stages NOT exercised, and why:
        Step 4 (Gatekeeper spend authorisation) and step 5 (Modal ESM-2 scoring) are skipped.
        `Orchestrator.run_stream()` hard-wires step 5 to `workers.esm_worker.score_variants`,
        and the loop learns from those ESM pseudo-likelihoods via the ledger's evidence table.
        There is no supported interface for substituting an experimental assay as the outcome
        source, so running the full loop would optimise ESM scores rather than measured fitness
        and would consume none of the experimental budget. Pricing an ESM call this benchmark
        never makes would also fabricate spend. The planning, selection and feedback stages that
        decide WHICH variants are measured are exercised in full.

    Evidence handed to production contains only measurements already revealed to this method.
    """

    def __init__(
        self,
        seed: int,
        wt_sequence: str,
        metadata: Any,
        use_live_muse: bool,
        pdb_path: Optional[Path],
        library: LibrarySpec,
        call_budget: Optional[int] = None,
        qualify_threshold: Optional[float] = None,
        probe_after_qualify: bool = False,
        ablation: str = "none",
        prior_spread: Optional[float] = None,
        ablate_components: Optional[Set[str]] = None,
    ) -> None:
        super().__init__(seed)
        # Ablation arm. 'none' is the shipped configuration; the others restore one
        # or both pre-fix behaviours so their effect can be measured on identical
        # seeds. Recorded in the arm name so a report cannot confuse two arms.
        self.ablation = ablation
        # Component-level ablation, for attributing the engine's advantage. Recorded
        # in the arm's notes so a run folder states what it measured - without this,
        # five ablation folders are distinguishable only by their timestamps.
        self.ablate_components = set(ablate_components or ())
        if self.ablate_components:
            self.notes.append(
                f"ABLATED_COMPONENTS={','.join(sorted(self.ablate_components))}"
            )
        if qualify_threshold is not None and probe_after_qualify:
            self.name = "panoptes_conditional_probe"
        elif qualify_threshold is not None:
            self.name = "panoptes_qualify_then_optimize"
        else:
            self.name = "panoptes_muse" if use_live_muse else "panoptes_default"
        # Which policy runs after qualification: deterministic conditional enumeration
        # (True) or greedy local frontier search (False).
        self.probe_after_qualify = bool(probe_after_qualify)
        # The ablation is NOT folded into self.name: the method registry and the
        # --methods filter both key on that name, so renaming the arm matched nothing
        # and silently ran zero seeds. It goes in the description instead.
        # The prior calibration is a module-level switch because raising its spread
        # cannot undo it (tanh saturates). Set per selector construction; seeds run
        # in separate processes under --parallel-seeds, so this does not race.
        import workers.math_worker as _mw

        # Honour the module default when no ablation is requested.
        #
        # This previously read `ablation not in ("raw-prior", "legacy")`, and `ablation`
        # is None for an ordinary run - so `None not in (...)` is True and EVERY default
        # run silently forced calibration ON, overriding math_worker's default of False.
        # That is the configuration docs/545.md §2 measured as worse and reverted: a
        # plain `--methods panoptes_default` produced 0.7165 / 33.83 against the
        # documented 0.7462 / 45.60, reproducing the reverted arm rather than the
        # shipped one. It also made `--prior-spread 8.0` look like an improvement when
        # it was only loose enough to approximate the uncalibrated behaviour.
        # Valid --prior-ablation values are ("none", "aux-weighted", "raw-prior",
        # "legacy"). Only the two that name the prior change it; "none" and
        # "aux-weighted" leave the shipped default in place.
        if ablation in ("raw-prior", "legacy"):
            _mw.PRIOR_CALIBRATION_ENABLED = False
        if prior_spread is not None:
            # An explicit spread means calibrate at that weight, overriding the
            # ablation default: this is the sweep that settles the prior's weight by
            # measurement instead of by argument.
            _mw.PRIOR_CALIBRATION_ENABLED = True
            _mw.PRIOR_FITNESS_SPREAD = float(prior_spread)
            self.notes.append(f"PRIOR_FITNESS_SPREAD={float(prior_spread)}")
        if ablation and ablation != "none":
            self.notes.append(
                f"ABLATION={ablation}: "
                f"prior_calibration={_mw.PRIOR_CALIBRATION_ENABLED}, "
                f"auxiliary_objective_weights="
                f'{ablation in ("aux-weighted", "legacy")}'
            )
        # When set, the selector switches to frontier search once a measured variant
        # reaches this fitness. None keeps the single-policy baseline.
        self.qualify_threshold = qualify_threshold
        self.qualified_at_round: Optional[int] = None
        self.wt_sequence = wt_sequence
        self.metadata = metadata
        # The library being run, held on the selector rather than read from module
        # state. Everything target-specific - the assayable positions handed to the
        # planner, the accession, the encoders, pool-restricted mode - comes from here.
        self.library = library
        self.use_live_muse = use_live_muse
        self.pdb_path = pdb_path
        self.call_budget = call_budget
        self.config = {
            "strategy_source": (
                "live Muse (Meta API)" if use_live_muse
                else "annotation-driven planner (metadata.binding_positions; no scripted per-target literals)"
            ),
            "selection": "workers.math_worker.select_next_micro_batch",
            "acquisition": "hybrid_thompson (production default)",
            "reflection_feedback": bool(use_live_muse),
            "stages_exercised": ["1_init", "2_strategy", "3_selection"] + (["6b_reflection"] if use_live_muse else []),
            "stages_skipped": ["4_authorization", "5_esm_execution"],
        }
        self._orchestrator = None
        self._mock_planner = None
        self.recorder = _LiveCallRecorder()
        self.live_calls = 0
        self.reflection_calls = 0
        self.last_batch_predictions: Dict[str, float] = {}

    # -- production handles -------------------------------------------------------------------
    def _planner(self):
        if self.use_live_muse:
            if self._orchestrator is None:
                from controller.orchestrator import Orchestrator, OrchestratorConfig

                self._orchestrator = Orchestrator(
                    OrchestratorConfig(
                        target_accession=self.library.accession,
                        use_mock_llm=False,
                        enable_live_llm=True,
                        enable_reflection=True,
                        random_seed=self.seed,
                    ),
                    # The saturated sites ARE this benchmark's ordered library:
                    # declaring them is experimental design, not a hint about which
                    # residues win. Without it the planner proposes distal allosteric
                    # residues that the library cannot measure, and whole rounds are
                    # discarded as not_in_eligible_universe.
                    assayable_positions=self.library.positions,
                )
            return self._orchestrator
        if self._mock_planner is None:
            self._mock_planner = _AnnotationDrivenPlanner()
        return self._mock_planner

    def _strategy(self, round_index: int, revealed: RevealedState):
        """Step 2, through whichever production planner this method is configured for."""
        top_evidence = [
            {"variant_sequence": self._sequence_for(v), "computed_score": score}
            for v, score in sorted(revealed.measured.items(), key=lambda kv: kv[1], reverse=True)[:5]
        ]
        planner = self._planner()
        if self.use_live_muse:
            if self.call_budget is not None and self.live_calls >= self.call_budget:
                raise BenchmarkError(
                    f"live model-call cap of {self.call_budget} reached (--max-model-calls)"
                )
            logging.getLogger("controller.orchestrator").addHandler(self.recorder)
            before = self.recorder.live_attempts
            try:
                strategy = planner._generate_strategy(
                    self.metadata, round_number=round_index + 1, top_evidence=top_evidence, current_hop=0,
                )
            finally:
                logging.getLogger("controller.orchestrator").removeHandler(self.recorder)
            if self.recorder.live_attempts > before:
                self.live_calls += 1
                self.model_calls += 1
            return strategy
        return planner.generate_strategy(
            self.metadata, round_number=round_index + 1, top_evidence=top_evidence, current_hop=0,
        )

    def _sequence_for(self, variant_id: str) -> str:
        from controller.gatekeeper import Variant

        return Variant(notation=variant_id).apply_to(self.wt_sequence)

    def _reflect(self, round_index: int, revealed_batch: Dict[str, float]) -> None:
        """Step 6b: production residual reflection, fed back exactly as run_stream() does."""
        from controller.prompts import ResidualOutlier

        outliers: List[ResidualOutlier] = []
        for variant_id, observed in revealed_batch.items():
            predicted = self.last_batch_predictions.get(variant_id, 0.0)
            outliers.append(ResidualOutlier(
                variant_notation=variant_id,
                predicted_score=round(float(predicted), 4),
                observed_score=round(float(observed), 4),
                residual=round(float(observed) - float(predicted), 4),
                mutation_details=variant_id.split(":"),
            ))
        outliers.sort(key=lambda o: o.residual, reverse=True)
        positive = [o for o in outliers if o.residual > 0.05][:3]
        negative = [o for o in reversed(outliers) if o.residual < -0.05][:3]
        if not positive and not negative:
            return
        if self.call_budget is not None and self.live_calls >= self.call_budget:
            self.notes.append("reflection skipped: live model-call cap reached")
            return
        orchestrator = self._planner()
        logging.getLogger("controller.orchestrator").addHandler(self.recorder)
        before = self.recorder.live_attempts
        try:
            reflection = orchestrator._generate_reflection(
                metadata=self.metadata,
                round_number=round_index + 1,
                positive_outliers=positive,
                negative_outliers=negative,
                pdb_id=getattr(self.metadata, "pdb_id", None),
            )
        finally:
            logging.getLogger("controller.orchestrator").removeHandler(self.recorder)
        # Production stores the diagnosis on the orchestrator; the next strategy merges it in.
        orchestrator.last_reflection = reflection
        self.reflection_calls += 1
        if self.recorder.live_attempts > before:
            self.live_calls += 1
            self.model_calls += 1

    # -- selection ----------------------------------------------------------------------------
    def _conditional_probe_batch(
        self, n: int, view: CandidateView, revealed: RevealedState,
    ) -> List[str]:
        """
        Deterministic conditional enumeration: fix the best marginal identity, then
        ENUMERATE a second position exhaustively.

        What this tests. A substitution whose benefit is conditional on another has a
        poor MARGINAL effect, so any surrogate fitted on marginals down-weights it, so
        it is never measured and the condition stays undiscoverable. Measured on TrpB:
        K227 mean marginal -0.021, better than parent in 38.9% of backgrounds, present
        in 100% of variants above 0.90 fitness.

        Why enumeration and not ranking. A previous version of this reserved batch slots
        for carriers of the best identity and then filled them BY ACQUISITION VALUE,
        which inherited the exact preference it was meant to bypass: fixing the good
        identity and ranking its carriers returns carriers without the conditional
        partner. Measured: K227 sampling unchanged at 1.00%, no tail hits, t=1.21. The
        acquisition function must be excluded from this block, not consulted by it.

        Mechanism, using only this campaign's own measurements:
          1. anchor   - the (position, residue) with the best mean observed fitness.
          2. scan     - a second position, the one with the fewest distinct residues
                        measured so far in the anchor's presence, so the block buys
                        coverage where the conditional structure is least constrained.
          3. enumerate- every residue at the scan position, holding the anchor fixed and
                        the remaining positions at the best observed background.

        This asserts no biology and names no residue: the anchor comes from measurements
        and the scan position from a coverage count. On a 4-site library a full scan is
        at most 20 variants, so it fits inside one batch.
        """
        positions = list(self.library.positions)

        # The whole background is the best MEASURED variant, not a per-identity mean.
        #
        # An earlier version picked the anchor as the (position, residue) with the best
        # mean observed fitness. With 25 starting measurements spread over 4 x 20 = 80
        # identities, each identity has one to three observations, so that mean is
        # noise: on a planted-epistasis landscape the probe anchored on an arbitrary
        # residue every round and never found the planted pair. Conditioning on an
        # actually-measured high-fitness variant needs no per-identity estimate at all.
        if not revealed.measured:
            return []
        best_id = max(revealed.measured, key=lambda v: revealed.measured[v])
        background = view.combos.get(best_id)
        if not background:
            return []

        # Enumerate the position with the fewest residues already measured in this
        # background's company — the one where the conditional structure is least
        # constrained by what we have seen.
        seen_at: Dict[int, Set[str]] = {s: set() for s in range(len(positions))}
        for vid in revealed.measured:
            combo = view.combos.get(vid)
            if not combo:
                continue
            # Count only variants that share most of this background, so coverage is
            # measured in the relevant neighbourhood rather than library-wide.
            #
            # Both genotypes are full length over the library's positions by
            # construction, and that is asserted rather than assumed: zipping a 3-site
            # genotype against a 4-site one compares the wrong sites and says nothing.
            if len(combo) != len(background):
                raise BenchmarkError(
                    f"genotype length mismatch: {combo!r} ({len(combo)}) vs background "
                    f"{background!r} ({len(background)}) on library {self.library.key!r}. "
                    f"Canonicalise to full length before comparing; never zip unequal "
                    f"genotypes."
                )
            shared = sum(1 for a, b in zip(combo, background) if a == b)
            if shared >= len(positions) - 2:
                for slot, residue in enumerate(combo):
                    seen_at[slot].add(residue)
        scan_slot = min(range(len(positions)), key=lambda s: len(seen_at[s]))

        proposals: List[str] = []
        ballot: List[str] = []
        if self.library.pool_restricted:
            # Same defect as _frontier_batch: enumerating a residue scan over a sparse
            # library produces genotypes nobody assayed, the batch comes back empty or
            # near-empty, and propose() falls through to the baseline policy. On a
            # pool-restricted library the scan is taken from the ENUMERATED POOL instead:
            # the untested members that share this background at every position except
            # the scan slot. That is the same conditional question - vary one position,
            # hold the rest - asked only of candidates that exist.
            for vid in view.untested(revealed):
                combo = view.combos.get(vid)
                if not combo or len(combo) != len(background):
                    continue
                differs = [s for s in range(len(background)) if combo[s] != background[s]]
                if differs == [scan_slot]:
                    ballot.append(vid)
                    if vid not in proposals:
                        proposals.append(vid)
                if len(proposals) >= n:
                    break
            self._note_ballot(ballot, "conditional_probe_pool_restricted")
            if not proposals:
                # The conditional probe is a policy that genuinely cannot be expressed
                # here: this background has no in-pool single-position variants at the
                # chosen slot. Say so rather than degrading into the baseline arm.
                raise BenchmarkError(
                    f"conditional probe cannot be expressed on pool-restricted library "
                    f"{self.library.key!r}: the untested pool holds no variant that differs "
                    f"from background {best_id!r} at position {positions[scan_slot]} alone. "
                    f"A sparse library does not contain the conditional slice this policy "
                    f"enumerates. Refusing to return an empty batch, which would silently "
                    f"degrade this arm into the baseline policy."
                )
            self.notes.append(
                f"conditional probe (pool-restricted): background {best_id or 'WT'} "
                f"(fitness {revealed.measured[best_id]:.4f}), position "
                f"{positions[scan_slot]} ({len(proposals)} in-pool variants)"
            )
            return proposals

        for residue in AA_ALPHABET:
            if residue == background[scan_slot]:
                continue
            combo = list(background)
            combo[scan_slot] = residue
            vid = _canonical_variant_id("".join(combo), self.library)
            if vid and vid in view.index_of and vid not in revealed.measured:
                ballot.append(vid)
                if vid not in proposals:
                    proposals.append(vid)
            if len(proposals) >= n:
                break
        self._note_ballot(ballot, "conditional_probe_scan")
        if proposals:
            self.notes.append(
                f"conditional probe: background {best_id or 'WT'} "
                f"(fitness {revealed.measured[best_id]:.4f}), enumerated position "
                f"{positions[scan_slot]} ({len(proposals)} variants)"
            )
        return proposals

    def _frontier_batch(self, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        """
        Greedy local search around the best measured variant, for use AFTER the
        campaign has already secured a qualifying lead.

        Enumerates untested Hamming-1 neighbours of the current best and ranks them
        with a ridge surrogate fitted on the revealed data. No diversity term, no
        Pareto weighting, no prior: deliberately the opposite of the baseline policy.

        The point is not that greedy search is better. It is that the two objectives -
        secure an acceptable lead, then improve it - are lexicographic, and a single
        acquisition function tuned for the first cannot serve the second. Our measured
        failure is specifically at the frontier: the baseline reaches >=0.70 in 47 of 48
        campaigns and >=0.90 in none, while a plain ridge regression reaches >=0.90 in
        about 17% by concentrating its sampling.
        """
        from sklearn.linear_model import Ridge

        best_id = max(revealed.measured, key=lambda v: revealed.measured[v])
        best_combo = view.combos.get(best_id)
        if not best_combo:
            return []

        # Where the candidates come from depends on the library's SELECTION MODE, exactly
        # as it does in propose(). propose() was made pool-restriction-aware and this was
        # not, which is the whole of defect 3: on HER2 the median library member has ZERO
        # in-library Hamming-1 neighbours (mean 0.18 of 190 possible), so this enumeration
        # returned an empty list, propose() fell through to the pool path, and in 41 of 48
        # seeds the frontier arm selected the IDENTICAL 125 variants as the reliability arm.
        # Its reported 3.062 against reliability's 3.125 was a budget penalty, not a policy
        # difference - the arm was not a distinct policy at all.
        if self.library.pool_restricted:
            # On a sparse library the design space IS the enumerated pool, so the frontier
            # policy ranks the untested pool rather than Hamming-1 neighbours nobody
            # assayed. The policy keeps its character - a greedy ridge surrogate with no
            # diversity term, no Pareto weighting and no prior, deliberately the opposite
            # of the baseline - over the candidates that actually exist.
            candidates = [
                (vid, view.combos[vid])
                for vid in view.untested(revealed)
                if vid in view.combos
            ]
            source = "frontier_pool_restricted"
        else:
            # Untested Hamming-1 neighbours of the incumbent: 4 positions x 19 residues.
            candidates = []
            for slot in range(self.library.n_sites):
                for aa in AA_ALPHABET:
                    if aa == best_combo[slot]:
                        continue
                    combo = best_combo[:slot] + aa + best_combo[slot + 1:]
                    vid = _canonical_variant_id(combo, self.library)
                    if vid in view.index_of and vid not in revealed.measured:
                        candidates.append((vid, combo))
            source = "frontier_hamming1"

        if not candidates:
            if self.library.pool_restricted:
                # The pool is enumerated and finite, so an empty pool here means the
                # library is exhausted, not that the policy is inexpressible. Fail loudly:
                # returning [] would send propose() back to the baseline policy and the arm
                # would silently stop being the frontier arm.
                raise BenchmarkError(
                    f"frontier policy has no candidates on pool-restricted library "
                    f"{self.library.key!r}: every eligible variant is already measured. "
                    f"Refusing to return an empty batch, which would silently degrade this "
                    f"arm into the baseline policy."
                )
            return []

        # Rank them by a ridge surrogate fitted on everything measured so far.
        rows, ys = [], []
        for vid, score in revealed.measured.items():
            combo = view.combos.get(vid)
            if combo:
                rows.append(_one_hot_matrix([combo], self.library)[0])
                ys.append(score)
        if len(rows) >= 10:
            model = Ridge(alpha=1.0).fit(np.asarray(rows), np.asarray(ys))
            X = _one_hot_matrix([c for _, c in candidates], self.library)
            order = np.argsort(-model.predict(X))
        else:
            order = np.arange(len(candidates))
        self._note_ballot([vid for vid, _ in candidates], source)
        return [candidates[int(i)][0] for i in order[:n]]

    def _note_ballot(self, candidates: Sequence[str], source: str) -> None:
        """
        Records a post-qualification ballot for coverage.

        These generators never reach `select_next_micro_batch`, so without this their
        candidates are invisible to the coverage diagnostic and the arm's ballot is
        under-counted - the same error as computing coverage from committed picks.
        Coverage only: there is no acquisition value here to trace per candidate.
        """
        sink = self.candidate_trace
        if sink is None:
            return
        try:
            sink.note_pool(candidates, source=source)
        except Exception as exc:
            logger.debug("ballot note skipped (%s): %s", source, exc)

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from workers.math_worker import select_next_micro_batch

        # What this round's measurements should be attributed to if execution reaches the
        # baseline path below. It is "baseline" for the baseline arms, and a named
        # fallback for a post-qualification arm whose own policy produced nothing.
        baseline_source = "baseline"

        # Qualify-then-optimise: run the baseline policy unchanged until a measured
        # variant clears the qualifying threshold, then spend the remaining budget on
        # frontier search. Under this benchmark's deterministic reveal, every trajectory
        # is identical to the baseline up to and including the round that qualifies, so
        # P(reach threshold) and time-to-threshold are preserved pathwise - the switch
        # can only affect what happens after the protected event. It does NOT preserve
        # top-1% yield or final best fitness, which is exactly what is being tested.
        if self.qualify_threshold is not None and revealed.measured:
            if max(revealed.measured.values()) >= self.qualify_threshold:
                # After qualification, spend the remaining budget on whichever
                # post-qualification policy this arm is configured for.
                second_phase = (
                    self._conditional_probe_batch(n, view, revealed)
                    if self.probe_after_qualify
                    else self._frontier_batch(n, view, revealed)
                )
                mode = "conditional-probe" if self.probe_after_qualify else "frontier"
                if second_phase:
                    if not self.qualified_at_round:
                        self.qualified_at_round = round_index
                        self.notes.append(
                            f"qualified before round {round_index}; {mode} mode"
                        )
                    self.last_batch_predictions = {}
                    self.last_batch_source = mode
                    return second_phase
                # The post-qualification policy produced nothing and this round is about
                # to run the BASELINE policy instead. In pool-restricted mode that cannot
                # happen - both second-phase generators now raise - so this is the
                # generative case, where the incumbent's Hamming-1 shell can genuinely be
                # exhausted late in a campaign. It is a real degradation of the arm into
                # the baseline, so it is counted and named rather than left invisible:
                # a reader who sees the frontier arm's number has to know how many of its
                # rounds were actually the baseline policy.
                self.notes.append(
                    f"POLICY_FALLBACK[{mode}->baseline] at round {round_index}: the "
                    f"post-qualification policy had no candidates, so this round ran the "
                    f"BASELINE policy. This round's measurements are not this arm's policy."
                )
                logger.warning(
                    "%s round %d: %s policy had no candidates; this round falls back to "
                    "the baseline policy.", self.name, round_index, mode,
                )
                baseline_source = f"{mode}_fallback_to_baseline"

        strategy = self._strategy(round_index, revealed)
        evidence = [
            {"variant_sequence": self._sequence_for(v), "computed_score": float(score)}
            for v, score in revealed.measured.items()
        ]
        # On a pool-restricted library the design space is the measured variants, so hand
        # the engine the untested ones instead of letting it generate Hamming neighbours
        # that nobody assayed. Everything after generation - prior, DKL, Pareto, DPP -
        # runs unchanged; only the source of candidates differs.
        candidate_pool = None
        if self.library.pool_restricted:
            from controller.gatekeeper import Variant

            candidate_pool = [Variant(notation=vid) for vid in view.untested(revealed)]
            if not candidate_pool:
                raise BenchmarkError(
                    "pool-restricted library exhausted: every eligible variant is "
                    "already measured"
                )

        # select_next_micro_batch draws its Thompson sample from numpy's global RNG; seeding it
        # here is a caller-side action that makes the production path reproducible.
        np.random.seed((self.seed * 7919 + round_index) % (2**31 - 1))
        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=self.wt_sequence,
            evidence=evidence,
            batch_size=n,
            catalytic_invariants=set(getattr(self.metadata, "catalytic_positions", []) or []),
            pdb_path=self.pdb_path,
            random_seed=(self.seed or 0) * 1000 + round_index,
            apply_auxiliary_objective_weights=self.ablation in ("aux-weighted", "legacy"),
            ablate=self.ablate_components or None,
            candidate_pool=candidate_pool,
            # Write-only diagnostics sink, attached by the harness. This is what
            # captures the engine's REAL ballot: the selection trace holds committed
            # picks plus harness rejections only, which cannot separate a narrow
            # generator from acquisition re-ranking the same shallow candidates.
            candidate_trace=self.candidate_trace,
        )
        self.last_batch_predictions = {}
        proposals: List[str] = []
        for variant, mean in zip(batch.variants, batch.predicted_means or [0.0] * len(batch.variants)):
            combo = _combo_from_notation(variant.notation, self.library)
            variant_id = _canonical_variant_id(combo, self.library) if combo else variant.notation
            self.last_batch_predictions[variant_id] = float(mean)
            proposals.append(variant_id)
        self.last_batch_source = baseline_source
        return proposals

    def observe(self, round_index: int, revealed_batch: Dict[str, float]) -> None:
        """Production's between-round feedback stage; only live-Muse runs implement one."""
        if self.use_live_muse and revealed_batch:
            try:
                self._reflect(round_index, revealed_batch)
            except Exception as exc:
                self.notes.append(f"reflection failed: {exc}")


# ==============================================================================================
# Protocol
# ==============================================================================================
@dataclass
class RoundMetrics:
    method: str
    seed: int
    round_index: int
    n_measurements: int
    best_fitness: float
    new_top1_hits: int
    cum_top1_excluding_initial: int
    cum_top1_including_initial: int
    found_top1: bool
    round_seconds: float
    cumulative_seconds: float
    model_calls: int
    token_usage: str
    requested: int
    accepted: int
    rejected: int
    status: str
    failure_reason: str
    #: How many extra rounds of candidate requests it took to fill this batch. 0 means the
    #: policy's first ranked list filled it outright.
    backfill_attempts: int = 0
    #: How many of `accepted` came from those extra requests rather than the first list.
    backfill_added: int = 0
    #: Which selection rule produced this round's measurements. "policy" is the arm's own
    #: rule; anything else names a secondary rule, so a reader can see which measurements
    #: did not come from the policy under test.
    selection_source: str = "policy"


@dataclass
class TraceRow:
    method: str
    seed: int
    round_index: int
    rank: int
    variant_id: str
    status: str
    reason: str
    revealed_fitness: Optional[float]


@dataclass
class RunOutcome:
    method: str
    seed: int
    status: str                       # completed | failed | unavailable | unsupported | skipped
    failure_reason: str = ""
    rounds: List[RoundMetrics] = field(default_factory=list)
    trace: List[TraceRow] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    model_calls: int = 0
    token_usage: Optional[int] = None


def run_single(
    selector: Selector,
    evaluator: Evaluator,
    view: CandidateView,
    initial_ids: Sequence[str],
    seed: int,
) -> RunOutcome:
    """
    Runs one method for one seed under the shared protocol.

    The initial sample is identical across methods within a seed; everything after it diverges,
    because each method learns only from the outcomes of its own choices.
    """
    outcome = RunOutcome(method=selector.name, seed=seed, status="completed")
    revealed = RevealedState()
    started = time.perf_counter()

    # Round 0: the shared starting sample.
    initial = evaluator.reveal(initial_ids)
    revealed.measured.update(initial)
    initial_hits = sum(1 for v in initial_ids if evaluator.is_top1(v))
    for rank, variant_id in enumerate(initial_ids):
        outcome.trace.append(TraceRow(
            method=selector.name, seed=seed, round_index=0, rank=rank, variant_id=variant_id,
            status="initial", reason="", revealed_fitness=initial[variant_id],
        ))
    elapsed = time.perf_counter() - started
    outcome.rounds.append(RoundMetrics(
        method=selector.name, seed=seed, round_index=0,
        n_measurements=len(revealed.measured), best_fitness=revealed.best(),
        new_top1_hits=initial_hits, cum_top1_excluding_initial=0,
        cum_top1_including_initial=initial_hits, found_top1=initial_hits > 0,
        round_seconds=round(elapsed, 4), cumulative_seconds=round(elapsed, 4),
        model_calls=0, token_usage="unknown", requested=len(initial_ids),
        accepted=len(initial_ids), rejected=0, status="completed", failure_reason="",
    ))

    # Harness-style methods (ALDE) run their own internal loop and return every round's picks
    # at once. The outcomes are still revealed round by round below, and the metrics, rejection
    # rules and trace are identical to a per-round method's.
    campaign: Optional[List[List[str]]] = None
    campaign_runner = getattr(selector, "run_campaign", None)
    if callable(campaign_runner):
        attach = getattr(selector, "attach_truth", None)
        if callable(attach):
            attach(evaluator.export_truth_table(selector.name))
        try:
            campaign = campaign_runner(view, initial_ids, view.n_rounds, view.batch_size)
        except Exception as exc:
            outcome.status = "failed"
            outcome.failure_reason = f"{type(exc).__name__}: {exc}"
            logger.error("%s seed %d campaign failed: %s", selector.name, seed, exc)
            outcome.model_calls = selector.model_calls
            outcome.notes.extend(selector.notes)
            return outcome

    cumulative_hits_excluding_initial = 0
    for round_index in range(1, view.n_rounds + 1):
        round_started = time.perf_counter()
        try:
            if campaign is not None:
                proposals = campaign[round_index - 1] if round_index - 1 < len(campaign) else []
            else:
                proposals = selector.propose(round_index, view.batch_size, view, revealed)
        except Exception as exc:
            outcome.status = "failed"
            outcome.failure_reason = f"round {round_index}: {type(exc).__name__}: {exc}"
            logger.error("%s seed %d failed at round %d: %s", selector.name, seed, round_index, exc)
            break

        accepted: List[str] = []
        seen_this_batch: Set[str] = set()
        # Every id this round has already passed judgement on, so a re-offered candidate
        # is not traced twice when the batch is backfilled below.
        judged: Set[str] = set()

        def _judge(candidates: Sequence[str]) -> int:
            """Applies the rejection rules to a proposal list. Returns how many it accepted."""
            before = len(accepted)
            for rank, variant_id in enumerate(candidates):
                if variant_id in judged:
                    # Already decided earlier this round. A longer ranked list repeats its
                    # own prefix, so this is a re-offer, not a new duplicate.
                    continue
                judged.add(variant_id)
                if len(accepted) >= view.batch_size:
                    reason = "over_batch_size"
                elif variant_id in seen_this_batch:
                    reason = "duplicate_within_batch"
                elif variant_id in revealed.measured:
                    reason = "already_observed"
                elif variant_id not in view.index_of:
                    reason = "not_in_eligible_universe"
                else:
                    reason = ""
                if reason:
                    outcome.trace.append(TraceRow(
                        method=selector.name, seed=seed, round_index=round_index, rank=rank,
                        variant_id=variant_id, status="rejected", reason=reason,
                        revealed_fitness=None,
                    ))
                    continue
                seen_this_batch.add(variant_id)
                accepted.append(variant_id)
            return len(accepted) - before

        _judge(proposals)
        first_pass_accepted = len(accepted)
        selection_source = getattr(selector, "last_batch_source", "policy")

        # Backfill a short batch from the SAME policy.
        #
        # The protocol requires exactly `batch_size` unique measurements per round, and
        # every published baseline measured all 125. Our arms did not: generated proposals
        # that are not in the eligible universe were rejected and never replaced, so GB1
        # reliability averaged 123.23 (min 120, under-filled in 42 of 48 seeds), GB1
        # frontier 124.02 (min 107) and TrpB frontier 120.23 (min 116). GB1 covers 149,360
        # of 160,000 possible 4-site combinations, so ~6.65% of free proposals do not
        # exist as measured variants. That is a protocol violation, and it biases AGAINST
        # our arms: they measured fewer variants than their competitors.
        #
        # The fill asks the policy for a LONGER ranked list and takes its next eligible
        # picks. That is the policy's own ordering continued, not a second rule: nothing is
        # padded with random draws. A campaign-style method (ALDE) precomputes every round
        # at once and is left alone - its rounds are not re-requestable, and those baselines
        # already fill their budget.
        backfill_attempts = 0
        if campaign is None and len(accepted) < view.batch_size:
            untested_remaining = len(view.variant_ids) - len(revealed.measured)
            overshoot = 0
            for backfill_attempts in range(1, BACKFILL_MAX_ATTEMPTS + 1):
                missing = view.batch_size - len(accepted)
                if missing <= 0:
                    backfill_attempts -= 1
                    break
                # Ask for the batch plus a growing margin, so the margin covers the
                # observed rejection rate with room to spare and is still bounded.
                overshoot = max(overshoot * 2, missing * 4, 8)
                request = min(view.batch_size + overshoot, max(untested_remaining, 1))
                try:
                    more = selector.propose(round_index, request, view, revealed)
                except Exception as exc:
                    logger.warning(
                        "%s seed %d round %d: backfill request for %d candidates failed "
                        "(%s: %s); the batch stays short.",
                        selector.name, seed, round_index, request, type(exc).__name__, exc,
                    )
                    break
                gained = _judge(more)
                if len(accepted) >= view.batch_size:
                    break
                if gained == 0 and not [v for v in more if v not in judged]:
                    # The policy has nothing new to offer; further attempts cannot help.
                    logger.warning(
                        "%s seed %d round %d: the policy offered no new candidates on "
                        "backfill attempt %d; the batch stays short at %d/%d.",
                        selector.name, seed, round_index, backfill_attempts,
                        len(accepted), view.batch_size,
                    )
                    break
        backfill_added = len(accepted) - first_pass_accepted
        if backfill_added:
            # The backfilled picks came from the same ranked policy list, so the source is
            # unchanged - but the fact that a fill was needed is recorded per round.
            selection_source = f"{selection_source}+backfilled({backfill_added})"

        # The batch is committed here; only now does the evaluator reveal its outcomes.
        batch_results = evaluator.reveal(accepted)
        revealed.measured.update(batch_results)
        new_hits = sum(1 for v in accepted if evaluator.is_top1(v))
        cumulative_hits_excluding_initial += new_hits
        for rank, variant_id in enumerate(accepted):
            outcome.trace.append(TraceRow(
                method=selector.name, seed=seed, round_index=round_index, rank=rank,
                variant_id=variant_id, status="accepted", reason="",
                revealed_fitness=batch_results[variant_id],
            ))

        observe = getattr(selector, "observe", None)
        if callable(observe):
            try:
                observe(round_index, batch_results)
            except Exception as exc:
                selector.notes.append(f"observe failed at round {round_index}: {exc}")

        round_seconds = time.perf_counter() - round_started
        # Counted from everything judged this round, not from the first proposal list:
        # with backfilling, `len(proposals)` is no longer the number of candidates offered.
        rejected = len(judged) - len(accepted)
        status = "completed" if len(accepted) == view.batch_size else "short_batch"
        outcome.rounds.append(RoundMetrics(
            method=selector.name, seed=seed, round_index=round_index,
            n_measurements=len(revealed.measured), best_fitness=revealed.best(),
            new_top1_hits=new_hits,
            cum_top1_excluding_initial=cumulative_hits_excluding_initial,
            cum_top1_including_initial=cumulative_hits_excluding_initial + initial_hits,
            found_top1=(cumulative_hits_excluding_initial + initial_hits) > 0,
            round_seconds=round(round_seconds, 4),
            cumulative_seconds=round(time.perf_counter() - started, 4),
            model_calls=selector.model_calls,
            token_usage="unknown" if selector.token_usage is None else str(selector.token_usage),
            requested=view.batch_size, accepted=len(accepted), rejected=rejected,
            status=status,
            failure_reason=(
                "" if status == "completed"
                else f"{rejected} selection(s) rejected; backfill exhausted after "
                     f"{backfill_attempts} attempt(s) at {len(accepted)}/{view.batch_size}"
            ),
            backfill_attempts=backfill_attempts,
            backfill_added=backfill_added,
            selection_source=selection_source,
        ))
        if status == "short_batch":
            logger.warning(
                "%s seed %d round %d: only %d/%d selections were usable AFTER %d backfill "
                "attempt(s). The arm has now spent less than the stated budget, so its "
                "result is not comparable with a baseline that measured the full batch.",
                selector.name, seed, round_index, len(accepted), view.batch_size,
                backfill_attempts,
            )

    outcome.model_calls = selector.model_calls
    outcome.token_usage = selector.token_usage
    outcome.notes.extend(selector.notes)
    return outcome


# ==============================================================================================
# Self-checks (run as part of --smoke-test)
# ==============================================================================================
class _ScriptedSelector(Selector):
    """A selector that returns a fixed list, used to prove the rejection rules actually bite."""

    name = "scripted"

    def __init__(self, seed: int, scripted: Sequence[str]) -> None:
        super().__init__(seed)
        self.scripted = list(scripted)

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        return list(self.scripted)


def run_self_checks(library: LibrarySpec) -> List[Tuple[str, bool, str]]:
    """
    Validates the benchmark's own mechanics against `library`: encoding round trips,
    genotype-length discipline, evaluator isolation, selection rejection rules and budget
    accounting. Returns (check, passed, detail) rows.
    """
    results: List[Tuple[str, bool, str]] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        results.append((name, bool(condition), detail))

    def raises(error: type, call) -> bool:
        try:
            call()
        except error:
            return True
        except Exception:
            return False
        return False

    # 1. Notation <-> combo round trip, including parent-residue slots.
    full_combo = library.parent_residues
    mutated = "".join(
        "A" if residue != "A" else "V" for residue in library.parent_residues
    )
    notation = _canonical_variant_id(mutated, library)
    check("notation round trip", _combo_from_notation(notation, library) == mutated,
          f"{mutated} -> {notation}")
    # A mutation string OMITS wild-type positions: this one names a single site and must
    # still canonicalise to a FULL-LENGTH genotype with the parent filled in elsewhere.
    partial_combo = mutated[0] + library.parent_residues[1:]
    partial = _canonical_variant_id(partial_combo, library)
    check("wild-type sites are omitted and refilled",
          _combo_from_notation(partial, library) == partial_combo
          and len(partial_combo) == library.n_sites,
          f"'{partial}' -> {partial_combo} ({library.n_sites} sites)")
    check("off-library notation rejected",
          _combo_from_notation("S65T", library) is None
          if 65 not in library.positions else True,
          "position outside the library")
    wrong_parent = f"{'A' if library.parent_residues[0] != 'A' else 'V'}{library.positions[0]}V"
    check("parent-mismatch rejected", _combo_from_notation(wrong_parent, library) is None,
          f"{wrong_parent}: wrong parent residue")
    # Unequal genotype lengths are refused, never zipped or truncated.
    check("short genotype refused, not truncated",
          raises(BenchmarkError, lambda: _canonical_variant_id(full_combo[:-1], library)),
          f"{library.n_sites - 1} of {library.n_sites} sites")
    check("over-long genotype refused",
          raises(BenchmarkError, lambda: _canonical_variant_id(full_combo + "A", library)),
          f"{library.n_sites + 1} of {library.n_sites} sites")

    universe = build_synthetic_universe(library, n_variants=600, seed=3)
    evaluator = Evaluator(universe)
    view = CandidateView(
        variant_ids=universe.variant_ids, one_hot=universe.one_hot, index_of=universe.index_of,
        combos=universe.combos, library=library, batch_size=5, n_rounds=1,
    )

    # 2. The selector-visible view carries no outcome information.
    check(
        "candidate view holds no fitness",
        not any("fitness" in field_name for field_name in vars(view)),
        f"fields: {sorted(vars(view))}",
    )

    # 3. Rejection rules: duplicate, already-observed, unknown, and over-budget selections.
    initial = universe.variant_ids[:5]
    fresh = universe.variant_ids[10:13]
    scripted = [fresh[0], fresh[0], initial[0], "V183Z", fresh[1], fresh[2]]
    outcome = run_single(_ScriptedSelector(0, scripted), evaluator, view, initial, seed=0)
    reasons = {row.reason for row in outcome.trace if row.status == "rejected"}
    accepted = [row.variant_id for row in outcome.trace if row.status == "accepted"]
    check("duplicate rejected", "duplicate_within_batch" in reasons, str(sorted(reasons)))
    check("already-observed rejected", "already_observed" in reasons, str(sorted(reasons)))
    check("unknown variant rejected", "not_in_eligible_universe" in reasons, str(sorted(reasons)))
    check("accepted only valid picks", accepted == fresh, f"accepted={accepted}")
    # A policy with nothing new to offer cannot be backfilled, so the batch stays short and
    # says so. What must NOT happen is a silent pad with random draws to hide the shortfall.
    check("unfillable short batch is reported short", outcome.rounds[-1].status == "short_batch",
          f"accepted {outcome.rounds[-1].accepted}/{outcome.rounds[-1].requested} after "
          f"{outcome.rounds[-1].backfill_attempts} backfill attempt(s)")
    check("an unfillable batch is not padded",
          outcome.rounds[-1].accepted == len(fresh),
          f"accepted {outcome.rounds[-1].accepted}, expected {len(fresh)}")

    # 4. Budget accounting: unique measurements == initial + accepted, never more.
    final = outcome.rounds[-1]
    check("budget counted from unique measurements",
          final.n_measurements == len(initial) + len(fresh),
          f"{final.n_measurements} == {len(initial)} + {len(fresh)}")

    # 5. A full-budget run consumes exactly the configured budget.
    view_full = CandidateView(
        variant_ids=universe.variant_ids, one_hot=universe.one_hot, index_of=universe.index_of,
        combos=universe.combos, library=library, batch_size=10, n_rounds=3,
    )
    full = run_single(RandomSelector(1), evaluator, view_full, universe.variant_ids[:10], seed=1)
    check("full budget consumed", full.rounds[-1].n_measurements == 10 + 3 * 10,
          f"{full.rounds[-1].n_measurements} measurements")
    check("top-1% threshold is evaluator-side",
          evaluator.top1_threshold > min(universe.fitness.values()),
          f"threshold {evaluator.top1_threshold:.4g}")

    # 6. The CLI's ablation vocabulary must match the engine's, or a name this parser
    # accepts fails inside the first round instead of at the prompt.
    try:
        from workers.math_worker import ABLATABLE_COMPONENTS as _engine_components

        check("ablation vocabulary matches the engine",
              set(ABLATABLE_COMPONENTS) == set(_engine_components),
              f"cli={sorted(ABLATABLE_COMPONENTS)} engine={sorted(_engine_components)}")
    except Exception as exc:
        check("ablation vocabulary matches the engine", False, f"{type(exc).__name__}: {exc}")

    # 7. Every legacy provenance stamp this benchmark can emit maps onto the explicit
    # vocabulary, and every member of that vocabulary is defined.
    check("legacy provenance stamps all map onto the vocabulary",
          all(isinstance(map_legacy_provenance_status(s), LabelProvenance)
              for s in LEGACY_PROVENANCE_STATUS_TO_ENUM),
          f"{len(LEGACY_PROVENANCE_STATUS_TO_ENUM)} mapped")
    check("every provenance kind is defined",
          all(m in LABEL_PROVENANCE_DEFINITIONS for m in LabelProvenance),
          f"{len(LABEL_PROVENANCE_DEFINITIONS)} definitions")

    # 8. A malformed ablation request is refused rather than running the full engine.
    refused = []
    for spec in (",", "  ", "dpp,,prior", "not_a_component"):
        try:
            parse_ablate_components(spec)
        except ConfigurationError:
            refused.append(spec)
    check("malformed --ablate-components values are refused", len(refused) == 4,
          f"refused {refused}")
    check("empty --ablate-components means no ablation", parse_ablate_components("") == (),
          "the parser default")

    # 9. The manifest's library identity is taken from the universe that was scored, not
    # from the library that was requested. With no active-library state there is nothing
    # to drift, but a manifest describing a different library from the one the universe
    # was built on would misattribute every number in the run.
    identity = _library_identity(universe.library)
    check("manifest library identity comes from the scored universe",
          identity["key"] == universe.library.key == library.key
          and identity["positions"] == list(library.positions)
          and identity["parent_residues"] == library.parent_residues,
          f"{identity['key']} / {identity['positions']} / {identity['parent_residues']}")
    return results


# ==============================================================================================
# The run manifest: the authoritative, machine-readable record of what was run
# ==============================================================================================
# Every path below must be present in a manifest. The list is the contract: a consumer can
# assert against it, and validate_manifest() refuses to write a manifest that breaks it.
#
# The split that matters: CONFIGURATION (what was run) must be complete and must contain
# nothing derived from the outcomes. The run that forced this rule reconstructed which of 62
# completed arms was which configuration by fingerprinting their RESULTS, which is
# scientifically indefensible - the identity of a configuration cannot depend on what it
# produced. Outcome-derived values live under `results`.
MANIFEST_REQUIRED_PATHS: Tuple[str, ...] = (
    "manifest_schema_version",
    "run_id",
    "configuration_fingerprint",
    "run.state",
    "run.state_reason",
    "run.abort_reason",
    "run.started_utc",
    "run.finished_utc",
    "run.timezone",
    "run.duration_seconds",
    "run.seeds_attempted",
    "run.seeds_completed",
    "run.arms_requested",
    "run.arms_completed",
    "run.failures",
    "run.warnings",
    "configuration.acquisition_policy.policy_id",
    "configuration.acquisition_policy.formula",
    "configuration.acquisition_policy.parameters",
    "configuration.components.ablatable",
    "configuration.components.enabled",
    "configuration.components.disabled",
    "configuration.components.states",
    "configuration.components.requested_spec",
    "configuration.model.surrogate",
    "configuration.model.hyperparameters",
    "configuration.dataset.name",
    "configuration.dataset.source",
    "configuration.dataset.version",
    "configuration.dataset.sha256",
    "configuration.dataset.unique_variants",
    "configuration.dataset.eligibility_rules",
    "configuration.dataset.excluded_records",
    "configuration.dataset.library",
    "configuration.dataset.library.key",
    "configuration.dataset.library.positions",
    "configuration.dataset.library.parent_residues",
    "configuration.dataset.library.labelled_stage",
    "configuration.dataset.library.pool_restricted",
    "configuration.dataset.library.selection_mode",
    "configuration.library_selection.requested",
    "configuration.library_selection.resolved",
    "configuration.assay.censoring",
    "configuration.assay.label_policy",
    "configuration.label_provenance.status",
    "configuration.label_provenance.legacy_status",
    "configuration.label_provenance.definition",
    "configuration.label_provenance.detail",
    "configuration.label_provenance.limitations",
    "configuration.label_provenance.open_blockers",
    "configuration.assay.fitness_direction",
    "configuration.assay.transformation",
    "configuration.assay.threshold_definition",
    "configuration.measurement_budget.initial_size",
    "configuration.measurement_budget.rounds",
    "configuration.measurement_budget.batch_size",
    "configuration.measurement_budget.total_unique_measurements",
    "configuration.measurement_budget.batch_schedule",
    "configuration.seeding.base_seed",
    "configuration.seeding.n_seeds",
    "configuration.seeding.seed_values",
    "configuration.seeding.initial_sets_sha256",
    "configuration.seeding.seed_block_sha256",
    "configuration.methods",
    "configuration.code.commit",
    "configuration.code.dirty_tree",
    "configuration.code.source_tree_sha256",
    "configuration.environment.python",
    "configuration.environment.platform",
    "configuration.environment.packages",
    "results.top1.threshold_value",
    "limitations",
)

# Keys that may never appear anywhere inside `configuration`, because they are outcomes.
_OUTCOME_KEY_DENYLIST = frozenset({
    "best_fitness", "mean_best_fitness", "std_best_fitness", "threshold_value", "global_best",
    "n_hits", "revealed_fitness", "found_top1", "top1_hits", "results", "summary_rows",
})


def _json_safe(value: Any) -> Any:
    """
    Converts a value into something `json.dumps` accepts without a `default=` fallback.

    A manifest that only serialises via `default=str` is not machine-readable: tuples,
    sets, Paths and enums all stringify into shapes no consumer can parse back.
    """
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = sorted(value, key=repr) if isinstance(value, (set, frozenset)) else list(value)
        return [_json_safe(v) for v in items]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return _json_safe(float(value))
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, datetime):
        return value.astimezone(timezone.utc).isoformat()
    return repr(value)


def configuration_fingerprint(configuration: Dict[str, Any]) -> str:
    """
    A stable identity for a configuration, computable from the manifest alone.

    Canonical JSON (sorted keys, no insignificant whitespace) over the whole
    `configuration` subtree. Two arms that differ in anything that was configured get
    different fingerprints; two runs of the same configuration that produced different
    results get the same one, because no outcome enters this subtree.
    """
    canonical = json.dumps(_json_safe(configuration), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def seed_block_identity(
    base_seed: int,
    n_seeds: int,
    seeds: Sequence[int],
    initial_sets: Dict[int, Sequence[str]],
    inline_limit: int = 5000,
) -> Dict[str, Any]:
    """
    Records the starting conditions precisely enough to reproduce or recognise them.

    The seed values alone are not the starting conditions: the shared initial sample is
    drawn from the eligible universe, so the same seed over a different eligibility filter
    is a different experiment. Each seed's initial set is hashed, and the whole block gets
    one hash, so two runs can be compared without anyone reading an outcome.
    """
    per_seed: Dict[str, str] = {}
    sizes: Dict[str, int] = {}
    for seed in seeds:
        ids = list(initial_sets.get(seed, ()))
        digest = hashlib.sha256()
        digest.update(f"seed={seed}\n".encode())
        for variant_id in ids:  # selection order is part of the identity
            digest.update(variant_id.encode())
            digest.update(b"\n")
        per_seed[str(seed)] = digest.hexdigest()
        sizes[str(seed)] = len(ids)

    block = hashlib.sha256()
    for seed in seeds:
        block.update(per_seed[str(seed)].encode())
    total = sum(sizes.values())
    identity: Dict[str, Any] = {
        "base_seed": base_seed,
        "n_seeds": n_seeds,
        "seed_values": list(seeds),
        "initial_set_sizes": sizes,
        "initial_sets_sha256": per_seed,
        "seed_block_sha256": block.hexdigest(),
        "sampler": "random.Random(seed).sample(eligible_variant_ids, initial_size); the "
                   "eligible list is sorted lexicographically, so the draw is reproducible",
        "shared_across_methods": True,
        "initial_sets_inlined": total <= inline_limit,
    }
    if identity["initial_sets_inlined"]:
        identity["initial_sets"] = {str(s): list(initial_sets.get(s, ())) for s in seeds}
    else:
        identity["initial_sets"] = None
        identity["initial_sets_omitted_reason"] = (
            f"{total} starting variants exceeds the {inline_limit}-id inline limit; the "
            "per-seed hashes above still identify the block"
        )
    return identity


def _code_identity(repository_state: Dict[str, Any]) -> Dict[str, Any]:
    """Commit, dirty-tree status and a content hash of the tree the run executed against."""
    git_available = bool(repository_state.get("git_available"))
    status_lines = [str(line) for line in (repository_state.get("status_lines") or [])]
    # An unknown tree state is recorded as None, never as False: 'we could not tell' must
    # not be readable as 'clean'.
    dirty: Optional[bool] = bool(status_lines) if git_available else None
    return {
        "commit": repository_state.get("head"),
        "git_available": git_available,
        "dirty_tree": dirty,
        "dirty_file_count": len(status_lines) if git_available else None,
        "dirty_files": status_lines[:250],
        "dirty_files_truncated": len(status_lines) > 250,
        "source_tree_sha256": repository_state.get("source_tree_sha256"),
        "files_hashed": repository_state.get("files_hashed"),
        "benchmark_version": BENCHMARK_VERSION,
        "benchmark_file": str(Path(__file__).resolve()),
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
    }


# The LibrarySpec attributes the manifest records. Named here so this module's whole
# read surface onto the library description is one list.
#
# It used to be read from `ACTIVE_LIBRARY` and the eight other module globals that
# `configure_library()` reassigned, and the manifest recorded those globals as evidence
# that the library had really been configured. There is no active-library state any more:
# the library is threaded explicitly and `Universe` carries the one it was built from,
# which is strictly better evidence. The manifest now attests to the library the universe
# was actually scored against, not to the last value a global happened to hold.
_LIBRARY_MANIFEST_ATTRIBUTES: Tuple[str, ...] = (
    "key", "label", "accession", "positions", "parent_residues", "dataset",
    "structure_pdb", "published_optimum", "published_measured_count",
    "published_imputed_count", "provenance_note", "labelled_stage", "pool_restricted",
    # Added by the datasets refactor; each one silently changes what was run.
    "unlabelled_placeholder", "allowed_aa_per_position", "fitness_direction",
    "assay_units", "assay_transformation", "detection_floor", "detection_ceiling",
    "replicate_count", "source_publication", "source_doi", "dataset_version",
    "duplicate_policy", "conflict_tolerance",
)


def _library_identity(library: LibrarySpec) -> Dict[str, Any]:
    """
    The library definition: which positions, which parent, which eligibility switches.

    A manifest that did not record this could not tell a TrpB arm from a GB1 arm, which
    is the multi-landscape version of the 62-arm reconstruction problem. `labelled_stage`
    and `pool_restricted` are load-bearing rather than descriptive: labelled_stage decides
    which rows are eligible at all (69% of the HER2 file is an unlabelled split), and
    pool_restricted changes where candidates come from, which makes the task strictly
    easier. `source` records how this spec reached the manifest, so a reader can tell the
    library the universe was built from apart from one that was merely requested.
    """
    identity: Dict[str, Any] = {
        attribute: _json_safe(getattr(library, attribute, None))
        for attribute in _LIBRARY_MANIFEST_ATTRIBUTES
    }
    identity["alphabet"] = AA_ALPHABET
    identity["n_sites"] = library.n_sites
    identity["one_hot_features"] = library.one_hot_width
    identity["selection_mode"] = library.selection_mode.value
    identity["metadata_json"] = str(library.metadata_json)
    identity["per_position_alphabets"] = [
        library.alphabet_for(slot) for slot in range(library.n_sites)
    ]
    identity["source"] = "Universe.library (the library this universe was built from)"
    identity["available_libraries"] = sorted(LIBRARIES)
    return identity


def _dataset_identity(universe: "Universe", args: argparse.Namespace) -> Dict[str, Any]:
    """
    Dataset name, source, version, hash and the exact eligibility rules applied.

    Name and source come from the library registry rather than a second table in this
    file, so there is one declaration per landscape. When --dataset points somewhere the
    active library does not own, that is recorded as an override instead of being
    described with the library's citation - scoring one library's variants against
    another's table must not be able to hide behind the right-looking name.
    """
    fingerprint = universe.fingerprint
    # The library the universe was BUILT from, not one read from module state: that is
    # what makes this record an attestation rather than a restatement of the request.
    library = universe.library
    owns_dataset = Path(fingerprint.path) == Path(library.dataset)

    if universe.synthetic:
        name = (
            f"SYNTHETIC additive + pairwise landscape over the "
            f"{library.n_sites}-site {library.key} geometry"
        )
        source = "generated in-process by build_synthetic_universe(); NOT biological"
        version = f"generator parameters in dataset.synthetic_generator; {fingerprint.path}"
    elif owns_dataset:
        name = library.label
        source = library.provenance_note or "UNDECLARED_BY_CALLER"
        version = (
            library.dataset_version
            or "prepared CSV as committed; content_sha256 below is the authoritative version"
        )
    else:
        name = Path(fingerprint.path).stem or "UNDECLARED"
        source = (
            f"UNDECLARED_BY_CALLER: --dataset points at a file library {library.key!r} "
            "does not declare, so that library's citation does NOT describe these labels"
        )
        version = "UNDECLARED_BY_CALLER"

    identity: Dict[str, Any] = {
        "name": name,
        "source": source,
        "version": {
            "declared": version,
            "content_sha256": fingerprint.sha256,
        },
        "path": fingerprint.path,
        "sha256": fingerprint.sha256,
        "size_bytes": fingerprint.size_bytes,
        "rows_parsed": fingerprint.rows_parsed,
        "unique_variants": fingerprint.unique_variants,
        "synthetic": bool(universe.synthetic),
        "declared_by_its_library": bool(owns_dataset),
        "library": _library_identity(library),
        "eligibility_rules": _eligibility_rules(library, args),
        "excluded_records": {
            "duplicates_dropped": universe.duplicates_dropped,
            "imputed_excluded": universe.excluded_imputed,
            "unmapped_dropped": universe.unmapped_dropped,
            "malformed_or_non_finite_dropped": universe.malformed_rows_dropped,
            # HER2's held-out split: 25,114 rows, 69% of the file. A manifest that did
            # not state this would describe the wrong universe.
            "unlabelled_excluded": universe.unlabelled_excluded,
            "subsampled_to": args.max_variants,
            "subsample_seed": args.subsample_seed if args.max_variants else None,
        },
        "published_split": {
            "measured": library.published_measured_count,
            "imputed_knn": library.published_imputed_count,
            "source": library.source_doi or library.provenance_note or "UNDECLARED",
        },
        "sequence_mapping": (
            f"variant id == position-ordered notation against the {library.accession} "
            f"parent; wild-type positions are OMITTED from the notation and refilled, so "
            f"every genotype is {library.n_sites} residue(s) long before any comparison"
        ),
    }
    if universe.synthetic:
        # The synthetic fingerprint hashes only seed and size, so the planted interaction
        # has to be stated separately or two arms on the same seed - one with a planted
        # sign-epistatic pair, one without - would share a configuration fingerprint.
        # This is the generator's ground truth, recorded AFTER the run; it is never in the
        # Universe a selector sees.
        identity["synthetic_generator"] = {
            "generator": "build_synthetic_universe",
            # Both recorded: `seed_requested` is what the CLI asked for and
            # `fingerprint_path` is what the generator actually stamped into the
            # universe. A caller that built the universe itself could desync them, and a
            # manifest must not quietly report one as the other.
            "seed_requested": getattr(args, "synthetic_seed", None),
            "fingerprint_path": fingerprint.path,
            "requested_size": getattr(args, "synthetic_size", None),
            "plant_sign_epistasis": bool(getattr(args, "plant_sign_epistasis", False)),
            "planted_interaction": _json_safe(universe.planted_interaction),
            "visibility": "ground truth for the test harness; never passed to a selector",
        }
    return identity


def _eligibility_rules(library: LibrarySpec, args: argparse.Namespace) -> List[str]:
    """The filters `load_universe` applies, including this library's own."""
    strict = not bool(getattr(args, "lenient_rows", False))
    rules = [
        (
            "a record must map onto the library's declared positions and agree with the "
            "parent residue at each of them; a record that does not is a REFUSAL naming "
            "the row, because a dataset that disagrees with the declared design is a "
            "configuration error, not noise"
            if strict else
            "--lenient-rows: a record that does not map onto the declared positions is "
            "DROPPED and counted in excluded_records.unmapped_dropped instead of refused"
        ),
        "the parent (zero-mutation) record is not a library member and is dropped",
        f"duplicate variant ids keep the first occurrence (policy "
        f"{library.duplicate_policy!r}, conflict tolerance {library.conflict_tolerance})",
        "blank, non-numeric and non-finite DMS_score values are dropped and counted",
        "records identified as imputed are dropped when a provenance list is available",
        "fewer than 1000 surviving variants aborts the run rather than benchmarking on "
        "too few records",
        "the surviving list is sorted lexicographically by variant id, which fixes the "
        "row order of the one-hot matrix and every tie-break",
    ]
    if library.labelled_stage is not None:
        rules.insert(0, (
            f"ONLY rows whose 'stage' column equals {library.labelled_stage!r} are "
            f"eligible: the others are an unlabelled split whose DMS_score is an absent "
            f"label, not a measurement. Each one is checked against the declared "
            f"placeholder {library.unlabelled_placeholder!r}, so a real measurement "
            f"cannot be discarded as an absent label, and admitting them would put an "
            f"invented value into every quantile"
        ))
    if library.pool_restricted:
        rules.append(
            "this library is POOL-RESTRICTED: it is a sparse sample of its combinatorial "
            "product, so the selector is handed the enumerated untested measured variants "
            "to rank rather than generating Hamming neighbours nobody assayed. Ranking an "
            "enumerated pool is strictly easier than generating under constraints"
        )
    if library.detection_floor is not None or library.detection_ceiling is not None:
        rules.append(
            f"censoring is RECORDED, never applied: a value at the declared floor "
            f"({library.detection_floor}) or ceiling ({library.detection_ceiling}) keeps "
            f"its measured value and is counted, so a block of ties is legible as a "
            f"detection limit rather than a real plateau"
        )
    return rules


def _label_provenance_block(universe: "Universe") -> Dict[str, Any]:
    """
    The explicit provenance claim for this run's labels, plus the limitations it carries.

    Three cases, kept apart because they are different claims. A positive published count
    (TrpB: 871) names it. A count of None (GB1: the prepared CSV records no split) names
    the library's own provenance note, because 'up to None imputed values' would read as
    a measured count. A recorded 0 (HER2) is a fact rather than an unknown, so it does not
    get 'up to 0 may be included' - but it is still UNVERIFIED, because the 0 is our
    registry's declaration and not an upstream flag, and laundering our own assertion into
    a provenance claim is exactly what this vocabulary exists to prevent.
    """
    library = universe.library
    imputed = library.published_imputed_count
    resolved = map_legacy_provenance_status(universe.provenance_status)
    limitations: List[str] = []
    if resolved is LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE:
        if imputed:
            limitations.append(
                f"Up to {imputed} published KNN-imputed values are present in the eligible "
                "universe and cannot be located, so no individual label in this run can be "
                "claimed as measured."
            )
        elif imputed == 0:
            limitations.append(
                f"Library {library.key!r} declares 0 imputed records, but that is THIS "
                "repository's declaration and not an upstream per-record flag, so the "
                "labels remain HISTORICAL_UNKNOWN_PROVENANCE. Verify against the source "
                "publication before describing them as measured. "
                + (library.provenance_note or "")
            )
        else:
            limitations.append(
                f"The measured/imputed split of library {library.key!r} is unverified and no "
                "published count is on record, so no individual label in this run can be "
                "claimed as measured. " + (library.provenance_note or "")
            )
        limitations.append(DIFFERENTIAL_BIAS_WARNING)
    if resolved is LabelProvenance.SIMULATED_NOISY_OUTCOME:
        limitations.append(
            "Labels are generated, not measured. Results validate mechanism only and say "
            "nothing about biological performance."
        )
    if resolved is LabelProvenance.HISTORICAL_MEASURED:
        limitations.append(
            "Labels are retrospective measurements from a published landscape. They support "
            "no prospective claim; only PROSPECTIVE_PHYSICAL_MEASUREMENT does."
        )
    return {
        "status": resolved.value,
        "legacy_status": universe.provenance_status,
        "library": library.key,
        "library_provenance_note": library.provenance_note or "",
        "published_imputed_count": imputed,
        "source_publication": library.source_publication,
        "source_doi": library.source_doi,
        "definition": LABEL_PROVENANCE_DEFINITIONS[resolved],
        "detail": universe.provenance_detail,
        "vocabulary": {m.value: LABEL_PROVENANCE_DEFINITIONS[m] for m in LabelProvenance},
        "provenance_mode_requested": None,  # filled by build_manifest from the CLI
        "limitations": limitations,
        "differential_bias": DIFFERENTIAL_BIAS_WARNING,
        "open_blockers": list(PROVENANCE_OPEN_BLOCKERS),
        "prospective": False,
    }


def _production_selection_identity(
    components_disabled: Sequence[str],
    ablation_arm: str,
    prior_spread: Optional[float],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Reads the production selector's effective acquisition policy and hyperparameters.

    Defaults are read from the live signature of `select_next_micro_batch` and the live
    module globals rather than restated here, so the manifest records what the run actually
    used instead of what this file once believed. The benchmark's own overrides are listed
    separately, because 'the default, except for these three arguments' is the only honest
    description of the call.
    """
    policy: Dict[str, Any] = {
        "policy_id": "hybrid_thompson",
        "formula": "acq(x) = EI(x; xi) + 0.3 * max(thompson_sample(x) - f_best, 0), where "
                   "thompson_sample(x) ~ N(mu(x), sigma(x)^2) and mu(x) = m(x) + mu_res(x)",
        "surrogate_mean": "m(x) = calibrated biophysical prior; the GP is fitted on residuals "
                          "y - m(x)",
        "batch_construction": "k-DPP diversity selection over the Pareto-filtered candidate "
                              "pool; candidates are enumerated from the strategy's clusters "
                              "and lead recombinants, capped at 3 mutations",
        "ablation_arm": ablation_arm,
        "components_disabled": list(components_disabled),
        "parameters": {},
        "engine": {},
    }
    model: Dict[str, Any] = {
        "surrogate": "workers.math_worker deep-kernel/GP on epistatic one-hot features with a "
                     "residual biophysical prior mean",
        "hyperparameters": {},
        "baselines": "per-method hyperparameters are recorded under configuration.methods",
    }

    try:
        import workers.math_worker as _mw
    except Exception as exc:  # pragma: no cover - only when the engine cannot be imported
        policy["engine"] = {"available": False, "error": f"{type(exc).__name__}: {exc}"}
        model["hyperparameters"] = {"available": False}
        return policy, model

    signature = inspect.signature(_mw.select_next_micro_batch)
    defaults: Dict[str, Any] = {}
    for name in ("exploration_xi", "acquisition_mode", "fantasy_variance",
                 "apply_auxiliary_objective_weights", "pareto_weights_config"):
        parameter = signature.parameters.get(name)
        if parameter is not None and parameter.default is not inspect.Parameter.empty:
            defaults[name] = _json_safe(parameter.default)
    prior_fn = signature.parameters.get("prior_mean_fn")
    defaults["prior_mean_fn"] = (
        getattr(prior_fn.default, "__qualname__", None) if prior_fn is not None else None
    )

    engine_file = Path(getattr(_mw, "__file__", "") or "")
    engine_sha = None
    if engine_file.is_file():
        engine_sha = _sha256_file(engine_file)[0]

    policy["policy_id"] = str(defaults.get("acquisition_mode") or "hybrid_thompson")
    policy["parameters"] = {
        "acquisition_mode": defaults.get("acquisition_mode"),
        "exploration_xi": defaults.get("exploration_xi"),
        "thompson_weight": 0.3,
        "fantasy_variance": defaults.get("fantasy_variance"),
        "pareto_weights_config": defaults.get("pareto_weights_config"),
        "apply_auxiliary_objective_weights": ablation_arm in ("aux-weighted", "legacy"),
        "prior_calibration_enabled": bool(_mw.PRIOR_CALIBRATION_ENABLED),
        "prior_fitness_spread": float(_mw.PRIOR_FITNESS_SPREAD),
        "prior_spread_cli_override": prior_spread,
        "prior_mean_fn": None if "prior" in set(components_disabled) else defaults["prior_mean_fn"],
    }
    policy["engine"] = {
        "available": True,
        "module": "workers.math_worker.select_next_micro_batch",
        "module_sha256": engine_sha,
        "signature_defaults": defaults,
        "overrides_passed_by_the_benchmark": {
            "batch_size": "--batch-size",
            "random_seed": "seed * 1000 + round_index",
            "numpy_global_seed": "(seed * 7919 + round_index) % (2**31 - 1), seeded by the "
                                 "caller because the Thompson draw uses numpy's global RNG",
            "catalytic_invariants": "metadata.catalytic_positions",
            "pdb_path": "workers/data/8EH0.pdb when present",
            "apply_auxiliary_objective_weights": ablation_arm in ("aux-weighted", "legacy"),
            "ablate": list(components_disabled),
        },
    }
    model["hyperparameters"] = dict(policy["parameters"])
    model["module_sha256"] = engine_sha
    return policy, model


def _assay_block(universe: "Universe") -> Dict[str, Any]:
    """
    How a raw assay number becomes the quantity being optimised, and how 'top' is defined.

    Read from the library's declared `AssayMeta` rather than asserted here: units,
    transformation, detection limits and replicate count are per-landscape facts, and the
    old module-level FITNESS_DIRECTION made one landscape's claim look universal.

    The threshold DEFINITION is configuration; the threshold VALUE is computed from the
    outcomes and therefore lives under `results`. Censoring COUNTS are outcome-adjacent
    but are properties of the eligible set, not of any method's behaviour, so they belong
    with the dataset: they say what kind of labels the run was scored against.
    """
    library = universe.library
    n = max(len(universe), 1)
    return {
        "fitness_direction": library.fitness_direction,
        "units": library.assay_units,
        "transformation": (
            f"declared {library.assay_transformation!r}; this benchmark applies none of its "
            "own - the prepared table's DMS_score column is the objective verbatim, never "
            "re-scaled, log-transformed, standardised or clipped here"
        ),
        "replicates": (
            f"{library.replicate_count} declared"
            if library.replicate_count is not None else
            "unknown: the prepared table carries one value per variant, so no replicate "
            "aggregation or measurement-noise model is applied"
        ),
        "censoring": {
            "detection_floor": library.detection_floor,
            "detection_ceiling": library.detection_ceiling,
            "at_floor": universe.censored_at_floor,
            "at_ceiling": universe.censored_at_ceiling,
            "fraction_at_floor": round(universe.censored_at_floor / n, 6),
            "policy": (
                "values at a declared limit keep their measured value and are COUNTED, not "
                "adjusted or dropped. A mass of ties at the floor is a detection limit, not "
                "a plateau, and a model treating them as exact continuous labels is "
                "mis-specified"
            ),
        },
        "label_policy": {
            "labelled_stage": library.labelled_stage,
            "unlabelled_placeholder": library.unlabelled_placeholder,
            "unlabelled_excluded": universe.unlabelled_excluded,
            "policy": (
                "absent label, measured zero and detection floor are three distinct states; "
                "only labelled rows enter the eligible universe, and each excluded row is "
                "checked against the declared placeholder so a real measurement cannot be "
                "discarded as an absent label"
            ),
        },
        "source_publication": library.source_publication,
        "source_doi": library.source_doi,
        "threshold_definition": {
            "statistic": "percentile of measured fitness over the eligible set of this run",
            "percentile": 99.0,
            "tie_handling": "fitness >= threshold counts as a hit; all ties count, so the hit "
                            "set can exceed 1% of the pool",
            "scope": "recomputed per run over the eligible set; a different eligibility "
                     "filter changes the threshold",
            "visibility": "evaluator-only; never passed to a selector or into candidate "
                          "metadata",
        },
    }


def _method_block(
    methods: Sequence[str],
    method_status: Dict[str, Dict[str, Any]],
    method_configs: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    block: Dict[str, Any] = {}
    for method in methods:
        status = dict(method_status.get(method, {}))
        probe = status.pop("probe", None)
        entry: Dict[str, Any] = {
            "status": status.get("status", "unknown"),
            "status_detail": status.get("detail", ""),
            "config": _json_safe(method_configs.get(method, {})),
        }
        if probe is not None:
            try:
                entry["availability_probe"] = json.loads(probe)
            except (TypeError, ValueError):
                entry["availability_probe"] = {"raw": str(probe)}
        block[method] = entry
    return block


def _per_arm_results(
    methods: Sequence[str],
    outcomes: Sequence[RunOutcome],
) -> Dict[str, Any]:
    """
    Outcome-derived aggregates, so a committed manifest carries the run's numbers too.

    Strictly under `results`: nothing here may be read back as part of a configuration's
    identity. summary.csv remains the full table.
    """
    per_arm: Dict[str, Any] = {}
    for method in methods:
        runs = [o for o in outcomes if o.method == method]
        finals = [o.rounds[-1] for o in runs if o.status == "completed" and o.rounds]
        per_arm[method] = {
            "seeds_attempted": len(runs),
            "seeds_completed": len(finals),
            "mean_best_fitness": (
                statistics.fmean(m.best_fitness for m in finals) if finals else None),
            "mean_top1_hits_including_initial": (
                statistics.fmean(m.cum_top1_including_initial for m in finals) if finals else None),
            "mean_unique_measurements": (
                statistics.fmean(m.n_measurements for m in finals) if finals else None),
            "total_model_calls": sum(o.model_calls for o in runs),
        }
    return per_arm


def build_manifest(
    *,
    run_id: str,
    args: argparse.Namespace,
    universe: "Universe",
    evaluator: Optional["Evaluator"],
    methods: Sequence[str],
    method_status: Dict[str, Dict[str, Any]],
    method_configs: Dict[str, Dict[str, Any]],
    components_disabled: Sequence[str],
    seeds: Sequence[int],
    initial_sets: Dict[int, Sequence[str]],
    outcomes: Sequence[RunOutcome],
    started: datetime,
    finished: datetime,
    repository_state_before: Dict[str, Any],
    repository_state_after: Optional[Dict[str, Any]] = None,
    self_checks: Sequence[Tuple[str, bool, str]] = (),
    extra_warnings: Sequence[str] = (),
    run_state: Optional[str] = None,
    abort_reason: str = "",
    limitations: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Builds the run manifest: the single authoritative, machine-readable record of a run.

    Two rules hold this together:
      1. `configuration` is complete enough to identify and re-run the arm on its own. The
         console transcript, the CLI string and the report are convenience, not record.
      2. `configuration` contains no outcome. Identifying an arm by what it produced is the
         practice this replaces, so the fingerprint deliberately cannot see results.
    """
    disabled = tuple(components_disabled)
    policy, model = _production_selection_identity(
        components_disabled=disabled,
        ablation_arm=str(getattr(args, "ablation", "none")),
        prior_spread=getattr(args, "prior_spread", None),
    )
    states = component_states(disabled)
    provenance = _label_provenance_block(universe)
    provenance["provenance_mode_requested"] = getattr(args, "provenance_mode", None)

    total_budget = args.initial_size + args.rounds * args.batch_size
    schedule = [{
        "round": 0,
        "n_measurements": args.initial_size,
        "kind": "shared_initial_sample",
        "selected_by": "random.Random(seed).sample over the eligible set, identical for every "
                       "method within a seed",
    }]
    for index in range(1, args.rounds + 1):
        schedule.append({
            "round": index,
            "n_measurements": args.batch_size,
            "kind": "selected",
            "selected_by": "the method under test; committed before any outcome is revealed",
        })

    configuration: Dict[str, Any] = {
        "benchmark_version": BENCHMARK_VERSION,
        "smoke_test": bool(getattr(args, "smoke_test", False)),
        "acquisition_policy": policy,
        "components": {
            "ablatable": list(ABLATABLE_COMPONENTS),
            "enabled": [name for name in ABLATABLE_COMPONENTS if states[name]],
            "disabled": list(disabled),
            "states": states,
            "requested_spec": getattr(args, "ablate_components", ""),
            "prior_ablation_arm": str(getattr(args, "ablation", "none")),
            "prior_spread_override": getattr(args, "prior_spread", None),
            # The qualification bar, in full. A run whose manifest does not state its
            # qualification bar is not admissible evidence, so all four facts are recorded:
            # what was asked for in each form, what it resolved to, and which form won.
            #
            # `qualify_threshold` is retained under its existing key for the readers and
            # greps that already use it, but it now carries the RESOLVED ABSOLUTE the run
            # was actually gated by, not the raw argument - which is None on a default run
            # and would have recorded "no bar" for every run that had one.
            "qualify_threshold": (
                float(evaluator.qualify_threshold_absolute) if evaluator is not None
                else getattr(args, "qualify_threshold", None)
            ),
            "qualify_threshold_absolute": (
                float(evaluator.qualify_threshold_absolute) if evaluator is not None else None
            ),
            "qualify_threshold_source": (
                evaluator.qualify_threshold_source if evaluator is not None else None
            ),
            "qualify_quantile": (
                evaluator.qualify_quantile if evaluator is not None else None
            ),
            "qualify_quantile_requested": getattr(args, "qualify_quantile", None),
            "qualify_threshold_requested": getattr(args, "qualify_threshold", None),
            # How attainable that bar is on THIS landscape. A quantile transfers the
            # rarity of the bar, not its attainability: a count of 1 means qualification
            # is equivalent to finding the global best.
            "qualify_qualifying_count": (
                int(evaluator.n_qualifying) if evaluator is not None else None
            ),
            "qualify_bar_below_top1": (
                bool(evaluator.qualify_threshold_absolute < evaluator.top1_threshold)
                if evaluator is not None else None
            ),
            "qualify_minimum_set_floor": {
                "enabled": False,
                "note": "No minimum-qualifying-set floor is applied. At the default "
                        "quantile the qualifying set is ONE variant on the CR6261 "
                        "landscapes, which makes qualification equivalent to finding the "
                        "global best; see qualify_qualifying_count. Applying a floor is a "
                        "protocol decision that must be taken in advance, not resolved "
                        "silently by the benchmark.",
            },
        },
        "model": model,
        "dataset": _dataset_identity(universe, args),
        "label_provenance": provenance,
        "assay": _assay_block(universe),
        "library_selection": {
            # Which library was asked for, which one the universe was built from, and
            # whether --dataset overrode its declared file. All three, because a mismatch
            # between them is the failure this records: one library's variants scored
            # against another library's table.
            "requested": getattr(args, "library", None),
            "resolved": universe.library.key,
            "dataset_override": (
                None if Path(str(getattr(args, "dataset", "") or "")) == Path(universe.library.dataset)
                else str(getattr(args, "dataset", "") or "")
            ),
            "lenient_rows": bool(getattr(args, "lenient_rows", False)),
        },
        "measurement_budget": {
            "initial_size": args.initial_size,
            "rounds": args.rounds,
            "batch_size": args.batch_size,
            "total_unique_measurements": total_budget,
            "batch_schedule": schedule,
            "accounting": "unique variants only; duplicate, already-observed, unknown and "
                          "over-budget selections are rejected, and a short batch is then "
                          "BACKFILLED from the same policy - the selector is asked for a "
                          "longer ranked list and its next eligible picks are taken. A "
                          "short batch is never padded with random draws. A batch that "
                          "still cannot fill is reported short, and the round metrics "
                          "carry backfill_attempts, backfill_added and selection_source "
                          "so a reader can see which measurements came from where",
            "backfill": {
                "enabled": True,
                "max_attempts": BACKFILL_MAX_ATTEMPTS,
                "source": "the arm's own policy, asked for a longer ranked list",
                "secondary_rule": None,
                "note": "Before this was implemented every baseline measured exactly 125 "
                        "unique variants while our arms did not (GB1 reliability mean "
                        "123.23 / min 120, GB1 frontier mean 124.02 / min 107, TrpB "
                        "frontier mean 120.23 / min 116), because proposals outside the "
                        "eligible universe were rejected and never replaced. That "
                        "under-filling biased AGAINST our arms. Per-round fill is recorded "
                        "in round_metrics.csv.",
            },
            "max_model_calls": getattr(args, "max_model_calls", None),
            "live_calls_permitted": bool(getattr(args, "enable_live_muse", False)),
        },
        "seeding": seed_block_identity(
            base_seed=args.base_seed,
            n_seeds=args.seeds,
            seeds=list(seeds),
            initial_sets=initial_sets,
        ),
        "methods": _method_block(methods, method_status, method_configs),
        "methods_requested": list(methods),
        "code": _code_identity(repository_state_before),
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "packages": {"numpy": np.__version__, **_package_versions()},
        },
    }

    # ---- Outcome-derived material. Nothing below may be referenced by the fingerprint. ----
    completed = [o for o in outcomes if o.status == "completed" and o.rounds]
    failures = [
        {"method": o.method, "seed": o.seed, "status": o.status, "reason": o.failure_reason}
        for o in outcomes if o.status != "completed"
    ]
    short_batches = sum(
        1 for o in outcomes for m in o.rounds if m.status == "short_batch"
    )

    warnings: List[str] = [str(w) for w in extra_warnings if str(w).strip()]
    for failure in failures:
        warnings.append(
            f"{failure['method']} seed {failure['seed']}: {failure['status']} - "
            f"{failure['reason']}"
        )
    if short_batches:
        warnings.append(
            f"{short_batches} round(s) committed a short batch: selections were rejected and "
            "the policy could not supply enough eligible replacements within "
            f"{BACKFILL_MAX_ATTEMPTS} backfill attempts, so those arms spent less than the "
            "stated budget."
        )
    for method in methods:
        attempted = [o for o in outcomes if o.method == method]
        if attempted and not [o for o in attempted if o.status == "completed" and o.rounds]:
            warnings.append(
                f"{method}: 0 of {len(attempted)} seeds completed - the arm did not run, which "
                "is not a null result."
            )
    excluded = configuration["dataset"]["excluded_records"]
    dropped = {
        key: excluded.get(key, 0) for key in
        ("duplicates_dropped", "imputed_excluded", "unmapped_dropped",
         "malformed_or_non_finite_dropped")
        if excluded.get(key, 0)
    }
    if dropped:
        warnings.append(
            f"records excluded while building the eligible universe: {dropped}. The eligible "
            "set defines the top-1% threshold, so an exclusion changes the scoring."
        )
    if excluded.get("subsampled_to"):
        warnings.append(
            f"the candidate universe was subsampled to {excluded['subsampled_to']} variants "
            f"(seed {excluded.get('subsample_seed')}): this is not the full published landscape."
        )
    if provenance["status"] != LabelProvenance.HISTORICAL_MEASURED.value:
        warnings.append(
            f"label provenance is {provenance['status']}: see "
            "configuration.label_provenance.limitations before quoting any number."
        )
    if evaluator is not None and evaluator.harness_exports:
        warnings.append(
            "the evaluator exported its complete fitness table to harness method(s) "
            f"{sorted(set(evaluator.harness_exports))}; see Evaluator.export_truth_table."
        )

    if run_state is None:
        if not outcomes:
            run_state = "FAILED"
        elif failures or short_batches:
            run_state = "PARTIAL"
        else:
            run_state = "SUCCESS"
    if run_state == "SUCCESS" and not completed:
        run_state = "FAILED"

    # PARTIAL must not be readable as 'crashed'. On an incompletely-saturated library
    # proposals outside the eligible universe are rejected; a short batch now means the
    # policy could not supply enough eligible replacements even after backfilling, and the
    # arm did then spend less than the stated budget, which a reader has to know. State the
    # reason rather than leaving one word to carry both meanings.
    if run_state == "ABORTED":
        state_reason = abort_reason or "aborted before any measurement"
    elif run_state == "FAILED":
        state_reason = (
            f"no run completed: {len(outcomes)} (method, seed) attempt(s), "
            f"{len(failures)} failed"
        )
    elif run_state == "PARTIAL":
        parts = []
        if failures:
            parts.append(f"{len(failures)} (method, seed) run(s) failed")
        if short_batches:
            parts.append(
                f"{short_batches} round(s) committed a short batch - rejected selections "
                f"are backfilled from the same policy, but it could not supply enough "
                f"eligible replacements within {BACKFILL_MAX_ATTEMPTS} attempts, so those "
                f"arms spent less than the stated budget. A short batch is never padded "
                f"with random draws. This can happen on a library that does not saturate "
                f"its combinatorial product"
            )
        state_reason = "; ".join(parts)
    else:
        state_reason = (
            "every requested arm completed every seed, and every round committed a full batch"
        )

    manifest: Dict[str, Any] = {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": run_id,
        "run": {
            "state": run_state,
            "state_reason": state_reason,
            "state_vocabulary": ["SUCCESS", "PARTIAL", "FAILED", "ABORTED"],
            "abort_reason": abort_reason,
            "started_utc": started.astimezone(timezone.utc).isoformat(),
            "finished_utc": finished.astimezone(timezone.utc).isoformat(),
            "timezone": "UTC",
            "started_epoch_seconds": started.timestamp(),
            "finished_epoch_seconds": finished.timestamp(),
            "duration_seconds": round((finished - started).total_seconds(), 3),
            "arms_requested": list(methods),
            "arms_completed": sorted({o.method for o in completed}),
            "seeds_attempted": len({(o.method, o.seed) for o in outcomes}),
            "seeds_completed": len(completed),
            "failures": failures,
            "warnings": warnings,
            "self_checks": [{"check": n, "passed": bool(p), "detail": d} for n, p, d in self_checks],
            "repository_unchanged": (
                None if repository_state_after is None else
                (repository_state_before.get("source_tree_sha256")
                 == repository_state_after.get("source_tree_sha256")
                 and repository_state_before.get("status_lines")
                 == repository_state_after.get("status_lines"))
            ),
            "repository_state_before": _json_safe(repository_state_before),
            "repository_state_after": _json_safe(repository_state_after),
        },
        "invocation": {
            "argv": list(sys.argv),
            "cwd": str(Path.cwd()),
            "resolved_args": {k: _json_safe(v) for k, v in sorted(vars(args).items())},
        },
        "configuration": _json_safe(configuration),
        "results": {
            "top1": {
                "threshold_value": float(evaluator.top1_threshold) if evaluator else None,
                "n_hits": int(evaluator.n_top1) if evaluator else None,
                "global_best": float(evaluator.global_best) if evaluator else None,
                "definition": "see configuration.assay.threshold_definition",
            },
            "per_arm": _per_arm_results(methods, outcomes),
            "notes": [f"{o.method}: {note}" for o in outcomes for note in o.notes],
        },
        "limitations": list(limitations or []),
    }
    manifest["configuration_fingerprint"] = configuration_fingerprint(manifest["configuration"])
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Dict[str, Any]) -> None:
    """
    Refuses a manifest that cannot identify its own run.

    Called on every manifest before it is written, so an incomplete record fails the run
    rather than being discovered months later by someone trying to reconstruct an arm.
    """
    missing: List[str] = []
    for path in MANIFEST_REQUIRED_PATHS:
        node: Any = manifest
        for key in path.split("."):
            if not isinstance(node, dict) or key not in node:
                missing.append(path)
                break
            node = node[key]
    if missing:
        raise ManifestIncompleteError(
            "the run manifest is missing required field(s): " + ", ".join(missing)
        )

    leaked: List[str] = []

    def _walk(node: Any, trail: str) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key in _OUTCOME_KEY_DENYLIST:
                    leaked.append(f"{trail}.{key}")
                _walk(value, f"{trail}.{key}")
        elif isinstance(node, list):
            for index, item in enumerate(node):
                _walk(item, f"{trail}[{index}]")

    _walk(manifest.get("configuration", {}), "configuration")
    if leaked:
        raise ManifestIncompleteError(
            "outcome-derived field(s) leaked into the manifest's configuration subtree, which "
            "would make a configuration's identity depend on its results: " + ", ".join(leaked)
        )


# ==============================================================================================
# Outputs
# ==============================================================================================
def _write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_trace_csv(path: Path, outcomes: Sequence[RunOutcome], stamp: str) -> None:
    _write_csv(
        path,
        ["provenance", "method", "seed", "round", "rank", "variant_id", "status", "reason", "revealed_fitness"],
        [
            [stamp, r.method, r.seed, r.round_index, r.rank, r.variant_id, r.status, r.reason,
             "" if r.revealed_fitness is None else f"{r.revealed_fitness:.9g}"]
            for outcome in outcomes for r in outcome.trace
        ],
    )


def write_round_metrics_csv(path: Path, outcomes: Sequence[RunOutcome], stamp: str) -> None:
    _write_csv(
        path,
        ["provenance", "method", "seed", "round", "n_unique_measurements", "best_measured_fitness",
         "new_top1_hits", "cum_top1_hits_excluding_initial", "cum_top1_hits_including_initial",
         "found_top1", "round_seconds", "cumulative_seconds", "model_calls", "token_usage",
         "requested", "accepted", "rejected", "status", "failure_reason",
         # Budget fill and attribution, per round: how hard the batch was to fill, and
         # which selection rule the measurements came from.
         "backfill_attempts", "backfill_added", "selection_source"],
        [
            [stamp, m.method, m.seed, m.round_index, m.n_measurements, f"{m.best_fitness:.9g}",
             m.new_top1_hits, m.cum_top1_excluding_initial, m.cum_top1_including_initial,
             int(m.found_top1), m.round_seconds, m.cumulative_seconds, m.model_calls, m.token_usage,
             m.requested, m.accepted, m.rejected, m.status, m.failure_reason,
             m.backfill_attempts, m.backfill_added, m.selection_source]
            for outcome in outcomes for m in outcome.rounds
        ],
    )


def summarise(outcomes: Sequence[RunOutcome], method_status: Dict[str, Dict[str, str]]) -> List[Dict[str, Any]]:
    """Aggregates across seeds. Unavailable and failed methods stay visible with zero successes."""
    by_method: Dict[str, List[RunOutcome]] = {}
    for outcome in outcomes:
        by_method.setdefault(outcome.method, []).append(outcome)

    rows: List[Dict[str, Any]] = []
    for method, status in method_status.items():
        runs = by_method.get(method, [])
        completed = [r for r in runs if r.status == "completed" and r.rounds]
        finals = [r.rounds[-1] for r in completed]
        best = [m.best_fitness for m in finals]
        hits_incl = [m.cum_top1_including_initial for m in finals]
        hits_excl = [m.cum_top1_excluding_initial for m in finals]
        measurements = [m.n_measurements for m in finals]
        runtimes = [m.cumulative_seconds for m in finals]

        def _mean(values: Sequence[float]) -> str:
            return f"{statistics.fmean(values):.6g}" if values else ""

        def _stdev(values: Sequence[float]) -> str:
            return f"{statistics.stdev(values):.6g}" if len(values) > 1 else ("0" if values else "")

        # An arm that attempted seeds and completed none did not 'run and find
        # nothing' - it failed to run. Reporting it as available with blank metrics
        # is indistinguishable from a real null result, so say so in the status.
        arm_status = status.get("status", "unknown")
        arm_detail = status.get("detail", "")
        if runs and not completed:
            arm_status = "FAILED"
            arm_detail = (
                f"0 of {len(runs)} seeds completed - the arm did not run. "
                f"See benchmark.log for the per-seed reason. "
            ) + arm_detail
            logger.error(
                "%s: 0 of %d seeds completed; reporting status=FAILED rather than a "
                "null result.", method, len(runs),
            )

        rows.append({
            "method": method,
            "status": arm_status,
            "status_detail": arm_detail,
            "seeds_attempted": len(runs),
            "seeds_completed": len(completed),
            "completion_rate": f"{len(completed) / len(runs):.3f}" if runs else "0.000",
            "mean_measurements_used": _mean(measurements),
            "mean_best_fitness": _mean(best),
            "std_best_fitness": _stdev(best),
            "mean_top1_hits_including_initial": _mean(hits_incl),
            "std_top1_hits_including_initial": _stdev(hits_incl),
            "mean_top1_hits_excluding_initial": _mean(hits_excl),
            "std_top1_hits_excluding_initial": _stdev(hits_excl),
            "success_rate_found_top1": f"{sum(1 for m in finals if m.found_top1) / len(finals):.3f}" if finals else "",
            "mean_runtime_seconds": _mean(runtimes),
            "total_model_calls": sum(r.model_calls for r in runs),
            "token_usage": "unknown",
            "estimated_cost_usd": "unknown",
        })
    return rows


def write_summary_csv(path: Path, rows: Sequence[Dict[str, Any]], stamp: str) -> None:
    if not rows:
        return
    header = ["provenance"] + list(rows[0].keys())
    _write_csv(path, header, [[stamp] + [row[k] for k in rows[0]] for row in rows])


def _plot_title(synthetic: bool, library_label: str) -> str:
    """
    Titles the figure from the library it was run on.

    The label is an argument, never module state: a GB1 or HER2 figure titled 'TrpB
    four-site landscape' is a mislabelled figure, and a figure outlives the run folder
    that would have explained it. A synthetic figure names the geometry it was generated
    over, because 'SYNTHETIC' alone does not say which library's sites were used.
    """
    if synthetic:
        return (
            f"SYNTHETIC landscape on the {library_label} geometry (NOT biological) - "
            "best value vs budget"
        )
    return f"{library_label} - best measured fitness vs budget"


def write_plot(
    path: Path,
    outcomes: Sequence[RunOutcome],
    stamp: str,
    synthetic: bool,
    cache_dir: Path,
    library_label: str,
) -> Optional[str]:
    """
    Best measured fitness vs total measurements, mean with a +/-1 s.d. band across seeds.

    `library_label` is required rather than defaulted: the title used to be the literal
    string "TrpB four-site landscape", so a GB1 or HER2 run produced a plot captioned
    with the wrong protein.
    """
    os.environ["MPLCONFIGDIR"] = str(cache_dir)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return f"plot skipped: matplotlib unavailable ({exc})"

    series: Dict[str, Dict[int, List[float]]] = {}
    for outcome in outcomes:
        if outcome.status != "completed":
            continue
        for metric in outcome.rounds:
            series.setdefault(outcome.method, {}).setdefault(metric.n_measurements, []).append(metric.best_fitness)
    if not series:
        return "plot skipped: no completed runs"

    figure, axes = plt.subplots(figsize=(9, 5.5))
    for method in sorted(series):
        points = sorted(series[method])
        means = np.array([statistics.fmean(series[method][p]) for p in points])
        devs = np.array([statistics.stdev(series[method][p]) if len(series[method][p]) > 1 else 0.0 for p in points])
        axes.plot(points, means, marker="o", label=method)
        axes.fill_between(points, means - devs, means + devs, alpha=0.15)

    axes.set_xlabel("Unique experimental measurements used")
    axes.set_ylabel("Best measured fitness found")
    axes.set_title(_plot_title(synthetic, library_label))
    axes.legend()
    axes.grid(alpha=0.25)
    figure.text(0.01, 0.01, stamp, fontsize=7, color="#666666")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return None


def write_report(
    path: Path,
    manifest: Dict[str, Any],
    summary_rows: Sequence[Dict[str, Any]],
    outcomes: Sequence[RunOutcome],
    stamp: str,
    plot_note: Optional[str],
) -> None:
    configuration = manifest["configuration"]
    dataset = configuration["dataset"]
    budget = configuration["measurement_budget"]
    components = configuration["components"]
    top1 = manifest["results"]["top1"]
    run = manifest["run"]

    lines: List[str] = []
    add = lines.append
    add("PanoptesBio isolated comparison benchmark")
    add("=" * 78)
    add(f"Run id:        {manifest['run_id']}")
    add(f"Run state:     {run['state']} - {run['state_reason']}")
    add(f"Started (UTC): {run['started_utc']}")
    add(f"Ended (UTC):   {run['finished_utc']}")
    add(f"Provenance:    {stamp}")
    add(f"Library:       {dataset['library'].get('key')} - {dataset['library'].get('label')}")
    add(f"               positions {dataset['library'].get('positions')}, parent "
        f"{dataset['library'].get('parent_residues')}"
        + (f", labelled stage {dataset['library'].get('labelled_stage')!r}"
           if dataset['library'].get('labelled_stage') else "")
        + (", POOL-RESTRICTED" if dataset['library'].get('pool_restricted') else ""))
    add(f"Dataset:       {dataset['path']}")
    add(f"               sha256 {dataset['sha256'][:16]}..., "
        f"{dataset['unique_variants']:,} eligible variants")
    add("")
    add("Configuration identity (manifest.json is authoritative; this is a summary)")
    add("-" * 78)
    add(f"  configuration_fingerprint: {manifest['configuration_fingerprint']}")
    add(f"  acquisition policy:        {configuration['acquisition_policy']['policy_id']}")
    add(f"  prior ablation arm:        {components['prior_ablation_arm']}")
    add(f"  components enabled:        {', '.join(components['enabled']) or '(none)'}")
    add(f"  components DISABLED:       {', '.join(components['disabled']) or '(none)'}")
    add(f"  code commit:               {configuration['code']['commit']} "
        f"(dirty tree: {configuration['code']['dirty_tree']})")
    add(f"  seed block:                {configuration['seeding']['seed_block_sha256'][:16]}... "
        f"over seeds {configuration['seeding']['seed_values']}")
    add(f"  label provenance:          {configuration['label_provenance']['status']} "
        f"(legacy {configuration['label_provenance']['legacy_status']})")
    add(f"  library:                   {configuration['library_selection']['resolved']} "
        f"(requested {configuration['library_selection']['requested']}), selection mode "
        f"{dataset['library']['selection_mode']}")
    if configuration["library_selection"]["dataset_override"]:
        add(f"  DATASET OVERRIDE:          {configuration['library_selection']['dataset_override']}")
    add("")
    add("Protocol")
    add("-" * 78)
    add(f"  {configuration['seeding']['n_seeds']} paired seeds; {budget['initial_size']} shared "
        f"starting variants, then {budget['rounds']} rounds x {budget['batch_size']}.")
    add(f"  Total budget per method per seed: {budget['total_unique_measurements']} unique "
        "measurements, including the start.")
    # The direction is the library's declaration, not a module constant: FITNESS_DIRECTION
    # no longer exists, and it was never the same claim for every landscape.
    add(f"  Fitness direction: {configuration['assay']['fitness_direction']}.")
    add(f"  Top-1% threshold: {top1['threshold_value']:.6g} (99th percentile over the eligible set).")
    add(f"  Ties: every variant at or above the threshold counts as a hit - "
        f"{top1['n_hits']:,} of {dataset['unique_variants']:,} variants qualify.")
    add(f"  Global optimum in the eligible set: {top1['global_best']:.6g}.")
    add("")
    add("Results")
    add("-" * 78)
    def _num(value: Any, places: int = 4) -> str:
        try:
            return f"{float(value):.{places}f}"
        except (TypeError, ValueError):
            return "-"

    add(f"  {'method':<17} {'status':<13} {'done':>6} {'best fit':>10} {'sd':>8} "
        f"{'top1':>7} {'sd':>7} {'secs':>8}")
    for row in summary_rows:
        add(f"  {row['method']:<17} {row['status']:<13} "
            f"{str(row['seeds_completed']) + '/' + str(row['seeds_attempted']):>6} "
            f"{_num(row['mean_best_fitness']):>10} {_num(row['std_best_fitness']):>8} "
            f"{_num(row['mean_top1_hits_including_initial'], 2):>7} "
            f"{_num(row['std_top1_hits_including_initial'], 2):>7} "
            f"{_num(row['mean_runtime_seconds'], 1):>8}")
    add("")
    add("  'best fit' is the best measured fitness found within the budget, averaged over seeds.")
    add("  'top1' counts top-1% variants found, INCLUDING those present in the shared starting sample")
    add("  (the per-round CSV reports the excluding-initial count separately).")
    add("")
    add("Method status")
    add("-" * 78)
    for row in summary_rows:
        add(f"  {row['method']}: {row['status']}")
        if row["status_detail"]:
            for chunk in str(row["status_detail"]).split("; "):
                add(f"      - {chunk}")
    add("")
    failures = [o for o in outcomes if o.status not in ("completed",)]
    if failures:
        add("Failed runs (not excluded from the denominators above)")
        add("-" * 78)
        for outcome in failures:
            add(f"  {outcome.method} seed {outcome.seed}: {outcome.status} - {outcome.failure_reason}")
        add("")
    notes = [(o.method, n) for o in outcomes for n in o.notes]
    if notes:
        add("Run notes")
        add("-" * 78)
        for method, note in notes:
            add(f"  {method}: {note}")
        add("")
    if run["warnings"]:
        add("Warnings recorded in the manifest")
        add("-" * 78)
        for warning in run["warnings"]:
            add(f"  - {warning}")
        add("")
    add("Limitations")
    add("-" * 78)
    for limitation in manifest["limitations"]:
        add(f"  - {limitation}")
    if plot_note:
        add(f"  - {plot_note}")
    add("")
    add("Cost")
    add("-" * 78)
    add("  No dollar figures are reported. Local baselines make no paid calls; the production")
    add("  adapters make no cloud calls in this benchmark; and the Meta API client in")
    add("  controller/orchestrator.py discards the provider's usage block, so token counts are")
    add("  genuinely unknown rather than estimated.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ==============================================================================================
# Integrity checks on the existing repository
# ==============================================================================================
def repository_state(ignore_prefixes: Sequence[str] = ()) -> Dict[str, Any]:
    """
    Records git status and a hash of tracked source/config files so drift is detectable.

    `ignore_prefixes` drops this benchmark's own outputs from the comparison: the run folder is
    created before the first snapshot and filled afterwards, which is expected, not drift.
    """
    state: Dict[str, Any] = {"git_available": False}

    def _is_ignored(line: str) -> bool:
        path = line[3:].strip().strip('"')
        return any(path.startswith(prefix) for prefix in ignore_prefixes)

    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
                              capture_output=True, text=True, timeout=30)
        status = subprocess.run(["git", "status", "--porcelain=v1", "--untracked-files=all"],
                                cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=60)
        if head.returncode == 0:
            state["git_available"] = True
            state["head"] = head.stdout.strip()
            state["status_lines"] = sorted(
                line for line in status.stdout.splitlines() if line.strip() and not _is_ignored(line)
            )
    except Exception as exc:
        state["error"] = str(exc)

    digest = hashlib.sha256()
    counted = 0
    for pattern in ("*.py", "*.txt", "*.json", "*.sql", "*.md"):
        for file_path in sorted(REPO_ROOT.glob(f"**/{pattern}")):
            parts = set(file_path.parts)
            if parts & {".git", "__pycache__", "node_modules", ".pytest_cache", "benchmark_comparison_results"}:
                continue
            if file_path.name == Path(__file__).name:
                continue
            relative = file_path.relative_to(REPO_ROOT).as_posix()
            if any(relative.startswith(prefix) for prefix in ignore_prefixes):
                continue
            try:
                digest.update(file_path.relative_to(REPO_ROOT).as_posix().encode())
                digest.update(file_path.read_bytes())
                counted += 1
            except OSError:
                continue
    state["source_tree_sha256"] = digest.hexdigest()
    state["files_hashed"] = counted
    return state


# ==============================================================================================
# Orchestration
# ==============================================================================================
ALL_METHODS = ("random", "ridge", "random_forest", "alde", "panoptes_default", "panoptes_muse",
               "panoptes_qualify_then_optimize", "panoptes_conditional_probe")
DEFAULT_METHODS = ("random", "ridge", "random_forest", "panoptes_default")


def _load_panoptes_context(
    library: LibrarySpec,
    pdb: bool = True,
) -> Tuple[Any, str, Optional[Path], LibraryDesign]:
    """
    Step 1 of the production loop: cached protein metadata, no network and no writes.

    This is also where the library's declared design meets the real wild-type sequence,
    so it returns the canonical `LibraryDesign`. Building it validates every declared
    position and parent residue against that sequence and refuses on a mismatch, which
    is the check that stops the dataset being aligned against the wrong reference.
    """
    from workers.rag_worker import fetch_protein_metadata

    metadata_json = library.metadata_json
    if not metadata_json.is_file():
        raise BenchmarkError(
            f"Parent metadata cache {metadata_json} is missing for library "
            f"{library.key!r}. Generate it with:\n"
            f"  python scripts/download_trpb_benchmark.py"
        )
    metadata = fetch_protein_metadata(library.accession, use_cache=True)
    sequence = metadata.sequence
    if max(library.positions) > len(sequence):
        raise BenchmarkError(
            f"library {library.key!r} declares position {max(library.positions)} but "
            f"{library.accession}'s sequence is {len(sequence)} residues long."
        )
    observed = "".join(sequence[p - 1] for p in library.positions)
    if observed != library.parent_residues:
        raise BenchmarkError(
            f"Parent residues at {list(library.positions)} are '{observed}', expected "
            f"'{library.parent_residues}'. Residue numbering does not match the dataset; "
            f"refusing to run."
        )
    # Raises WildTypeMismatchError / AlphabetError on any disagreement, rather than
    # quietly tolerating a library description the sequence does not support.
    design = library.to_design(sequence)
    structure = library.structure_pdb
    pdb_path = structure if (pdb and structure is not None and structure.is_file()) else None
    return metadata, sequence, pdb_path, design


def build_selector_factories(
    args: argparse.Namespace,
    methods: Sequence[str],
    library: LibrarySpec,
    qualify_absolute: Optional[float] = None,
) -> Tuple[Dict[str, Callable[[int], Selector]], Dict[str, Dict[str, str]]]:
    """
    Resolves each requested method into a factory, or a visible unavailable/unsupported status.

    `library` is threaded in rather than read from module state, so the factories a run
    builds are tied to the library that run loaded.

    Optional integrations are resolved defensively: a failure here marks that one method and
    never prevents the local baselines from running.
    """
    factories: Dict[str, Callable[[int], Selector]] = {}
    status: Dict[str, Dict[str, str]] = {}

    for method in methods:
        if method == "random":
            factories[method] = lambda seed: RandomSelector(seed)
            status[method] = {"status": "available", "detail": "uniform random, no model"}
        elif method == "ridge":
            factories[method] = lambda seed: RidgeSelector(seed)
            status[method] = {
                "status": "available",
                "detail": (
                    f"sklearn Ridge(alpha=1.0) on "
                    f"{library.n_sites * len(AA_ALPHABET)} one-hot features"
                ),
            }
        elif method == "random_forest":
            factories[method] = lambda seed: RandomForestSelector(seed)
            status[method] = {"status": "available", "detail": "sklearn RandomForestRegressor(n_estimators=300)"}
        elif method == "alde":
            # Resolution and probing are separate: the resolver records every path it
            # looked at so an UNAVAILABLE reason can name them, which is what turns an
            # environment problem into something the next person can act on.
            resolution = resolve_alde_paths(args)
            checkout = resolution["checkout"]
            python_exe = resolution["python_exe"]
            probe = probe_alde(checkout, python_exe, resolution)
            if probe["available"]:
                workspace = (Path(args.alde_workspace).expanduser().resolve()
                             if args.alde_workspace else checkout.parent / "benchmark_runs")
                workspace.mkdir(parents=True, exist_ok=True)
                factories[method] = (
                    lambda seed, c=checkout, pe=python_exe, w=workspace, p=probe:
                    AldeSelector(seed, c, pe, w, p, timeout_seconds=args.alde_timeout,
                                 train_iter=args.alde_train_iter)
                )
                versions = probe.get("dependencies", {})
                status[method] = {
                    "status": "available",
                    "detail": (
                        f"checkout {probe['checkout_path']} @ {(probe['commit'] or 'unknown')[:12]}; "
                        f"onehot + DNN_ENSEMBLE + TS, run unmodified in {python_exe}; "
                        f"torch {versions.get('torch')}, botorch {versions.get('botorch')}, "
                        f"gpytorch {versions.get('gpytorch')}, numpy {versions.get('numpy')}; "
                        f"outputs under {workspace}"
                    ),
                }
            else:
                status[method] = {"status": "UNAVAILABLE", "detail": probe["reason"]}
            status[method]["probe"] = json.dumps(probe)
        elif method in ("panoptes_default", "panoptes_muse",
                        "panoptes_qualify_then_optimize", "panoptes_conditional_probe"):
            live = method == "panoptes_muse"
            # The ABSOLUTE the Evaluator resolved from the eligible universe - not
            # args.qualify_threshold, which is now an optional override and is None on a
            # default run. Reading the raw arg here is what let the bar never transfer.
            qualify = (
                qualify_absolute
                if method in ("panoptes_qualify_then_optimize", "panoptes_conditional_probe")
                else None
            )
            if qualify is None and method in (
                "panoptes_qualify_then_optimize", "panoptes_conditional_probe"
            ):
                raise BenchmarkError(
                    f"{method} needs a resolved qualification bar, but none was supplied. "
                    f"Pass the Evaluator's qualify_threshold_absolute: an arm whose bar is "
                    f"unset silently becomes the single-policy baseline."
                )
            probe = method == "panoptes_conditional_probe"
            if live and not args.enable_live_muse:
                status[method] = {
                    "status": "SKIPPED",
                    "detail": "live Muse comparison requires explicit opt-in: pass --enable-live-muse "
                              "(and set LLAMA_API_KEY or META_API_KEY)",
                }
                continue
            if live and not (os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")):
                status[method] = {
                    "status": "UNAVAILABLE",
                    "detail": "--enable-live-muse was passed but neither LLAMA_API_KEY nor META_API_KEY is "
                              "set, so no live call can be made. Refusing to run the mock planner and label "
                              "the result 'live'.",
                }
                continue
            try:
                metadata, sequence, pdb_path, _design = _load_panoptes_context(library)
            except Exception as exc:
                status[method] = {"status": "UNSUPPORTED", "detail": f"{type(exc).__name__}: {exc}"}
                continue
            # Parsed once, up front: an invalid --ablate-components value must not be
            # discovered inside the first round of the first seed.
            components = resolve_ablate_components(args)
            factories[method] = (
                lambda seed, m=metadata, s=sequence, p=pdb_path, l=live, q=qualify, pr=probe,
                lib=library, c=components:
                PanoptesSelector(seed, s, m, l, p, lib, call_budget=args.max_model_calls,
                                 ablation=args.ablation,
                                 prior_spread=args.prior_spread,
                                 ablate_components=set(c),
                                 qualify_threshold=q, probe_after_qualify=pr)
            )
            detail = (
                "production path: fetch_protein_metadata -> "
                + ("Orchestrator._generate_strategy (live Muse)" if live
                   else "_AnnotationDrivenPlanner.generate_strategy (search space from annotated sites)")
                + " -> select_next_micro_batch"
                + (" -> Orchestrator._generate_reflection" if live else "")
                + "; steps 4 (gatekeeper) and 5 (ESM execution) are not exercised - "
                  "run_stream() offers no interface to substitute the experimental assay for ESM scoring"
                # The arm's identity, so a report row cannot be confused with another arm's.
                # manifest.json remains the authoritative, machine-readable version.
                + f"; prior ablation arm '{args.ablation}'"
                + (f"; components disabled: {', '.join(components)}" if components
                   else "; all selector components enabled")
            )
            status[method] = {"status": "available", "detail": detail}
        else:
            status[method] = {"status": "UNSUPPORTED", "detail": f"unknown method '{method}'"}
    return factories, status


def execute(args: argparse.Namespace, library: LibrarySpec) -> int:
    """
    Runs one benchmark against one library. The library is an argument, not state.

    This frame only creates the run folder and attaches the run log; `_execute_run` does
    the work, so the log handler can be detached in a `finally` rather than leaking for
    the life of the process.
    """
    started = datetime.now(timezone.utc)
    run_id = f"run_{started.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    output_root = Path(args.output_dir).expanduser().resolve()
    run_dir = output_root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)  # never overwrites an existing run
    cache_dir = run_dir / "cache"
    cache_dir.mkdir()

    file_handler = logging.FileHandler(run_dir / "benchmark.log", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logging.getLogger().addHandler(file_handler)
    logger.info("Run %s -> %s", run_id, run_dir)
    try:
        return _execute_run(args, library, run_id=run_id, run_dir=run_dir,
                            cache_dir=cache_dir, output_root=output_root, started=started)
    finally:
        # Detached and closed even when the run aborts: otherwise the log file stays open
        # for the lifetime of the process and every later run in it logs to this folder too.
        logging.getLogger().removeHandler(file_handler)
        file_handler.close()


def _execute_run(
    args: argparse.Namespace,
    library: LibrarySpec,
    *,
    run_id: str,
    run_dir: Path,
    cache_dir: Path,
    output_root: Path,
    started: datetime,
) -> int:
    try:
        ignore_prefixes = (output_root.relative_to(REPO_ROOT).as_posix(),)
    except ValueError:
        ignore_prefixes = ()
    before_state = repository_state(ignore_prefixes)

    # ---- Dataset ---------------------------------------------------------------------------
    self_checks: List[Tuple[str, bool, str]] = []
    # A planted-epistasis run is a synthetic-landscape run by definition: there is no
    # real dataset with an interaction we chose. Treating the flag as implying synthetic
    # avoids the alternative, which is loading the real dataset and silently ignoring
    # the plant.
    use_synthetic = args.smoke_test or args.plant_sign_epistasis

    if use_synthetic:
        self_checks = run_self_checks(library)
        failed = [name for name, passed, _ in self_checks if not passed]
        print("\nSelf-checks")
        print("-" * 78)
        for name, passed, detail in self_checks:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        if failed:
            raise BenchmarkError(f"self-checks failed: {failed}")
        universe = build_synthetic_universe(
            library,
            n_variants=args.synthetic_size,
            seed=args.synthetic_seed,
            plant_sign_epistasis=args.plant_sign_epistasis,
        )
        logger.warning("SMOKE TEST: using a SYNTHETIC landscape. These numbers are not biology.")
    else:
        universe = load_universe(
            library=library,
            dataset_path=Path(args.dataset).expanduser().resolve(),
            provenance_mode=args.provenance_mode,
            imputed_ids_path=Path(args.imputed_ids).expanduser().resolve() if args.imputed_ids else None,
            max_variants=args.max_variants,
            subsample_seed=args.subsample_seed,
            strict_rows=not args.lenient_rows,
        )

    # Carries the explicit vocabulary first and the legacy string after it, so old grep
    # patterns over existing run folders keep working.
    stamp = provenance_stamp(universe)
    logger.info(
        "Eligible universe: %d variants (%s)", len(universe), stamp,
    )

    # The qualification bar is resolved inside the Evaluator, against the same eligible
    # universe that defines the top-1% threshold. argparse has already guaranteed that at
    # most one of the two forms was given.
    evaluator = Evaluator(
        universe,
        qualify_quantile=getattr(args, "qualify_quantile", None),
        qualify_threshold=getattr(args, "qualify_threshold", None),
    )
    logger.info(
        "Qualification bar: %.4f (%s%s); %d of %d eligible variants qualify. Top-1%% "
        "threshold is %.4f.",
        evaluator.qualify_threshold_absolute,
        evaluator.qualify_threshold_source,
        f" q={evaluator.qualify_quantile}" if evaluator.qualify_quantile is not None else "",
        evaluator.n_qualifying,
        len(universe),
        evaluator.top1_threshold,
    )
    if evaluator.n_qualifying <= 1:
        logger.warning(
            "Only %d eligible variant(s) meet the qualification bar %.4f, so qualifying is "
            "equivalent to finding the global best on this landscape. A quantile transfers "
            "the rarity of the bar, not its attainability.",
            evaluator.n_qualifying, evaluator.qualify_threshold_absolute,
        )
    if evaluator.qualify_threshold_absolute < evaluator.top1_threshold:
        logger.warning(
            "The qualification bar %.4f sits BELOW this landscape's own top-1%% threshold "
            "%.4f, so 'qualifying' is an easier bar than being a hit, which inverts the "
            "stage's purpose.",
            evaluator.qualify_threshold_absolute, evaluator.top1_threshold,
        )
    view = CandidateView(
        variant_ids=universe.variant_ids,
        one_hot=universe.one_hot,
        index_of=universe.index_of,
        combos=universe.combos,
        library=universe.library,
        batch_size=args.batch_size,
        n_rounds=args.rounds,
    )

    total_budget = args.initial_size + args.rounds * args.batch_size
    if total_budget > len(universe):
        raise BenchmarkError(
            f"Budget of {total_budget} measurements exceeds the eligible universe of {len(universe)}."
        )

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    unknown = [m for m in methods if m not in ALL_METHODS]
    if unknown:
        raise BenchmarkError(f"Unknown method(s) {unknown}. Choose from: {', '.join(ALL_METHODS)}")
    components_disabled = resolve_ablate_components(args)
    factories, method_status = build_selector_factories(
        args, methods, library, qualify_absolute=evaluator.qualify_threshold_absolute,
    )

    # The starting samples are drawn before anything runs: they are part of the
    # configuration, not a by-product of it, so the manifest can record the seed block
    # even for a run that aborts before its first measurement.
    seeds = list(range(args.base_seed, args.base_seed + args.seeds))
    initial_sets: Dict[int, List[str]] = {
        seed: random.Random(seed).sample(universe.variant_ids, args.initial_size)
        for seed in seeds
    }

    def _manifest(
        outcomes: Sequence[RunOutcome],
        finished: datetime,
        after_state: Optional[Dict[str, Any]] = None,
        run_state: Optional[str] = None,
        abort_reason: str = "",
        extra_warnings: Sequence[str] = (),
        limitations: Sequence[str] = (),
    ) -> Dict[str, Any]:
        return build_manifest(
            run_id=run_id,
            args=args,
            universe=universe,
            evaluator=evaluator,
            methods=methods,
            method_status=method_status,
            method_configs=_method_configs(factories, methods),
            components_disabled=components_disabled,
            seeds=seeds,
            initial_sets=initial_sets,
            outcomes=outcomes,
            started=started,
            finished=finished,
            repository_state_before=before_state,
            repository_state_after=after_state,
            self_checks=self_checks,
            extra_warnings=extra_warnings,
            run_state=run_state,
            abort_reason=abort_reason,
            limitations=limitations,
        )

    # A requested arm that cannot run aborts here - before a single measurement is spent -
    # and the aborted run still leaves a manifest saying what was missing and why.
    try:
        enforce_requested_methods_available(
            methods, method_status,
            allow_unavailable=bool(getattr(args, "allow_unavailable_baselines", False)),
        )
    except BaselineUnavailableError as exc:
        manifest = _manifest(
            outcomes=[],
            finished=datetime.now(timezone.utc),
            run_state="ABORTED",
            abort_reason=str(exc),
            extra_warnings=[f"run ABORTED before any measurement: {exc}"],
            limitations=["ABORTED RUN: no measurements were taken. See run.abort_reason."],
        )
        (run_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8")
        logger.error("Run %s ABORTED: %s", run_id, exc)
        print(f"Aborted-run manifest: {run_dir / 'manifest.json'}", file=sys.stderr)
        raise

    # ---- Run -------------------------------------------------------------------------------
    # Captures each arm's REAL candidate pool, so coverage is not computed from the
    # selection trace alone (which holds committed picks plus harness rejections and is
    # a lower bound on the ballot). Write-only from the selector's side.
    from veritas.diagnostics.benchmark_hook import CandidatePoolRecorder

    pool_recorder = CandidatePoolRecorder()

    outcomes: List[RunOutcome] = []
    for seed in seeds:
        initial_ids = initial_sets[seed]  # shared across methods within a seed
        for method in methods:
            factory = factories.get(method)
            if factory is None:
                continue
            logger.info("seed %d | %s", seed, method)
            try:
                selector = factory(seed)
                pool_recorder.attach(selector, method=method, seed=seed)
                outcome = run_single(selector, evaluator, view, initial_ids, seed)
            except Exception as exc:
                logger.error("%s seed %d could not start: %s", method, seed, exc)
                outcome = RunOutcome(method=method, seed=seed, status="failed",
                                     failure_reason=f"{type(exc).__name__}: {exc}")
            outcomes.append(outcome)
            if outcome.status == "completed" and outcome.rounds:
                final = outcome.rounds[-1]
                logger.info(
                    "  -> best %.6g | top1 %d (incl. start) | %d measurements | %.1fs",
                    final.best_fitness, final.cum_top1_including_initial,
                    final.n_measurements, final.cumulative_seconds,
                )

    for method in methods:
        if method in factories and any(o.method == method and o.status == "failed" for o in outcomes):
            method_status[method]["detail"] += "; one or more seeds failed (see report)"

    # ---- Outputs ---------------------------------------------------------------------------
    summary_rows = summarise(outcomes, method_status)
    write_trace_csv(run_dir / "selection_trace.csv", outcomes, stamp)
    write_round_metrics_csv(run_dir / "round_metrics.csv", outcomes, stamp)
    write_summary_csv(run_dir / "summary.csv", summary_rows, stamp)
    # Candidate-GENERATION diagnostics. One hook: it derives the library from the
    # eligible ids, writes generation_coverage.{csv,json} and depth_returns.csv beside
    # the files above, warns on an arm whose generator collapsed, and never raises.
    # `evaluator.is_top1` is consulted here, after the run, over variants the campaigns
    # already paid to measure; it is never handed to a selector.
    from veritas.diagnostics.benchmark_hook import emit_coverage_artifacts

    emit_coverage_artifacts(run_dir, outcomes, eligible_ids=universe.variant_ids,
                            stamp=stamp, top1_predicate=evaluator.is_top1, run_id=run_id,
                            proposed_pools=pool_recorder)

    plot_note = write_plot(run_dir / "best_fitness_vs_measurements.png", outcomes, stamp,
                           universe.synthetic, cache_dir, library.label)
    if plot_note:
        logger.warning(plot_note)

    after_state = repository_state(ignore_prefixes)
    unchanged = (
        before_state.get("source_tree_sha256") == after_state.get("source_tree_sha256")
        and before_state.get("status_lines") == after_state.get("status_lines")
    )

    limitations = [
        "Methods diverge after the shared starting sample, so per-round comparisons are paired "
        "only at round 0.",
        "panoptes_default / panoptes_muse exercise production planning, selection and (live only) "
        "reflection. Gatekeeper authorisation and Modal ESM-2 execution are not exercised, because "
        "Orchestrator.run_stream() has no interface for substituting an experimental assay for ESM "
        "scoring in step 5.",
        "panoptes_default plans from the target's annotated active sites "
        "(metadata.binding_positions), not from a per-target script. It deliberately does NOT "
        "use workers.mock_workers.MockLLMWorker, whose TrpB branch hardcodes positions "
        "183/184/227/228 and suggests three of the four global-optimum residues - that would "
        "encode the answer this benchmark measures the discovery of.",
        "panoptes_muse is given this library's saturated positions as its assayable scope, "
        "which is enforced deterministically on the planner's output. Without it the planner "
        "proposed distal allosteric residues the four-site library cannot measure and whole "
        "rounds were discarded as not_in_eligible_universe. Check the accepted/requested "
        "columns in round_metrics.csv: a short batch still means proposals were rejected.",
        "Token usage is unknown: controller/orchestrator.py discards the provider's usage block.",
        "Top-1% membership is defined over the eligible set in this run; a different eligibility "
        "filter changes the threshold.",
    ]
    if universe.synthetic:
        limitations.insert(0, "SYNTHETIC RUN: values are generated, not measured. These results say "
                              "nothing about biological performance.")
    if universe.unlabelled_excluded:
        limitations.insert(0, (
            f"{universe.unlabelled_excluded:,} row(s) of {universe.fingerprint.path} carry no "
            f"label (stage != {library.labelled_stage!r}) and are excluded. Their "
            f"DMS_score placeholder of {library.unlabelled_placeholder!r} is an ABSENT "
            f"LABEL, not a measurement; the eligible universe is the labelled split only."
        ))
    if universe.censored_at_floor:
        limitations.insert(0, (
            f"{universe.censored_at_floor:,} of {len(universe):,} eligible measurements "
            f"({100.0 * universe.censored_at_floor / max(len(universe), 1):.1f}%) sit at the "
            f"declared assay floor of {library.detection_floor} and are LEFT-CENSORED, not "
            f"exact values. This library is censored from below: a large mass of ties at the "
            f"floor is expected, and a model treating them as exact continuous labels is "
            f"mis-specified."
        ))
    if library.pool_restricted:
        limitations.insert(0, (
            f"POOL-RESTRICTED run: library {library.key!r} is a sparse sample of its "
            f"combinatorial product, so the selector was handed the enumerated untested "
            f"measured variants and RANKED them. That is a strictly easier problem than "
            f"generating under mutation constraints, and success here does not validate "
            f"unconstrained generation."
        ))
    # The provenance limitations are generated from the resolved vocabulary - including
    # the library-specific imputed-count text, which is None for GB1 and a recorded 0 for
    # HER2 - so the report and the manifest cannot disagree about what the labels are.
    for limitation in reversed(_label_provenance_block(universe)["limitations"]):
        limitations.insert(0, limitation)
    if not unchanged:
        limitations.insert(0, "WARNING: the repository source hash changed during this run - see manifest.")

    extra_warnings: List[str] = []
    if plot_note:
        extra_warnings.append(plot_note)
    if not unchanged:
        extra_warnings.append(
            "the repository source hash or git status changed during this run; see "
            "run.repository_state_before / run.repository_state_after"
        )

    manifest = _manifest(
        outcomes=outcomes,
        finished=datetime.now(timezone.utc),
        after_state=after_state,
        extra_warnings=extra_warnings,
        limitations=limitations,
    )
    # No `default=str` fallback: a manifest that needs one is not machine-readable.
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    write_report(run_dir / "report.txt", manifest, summary_rows, outcomes, stamp, plot_note)

    print()
    print((run_dir / "report.txt").read_text(encoding="utf-8"))
    print(f"Artifacts: {run_dir}")
    print(f"Run state: {manifest['run']['state']} | "
          f"configuration fingerprint: {manifest['configuration_fingerprint'][:16]}")
    if not unchanged:
        print("WARNING: repository files changed during this run. Inspect manifest.json "
              "(run.repository_state_before / run.repository_state_after).", file=sys.stderr)
        return EXIT_REPOSITORY_CHANGED
    return EXIT_OK


def _method_configs(
    factories: Dict[str, Callable[[int], Selector]],
    methods: Sequence[str],
) -> Dict[str, Dict[str, Any]]:
    """
    Collects each available method's own hyperparameter block for the manifest.

    A construction failure here must not cost the manifest: it is recorded as the config
    for that method instead of propagating.
    """
    configs: Dict[str, Dict[str, Any]] = {}
    for method in methods:
        factory = factories.get(method)
        if factory is None:
            configs[method] = {}
            continue
        try:
            configs[method] = dict(factory(0).config)
        except Exception as exc:
            configs[method] = {"config_unavailable": f"{type(exc).__name__}: {exc}"}
    return configs


def _package_versions() -> Dict[str, Optional[str]]:
    versions: Dict[str, Optional[str]] = {}
    for name in ("sklearn", "scipy", "matplotlib", "torch", "gpytorch", "botorch", "pydantic"):
        try:
            versions[name] = getattr(__import__(name), "__version__", "unknown")
        except Exception:
            versions[name] = None
    return versions


# ==============================================================================================
# CLI
# ==============================================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="benchmark_panoptes_isolated.py",
        description=(
            "Isolated comparison benchmark for PanoptesBio on a registered combinatorial "
            "library (--library). Writes only into a fresh folder under "
            "benchmark_comparison_results/ and never modifies the existing project."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # mechanism validation on synthetic data, a few seconds, no project data needed
  python benchmark_panoptes_isolated.py --smoke-test

  # real data, local baselines plus the production path (no paid calls)
  python benchmark_panoptes_isolated.py --imputed-ids data/trpb_imputed.txt

  # optional: the published ALDE implementation from a local checkout
  python benchmark_panoptes_isolated.py --methods alde --alde-path ~/src/ALDE

  # optional: live Muse planning (makes paid Meta API calls; capped)
  python benchmark_panoptes_isolated.py --methods panoptes_muse --enable-live-muse --max-model-calls 20
""",
    )
    parser.add_argument("--smoke-test", action="store_true",
                        help="Run a small bounded CPU test on a clearly-labelled SYNTHETIC landscape.")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET),
                        help="Path to the prepared CSV for the chosen --library. Left at this "
                             "default while --library names another library, it is replaced by "
                             "that library's own dataset (default: %(default)s).")
    parser.add_argument("--methods", default=",".join(DEFAULT_METHODS),
                        help=f"Comma-separated subset of: {', '.join(ALL_METHODS)} (default: %(default)s).")
    parser.add_argument("--seeds", type=int, default=20, help="Number of paired seeds (default: %(default)s).")
    parser.add_argument("--base-seed", type=int, default=0, help="First seed value (default: %(default)s).")
    parser.add_argument("--initial-size", type=int, default=25,
                        help="Shared starting variants per seed (default: %(default)s).")
    parser.add_argument("--rounds", type=int, default=4,
                        help="Selection rounds after the start (default: %(default)s).")
    parser.add_argument("--batch-size", type=int, default=25,
                        help="New variants per round (default: %(default)s).")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Parent folder for run folders (default: %(default)s).")
    parser.add_argument(
        "--ablation",
        choices=("none", "aux-weighted", "raw-prior", "legacy"),
        default="none",
        help=(
            "Restore a pre-fix behaviour for a paired comparison. none = shipped "
            "configuration; aux-weighted = re-apply the synthesizability and "
            "developability multipliers to acquisition value; raw-prior = bypass the "
            "prior fitness-scale calibration; legacy = both. The arm name in every "
            "report carries the setting."
        ),
    )
    parser.add_argument(
        "--prior-spread", type=float, default=None,
        help=(
            "Calibrate the biophysical prior to this magnitude in fitness units and "
            "measure the result. Omitted leaves the shipped default (calibration off, "
            "raw magnitudes), which is the best configuration measured so far. 0.0 "
            "ablates the prior entirely."
        ),
    )
    parser.add_argument(
        "--ablate-components", default="",
        help=(
            "Comma-separated selector components to disable, for attributing the "
            f"engine advantage: {', '.join(ABLATABLE_COMPONENTS)}. Omitting the flag "
            "disables none. A value that names no component (',' or whitespace), an "
            "empty entry ('dpp,,prior') or an unrecognised name is rejected before the "
            "run starts, because each of those used to produce an unablated run in a "
            "folder labelled as an ablation arm."
        ),
    )
    parser.add_argument(
        "--allow-unavailable-baselines", action="store_true",
        help=(
            "Continue when a method named in --methods cannot run, reporting it as a row "
            "instead of failing. OFF by default: a requested ALDE baseline was once "
            "silently marked UNAVAILABLE because ALDE_PATH was unset, and two complete "
            "benchmark runs were quoted from that report. Use this only to probe what an "
            "environment can run - a run made with it is not a complete comparison."
        ),
    )
    parser.add_argument("--provenance-mode", choices=("require", "acknowledge-unflagged"), default="require",
                        help="'require' (default) refuses to run unless measured records can be "
                             "distinguished from the published KNN-imputed ones.")
    parser.add_argument("--imputed-ids", default=None,
                        help="File listing the imputed (non-measured) variants to exclude.")
    parser.add_argument("--lenient-rows", action="store_true",
                        help="Drop and COUNT dataset rows that do not map onto the library's "
                             "declared positions, instead of refusing to run. Off by default: a "
                             "dataset whose rows disagree with the declared design is a "
                             "configuration error, and dropping them silently shrinks the "
                             "eligible universe without changing any printed number.")
    parser.add_argument("--max-variants", type=int, default=None,
                        help="Optionally subsample the candidate universe (seeded) for a faster run.")
    parser.add_argument("--subsample-seed", type=int, default=0, help="Seed for --max-variants.")
    parser.add_argument("--synthetic-size", type=int, default=0,
                        help="Variant count for the synthetic smoke-test landscape; 0 (the default) "
                             "covers the complete 20^4-1 combo space so the production adapter's "
                             "selections are in-universe.")
    parser.add_argument("--alde-path", default=None,
                        help="Path to an existing clone of https://github.com/jsunn-y/ALDE.")
    parser.add_argument("--alde-python", default=None,
                        help="Python executable of ALDE's own isolated environment. Defaults to "
                             "<alde-path>/../.venv-alde/Scripts/python.exe (or bin/python).")
    parser.add_argument("--alde-workspace", default=None,
                        help="Folder for ALDE's generated data, tensors and logs. Defaults to "
                             "<alde-path>/../benchmark_runs. Never inside this repository.")
    parser.add_argument("--alde-timeout", type=int, default=7200,
                        help="Per-seed timeout in seconds for the ALDE subprocess (default: %(default)s).")
    parser.add_argument("--alde-train-iter", type=int, default=300,
                        help="ALDE train_iter; 300 is the authors' default (default: %(default)s).")
    parser.add_argument("--enable-live-muse", action="store_true",
                        help="Explicitly permit live Meta/Muse API calls. Off by default; without it the "
                             "live comparison is reported as SKIPPED rather than silently mocked.")
    parser.add_argument("--max-model-calls", type=int, default=50,
                        help="Hard cap on live model calls across the whole run (default: %(default)s).")
    parser.add_argument("--library", choices=tuple(sorted(LIBRARIES)),
                        default=DEFAULT_LIBRARY_KEY,
                        help="Which registered library to run. 'trpb' is the development "
                             "landscape and has been extensively analysed, so results on it "
                             "are development evidence; 'gb1' is an independent saturated "
                             "library and is where a generalisation claim has to be earned; "
                             "'her2' is a sparse ten-site library run POOL-RESTRICTED "
                             "(default: %(default)s).")
    parser.add_argument("--synthetic-seed", type=int, default=7,
                        help="RNG seed for the synthetic landscape. Varying it gives "
                             "independent landscapes, which is how a policy is tested "
                             "without being tuned to one interaction (default: %(default)s).")
    parser.add_argument("--plant-sign-epistasis", action="store_true",
                        help="Plant a sign-epistatic residue pair in the synthetic "
                             "landscape: each residue alone is deleterious, together "
                             "strongly beneficial. Reproduces the structure that defeats "
                             "marginal-fitted surrogates on TrpB, at a pair chosen by the "
                             "generator's RNG rather than one we already know. Requires "
                             "--synthetic-size.")
    # The qualification bar, stated EITHER as a quantile of the eligible universe (the
    # default, because it is the only form that transfers between landscapes) or as an
    # absolute fitness (an explicit override). Mutually exclusive: a run gated by two
    # different bars cannot say which one it ran under.
    qualify = parser.add_mutually_exclusive_group()
    qualify.add_argument("--qualify-quantile", type=float, default=None,
                         help=f"Quantile (0-100) of the eligible universe at which "
                              f"panoptes_qualify_then_optimize switches from the baseline "
                              f"policy to frontier search (default: "
                              f"{DEFAULT_QUALIFY_QUANTILE}). Resolved to an absolute fitness "
                              f"inside the Evaluator, from the same values array that defines "
                              f"the top-1% threshold, so the two bars cannot drift apart. The "
                              f"default reproduces the old absolute 0.70 on TrpB exactly "
                              f"(0.6998, the identical 64 variants). Mutually exclusive with "
                              f"--qualify-threshold.")
    qualify.add_argument("--qualify-threshold", type=float, default=None,
                         help=f"Absolute fitness at which panoptes_qualify_then_optimize "
                              f"switches to frontier search - an EXPLICIT OVERRIDE of "
                              f"--qualify-quantile. There is no default: an absolute is not "
                              f"comparable across landscapes (the pre-fix default "
                              f"{LEGACY_QUALIFY_ABSOLUTE} admitted 100% of three of them), so "
                              f"it has to be asked for. Under deterministic reveal, P(reach "
                              f"threshold) and time-to-threshold are preserved pathwise, "
                              f"because every trajectory is identical to the baseline until "
                              f"the protected event has already occurred.")
    parser.add_argument("--verbose", action="store_true", help="Debug-level console logging.")
    return parser


def _load_repo_dotenv() -> None:
    """
    Loads REPO_ROOT/.env before any credential is read.

    `controller.orchestrator` calls load_dotenv() at import time, but it is imported
    lazily inside PanoptesSelector._planner(), which runs long after the live-Muse
    availability gate has already consulted os.environ. Without this, a key present
    in .env makes the gate report panoptes_muse as UNAVAILABLE and the live
    comparison is silently skipped. Deliberate and explicit rather than relying on
    another module's import side effect.
    """
    env_path = REPO_ROOT / ".env"
    if env_path.is_file():
        try:
            from dotenv import load_dotenv

            load_dotenv(env_path, override=False)
        except ImportError:
            logger.debug("python-dotenv unavailable; relying on the ambient environment")

    # After .env is loaded, so VERITAS_USE_OS_TRUSTSTORE can be set there, and before
    # any live planner or Modal call. On an intercepted network, omitting this makes
    # every live arm fail over to a mock and report as though it had run live.
    try:
        from controller.tls_trust import enable_os_trust_store

        if enable_os_trust_store():
            logger.info("TLS verification is using the OS certificate store")
    except ImportError:
        pass


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    _load_repo_dotenv()

    logging.basicConfig(  # configured before resolve_library so its log line is visible
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("workers").setLevel(logging.WARNING)
    logging.getLogger("controller").setLevel(logging.WARNING)

    # The library is resolved once and then passed down every call path. Nothing is
    # reassigned at module scope, so this is ordering for the log line only.
    spec = resolve_library(args.library)
    if args.dataset == str(LIBRARIES[DEFAULT_LIBRARY_KEY].dataset) and args.library != DEFAULT_LIBRARY_KEY:
        # --dataset still carries the argparse default, so the user selected a library
        # without naming a file. Point it at that library's dataset instead of silently
        # scoring one library's variants against another's table.
        args.dataset = str(spec.dataset)

    if args.smoke_test:
        # A bounded CPU run: small landscape, few seeds, short rounds, no paid calls.
        args.seeds = min(args.seeds, 2)
        args.initial_size = min(args.initial_size, 10)
        args.rounds = min(args.rounds, 2)
        args.batch_size = min(args.batch_size, 10)
        args.enable_live_muse = False

    for name, value in (("--seeds", args.seeds), ("--rounds", args.rounds),
                        ("--batch-size", args.batch_size), ("--initial-size", args.initial_size)):
        if value < 1:
            parser.error(f"{name} must be >= 1")

    try:
        # Validated here, before any run folder exists: a malformed ablation request must
        # not get as far as creating an output folder, let alone a round.
        resolve_ablate_components(args)
    except ConfigurationError as exc:
        print(f"\nCONFIGURATION ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return EXIT_ERROR

    try:
        return execute(args, spec)
    except BaselineUnavailableError as exc:
        print(f"\nREQUESTED BASELINE UNAVAILABLE\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return EXIT_BASELINE_UNAVAILABLE
    except ProvenanceError as exc:
        print(f"\nPROVENANCE ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return EXIT_PROVENANCE
    except ConfigurationError as exc:
        print(f"\nCONFIGURATION ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return EXIT_ERROR
    except BenchmarkError as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    raise SystemExit(main())


