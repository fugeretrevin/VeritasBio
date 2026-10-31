"""
veritas.datasets - the library/dataset concern, extracted from the benchmark.

A new combinatorial protein library becomes benchmarkable by describing it
(`LibraryDesign`) and pointing an `IngestConfig` at its prepared table. No optimisation
code changes, and nothing in here imports the optimiser.

    from veritas.datasets import (
        AssayMeta, FitnessDirection, IngestConfig, LibraryDesign, SelectionMode,
        VariablePosition, ingest_dataset,
    )

    design = LibraryDesign(
        target_id="SPG1_STRSG",
        assay_id="gb1_four_site_wu_2016",
        wild_type_sequence=wt,
        positions=tuple(
            VariablePosition(position=p, wild_type_aa=wt[p - 1], allowed_aa=tuple(AA_ALPHABET))
            for p in (39, 40, 41, 54)
        ),
        fitness_direction=FitnessDirection.HIGHER_IS_BETTER,
        selection_mode=SelectionMode.POOL_RESTRICTED,
    )
    data = ingest_dataset(IngestConfig(design=design, source=csv_path))

Design notes that are load-bearing rather than stylistic:

  * No module-level "active library". Every encoder takes its design explicitly, so two
    libraries are usable in one process - the prerequisite for a cross-protein study.
  * A mutation string OMITS wild-type positions. Everything canonicalises to full length
    over the design's positions before comparing genotypes; unequal lengths raise.
  * Absent label, measured zero and detection floor are three distinct states.
"""

from .errors import (
    AlphabetError,
    ConflictingMeasurementError,
    DatasetError,
    DuplicateGenotypeError,
    GenotypeLengthError,
    IndelError,
    LabelPolicyError,
    PositionOutsideLibraryError,
    SchemaError,
    WildTypeMismatchError,
)
from .ingest import (
    ColumnMap,
    IngestConfig,
    IngestedDataset,
    IngestReport,
    VariantRecord,
    ingest_dataset,
)
from .schema import (
    AA_ALPHABET,
    AssayMeta,
    Censoring,
    FitnessDirection,
    LabelStatus,
    LibraryDesign,
    SelectionMode,
    VariablePosition,
)

__all__ = [
    "AA_ALPHABET",
    "AlphabetError",
    "AssayMeta",
    "Censoring",
    "ColumnMap",
    "ConflictingMeasurementError",
    "DatasetError",
    "DuplicateGenotypeError",
    "FitnessDirection",
    "GenotypeLengthError",
    "IndelError",
    "IngestConfig",
    "IngestReport",
    "IngestedDataset",
    "LabelPolicyError",
    "LabelStatus",
    "LibraryDesign",
    "PositionOutsideLibraryError",
    "SchemaError",
    "SelectionMode",
    "VariablePosition",
    "VariantRecord",
    "WildTypeMismatchError",
    "ingest_dataset",
]
