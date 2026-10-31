"""
Ingestion: a prepared table plus a `LibraryDesign` becomes validated records.

What this guarantees, and why each one is here:

  * Every genotype is FULL LENGTH over the library's positions. Mutation strings omit
    wild-type sites, so `V183A:F184A:V227A` is canonicalised to `AAAS` before anything
    compares it to another genotype.
  * A missing label, a true measured zero and a detection floor are three different
    states. `fitness` is `None` when no label exists, so an absent label cannot be read
    as a measurement of zero; a value at the declared floor is `LEFT_CENSORED` and keeps
    its censoring flag; a measured 0.0 with no floor declared stays `EXACT`.
  * Duplicate genotypes are refused by default. Aggregating them requires naming a
    deterministic policy AND a tolerance, and replicates that disagree beyond that
    tolerance are refused rather than averaged into a plausible number.
  * Validation rejects. A wrong wild-type reference, a residue outside the declared
    alphabet, an indel, an inconsistent genotype length, a position the design does not
    declare, a non-numeric score and a measurement below the declared detection limit
    are all errors that name the row.
  * Nothing infers the permissible alphabet from the observed sequences. An assay that
    covers 11,277 of 20^10 genotypes would otherwise "declare" a library of whatever it
    happened to sample.
"""

from __future__ import annotations

import csv
import hashlib
import math
import random
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .errors import (
    ConflictingMeasurementError,
    DuplicateGenotypeError,
    LabelPolicyError,
    SchemaError,
)
from .schema import (
    Censoring,
    FitnessDirection,
    LabelStatus,
    LibraryDesign,
    SelectionMode,
)

#: Value transformations that may be applied at ingestion. A transformation is RECORDED
#: on the design's assay metadata, so a downstream reader can tell a log-scale table from
#: a linear one instead of inferring it from the value range.
_TRANSFORMS: Dict[str, Any] = {
    "none": lambda v: v,
    "identity": lambda v: v,
    "log10": math.log10,
    "log2": math.log2,
    "ln": math.log,
    "negate": lambda v: -v,
}
#: Transformations that reverse which end of the scale is better.
_DIRECTION_FLIPPING = ("negate",)
_POSITIVE_ONLY = ("log10", "log2", "ln")

_DUPLICATE_POLICIES = ("reject", "first", "mean", "median")
_FLOAT_TOLERANCE = 1e-9


@dataclass(frozen=True)
class ColumnMap:
    """Which columns of the prepared table hold what. Renaming is configuration."""

    variant: str = "mutant"
    fitness: str = "DMS_score"
    #: Optional. When the column is present it is cross-checked against the mutation
    #: string, which is how an indel in the sequence is caught even if the notation hides it.
    sequence: Optional[str] = "mutated_sequence"
    #: Optional. Names the split column used by `labelled_stage`.
    stage: Optional[str] = "stage"
    #: Optional. A column flagging model-imputed rather than measured values.
    provenance: Optional[str] = None
    #: Optional. A column holding the replicate count behind each value.
    n_replicates: Optional[str] = None


@dataclass(frozen=True)
class IngestConfig:
    """Everything that is configuration rather than code for one dataset."""

    design: LibraryDesign
    source: Path
    column_map: ColumnMap = field(default_factory=ColumnMap)
    #: When set, ONLY rows carrying this value in the stage column hold labels. Every
    #: other row is recorded as an absent label and never as a measurement.
    labelled_stage: Optional[str] = None
    #: The numeric placeholder the unlabelled rows carry. When given, every unlabelled
    #: row must carry exactly it; a different value means a real measurement is about to
    #: be discarded, which is refused.
    unlabelled_placeholder: Optional[float] = None
    duplicate_policy: str = "reject"
    #: The spread above which repeated measurements of one genotype are a conflict rather
    #: than replicates. Only consulted by an aggregating policy.
    conflict_tolerance: float = 0.0
    transform: str = "none"
    #: Cross-check the sequence column against the mutation string when both are present.
    verify_sequence_column: bool = True
    max_records: Optional[int] = None
    subsample_seed: int = 0

    def __post_init__(self) -> None:
        if self.duplicate_policy not in _DUPLICATE_POLICIES:
            raise SchemaError(
                f"duplicate_policy {self.duplicate_policy!r} is not one of "
                f"{list(_DUPLICATE_POLICIES)}. 'reject' refuses any duplicate genotype; the "
                f"others aggregate deterministically and are recorded in the report."
            )
        if self.transform not in _TRANSFORMS:
            raise SchemaError(
                f"transform {self.transform!r} is not one of {sorted(_TRANSFORMS)}"
            )
        if self.conflict_tolerance < 0:
            raise SchemaError("conflict_tolerance must be >= 0")
        if self.max_records is not None and self.max_records < 1:
            raise SchemaError("max_records must be >= 1 when given")


@dataclass(frozen=True)
class VariantRecord:
    """One genotype of the library, with its label state made explicit."""

    variant_id: str                 # canonical notation, wild-type sites omitted
    combo: str                      # FULL length over design.positions
    #: None when no label exists. There is no sentinel number for "unmeasured".
    fitness: Optional[float]
    label_status: LabelStatus
    censoring: Censoring
    #: The value before `IngestConfig.transform`, kept so a transformation is auditable.
    raw_fitness: Optional[float] = None
    n_replicates: Optional[int] = None
    replicate_values: Tuple[float, ...] = ()
    stage: Optional[str] = None
    is_imputed: bool = False

    @property
    def is_measured(self) -> bool:
        return self.label_status is LabelStatus.MEASURED

    @property
    def is_censored(self) -> bool:
        return self.censoring in (Censoring.LEFT_CENSORED, Censoring.RIGHT_CENSORED)


@dataclass(frozen=True)
class IngestReport:
    """What ingestion did, in numbers a manifest can carry."""

    rows_read: int
    rows_labelled: int
    rows_absent_label: int
    wild_type_rows: int
    duplicate_groups: int
    duplicate_rows_collapsed: int
    duplicate_aggregation: str
    conflict_tolerance: float
    exact: int
    left_censored: int
    right_censored: int
    imputed_rows: int
    transform: str
    selection_mode: str
    fitness_direction: str
    combinatorial_size: int
    combinations_present: int
    combinations_missing: int
    coverage_fraction: float
    subsampled_to: Optional[int] = None

    @property
    def library_is_complete(self) -> bool:
        """True when every non-wild-type genotype the design permits carries a label."""
        return self.combinations_missing == 0

    def to_dict(self) -> Dict[str, Any]:
        payload = {
            "rows_read": self.rows_read,
            "rows_labelled": self.rows_labelled,
            "rows_absent_label": self.rows_absent_label,
            "wild_type_rows": self.wild_type_rows,
            "duplicate_groups": self.duplicate_groups,
            "duplicate_rows_collapsed": self.duplicate_rows_collapsed,
            "duplicate_aggregation": self.duplicate_aggregation,
            "conflict_tolerance": self.conflict_tolerance,
            "exact": self.exact,
            "left_censored": self.left_censored,
            "right_censored": self.right_censored,
            "imputed_rows": self.imputed_rows,
            "transform": self.transform,
            "selection_mode": self.selection_mode,
            "fitness_direction": self.fitness_direction,
            "combinatorial_size": self.combinatorial_size,
            "combinations_present": self.combinations_present,
            "combinations_missing": self.combinations_missing,
            "coverage_fraction": self.coverage_fraction,
            "library_is_complete": self.library_is_complete,
        }
        if self.subsampled_to is not None:
            payload["subsampled_to"] = self.subsampled_to
        return payload


@dataclass(frozen=True)
class IngestedDataset:
    """The validated result: labelled records, absent labels, and what happened."""

    design: LibraryDesign
    records: Tuple[VariantRecord, ...]
    #: Rows that exist in the file but carry NO measurement. Kept separate from
    #: `records` so no aggregate can accidentally include them.
    absent_labels: Tuple[VariantRecord, ...]
    report: IngestReport
    source_path: str
    source_sha256: str
    source_size_bytes: int

    @property
    def selection_mode(self) -> SelectionMode:
        return self.design.selection_mode

    def __len__(self) -> int:
        return len(self.records)

    def variant_ids(self) -> Tuple[str, ...]:
        return tuple(r.variant_id for r in self.records)

    def combos(self) -> Dict[str, str]:
        return {r.variant_id: r.combo for r in self.records}

    def fitness_table(self) -> Dict[str, float]:
        """Measured values only. An absent label is absent from this mapping."""
        return {r.variant_id: float(r.fitness) for r in self.records if r.fitness is not None}

    def censored_ids(self) -> Tuple[str, ...]:
        return tuple(r.variant_id for r in self.records if r.is_censored)

    def censoring_table(self) -> Dict[str, str]:
        """Censoring per measured genotype, so a model can honour it rather than guess."""
        return {r.variant_id: r.censoring.value for r in self.records}

    def one_hot(self) -> np.ndarray:
        return self.design.one_hot([r.combo for r in self.records])

    def candidate_pool(self, exclude: Sequence[str] = ()) -> Tuple[str, ...]:
        """
        The enumerated feasible library, for POOL_RESTRICTED selection.

        Refused in GENERATIVE mode on purpose: handing a generative policy the measured
        set would turn generation into ranking and quietly convert the harder problem
        into the easier one.
        """
        if not self.design.is_pool_restricted:
            raise SchemaError(
                f"library {self.design.assay_id} declares selection_mode="
                f"{self.design.selection_mode.value}; an enumerated candidate pool is only "
                f"defined for pool_restricted mode. Use mutation_constraints() to generate "
                f"under the library's constraints instead - success at ranking an enumerated "
                f"pool does not validate unconstrained generation."
            )
        skip = set(exclude)
        return tuple(r.variant_id for r in self.records if r.variant_id not in skip)

    def mutation_constraints(self) -> Dict[str, Any]:
        return self.design.mutation_constraints()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "design": self.design.to_dict(),
            "source": {
                "path": self.source_path,
                "sha256": self.source_sha256,
                "size_bytes": self.source_size_bytes,
            },
            "ingestion": self.report.to_dict(),
        }


def _sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _parse_float(text: str, where: str) -> float:
    try:
        value = float(text)
    except (TypeError, ValueError):
        raise SchemaError(
            f"{where}: {text!r} is not a number. A non-numeric score is refused rather "
            f"than dropped, so a column or delimiter mistake cannot shrink the dataset "
            f"quietly."
        ) from None
    if not math.isfinite(value):
        raise SchemaError(f"{where}: {text!r} is not finite")
    return value


def _is_truthy(text: Optional[str]) -> bool:
    return (text or "").strip().lower() in ("1", "true", "yes", "y", "imputed")


def _classify(value: float, design: LibraryDesign, where: str) -> Censoring:
    floor = design.assay.detection_floor
    ceiling = design.assay.detection_ceiling
    if floor is not None:
        if value < floor - _FLOAT_TOLERANCE:
            raise SchemaError(
                f"{where}: value {value!r} is below the declared detection floor {floor!r}. "
                f"Either the floor is wrong or the value is, and clamping it would invent a "
                f"measurement; refusing the row."
            )
        if value <= floor + _FLOAT_TOLERANCE:
            return Censoring.LEFT_CENSORED
    if ceiling is not None and value >= ceiling - _FLOAT_TOLERANCE:
        return Censoring.RIGHT_CENSORED
    return Censoring.EXACT


def _apply_transform(value: float, name: str, where: str) -> float:
    if name in _POSITIVE_ONLY and value <= 0:
        raise SchemaError(
            f"{where}: {name} is undefined for {value!r}. A transformation that cannot be "
            f"applied is an error, not a NaN to be carried forward."
        )
    result = _TRANSFORMS[name](value)
    if not math.isfinite(result):
        raise SchemaError(f"{where}: {name}({value!r}) is not finite")
    return float(result)


@dataclass
class _RawRow:
    variant_id: str
    combo: str
    value: Optional[float]
    stage: Optional[str]
    is_imputed: bool
    line: int
    n_replicates: Optional[int]


def ingest_dataset(config: IngestConfig) -> IngestedDataset:
    """
    Reads `config.source` into validated `VariantRecord`s under `config.design`.

    Raises a subclass of `veritas.datasets.errors.DatasetError` naming the offending row
    rather than returning a partially repaired dataset.
    """
    design = config.design
    path = Path(config.source)
    if not path.is_file():
        raise SchemaError(f"dataset not found: {path}")
    sha, size = _sha256_file(path)
    columns = config.column_map

    raw_rows: List[_RawRow] = []
    absent: List[VariantRecord] = []
    wild_type_rows = 0
    rows_read = 0
    imputed_rows = 0

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]
        for required in (columns.variant, columns.fitness):
            if required not in fieldnames:
                raise SchemaError(
                    f"{path.name} has no {required!r} column (found {fieldnames}). Point "
                    f"IngestConfig.column_map at the right columns rather than renaming the data."
                )
        stage_column = columns.stage if (columns.stage and columns.stage in fieldnames) else None
        if config.labelled_stage is not None and stage_column is None:
            raise SchemaError(
                f"library {design.assay_id} declares labelled_stage="
                f"{config.labelled_stage!r}, but {path.name} has no "
                f"{columns.stage!r} column (found {fieldnames}). Refusing to admit rows as "
                f"measurements without knowing which split they came from."
            )
        sequence_column = (
            columns.sequence if (columns.sequence and columns.sequence in fieldnames) else None
        )
        provenance_column = (
            columns.provenance if (columns.provenance and columns.provenance in fieldnames) else None
        )
        replicate_column = (
            columns.n_replicates
            if (columns.n_replicates and columns.n_replicates in fieldnames)
            else None
        )

        for offset, row in enumerate(reader):
            line = offset + 2  # 1-indexed, header included
            rows_read += 1
            where = f"{path.name} line {line}"
            notation = (row.get(columns.variant) or "").strip()
            combo = design.combo_from_mutations(notation)
            variant_id = design.canonical_variant_id(combo)

            if sequence_column is not None and config.verify_sequence_column:
                sequence = (row.get(sequence_column) or "").strip()
                if sequence:
                    from_sequence = design.combo_from_sequence(sequence)
                    if from_sequence != combo:
                        raise SchemaError(
                            f"{where}: the mutation string {notation!r} gives genotype {combo!r} "
                            f"but the {sequence_column!r} column gives {from_sequence!r}. The two "
                            f"descriptions of one row disagree; refusing to pick one."
                        )

            stage = (row.get(stage_column) or "").strip() if stage_column else None
            is_imputed = _is_truthy(row.get(provenance_column)) if provenance_column else False
            if is_imputed:
                imputed_rows += 1
            n_replicates: Optional[int] = None
            if replicate_column is not None:
                text = (row.get(replicate_column) or "").strip()
                if text:
                    n_replicates = int(_parse_float(text, f"{where} ({replicate_column})"))

            unlabelled = config.labelled_stage is not None and stage != config.labelled_stage
            raw_text = (row.get(columns.fitness) or "").strip()

            if unlabelled:
                # The row exists but holds no measurement. When a placeholder is declared,
                # every unlabelled row must carry exactly it - otherwise a real
                # measurement is about to be thrown away under the label policy.
                if config.unlabelled_placeholder is not None:
                    observed = _parse_float(raw_text, f"{where} ({columns.fitness})")
                    if abs(observed - config.unlabelled_placeholder) > _FLOAT_TOLERANCE:
                        raise LabelPolicyError(
                            f"{where}: row is outside the labelled stage "
                            f"{config.labelled_stage!r} but carries "
                            f"{columns.fitness}={observed!r}, not the declared placeholder "
                            f"{config.unlabelled_placeholder!r}. That looks like a real "
                            f"measurement in the held-out split; refusing to discard it as an "
                            f"absent label."
                        )
                if not variant_id:
                    wild_type_rows += 1
                    continue
                absent.append(VariantRecord(
                    variant_id=variant_id,
                    combo=combo,
                    fitness=None,
                    label_status=LabelStatus.ABSENT,
                    censoring=Censoring.NOT_APPLICABLE,
                    raw_fitness=None,
                    stage=stage,
                    is_imputed=is_imputed,
                ))
                continue

            if not variant_id:
                # The unmutated reference. It carries no mutations, so it is not a member
                # of the library; its value is a useful anchor but not a candidate.
                wild_type_rows += 1
                continue

            value = _parse_float(raw_text, f"{where} ({columns.fitness})")
            _classify(value, design, where)  # rejects a value below the declared floor
            raw_rows.append(_RawRow(
                variant_id=variant_id, combo=combo, value=value, stage=stage,
                is_imputed=is_imputed, line=line, n_replicates=n_replicates,
            ))

    if not raw_rows:
        raise SchemaError(
            f"{path.name}: no labelled rows survived ingestion for library {design.assay_id}. "
            f"Check the column map and the labelled_stage policy."
        )

    # ---- duplicates and replicate aggregation --------------------------------------
    grouped: Dict[str, List[_RawRow]] = {}
    for item in raw_rows:
        grouped.setdefault(item.variant_id, []).append(item)

    duplicate_groups = sum(1 for rows in grouped.values() if len(rows) > 1)
    duplicate_rows_collapsed = sum(len(rows) - 1 for rows in grouped.values() if len(rows) > 1)
    if duplicate_groups and config.duplicate_policy == "reject":
        offenders = sorted(v for v, rows in grouped.items() if len(rows) > 1)[:5]
        raise DuplicateGenotypeError(
            f"{path.name}: {duplicate_groups} genotype(s) appear more than once, e.g. "
            f"{offenders}. Two spellings of one genotype are the same genotype once "
            f"wild-type sites are filled in. Choose a duplicate_policy "
            f"('first', 'mean' or 'median') and a conflict_tolerance to aggregate them "
            f"deterministically; they are not collapsed silently."
        )

    records: List[VariantRecord] = []
    exact = left_censored = right_censored = 0
    for variant_id in sorted(grouped):
        rows = grouped[variant_id]
        values = tuple(float(r.value) for r in rows)  # file order, deterministic
        if len(values) > 1:
            spread = max(values) - min(values)
            if spread > config.conflict_tolerance + _FLOAT_TOLERANCE:
                raise ConflictingMeasurementError(
                    f"{path.name}: genotype {variant_id!r} is measured {len(values)} times with "
                    f"values {list(values)} (spread {spread:.6g} > conflict_tolerance "
                    f"{config.conflict_tolerance:.6g}). Averaging measurements that disagree this "
                    f"much would invent a value; raise conflict_tolerance deliberately or fix the "
                    f"data. Rows: {[r.line for r in rows]}."
                )
        if config.duplicate_policy == "median":
            aggregated = float(statistics.median(values))
        elif config.duplicate_policy == "mean":
            aggregated = float(statistics.fmean(values))
        else:  # 'first', and the single-row case under every policy
            aggregated = values[0]

        censoring = _classify(aggregated, design, f"{path.name} genotype {variant_id}")
        if censoring is Censoring.LEFT_CENSORED:
            left_censored += 1
        elif censoring is Censoring.RIGHT_CENSORED:
            right_censored += 1
        else:
            exact += 1

        where = f"{path.name} genotype {variant_id}"
        records.append(VariantRecord(
            variant_id=variant_id,
            combo=rows[0].combo,
            fitness=_apply_transform(aggregated, config.transform, where),
            label_status=LabelStatus.MEASURED,
            censoring=censoring,
            raw_fitness=aggregated,
            n_replicates=(
                rows[0].n_replicates if (len(rows) == 1 and rows[0].n_replicates is not None)
                else len(values)
            ),
            replicate_values=values,
            stage=rows[0].stage,
            is_imputed=any(r.is_imputed for r in rows),
        ))

    subsampled_to: Optional[int] = None
    if config.max_records is not None and config.max_records < len(records):
        keep = sorted(random.Random(config.subsample_seed).sample(range(len(records)), config.max_records))
        records = [records[i] for i in keep]
        subsampled_to = config.max_records
        exact = sum(1 for r in records if r.censoring is Censoring.EXACT)
        left_censored = sum(1 for r in records if r.censoring is Censoring.LEFT_CENSORED)
        right_censored = sum(1 for r in records if r.censoring is Censoring.RIGHT_CENSORED)

    # The transformation and any direction flip belong on the design that is handed
    # downstream, so a reader never has to infer the scale from the values.
    out_design = design
    if config.transform not in ("none", "identity"):
        out_design = out_design.with_assay(
            type(design.assay)(**{**design.assay.to_dict(), "transformation": config.transform})
        )
        if config.transform in _DIRECTION_FLIPPING:
            flipped = (
                FitnessDirection.HIGHER_IS_BETTER
                if design.fitness_direction is FitnessDirection.LOWER_IS_BETTER
                else FitnessDirection.LOWER_IS_BETTER
            )
            out_design = out_design.with_fitness_direction(flipped)

    # Coverage excludes the wild type from both sides: it is a permitted genotype but not
    # a library member, so counting it as "missing" would make a complete library look
    # incomplete by exactly one.
    non_wild_type_space = design.combinatorial_size - 1
    present = len(records)
    report = IngestReport(
        rows_read=rows_read,
        rows_labelled=len(raw_rows),
        rows_absent_label=len(absent),
        wild_type_rows=wild_type_rows,
        duplicate_groups=duplicate_groups,
        duplicate_rows_collapsed=duplicate_rows_collapsed,
        duplicate_aggregation=config.duplicate_policy,
        conflict_tolerance=config.conflict_tolerance,
        exact=exact,
        left_censored=left_censored,
        right_censored=right_censored,
        imputed_rows=imputed_rows,
        transform=config.transform,
        selection_mode=design.selection_mode.value,
        fitness_direction=out_design.fitness_direction.value,
        combinatorial_size=design.combinatorial_size,
        combinations_present=present,
        combinations_missing=max(non_wild_type_space - present, 0),
        coverage_fraction=(present / non_wild_type_space) if non_wild_type_space else 1.0,
        subsampled_to=subsampled_to,
    )

    return IngestedDataset(
        design=out_design,
        records=tuple(records),
        absent_labels=tuple(sorted(absent, key=lambda r: r.variant_id)),
        report=report,
        source_path=str(path),
        source_sha256=sha,
        source_size_bytes=size,
    )


__all__ = [
    "ColumnMap",
    "IngestConfig",
    "IngestReport",
    "IngestedDataset",
    "VariantRecord",
    "ingest_dataset",
]
