"""
Ingestion errors.

Every one of these is a REFUSAL, not a warning. The alternative - repairing a record so
the run can continue - is what turns a dataset defect into a published number: an indel
silently read as a substitution, a conflicting replicate silently averaged, an absent
label silently read as a measured zero. Each error names the offending row and what was
expected, so the fix is a change to the data or the configuration rather than a guess.
"""

from __future__ import annotations


class DatasetError(ValueError):
    """Base class for every dataset schema or ingestion refusal."""


class SchemaError(DatasetError):
    """The design, the configuration or a row violates the canonical schema."""


class WildTypeMismatchError(DatasetError):
    """A declared or stated wild-type residue disagrees with the reference sequence."""


class AlphabetError(DatasetError):
    """A residue is outside the alphabet DECLARED for that position.

    Never raised on the basis of what happens to be observed: an incomplete assay does
    not get to narrow the library's permissible alphabet.
    """


class IndelError(DatasetError):
    """An insertion, deletion, frameshift or truncation was found.

    These are length-changing events. Treating one as a substitution corrupts every
    position index after it, so they are rejected rather than coerced.
    """


class GenotypeLengthError(DatasetError):
    """A genotype is not full length over the library's variable positions.

    Raised instead of zipping or truncating two genotypes of unequal length, which
    silently compares site i of one against site i of the other.
    """


class PositionOutsideLibraryError(DatasetError):
    """A substitution sits at a position the library design does not declare."""


class DuplicateGenotypeError(DatasetError):
    """The same genotype appears more than once and no aggregation policy was chosen."""


class ConflictingMeasurementError(DatasetError):
    """Repeated measurements of one genotype disagree beyond the declared tolerance."""


class LabelPolicyError(DatasetError):
    """A row's label availability contradicts the declared label policy."""


__all__ = [
    "DatasetError",
    "SchemaError",
    "WildTypeMismatchError",
    "AlphabetError",
    "IndelError",
    "GenotypeLengthError",
    "PositionOutsideLibraryError",
    "DuplicateGenotypeError",
    "ConflictingMeasurementError",
    "LabelPolicyError",
]
