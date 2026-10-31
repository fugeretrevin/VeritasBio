"""
The canonical schema for a combinatorial protein library.

A `LibraryDesign` is the single description a new library needs in order to be
benchmarked: the target and assay it belongs to, the wild-type reference it is defined
against, which positions vary, what each position is PERMITTED to carry, which direction
of fitness is better, and whether candidates come from an enumerated pool or are
generated under mutation constraints.

REQUIRED (a library cannot be ingested without these):
    target_id, assay_id, wild_type_sequence, positions (position + wild-type residue +
    allowed alphabet per position), fitness_direction, selection_mode.

OPTIONAL (recorded when known, never guessed - see `AssayMeta`):
    assay units and transformation, detection floor/ceiling, replicate count, source
    publication / DOI / dataset version, measurement provenance, free-text notes.

Every encoder here takes its library explicitly. There is no module-level "active
library": two designs are usable in one process, which is what a cross-protein study
needs and what a set of module globals cannot give.

THE PARSING RULE that everything else depends on: a mutation string OMITS wild-type
positions. `V183A:F184A:V227A` is a FOUR-site genotype whose fourth site is wild type.
Every genotype is canonicalised to full length over the library's positions, filling the
wild-type residue where a site is omitted, BEFORE anything compares two of them. Unequal
lengths are an error, never a zip or a truncation.
"""

from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .errors import (
    AlphabetError,
    GenotypeLengthError,
    IndelError,
    PositionOutsideLibraryError,
    SchemaError,
    WildTypeMismatchError,
)

#: The 20 proteinogenic amino acids. A position's permissible alphabet is a subset of
#: this, DECLARED by the library design - it is never inferred from observed sequences,
#: because an incomplete assay would then silently shrink the library.
AA_ALPHABET: str = "ACDEFGHIKLMNPQRSTVWY"

_SUBSTITUTION = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")
_TRUNCATION = re.compile(r"^([A-Za-z])(\d+)$")
_INDEL_MARKERS = ("INS", "DEL", "DUP", "FS", "EXT", "*", "-", "_", ">")
#: Spellings of "this row is the unmutated reference", not a library member.
WILD_TYPE_TOKENS = ("", "WT", "WILDTYPE", "WILD_TYPE", "NONE")


class FitnessDirection(str, Enum):
    """Which direction of the measured value is better. Declared, never inferred."""

    HIGHER_IS_BETTER = "higher_is_better"
    LOWER_IS_BETTER = "lower_is_better"


class SelectionMode(str, Enum):
    """
    How candidates are produced, declared per library and kept distinct.

    POOL_RESTRICTED selects from an enumerated feasible library. It is the default for
    retrospective benchmarking, because a sparse measured sample of a combinatorial
    product cannot score a proposal nobody assayed.

    GENERATIVE proposes variants under the library's mutation constraints. Success in
    pool mode does NOT validate unconstrained generation: in pool mode the engine is
    handed the feasible set and only has to rank it, which is a strictly easier problem.
    """

    POOL_RESTRICTED = "pool_restricted"
    GENERATIVE = "generative"


class LabelStatus(str, Enum):
    """Whether a row carries a measurement at all."""

    MEASURED = "measured"
    #: No label exists. The file may still hold a numeric placeholder (the HER2 held-out
    #: split carries DMS_score exactly 0.0 on all 25,114 of its rows); that number is
    #: not a measurement and must never enter a fitness table.
    ABSENT = "absent"


class Censoring(str, Enum):
    """How a measured value relates to the assay's detection limits."""

    EXACT = "exact"
    #: At or below the detection floor: the true value is somewhere at or below this.
    LEFT_CENSORED = "left_censored"
    #: At or above the saturation ceiling.
    RIGHT_CENSORED = "right_censored"
    #: There is no measurement to censor.
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True)
class VariablePosition:
    """One variable site: where it is, what the wild type carries, what it may carry."""

    position: int            # 1-indexed into the wild-type reference sequence
    wild_type_aa: str
    allowed_aa: Tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.position, int) or self.position < 1:
            raise SchemaError(f"position must be a 1-indexed integer, got {self.position!r}")
        if len(self.wild_type_aa) != 1:
            raise SchemaError(
                f"position {self.position}: wild_type_aa must be one residue, "
                f"got {self.wild_type_aa!r}"
            )
        if self.wild_type_aa not in AA_ALPHABET:
            raise AlphabetError(
                f"position {self.position}: wild-type residue {self.wild_type_aa!r} is not "
                f"one of the 20 proteinogenic amino acids"
            )
        if not self.allowed_aa:
            raise SchemaError(
                f"position {self.position}: allowed_aa is empty. The permissible alphabet is "
                f"part of the library design and must be declared; it is never inferred from "
                f"the observed sequences, because an incomplete assay would shrink it silently."
            )
        if len(set(self.allowed_aa)) != len(self.allowed_aa):
            raise SchemaError(f"position {self.position}: allowed_aa has repeated residues")
        illegal = [aa for aa in self.allowed_aa if aa not in AA_ALPHABET or len(aa) != 1]
        if illegal:
            raise AlphabetError(
                f"position {self.position}: allowed_aa contains non-amino-acid entries {illegal}"
            )


@dataclass(frozen=True)
class AssayMeta:
    """
    OPTIONAL metadata. Every field defaults to unknown rather than to a plausible value:
    a guessed detection floor or a guessed unit is worse than an admitted gap, because
    it is indistinguishable from a recorded fact downstream.
    """

    units: Optional[str] = None
    #: A transformation already applied to the stored values, e.g. "log10".
    transformation: Optional[str] = None
    #: Values at or below this are LEFT_CENSORED; values strictly below it contradict
    #: the declared limit and are rejected.
    detection_floor: Optional[float] = None
    #: Values at or above this are RIGHT_CENSORED.
    detection_ceiling: Optional[float] = None
    replicate_count: Optional[int] = None
    measurement_provenance: str = "unspecified"
    source_publication: Optional[str] = None
    source_doi: Optional[str] = None
    dataset_version: Optional[str] = None
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass(frozen=True)
class LibraryDesign:
    """A combinatorial library, described completely enough to ingest and benchmark."""

    # ---- REQUIRED ----------------------------------------------------------------
    target_id: str
    assay_id: str
    wild_type_sequence: str
    positions: Tuple[VariablePosition, ...]
    fitness_direction: FitnessDirection
    selection_mode: SelectionMode
    # ---- OPTIONAL ----------------------------------------------------------------
    assay: AssayMeta = field(default_factory=AssayMeta)
    #: Only meaningful in GENERATIVE mode; defaults to the number of variable sites.
    max_mutation_depth: Optional[int] = None
    label: str = ""

    # -- construction-time validation ---------------------------------------------
    def __post_init__(self) -> None:
        if not self.target_id or not self.assay_id:
            raise SchemaError("target_id and assay_id are required")
        if not self.wild_type_sequence:
            raise SchemaError("wild_type_sequence is required")
        bad = sorted({c for c in self.wild_type_sequence if c not in AA_ALPHABET})
        if bad:
            raise AlphabetError(
                f"wild_type_sequence for {self.target_id} contains non-amino-acid characters {bad}"
            )
        if not self.positions:
            raise SchemaError(f"{self.target_id}: a library must declare at least one variable position")
        if not isinstance(self.fitness_direction, FitnessDirection):
            raise SchemaError(f"fitness_direction must be a FitnessDirection, got {self.fitness_direction!r}")
        if not isinstance(self.selection_mode, SelectionMode):
            raise SchemaError(f"selection_mode must be a SelectionMode, got {self.selection_mode!r}")

        seen: set = set()
        previous = 0
        for site in self.positions:
            if site.position in seen:
                raise SchemaError(f"{self.target_id}: position {site.position} is declared twice")
            if site.position <= previous:
                raise SchemaError(
                    f"{self.target_id}: positions must be declared in strictly increasing order "
                    f"(saw {site.position} after {previous}); the order fixes the genotype's slot "
                    f"order, so it cannot be left to the caller"
                )
            seen.add(site.position)
            previous = site.position
            if site.position > len(self.wild_type_sequence):
                raise SchemaError(
                    f"{self.target_id}: position {site.position} is outside the "
                    f"{len(self.wild_type_sequence)}-residue wild-type reference"
                )
            reference = self.wild_type_sequence[site.position - 1]
            if site.wild_type_aa != reference:
                raise WildTypeMismatchError(
                    f"{self.target_id}: the design declares {site.wild_type_aa}{site.position} but "
                    f"the wild-type reference carries {reference} there. Either the reference "
                    f"sequence or the residue numbering is wrong; refusing to align against it."
                )
            if site.wild_type_aa not in site.allowed_aa:
                raise SchemaError(
                    f"{self.target_id}: position {site.position} allows {list(site.allowed_aa)} but "
                    f"not its own wild-type residue {site.wild_type_aa!r}. A mutation string omits "
                    f"wild-type sites, so the wild-type residue must be representable at every site."
                )
        if self.max_mutation_depth is not None and not (1 <= self.max_mutation_depth <= len(self.positions)):
            raise SchemaError(
                f"{self.target_id}: max_mutation_depth {self.max_mutation_depth} must be between 1 "
                f"and the {len(self.positions)} declared positions"
            )

        # Derived lookups, built once. Cached on a frozen instance deliberately: they are
        # pure functions of the declared fields.
        object.__setattr__(self, "_slot_of", {s.position: i for i, s in enumerate(self.positions)})
        object.__setattr__(self, "_allowed_sets", tuple(frozenset(s.allowed_aa) for s in self.positions))
        object.__setattr__(
            self, "_encoding_index",
            tuple({aa: i for i, aa in enumerate(s.allowed_aa)} for s in self.positions),
        )
        object.__setattr__(self, "_wild_type_combo", "".join(s.wild_type_aa for s in self.positions))

    # -- shape --------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        return len(self.positions)

    @property
    def site_positions(self) -> Tuple[int, ...]:
        return tuple(s.position for s in self.positions)

    @property
    def wild_type_combo(self) -> str:
        """The full-length genotype of the unmutated reference over the variable sites."""
        return getattr(self, "_wild_type_combo")

    @property
    def one_hot_width(self) -> int:
        return sum(len(s.allowed_aa) for s in self.positions)

    @property
    def combinatorial_size(self) -> int:
        """How many genotypes the declared design permits, wild type included."""
        total = 1
        for site in self.positions:
            total *= len(site.allowed_aa)
        return total

    @property
    def is_saturated_alphabet(self) -> bool:
        return all(len(s.allowed_aa) == len(AA_ALPHABET) for s in self.positions)

    @property
    def is_pool_restricted(self) -> bool:
        return self.selection_mode is SelectionMode.POOL_RESTRICTED

    @property
    def is_generative(self) -> bool:
        return self.selection_mode is SelectionMode.GENERATIVE

    @property
    def effective_max_mutation_depth(self) -> int:
        return self.max_mutation_depth if self.max_mutation_depth is not None else self.n_sites

    def mutation_constraints(self) -> Dict[str, Any]:
        """The constraints a GENERATIVE proposer must respect. Also useful as a record."""
        return {
            "positions": list(self.site_positions),
            "allowed_aa": {str(s.position): list(s.allowed_aa) for s in self.positions},
            "max_mutation_depth": self.effective_max_mutation_depth,
            "wild_type_combo": self.wild_type_combo,
            "selection_mode": self.selection_mode.value,
        }

    # -- derivation ---------------------------------------------------------------
    def with_assay(self, assay: AssayMeta) -> "LibraryDesign":
        return dataclasses.replace(self, assay=assay)

    def with_fitness_direction(self, direction: FitnessDirection) -> "LibraryDesign":
        return dataclasses.replace(self, fitness_direction=direction)

    # -- genotype handling --------------------------------------------------------
    def assert_full_length(self, combo: str, context: str = "genotype") -> None:
        """
        The guard that replaces zipping two genotypes of different length.

        Comparing site i of a 3-character genotype against site i of a 4-character one
        is always wrong and always silent, so unequal length is raised here instead.
        """
        if not isinstance(combo, str):
            raise GenotypeLengthError(f"{context} must be a string, got {type(combo).__name__}")
        if len(combo) != self.n_sites:
            raise GenotypeLengthError(
                f"{context} {combo!r} has {len(combo)} site(s) but library {self.assay_id} declares "
                f"{self.n_sites} ({list(self.site_positions)}). A mutation string omits wild-type "
                f"positions, so canonicalise it to full length with combo_from_mutations() before "
                f"comparing; genotypes are never zipped or truncated to match."
            )

    def _assert_alphabet(self, combo: str, context: str = "genotype") -> None:
        allowed: Tuple[frozenset, ...] = getattr(self, "_allowed_sets")
        for slot, residue in enumerate(combo):
            if residue not in allowed[slot]:
                site = self.positions[slot]
                raise AlphabetError(
                    f"{context} {combo!r}: residue {residue!r} is not permitted at position "
                    f"{site.position}, which allows {sorted(site.allowed_aa)}. The permissible "
                    f"alphabet comes from the library design, not from the observed data."
                )

    def validate_combo(self, combo: str, context: str = "genotype") -> str:
        self.assert_full_length(combo, context)
        self._assert_alphabet(combo, context)
        return combo

    def canonical_variant_id(self, combo: str) -> str:
        """
        Position-ordered mutation notation for a full-length genotype, wild-type sites
        omitted. The wild-type genotype maps to the empty string: it carries no mutations
        and is not a library member.
        """
        self.validate_combo(combo, "genotype")
        return ":".join(
            f"{site.wild_type_aa}{site.position}{combo[slot]}"
            for slot, site in enumerate(self.positions)
            if combo[slot] != site.wild_type_aa
        )

    def combo_from_mutations(self, notation: str) -> str:
        """
        Canonicalises a mutation string to a FULL-LENGTH genotype over this library's
        positions, filling the wild-type residue at every site the string omits.

        Rejects, rather than repairs: an indel or truncation (`IndelError`), a position
        outside the design (`PositionOutsideLibraryError`), a stated wild-type residue
        that disagrees with the reference (`WildTypeMismatchError`), a residue outside
        the declared alphabet (`AlphabetError`), and two tokens at one position or an
        unparseable token (`SchemaError`).
        """
        if notation is None:
            raise SchemaError("mutation string is missing")
        text = notation.strip()
        if text.upper() in WILD_TYPE_TOKENS:
            return self.wild_type_combo

        slot_of: Mapping[int, int] = getattr(self, "_slot_of")
        allowed: Tuple[frozenset, ...] = getattr(self, "_allowed_sets")
        residues = list(self.wild_type_combo)
        touched: Dict[int, str] = {}

        for raw in text.split(":"):
            token = raw.strip()
            if not token:
                raise SchemaError(f"{notation!r}: empty mutation token")
            match = _SUBSTITUTION.match(token)
            if match is None:
                upper = token.upper()
                if any(marker in upper for marker in _INDEL_MARKERS) or _TRUNCATION.match(token):
                    raise IndelError(
                        f"{notation!r}: token {token!r} is an insertion, deletion, frameshift or "
                        f"truncation, not a substitution. Length-changing events shift every "
                        f"position index after them, so library {self.assay_id} refuses the row "
                        f"instead of reading it as a substitution."
                    )
                raise SchemaError(
                    f"{notation!r}: token {token!r} is not substitution notation "
                    f"(expected e.g. 'V183A')"
                )
            wt_aa, position_text, mut_aa = match.group(1).upper(), match.group(2), match.group(3).upper()
            position = int(position_text)
            if position not in slot_of:
                raise PositionOutsideLibraryError(
                    f"{notation!r}: position {position} is not one of library {self.assay_id}'s "
                    f"declared variable positions {list(self.site_positions)}. The inferred "
                    f"mutation positions must agree with the declared library design."
                )
            slot = slot_of[position]
            site = self.positions[slot]
            if wt_aa != site.wild_type_aa:
                raise WildTypeMismatchError(
                    f"{notation!r}: token {token!r} states wild-type {wt_aa} at position "
                    f"{position}, but {self.target_id}'s reference carries {site.wild_type_aa}. "
                    f"Refusing to align against the wrong wild-type sequence."
                )
            if mut_aa not in allowed[slot]:
                raise AlphabetError(
                    f"{notation!r}: residue {mut_aa!r} is not permitted at position {position}, "
                    f"which allows {sorted(site.allowed_aa)}."
                )
            if slot in touched and touched[slot] != mut_aa:
                raise SchemaError(
                    f"{notation!r}: position {position} is substituted twice, to {touched[slot]!r} "
                    f"and {mut_aa!r}. One genotype cannot carry two residues at one site."
                )
            touched[slot] = mut_aa
            residues[slot] = mut_aa
        return "".join(residues)

    def combo_from_sequence(self, sequence: str) -> str:
        """
        Derives a full-length genotype by aligning a mutated sequence against the
        wild-type reference.

        A length difference is an INDEL and is refused: it is not a set of
        substitutions, and reading it as one misassigns every downstream position.
        """
        if sequence is None:
            raise SchemaError("mutated sequence is missing")
        text = sequence.strip().upper()
        if len(text) != len(self.wild_type_sequence):
            raise IndelError(
                f"mutated sequence is {len(text)} residues but {self.target_id}'s wild-type "
                f"reference is {len(self.wild_type_sequence)}. A length difference is an "
                f"insertion or deletion, not a substitution set; refusing to align it."
            )
        slot_of: Mapping[int, int] = getattr(self, "_slot_of")
        residues = list(self.wild_type_combo)
        for index, (observed, reference) in enumerate(zip(text, self.wild_type_sequence)):
            if observed == reference:
                continue
            position = index + 1
            if position not in slot_of:
                raise PositionOutsideLibraryError(
                    f"mutated sequence substitutes {reference}{position}{observed}, which is "
                    f"outside library {self.assay_id}'s declared positions "
                    f"{list(self.site_positions)}."
                )
            residues[slot_of[position]] = observed
        combo = "".join(residues)
        self._assert_alphabet(combo, "genotype from sequence")
        return combo

    def apply_to_sequence(self, combo: str) -> str:
        """The full-length protein sequence for a genotype."""
        self.validate_combo(combo)
        residues = list(self.wild_type_sequence)
        for slot, site in enumerate(self.positions):
            residues[site.position - 1] = combo[slot]
        return "".join(residues)

    def mutation_depth(self, combo: str) -> int:
        """How many sites differ from wild type. Full length is enforced first."""
        self.assert_full_length(combo)
        return self.hamming_distance(combo, self.wild_type_combo)

    def hamming_distance(self, left: str, right: str) -> int:
        """Site-wise distance between two FULL-LENGTH genotypes."""
        self.assert_full_length(left, "left genotype")
        self.assert_full_length(right, "right genotype")
        return sum(1 for a, b in zip(left, right) if a != b)

    def one_hot(self, combos: Sequence[str]) -> np.ndarray:
        """
        One-hot encodes genotypes as sum(len(allowed_aa)) binary features, in the
        declared position order. 80 columns for a saturated four-site library, 200 for a
        saturated ten-site one, 4 for a two-site binary one.
        """
        encoding: Tuple[Dict[str, int], ...] = getattr(self, "_encoding_index")
        offsets: List[int] = []
        running = 0
        for site in self.positions:
            offsets.append(running)
            running += len(site.allowed_aa)
        matrix = np.zeros((len(combos), running), dtype=np.float32)
        for row, combo in enumerate(combos):
            self.assert_full_length(combo, f"genotype at row {row}")
            for slot, residue in enumerate(combo):
                column = encoding[slot].get(residue)
                if column is None:
                    site = self.positions[slot]
                    raise AlphabetError(
                        f"row {row}: residue {residue!r} is not permitted at position "
                        f"{site.position}, which allows {sorted(site.allowed_aa)}."
                    )
                matrix[row, offsets[slot] + column] = 1.0
        return matrix

    def enumerate_combos(self) -> Iterable[str]:
        """
        Every genotype the declared design permits, in lexicographic slot order.

        Only usable on a small design: `combinatorial_size` is 20^10 for the HER2 CDR-H3
        library, which is why that library is POOL_RESTRICTED.
        """
        import itertools

        for parts in itertools.product(*(site.allowed_aa for site in self.positions)):
            yield "".join(parts)

    def to_dict(self) -> Dict[str, Any]:
        """A record of the design, for a run manifest."""
        return {
            "target_id": self.target_id,
            "assay_id": self.assay_id,
            "label": self.label,
            "wild_type_length": len(self.wild_type_sequence),
            "positions": [
                {
                    "position": s.position,
                    "wild_type_aa": s.wild_type_aa,
                    "allowed_aa": list(s.allowed_aa),
                }
                for s in self.positions
            ],
            "wild_type_combo": self.wild_type_combo,
            "fitness_direction": self.fitness_direction.value,
            "selection_mode": self.selection_mode.value,
            "max_mutation_depth": self.effective_max_mutation_depth,
            "combinatorial_size": self.combinatorial_size,
            "one_hot_width": self.one_hot_width,
            "assay": self.assay.to_dict(),
        }


__all__ = [
    "AA_ALPHABET",
    "AssayMeta",
    "Censoring",
    "FitnessDirection",
    "LabelStatus",
    "LibraryDesign",
    "SelectionMode",
    "VariablePosition",
    "WILD_TYPE_TOKENS",
]
