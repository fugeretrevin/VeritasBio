"""
veritas/diagnostics/genotype.py - canonicalising a mutation string onto a library.

Why this module exists at all
----------------------------
Mutation strings OMIT wild-type positions. `V183A:F184A:V227A` over the four-site
TrpB library is a FOUR-site genotype whose fourth site (228) is wild type - it is
NOT a three-site genotype, and it is not comparable to a three-token string from a
different library by position order. Getting that wrong once already cost us an
analysis: zipping two notations of unequal token count silently compared position
183 against 184 and reported a mutation depth that was too low.

So every depth or coverage number in this package goes through
`LibrarySpec.canonical_genotype`, which returns a string of exactly
`len(positions)` residues with the parent residue filled in wherever the notation
is silent. Unequal lengths are refused, never zipped and never truncated.

Mutation depth is defined as the number of LIBRARY POSITIONS whose residue differs
from the parent (wild type). It is a property of the genotype, not of the notation's
token count.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Tuple

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_SET = frozenset(AA_ALPHABET)

WILD_TYPE_NOTATIONS = frozenset({"", "WT", "wt", "Wt", "wild_type", "wildtype"})


class GenotypeParseError(ValueError):
    """A notation or genotype could not be placed on the library without guessing."""


@dataclass(frozen=True, slots=True)
class LibrarySpec:
    """
    The ordered set of saturated positions of a combinatorial library, with the parent
    (wild-type) residue at each one.

    `positions` is strictly increasing so that a canonical genotype has one, stable
    character order regardless of the order tokens appear in a notation.
    """

    positions: Tuple[int, ...]
    parent_residues: str

    def __post_init__(self) -> None:
        if not self.positions:
            raise ValueError("a library must declare at least one position")
        if len(self.positions) != len(self.parent_residues):
            raise ValueError(
                f"{len(self.positions)} position(s) but {len(self.parent_residues)} parent "
                "residue(s); the two must correspond one to one"
            )
        if list(self.positions) != sorted(set(self.positions)):
            raise ValueError(
                f"positions must be unique and strictly increasing, got {self.positions}"
            )
        bad = sorted(set(self.parent_residues) - _AA_SET)
        if bad:
            raise ValueError(f"parent residues {bad} are not standard amino acids")

    # -- identity --------------------------------------------------------------------
    @property
    def width(self) -> int:
        """Number of library positions; the canonical genotype length."""
        return len(self.positions)

    @property
    def max_depth(self) -> int:
        """The deepest genotype the library can express."""
        return len(self.positions)

    def slot_of(self, position: int) -> int:
        try:
            return self.positions.index(position)
        except ValueError:
            raise GenotypeParseError(
                f"position {position} is outside the library {list(self.positions)}"
            ) from None

    # -- parsing ---------------------------------------------------------------------
    def canonical_genotype(self, notation: str) -> str:
        """
        Full-length genotype for a mutation notation, parent residue filled in wherever
        the notation omits a position.

        Raises GenotypeParseError for anything that cannot be placed without guessing:
        a malformed token, a position outside the library, a wild-type residue that
        disagrees with the parent, a non-standard substituted residue, or the same
        position named twice.
        """
        if not isinstance(notation, str):
            raise GenotypeParseError(
                f"expected a mutation notation string, got {type(notation).__name__}"
            )
        text = notation.strip()
        if text in WILD_TYPE_NOTATIONS:
            return self.parent_residues

        residues = list(self.parent_residues)
        touched: Dict[int, str] = {}
        for raw_token in text.split(":"):
            token = raw_token.strip()
            if len(token) < 3:
                raise GenotypeParseError(f"token {token!r} in {notation!r} is not a mutation")
            wt_aa, mut_aa, position_text = token[0], token[-1], token[1:-1]
            if not position_text.isdigit():
                raise GenotypeParseError(
                    f"token {token!r} in {notation!r} has no integer position"
                )
            position = int(position_text)
            slot = self.slot_of(position)
            if wt_aa != self.parent_residues[slot]:
                raise GenotypeParseError(
                    f"token {token!r} declares wild type {wt_aa!r} at position {position}, "
                    f"but the library's parent residue there is "
                    f"{self.parent_residues[slot]!r}"
                )
            if mut_aa not in _AA_SET:
                raise GenotypeParseError(
                    f"token {token!r} substitutes {mut_aa!r}, which is not a standard "
                    "amino acid"
                )
            if position in touched:
                raise GenotypeParseError(
                    f"{notation!r} names position {position} twice "
                    f"({touched[position]!r} and {mut_aa!r})"
                )
            touched[position] = mut_aa
            residues[slot] = mut_aa

        genotype = "".join(residues)
        if len(genotype) != self.width:
            # Unreachable by construction; asserted because a silent length drift here
            # is exactly the bug this module exists to prevent.
            raise GenotypeParseError(
                f"canonicalisation produced {len(genotype)} residues for a "
                f"{self.width}-site library"
            )
        return genotype

    def depth_of_genotype(self, genotype: str) -> int:
        """
        Mutation depth of an already-canonical genotype.

        Refuses an unequal-length genotype instead of zipping it against the parent,
        which would silently compare mismatched positions.
        """
        if not isinstance(genotype, str):
            raise GenotypeParseError(
                f"expected a genotype string, got {type(genotype).__name__}"
            )
        if len(genotype) != self.width:
            raise GenotypeParseError(
                f"genotype {genotype!r} has {len(genotype)} residues but the library has "
                f"{self.width} positions; refusing to compare unequal lengths"
            )
        return sum(1 for a, b in zip(genotype, self.parent_residues) if a != b)

    def mutation_depth(self, notation: str) -> int:
        """Number of library positions at which `notation` differs from the parent."""
        return self.depth_of_genotype(self.canonical_genotype(notation))

    def to_dict(self) -> Mapping[str, object]:
        return {
            "positions": list(self.positions),
            "parent_residues": self.parent_residues,
            "width": self.width,
        }


def infer_library_from_variant_ids(variant_ids: Iterable[str]) -> LibrarySpec:
    """
    Derives the library definition from the variant ids themselves.

    Every token carries its own wild-type residue (`V183A` says position 183 is V in
    the parent), so a set of ids over one saturated library determines both the
    positions and the parent residues. This keeps the diagnostics free of any
    dependence on module-level library constants, which are being moved elsewhere.

    Raises GenotypeParseError when two ids disagree about the parent residue at a
    position, or when no id names any position at all.
    """
    parents: Dict[int, str] = {}
    for variant_id in variant_ids:
        if not isinstance(variant_id, str):
            raise GenotypeParseError(
                f"expected variant id strings, got {type(variant_id).__name__}"
            )
        text = variant_id.strip()
        if text in WILD_TYPE_NOTATIONS:
            continue
        for raw_token in text.split(":"):
            token = raw_token.strip()
            if len(token) < 3 or not token[1:-1].isdigit():
                raise GenotypeParseError(f"token {token!r} in {variant_id!r} is not a mutation")
            position, wt_aa = int(token[1:-1]), token[0]
            if wt_aa not in _AA_SET:
                raise GenotypeParseError(
                    f"token {token!r} declares non-standard wild type {wt_aa!r}"
                )
            known = parents.get(position)
            if known is not None and known != wt_aa:
                raise GenotypeParseError(
                    f"conflicting parent residues at position {position}: "
                    f"{known!r} and {wt_aa!r}"
                )
            parents[position] = wt_aa

    if not parents:
        raise GenotypeParseError("no mutation tokens found; cannot infer a library")
    ordered = tuple(sorted(parents))
    return LibrarySpec(
        positions=ordered,
        parent_residues="".join(parents[p] for p in ordered),
    )
