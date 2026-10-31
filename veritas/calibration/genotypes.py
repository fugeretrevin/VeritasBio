"""Full-length genotype canonicalisation over a combinatorial library's positions.

A mutation string OMITS wild-type positions, so `V183A:F184A:V227A` is a four-site
genotype whose fourth site is wild type. Comparing such strings directly - or zipping
genotypes of unequal length - understates Hamming distances and has already inverted a
diversity conclusion in this repository (`scripts/analysis/README.md`, "Genotype
canonicalisation"). Everything here therefore canonicalises to a fixed-length tuple
over the library's declared positions and `hamming` refuses unequal lengths rather than
truncating.

This is the same contract as `scripts/analysis/diversity.py` and
`benchmark_panoptes_isolated.py:_combo_from_notation`, lifted into an importable module
so a third implementation of it does not have to be written per analysis.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

Combo = Tuple[str, ...]
Mutation = Tuple[int, str, str]


class GenotypeParseError(ValueError):
    """A variant identifier does not describe a member of the declared library."""


def parse_mutations(variant_id: str) -> Optional[List[Mutation]]:
    """`'V183N:F184R'` -> `[(183, 'V', 'N'), (184, 'F', 'R')]`.

    The wild type - written `WT` or as an empty identifier - is the empty list, which is
    a genotype, not a failure. Returns None for anything that is not substitution
    notation, so a caller can count malformed rows instead of crashing on them.
    """
    s = (variant_id or "").strip()
    if not s:
        return []
    if s.upper() == "WT":
        return []
    out: List[Mutation] = []
    for token in s.split(":"):
        if not token:
            continue
        if len(token) < 3:
            return None
        digits = "".join(c for c in token if c.isdigit())
        if not digits:
            return None
        wt_residue, mutant_residue = token[0], token[-1]
        if not (wt_residue.isalpha() and mutant_residue.isalpha()):
            return None
        # 'V183N' is wt + position + mutant; anything else in between is not notation.
        if token != f"{wt_residue}{digits}{mutant_residue}":
            return None
        out.append((int(digits), wt_residue, mutant_residue))
    return out


def canonical_combo(
    variant_id: str,
    positions: Sequence[int],
    parent_residues: str,
    *,
    strict: bool = False,
) -> Optional[Combo]:
    """The residue at EVERY library position, parent residue where the string is silent.

    The result always has `len(positions)` entries, which is what makes Hamming
    distances between two genotypes of the same library well defined.

    Returns None when the identifier touches a position outside the library, disagrees
    with the declared parent residue at a position, repeats a position, or is not
    substitution notation. With `strict=True` those cases raise `GenotypeParseError`
    instead - a dataset whose rows do not match the declared design is a configuration
    error, not noise.
    """
    if len(positions) != len(parent_residues):
        raise GenotypeParseError(
            f"library declares {len(positions)} positions but "
            f"{len(parent_residues)} parent residues"
        )
    parent: Dict[int, str] = dict(zip(positions, parent_residues))

    def reject(why: str) -> Optional[Combo]:
        if strict:
            raise GenotypeParseError(f"{variant_id!r}: {why}")
        return None

    mutations = parse_mutations(variant_id)
    if mutations is None:
        return reject("not substitution notation")
    seen: Dict[int, str] = {}
    for pos, wt_residue, mutant_residue in mutations:
        if pos not in parent:
            return reject(f"position {pos} is outside the library")
        if parent[pos] != wt_residue:
            return reject(
                f"position {pos} declares parent {parent[pos]!r} "
                f"but the row says {wt_residue!r}"
            )
        if pos in seen:
            return reject(f"position {pos} appears twice")
        seen[pos] = mutant_residue
    combo = tuple(seen.get(pos, parent[pos]) for pos in positions)
    assert len(combo) == len(positions), "canonical combo must be full length"
    return combo


def hamming(a: Sequence[str], b: Sequence[str]) -> int:
    """Sites at which two genotypes of the same library differ.

    Raises on unequal lengths. Never zip or truncate: that is the bug this module
    exists to prevent.
    """
    if len(a) != len(b):
        raise ValueError(f"hamming on unequal lengths: {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def depth(combo: Sequence[str], parent_residues: str) -> int:
    """Library positions at which `combo` differs from the parent. Depth 0 is the WT."""
    return hamming(combo, tuple(parent_residues))


def min_pairwise_hamming(combos: Iterable[Combo]) -> Optional[int]:
    """The closest two genotypes in a set, or None when there are fewer than two.

    Reported alongside a diversity ceiling because it says WHY the ceiling binds (or
    does not): a hit set whose minimum pairwise distance is already >= d cannot be
    thinned by a d-separation filter at all - which is exactly HER2's situation.
    """
    hist = pair_distance_histogram(combos)
    return min(hist) if hist else None


def pair_distance_histogram(combos: Iterable[Combo]) -> Dict[int, int]:
    """Counts of unordered pairs at each Hamming distance over a genotype set.

    Vectorised in blocks because a 1,655-hit tier is 1.37 million pairs. Lengths are
    checked up front and a ragged set raises rather than being padded or truncated.
    """
    items = list(combos)
    if len(items) < 2:
        return {}
    width = len(items[0])
    if any(len(c) != width for c in items):
        raise ValueError("pair_distance_histogram on genotypes of unequal length")

    matrix = np.array([list(c) for c in items], dtype="<U1")
    n = matrix.shape[0]
    counts = np.zeros(width + 1, dtype=np.int64)
    block = max(1, min(n, 4_000_000 // max(1, n * width)))
    for start in range(0, n, block):
        stop = min(n, start + block)
        # Distances from rows [start:stop] to every row at a strictly greater index.
        d = (matrix[start:stop, None, :] != matrix[None, :, :]).sum(axis=2)
        rows = np.arange(start, stop)[:, None]
        cols = np.arange(n)[None, :]
        counts += np.bincount(d[cols > rows].ravel(), minlength=width + 1)
    return {int(k): int(v) for k, v in enumerate(counts) if v}


def pairs_closer_than(combos: Iterable[Combo], d: int) -> int:
    """Unordered pairs at Hamming distance < d: what a d-separation filter can act on.

    Zero means the greedy d-separation rule removes nothing from the set, so the
    diversity-adjusted ceiling equals the raw count and the endpoint is degenerate.
    """
    hist = pair_distance_histogram(combos)
    return sum(v for k, v in hist.items() if k < d)
