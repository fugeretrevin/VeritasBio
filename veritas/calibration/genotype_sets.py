"""Do two landscapes of one protein actually share a genotype set?

H2 (§6.7) needs a within-genotype-set paired contrast: seed s must draw the SAME 25
initial genotypes on both readouts, which is only possible when the two tables hold the
same genotypes. eqFP611 blue and red do, row for row. CR6261 H1 and H9 hold 1,886 and
1,841 variants, so they cannot - protocol open item O-11 asks how large the shared set
actually is, and this module measures it instead of assuming it.

Comparison is over canonical full-length combos, never over mutation strings.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import AbstractSet, Dict, Iterable

from veritas.calibration.genotypes import Combo, depth


@dataclass(frozen=True)
class SetComparison:
    """Sizes of two genotype sets, their intersection and both remainders."""

    n_left: int
    n_right: int
    n_shared: int
    n_left_only: int
    n_right_only: int

    @property
    def identical(self) -> bool:
        return self.n_left_only == 0 and self.n_right_only == 0 and self.n_left == self.n_right

    @property
    def jaccard(self) -> float:
        union = self.n_shared + self.n_left_only + self.n_right_only
        return self.n_shared / union if union else 0.0

    @property
    def shared_share_left(self) -> float:
        return self.n_shared / self.n_left if self.n_left else 0.0

    @property
    def shared_share_right(self) -> float:
        return self.n_shared / self.n_right if self.n_right else 0.0


def compare_genotype_sets(left: AbstractSet[Combo], right: AbstractSet[Combo]) -> SetComparison:
    """Intersection size and both exclusive remainders for two genotype sets."""
    left_set, right_set = set(left), set(right)
    shared = left_set & right_set
    return SetComparison(
        n_left=len(left_set),
        n_right=len(right_set),
        n_shared=len(shared),
        n_left_only=len(left_set - right_set),
        n_right_only=len(right_set - left_set),
    )


@dataclass(frozen=True)
class SharedSetViability:
    """Is the shared set large enough, and unbiased enough, to pair on?"""

    comparison: SetComparison
    parent_residues: str
    depth_profile_shared: Dict[int, int]
    depth_profile_left_only: Dict[int, int]
    depth_profile_right_only: Dict[int, int]
    budget: int
    tier_fraction: float

    @property
    def tier_on_shared(self) -> int:
        """Hits in the shared set at the target's tier fraction, by construction."""
        return int(math.ceil(self.tier_fraction * self.comparison.n_shared))

    @property
    def budget_share_of_shared(self) -> float:
        return self.budget / self.comparison.n_shared if self.comparison.n_shared else float("inf")

    @property
    def budget_exceeds_tier(self) -> bool:
        return self.budget >= self.tier_on_shared


def depth_profile(combos: Iterable[Combo], parent_residues: str) -> Dict[int, int]:
    """`{depth: count}` over a genotype set, for checking whether a remainder is biased."""
    out: Dict[int, int] = {}
    for combo in combos:
        d = depth(combo, parent_residues)
        out[d] = out.get(d, 0) + 1
    return dict(sorted(out.items()))


def assess_shared_set(
    left: AbstractSet[Combo],
    right: AbstractSet[Combo],
    *,
    parent_residues: str,
    budget: int = 125,
    tier_fraction: float = 0.01,
) -> SharedSetViability:
    """The O-11 answer: intersection size, both remainders, and their depth profiles.

    The depth profiles matter because an intersection can be large and still unusable:
    if the variants missing from one side sit at particular depths, the shared set is a
    depth-biased subsample of the library rather than a smaller copy of it.
    """
    left_set, right_set = set(left), set(right)
    shared = left_set & right_set
    return SharedSetViability(
        comparison=compare_genotype_sets(left_set, right_set),
        parent_residues=parent_residues,
        depth_profile_shared=depth_profile(shared, parent_residues),
        depth_profile_left_only=depth_profile(left_set - right_set, parent_residues),
        depth_profile_right_only=depth_profile(right_set - left_set, parent_residues),
        budget=budget,
        tier_fraction=tier_fraction,
    )
