"""Is diversity-adjusted yield at Hamming d >= 2 a usable endpoint on this target?

The pre-registered primary endpoint (§2 of docs/preregistration-cross-protein.md) sorts
hits by fitness descending and keeps a hit only when it is at least `d` from everything
already kept. On the saturated four-site libraries that adjustment is the MORE
discriminating endpoint: TrpB collapses 1,600 raw hits to 166 at d >= 2. On a sparse
library it can collapse nothing at all - HER2's 113 raw hits give a d >= 2 ceiling of
113, numerically identical to the raw count for every arm, because the minimum pairwise
Hamming distance across its whole top-1% set is already 2.

An endpoint that cannot separate two arms is not a conservative choice; it silently
substitutes the secondary endpoint for the primary one. This module is the precondition
check: compute the ceilings, compare them to the raw count, and return a verdict per
target BEFORE the target is used as a primary endpoint.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from veritas.calibration.genotypes import Combo, hamming

DEGENERATE = "DEGENERATE"
NEAR_DEGENERATE = "NEAR-DEGENERATE"
USABLE = "USABLE"

#: A diversity adjustment that removes fewer than this share of the hit set carries no
#: usable information about how concentrated an arm's hits are.
INFORMATIVENESS_FLOOR = 0.05

Hit = Tuple[float, Combo]


def greedy_diverse(hits: Iterable[Hit], d: int) -> int:
    """The protocol's rule: fitness-descending, keep a hit only if >= d from all kept.

    Applied strictly after measurement, so it can never advantage an arm that knew
    about it. Input need not be sorted; it is sorted here so the rule cannot depend on
    the caller's ordering.
    """
    kept: List[Combo] = []
    for _fitness, combo in sorted(hits, key=lambda p: -p[0]):
        if all(hamming(combo, k) >= d for k in kept):
            kept.append(combo)
    return len(kept)


def diversity_ceiling(hits: Iterable[Hit], distances: Sequence[int] = (2, 3, 4)) -> Dict[int, int]:
    """The same greedy rule on the WHOLE hit set: the most distinct hits anybody can report."""
    materialised = list(hits)
    return {int(d): greedy_diverse(materialised, int(d)) for d in distances}


@dataclass(frozen=True)
class EndpointVerdict:
    """Whether the diversity adjustment carries information on one target at one d."""

    raw: int
    ceiling: int
    informativeness: float
    verdict: str
    threshold: float = INFORMATIVENESS_FLOOR

    @property
    def usable_as_primary(self) -> bool:
        return self.verdict == USABLE


def endpoint_verdict(*, raw: int, ceiling: int, floor: float = INFORMATIVENESS_FLOOR) -> EndpointVerdict:
    """Classify one (raw count, diversity ceiling) pair.

    `informativeness` is the share of the hit set the adjustment removes at the ceiling:
    `1 - ceiling/raw`. Zero means the adjusted endpoint is numerically identical to the
    raw count for every arm - the adjustment is a no-op, not a stricter test.
    """
    if ceiling > raw:
        raise ValueError(f"ceiling {ceiling} exceeds the raw hit count {raw}")
    if raw <= 0:
        return EndpointVerdict(raw=raw, ceiling=ceiling, informativeness=0.0,
                               verdict=DEGENERATE, threshold=floor)
    informativeness = 1.0 - ceiling / raw
    if ceiling == raw:
        verdict = DEGENERATE
    elif informativeness < floor:
        verdict = NEAR_DEGENERATE
    else:
        verdict = USABLE
    return EndpointVerdict(raw=raw, ceiling=ceiling, informativeness=informativeness,
                           verdict=verdict, threshold=floor)
