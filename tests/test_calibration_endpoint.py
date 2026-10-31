"""Defect 2: the primary endpoint (diversity-adjusted yield at d>=2) can be degenerate.

On HER2 the d>=2 ceiling equals the raw top-1% ceiling (113 = 113), so the diversity
adjustment carries no information and cannot be the primary endpoint there. These
tests pin the greedy rule, the ceiling, and the precondition verdict.
"""
from __future__ import annotations

import pytest

from veritas.calibration.endpoint import (
    DEGENERATE,
    NEAR_DEGENERATE,
    USABLE,
    diversity_ceiling,
    endpoint_verdict,
    greedy_diverse,
)


def hits(*rows):
    """(fitness, combo) pairs, highest fitness first."""
    return [(f, tuple(c)) for f, c in rows]


def test_greedy_diverse_keeps_the_fittest_and_drops_close_neighbours():
    h = hits((3.0, "AAAA"), (2.0, "AAAB"), (1.0, "BBBB"))
    assert greedy_diverse(h, 2) == 2        # AAAB is 1 away from AAAA, dropped
    assert greedy_diverse(h, 1) == 3        # every distinct genotype is >= 1 apart
    # AAAA and BBBB differ at all four sites, so both survive "differs everywhere".
    assert greedy_diverse(h, 4) == 2


def test_greedy_diverse_walks_in_fitness_order_not_input_order():
    unsorted = hits((1.0, "AAAB"), (3.0, "AAAA"), (2.0, "AABB"))
    # Sorted: AAAA(3) kept, AABB(2) is 2 away kept, AAAB(1) is 1 from AAAA dropped.
    assert greedy_diverse(unsorted, 2) == 2


def test_greedy_diverse_on_an_empty_hit_set():
    assert greedy_diverse([], 2) == 0


def test_diversity_ceiling_reports_every_requested_distance():
    h = hits((3.0, "AAAA"), (2.0, "AAAB"), (1.0, "BBBB"))
    assert diversity_ceiling(h, (2, 3, 4)) == {2: 2, 3: 2, 4: 2}


def test_verdict_is_degenerate_when_the_ceiling_equals_the_raw_count():
    # HER2 at d>=2: 113 of 113.
    v = endpoint_verdict(raw=113, ceiling=113)
    assert v.verdict == DEGENERATE
    assert v.informativeness == pytest.approx(0.0)


def test_verdict_is_near_degenerate_just_below_the_raw_count():
    # HER2 at d>=3: 111 of 113 -- a 1.8% reduction carries no usable signal.
    v = endpoint_verdict(raw=113, ceiling=111)
    assert v.verdict == NEAR_DEGENERATE
    assert v.informativeness == pytest.approx(2 / 113, abs=1e-9)


def test_verdict_is_usable_where_the_adjustment_binds():
    # eqFP611 blue at d>=2: 39 of 82; TrpB: 166 of 1,600.
    assert endpoint_verdict(raw=82, ceiling=39).verdict == USABLE
    assert endpoint_verdict(raw=1600, ceiling=166).verdict == USABLE


def test_verdict_threshold_is_five_percent_and_is_reported():
    assert endpoint_verdict(raw=1000, ceiling=951).verdict == NEAR_DEGENERATE
    assert endpoint_verdict(raw=1000, ceiling=949).verdict == USABLE
    assert endpoint_verdict(raw=1000, ceiling=949).threshold == pytest.approx(0.05)


def test_verdict_on_an_empty_hit_set_is_degenerate_rather_than_a_zero_division():
    assert endpoint_verdict(raw=0, ceiling=0).verdict == DEGENERATE
