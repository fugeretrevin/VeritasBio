"""Defect 1: --qualify-threshold is an absolute constant on incomparable scales.

A quantile of the eligible universe transfers; an absolute fitness does not. These
tests pin the quantile <-> absolute round trip, the share of the universe that
qualifies, and the probability that 25 random initial draws qualify immediately --
which is the quantity that showed the bug (every HER2 seed logged "qualified before
round 1").
"""
from __future__ import annotations

import numpy as np
import pytest

from veritas.calibration.qualification import (
    absolute_for_quantile,
    equivalent_quantile,
    p_qualify_in_draw,
    qualifying_count,
    qualifying_fraction,
)


def test_absolute_for_quantile_matches_the_benchmarks_own_percentile_rule():
    # top1.threshold is np.percentile(values, 99.0); the qualification threshold must
    # be defined by the same machinery so the two are comparable.
    values = np.arange(0.0, 1000.0)
    assert absolute_for_quantile(values, 99.0) == pytest.approx(
        float(np.percentile(values, 99.0))
    )


@pytest.mark.parametrize("q", [50.0, 90.0, 99.0, 99.5, 99.9])
def test_equivalent_quantile_round_trips_through_absolute_for_quantile(q):
    rng = np.random.default_rng(0)
    values = rng.normal(size=5000)
    absolute = absolute_for_quantile(values, q)
    assert equivalent_quantile(values, absolute) == pytest.approx(q, abs=1e-6)


def test_equivalent_quantile_is_clamped_outside_the_observed_range():
    values = np.array([0.0, 1.0, 2.0])
    assert equivalent_quantile(values, -5.0) == 0.0
    assert equivalent_quantile(values, 99.0) == 100.0


def test_qualifying_fraction_uses_at_or_above_semantics():
    # The selector qualifies on `max(measured) >= threshold`, so the boundary counts.
    values = np.array([0.0, 1.0, 2.0, 3.0])
    assert qualifying_count(values, 2.0) == 2
    assert qualifying_fraction(values, 2.0) == pytest.approx(0.5)
    assert qualifying_fraction(values, 3.0) == pytest.approx(0.25)
    assert qualifying_fraction(values, 4.0) == 0.0


def test_p_qualify_in_draw_is_exact_hypergeometric_not_binomial():
    # 25 draws without replacement from a universe of 100 holding 10 qualifiers.
    # 1 - C(90,25)/C(100,25) = 0.95211; the with-replacement approximation would give
    # 1 - 0.9**25 = 0.92821, so the two are distinguishable at this size.
    p = p_qualify_in_draw(n_universe=100, n_qualifying=10, n_draws=25)
    assert p == pytest.approx(0.9521134438156095, rel=1e-12)


def test_p_qualify_in_draw_edge_cases():
    assert p_qualify_in_draw(n_universe=100, n_qualifying=0, n_draws=25) == 0.0
    assert p_qualify_in_draw(n_universe=100, n_qualifying=100, n_draws=25) == 1.0
    # More draws than non-qualifiers: qualification is certain.
    assert p_qualify_in_draw(n_universe=30, n_qualifying=10, n_draws=25) == 1.0


def test_a_threshold_above_the_global_best_can_never_qualify():
    values = np.array([0.1, 0.5, 0.9])
    assert qualifying_fraction(values, 1.5) == 0.0
    assert p_qualify_in_draw(n_universe=3, n_qualifying=0, n_draws=3) == 0.0
