"""A qualification threshold that transfers between landscapes.

`--qualify-threshold` is an absolute fitness and defaults to 0.70. Fitness scales are
not comparable across this corpus: TrpB runs 0 to 1.0, GB1's global best is 8.762,
HER2's is 15.084, and the three bnAb landscapes are -log10(Kd) with an assay floor at
7.0 - so an absolute 0.70 sits *below the floor* there and qualifies every variant
instantly, while on TrpB it is a bar in the extreme upper tail. The frontier policy's
qualification stage therefore gates nothing outside TrpB.

The fix expressed here: state the bar as a QUANTILE of the eligible universe's fitness
distribution and convert it to an absolute per landscape. The conversion uses
`numpy.percentile`, which is the same machinery the benchmark already uses to define
`top1.threshold` (`np.percentile(values, 99.0)`), so the two thresholds are expressed
on one scale and can be read against each other.

Nothing here changes the selector. It computes the numbers a protocol amendment needs.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np


def absolute_for_quantile(values: Sequence[float], quantile: float) -> float:
    """The absolute fitness at `quantile` (0-100) of the eligible universe.

    Linear interpolation, i.e. `np.percentile`'s default and the benchmark's own
    convention for `top1.threshold`.
    """
    if not 0.0 <= quantile <= 100.0:
        raise ValueError(f"quantile must be in [0, 100], got {quantile}")
    return float(np.percentile(np.asarray(values, dtype=np.float64), quantile))


def equivalent_quantile(values: Sequence[float], absolute: float) -> float:
    """The quantile whose absolute equivalent is `absolute` - the inverse of the above.

    This is what answers "what quantile reproduces 0.70's behaviour on TrpB": invert
    the threshold the parameter was tuned on, then read that quantile off every other
    landscape.

    Clamped to [0, 100] for a threshold outside the observed range. Where the
    distribution is flat at `absolute` (a censored mass, e.g. HER2's 51.7% at 0.1 or
    the bnAbs' assay floor at 7.0) the lowest quantile that attains the value is
    returned, which is the conservative reading: it is the first quantile at which the
    bar is met.
    """
    v = np.sort(np.asarray(values, dtype=np.float64))
    n = v.size
    if n == 0:
        raise ValueError("cannot invert a percentile over an empty universe")
    if absolute <= v[0]:
        return 0.0
    if absolute >= v[-1]:
        return 100.0
    # np.percentile(q) interpolates between v[i] and v[i+1] at index q/100*(n-1).
    # Take the FIRST index whose value reaches `absolute`, then interpolate backwards
    # within the preceding step.
    i = int(np.searchsorted(v, absolute, side="left"))
    lo, hi = v[i - 1], v[i]
    frac = 0.0 if hi == lo else (absolute - lo) / (hi - lo)
    return float(100.0 * (i - 1 + frac) / (n - 1))


def qualifying_count(values: Sequence[float], threshold: float) -> int:
    """Variants a run would accept as qualifying: `fitness >= threshold`.

    At-or-above semantics, matching the selector
    (`max(revealed.measured.values()) >= self.qualify_threshold`).
    """
    v = np.asarray(values, dtype=np.float64)
    return int(np.count_nonzero(v >= threshold))


def qualifying_fraction(values: Sequence[float], threshold: float) -> float:
    """Share of the eligible universe at or above `threshold`."""
    v = np.asarray(values, dtype=np.float64)
    if v.size == 0:
        return 0.0
    return float(qualifying_count(v, threshold) / v.size)


def p_qualify_in_draw(*, n_universe: int, n_qualifying: int, n_draws: int) -> float:
    """P(at least one qualifier) in `n_draws` distinct variants drawn uniformly.

    Exact hypergeometric, because the protocol's 25 initial observations are drawn
    WITHOUT replacement; the with-replacement approximation `1-(1-f)**n` is materially
    different on the small libraries (CR6261 holds 1,886 variants).

    This is the diagnostic that exposed Defect 1: on HER2 it is ~1.0 at the absolute
    0.70 default, so every seed logged "qualified before round 1" and the frontier
    policy's qualification stage never gated anything.
    """
    if n_universe <= 0 or n_draws <= 0 or n_qualifying <= 0:
        return 0.0
    if n_qualifying >= n_universe:
        return 1.0
    n_other = n_universe - n_qualifying
    if n_draws > n_other:
        return 1.0
    # P(none) = C(n_other, n_draws) / C(n_universe, n_draws)
    log_p_none = (
        math.lgamma(n_other + 1)
        - math.lgamma(n_other - n_draws + 1)
        - math.lgamma(n_universe + 1)
        + math.lgamma(n_universe - n_draws + 1)
    )
    return float(1.0 - math.exp(log_p_none))


@dataclass(frozen=True)
class QualificationCalibration:
    """What one landscape's qualification bar is, under both definitions."""

    key: str
    n_universe: int
    quantile: float
    absolute_at_quantile: float
    n_qualifying_at_quantile: int
    p_qualify_initial_at_quantile: float
    legacy_absolute: float
    n_qualifying_legacy: int
    p_qualify_initial_legacy: float
    fitness_min: float
    fitness_max: float
    top1_threshold: float

    @property
    def fraction_at_quantile(self) -> float:
        return self.n_qualifying_at_quantile / self.n_universe if self.n_universe else 0.0

    @property
    def fraction_legacy(self) -> float:
        return self.n_qualifying_legacy / self.n_universe if self.n_universe else 0.0

    @property
    def behaviour_changes(self) -> bool:
        """Does the quantile rule accept a different set than the legacy constant?"""
        return self.n_qualifying_at_quantile != self.n_qualifying_legacy


def calibrate(
    key: str,
    values: Sequence[float],
    *,
    quantile: float,
    legacy_absolute: float = 0.70,
    n_initial: int = 25,
    top1_threshold: float,
) -> QualificationCalibration:
    """Both qualification bars on one landscape, with the share and the draw risk."""
    v = np.asarray(values, dtype=np.float64)
    absolute = absolute_for_quantile(v, quantile)
    n_q = qualifying_count(v, absolute)
    n_legacy = qualifying_count(v, legacy_absolute)
    n = int(v.size)
    return QualificationCalibration(
        key=key,
        n_universe=n,
        quantile=quantile,
        absolute_at_quantile=absolute,
        n_qualifying_at_quantile=n_q,
        p_qualify_initial_at_quantile=p_qualify_in_draw(
            n_universe=n, n_qualifying=n_q, n_draws=n_initial
        ),
        legacy_absolute=float(legacy_absolute),
        n_qualifying_legacy=n_legacy,
        p_qualify_initial_legacy=p_qualify_in_draw(
            n_universe=n, n_qualifying=n_legacy, n_draws=n_initial
        ),
        fitness_min=float(v.min()),
        fitness_max=float(v.max()),
        top1_threshold=float(top1_threshold),
    )
