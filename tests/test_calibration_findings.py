"""Regression pins for the three measured findings.

Unlike the anchors in `test_calibration_landscapes.py`, these numbers were produced by
this code rather than pre-confirmed elsewhere. They are pinned so that a later change
to the loader, the eligibility rules or the parsing cannot move a figure that a protocol
amendment cites without a test going red.

Reproduce the full tables with:

    .venv\\Scripts\\python.exe -X utf8 -m veritas.calibration.report
"""
from __future__ import annotations

import pytest

from veritas.calibration.endpoint import DEGENERATE, NEAR_DEGENERATE, endpoint_verdict
from veritas.calibration.genotype_sets import compare_genotype_sets
from veritas.calibration.genotypes import min_pairwise_hamming, pairs_closer_than
from veritas.calibration.landscapes import load_landscape
from veritas.calibration.qualification import (
    absolute_for_quantile,
    equivalent_quantile,
    p_qualify_in_draw,
    qualifying_count,
)
from veritas.calibration.report import (
    CALIBRATION_LANDSCAPE,
    LEGACY_QUALIFY_ABSOLUTE,
    N_INITIAL,
    build_report,
)

#: The proposed replacement default, derived below from TrpB.
PROPOSED_QUANTILE = 99.96


# --- Defect 1 -----------------------------------------------------------------------
def test_the_legacy_absolute_is_an_extreme_tail_quantile_on_trpb():
    trpb = load_landscape(CALIBRATION_LANDSCAPE)
    q = equivalent_quantile(trpb.fitness, LEGACY_QUALIFY_ABSOLUTE)
    assert q == pytest.approx(99.9602, abs=5e-4)
    # 64 of 159,999 variants, and the bar sits ABOVE the top-1% threshold.
    assert qualifying_count(trpb.fitness, LEGACY_QUALIFY_ABSOLUTE) == 64
    assert LEGACY_QUALIFY_ABSOLUTE > trpb.threshold


def test_the_proposed_quantile_preserves_trpbs_behaviour():
    trpb = load_landscape(CALIBRATION_LANDSCAPE)
    absolute = absolute_for_quantile(trpb.fitness, PROPOSED_QUANTILE)
    assert absolute == pytest.approx(0.70, abs=5e-4)
    # Same accepted set as the constant it replaces: no behaviour change on TrpB.
    assert qualifying_count(trpb.fitness, absolute) == 64


@pytest.mark.parametrize(
    "key,expected_qualifying",
    [
        ("her2", 1_027),          # ~9% of the library clears an absolute 0.70
        ("cr6261_h1", 1_886),     # the whole library: 0.70 is below the 7.0 assay floor
        ("cr6261_h9", 1_841),
        ("cr9114_h1", 48_840),
        ("creilov", 165_427),
    ],
)
def test_the_legacy_absolute_gates_nothing_outside_trpb(key, expected_qualifying):
    ls = load_landscape(key)
    n_q = qualifying_count(ls.fitness, LEGACY_QUALIFY_ABSOLUTE)
    assert n_q == expected_qualifying
    # Which is why every seed logged "qualified before round 1".
    p = p_qualify_in_draw(n_universe=ls.n_universe, n_qualifying=n_q, n_draws=N_INITIAL)
    assert p > 0.9


def test_the_quantile_rule_equalises_the_qualification_bar_across_landscapes():
    """The property an absolute threshold does not have.

    Under the quantile rule, P(a random 25-variant start already qualifies) is within a
    factor of two across all nine landscapes; under the absolute 0.70 it ranges from
    0.01 on TrpB to 1.00 on four of them.
    """
    keys = ["trpb", "gb1", "her2", "eqfp611_blue", "eqfp611_red",
            "creilov", "cr9114_h1", "cr6261_h1", "cr6261_h9"]
    ps = []
    for key in keys:
        ls = load_landscape(key)
        absolute = absolute_for_quantile(ls.fitness, PROPOSED_QUANTILE)
        n_q = qualifying_count(ls.fitness, absolute)
        ps.append(p_qualify_in_draw(
            n_universe=ls.n_universe, n_qualifying=n_q, n_draws=N_INITIAL))
    assert max(ps) / min(ps) < 2.0
    assert max(ps) < 0.02


# --- Defect 2 -----------------------------------------------------------------------
def test_her2_primary_endpoint_is_degenerate_and_why():
    her2 = load_landscape("her2")
    ceilings = her2.ceilings((2, 3, 4))
    assert (her2.n_hits, ceilings[2], ceilings[3], ceilings[4]) == (113, 113, 111, 94)
    assert endpoint_verdict(raw=her2.n_hits, ceiling=ceilings[2]).verdict == DEGENERATE
    assert endpoint_verdict(raw=her2.n_hits, ceiling=ceilings[3]).verdict == NEAR_DEGENERATE
    # The mechanism: nothing in the hit set is closer than 2, and only 2 of the
    # C(113,2) = 6,328 pairs sit at that distance.
    combos = her2.hit_combos
    assert min_pairwise_hamming(combos) == 2
    assert pairs_closer_than(combos, 2) == 0
    assert pairs_closer_than(combos, 3) == 2


@pytest.mark.parametrize(
    "key,raw,d2,d3",
    [
        ("trpb", 1_600, 166, 20),
        ("gb1", 1_494, 146, 24),
        ("eqfp611_blue", 82, 39, 14),
        ("eqfp611_red", 82, 28, 10),
        ("cr6261_h1", 19, 10, 6),
        ("cr6261_h9", 19, 13, 4),
    ],
)
def test_the_adjustment_binds_everywhere_except_her2(key, raw, d2, d3):
    ls = load_landscape(key)
    ceilings = ls.ceilings((2, 3))
    assert (ls.n_hits, ceilings[2], ceilings[3]) == (raw, d2, d3)
    assert endpoint_verdict(raw=raw, ceiling=d2).usable_as_primary


# --- Defect 3 -----------------------------------------------------------------------
def test_the_cr6261_intersection_size():
    h1, h9 = load_landscape("cr6261_h1"), load_landscape("cr6261_h9")
    cmp = compare_genotype_sets(set(h1.combos), set(h9.combos))
    assert (cmp.n_left, cmp.n_right) == (1_886, 1_841)
    assert cmp.n_shared == 1_811
    assert (cmp.n_left_only, cmp.n_right_only) == (75, 30)
    assert cmp.shared_share_left == pytest.approx(0.9602, abs=5e-4)
    assert cmp.shared_share_right == pytest.approx(0.9837, abs=5e-4)


def test_the_cr6261_hit_sets_do_not_overlap_at_all():
    """Two antigens, one library: 17 of H1's 19 hits and all 19 of H9's lie in the
    shared set, and not one genotype is a hit on both."""
    h1, h9 = load_landscape("cr6261_h1"), load_landscape("cr6261_h9")
    shared = set(h1.combos) & set(h9.combos)
    h1_hits = {c for c in h1.hit_combos if c in shared}
    h9_hits = {c for c in h9.hit_combos if c in shared}
    assert (len(h1_hits), len(h9_hits)) == (17, 19)
    assert h1_hits & h9_hits == set()


# --- the report itself ---------------------------------------------------------------
def test_the_report_runs_on_a_subset_and_reproduces_the_anchors():
    text, payload = build_report(("her2", "cr6261_h1", "cr6261_h9"))
    assert payload["anchors_reproduced"] is True
    assert "ANCHORS REPRODUCED" in text
    assert payload["endpoint"]["degenerate_d2"] == ["her2"]
    assert payload["genotype_sets"]["cr6261"]["n_shared"] == 1_811
    # A subset must not silently report landscapes it was not asked for.
    assert {r["landscape"] for r in payload["endpoint"]["rows"]} == {
        "trpb", "her2", "cr6261_h1", "cr6261_h9"}
