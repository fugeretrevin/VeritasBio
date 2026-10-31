"""Defect 3 / protocol open item O-11: is the CR6261 H1/H9 genotype set shared?

H2's replication instance needs a within-genotype-set paired contrast: seed s must be
able to draw the SAME 25 initial genotypes on both readouts. eqFP611 blue/red satisfy
that by construction (an identical 8,191-genotype set, row for row). CR6261 H1 and H9
hold 1,886 and 1,841 variants and the sets are NOT identical, so the shared set has to
be measured rather than assumed.
"""
from __future__ import annotations

import pytest

from veritas.calibration.genotype_sets import compare_genotype_sets


def test_compare_counts_the_intersection_and_both_exclusive_remainders():
    cmp = compare_genotype_sets({("A",), ("B",), ("C",)}, {("B",), ("C",), ("D",)})
    assert cmp.n_left == 3
    assert cmp.n_right == 3
    assert cmp.n_shared == 2
    assert cmp.n_left_only == 1
    assert cmp.n_right_only == 1
    assert cmp.jaccard == pytest.approx(0.5)
    assert not cmp.identical


def test_identical_sets_are_reported_as_identical():
    cmp = compare_genotype_sets({("A",), ("B",)}, {("B",), ("A",)})
    assert cmp.identical
    assert cmp.jaccard == 1.0
    assert cmp.n_left_only == cmp.n_right_only == 0


def test_shared_coverage_is_reported_per_side():
    cmp = compare_genotype_sets({("A",), ("B",), ("C",), ("D",)}, {("A",), ("B",)})
    assert cmp.shared_share_left == pytest.approx(0.5)
    assert cmp.shared_share_right == pytest.approx(1.0)


def test_eqfp611_channels_share_their_genotype_set_row_for_row():
    from veritas.calibration.landscapes import load_landscape

    blue, red = load_landscape("eqfp611_blue"), load_landscape("eqfp611_red")
    cmp = compare_genotype_sets(set(blue.combos), set(red.combos))
    assert cmp.identical
    assert cmp.n_shared == 8_191


def test_cr6261_channels_do_not_share_their_genotype_set():
    from veritas.calibration.landscapes import load_landscape

    h1, h9 = load_landscape("cr6261_h1"), load_landscape("cr6261_h9")
    cmp = compare_genotype_sets(set(h1.combos), set(h9.combos))
    assert cmp.n_left == 1_886
    assert cmp.n_right == 1_841
    assert not cmp.identical
    # The quantity O-11 asks for: a non-empty shared set smaller than either side.
    assert 0 < cmp.n_shared < cmp.n_right
    assert cmp.n_left_only + cmp.n_shared == cmp.n_left
    assert cmp.n_right_only + cmp.n_shared == cmp.n_right
