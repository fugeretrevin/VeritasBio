"""Genotype canonicalisation: the bug that inverted a conclusion must not come back.

Mutation strings omit wild-type positions, so `V183A:F184A:V227A` is a FOUR-site
genotype whose fourth site is wild type. An implementation that compares the raw
strings, or zips genotypes of unequal length, understates Hamming distances and
inverts diversity conclusions (scripts/analysis/README.md, "Two corrections worth
knowing about"). These tests pin the contract.
"""
from __future__ import annotations

import pytest

from veritas.calibration.genotypes import (
    GenotypeParseError,
    canonical_combo,
    depth,
    hamming,
    min_pairwise_hamming,
    pair_distance_histogram,
    pairs_closer_than,
    parse_mutations,
)

TRPB_POS = (183, 184, 227, 228)
TRPB_WT = "VFVS"


def test_parse_mutations_returns_position_wt_mutant_triples():
    assert parse_mutations("V183N:F184R") == [(183, "V", "N"), (184, "F", "R")]


def test_parse_mutations_treats_literal_wt_as_no_substitutions():
    assert parse_mutations("WT") == []
    assert parse_mutations("wt") == []
    assert parse_mutations("") == []


@pytest.mark.parametrize("bad", ["V", "183", "V1", "x:y", "V183N:Q"])
def test_parse_mutations_rejects_non_substitution_notation(bad):
    assert parse_mutations(bad) is None


def test_canonical_combo_fills_wild_type_at_omitted_positions():
    # Three tokens, FOUR sites: position 228 is silent and therefore wild type.
    combo = canonical_combo("V183A:F184A:V227A", TRPB_POS, TRPB_WT)
    assert combo == ("A", "A", "A", "S")
    assert len(combo) == len(TRPB_POS)


def test_canonical_combo_of_wild_type_is_the_parent_combo():
    assert canonical_combo("WT", TRPB_POS, TRPB_WT) == ("V", "F", "V", "S")


def test_canonical_combo_rejects_a_position_outside_the_library():
    assert canonical_combo("V183A:Q999K", TRPB_POS, TRPB_WT) is None


def test_canonical_combo_rejects_a_token_disagreeing_with_the_parent_residue():
    # Position 183 is V in the parent; a row claiming Q183A describes another design.
    assert canonical_combo("Q183A", TRPB_POS, TRPB_WT) is None


def test_canonical_combo_strict_mode_raises_instead_of_dropping():
    with pytest.raises(GenotypeParseError):
        canonical_combo("V183A:Q999K", TRPB_POS, TRPB_WT, strict=True)


def test_the_truncation_bug_regression():
    """A 3-token and a 4-token string that differ at exactly one site.

    zip() over ('A','A','A') and ('A','A','A','A') compares three sites and returns
    0, which is the bug. Canonical full-length genotypes return 1.
    """
    short = canonical_combo("V183A:F184A:V227A", TRPB_POS, TRPB_WT)
    long = canonical_combo("V183A:F184A:V227A:S228A", TRPB_POS, TRPB_WT)
    assert len(short) == len(long) == 4
    assert hamming(short, long) == 1


def test_hamming_refuses_unequal_lengths():
    with pytest.raises(ValueError):
        hamming(("A", "A", "A"), ("A", "A", "A", "A"))


def test_depth_counts_library_positions_differing_from_the_parent():
    assert depth(canonical_combo("WT", TRPB_POS, TRPB_WT), TRPB_WT) == 0
    assert depth(canonical_combo("V183A", TRPB_POS, TRPB_WT), TRPB_WT) == 1
    assert depth(canonical_combo("V183A:F184A:V227A", TRPB_POS, TRPB_WT), TRPB_WT) == 3


def test_min_pairwise_hamming_over_a_hit_set():
    combos = [("A", "A"), ("A", "B"), ("C", "D")]
    assert min_pairwise_hamming(combos) == 1
    assert min_pairwise_hamming([("A", "A")]) is None


def test_pair_distance_histogram_counts_unordered_pairs():
    combos = [("A", "A"), ("A", "B"), ("C", "D")]
    # AA-AB: 1, AA-CD: 2, AB-CD: 2
    assert pair_distance_histogram(combos) == {1: 1, 2: 2}
    assert pair_distance_histogram([("A", "A")]) == {}


def test_pair_distance_histogram_refuses_ragged_genotypes():
    with pytest.raises(ValueError):
        pair_distance_histogram([("A", "A"), ("A", "A", "A")])


def test_pairs_closer_than_is_what_a_d_separation_filter_can_act_on():
    combos = [("A", "A"), ("A", "B"), ("C", "D")]
    assert pairs_closer_than(combos, 2) == 1      # the one pair at distance 1
    assert pairs_closer_than(combos, 3) == 3      # every pair is closer than 3
    assert pairs_closer_than(combos, 1) == 0      # no duplicates


def test_pair_distance_histogram_blocks_agree_with_the_naive_count():
    # Exercises the blocked path against a direct count on a set big enough to block.
    rng = __import__("random").Random(7)
    combos = [tuple(rng.choice("AB") for _ in range(8)) for _ in range(200)]
    hist = pair_distance_histogram(combos)
    assert sum(hist.values()) == 200 * 199 // 2
    naive = {}
    for i in range(len(combos)):
        for j in range(i + 1, len(combos)):
            d = hamming(combos[i], combos[j])
            naive[d] = naive.get(d, 0) + 1
    assert hist == dict(sorted(naive.items()))
