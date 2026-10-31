"""The loader must reproduce the pre-confirmed figures before it is trusted on anything new.

Anchors (docs/preregistration-cross-protein.md §6.1, §6.3; docs/dataset-inventory.md §4):

    TrpB shallow (depth<=2)  168/2,242   = 7.4933%    deep  1,432/157,757 = 0.9077%
    GB1  shallow (depth<=2)   72/2,167   = 3.3226%    deep  1,422/147,193 = 0.9661%

Convention: depth 0 (the wild type) is excluded from both shells; it is not a library
member and does not enter the eligible universe either, which is what makes GB1 149,360
rather than 149,361 and is what every run manifest records.

These tests read the prepared CSVs under scripts/data/. The two largest landscapes
(CreiLOV, CR9114) are exercised by the report rather than here, to keep the suite quick.
"""
from __future__ import annotations

import pytest

from veritas.calibration.landscapes import REGISTRY, load_landscape


@pytest.fixture(scope="module")
def trpb():
    return load_landscape("trpb")


@pytest.fixture(scope="module")
def gb1():
    return load_landscape("gb1")


def test_the_registry_holds_the_nine_study_landscapes():
    assert set(REGISTRY) == {
        "trpb", "gb1", "her2",
        "eqfp611_blue", "eqfp611_red", "creilov",
        "cr9114_h1", "cr6261_h1", "cr6261_h9",
    }


def test_trpb_eligible_universe_and_hit_set(trpb):
    assert trpb.n_universe == 159_999
    assert trpb.n_hits == 1_600
    assert trpb.threshold_source.startswith("manifest")


def test_trpb_shallow_shell_anchor(trpb):
    shell = trpb.shell(depth_cut=2)
    assert (shell.n_hits, shell.n_variants) == (168, 2_242)
    assert shell.rate * 100 == pytest.approx(7.4933, abs=5e-4)


def test_trpb_deep_shell_anchor(trpb):
    deep = trpb.deep_shell(depth_cut=2)
    assert (deep.n_hits, deep.n_variants) == (1_432, 157_757)
    assert deep.rate * 100 == pytest.approx(0.9077, abs=5e-4)


def test_gb1_eligible_universe_excludes_the_wild_type_row(gb1):
    # 149,361 rows in the CSV, one of which is the literal WT.
    assert gb1.n_universe == 149_360
    assert gb1.n_hits == 1_494


def test_gb1_shallow_shell_anchor(gb1):
    shell = gb1.shell(depth_cut=2)
    assert (shell.n_hits, shell.n_variants) == (72, 2_167)
    assert shell.rate * 100 == pytest.approx(3.3226, abs=5e-4)


def test_gb1_deep_shell_anchor(gb1):
    deep = gb1.deep_shell(depth_cut=2)
    assert (deep.n_hits, deep.n_variants) == (1_422, 147_193)
    assert deep.rate * 100 == pytest.approx(0.9661, abs=5e-4)


def test_shells_partition_the_mutant_universe(trpb, gb1):
    for ls in (trpb, gb1):
        shell, deep = ls.shell(depth_cut=2), ls.deep_shell(depth_cut=2)
        assert shell.n_variants + deep.n_variants == ls.n_universe
        assert shell.n_hits + deep.n_hits == ls.n_hits
        assert min(ls.depths) >= 1          # no depth-0 row survives ingestion


def test_diversity_ceilings_match_the_published_inventory(trpb, gb1):
    assert trpb.ceilings((2, 3, 4)) == {2: 166, 3: 20, 4: 4}
    assert gb1.ceilings((2, 3, 4)) == {2: 146, 3: 24, 4: 5}


def test_her2_universe_hit_set_and_degenerate_ceiling():
    her2 = load_landscape("her2")
    assert her2.n_universe == 11_276          # 11,277 train rows less the literal WT
    assert her2.n_hits == 113
    # The defect: the diversity adjustment removes nothing on HER2.
    assert her2.ceilings((2,))[2] == 113


def test_cr6261_universes_and_ceilings():
    h1, h9 = load_landscape("cr6261_h1"), load_landscape("cr6261_h9")
    assert (h1.n_universe, h1.n_hits) == (1_886, 19)
    assert (h9.n_universe, h9.n_hits) == (1_841, 19)
    assert h1.ceilings((2, 3, 4)) == {2: 10, 3: 6, 4: 4}
    assert h9.ceilings((2, 3, 4)) == {2: 13, 3: 4, 4: 2}


def test_eqfp611_channels_share_a_universe_size_and_hit_tier():
    blue, red = load_landscape("eqfp611_blue"), load_landscape("eqfp611_red")
    assert blue.n_universe == red.n_universe == 8_191
    assert blue.n_hits == red.n_hits == 82
    assert blue.ceilings((2, 3))[2] == 39
    assert red.ceilings((2, 3))[2] == 28


def test_computed_percentile_agrees_with_every_recorded_manifest_threshold():
    """Where a manifest exists, the loader's own 99th percentile must reproduce it."""
    for key in ("trpb", "gb1", "her2"):
        ls = load_landscape(key)
        assert ls.threshold == pytest.approx(ls.computed_threshold, rel=1e-12)
