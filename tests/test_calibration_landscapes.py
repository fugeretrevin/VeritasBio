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
    """UPDATED 2026-11-01, and the old literals were pinning a defect.

    This asserted 19 hits on each CR6261 landscape, with d>=2 ceilings of 10 and 13.
    Those numbers required `_load_cached`'s hard-coded `np.percentile(fitness, 99.0)` —
    a 1% bar — while **both landscapes declare `tier_fraction=0.02`** under §4.5's tier
    ladder. The declared tier implies 38 and 37 members, and `tier_size` said so all
    along, so `n_hits` and `tier_size` disagreed on these two landscapes and nowhere else.

    The change was measured before it was made. The blast radius is exactly these two:
    the other seven read their threshold from a run manifest, which rule 2.3 makes the
    authority, and neither CR6261 landscape has ever been scored — 0 of 170 run manifests
    mention them — so no measured threshold was overridden. The figures that move are
    documentation only, and `docs/HANDOFF.md` records the move.

    Kept as literals rather than relaxed, because the point of this test is to notice if
    they move again.
    """
    h1, h9 = load_landscape("cr6261_h1"), load_landscape("cr6261_h9")
    assert (h1.n_universe, h1.n_hits) == (1_886, 38)
    assert (h9.n_universe, h9.n_hits) == (1_841, 37)
    assert h1.ceilings((2, 3, 4)) == {2: 18, 3: 9, 4: 6}
    assert h9.ceilings((2, 3, 4)) == {2: 21, 3: 8, 4: 2}


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


# ======================================================================================
# Manifest threshold reading (rule 2.3: read it, never recompute it)
#
# `_manifest_threshold` knew only the legacy flat schema (`dataset.*`, `top1.threshold`)
# while today's runs are `manifest_schema_version` 2.x and record
# `configuration.dataset.library` and `results.top1.threshold_value`. Every modern run
# was therefore skipped, and on no match the loader RECOMPUTES the 99th percentile -
# which rule 2.3 forbids outright. It went unnoticed because 160 legacy run directories
# survive in the main checkout and satisfied the legacy path, and because the recomputed
# value happens to agree numerically on all nine landscapes today. Latent, not harmless:
# it would have begun failing silently the moment the legacy directories were pruned.
# ======================================================================================
import json  # noqa: E402

from veritas.calibration import landscapes as lsmod  # noqa: E402

_SPEC = REGISTRY["creilov"]


def _write_manifest(directory, name, payload):
    run = directory / name
    run.mkdir(parents=True, exist_ok=True)
    (run / "manifest.json").write_text(json.dumps(payload), encoding="utf-8")
    return run


def _current_schema(threshold, *, key="creilov"):
    """A manifest in the shape runs actually write today."""
    return {
        "manifest_schema_version": "2.1",
        "configuration": {"dataset": {
            "synthetic": False,
            "library": {"key": key, "positions": list(_SPEC.positions),
                        "parent_residues": _SPEC.parent_residues},
        }},
        "results": {"top1": {"threshold_value": threshold, "n_hits": 1655}},
    }


def _legacy_schema(threshold):
    return {
        "dataset": {"synthetic": False,
                    "active_site_positions": list(_SPEC.positions),
                    "parent_residues": _SPEC.parent_residues},
        "top1": {"threshold": threshold},
    }


def test_the_current_manifest_schema_is_read(tmp_path, monkeypatch):
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    _write_manifest(tmp_path, "run_a", _current_schema(4.125))
    assert lsmod._manifest_threshold(_SPEC) == (
        4.125, "manifest (1 run(s), e.g. run_a)")


def test_the_legacy_manifest_schema_is_still_read(tmp_path, monkeypatch):
    """Both shapes, because 160 legacy directories are live evidence."""
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    _write_manifest(tmp_path, "run_old", _legacy_schema(4.125))
    threshold, source = lsmod._manifest_threshold(_SPEC)
    assert (threshold, source.startswith("manifest")) == (4.125, True)


def test_a_matched_run_with_no_readable_threshold_raises(tmp_path, monkeypatch):
    """The decisive test. Rule 2.3 forbids recomputing, so a run that scored this
    landscape and whose threshold cannot be read must NOT fall through to a percentile.
    Falling through is what made the schema drift invisible."""
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    payload = _current_schema(4.125)
    del payload["results"]["top1"]["threshold_value"]
    _write_manifest(tmp_path, "run_broken", payload)
    with pytest.raises(ValueError, match="record no threshold at any known path"):
        lsmod._manifest_threshold(_SPEC)


def test_manifests_that_disagree_still_raise(tmp_path, monkeypatch):
    """Two thresholds for one landscape means the eligible universe moved."""
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    _write_manifest(tmp_path, "run_a", _current_schema(4.125))
    _write_manifest(tmp_path, "run_b", _current_schema(4.130))
    with pytest.raises(ValueError, match="disagree on top1.threshold"):
        lsmod._manifest_threshold(_SPEC)


def test_a_run_for_another_landscape_is_not_mistaken_for_this_one(tmp_path, monkeypatch):
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    other = dict(_current_schema(9.6))
    other["configuration"]["dataset"]["library"] = {
        "key": "trpb", "positions": list(REGISTRY["trpb"].positions),
        "parent_residues": REGISTRY["trpb"].parent_residues}
    _write_manifest(tmp_path, "run_trpb", other)
    assert lsmod._manifest_threshold(_SPEC) is None


def test_a_renamed_key_on_matching_geometry_is_refused(tmp_path, monkeypatch):
    """Matching on geometry AND key, so neither a rename nor a reused key passes."""
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    _write_manifest(tmp_path, "run_renamed", _current_schema(4.125, key="creilov_v2"))
    assert lsmod._manifest_threshold(_SPEC) is None


def test_a_synthetic_run_is_skipped_in_both_schemas(tmp_path, monkeypatch):
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    payload = _current_schema(4.125)
    payload["configuration"]["dataset"]["synthetic"] = True
    _write_manifest(tmp_path, "run_synth", payload)
    assert lsmod._manifest_threshold(_SPEC) is None


def test_the_fallback_reason_does_not_claim_no_manifest_exists(tmp_path, monkeypatch):
    """The old source string was the single literal "computed (no run manifest exists)",
    which was FALSE whenever manifests existed but none were recognised - exactly the
    condition the schema drift created. A provenance field that lies about its own
    provenance is worse than one that is merely unhelpful."""
    monkeypatch.setattr(lsmod, "RUNS_DIR", tmp_path)
    assert lsmod._no_manifest_reason(_SPEC) == "computed (no run manifest exists)"

    other = dict(_current_schema(9.6))
    other["configuration"]["dataset"]["library"] = {
        "key": "trpb", "positions": list(REGISTRY["trpb"].positions),
        "parent_residues": REGISTRY["trpb"].parent_residues}
    _write_manifest(tmp_path, "run_trpb", other)
    reason = lsmod._no_manifest_reason(_SPEC)
    assert "none scored creilov" in reason
    assert "NOT a measured threshold" in reason
    assert "no run manifest exists" not in reason


def test_the_computed_percentile_follows_the_declared_tier():
    """A landscape's bootstrap threshold must describe the tier it declares.

    `_load_cached` used a hard-coded `np.percentile(fitness, 99.0)` for every landscape
    while the `tier_size` property used `tier_fraction`, so any landscape declaring a
    tier other than 1% had a threshold and a tier size describing different tiers. It bit
    exactly two: CR6261-H1 and CR6261-H9 declare `tier_fraction=0.02` under the protocol's
    tier ladder and were handed a 1% bar - 19 hits where their own tier implies 38 and 37.

    Asserted as agreement between the two quantities rather than as the literals 38/37,
    so the test states the invariant and does not have to be edited if a tier changes.
    """
    for key in sorted(REGISTRY):
        ls = load_landscape(key)
        if ls.threshold_source.startswith("manifest"):
            continue  # rule 2.3: the manifest is the authority, not the percentile
        assert ls.n_hits == ls.tier_size, (
            f"{key}: declares a {REGISTRY[key].tier_fraction:.0%} tier implying "
            f"{ls.tier_size} members, but its computed threshold yields {ls.n_hits}")


def test_cr6261_gets_its_registered_two_percent_tier():
    """The two landscapes the hard-coded 99.0 actually broke.

    Pinned as literals because these are the figures the fix moved, and a silent
    regression to 19 would restore a tier the protocol does not declare.
    """
    assert (load_landscape("cr6261_h1").n_hits, load_landscape("cr6261_h9").n_hits) == (38, 37)


def test_reading_the_current_schema_moved_no_hit_set():
    """The safety check for the fix itself.

    Seven of nine landscapes were silently recomputing before it. If the manifest
    threshold differed from the recomputed percentile anywhere, every figure on that
    landscape would move. Measured: it differs nowhere, so the defect was latent and no
    published figure changes. This test is what makes that claim re-checkable rather
    than a one-off assertion in a commit message.
    """
    import numpy as np

    for key in sorted(REGISTRY):
        ls = load_landscape(key)
        recomputed = int(np.count_nonzero(ls.fitness >= ls.computed_threshold))
        assert ls.n_hits == recomputed, (
            f"{key}: manifest threshold yields {ls.n_hits} hits, recomputed percentile "
            f"yields {recomputed} - a figure on this landscape has moved")
