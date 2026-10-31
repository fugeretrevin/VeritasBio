"""The pre-flight gate: six assertions over artefacts we already write.

Every check here exists because a defect was recorded in one of our own artefacts and
nothing asserted on it. So every test comes in a pair: the check must FAIL on an input
constructed to carry the defect, and PASS on the same input without it. A check that
cannot catch the defect it was written for is worthless, so the bad inputs reuse the
real figures - ALDE's 26 of 125, HER2's 113 == 113, CR6261's 198 of 220 constant
one-hot columns.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from veritas import preflight as pf


# ======================================================================================
# The result type
# ======================================================================================
def test_a_check_result_refuses_to_be_read_as_a_bare_boolean():
    """`if result:` is the bug this gate exists to prevent, so it raises."""
    result = pf.CheckResult(
        check="x", status=pf.FAIL, measured=26, expected=125, summary="short"
    )
    with pytest.raises(TypeError, match="bare boolean"):
        bool(result)


def test_a_failure_states_the_measured_value_against_the_expectation():
    result = pf.CheckResult(
        check="budget", status=pf.FAIL, measured=26, expected=125, summary="short budget"
    )
    rendered = result.render()
    assert "26" in rendered and "125" in rendered
    assert result.failed is True
    assert result.to_dict()["measured"] == 26


def test_not_applicable_must_carry_a_reason():
    with pytest.raises(ValueError, match="reason"):
        pf.CheckResult(check="x", status=pf.NOT_APPLICABLE, measured=None,
                       expected=None, summary="n/a")
    ok = pf.CheckResult(check="x", status=pf.NOT_APPLICABLE, measured=None, expected=None,
                        summary="n/a", reason="no run directory supplied")
    assert ok.failed is False
    assert ok.applicable is False


def test_a_report_exits_non_zero_when_any_check_failed():
    good = pf.CheckResult(check="a", status=pf.PASS, measured=1, expected=1, summary="ok")
    bad = pf.CheckResult(check="b", status=pf.FAIL, measured=0, expected=1, summary="no")
    assert pf.PreflightReport(library="x", arms=("a",), results=(good,)).exit_code == 0
    assert pf.PreflightReport(library="x", arms=("a",), results=(good, bad)).exit_code == 1


def test_a_report_renders_every_check_and_names_the_failures():
    good = pf.CheckResult(check="a", status=pf.PASS, measured=1, expected=1, summary="ok")
    bad = pf.CheckResult(check="b", status=pf.FAIL, measured=0, expected=1, summary="no")
    text = pf.PreflightReport(library="her2", arms=("alde",), results=(good, bad)).render()
    assert "her2" in text
    assert "a" in text and "b" in text
    assert pf.FAIL in text


# ======================================================================================
# 1. Corpus integrity
# ======================================================================================
def test_the_corpus_lock_check_delegates_to_the_freeze_script():
    """The status mirrors `freeze_cross_protein_corpus.py verify`'s own exit code.

    Deliberately NOT asserted here: that the corpus on disk verifies right now. That is
    the gate's answer, not a property of this code, and it moves whenever the pin is
    being updated.
    """
    lock = pf.result_for(pf.check_corpus_integrity(), "corpus.lock")
    assert lock.detail["delegated_to"].endswith("verify")
    assert lock.detail["verify_output"]
    assert lock.status == (pf.PASS if lock.detail["verify_exit_code"] == 0 else pf.FAIL)
    assert lock.measured == lock.detail["verify_exit_code"]


def test_the_corpus_lock_check_fails_when_the_corpus_is_absent(tmp_path):
    results = pf.check_corpus_integrity(corpus_dir=tmp_path / "nothing-here")
    lock = pf.result_for(results, "corpus.lock")
    assert lock.status == pf.FAIL
    assert lock.detail["verify_exit_code"] != 0
    # The failure names the pinned paths it could not find, rather than only "failed".
    assert "Data/" in lock.detail["verify_output"]


def test_a_clean_table_read_while_a_raw_table_exists_is_flagged():
    """CombinGym's clean CR9114 workbook dropped a whole quadrant. See docs/cr9114-reprep.md."""
    results = pf.check_corpus_integrity(
        landscape_files={"cr6261_h1": "Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx"},
    )
    derived = pf.result_for(results, "corpus.derived_table")
    assert derived.status == pf.FAIL
    assert derived.measured == 1
    flagged = derived.detail["flagged"][0]
    assert flagged["landscape"] == "cr6261_h1"
    assert "Raw" in flagged["raw_alternative"]


def test_a_landscape_reading_the_raw_table_is_not_flagged():
    results = pf.check_corpus_integrity(
        landscape_files={"cr9114_h1": "Data/DMS/Raw/bnAbs_CR9114_raw.csv"},
    )
    derived = pf.result_for(results, "corpus.derived_table")
    assert derived.status == pf.PASS
    assert derived.measured == 0


def test_an_acknowledged_derived_table_is_recorded_rather_than_silently_ignored():
    results = pf.check_corpus_integrity(
        landscape_files={"cr6261_h1": "Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx"},
        acknowledged_derived={"cr6261_h1": "raw table has no H1/H9 split; reviewed 2026-10-31"},
    )
    derived = pf.result_for(results, "corpus.derived_table")
    assert derived.status == pf.PASS
    assert derived.detail["acknowledged"]["cr6261_h1"].startswith("raw table has no")


def test_the_default_scope_is_the_landscapes_this_study_reads_and_cr9114_stays_on_raw():
    """If cr9114_h1 is ever repointed at the clean workbook, this fails."""
    derived = pf.result_for(pf.check_corpus_integrity(), "corpus.derived_table")
    examined = set(derived.detail["landscapes_examined"])
    assert "cr9114_h1" in examined
    assert "hiv1_protease" not in examined  # a reserve landscape nothing scores against
    assert "cr9114_h1" not in {f["landscape"] for f in derived.detail["flagged"]}


def test_the_prepared_table_digest_pin_catches_a_changed_prepared_table():
    pinned = pf.PREPARED_TABLE_DIGESTS["cr9114_h1"]
    good = pf.result_for(
        pf.check_corpus_integrity(prepared_digests={"cr9114_h1": pinned}), "corpus.prepared")
    assert good.status == pf.PASS

    bad = pf.result_for(
        pf.check_corpus_integrity(prepared_digests={"cr9114_h1": "0" * 64}), "corpus.prepared")
    assert bad.status == pf.FAIL
    assert pinned in str(bad.detail["mismatched"])


def test_the_prepared_table_check_records_the_lock_gap_instead_of_passing_silently():
    results = pf.check_corpus_integrity(prepared_digests={})
    prepared = pf.result_for(results, "corpus.prepared")
    assert prepared.status == pf.NOT_APPLICABLE
    assert "lock" in prepared.reason


# ======================================================================================
# 2. Encoding compatibility
# ======================================================================================
def test_constant_columns_are_detected_and_counted():
    varying = np.array([[0.0, 1.0], [1.0, 0.0]])
    assert pf.constant_columns(varying).size == 0
    constant = np.array([[0.0, 0.0], [1.0, 0.0]])
    assert pf.constant_columns(constant).tolist() == [1]


def test_the_encoding_check_fails_at_her2_scale_and_reports_the_count():
    """HER2: 55 of 200 one-hot columns are constant - a residue never observed there."""
    rng = np.random.default_rng(0)
    X = rng.integers(0, 2, size=(64, 200)).astype(float)
    X[:, :55] = 0.0
    result = pf.check_encoding_compatibility(arm="alde", X=X, library="her2")
    assert result.status == pf.FAIL
    assert result.measured == 55
    assert result.expected == 0
    assert "55" in result.render() and "200" in result.render()


def test_the_encoding_check_passes_on_a_saturated_library():
    """TrpB and GB1 have 0 constant columns: every amino acid appears at every site."""
    X = np.eye(80)[np.arange(80) % 80]
    result = pf.check_encoding_compatibility(arm="alde", X=X, library="trpb")
    assert result.status == pf.PASS
    assert result.measured == 0


def test_availability_is_not_compatibility():
    """The existing probe reports ALDE available on all nine landscapes. It is not."""
    X = np.zeros((8, 4))
    result = pf.check_encoding_compatibility(arm="alde", X=X, library="cr9114_h1",
                                             available=True)
    assert result.status == pf.FAIL
    assert result.detail["available"] is True
    assert result.detail["compatible"] is False


def test_the_encoding_check_is_not_applicable_to_an_arm_without_an_encoded_design_space():
    X = np.zeros((8, 4))
    result = pf.check_encoding_compatibility(arm="ridge", X=X, library="her2")
    assert result.status == pf.NOT_APPLICABLE
    assert "encoded design space" in result.reason


def test_one_hot_refuses_ragged_genotypes_rather_than_zipping_them():
    with pytest.raises(ValueError, match="unequal length"):
        pf.one_hot([("A", "C"), ("A",)])


def test_one_hot_allocates_twenty_columns_per_site_like_aldes_encoder():
    X = pf.one_hot([("A", "C"), ("C", "A")])
    assert X.shape == (2, 2 * len(pf.ALL_AAS))
    assert X.sum(axis=1).tolist() == [2.0, 2.0]


@pytest.mark.parametrize("key,expected", [
    ("cr6261_h1", (198, 220)),
    ("cr6261_h9", (198, 220)),
    ("eqfp611_blue", (234, 260)),
])
def test_the_computed_constant_column_counts_reproduce_the_recorded_probe(key, expected):
    """docs/alde-compatibility.md §1, measured against ALDE's own code. If this check
    cannot reproduce those counts it cannot catch the defect."""
    X = pf.encode_landscape(key)
    assert (pf.constant_columns(X).size, X.shape[1]) == expected
    assert pf.REFERENCE_CONSTANT_COLUMNS[key] == expected


# ======================================================================================
# 3. Budget
# ======================================================================================
def _obs(**kwargs):
    base = dict(method="alde", seed=0, n_unique_measurements=125)
    base.update(kwargs)
    return pf.BudgetObservation(**base)


def test_the_budget_check_passes_when_distinct_measurements_equal_the_budget():
    results = pf.check_budget([_obs()], declared_budget=125)
    assert pf.result_for(results, "budget.alde").status == pf.PASS


def test_her2_alde_26_of_125_is_degenerate():
    """26 = the 25 shared starting variants plus one variant measured 100 times."""
    results = pf.check_budget(
        [_obs(n_unique_measurements=26, rejected_duplicate=99)],
        declared_budget=125, initial_size=25,
    )
    result = pf.result_for(results, "budget.alde")
    assert result.status == pf.FAIL
    assert result.measured == 26
    assert result.expected == 125
    assert result.detail["classification"] == pf.BUDGET_DEGENERATE
    assert "unrecoverable" in result.summary


def test_gb1_reliability_120_of_125_is_recoverable_not_degenerate():
    """Our own arms measured 120-124 where baselines measured 125: ineligible proposals."""
    results = pf.check_budget(
        [_obs(method="panoptes_default", n_unique_measurements=120, rejected_ineligible=5)],
        declared_budget=125,
    )
    result = pf.result_for(results, "budget.panoptes_default")
    assert result.status == pf.FAIL
    assert result.detail["classification"] == pf.BUDGET_RECOVERABLE
    assert "unrecoverable" not in result.summary


def test_trpb_frontier_116_of_125_is_recoverable():
    results = pf.check_budget(
        [_obs(method="panoptes_qualify_then_optimize", n_unique_measurements=116,
              rejected_ineligible=9)],
        declared_budget=125,
    )
    result = pf.result_for(results, "budget.panoptes_qualify_then_optimize")
    assert result.detail["classification"] == pf.BUDGET_RECOVERABLE


def test_degenerate_and_recoverable_are_not_collapsed_into_one_failure():
    results = pf.check_budget(
        [_obs(n_unique_measurements=26, rejected_duplicate=99),
         _obs(method="panoptes_default", n_unique_measurements=120, rejected_ineligible=5)],
        declared_budget=125,
    )
    classifications = {r.check: r.detail["classification"] for r in results}
    assert classifications["budget.alde"] == pf.BUDGET_DEGENERATE
    assert classifications["budget.panoptes_default"] == pf.BUDGET_RECOVERABLE


def test_an_arm_that_offered_fewer_candidates_than_requested_is_recoverable():
    """The frontier arm's shortfall has no rejections at all: the eligible pool ran out.

    GB1 reliability (mean 123.23, min 120) and TrpB frontier (mean 120.23, min 116) are
    this case, and `round_metrics.csv` already records requested/accepted/rejected.
    """
    results = pf.check_budget(
        [_obs(method="panoptes_qualify_then_optimize", n_unique_measurements=121,
              requested=125, accepted=121, rejected=0)],
        declared_budget=125,
    )
    result = pf.result_for(results, "budget.panoptes_qualify_then_optimize")
    assert result.status == pf.FAIL
    assert result.detail["classification"] == pf.BUDGET_RECOVERABLE
    assert "0 rejections" in result.detail["reasons"][0]


def test_an_arm_whose_rejections_collapse_onto_one_candidate_stays_degenerate():
    """Same distinct count, opposite cause: duplicates, not an exhausted pool."""
    results = pf.check_budget(
        [_obs(n_unique_measurements=121, requested=125, accepted=96, rejected=29,
              rejected_duplicate=29)],
        declared_budget=125,
    )
    assert pf.result_for(results, "budget.alde").detail["classification"] == (
        pf.BUDGET_DEGENERATE)


def test_the_newer_round_metrics_columns_classify_exactly_rather_than_by_inference(tmp_path):
    """`round_metrics.csv` now records `backfill_added`, `distinct_proposed` and
    `degenerate_proposals`. Those are the run's own determination, so the gate reuses
    them instead of inferring from rejection counts."""
    path = tmp_path / "round_metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["provenance", "method", "seed", "round", "n_unique_measurements",
                    "requested", "accepted", "rejected", "status",
                    "backfill_added", "distinct_proposed", "degenerate_proposals"])
        # A collapsed surrogate: the run itself flagged the proposals as degenerate.
        w.writerow(["s", "alde", "0", "1", "40", "125", "40", "85", "short_batch",
                    "0", "1", "1"])
        # Ineligible proposals, filled from elsewhere: the arm searched.
        w.writerow(["s", "frontier", "0", "1", "121", "125", "121", "4", "completed",
                    "4", "121", "0"])
    observations = {o.method: o for o in pf.budget_observations_from_run(tmp_path)}
    assert observations["alde"].degenerate_proposals == 1
    assert observations["frontier"].backfill_added == 4
    assert observations["frontier"].distinct_proposed == 121

    results = pf.check_budget(observations.values(), declared_budget=125)
    assert pf.result_for(results, "budget.alde").detail["classification"] == (
        pf.BUDGET_DEGENERATE)
    assert pf.result_for(results, "budget.frontier").detail["classification"] == (
        pf.BUDGET_RECOVERABLE)


def test_a_shortfall_with_no_reason_breakdown_is_unclassified_and_still_fails():
    results = pf.check_budget([_obs(n_unique_measurements=100)], declared_budget=125)
    result = pf.result_for(results, "budget.alde")
    assert result.status == pf.FAIL
    assert result.detail["classification"] == pf.BUDGET_UNCLASSIFIED
    assert "selection_trace.csv" in result.detail["missing_evidence"]


def test_budget_observations_are_read_from_round_metrics_and_the_selection_trace(tmp_path):
    rounds = tmp_path / "round_metrics.csv"
    with rounds.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["provenance", "method", "seed", "round", "n_unique_measurements",
                    "requested", "accepted", "rejected", "status"])
        w.writerow(["s", "alde", "0", "0", "25", "25", "25", "0", "completed"])
        w.writerow(["s", "alde", "0", "1", "26", "25", "1", "24", "short_batch"])
        w.writerow(["s", "frontier", "0", "0", "25", "25", "25", "0", "completed"])
        w.writerow(["s", "frontier", "0", "1", "46", "25", "21", "0", "short_batch"])
    trace = tmp_path / "selection_trace.csv"
    with trace.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["provenance", "method", "seed", "round", "rank", "variant_id",
                    "status", "reason", "revealed_fitness"])
        for _ in range(24):
            w.writerow(["s", "alde", "0", "1", "1", "A106T", "rejected",
                        "duplicate_within_batch", ""])

    observations = pf.budget_observations_from_run(tmp_path)
    assert len(observations) == 2
    alde = next(o for o in observations if o.method == "alde")
    assert alde.n_unique_measurements == 26
    assert alde.rejected_duplicate == 24
    assert alde.rejected_ineligible == 0
    assert (alde.requested, alde.accepted, alde.rejected) == (50, 26, 24)

    results = pf.check_budget(observations, declared_budget=125)
    assert pf.result_for(results, "budget.alde").detail["classification"] == (
        pf.BUDGET_DEGENERATE)
    # Same artefacts, the other arm, the other classification: nothing is collapsed.
    assert pf.result_for(results, "budget.frontier").detail["classification"] == (
        pf.BUDGET_RECOVERABLE)


# ======================================================================================
# 4. Endpoint viability
# ======================================================================================
def test_her2_is_degenerate_at_the_intended_distance():
    """113 == 113: the diversity adjustment removes nothing, so it is not a stricter test."""
    result = pf.check_endpoint_viability(library="her2", raw=113, ceiling=113,
                                         distance=2, budget=125)
    assert result.status == pf.FAIL
    assert result.detail["verdict"] == "DEGENERATE"
    assert result.measured["ceiling"] == 113


def test_trpb_is_usable_at_the_intended_distance():
    result = pf.check_endpoint_viability(library="trpb", raw=1600, ceiling=166,
                                         distance=2, budget=125)
    assert result.status == pf.PASS
    assert result.detail["verdict"] == "USABLE"


def test_budget_saturation_fails_even_where_the_adjustment_binds():
    """eqFP611: the d>=2 adjustment binds (39 of 82) but 125 measurements >= an 82-hit tier."""
    result = pf.check_endpoint_viability(library="eqfp611_blue", raw=82, ceiling=39,
                                         distance=2, budget=125)
    assert result.status == pf.FAIL
    assert result.detail["verdict"] == "USABLE"
    assert result.detail["budget_saturates_tier"] is True


def test_cr6261_nineteen_hit_tier_is_saturated_by_the_budget():
    result = pf.check_endpoint_viability(library="cr6261_h1", raw=19, ceiling=10,
                                         distance=2, budget=125)
    assert result.status == pf.FAIL
    assert result.detail["budget_saturates_tier"] is True


def test_the_endpoint_check_reuses_the_calibration_verdict_including_near_degenerate():
    result = pf.check_endpoint_viability(library="her2", raw=113, ceiling=111,
                                         distance=3, budget=50)
    assert result.detail["verdict"] == "NEAR-DEGENERATE"
    assert result.status == pf.FAIL


def test_an_unsaturated_budget_under_a_binding_adjustment_passes():
    result = pf.check_endpoint_viability(library="creilov", raw=1655, ceiling=400,
                                         distance=2, budget=125)
    assert result.status == pf.PASS
    assert result.detail["budget_saturates_tier"] is False


# ======================================================================================
# 5. Threshold transfer
# ======================================================================================
TRPB_LIKE = np.linspace(0.0, 1.0, 10_001)
BNAB_LIKE = np.linspace(7.0, 10.15, 2_001)


def _manifest(values, quantile, *, source="eligible universe fitness", absolute=None,
              top1=None):
    from veritas.calibration.qualification import absolute_for_quantile
    return {
        "qualification": {
            "quantile": quantile,
            "absolute": absolute_for_quantile(values, quantile) if absolute is None else absolute,
            "source": source,
        },
        "top1": {"threshold": float(np.percentile(values, 99.0))},
    }


def test_a_quantile_bar_resolved_from_the_same_values_array_passes():
    manifest = _manifest(TRPB_LIKE, 99.0)
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE, manifest=manifest)
    assert result.status == pf.PASS, result.render()


def test_an_absolute_constant_bar_fails_because_it_does_not_transfer():
    """0.70 is TrpB's 99.96th percentile and below the global minimum on four landscapes."""
    manifest = {"qualification": {"absolute": 0.70}, "top1": {"threshold": 0.99}}
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE, manifest=manifest)
    assert result.status == pf.FAIL
    assert "quantile" in result.detail["violations"][0]


def test_a_bar_below_the_global_minimum_gates_nothing_and_fails():
    manifest = {
        "qualification": {"quantile": 0.0, "absolute": 0.70, "source": "legacy constant"},
        "top1": {"threshold": float(np.percentile(BNAB_LIKE, 99.0))},
    }
    result = pf.check_threshold_transfer(library="cr6261_h1", values=BNAB_LIKE,
                                         manifest=manifest)
    assert result.status == pf.FAIL
    assert any("below the global minimum" in v for v in result.detail["violations"])
    assert result.detail["p_qualify_in_initial_draw"] == pytest.approx(1.0)


def test_a_bar_the_initial_draw_already_clears_gates_nothing():
    """HER2 at the real figures: 1,027 of 11,276 variants clear the absolute 0.70, so a
    25-variant draw already qualifies with probability 0.908. Every HER2 seed logged
    "qualified before round 1" and the qualification stage never ran. A `p == 1.0` test
    would not have caught this, which is why GATING_P_MAX is 0.5 and not 1.0."""
    her2_like = np.concatenate([np.full(11_276 - 1_027, 0.1),
                                np.linspace(0.8, 15.08, 1_027)])
    manifest = {
        "qualification": {"quantile": 0.0, "absolute": 0.70, "source": "legacy constant"},
        "top1": {"threshold": float(np.percentile(her2_like, 99.0))},
    }
    result = pf.check_threshold_transfer(library="her2", values=her2_like,
                                         manifest=manifest)
    assert result.status == pf.FAIL
    assert result.detail["p_qualify_in_initial_draw"] == pytest.approx(0.908, abs=0.01)
    assert any("gates nothing" in v for v in result.detail["violations"])


def test_a_bar_that_binds_on_most_seeds_is_not_flagged_as_ungating():
    """The companion to the test above: a bar the starting sample usually does NOT clear
    is a working gate, and must not be reported as one that gates nothing."""
    manifest = _manifest(TRPB_LIKE, 99.0)
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE,
                                         manifest=manifest)
    assert result.status == pf.PASS
    assert result.detail["p_qualify_in_initial_draw"] < pf.GATING_P_MAX


def test_a_bar_resolved_from_a_different_values_array_fails():
    manifest = _manifest(TRPB_LIKE, 99.0)
    # Same quantile, different universe: the recorded absolute no longer resolves.
    result = pf.check_threshold_transfer(library="trpb", values=BNAB_LIKE, manifest=manifest)
    assert result.status == pf.FAIL
    assert any("does not resolve" in v for v in result.detail["violations"])


def test_the_manifest_must_record_the_quantile_the_absolute_and_the_source():
    manifest = _manifest(TRPB_LIKE, 99.0)
    del manifest["qualification"]["source"]
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE, manifest=manifest)
    assert result.status == pf.FAIL
    assert any("source" in v for v in result.detail["violations"])


def test_the_bar_is_read_from_the_shape_real_manifests_actually_carry():
    """A real run records `configuration.components.qualify_threshold: 0.7` and
    `results.top1.threshold_value`, and nothing else - no quantile, no source."""
    manifest = {
        "configuration": {"components": {"qualify_threshold": 0.7}},
        "results": {"top1": {"threshold_value": float(np.percentile(TRPB_LIKE, 99.0))}},
    }
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE,
                                         manifest=manifest)
    assert result.status == pf.FAIL
    assert result.detail["read_from"] == "configuration.components.qualify_threshold"
    assert result.detail["top1_read_from"] == "results.top1.threshold_value"
    assert result.measured["absolute"] == 0.7
    assert "quantile" in result.detail["violations"][0]


def test_a_manifest_with_no_qualification_bar_names_every_path_it_checked():
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE, manifest={})
    assert result.status == pf.FAIL
    assert "configuration.components.qualify_threshold" in result.detail["violations"][0]
    assert result.detail["manifest_paths_checked"]


def test_the_bar_must_come_from_the_array_that_defines_top1_threshold():
    manifest = _manifest(TRPB_LIKE, 99.0)
    manifest["top1"]["threshold"] = 12345.0
    result = pf.check_threshold_transfer(library="trpb", values=TRPB_LIKE, manifest=manifest)
    assert result.status == pf.FAIL
    assert any("top1.threshold" in v for v in result.detail["violations"])


# ======================================================================================
# 6. Row-order invariance
# ======================================================================================
ROWS = [{"id": f"v{i}", "fitness": 1.0 - i / 100.0} for i in range(100)]


def _invariant_runner(rows, seed):
    """Any correct method is invariant to input row order."""
    return sorted(r["id"] for r in sorted(rows, key=lambda r: -r["fitness"])[:5])


def _index_zero_runner(rows, seed):
    """ALDE with an all-NaN surrogate: np.argmax over all-NaN returns 0, forever."""
    return [rows[0]["id"]] * 5


def test_row_order_invariance_passes_for_an_order_invariant_runner():
    result = pf.check_row_order_invariance(rows=ROWS, runner=_invariant_runner, seed=0)
    assert result.status == pf.PASS
    assert result.detail["shuffle_seed"] == pf.DEFAULT_SHUFFLE_SEED


def test_row_order_invariance_fails_for_a_runner_that_returns_index_zero_forever():
    result = pf.check_row_order_invariance(rows=ROWS, runner=_index_zero_runner, seed=0)
    assert result.status == pf.FAIL
    assert result.detail["baseline_fingerprint"] != result.detail["shuffled_fingerprint"]


def test_a_failure_records_the_shuffle_seed_so_it_reproduces():
    result = pf.check_row_order_invariance(rows=ROWS, runner=_index_zero_runner, seed=3,
                                           shuffle_seed=12345)
    assert result.detail["shuffle_seed"] == 12345
    assert result.detail["run_seed"] == 3
    assert "12345" in result.render()
    # Same seed, same permutation: the failure is reproducible.
    assert (pf.shuffled_order(len(ROWS), 12345) == pf.shuffled_order(len(ROWS), 12345))
    assert pf.shuffled_order(len(ROWS), 12345) != pf.shuffled_order(len(ROWS), 999)


def test_row_order_invariance_without_a_runner_is_not_applicable_and_says_why():
    result = pf.check_row_order_invariance(rows=ROWS)
    assert result.status == pf.NOT_APPLICABLE
    assert "run" in result.reason


def test_the_row_zero_hazard_is_reported_for_a_table_sorted_by_fitness_descending():
    """CR9114's row 0 IS the global optimum, so a wholly broken arm posts a perfect result."""
    fitness = np.array([r["fitness"] for r in ROWS])
    result = pf.check_row_order_invariance(rows=ROWS, fitness=fitness)
    hazard = result.detail["row_zero_hazard"]
    assert hazard["row_zero_is_global_optimum"] is True
    assert hazard["percentile"] == pytest.approx(100.0)

    shuffled_fitness = np.array([0.5, 1.0, 0.2])
    other = pf.check_row_order_invariance(rows=ROWS[:3], fitness=shuffled_fitness)
    assert other.detail["row_zero_hazard"]["row_zero_is_global_optimum"] is False


def order_invariant_combo_runner(rows, seed):
    """A runner over a landscape's canonical combos that does not depend on row order."""
    return {"n": len(rows), "distinct": len(set(rows))}


def order_dependent_combo_runner(rows, seed):
    """The ALDE failure mode in miniature: whatever happens to be at index 0."""
    return {"picked": [rows[0]] * 5}


def test_a_single_row_table_cannot_be_shuffled_and_is_not_applicable():
    result = pf.check_row_order_invariance(rows=ROWS[:1], runner=_invariant_runner)
    assert result.status == pf.NOT_APPLICABLE


# ======================================================================================
# The entry point and the CLI
# ======================================================================================
def test_run_preflight_returns_a_result_for_every_check():
    report = pf.run_preflight(library="cr6261_h1", arms=("alde",),
                              only={pf.CHECK_ENCODING})
    assert {r.check.split(".")[0] for r in report.results} == set(pf.ALL_CHECKS)
    assert report.exit_code == 1  # 198 of 220 constant columns


def test_a_skipped_check_is_recorded_rather_than_silent():
    report = pf.run_preflight(library="cr6261_h1", arms=("alde",),
                              only={pf.CHECK_ENCODING})
    skipped = [r for r in report.results if r.status == pf.SKIPPED]
    assert {r.check for r in skipped} == set(pf.ALL_CHECKS) - {pf.CHECK_ENCODING}
    assert all(r.reason for r in skipped)


def test_skip_and_only_are_both_honoured():
    report = pf.run_preflight(library="cr6261_h1", arms=("alde",),
                              skip=set(pf.ALL_CHECKS))
    assert all(r.status == pf.SKIPPED for r in report.results)
    assert report.exit_code == 0


def test_the_report_serialises_to_json():
    report = pf.run_preflight(library="cr6261_h1", arms=("alde",), only={pf.CHECK_ENCODING})
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["library"] == "cr6261_h1"
    assert payload["exit_code"] == 1


def test_the_cli_exits_non_zero_on_a_failing_check(capsys):
    code = pf.main(["--library", "cr6261_h1", "--arms", "alde", "--only", "encoding"])
    assert code == 1
    out = capsys.readouterr().out
    assert "198" in out


def test_the_cli_exits_zero_when_nothing_applicable_fails(capsys):
    code = pf.main(["--library", "cr6261_h1", "--arms", "ridge", "--only", "encoding"])
    assert code == 0
    assert pf.NOT_APPLICABLE in capsys.readouterr().out


def test_the_cli_invokes_a_row_order_runner_against_a_real_landscape(capsys):
    """`--invariance-runner module:function` is how the check is invoked for real."""
    bad = pf.main(["--library", "cr6261_h1", "--arms", "alde", "--only", "row_order",
                   "--invariance-runner",
                   "tests.test_preflight:order_dependent_combo_runner"])
    assert bad == 1
    assert "shuffled" in capsys.readouterr().out

    good = pf.main(["--library", "cr6261_h1", "--arms", "alde", "--only", "row_order",
                    "--invariance-runner",
                    "tests.test_preflight:order_invariant_combo_runner"])
    assert good == 0


def test_a_bad_runner_spec_is_an_error_not_a_pass(capsys):
    assert pf.main(["--library", "cr6261_h1", "--arms", "alde", "--only", "row_order",
                    "--invariance-runner", "not-a-module-spec"]) == 2


def test_the_cli_can_emit_json(capsys):
    code = pf.main(["--library", "cr6261_h1", "--arms", "alde", "--only", "encoding",
                    "--json"])
    assert code == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["library"] == "cr6261_h1"


def test_the_docs_describe_every_check():
    text = (Path(pf.__file__).resolve().parents[1] / "docs" / "preflight.md").read_text(
        encoding="utf-8")
    for check in pf.ALL_CHECKS:
        assert check in text
