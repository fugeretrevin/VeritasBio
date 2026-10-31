"""
tests/test_benchmark_manifest_provenance.py - Regression tests for the three silent-failure
defects found while reconstructing the 545 ablation round.

Each class below pins one real incident:

  1. TestRequestedBaselineAvailability
     A requested ALDE baseline was reported "UNAVAILABLE" in a report row because ALDE_PATH
     was unset. Two full benchmark runs completed before anyone read the row. A baseline that
     was explicitly requested and cannot run must abort the run.

  2. TestManifestIsAuthoritative / TestManifestIdentifiesConfiguration
     Run manifests did not record which components were ablated, so recovering which of 62
     completed arms was which configuration required fingerprinting their OUTCOMES and
     matching them back. Outcomes must never be used to identify a configuration: the
     manifest has to say what was run, in machine-readable form.

  3. TestAblateComponentsParsing
     `--ablate-components` with an empty or separator-only value silently ran the full
     engine, and an unknown component name failed only once the first round was already
     under way.

  4. TestLabelProvenanceVocabulary
     `provenance=UNVERIFIED_INCLUDES_IMPUTED` is a single opaque string. The explicit
     vocabulary keeps measured, imputed, unknown-provenance, simulated, prospective and
     predicted labels non-interchangeable, and records the differential-bias limitation
     instead of assuming a shared contaminated universe cancels out.
"""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple
from unittest import mock

import benchmark_panoptes_isolated as bench


# ----------------------------------------------------------------------------------------------
# Shared fixtures
# ----------------------------------------------------------------------------------------------
_UNIVERSE_CACHE: Dict[Tuple[str, int], bench.Universe] = {}

# The synthetic generator enumerates a 4-character combo space, so it supports the
# four-site libraries. HER2 is ten-site and pool-restricted; its spec is read from the
# registry rather than simulated.
FOUR_SITE_LIBRARIES = ("trpb", "gb1")


def _library(key: str = "trpb") -> bench.LibrarySpec:
    """
    The library, resolved explicitly. There is no active-library state to configure: the
    spec is threaded into every call and `Universe` carries the one it was built from.
    """
    return bench.resolve_library(key)


def _tiny_universe(seed: int = 3, library: str = "trpb") -> bench.Universe:
    """A small, fast, clearly-synthetic universe. Values are generated, not biological."""
    key = (library, seed)
    if key not in _UNIVERSE_CACHE:
        _UNIVERSE_CACHE[key] = bench.build_synthetic_universe(
            _library(library), n_variants=300, seed=seed,
        )
    return _UNIVERSE_CACHE[key]


def _args(extra: Sequence[str] = ()) -> Any:
    # --synthetic-seed matches the seed _tiny_universe() generates with, so the manifest's
    # requested seed and the universe's own fingerprint agree as they do in a real run.
    return bench.build_parser().parse_args(
        ["--smoke-test", "--synthetic-size", "300", "--synthetic-seed", "3", "--seeds", "2",
         "--initial-size", "5", "--rounds", "2", "--batch-size", "5", *extra]
    )


_FAKE_REPO_STATE = {
    "git_available": True,
    "head": "0123456789abcdef0123456789abcdef01234567",
    "status_lines": [" M benchmark_panoptes_isolated.py"],
    "source_tree_sha256": "f" * 64,
    "files_hashed": 123,
}


def _build_manifest(
    *,
    components_disabled: Sequence[str] = (),
    outcomes: Sequence[bench.RunOutcome] = (),
    universe: bench.Universe = None,
    run_state: str = None,
    args: Any = None,
) -> Dict[str, Any]:
    universe = universe or _tiny_universe()
    args = args or _args()
    evaluator = bench.Evaluator(universe)
    seeds = [0, 1]
    initial_sets = {s: universe.variant_ids[s * 5:s * 5 + 5] for s in seeds}
    started = datetime(2026, 10, 31, 12, 0, 0, tzinfo=timezone.utc)
    return bench.build_manifest(
        run_id="run_test_0001",
        args=args,
        universe=universe,
        evaluator=evaluator,
        methods=["random"],
        method_status={"random": {"status": "available", "detail": "uniform random"}},
        method_configs={"random": {"strategy": "uniform without replacement"}},
        components_disabled=components_disabled,
        seeds=seeds,
        initial_sets=initial_sets,
        outcomes=outcomes,
        started=started,
        finished=started + timedelta(seconds=42),
        repository_state_before=_FAKE_REPO_STATE,
        repository_state_after=_FAKE_REPO_STATE,
        run_state=run_state,
    )


def _completed_outcome(method: str = "random", seed: int = 0, best: float = 0.5) -> bench.RunOutcome:
    """A minimal completed RunOutcome, used to vary outcomes while holding configuration fixed."""
    return bench.RunOutcome(
        method=method,
        seed=seed,
        status="completed",
        rounds=[bench.RoundMetrics(
            method=method, seed=seed, round_index=1, n_measurements=10, best_fitness=best,
            new_top1_hits=1, cum_top1_excluding_initial=1, cum_top1_including_initial=1,
            found_top1=True, round_seconds=0.1, cumulative_seconds=0.2, model_calls=1,
            token_usage="unknown", requested=5, accepted=5, rejected=0,
            status="completed", failure_reason="",
        )],
    )


# ==============================================================================================
# 3. --ablate-components handling
# ==============================================================================================
class TestAblateComponentsParsing(unittest.TestCase):
    def test_default_empty_string_means_no_ablation(self) -> None:
        self.assertEqual(bench.parse_ablate_components(""), ())
        self.assertEqual(bench.parse_ablate_components(None), ())

    def test_whitespace_only_value_is_rejected(self) -> None:
        """`--ablate-components " "` asked for an ablation and got the full engine."""
        with self.assertRaises(bench.ConfigurationError):
            bench.parse_ablate_components("   ")

    def test_separator_only_value_is_rejected(self) -> None:
        for spec in (",", ",,", " , "):
            with self.subTest(spec=spec):
                with self.assertRaises(bench.ConfigurationError):
                    bench.parse_ablate_components(spec)

    def test_empty_token_inside_a_list_is_rejected(self) -> None:
        for spec in ("dpp,,gflownet", "dpp,", ",dpp"):
            with self.subTest(spec=spec):
                with self.assertRaises(bench.ConfigurationError):
                    bench.parse_ablate_components(spec)

    def test_unknown_component_is_rejected_at_parse_time(self) -> None:
        with self.assertRaises(bench.ConfigurationError) as caught:
            bench.parse_ablate_components("dpp,gflownnet")
        message = str(caught.exception)
        self.assertIn("gflownnet", message)
        for known in bench.ABLATABLE_COMPONENTS:
            self.assertIn(known, message)

    def test_known_components_are_normalised(self) -> None:
        self.assertEqual(bench.parse_ablate_components("gflownet, dpp ,dpp"), ("dpp", "gflownet"))

    def test_component_vocabulary_matches_the_production_engine(self) -> None:
        """A name the benchmark accepts but math_worker rejects would fail mid-run."""
        from workers.math_worker import ABLATABLE_COMPONENTS as production

        self.assertEqual(set(bench.ABLATABLE_COMPONENTS), set(production))

    def test_component_states_report_enabled_and_disabled(self) -> None:
        states = bench.component_states(("dpp",))
        self.assertFalse(states["dpp"])
        self.assertTrue(states["gflownet"])
        self.assertEqual(set(states), set(bench.ABLATABLE_COMPONENTS))


# ==============================================================================================
# 1. A requested-but-unavailable baseline is fatal
# ==============================================================================================
class TestRequestedBaselineAvailability(unittest.TestCase):
    def _status_without_alde_path(self, extra: Sequence[str] = ()) -> Dict[str, Dict[str, Any]]:
        args = _args(["--methods", "alde", *extra])
        # ALDE lives OUTSIDE this repository, so the search also looks at sibling
        # directories. Point the sibling search somewhere empty and unset the variable,
        # so this test exercises the genuinely-missing case on any machine.
        environment = {k: v for k, v in os.environ.items()
                       if k not in ("ALDE_PATH", "ALDE_PYTHON")}
        with mock.patch.dict(os.environ, environment, clear=True), \
             mock.patch.object(bench, "_repository_family_roots",
                               lambda: [Path(tempfile.gettempdir()) / "no-such-family-root"]):
            _, status = bench.build_selector_factories(args, ["alde"], _library())
        return status

    def test_unset_alde_path_makes_the_requested_baseline_fatal(self) -> None:
        status = self._status_without_alde_path()
        self.assertEqual(status["alde"]["status"], "UNAVAILABLE")
        with self.assertRaises(bench.BaselineUnavailableError) as caught:
            bench.enforce_requested_methods_available(["alde"], status)
        message = str(caught.exception)
        self.assertIn("alde", message)
        self.assertIn("ALDE_PATH", message)

    def test_unavailable_reason_names_every_path_it_searched(self) -> None:
        """
        An overnight agent reported ALDE "not on disk anywhere" and skipped the arm. The
        refusal has to say where it looked, or the environment problem is unactionable.
        """
        status = self._status_without_alde_path()
        reason = status["alde"]["detail"]
        self.assertIn("Looked in", reason)
        self.assertIn("--alde-path", reason)
        self.assertIn("ALDE_PATH", reason)
        self.assertIn("sibling of the repository", reason)
        self.assertIn("NOT vendored", json.loads(status["alde"]["probe"])["note"])
        probe = json.loads(status["alde"]["probe"])
        self.assertTrue(probe["searched_checkouts"])
        for entry in probe["searched_checkouts"]:
            self.assertIn("source", entry)
            self.assertIn("exists", entry)

    def test_an_explicit_path_is_never_replaced_by_a_discovered_one(self) -> None:
        """
        A missing --alde-path must NOT fall through to a sibling clone. Running a
        different checkout from the one the operator named, and reporting it as theirs,
        is the same class of defect as the silent UNAVAILABLE row.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            real = Path(tmp_dir) / "family" / "ALDE"
            real.mkdir(parents=True)
            missing = Path(tmp_dir) / "not-the-checkout"
            args = _args(["--methods", "alde", "--alde-path", str(missing)])
            with mock.patch.object(bench, "_repository_family_roots",
                                   lambda: [Path(tmp_dir) / "family"]):
                resolution = bench.resolve_alde_paths(args)
            # The sibling clone exists and would otherwise have been selected.
            self.assertTrue(real.is_dir())
            self.assertIsNone(resolution["checkout"])
            self.assertIn("refusing to fall back", resolution["blocked"])
            self.assertIn(str(missing), resolution["blocked"])
            probe = bench.probe_alde(None, None, resolution)
            self.assertFalse(probe["available"])
            self.assertIn("refusing to fall back", probe["reason"])

    def test_an_explicit_interpreter_is_never_replaced_by_a_convention(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            checkout = Path(tmp_dir) / "family" / "ALDE"
            checkout.mkdir(parents=True)
            convention = Path(tmp_dir) / "family" / ".venv-alde" / "Scripts"
            convention.mkdir(parents=True)
            (convention / "python.exe").write_text("", encoding="utf-8")
            args = _args(["--methods", "alde",
                          "--alde-python", str(Path(tmp_dir) / "no-such-python.exe")])
            with mock.patch.object(bench, "_repository_family_roots",
                                   lambda: [Path(tmp_dir) / "family"]):
                resolution = bench.resolve_alde_paths(args)
            self.assertEqual(resolution["checkout"], checkout.resolve())
            self.assertIsNone(resolution["python_exe"])
            self.assertIn("refusing to fall back", resolution["blocked"])

    def test_worktree_search_reaches_outside_the_worktree(self) -> None:
        """
        A worktree sits at <repo>/.claude/worktrees/<name>, so the parent of the running
        tree is the worktrees folder. The search must also walk back to the real
        checkout's parent, which is where a sibling ALDE clone actually lives.
        """
        roots = bench._repository_family_roots()
        self.assertTrue(roots)
        parts = bench.REPO_ROOT.parts
        if ".claude" in parts and "worktrees" in parts:
            real_root_parent = Path(*parts[:parts.index(".claude")]).parent
            self.assertIn(real_root_parent, roots)
            self.assertNotEqual(roots[0], real_root_parent)  # both layouts are covered

    def test_available_baselines_are_accepted(self) -> None:
        args = _args(["--methods", "random,ridge"])
        _, status = bench.build_selector_factories(args, ["random", "ridge"], _library())
        bench.enforce_requested_methods_available(["random", "ridge"], status)

    def test_skipped_live_arm_is_also_fatal_when_requested(self) -> None:
        """Requesting panoptes_muse without the opt-in used to produce a quiet SKIPPED row."""
        args = _args(["--methods", "panoptes_muse"])
        _, status = bench.build_selector_factories(args, ["panoptes_muse"], _library())
        self.assertNotEqual(status["panoptes_muse"]["status"], "available")
        with self.assertRaises(bench.BaselineUnavailableError):
            bench.enforce_requested_methods_available(["panoptes_muse"], status)

    def test_override_keeps_the_report_only_behaviour(self) -> None:
        status = self._status_without_alde_path()
        bench.enforce_requested_methods_available(
            ["alde"], status, allow_unavailable=True
        )  # must not raise

    def test_cli_run_aborts_before_any_seed_and_records_why(self) -> None:
        """End to end: the run must not complete, and the folder must say what was missing."""
        environment = {k: v for k, v in os.environ.items()
                       if k not in ("ALDE_PATH", "ALDE_PYTHON")}
        with tempfile.TemporaryDirectory() as tmp_dir:
            with mock.patch.dict(os.environ, environment, clear=True), \
                 mock.patch.object(bench, "_repository_family_roots",
                                   lambda: [Path(tmp_dir) / "no-such-family-root"]), \
                 mock.patch.object(bench, "_load_repo_dotenv", lambda: None):
                code = bench.main([
                    "--smoke-test", "--synthetic-size", "300", "--methods", "random,alde",
                    "--seeds", "1", "--initial-size", "5", "--rounds", "1", "--batch-size", "5",
                    "--output-dir", tmp_dir,
                ])
            self.assertEqual(code, bench.EXIT_BASELINE_UNAVAILABLE)
            manifests = list(Path(tmp_dir).glob("*/manifest.json"))
            self.assertEqual(len(manifests), 1, "an aborted run must still leave a manifest")
            manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
            self.assertEqual(manifest["run"]["state"], "ABORTED")
            self.assertIn("alde", manifest["run"]["abort_reason"])
            self.assertEqual(manifest["run"]["seeds_attempted"], 0)
            self.assertFalse(list(manifests[0].parent.glob("round_metrics.csv")))


# ==============================================================================================
# 4. The explicit label-provenance vocabulary
# ==============================================================================================
class TestLabelProvenanceVocabulary(unittest.TestCase):
    EXPECTED = (
        "HISTORICAL_MEASURED",
        "HISTORICAL_IMPUTED",
        "HISTORICAL_UNKNOWN_PROVENANCE",
        "SIMULATED_NOISY_OUTCOME",
        "PROSPECTIVE_PHYSICAL_MEASUREMENT",
        "MODEL_PREDICTION",
    )

    def test_vocabulary_is_exactly_the_six_declared_kinds(self) -> None:
        self.assertEqual(
            tuple(member.name for member in bench.LabelProvenance), self.EXPECTED
        )
        self.assertEqual(
            {member.value for member in bench.LabelProvenance}, set(self.EXPECTED)
        )

    def test_members_are_not_interchangeable(self) -> None:
        members = list(bench.LabelProvenance)
        for i, left in enumerate(members):
            for right in members[i + 1:]:
                self.assertNotEqual(left, right)
        # Not a str subclass: an imputed label cannot silently compare equal to the
        # string a CSV happens to carry, which is how the old single stamp spread.
        self.assertNotIsInstance(bench.LabelProvenance.HISTORICAL_MEASURED, str)
        self.assertNotEqual(bench.LabelProvenance.HISTORICAL_MEASURED, "HISTORICAL_MEASURED")

    def test_every_member_is_defined(self) -> None:
        for member in bench.LabelProvenance:
            self.assertIn(member, bench.LABEL_PROVENANCE_DEFINITIONS)
            self.assertTrue(bench.LABEL_PROVENANCE_DEFINITIONS[member].strip())

    def test_legacy_statuses_map_onto_the_vocabulary(self) -> None:
        cases = {
            "VERIFIED_COLUMN": bench.LabelProvenance.HISTORICAL_MEASURED,
            "VERIFIED_EXTERNAL_LIST": bench.LabelProvenance.HISTORICAL_MEASURED,
            "UNVERIFIED_INCLUDES_IMPUTED": bench.LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE,
            "SYNTHETIC_NOT_BIOLOGICAL": bench.LabelProvenance.SIMULATED_NOISY_OUTCOME,
        }
        for legacy, expected in cases.items():
            with self.subTest(legacy=legacy):
                self.assertIs(bench.map_legacy_provenance_status(legacy), expected)

    def test_unknown_legacy_status_is_not_guessed(self) -> None:
        with self.assertRaises(bench.ProvenanceError):
            bench.map_legacy_provenance_status("PROBABLY_FINE")

    def test_manifest_keeps_the_legacy_stamp_and_the_differential_bias_limitation(self) -> None:
        manifest = _build_manifest()
        block = manifest["configuration"]["label_provenance"]
        self.assertEqual(block["status"], bench.LabelProvenance.SIMULATED_NOISY_OUTCOME.value)
        self.assertEqual(block["legacy_status"], "SYNTHETIC_NOT_BIOLOGICAL")
        joined = " ".join(block["limitations"]) + " " + block["differential_bias"] + " " + \
            json.dumps(block["open_blockers"])
        self.assertIn("differential", joined.lower())
        self.assertIn("TrpB", joined)
        self.assertIn("GB1", joined)

    def test_stamp_carries_both_vocabularies(self) -> None:
        universe = _tiny_universe()
        stamp = bench.provenance_stamp(universe)
        self.assertIn("SIMULATED_NOISY_OUTCOME", stamp)
        self.assertIn("SYNTHETIC_NOT_BIOLOGICAL", stamp)


# ==============================================================================================
# 2. The manifest is authoritative and machine-readable
# ==============================================================================================
class TestManifestIsAuthoritative(unittest.TestCase):
    def test_every_required_field_is_present(self) -> None:
        manifest = _build_manifest()
        bench.validate_manifest(manifest)  # must not raise
        for path in bench.MANIFEST_REQUIRED_PATHS:
            with self.subTest(path=path):
                node: Any = manifest
                for key in path.split("."):
                    self.assertIsInstance(node, dict, f"{path} is not reachable")
                    self.assertIn(key, node, f"{path} missing at '{key}'")
                    node = node[key]

    def test_validation_rejects_an_incomplete_manifest(self) -> None:
        manifest = _build_manifest()
        del manifest["configuration"]["components"]["disabled"]
        with self.assertRaises(bench.ManifestIncompleteError):
            bench.validate_manifest(manifest)

    def test_manifest_is_json_serialisable_without_a_fallback(self) -> None:
        """A manifest that only survives `default=str` is not machine-readable."""
        manifest = _build_manifest(outcomes=[_completed_outcome()])
        reloaded = json.loads(json.dumps(manifest))
        self.assertEqual(reloaded["configuration_fingerprint"], manifest["configuration_fingerprint"])

    def test_times_are_unambiguous_utc(self) -> None:
        manifest = _build_manifest()
        self.assertEqual(manifest["run"]["timezone"], "UTC")
        for key in ("started_utc", "finished_utc"):
            value = manifest["run"][key]
            self.assertTrue(value.endswith("+00:00"), value)
            self.assertEqual(datetime.fromisoformat(value).tzinfo, timezone.utc)
        self.assertEqual(manifest["run"]["duration_seconds"], 42.0)

    def test_run_state_is_recorded(self) -> None:
        failed = _build_manifest(outcomes=[])
        self.assertEqual(failed["run"]["state"], "FAILED")
        succeeded = _build_manifest(outcomes=[_completed_outcome(seed=0), _completed_outcome(seed=1)])
        self.assertEqual(succeeded["run"]["state"], "SUCCESS")
        self.assertEqual(succeeded["run"]["seeds_completed"], 2)
        partial = _build_manifest(outcomes=[
            _completed_outcome(seed=0),
            bench.RunOutcome(method="random", seed=1, status="failed", failure_reason="boom"),
        ])
        self.assertEqual(partial["run"]["state"], "PARTIAL")
        self.assertTrue(any("boom" in w for w in partial["run"]["warnings"]))

    def test_components_enabled_and_disabled_are_both_listed(self) -> None:
        manifest = _build_manifest(components_disabled=("dpp", "prior"))
        components = manifest["configuration"]["components"]
        self.assertEqual(components["disabled"], ["dpp", "prior"])
        self.assertEqual(
            set(components["enabled"]) | set(components["disabled"]),
            set(bench.ABLATABLE_COMPONENTS),
        )
        self.assertFalse(set(components["enabled"]) & set(components["disabled"]))

    def test_budget_and_batch_schedule_are_explicit(self) -> None:
        manifest = _build_manifest()
        budget = manifest["configuration"]["measurement_budget"]
        self.assertEqual(budget["initial_size"], 5)
        self.assertEqual(budget["rounds"], 2)
        self.assertEqual(budget["batch_size"], 5)
        self.assertEqual(budget["total_unique_measurements"], 15)
        schedule = budget["batch_schedule"]
        self.assertEqual([entry["round"] for entry in schedule], [0, 1, 2])
        self.assertEqual(sum(entry["n_measurements"] for entry in schedule), 15)

    def test_seed_block_identity_is_recorded(self) -> None:
        manifest = _build_manifest()
        seeding = manifest["configuration"]["seeding"]
        self.assertEqual(seeding["seed_values"], [0, 1])
        self.assertEqual(sorted(seeding["initial_sets_sha256"]), ["0", "1"])
        self.assertTrue(seeding["seed_block_sha256"])
        # Two different starting blocks must not share an identity.
        universe = _tiny_universe()
        other = bench.seed_block_identity(
            base_seed=0, n_seeds=2, seeds=[0, 1],
            initial_sets={0: universe.variant_ids[20:25], 1: universe.variant_ids[25:30]},
        )
        self.assertNotEqual(other["seed_block_sha256"], seeding["seed_block_sha256"])

    def test_code_commit_and_dirty_state_are_recorded(self) -> None:
        manifest = _build_manifest()
        code = manifest["configuration"]["code"]
        self.assertEqual(code["commit"], _FAKE_REPO_STATE["head"])
        self.assertTrue(code["dirty_tree"])
        self.assertEqual(code["dirty_file_count"], 1)
        clean = bench._code_identity({**_FAKE_REPO_STATE, "status_lines": []})
        self.assertFalse(clean["dirty_tree"])
        unknown = bench._code_identity({"git_available": False})
        self.assertIsNone(unknown["dirty_tree"], "an unknown tree state must not read as clean")

    def test_dataset_identity_and_eligibility_rules_are_recorded(self) -> None:
        manifest = _build_manifest()
        dataset = manifest["configuration"]["dataset"]
        for key in ("name", "source", "version", "sha256", "unique_variants",
                    "eligibility_rules", "excluded_records"):
            self.assertIn(key, dataset)
        self.assertTrue(dataset["eligibility_rules"])
        for key in ("duplicates_dropped", "imputed_excluded", "unmapped_dropped",
                    "malformed_or_non_finite_dropped", "unlabelled_excluded",
                    "subsampled_to"):
            self.assertIn(key, dataset["excluded_records"])

    def test_assay_transformation_and_threshold_definition_are_recorded(self) -> None:
        manifest = _build_manifest()
        assay = manifest["configuration"]["assay"]
        # The direction is the library's declaration, not a module constant: the old
        # FITNESS_DIRECTION global made one landscape's claim look universal.
        self.assertEqual(assay["fitness_direction"],
                         manifest["configuration"]["dataset"]["library"]["fitness_direction"])
        self.assertTrue(assay["transformation"])
        definition = assay["threshold_definition"]
        self.assertEqual(definition["percentile"], 99.0)
        self.assertIn("tie_handling", definition)
        # The threshold VALUE is outcome-derived, so it belongs with the results.
        self.assertNotIn("threshold_value", json.dumps(manifest["configuration"]))
        self.assertIsInstance(manifest["results"]["top1"]["threshold_value"], float)

    def test_acquisition_policy_and_hyperparameters_are_recorded(self) -> None:
        manifest = _build_manifest()
        policy = manifest["configuration"]["acquisition_policy"]
        self.assertEqual(policy["policy_id"], "hybrid_thompson")
        self.assertIn("formula", policy)
        self.assertIn("parameters", policy)
        self.assertIn("hyperparameters", manifest["configuration"]["model"])


class TestManifestIdentifiesConfiguration(unittest.TestCase):
    """
    The 62-arm reconstruction incident: a manifest must identify its configuration without
    anyone looking at what the run produced.
    """

    def test_two_arms_differing_only_in_ablation_have_distinguishable_manifests(self) -> None:
        identical_outcomes = [_completed_outcome(best=0.5)]
        full = _build_manifest(components_disabled=(), outcomes=identical_outcomes)
        ablated = _build_manifest(components_disabled=("dpp",), outcomes=identical_outcomes)

        self.assertNotEqual(full["configuration_fingerprint"], ablated["configuration_fingerprint"])
        self.assertNotEqual(
            full["configuration"]["components"]["disabled"],
            ablated["configuration"]["components"]["disabled"],
        )

    def test_fingerprint_ignores_outcomes_entirely(self) -> None:
        """Same configuration, different results: the configuration identity must not move."""
        quiet = _build_manifest(outcomes=[_completed_outcome(best=0.11)])
        loud = _build_manifest(outcomes=[_completed_outcome(best=0.99), _completed_outcome(seed=1)])
        self.assertEqual(quiet["configuration_fingerprint"], loud["configuration_fingerprint"])
        self.assertNotEqual(quiet["results"], loud["results"])

    def test_configuration_subtree_holds_no_outcome_fields(self) -> None:
        manifest = _build_manifest(outcomes=[_completed_outcome(best=0.99)])
        forbidden = {
            "best_fitness", "mean_best_fitness", "threshold_value", "global_best", "n_hits",
            "revealed_fitness", "top1_hits", "found_top1", "summary_rows", "results",
        }
        seen: List[str] = []

        def walk(node: Any) -> None:
            if isinstance(node, dict):
                for key, value in node.items():
                    if key in forbidden:
                        seen.append(key)
                    walk(value)
            elif isinstance(node, list):
                for item in node:
                    walk(item)

        walk(manifest["configuration"])
        self.assertEqual(seen, [], f"outcome field(s) leaked into the configuration: {seen}")

    def test_fingerprint_is_recomputable_from_the_manifest_alone(self) -> None:
        manifest = _build_manifest(components_disabled=("gflownet",))
        reloaded = json.loads(json.dumps(manifest))
        self.assertEqual(
            bench.configuration_fingerprint(reloaded["configuration"]),
            reloaded["configuration_fingerprint"],
        )

    def test_every_ablation_arm_of_a_sweep_is_uniquely_identified(self) -> None:
        """The 2^5 joint design: 32 arms, 32 distinct manifests, no outcomes consulted."""
        import itertools

        fingerprints = {}
        for size in range(len(bench.ABLATABLE_COMPONENTS) + 1):
            for combo in itertools.combinations(sorted(bench.ABLATABLE_COMPONENTS), size):
                manifest = _build_manifest(components_disabled=combo)
                fingerprint = manifest["configuration_fingerprint"]
                self.assertNotIn(
                    fingerprint, fingerprints,
                    f"arm {combo} collides with {fingerprints.get(fingerprint)}",
                )
                fingerprints[fingerprint] = combo
        self.assertEqual(len(fingerprints), 32)

    def test_differing_prior_spread_is_also_distinguishable(self) -> None:
        base = _build_manifest(args=_args())
        swept = _build_manifest(args=_args(["--prior-spread", "8.0"]))
        self.assertNotEqual(base["configuration_fingerprint"], swept["configuration_fingerprint"])

    def test_planted_epistasis_is_distinguishable_from_an_unplanted_landscape(self) -> None:
        """
        The planted pair is generator configuration, not a result. Without it in the
        manifest, a planted and an unplanted arm on the same seed would be indistinguishable
        except by their outcomes.
        """
        plain = bench.build_synthetic_universe(_library(), n_variants=300, seed=5)
        planted = bench.build_synthetic_universe(
            _library(), n_variants=300, seed=5, plant_sign_epistasis=True)
        plain_manifest = _build_manifest(universe=plain, args=_args(["--synthetic-seed", "5"]))
        planted_manifest = _build_manifest(
            universe=planted,
            args=_args(["--synthetic-seed", "5", "--plant-sign-epistasis"]),
        )
        self.assertNotEqual(
            plain_manifest["configuration_fingerprint"],
            planted_manifest["configuration_fingerprint"],
        )
        generator = planted_manifest["configuration"]["dataset"]["synthetic_generator"]
        self.assertTrue(generator["plant_sign_epistasis"])
        self.assertIsNotNone(generator["planted_interaction"])
        self.assertIn("position_a", generator["planted_interaction"])
        self.assertIsNone(
            plain_manifest["configuration"]["dataset"]["synthetic_generator"]
            ["planted_interaction"]
        )


class TestManifestRecordsTheLibrary(unittest.TestCase):
    """
    Multi-landscape runs are the same problem as the 62 ablation arms: a manifest that
    does not name its library cannot tell a TrpB arm from a GB1 one.

    Since the datasets refactor there is no active-library state to configure, so the
    manifest's library identity is taken from `Universe.library` - the library the
    universe was actually built from, which is stronger evidence than any global was.
    """

    def test_library_identity_reads_the_real_library_spec(self) -> None:
        for key in sorted(bench.LIBRARIES):
            with self.subTest(library=key):
                spec = bench.resolve_library(key)
                identity = bench._library_identity(spec)
                self.assertIn("Universe.library", identity["source"])
                self.assertEqual(identity["key"], spec.key)
                self.assertEqual(identity["label"], spec.label)
                self.assertEqual(identity["accession"], spec.accession)
                self.assertEqual(identity["positions"], list(spec.positions))
                self.assertEqual(identity["parent_residues"], spec.parent_residues)
                self.assertEqual(identity["n_sites"], spec.n_sites)
                self.assertEqual(identity["one_hot_features"], spec.one_hot_width)
                self.assertEqual(identity["selection_mode"], spec.selection_mode.value)
                # The switches that silently change what was run (38e98be, b681c57).
                self.assertEqual(identity["labelled_stage"], spec.labelled_stage)
                self.assertEqual(identity["pool_restricted"], spec.pool_restricted)
                self.assertEqual(identity["unlabelled_placeholder"],
                                 spec.unlabelled_placeholder)
                self.assertEqual(identity["fitness_direction"], spec.fitness_direction)
                self.assertEqual(identity["detection_floor"], spec.detection_floor)

    def test_manifest_library_comes_from_the_scored_universe(self) -> None:
        """
        The replacement for the old `registry == ACTIVE_LIBRARY` check: the recorded
        library must be the one the universe was built from, not one passed alongside it.
        A manifest describing a different library would misattribute every number.
        """
        universe = _tiny_universe(library="gb1")
        manifest = _build_manifest(universe=universe, args=_args(["--library", "gb1"]))
        recorded = manifest["configuration"]["dataset"]["library"]
        self.assertEqual(recorded["key"], universe.library.key)
        self.assertEqual(recorded["positions"], list(universe.library.positions))
        self.assertEqual(manifest["configuration"]["library_selection"]["resolved"],
                         universe.library.key)
        self.assertEqual(manifest["configuration"]["library_selection"]["requested"], "gb1")

    def test_a_library_mismatch_is_visible_in_the_manifest(self) -> None:
        """
        Asking for one library and scoring another's universe must be legible, not silent:
        `requested` and `resolved` disagree rather than both reporting the request.
        """
        universe = _tiny_universe(library="gb1")
        manifest = _build_manifest(universe=universe, args=_args(["--library", "trpb"]))
        selection = manifest["configuration"]["library_selection"]
        self.assertEqual(selection["requested"], "trpb")
        self.assertEqual(selection["resolved"], "gb1")
        self.assertNotEqual(selection["requested"], selection["resolved"])

    def test_her2_sparse_library_switches_are_recorded_and_not_defaults(self) -> None:
        """HER2 is the library where every switch is on; they must not read as defaults."""
        spec = bench.resolve_library("her2")
        identity = bench._library_identity(spec)
        self.assertEqual(identity["labelled_stage"], "train")
        self.assertTrue(identity["pool_restricted"])
        self.assertEqual(identity["n_sites"], 10)
        self.assertEqual(identity["one_hot_features"], 200)
        self.assertEqual(identity["selection_mode"], "pool_restricted")
        self.assertIsNotNone(identity["unlabelled_placeholder"])
        rules = bench._eligibility_rules(spec, _args())
        joined = " ".join(rules)
        self.assertIn("'train'", joined)
        self.assertIn("POOL-RESTRICTED", joined)
        self.assertIn("placeholder", joined)

    def test_censoring_and_label_policy_are_recorded_for_a_censored_library(self) -> None:
        """
        HER2 is censored from below (51.7% of labels at the assay floor) and 69% of its
        file is an unlabelled split. Both are properties of what the run was scored
        against, so both belong in the manifest rather than only in a log line.
        """
        universe = _tiny_universe()
        object.__setattr__(universe, "library", bench.resolve_library("her2"))
        object.__setattr__(universe, "censored_at_floor", 5_833)
        object.__setattr__(universe, "unlabelled_excluded", 25_114)
        try:
            assay = bench._assay_block(universe)
            self.assertEqual(assay["censoring"]["at_floor"], 5_833)
            self.assertEqual(assay["censoring"]["detection_floor"],
                             bench.LIBRARIES["her2"].detection_floor)
            self.assertIn("COUNTED", assay["censoring"]["policy"])
            self.assertEqual(assay["label_policy"]["unlabelled_excluded"], 25_114)
            self.assertEqual(assay["label_policy"]["labelled_stage"], "train")
            self.assertIn("three distinct states", assay["label_policy"]["policy"])
        finally:
            object.__setattr__(universe, "library", bench.resolve_library("trpb"))
            object.__setattr__(universe, "censored_at_floor", 0)
            object.__setattr__(universe, "unlabelled_excluded", 0)

    def test_two_libraries_have_distinguishable_manifests(self) -> None:
        fingerprints: Dict[str, str] = {}
        libraries: Dict[str, Any] = {}
        for key in FOUR_SITE_LIBRARIES:
            manifest = _build_manifest(
                universe=_tiny_universe(library=key), args=_args(["--library", key]),
            )
            fingerprints[key] = manifest["configuration_fingerprint"]
            libraries[key] = manifest["configuration"]["dataset"]["library"]
        self.assertNotEqual(fingerprints["trpb"], fingerprints["gb1"])
        self.assertNotEqual(libraries["trpb"]["positions"], libraries["gb1"]["positions"])
        self.assertNotEqual(libraries["trpb"]["parent_residues"],
                            libraries["gb1"]["parent_residues"])

    def test_unknown_imputed_count_does_not_print_as_a_number(self) -> None:
        """GB1 records no published imputed count; 'up to None imputed' must never appear."""
        universe = _tiny_universe(library="gb1")
        object.__setattr__(universe, "provenance_status", "UNVERIFIED_INCLUDES_IMPUTED")
        try:
            block = bench._label_provenance_block(universe)
            text = " ".join(block["limitations"])
            self.assertNotIn("None", text)
            self.assertIn("unverified", text.lower())
            self.assertEqual(block["library"], "gb1")
        finally:
            object.__setattr__(universe, "provenance_status", "SYNTHETIC_NOT_BIOLOGICAL")

    def test_declared_zero_imputed_is_still_unverified(self) -> None:
        """
        HER2 declares 0 imputed records. That is OUR declaration, not an upstream flag, so
        the labels stay HISTORICAL_UNKNOWN_PROVENANCE and the limitation says why. Reading
        our own registry back as a provenance claim is what the vocabulary exists to stop.
        """
        universe = _tiny_universe()
        object.__setattr__(universe, "library", bench.resolve_library("her2"))
        object.__setattr__(universe, "provenance_status", "UNVERIFIED_INCLUDES_IMPUTED")
        try:
            block = bench._label_provenance_block(universe)
            self.assertEqual(block["status"],
                             bench.LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE.value)
            self.assertEqual(block["published_imputed_count"], 0)
            text = " ".join(block["limitations"])
            self.assertIn("declares 0 imputed records", text)
            self.assertIn("not an upstream per-record flag", text)
            self.assertNotIn("up to 0", text)
        finally:
            object.__setattr__(universe, "library", bench.resolve_library("trpb"))
            object.__setattr__(universe, "provenance_status", "SYNTHETIC_NOT_BIOLOGICAL")

    def test_plot_title_names_the_library_it_ran_on(self) -> None:
        """A GB1 figure titled 'TrpB' outlives the folder that would have explained it."""
        trpb = bench.LIBRARIES["trpb"].label
        gb1 = bench.LIBRARIES["gb1"].label
        self.assertIn("TrpB", bench._plot_title(False, trpb))
        gb1_title = bench._plot_title(False, gb1)
        self.assertIn("GB1", gb1_title)
        self.assertNotIn("TrpB", gb1_title)
        synthetic_title = bench._plot_title(True, gb1)
        self.assertIn("SYNTHETIC", synthetic_title)
        self.assertIn("NOT biological", synthetic_title)
        self.assertIn("GB1", synthetic_title)  # names the geometry it was generated over


if __name__ == "__main__":
    unittest.main()
