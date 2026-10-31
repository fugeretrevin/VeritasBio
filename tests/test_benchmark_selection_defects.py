"""Regression tests for three defects that made a benchmark result non-comparable.

Each one was measured on real runs, and each test here fails against the code as it stood
before the corresponding fix.

Defect 1 - THE QUALIFICATION BAR NEVER TRANSFERRED.
    `--qualify-threshold` was an absolute fitness defaulting to 0.70. On TrpB, the only
    landscape it was tuned on, that is the 99.9602nd percentile of the eligible universe
    (64 of 159,999 variants). Fitness scales are not comparable across the corpus, so the
    same constant admitted 31.6% of GB1, 9.11% of HER2 and 100% of three bnAb landscapes
    whose assay floor is 7.0 on a -log10(Kd) scale. On 8 of 9 landscapes it sat BELOW that
    target's own top-1% threshold, making "qualifying" easier than being a hit.

Defect 2 - SHORT BATCHES MEANT THE BUDGET WAS NEVER MATCHED.
    Every baseline measured exactly 125 unique variants; our arms averaged as low as
    120.23 with a minimum of 107, because proposals outside the eligible universe were
    rejected and never replaced. The protocol requires exactly 125, so this was a protocol
    violation that biased AGAINST our arms.

Defect 3 - `_frontier_batch` WAS NOT POOL-RESTRICTION-AWARE.
    `propose()` was made pool-restriction-aware and `_frontier_batch` never was. On HER2
    the median library member has zero in-library Hamming-1 neighbours, so the frontier
    batch came back empty and `propose()` fell through to the baseline pool path: in 41 of
    48 seeds the frontier arm selected the IDENTICAL 125 variants as the reliability arm.
"""
from __future__ import annotations

import sys
import unittest
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

_UTC = datetime(2026, 10, 31, 12, 0, 0, tzinfo=timezone.utc)

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import benchmark_panoptes_isolated as bench  # noqa: E402


# ==============================================================================================
# Shared fixtures
# ==============================================================================================
_UNIVERSE_CACHE: Dict[str, bench.Universe] = {}


def _real_universe(key: str) -> bench.Universe:
    """A real landscape, loaded once. These tests assert on published fitness values."""
    if key not in _UNIVERSE_CACHE:
        library = bench.resolve_library(key)
        _UNIVERSE_CACHE[key] = bench.load_universe(
            library, library.dataset, "acknowledge-unflagged", None,
        )
    return _UNIVERSE_CACHE[key]


def _synthetic(key: str = "trpb", n: int = 400, seed: int = 5) -> bench.Universe:
    return bench.build_synthetic_universe(bench.resolve_library(key), n_variants=n, seed=seed)


def _view(universe: bench.Universe, batch_size: int = 5, n_rounds: int = 2) -> bench.CandidateView:
    return bench.CandidateView(
        variant_ids=universe.variant_ids,
        one_hot=universe.one_hot,
        index_of=universe.index_of,
        combos=universe.combos,
        library=universe.library,
        batch_size=batch_size,
        n_rounds=n_rounds,
    )


# ==============================================================================================
# Defect 1: the qualification bar must be a quantile, resolved against the eligible universe
# ==============================================================================================
class TestQualificationBarTransfers(unittest.TestCase):
    """The bar is stated as a quantile and resolved per landscape, from one values array."""

    def test_trpb_quantile_reproduces_the_legacy_absolute_exactly(self) -> None:
        """The development result must not move: q=99.96 is 0.6998 and the same 64 variants.

        This is the whole justification for the chosen default. If this fails, the fix has
        silently changed the landscape the parameter was tuned on.
        """
        universe = _real_universe("trpb")
        evaluator = bench.Evaluator(universe)

        self.assertEqual(evaluator.qualify_threshold_source, "quantile")
        self.assertEqual(evaluator.qualify_quantile, 99.96)
        self.assertAlmostEqual(evaluator.qualify_threshold_absolute, 0.6998, places=4)

        resolved = {
            v for v in universe.variant_ids
            if universe.fitness[v] >= evaluator.qualify_threshold_absolute
        }
        legacy = {
            v for v in universe.variant_ids
            if universe.fitness[v] >= bench.LEGACY_QUALIFY_ABSOLUTE
        }
        self.assertEqual(len(legacy), 64)
        # The IDENTICAL set, not merely the same count.
        self.assertEqual(resolved, legacy)
        self.assertEqual(evaluator.n_qualifying, 64)

    def test_the_resolved_absolute_transfers_to_each_landscape(self) -> None:
        """One quantile, three landscapes, three different absolutes."""
        expected = {"trpb": 0.6998, "gb1": 5.362, "her2": 6.379}
        for key, absolute in expected.items():
            with self.subTest(library=key):
                evaluator = bench.Evaluator(_real_universe(key))
                self.assertAlmostEqual(
                    evaluator.qualify_threshold_absolute, absolute, places=3,
                )

    def test_the_legacy_absolute_did_not_gate_gb1_or_her2(self) -> None:
        """Why the constant had to go: it admitted 3.16% of GB1 and 9.11% of HER2.

        HER2's 1,027 is why every HER2 seed logged "qualified before round 1".

        NOTE on the GB1 figure. The defect report quoted 31.6% for GB1. The measured share
        is 4,718 of 149,360 = 3.159%, i.e. the quoted figure is out by a factor of ten.
        That does not weaken the defect - 3.16% is still 316x the 0.0100 rarity the bar had
        on TrpB, and 4,718 qualifiers still make qualification near-certain from a
        25-variant start - but the number itself is pinned here as measured rather than as
        quoted, so a reader does not carry the wrong magnitude forward.
        """
        gb1 = _real_universe("gb1")
        n_gb1 = sum(
            1 for v in gb1.variant_ids
            if gb1.fitness[v] >= bench.LEGACY_QUALIFY_ABSOLUTE
        )
        self.assertEqual(n_gb1, 4718)
        self.assertEqual(len(gb1.variant_ids), 149360)
        self.assertAlmostEqual(n_gb1 / len(gb1.variant_ids), 0.0316, places=3)

        her2 = _real_universe("her2")
        n_her2 = sum(
            1 for v in her2.variant_ids
            if her2.fitness[v] >= bench.LEGACY_QUALIFY_ABSOLUTE
        )
        self.assertEqual(n_her2, 1027)
        self.assertEqual(len(her2.variant_ids), 11276)
        # And the quantile bar does gate: a tiny tail instead of 9.11%.
        evaluator = bench.Evaluator(her2)
        self.assertLess(evaluator.n_qualifying, 20)

    def test_the_bar_is_resolved_from_the_same_array_as_the_top1_threshold(self) -> None:
        """One definition of the eligible universe, so the two bars cannot drift apart.

        HER2 is the case that matters: its `labelled_stage="train"` exclusion removes
        25,114 unlabelled rows. If the qualification bar were computed over a different
        array it would be a bar on a different universe than the one being scored.
        """
        for key in ("trpb", "gb1", "her2"):
            with self.subTest(library=key):
                universe = _real_universe(key)
                evaluator = bench.Evaluator(universe)
                values = np.array(
                    [universe.fitness[v] for v in universe.variant_ids], dtype=np.float64,
                )
                self.assertEqual(
                    evaluator.qualify_threshold_absolute,
                    float(np.percentile(values, 99.96)),
                )
                self.assertEqual(evaluator.top1_threshold, float(np.percentile(values, 99.0)))
                # On every one of these the quantile bar is above the top-1% bar, which is
                # what makes qualifying a HARDER bar than being a hit, as intended.
                self.assertGreater(
                    evaluator.qualify_threshold_absolute, evaluator.top1_threshold,
                )

    def test_an_explicit_absolute_overrides_the_quantile_and_is_labelled(self) -> None:
        evaluator = bench.Evaluator(_synthetic(), qualify_threshold=0.70)
        self.assertEqual(evaluator.qualify_threshold_source, "explicit-override")
        self.assertEqual(evaluator.qualify_threshold_absolute, 0.70)
        self.assertIsNone(evaluator.qualify_quantile)

    def test_a_custom_quantile_is_honoured(self) -> None:
        universe = _synthetic()
        evaluator = bench.Evaluator(universe, qualify_quantile=50.0)
        values = np.array([universe.fitness[v] for v in universe.variant_ids])
        self.assertEqual(evaluator.qualify_quantile, 50.0)
        self.assertEqual(
            evaluator.qualify_threshold_absolute, float(np.percentile(values, 50.0)),
        )

    def test_the_two_forms_are_mutually_exclusive(self) -> None:
        with self.assertRaises(bench.BenchmarkError):
            bench.Evaluator(_synthetic(), qualify_quantile=99.0, qualify_threshold=0.7)

    def test_the_cli_refuses_both_forms_at_once(self) -> None:
        parser = bench.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(["--qualify-quantile", "99.9", "--qualify-threshold", "0.7"])

    def test_the_cli_default_is_the_quantile_and_the_absolute_has_no_default(self) -> None:
        """The absolute must not have a default: that default was the defect."""
        args = bench.build_parser().parse_args([])
        self.assertIsNone(args.qualify_threshold)
        self.assertIsNone(args.qualify_quantile)
        # Unset means the module default, which is a quantile.
        self.assertEqual(bench.DEFAULT_QUALIFY_QUANTILE, 99.96)

    def test_an_out_of_range_quantile_is_refused(self) -> None:
        with self.assertRaises(bench.BenchmarkError):
            bench.Evaluator(_synthetic(), qualify_quantile=101.0)

    def test_the_manifest_states_the_qualification_bar(self) -> None:
        """A run whose manifest does not state its qualification bar is not evidence."""
        universe = _synthetic()
        args = bench.build_parser().parse_args(
            ["--smoke-test", "--synthetic-size", "400", "--synthetic-seed", "5",
             "--seeds", "1", "--initial-size", "5", "--rounds", "1", "--batch-size", "5"]
        )
        evaluator = bench.Evaluator(
            universe,
            qualify_quantile=args.qualify_quantile,
            qualify_threshold=args.qualify_threshold,
        )
        manifest = _build_manifest(args=args, universe=universe, evaluator=evaluator)
        config = manifest["configuration"]["components"]

        self.assertEqual(config["qualify_quantile"], 99.96)
        self.assertEqual(config["qualify_threshold_source"], "quantile")
        self.assertAlmostEqual(
            config["qualify_threshold_absolute"], evaluator.qualify_threshold_absolute,
        )
        # The pre-existing key is reconciled, not duplicated with a different meaning: it
        # carries the resolved absolute the run was actually gated by.
        self.assertEqual(config["qualify_threshold"], config["qualify_threshold_absolute"])
        # The attainability caveat is recorded rather than silently resolved.
        self.assertEqual(config["qualify_qualifying_count"], evaluator.n_qualifying)
        self.assertFalse(config["qualify_minimum_set_floor"]["enabled"])

    def test_the_manifest_records_an_explicit_override_as_such(self) -> None:
        universe = _synthetic()
        args = bench.build_parser().parse_args(
            ["--smoke-test", "--synthetic-size", "400", "--synthetic-seed", "5",
             "--seeds", "1", "--initial-size", "5", "--rounds", "1", "--batch-size", "5",
             "--qualify-threshold", "0.7"]
        )
        evaluator = bench.Evaluator(
            universe,
            qualify_quantile=args.qualify_quantile,
            qualify_threshold=args.qualify_threshold,
        )
        config = _build_manifest(
            args=args, universe=universe, evaluator=evaluator,
        )["configuration"]["components"]
        self.assertEqual(config["qualify_threshold_source"], "explicit-override")
        self.assertIsNone(config["qualify_quantile"])
        self.assertEqual(config["qualify_threshold_absolute"], 0.70)

    def test_a_qualify_arm_cannot_be_built_without_a_resolved_bar(self) -> None:
        """An arm whose bar is unset would silently become the single-policy baseline."""
        args = bench.build_parser().parse_args([])
        with self.assertRaises(bench.BenchmarkError):
            bench.build_selector_factories(
                args, ["panoptes_qualify_then_optimize"], bench.resolve_library("trpb"),
                qualify_absolute=None,
            )


def _build_manifest(*, args, universe, evaluator):
    seeds = [0]
    return bench.build_manifest(
        run_id="run_defects_0001",
        args=args,
        universe=universe,
        evaluator=evaluator,
        methods=["random"],
        method_status={"random": {"status": "available", "detail": "uniform random"}},
        method_configs={"random": {"strategy": "uniform without replacement"}},
        components_disabled=(),
        seeds=seeds,
        initial_sets={s: universe.variant_ids[:5] for s in seeds},
        outcomes=[],
        started=_UTC,
        finished=_UTC,
        repository_state_before={
            "git_available": True,
            "head": "0" * 40,
            "status_lines": [],
            "source_tree_sha256": "f" * 64,
            "files_hashed": 1,
        },
        repository_state_after=None,
    )


# ==============================================================================================
# Defect 2: a short batch must be backfilled from the same policy
# ==============================================================================================
class _IneligibleThenEligible(bench.Selector):
    """A selector whose ranked list is mostly variants that do not exist.

    This reproduces GB1's real failure mode: ~6.65% of free proposals are not in the
    eligible universe, so the committed batch came up short. Here the rate is made extreme
    so the shortfall is unmistakable, and the ranked list is a stable PREFIX-EXTENDED list:
    asking for more returns the same order with more entries, exactly as a real ranking
    policy does.
    """

    name = "ineligible_then_eligible"

    def __init__(self, seed: int, pool: Sequence[str], bad_every: int = 2) -> None:
        super().__init__(seed)
        self.pool = list(pool)
        self.bad_every = bad_every
        self.requests: List[int] = []

    def propose(self, round_index, n, view, revealed):
        self.requests.append(n)
        out: List[str] = []
        pool_index = 0
        i = 0
        while len(out) < n and pool_index < len(self.pool):
            if i % self.bad_every == 0:
                # A genotype nobody assayed: rejected as not_in_eligible_universe.
                out.append(f"NOSUCHVARIANT{i:05d}")
            else:
                out.append(self.pool[pool_index])
                pool_index += 1
            i += 1
        return out


class _FixedShortList(bench.Selector):
    """A selector with nothing new to offer, however much it is asked for."""

    name = "fixed_short_list"

    def __init__(self, seed: int, scripted: Sequence[str]) -> None:
        super().__init__(seed)
        self.scripted = list(scripted)
        self.requests: List[int] = []

    def propose(self, round_index, n, view, revealed):
        self.requests.append(n)
        return list(self.scripted)


class TestShortBatchesAreBackfilled(unittest.TestCase):

    def test_a_short_batch_is_filled_to_the_full_budget(self) -> None:
        """The protocol requires exactly batch_size unique measurements per round."""
        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=10, n_rounds=2)
        initial = universe.variant_ids[:5]
        pool = [v for v in universe.variant_ids if v not in set(initial)]

        selector = _IneligibleThenEligible(0, pool)
        outcome = bench.run_single(selector, evaluator, view, initial, seed=0)

        self.assertEqual(outcome.status, "completed")
        for metrics in outcome.rounds[1:]:
            with self.subTest(round=metrics.round_index):
                # Without backfilling this was 5 of 10, every round.
                self.assertEqual(metrics.accepted, view.batch_size)
                self.assertEqual(metrics.status, "completed")
                self.assertGreater(metrics.backfill_attempts, 0)
                self.assertGreater(metrics.backfill_added, 0)
        # The full budget was spent.
        self.assertEqual(
            outcome.rounds[-1].n_measurements,
            len(initial) + view.n_rounds * view.batch_size,
        )

    def test_the_fill_asks_the_policy_for_more_rather_than_drawing_at_random(self) -> None:
        """A backfill must be the policy's own next choices, never a random pad."""
        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=10, n_rounds=1)
        initial = universe.variant_ids[:5]
        pool = [v for v in universe.variant_ids if v not in set(initial)]

        selector = _IneligibleThenEligible(0, pool)
        outcome = bench.run_single(selector, evaluator, view, initial, seed=0)

        # The selector was asked again, with a LARGER request than the batch size.
        self.assertGreater(len(selector.requests), 1)
        self.assertEqual(selector.requests[0], view.batch_size)
        self.assertGreater(selector.requests[1], view.batch_size)
        # Every accepted variant came from the policy's own ranked pool.
        accepted = [r.variant_id for r in outcome.trace if r.status == "accepted"
                    and r.round_index > 0]
        self.assertTrue(set(accepted).issubset(set(pool)))

    def test_backfilling_is_bounded(self) -> None:
        self.assertEqual(bench.BACKFILL_MAX_ATTEMPTS, 4)
        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=5, n_rounds=1)
        initial = universe.variant_ids[:5]

        # A policy that can never fill the batch: one eligible pick and nothing else.
        selector = _FixedShortList(0, [universe.variant_ids[20], "NOSUCHVARIANT"])
        outcome = bench.run_single(selector, evaluator, view, initial, seed=0)

        metrics = outcome.rounds[-1]
        self.assertEqual(metrics.status, "short_batch")
        self.assertLessEqual(metrics.backfill_attempts, bench.BACKFILL_MAX_ATTEMPTS)
        # And it is loud about it rather than padding to hide the shortfall.
        self.assertEqual(metrics.accepted, 1)
        self.assertIn("backfill", metrics.failure_reason)

    def test_a_re_offered_candidate_is_not_traced_twice(self) -> None:
        """A longer ranked list repeats its own prefix; that is not a new duplicate."""
        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=10, n_rounds=1)
        initial = universe.variant_ids[:5]
        pool = [v for v in universe.variant_ids if v not in set(initial)]

        outcome = bench.run_single(
            _IneligibleThenEligible(0, pool), evaluator, view, initial, seed=0,
        )
        rows = [r for r in outcome.trace if r.round_index == 1]
        ids = [r.variant_id for r in rows]
        self.assertEqual(len(ids), len(set(ids)), "a candidate was judged twice")

    def test_the_round_metrics_name_where_the_measurements_came_from(self) -> None:
        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=10, n_rounds=1)
        initial = universe.variant_ids[:5]
        pool = [v for v in universe.variant_ids if v not in set(initial)]

        outcome = bench.run_single(
            _IneligibleThenEligible(0, pool), evaluator, view, initial, seed=0,
        )
        source = outcome.rounds[-1].selection_source
        self.assertIn("backfilled", source)

    def test_the_round_metrics_csv_carries_the_fill_columns(self) -> None:
        import csv
        import tempfile

        universe = _synthetic(n=400)
        evaluator = bench.Evaluator(universe)
        view = _view(universe, batch_size=10, n_rounds=1)
        initial = universe.variant_ids[:5]
        pool = [v for v in universe.variant_ids if v not in set(initial)]
        outcome = bench.run_single(
            _IneligibleThenEligible(0, pool), evaluator, view, initial, seed=0,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "round_metrics.csv"
            bench.write_round_metrics_csv(path, [outcome], "stamp")
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        for column in ("backfill_attempts", "backfill_added", "selection_source"):
            self.assertIn(column, rows[0])

    def test_the_manifest_records_the_backfill_rule(self) -> None:
        universe = _synthetic()
        args = bench.build_parser().parse_args(
            ["--smoke-test", "--synthetic-size", "400", "--synthetic-seed", "5",
             "--seeds", "1", "--initial-size", "5", "--rounds", "1", "--batch-size", "5"]
        )
        manifest = _build_manifest(
            args=args, universe=universe, evaluator=bench.Evaluator(universe),
        )
        backfill = manifest["configuration"]["measurement_budget"]["backfill"]
        self.assertTrue(backfill["enabled"])
        self.assertEqual(backfill["max_attempts"], bench.BACKFILL_MAX_ATTEMPTS)
        # No unnamed secondary rule is in play.
        self.assertIsNone(backfill["secondary_rule"])


# ==============================================================================================
# Defect 3: the post-qualification policies must be pool-restriction-aware
# ==============================================================================================
def _her2_like_pool_universe() -> bench.Universe:
    """A sparse pool-restricted universe whose members have no Hamming-1 neighbours.

    This is HER2's shape, built small so the test is fast: every genotype differs from
    every other at two or more positions, so the Hamming-1 enumeration that
    `_frontier_batch` used to perform finds nothing in the library.
    """
    library = bench.resolve_library("her2")
    n_sites = library.n_sites
    alphabets = library.allowed_aa_per_position or tuple("ACDEFGHIKLMNPQRSTVWY"
                                                         for _ in range(n_sites))
    rng = np.random.default_rng(11)

    combos: List[str] = []
    for index in range(60):
        # Vary EVERY position per member, so no two members are Hamming-1 apart.
        combo = "".join(
            alphabets[s][int(rng.integers(len(alphabets[s])))] for s in range(n_sites)
        )
        combos.append(combo)
    combos = sorted(set(combos))

    ids, kept, scores = [], [], []
    for i, combo in enumerate(combos):
        vid = bench._canonical_variant_id(combo, library)
        if not vid or vid in ids:
            continue
        ids.append(vid)
        kept.append(combo)
        scores.append(float(i) / len(combos))

    return bench.Universe(
        variant_ids=ids,
        combos=dict(zip(ids, kept)),
        fitness=dict(zip(ids, scores)),
        one_hot=bench._one_hot_matrix(kept, library),
        index_of={v: i for i, v in enumerate(ids)},
        library=library,
        fingerprint=bench.DatasetFingerprint(
            path="<synthetic:her2-like>", sha256="0" * 64, size_bytes=0,
            rows_parsed=len(ids), unique_variants=len(ids),
        ),
        provenance_status="SYNTHETIC_NOT_BIOLOGICAL",
        provenance_detail="sparse pool-restricted fixture for defect 3",
        planted_interaction=None,
        excluded_imputed=0,
        duplicates_dropped=0,
        synthetic=True,
    )


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


class TestFrontierIsPoolRestrictionAware(unittest.TestCase):

    def setUp(self) -> None:
        self.universe = _her2_like_pool_universe()
        self.library = self.universe.library
        self.assertTrue(self.library.pool_restricted, "fixture must be pool-restricted")
        self.view = _view(self.universe, batch_size=5, n_rounds=2)

    def test_the_fixture_really_has_no_in_library_hamming1_neighbours(self) -> None:
        """Establishes the precondition: this is HER2's sparsity, not a contrived empty pool."""
        combos = list(self.universe.combos.values())
        pairs = sum(
            1 for i, a in enumerate(combos) for b in combos[i + 1:] if _hamming(a, b) == 1
        )
        self.assertEqual(pairs, 0)

    def _selector(self, probe: bool = False) -> bench.PanoptesSelector:
        return bench.PanoptesSelector(
            seed=0,
            wt_sequence="M" * 300,
            metadata=None,
            use_live_muse=False,
            pdb_path=None,
            library=self.library,
            qualify_threshold=0.0,          # already qualified, so the second phase runs
            probe_after_qualify=probe,
        )

    def test_frontier_selects_from_the_enumerated_pool(self) -> None:
        """Before the fix this returned [] and the arm became the reliability arm."""
        revealed = bench.RevealedState(
            measured={v: self.universe.fitness[v] for v in self.universe.variant_ids[:12]}
        )
        batch = self._selector()._frontier_batch(5, self.view, revealed)

        self.assertEqual(len(batch), 5)
        untested = set(self.view.untested(revealed))
        self.assertTrue(set(batch).issubset(untested))

    def test_frontier_fills_the_requested_batch_on_a_sparse_library(self) -> None:
        """The 1-variant round-1 batch that cost 24 of 125 measurements must not recur."""
        revealed = bench.RevealedState(
            measured={v: self.universe.fitness[v] for v in self.universe.variant_ids[:12]}
        )
        for n in (1, 5, 10):
            with self.subTest(n=n):
                self.assertEqual(len(self._selector()._frontier_batch(n, self.view, revealed)), n)

    def test_frontier_does_not_fall_through_to_the_baseline_pool_path(self) -> None:
        """propose() must return the frontier policy's batch, not the baseline's.

        This is the measured defect: in 41 of 48 HER2 seeds the frontier arm selected the
        identical 125 variants as the reliability arm because its batch came back empty.
        """
        revealed = bench.RevealedState(
            measured={v: self.universe.fitness[v] for v in self.universe.variant_ids[:12]}
        )
        selector = self._selector()
        batch = selector.propose(1, 5, self.view, revealed)

        self.assertEqual(len(batch), 5)
        self.assertEqual(selector.qualified_at_round, 1)
        self.assertEqual(selector.last_batch_source, "frontier")
        self.assertNotIn("baseline", selector.last_batch_source)
        self.assertFalse(
            [n for n in selector.notes if "POLICY_FALLBACK" in n],
            "the frontier arm silently fell back to the baseline policy",
        )

    def test_an_exhausted_pool_fails_loudly_instead_of_degrading(self) -> None:
        """A policy that cannot be expressed must say so, not become another arm."""
        revealed = bench.RevealedState(
            measured={v: self.universe.fitness[v] for v in self.universe.variant_ids}
        )
        with self.assertRaises(bench.BenchmarkError) as caught:
            self._selector()._frontier_batch(5, self.view, revealed)
        self.assertIn("pool-restricted", str(caught.exception))

    def test_the_conditional_probe_is_pool_restriction_aware_too(self) -> None:
        """`_conditional_probe_batch` has the same shape and the same defect.

        On a sparse library its residue scan enumerates genotypes nobody assayed. In
        pool-restricted mode it must either find the conditional slice in the pool or fail
        loudly - never return an empty batch that sends propose() to the baseline.
        """
        revealed = bench.RevealedState(
            measured={v: self.universe.fitness[v] for v in self.universe.variant_ids[:12]}
        )
        selector = self._selector(probe=True)
        try:
            batch = selector._conditional_probe_batch(5, self.view, revealed)
        except bench.BenchmarkError as exc:
            # The honest outcome on this fixture: no in-pool single-position variants
            # exist, so the policy genuinely cannot be expressed and says so.
            self.assertIn("conditional probe", str(exc))
            self.assertIn("pool-restricted", str(exc))
        else:
            self.assertTrue(batch)
            self.assertTrue(set(batch).issubset(set(self.view.untested(revealed))))

    def test_a_generative_library_keeps_the_hamming1_frontier(self) -> None:
        """The fix must not change the saturated-library behaviour the arm was built for."""
        universe = _synthetic("trpb", n=400, seed=5)
        self.assertFalse(universe.library.pool_restricted)
        view = _view(universe, batch_size=5)
        revealed = bench.RevealedState(
            measured={v: universe.fitness[v] for v in universe.variant_ids[:12]}
        )
        selector = bench.PanoptesSelector(
            seed=0, wt_sequence="M" * 300, metadata=None, use_live_muse=False,
            pdb_path=None, library=universe.library, qualify_threshold=0.0,
        )
        batch = selector._frontier_batch(5, view, revealed)
        best = max(revealed.measured, key=lambda v: revealed.measured[v])
        best_combo = universe.combos[best]
        for vid in batch:
            with self.subTest(variant=vid):
                self.assertEqual(_hamming(universe.combos[vid], best_combo), 1)


if __name__ == "__main__":
    unittest.main()
