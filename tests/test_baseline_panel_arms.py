"""
The baseline panel's arms, and the guard that stops one of them becoming the reference.

Every comparative claim in this project was against ONE baseline, `ridge`. The arms added
on 2026-11-01 (`gbm_onehot`, `knn_hamming`, `gp_onehot`, `greedy_measured`) measure how
strong ridge is relative to its neighbours. Pre-registration:
docs/prereg-baseline-panel.md.

Two kinds of test here, and the second kind is the one that matters most.

1. CONTRACT. A panel arm whose budget accounting or initial-sample handling differs from
   `ridge`'s is not an additional arm in the same comparison - it is a different protocol
   wearing the same budget, and the paired test over seeds is then fictional. So the arms
   are tested for the properties the pairing rests on, against `ridge` as the oracle for
   what those properties are, rather than against a hand-written expectation.

2. NO PROMOTION. docs/HANDOFF.md 5.4 makes promoting a second model to the status of
   reference baseline an amendment to a frozen hypothesis, which requires an explicit dated
   amendment record. A tuple edit is not an amendment record, and the edit is a one-liner
   that would pass review unnoticed, so `DEFAULT_METHODS` is pinned here. If this test
   fails, either an amendment has been written and this test should be updated to cite it,
   or a frozen hypothesis has just been amended by accident.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from typing import Dict, List

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import benchmark_panoptes_isolated as bench  # noqa: E402

#: The arms added for the panel. `random`, `ridge` and `random_forest` pre-date it.
NEW_ARMS = ("gbm_onehot", "knn_hamming", "gp_onehot", "greedy_measured")


def _universe(key: str = "trpb", n: int = 300, seed: int = 11) -> bench.Universe:
    return bench.build_synthetic_universe(bench.resolve_library(key), n_variants=n, seed=seed)


def _view(universe: bench.Universe, batch_size: int = 8, n_rounds: int = 2) -> bench.CandidateView:
    return bench.CandidateView(
        variant_ids=universe.variant_ids,
        one_hot=universe.one_hot,
        index_of=universe.index_of,
        combos=universe.combos,
        library=universe.library,
        batch_size=batch_size,
        n_rounds=n_rounds,
    )


def _revealed(universe: bench.Universe, n: int, seed: int = 3) -> bench.RevealedState:
    rng = np.random.default_rng(seed)
    picks = rng.choice(len(universe.variant_ids), size=n, replace=False)
    return bench.RevealedState(
        measured={universe.variant_ids[int(i)]: float(universe.fitness[universe.variant_ids[int(i)]])
                  for i in picks}
    )


def _selector(name: str, seed: int) -> bench.Selector:
    classes = {
        "random": bench.RandomSelector,
        "ridge": bench.RidgeSelector,
        "random_forest": bench.RandomForestSelector,
        "gbm_onehot": bench.GbmOneHotSelector,
        "knn_hamming": bench.KnnHammingSelector,
        "gp_onehot": bench.GpOneHotSelector,
        "greedy_measured": bench.GreedyMeasuredSelector,
    }
    return classes[name](seed)


def _hamming(a: str, b: str) -> int:
    """Hamming distance over FULL-LENGTH per-site genotypes.

    Rule 2.1: mutation strings omit wild-type positions, and comparing one of those by
    position silently truncates the comparison - a bug that already inverted a published
    conclusion about hit diversity here. `CandidateView.combos` carries full-length
    genotypes, so equal length is ASSERTED rather than hoped for, and nothing is zipped.
    """
    if len(a) != len(b):
        raise AssertionError(
            f"genotypes {a!r} and {b!r} differ in length ({len(a)} vs {len(b)}); "
            f"wild-type positions must be present, not omitted"
        )
    return sum(1 for i in range(len(a)) if a[i] != b[i])


class PanelArmsShareRidgesContract(unittest.TestCase):
    """The properties the paired test over seeds rests on, with `ridge` as the oracle."""

    def test_round_one_with_nothing_revealed_is_ridges_own_batch(self):
        """
        Identical first batch, not merely an identically sized one.

        Every panel arm delegates round 1 to `RandomSelector(self.seed)` exactly as `ridge`
        does. If one of them drew its own first batch instead, the arms would no longer
        share the shared starting sample's successor and the pairing would break at the
        first round rather than visibly.
        """
        universe = _universe()
        view = _view(universe)
        empty = bench.RevealedState()
        expected = _selector("ridge", 7).propose(1, 8, view, empty)
        for name in NEW_ARMS:
            with self.subTest(arm=name):
                self.assertEqual(_selector(name, 7).propose(1, 8, view, empty), expected)

    def test_budget_is_matched_per_arm(self):
        """Exactly n distinct untested variants, like ridge - no short batch, no duplicate."""
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 40)
        for name in ("ridge",) + NEW_ARMS:
            with self.subTest(arm=name):
                batch = _selector(name, 2).propose(2, 25, view, revealed)
                self.assertEqual(len(batch), 25)
                self.assertEqual(len(set(batch)), 25, "an arm proposed a duplicate")
                self.assertFalse(set(batch) & set(revealed.measured),
                                 "an arm re-proposed an already-measured variant")

    def test_batch_is_capped_by_the_pool_not_by_the_request(self):
        """`min(n, len(pool))`, as ridge does: a request larger than the pool is not an error."""
        universe = _universe(n=60)
        view = _view(universe)
        revealed = _revealed(universe, 50, seed=4)
        pool = len(view.untested(revealed))
        for name in ("ridge",) + NEW_ARMS:
            with self.subTest(arm=name):
                batch = _selector(name, 1).propose(2, pool + 15, view, revealed)
                self.assertEqual(len(batch), pool)

    def test_arms_are_deterministic_given_the_seed(self):
        """Two constructions at one seed must agree, or no run is reproducible."""
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 40)
        for name in NEW_ARMS:
            with self.subTest(arm=name):
                first = _selector(name, 5).propose(2, 12, view, revealed)
                second = _selector(name, 5).propose(2, 12, view, revealed)
                self.assertEqual(first, second)

    def test_no_arm_consults_the_evaluator(self):
        """
        A panel arm sees only `CandidateView` and `RevealedState`.

        Asserted structurally: `propose` is called with no evaluator in scope, so an arm
        that reached for true fitness would raise rather than quietly score itself against
        outcomes (rule 2.9).
        """
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 32)
        for name in NEW_ARMS:
            with self.subTest(arm=name):
                selector = _selector(name, 0)
                batch = selector.propose(2, 10, view, revealed)
                self.assertTrue(set(batch).issubset(set(view.variant_ids)))


class GreedyMeasuredIsAPureLocalExploit(unittest.TestCase):
    """`greedy_measured` has no model, so its rule is checkable exactly."""

    def test_it_picks_the_best_measured_variants_neighbourhood_nearest_first(self):
        """
        The batch is ordered by Hamming distance to the best measured variant, and no
        candidate at a larger distance is taken before an available nearer one.
        """
        universe = _universe(n=500, seed=21)
        view = _view(universe)
        revealed = _revealed(universe, 30, seed=9)
        best = max(revealed.measured, key=lambda v: revealed.measured[v])
        batch = bench.GreedyMeasuredSelector(0).propose(2, 20, view, revealed)

        distances = [_hamming(view.combos[v], view.combos[best]) for v in batch]
        self.assertEqual(distances, sorted(distances),
                         "the batch is not ordered nearest-first from the best lead")

        taken = max(distances)
        pool = set(view.untested(revealed))
        nearer_available = {
            v for v in pool
            if _hamming(view.combos[v], view.combos[best]) < taken
        }
        self.assertTrue(nearer_available.issubset(set(batch)),
                        "a farther candidate was taken while a nearer one was still free")

    def test_it_makes_no_model_call(self):
        """The control for 'is the surrogate doing anything?' must have no surrogate."""
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 30)
        selector = bench.GreedyMeasuredSelector(0)
        selector.propose(2, 10, view, revealed)
        self.assertEqual(selector.model_calls, 0)

    def test_a_wild_type_omitting_genotype_is_refused_not_truncated(self):
        """
        Rule 2.1, enforced rather than documented.

        A genotype shorter than the library's site count is a mutation string with wild-type
        positions omitted. Accepting one would understate every distance this arm computes,
        which is the shape of the bug that already inverted a hit-diversity conclusion here.
        """
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 30)
        best = max(revealed.measured, key=lambda v: revealed.measured[v])
        combos: Dict[str, str] = dict(view.combos)
        combos[best] = combos[best][:-1]          # drop a site, as a mutation string would
        broken = bench.CandidateView(
            variant_ids=view.variant_ids, one_hot=view.one_hot, index_of=view.index_of,
            combos=combos, library=view.library, batch_size=view.batch_size,
            n_rounds=view.n_rounds,
        )
        with self.assertRaises(bench.BenchmarkError) as caught:
            bench.GreedyMeasuredSelector(0).propose(2, 10, broken, revealed)
        self.assertIn("Wild-type positions must be PRESENT", str(caught.exception))


class GpArmFitsItsNoiseTerm(unittest.TestCase):
    def test_the_noise_level_is_learned_and_not_asserted(self):
        """
        A noise-free GP interpolates the measurements exactly and its mean surface is
        dominated by replicate noise, which is the whole reason this arm carries a
        WhiteKernel. The fitted level is recorded per round, so a future reader can see
        whether the GP concluded the assay was noisy.
        """
        universe = _universe()
        view = _view(universe)
        revealed = _revealed(universe, 40)
        selector = bench.GpOneHotSelector(0)
        selector.propose(2, 10, view, revealed)
        self.assertTrue(any("WhiteKernel" in note for note in selector.notes),
                        f"no fitted kernel recorded; notes were {selector.notes}")


class NoBaselineIsPromoted(unittest.TestCase):
    """
    The guard. HANDOFF 5.4 makes promoting a second model to reference-baseline status an
    amendment to a frozen hypothesis, which needs an explicit dated amendment record.
    """

    def test_the_panel_arms_are_registered_but_not_default(self):
        for name in NEW_ARMS:
            with self.subTest(arm=name):
                self.assertIn(name, bench.ALL_METHODS, "arm is not runnable")
                self.assertNotIn(
                    name, bench.DEFAULT_METHODS,
                    f"{name} has entered DEFAULT_METHODS. That promotes it toward reference "
                    f"status and amends a frozen hypothesis (HANDOFF 5.4). An amendment "
                    f"needs a dated record, not a tuple edit.",
                )

    def test_default_methods_is_exactly_what_the_frozen_study_registered(self):
        self.assertEqual(
            bench.DEFAULT_METHODS,
            ("random", "ridge", "random_forest", "panoptes_default"),
        )

    def test_ridge_is_still_the_reference_in_the_panel_order(self):
        """
        `random` first, then `ridge`. The order is load-bearing in two ways: a panel without
        `random` cannot calibrate the others - nothing then says whether a landscape's
        top-1% count can discriminate methods at all, which is the HER2 informativeness
        lesson - and a panel whose second entry is not `ridge` is no longer comparable to
        any figure already published in this project.
        """
        self.assertEqual(bench.BASELINE_PANEL_METHODS[:2], ("random", "ridge"))
        for name in NEW_ARMS:
            self.assertIn(name, bench.BASELINE_PANEL_METHODS)
        self.assertIn("panoptes_default", bench.BASELINE_PANEL_METHODS)

    def test_every_panel_arm_resolves_to_a_factory(self):
        """A silently unavailable arm is how our strongest baseline once left a 48-seed run."""
        library = bench.resolve_library("trpb")
        args = bench.build_parser().parse_args(["--library", "trpb"])
        methods: List[str] = [m for m in bench.BASELINE_PANEL_METHODS if m != "panoptes_default"]
        factories, status = bench.build_selector_factories(args, methods, library)
        for name in methods:
            with self.subTest(arm=name):
                self.assertIn(name, factories, f"{name}: {status.get(name)}")
                self.assertEqual(status[name]["status"], "available")


if __name__ == "__main__":
    unittest.main()
