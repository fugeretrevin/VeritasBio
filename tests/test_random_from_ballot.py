"""The coin must actually be a coin, and the record of it must read the code path.

Why this exists. `scripts/analysis/ranking_decomposition.py` measured, over all 48 CreiLOV
campaigns, that every stage after the GP residual makes the ranking's Spearman rho against
true fitness WORSE - GP residual +0.148, composite mu(x) -0.032, post-acquisition -0.104,
post-Pareto -0.121 - and that the committed batch scores 0.81x its own ballot's base rate.
That is a claim about the RANKING STACK, and no existing arm tests it: `--ablate-components`
can remove stages one at a time but cannot remove the stack, because something still has to
order the ballot. `panoptes_random_from_ballot` removes it entirely by drawing uniformly from
the post-veto ballot.

These tests assert the MECHANISM - that the coin reaches the selection, that nothing the
ranking stack computes can move it, and that the manifest's record of what ran is read from
the engine's return value rather than from the flag that asked for it (HANDOFF rule 2.12).
They never assert a candidate count, a hit count or a pass rate.

`tests/test_spatial_filter_ablation.py` is the worked example, including its warning: its
first version passed `evidence_records=[], strategy=None, wt_sequence=""`, which raised
TypeError before the stage under test was reached, so BOTH arms recorded zero and the test
passed while proving nothing. Every behavioural test below is therefore paired with a
NON-VACUITY guard that fails if the invocation stops exercising the path, or if the control
arm stops responding to the perturbation the test applies.
"""
from __future__ import annotations

import unittest
from unittest import mock

import numpy as np

from workers import math_worker

_BENCH = None


def _benchmark_module():
    """Import `benchmark_panoptes_isolated` by path, once.

    It must be registered in `sys.modules` BEFORE `exec_module`, or dataclass field
    resolution fails with `'NoneType' object has no attribute '__dict__'` - `dataclasses`
    looks the defining module up by name while the class body executes.
    """
    global _BENCH
    if _BENCH is None:
        import importlib.util
        import sys

        spec = importlib.util.spec_from_file_location(
            "benchmark_panoptes_isolated", "benchmark_panoptes_isolated.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules["benchmark_panoptes_isolated"] = module
        spec.loader.exec_module(module)
        _BENCH = module
    return _BENCH


def _trpb_selection_inputs():
    """A realistic invocation that REACHES the batch-construction branch.

    TrpB rather than CreiLOV because it is the saturated 4-site library the self-checks
    already anchor on, and because a test that needs a 15-site library to run is a test
    nobody runs. The pool is supplied, which is the pool-restricted path: it skips the
    GFlowNet sampler and therefore makes the ballot deterministic and nameable, which is
    what lets these tests compare two selections for equality at all.
    """
    from controller.gatekeeper import Variant
    from controller.prompts import ResidueCluster, SearchStrategy

    bench = _benchmark_module()
    library = bench.LIBRARIES["trpb"]
    wt = bench.resolve_library_metadata(library).sequence
    strategy = SearchStrategy(
        target_protein="TrpB", target_accession="P0A879",
        hypothesis="exercise the batch-construction branch",
        target_clusters=[ResidueCluster(
            cluster_id="sites", positions=list(library.positions),
            rationale="the ordered library")],
        max_mutation_depth=3,
    )
    # The pool is EVERY single substitution the library declares, plus a handful of
    # multi-mutants. Deliberately large: the first version of this fixture listed thirteen
    # variants, of which the novelty filter removed four and the spatial viability veto a
    # further five, leaving a post-veto ballot of FOUR. With a batch of four the ranked arm
    # then commits the entire ballot, so "the coin picked a different set" and "the ranking
    # moved under a perturbation" were both unmeasurable - a ballot smaller than the batch
    # makes every selection rule identical. Singles cannot be vetoed (the veto scores
    # pairwise distance between mutated positions, so a depth-1 candidate has no pair), so
    # they guarantee the ballot stays well above any batch size used here.
    pool = []
    for slot, position in enumerate(library.positions):
        parent = wt[position - 1]
        for residue in library.alphabet_for(slot):
            if residue != parent:
                pool.append(Variant(notation=f"{parent}{position}{residue}"))
    pool += [Variant(notation=n) for n in (
        "V183L:F184I", "V183L:S228G", "F184I:V227K", "V227K:S228G",
        "V183L:F184I:S228G", "F184I:V227K:S228G", "V183A:F184L:V227K",
    )]
    evidence = [
        {"variant_sequence": Variant(notation="V183L").apply_to(wt),
         "computed_score": 0.40},
        {"variant_sequence": Variant(notation="F184I").apply_to(wt),
         "computed_score": 0.55},
        {"variant_sequence": Variant(notation="S228G").apply_to(wt),
         "computed_score": 0.31},
        {"variant_sequence": Variant(notation="V227K").apply_to(wt),
         "computed_score": 0.62},
        {"variant_sequence": Variant(notation="V183L:F184I").apply_to(wt),
         "computed_score": 0.48},
        {"variant_sequence": Variant(notation="V183I").apply_to(wt),
         "computed_score": 0.22},
    ]
    return strategy, wt, pool, evidence


def _select(selection_rule, *, ablate=(), quality_perturbation=None, batch_size=4):
    """One real selection round, returning (committed notations, applied rule).

    `quality_perturbation`, when given, is a factor applied to every Expected Improvement
    the engine computes. EI is the sole input to `acquisition_scores` under both acquisition
    branches and it propagates through the Pareto multiplier into `quality_scores`, so
    negating it reverses the ranking the stack would commit. Nothing else about the call
    changes: same seed, same pool, same evidence, same surrogate.
    """
    strategy, wt, pool, evidence = _trpb_selection_inputs()
    real_ei = math_worker.compute_expected_improvement

    def perturbed(mean, std, f_best, xi=0.01):
        return real_ei(mean, std, f_best, xi) * quality_perturbation

    patch = (mock.patch.object(math_worker, "compute_expected_improvement", perturbed)
             if quality_perturbation is not None
             else contextlib_nullcontext())
    with patch:
        batch = math_worker.select_next_micro_batch(
            strategy=strategy, wt_sequence=wt, evidence=evidence,
            batch_size=batch_size, candidate_pool=pool,
            ablate=set(ablate) or None, random_seed=1234,
            selection_rule=selection_rule,
        )
    return [v.notation for v in batch.variants], batch.selection_rule_applied


def _live_selector(selection_rule):
    """A real `PanoptesSelector` on TrpB, with a revealed state it can plan from.

    Returns (selector, revealed, strategy). The caller patches `_strategy` to return that
    strategy, which is the only part of the production path that needs network-free
    metadata resolution it has already done once; everything else - generation, the
    GFlowNet sampler, the veto, the surrogate, the batch branch and the read-back - runs
    for real. `view` is unused on this path because TrpB is neither pool-restricted nor
    run with --production-candidate-source enumerate, so `propose` never calls
    `view.untested`; passing None would fail loudly if that ever changed.
    """
    from controller.prompts import ResidueCluster, SearchStrategy

    bench = _benchmark_module()
    library = bench.LIBRARIES["trpb"]
    metadata = bench.resolve_library_metadata(library)
    selector = bench.PanoptesSelector(
        0, metadata.sequence, metadata, False, None, library,
        selection_rule=selection_rule)
    revealed = bench.RevealedState()
    revealed.measured.update({
        "V183L": 0.40, "F184I": 0.55, "S228G": 0.31, "V227K": 0.62,
        "V183L:F184I": 0.48, "V183I": 0.22,
    })
    strategy = SearchStrategy(
        target_protein="TrpB", target_accession="P0A879",
        hypothesis="exercise the production propose() path",
        target_clusters=[ResidueCluster(
            cluster_id="sites", positions=list(library.positions),
            rationale="the ordered library")],
        max_mutation_depth=3,
    )
    return selector, revealed, strategy


class contextlib_nullcontext:
    """A no-op context manager, so `_select` has one code path for both cases."""

    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


class TestTheRuleIsRegistered(unittest.TestCase):
    def test_the_engine_knows_the_rule(self) -> None:
        self.assertIn("random_from_ballot", math_worker.SELECTION_RULES)
        self.assertIn("ranked", math_worker.SELECTION_RULES)

    def test_ranked_is_the_default_so_no_existing_figure_moves(self) -> None:
        """The shipped engine must be what a caller that passes nothing gets.

        If this default ever flips, every number measured before the flip and every number
        measured after it are different arms under one name - the hazard HANDOFF §1 records
        for `DPP_DIVERSITY_ENABLED` and `TWO_STAGE_GP_ENABLED`.
        """
        import inspect

        signature = inspect.signature(math_worker.select_next_micro_batch)
        self.assertEqual(signature.parameters["selection_rule"].default, "ranked")

    def test_an_unrecognised_rule_is_refused_rather_than_defaulted(self) -> None:
        """A typo must not run the ranking stack inside a folder labelled as the coin.

        This is the same defect `parse_ablate_components` was hardened against: an
        unrecognised value that falls through produces an unablated run reported as an
        ablation arm.
        """
        strategy, wt, pool, evidence = _trpb_selection_inputs()
        with self.assertRaises(ValueError) as caught:
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=evidence,
                batch_size=2, candidate_pool=pool,
                selection_rule="random-from-ballot",  # hyphens, not underscores
            )
        self.assertIn("random-from-ballot", str(caught.exception))

    def test_the_benchmark_registers_it_as_a_method(self) -> None:
        bench = _benchmark_module()
        self.assertIn("panoptes_random_from_ballot", bench.ALL_METHODS)

    def test_it_is_not_promoted_into_the_defaults(self) -> None:
        """This changes no default. DEFAULT_METHODS is what an unqualified run measures."""
        bench = _benchmark_module()
        self.assertNotIn("panoptes_random_from_ballot", bench.DEFAULT_METHODS)

    def test_the_arm_carries_its_own_name(self) -> None:
        """The method registry and the --methods filter key on `name`.

        An arm renamed only in its description silently ran zero seeds once; an arm NOT
        renamed would be tabulated as `panoptes_default` and the two could not be paired
        apart.
        """
        bench = _benchmark_module()
        library = bench.LIBRARIES["trpb"]
        wt = bench.resolve_library_metadata(library).sequence
        selector = bench.PanoptesSelector(
            0, wt, bench.resolve_library_metadata(library), False, None, library,
            selection_rule="random_from_ballot")
        self.assertEqual(selector.name, "panoptes_random_from_ballot")
        default = bench.PanoptesSelector(
            0, wt, bench.resolve_library_metadata(library), False, None, library)
        self.assertEqual(default.name, "panoptes_default")


class TestTheCoinIgnoresTheRankingStack(unittest.TestCase):
    """The decisive property: the committed batch is a function of the ballot and the RNG.

    Each test here is paired with the control arm under the identical perturbation. If the
    control stops responding, the perturbation no longer reaches the ranking stack and these
    tests prove nothing - so the control assertion is a test in its own right, not a
    comment.
    """

    def test_the_ranked_arm_DOES_respond_to_the_perturbation(self) -> None:
        """NON-VACUITY GUARD. Reversing EI must change what the ranking stack commits.

        If this fails, `quality_perturbation` is not reaching the scores and the invariance
        test below cannot distinguish a working coin from a dead perturbation.
        """
        baseline, rule = _select("ranked")
        self.assertEqual(rule, "top_n_quality")
        reversed_scores, _ = _select("ranked", quality_perturbation=-1.0)
        self.assertNotEqual(
            baseline, reversed_scores,
            "reversing every Expected Improvement did not change the ranked arm's batch, "
            "so the perturbation does not reach the ranking stack and the invariance test "
            "below is vacuous")

    def test_the_coin_is_invariant_to_the_same_perturbation(self) -> None:
        baseline, rule = _select("random_from_ballot")
        self.assertEqual(rule, "random_from_ballot")
        perturbed, _ = _select("random_from_ballot", quality_perturbation=-1.0)
        self.assertEqual(
            baseline, perturbed,
            "the coin's batch moved when the acquisition scores were reversed, so the "
            "ranking stack still reaches the selection")

    def test_the_coin_is_invariant_to_the_pareto_ranks_too(self) -> None:
        """Pareto enters as a MULTIPLIER on acquisition value, by a different route from EI.

        Checked separately because an implementation that drew uniformly from the
        Pareto-filtered subset instead of from the whole ballot would pass the EI test and
        fail this one.
        """
        strategy, wt, pool, evidence = _trpb_selection_inputs()
        real_ranks = math_worker.compute_non_dominated_pareto_ranks

        def forced(objectives: np.ndarray) -> np.ndarray:
            """Ranks that penalise by ballot position, hard.

            Not an inversion of the real objectives: with this pool the real ranks and
            their inversion gave the ranked arm the same top batch, so the non-vacuity
            guard below failed and the test proved nothing. A rank vector of 0..n-1 makes
            the Pareto multiplier 1/(1 + 0.3*i), a factor of ~23 across the ballot, which
            is large enough to reorder any ranking that reads it at all.
            """
            return np.arange(len(np.asarray(objectives)), dtype=int)

        def run(rule, patched):
            target = forced if patched else real_ranks
            with mock.patch.object(
                math_worker, "compute_non_dominated_pareto_ranks", target
            ):
                batch = math_worker.select_next_micro_batch(
                    strategy=strategy, wt_sequence=wt, evidence=evidence,
                    batch_size=4, candidate_pool=pool, random_seed=1234,
                    selection_rule=rule,
                )
            return [v.notation for v in batch.variants]

        # NON-VACUITY GUARD first: the ranked arm must move.
        self.assertNotEqual(
            run("ranked", False), run("ranked", True),
            "inverting the Pareto objectives did not change the ranked arm's batch, so "
            "this perturbation does not reach the stack")
        self.assertEqual(
            run("random_from_ballot", False), run("random_from_ballot", True),
            "the coin's batch moved when the Pareto objectives were inverted")

    def test_the_coin_does_not_simply_reproduce_the_ranked_batch(self) -> None:
        """A coin that happens to agree with the ranking stack would pass every invariance
        test above while measuring nothing. On a ballot of 13 with a batch of 4 the chance
        of exact agreement is small but not zero, so this pins the implementations apart
        rather than the outcome: it is the one assertion here that could in principle be a
        coincidence, and it is kept because a silent aliasing of the two rules is the more
        likely failure."""
        self.assertNotEqual(_select("ranked")[0], _select("random_from_ballot")[0])


class TestTheCoinDrawsFromTheEnginesOwnBallot(unittest.TestCase):
    def test_every_committed_variant_is_a_ballot_member(self) -> None:
        """The arm must select FROM the engine's pool, not propose a fresh one.

        Read off the candidate trace, which `select_next_micro_batch` hands the post-veto
        ballot and the committed indices in one call - so this checks the engine's own
        account of its ballot, not a reconstruction.
        """
        strategy, wt, pool, evidence = _trpb_selection_inputs()

        class _Sink:
            def __init__(self):
                self.rounds = []

            def record_round(self, **kw):
                self.rounds.append(kw)

        sink = _Sink()
        batch = math_worker.select_next_micro_batch(
            strategy=strategy, wt_sequence=wt, evidence=evidence,
            batch_size=4, candidate_pool=pool, random_seed=1234,
            selection_rule="random_from_ballot", candidate_trace=sink,
        )
        self.assertEqual(len(sink.rounds), 1, "the trace sink was never called, so this "
                                              "test is reading nothing")
        recorded = sink.rounds[0]
        ballot = list(recorded["candidates"])
        self.assertGreater(len(ballot), 0)
        for variant in batch.variants:
            self.assertIn(variant.notation, ballot)
        self.assertEqual(
            sorted(recorded["selected_indices"]),
            sorted(ballot.index(v.notation) for v in batch.variants))

    def test_the_ranking_is_still_computed_and_recorded(self) -> None:
        """The stack is bypassed, NOT removed.

        The acquisition values and Pareto ranks must still reach the trace, because the
        whole point of the arm is to score the committed picks against the ranking they did
        not use. An implementation that short-circuited the scoring would make that
        comparison impossible and would silently also change the ballot, since the GFlowNet
        reward reads the same surrogate.
        """
        strategy, wt, pool, evidence = _trpb_selection_inputs()

        class _Sink:
            def __init__(self):
                self.rounds = []

            def record_round(self, **kw):
                self.rounds.append(kw)

        sink = _Sink()
        math_worker.select_next_micro_batch(
            strategy=strategy, wt_sequence=wt, evidence=evidence,
            batch_size=4, candidate_pool=pool, random_seed=1234,
            selection_rule="random_from_ballot", candidate_trace=sink,
        )
        recorded = sink.rounds[0]
        n = len(recorded["candidates"])
        for field in ("acquisition_score", "pareto_rank", "quality_score",
                      "predicted_fitness", "posterior_uncertainty"):
            self.assertEqual(len(recorded[field]), n, f"{field} was not recorded per "
                                                      f"ballot member")
        self.assertGreater(
            len(set(float(v) for v in recorded["acquisition_score"])), 1,
            "every acquisition score is identical, so the ranking stack did not actually "
            "run and the trace cannot be used to score the coin's picks against it")

    def test_it_fills_the_batch_with_distinct_substitution_sets(self) -> None:
        """Budget matching is the entire basis of the paired test.

        Distinct candidate NOTATIONS can collapse to one substitution set; the harness
        rejects the duplicate, and an arm whose duplicates are not deduplicated here spends
        a smaller budget than the arm it is paired against. That confounds 'ranking removed'
        with 'fewer measurements', which is the exact confound the top-n branch's own
        comment records.
        """
        notations, _ = _select("random_from_ballot", batch_size=6)
        self.assertEqual(len(notations), 6)
        identities = {frozenset(n.split(":")) for n in notations}
        self.assertEqual(len(identities), len(notations))

    def test_a_longer_request_extends_the_same_draw(self) -> None:
        """The harness's backfill asks the policy for a LONGER list and takes its next
        eligible picks, on the documented understanding that 'a longer ranked list repeats
        its own prefix'. A coin reseeded per request would violate that and the backfilled
        picks would be a second, different draw - a secondary rule the manifest does not
        name."""
        short, _ = _select("random_from_ballot", batch_size=3)
        long, _ = _select("random_from_ballot", batch_size=7)
        self.assertEqual(long[:len(short)], short)


class TestTheRecordReadsTheCodePathNotTheFlag(unittest.TestCase):
    """HANDOFF rule 2.12: a manifest's component state is only as honest as what it reads.

    `component_states()` once derived the record from `--ablate-components` alone, so `dpp`
    and `dkl` - gated by module globals - were recorded as ENABLED in every run where they
    were off, and the HER2 null control's t=-8.29 movement could not be attributed from the
    artifacts because both manifests claimed the same state.
    """

    def test_the_engine_reports_which_branch_ran(self) -> None:
        self.assertEqual(_select("random_from_ballot")[1], "random_from_ballot")
        self.assertEqual(_select("ranked")[1], "top_n_quality")

    def test_the_reported_value_is_not_an_echo_of_the_argument(self) -> None:
        """NON-VACUITY GUARD for the whole class.

        If `selection_rule_applied` merely echoed `selection_rule`, it could never
        distinguish the two shipped ranking branches - and 'ranked' is what the caller
        asks for while `top_n_quality` is what runs, because DPP_DIVERSITY_ENABLED
        defaults False. That difference is the proof the value is read from the branch.
        """
        requested, applied = "ranked", _select("ranked")[1]
        self.assertNotEqual(
            requested, applied,
            "the engine echoed the requested rule back, so the record is reading the "
            "argument and not the code path")
        self.assertIn(applied, ("kdpp", "top_n_quality"))

    def test_the_kdpp_branch_is_named_separately(self) -> None:
        """With k-DPP restored the SAME request must report a different applied rule.

        Set on the module global, which is where the behaviour is gated - the same place
        `component_states` had to learn to read.
        """
        original = math_worker.DPP_DIVERSITY_ENABLED
        try:
            math_worker.DPP_DIVERSITY_ENABLED = True
            self.assertEqual(_select("ranked")[1], "kdpp")
        finally:
            math_worker.DPP_DIVERSITY_ENABLED = original
        self.assertEqual(_select("ranked")[1], "top_n_quality")

    def test_the_selector_refuses_a_round_whose_rule_cannot_be_verified(self) -> None:
        """If the switch stops reaching the call site, the arm must FAIL, not score.

        Simulated by making the engine report the other rule, which is what a regression
        that dropped the keyword would produce. An arm that scored a round here would
        produce a run folder labelled as the coin arm whose batches came from the ranking
        stack - the 2.12 defect, in the direction that inverts a conclusion.
        """
        bench = _benchmark_module()
        selector, revealed, strategy = _live_selector("random_from_ballot")
        real = math_worker.select_next_micro_batch

        def mislabelled(*args, **kwargs):
            batch = real(*args, **kwargs)
            batch.selection_rule_applied = "top_n_quality"
            return batch

        with mock.patch.object(type(selector), "_strategy",
                               lambda self, ri, rev: strategy):
            with mock.patch.object(math_worker, "select_next_micro_batch", mislabelled):
                with self.assertRaises(bench.BenchmarkError) as caught:
                    selector.propose(1, 4, None, revealed)
        self.assertIn("2.12", str(caught.exception))

    def test_the_arm_records_what_it_did_in_its_notes(self) -> None:
        """`results.notes` in the manifest is built from the selector's notes, so the
        engine-reported rule lands in the committed artifact rather than only in a log.

        This runs the REAL generative path - no supplied pool - because that is the path
        every CreiLOV, TrpB and GB1 campaign takes, and a read-back that only works in the
        pool-restricted regime would be verified on a path the campaigns never use.
        """
        selector, revealed, strategy = _live_selector("random_from_ballot")
        with mock.patch.object(type(selector), "_strategy",
                               lambda self, ri, rev: strategy):
            proposals = selector.propose(1, 4, None, revealed)
        self.assertEqual(len(proposals), 4)
        self.assertEqual(selector.selection_rule_observed, "random_from_ballot")
        self.assertTrue(
            any("SELECTION_RULE_APPLIED=random_from_ballot" in n for n in selector.notes),
            f"the engine-reported rule is not in the arm's notes: {selector.notes}")

    def test_the_default_arm_on_the_same_path_reports_the_ranking_branch(self) -> None:
        """NON-VACUITY GUARD for the test above: the generative path must distinguish the
        two arms, or the read-back is only echoing a constant."""
        selector, revealed, strategy = _live_selector("ranked")
        with mock.patch.object(type(selector), "_strategy",
                               lambda self, ri, rev: strategy):
            selector.propose(1, 4, None, revealed)
        self.assertEqual(selector.selection_rule_observed, "top_n_quality")


class TestTheWildTypeAnchorIsAblatable(unittest.TestCase):
    """The synthetic `y = 0.0` wild-type row, which every CreiLOV figure is measured
    through: 15.68 label standard deviations from the centre of that landscape's real
    labels, asserted at variance 1e-4 against a real measurement's 1e-2.

    It is ablatable so that a random-from-ballot win can be separated into 'the ranking
    stack was harmful' and 'the surrogate was poisoned', because a coin is immune to a
    poisoned surrogate and the default arm is not.
    """

    def test_it_is_in_the_ablation_vocabulary(self) -> None:
        self.assertIn("wt_anchor", math_worker.ABLATABLE_COMPONENTS)

    def test_the_benchmark_accepts_the_name(self) -> None:
        bench = _benchmark_module()
        self.assertIn("wt_anchor", bench.ABLATABLE_COMPONENTS)
        self.assertEqual(bench.parse_ablate_components("wt_anchor"), ("wt_anchor",))

    def _training_rows(self, ablate) -> int:
        """How many rows the surrogate was fitted on.

        Spied at `GaussianProcessRegressor.fit`, which is where the row count becomes a
        fact. Patching the featuriser instead would have counted candidate featurisations
        too and the two arms would have differed for the wrong reason.
        """
        from sklearn.gaussian_process import GaussianProcessRegressor

        strategy, wt, pool, evidence = _trpb_selection_inputs()
        seen = []
        real_fit = GaussianProcessRegressor.fit

        def spy(self, X, y, **kw):
            seen.append(np.asarray(X).shape[0])
            return real_fit(self, X, y, **kw)

        with mock.patch.object(GaussianProcessRegressor, "fit", spy):
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=evidence,
                batch_size=2, candidate_pool=pool, random_seed=7,
                ablate=set(ablate) or None,
            )
        self.assertTrue(seen, "the surrogate was never fitted, so this measures nothing")
        return seen[0]

    def test_ablating_it_removes_exactly_one_training_row(self) -> None:
        """A difference of ONE, not a literal count.

        Asserting '7 rows' would pin the fixture's evidence length into the test, which is
        the mistake `test_spatial_filter_ablation` records twice over. The mechanism is
        'one synthetic row, present or absent'.
        """
        with_anchor = self._training_rows(())
        without = self._training_rows(("wt_anchor",))
        self.assertEqual(
            with_anchor - without, 1,
            f"ablating wt_anchor changed the training set by "
            f"{with_anchor - without} rows, not 1")

    def test_the_default_still_carries_the_row(self) -> None:
        """NON-VACUITY GUARD: the row must be there to begin with, or the test above
        passes on an engine that never injected it."""
        strategy, wt, pool, evidence = _trpb_selection_inputs()
        self.assertEqual(self._training_rows(()), len(evidence) + 1)

    def test_no_evidence_and_no_anchor_is_refused_by_name(self) -> None:
        """With the anchor gone the only unconditional row is gone, so a call with no
        evidence has nothing to fit. It must say that, rather than surfacing as an
        IndexError on a feature-matrix shape."""
        strategy, wt, pool, _ = _trpb_selection_inputs()
        with self.assertRaises(ValueError) as caught:
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=[],
                batch_size=2, candidate_pool=pool, ablate={"wt_anchor"},
            )
        self.assertIn("wt_anchor", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
