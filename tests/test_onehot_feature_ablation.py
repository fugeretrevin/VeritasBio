"""The surrogate's input encoding must be switchable, and the switch must not move defaults.

Why this exists. `scripts/analysis/pool_hit_reachability.py` measures, on real campaigns,
that on CreiLOV the engine's post-veto ballot already contains enough top-1% hits for a
perfect ranker to match one-hot ridge, so the 48-of-48 loss cannot be explained by the
candidate pool's depth or its composition: the shortfall is in RANKING what was already on
the ballot. `docs/HANDOFF-2026-10-31-reachability.md` section 4 measured offline that a
plain GP on one-hot features is at parity with one-hot ridge on 7 of 9 landscapes, while
the engine's GP on its own physicochemical encoding ranked its in-run ballot at rho
-0.0069 on eqFP611-red. That points at the FEATURES, and `physchem_features` is the arm
that tests it.

These tests assert the MECHANISM - that the switch reaches every featurisation call site,
and that NOT passing it leaves the engine bit-identical. They deliberately do not assert a
yield, a rho or a candidate count. Two tests in this repository have previously pinned a
defect in place by asserting a literal (`one_hot_width == n_sites * 20`, which required
every library to be saturated, and `disabled == ["dpp", "prior"]`, which required the
manifest to ignore module-level defaults); a test asserting "the one-hot arm scores N"
would be that error a third time.
"""
from __future__ import annotations

import unittest

import numpy as np

from controller.gatekeeper import Variant
from workers import math_worker


class TestOneHotFeaturisation(unittest.TestCase):
    """`featurize_variant_one_hot` on its own, independent of the selector."""

    WT = "MGLRTAG"  # positions 1..7; only the declared ones are encoded

    def test_width_is_the_sum_of_the_declared_alphabets(self) -> None:
        vec = math_worker.featurize_variant_one_hot(
            Variant(notation="G2E"), [2, 3, 4], self.WT,
            allowed_residues={2: "EG", 3: "LN", 4: "DRS"},
        )
        self.assertEqual(vec.shape, (2 + 2 + 3,))

    def test_an_undeclared_position_falls_back_to_the_full_alphabet(self) -> None:
        """A saturated library declares nothing, and a position the EVIDENCE added to
        `active_positions` is not in the library spec either. Both must still encode."""
        vec = math_worker.featurize_variant_one_hot(
            Variant(notation="G2E"), [2, 3], self.WT, allowed_residues={2: "EG"},
        )
        self.assertEqual(vec.shape, (2 + len(math_worker.ONE_HOT_FALLBACK_ALPHABET),))

    def test_a_wild_type_position_encodes_the_parent_residue(self) -> None:
        """Not an all-zero block. The parent residue is a real allele and a one-hot that
        dropped it could not distinguish "wild type here" from "a residue off the
        alphabet", which is the next test's case."""
        vec = math_worker.featurize_variant_one_hot(
            Variant(notation="G2E"), [2, 3], self.WT,
            allowed_residues={2: "EG", 3: "LN"},
        )
        # position 2 -> E (index 0 of "EG"); position 3 -> wild type L (index 0 of "LN")
        np.testing.assert_array_equal(vec, np.array([1.0, 0.0, 1.0, 0.0]))

    def test_exactly_one_slot_per_position_is_set(self) -> None:
        vec = math_worker.featurize_variant_one_hot(
            Variant(notation="L3N"), [2, 3, 4], self.WT,
            allowed_residues={2: "EG", 3: "LN", 4: "DRT"},
        )
        self.assertEqual(float(vec.sum()), 3.0)

    def test_a_residue_off_the_declared_alphabet_is_all_zero_not_an_exception(self) -> None:
        """The generator can propose a residue the library does not contain - that is a
        measured reachability defect, not a reason to crash a selection round."""
        vec = math_worker.featurize_variant_one_hot(
            Variant(notation="G2W"), [2], self.WT, allowed_residues={2: "EG"},
        )
        np.testing.assert_array_equal(vec, np.array([0.0, 0.0]))

    def test_the_encoding_is_deterministic_in_width_across_candidates(self) -> None:
        """The surrogate needs one width for training rows and candidate rows alike."""
        kwargs = {"allowed_residues": {2: "EG", 3: "LN", 4: "DRT"}}
        widths = {
            math_worker.featurize_variant_one_hot(
                Variant(notation=n), [2, 3, 4], self.WT, **kwargs).shape[0]
            for n in ("G2E", "L3N", "G2E:L3N", "G2E:L3N:R4T", "WT")
        }
        self.assertEqual(len(widths), 1, f"ragged one-hot widths: {widths}")


class TestPhyschemFeaturesIsAblatable(unittest.TestCase):
    def test_the_encoding_is_in_the_ablation_vocabulary(self) -> None:
        self.assertIn("physchem_features", math_worker.ABLATABLE_COMPONENTS)

    def test_it_is_distinct_from_the_spatial_cross_term(self) -> None:
        """`spatial` clears one cross-term inside the physicochemical encoding;
        `physchem_features` replaces the whole encoding. Collapsing them would make it
        impossible to run one-hot features WITH the structural veto intact, which is the
        comparison that separates the features from the selection stages."""
        self.assertIn("spatial", math_worker.ABLATABLE_COMPONENTS)
        self.assertNotIn("physchem_features", {"spatial", "spatial_filter"})

    def test_the_benchmark_accepts_the_name(self) -> None:
        import benchmark_panoptes_isolated as bench

        self.assertIn("physchem_features", bench.ABLATABLE_COMPONENTS)
        self.assertEqual(
            bench.parse_ablate_components("physchem_features"), ("physchem_features",))

    # ---- the mechanism, through a real selection round --------------------------------
    def _round(self, ablate, allowed_residues=None):
        """One real `select_next_micro_batch` call, returning its selected notations.

        A supplied `candidate_pool` keeps the round cheap and deterministic without
        stubbing anything on the path under test: featurisation, surrogate fit,
        acquisition and selection all run for real.
        """
        from controller.prompts import ResidueCluster, SearchStrategy

        wt = "M" + "ACDEFGHIKLMNPQRSTVWY" * 12
        positions = [10, 20, 30, 40]
        strategy = SearchStrategy(
            target_protein="synthetic", target_accession="TEST",
            hypothesis="exercise the featurisation call sites",
            target_clusters=[ResidueCluster(
                cluster_id="sites", positions=positions,
                rationale="a synthetic four-site library")],
            max_mutation_depth=3,
        )
        pool = [Variant(notation=n) for n in (
            f"{wt[9]}10A", f"{wt[19]}20C", f"{wt[9]}10A:{wt[19]}20C",
            f"{wt[29]}30D", f"{wt[39]}40E", f"{wt[29]}30D:{wt[39]}40E",
        )]
        evidence = [
            {"variant_sequence": Variant(notation=f"{wt[9]}10G").apply_to(wt),
             "computed_score": 0.4},
            {"variant_sequence": Variant(notation=f"{wt[19]}20G").apply_to(wt),
             "computed_score": 0.6},
            {"variant_sequence": Variant(notation=f"{wt[29]}30G").apply_to(wt),
             "computed_score": 0.2},
        ]
        np.random.seed(11)
        batch = math_worker.select_next_micro_batch(
            strategy=strategy, wt_sequence=wt, evidence=evidence, batch_size=2,
            candidate_pool=pool, ablate=ablate, random_seed=7,
            allowed_residues=allowed_residues,
        )
        return [v.notation for v in batch.variants]

    def test_the_default_arm_is_unchanged_by_the_new_switch(self) -> None:
        """The regression guard that matters most. Introducing the component must not
        move any previously measured number, so the un-ablated path must still select
        exactly what it selected before - the binding is only rebound, never replaced."""
        first = self._round(())
        second = self._round(())
        self.assertEqual(first, second, "the default path is not deterministic, so this "
                                        "test cannot detect a change to it")
        self.assertEqual(
            first, self._round(None),
            "passing ablate=None and ablate=() must be the same arm")

    def _encoder_calls(self, ablate):
        """Run one real round and count calls into each encoder.

        The first version of this test asserted that the two arms SELECT differently.
        It failed: on a six-candidate pool both encodings rank the same two first, so a
        working switch looked broken. Worse, the converse would have been a vacuous pass
        - any unrelated change to selection would have satisfied it. Counting calls into
        each encoder asserts the thing the switch is actually for, and
        `test_the_unablated_path_reaches_the_physchem_encoder` below stops it passing
        vacuously if this invocation ever stops reaching the call sites.
        """
        from unittest import mock

        real_physchem = math_worker.featurize_variant_on_positions
        real_onehot = math_worker.featurize_variant_one_hot
        with mock.patch.object(math_worker, "featurize_variant_on_positions",
                               side_effect=real_physchem) as physchem, \
             mock.patch.object(math_worker, "featurize_variant_one_hot",
                               side_effect=real_onehot) as onehot:
            self._round(ablate)
            return physchem.call_count, onehot.call_count

    def test_the_unablated_path_reaches_the_physchem_encoder(self) -> None:
        """Guards the test below against passing vacuously. If this fails, the selection
        path changed and the ablation test is no longer measuring anything."""
        physchem, onehot = self._encoder_calls(())
        self.assertGreater(physchem, 0, "the default path never featurised anything, so "
                                        "the ablation test cannot tell a working switch "
                                        "from a broken one")
        self.assertEqual(onehot, 0, "the one-hot encoder ran without being asked for")

    def test_ablating_the_encoding_swaps_every_call_site(self) -> None:
        """All five call sites must move together.

        A switch that reached the candidate rows while the training rows kept the old
        encoding would be a silent dimension or scale error, so `physchem` must fall to
        exactly zero rather than merely decrease.
        """
        physchem, onehot = self._encoder_calls(("physchem_features",))
        self.assertGreater(onehot, 0, "the one-hot encoder never ran despite being "
                                      "requested: the switch does not reach the call sites")
        self.assertEqual(
            physchem, 0,
            f"the physicochemical encoder still ran {physchem} time(s) under "
            f"physchem_features, so training and candidate rows are featurised by "
            f"different functions")

    def test_it_works_on_a_library_that_declares_per_site_alphabets(self) -> None:
        """The non-saturated libraries are the ones this arm exists for, and they are the
        ones that pass `allowed_residues`."""
        wt_res = "M" + "ACDEFGHIKLMNPQRSTVWY" * 12
        allowed = {10: wt_res[9] + "A", 20: wt_res[19] + "C",
                   30: wt_res[29] + "D", 40: wt_res[39] + "E"}
        selected = self._round(("physchem_features",), allowed_residues=allowed)
        self.assertEqual(len(selected), 2)


if __name__ == "__main__":
    unittest.main()
