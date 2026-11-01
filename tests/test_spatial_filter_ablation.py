"""The spatial viability veto must be switchable, and the switch must actually switch it.

Why this exists. `filter_spatially_viable_variants` rejects any multi-mutant whose mutated
positions have mean pairwise C-alpha distance > 16 A or max > 24 A. That rule is monotone
in mutation depth, so it removes a growing share of every library as depth rises. Measured
over all nine libraries' full universes by `scripts/analysis/spatial_veto_profile.py`, the
fraction of each landscape's top-1% hits that can pass it at all runs from 62.6% (TrpB)
down to **0 of 82 on eqFP611-red** - a ceiling on achievable recall that applies before
any ranking happens.

Until now no flag could disable it. `--ablate-components spatial` clears the 3D cross-term
from the FEATURES and leaves the veto running, which is a different stage. So an
experiment asking "was the candidate pool the bottleneck?" could not distinguish that from
"was the filter the bottleneck?", and a null result would have been attributed to the
pool.

These tests assert the MECHANISM - that the switch reaches the call - and never a
candidate count or a pass rate. Two tests in this repository have previously pinned a
defect in place by asserting a literal (`one_hot_width == n_sites * 20`, which required
every library to be saturated, and `disabled == ["dpp", "prior"]`, which required the
manifest to ignore module-level defaults). A test that pins "the filter removes N
candidates" would be the same error a third time.
"""
from __future__ import annotations

import unittest
from unittest import mock

from workers import math_worker

_BENCH = None


def _benchmark_module():
    """Import `benchmark_panoptes_isolated` by path, once.

    It must be registered in `sys.modules` BEFORE `exec_module`, or dataclass field
    resolution fails with `'NoneType' object has no attribute '__dict__'` -
    `dataclasses` looks the defining module up by name while the class body executes.
    Cached because loading it is slow and every call here wants the same registry.
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


class TestSpatialFilterIsAblatable(unittest.TestCase):
    def test_the_veto_is_in_the_ablation_vocabulary(self) -> None:
        self.assertIn("spatial_filter", math_worker.ABLATABLE_COMPONENTS)

    def test_it_is_distinct_from_the_feature_cross_term(self) -> None:
        """`spatial` and `spatial_filter` are different stages and both must exist.

        `spatial` removes the 3D term from `featurize_variant_on_positions`; this one
        disables the viability veto. Collapsing them would make it impossible to ablate
        the veto while keeping the structural features, which is exactly the comparison
        the veto measurement calls for.
        """
        self.assertIn("spatial", math_worker.ABLATABLE_COMPONENTS)
        self.assertNotEqual("spatial", "spatial_filter")

    def test_the_benchmark_accepts_the_name(self) -> None:
        """The CLI vocabulary is imported from the engine, so a name the engine knows
        must be a name the benchmark accepts. `run_self_checks` asserts the two agree;
        this pins the specific name so a rename cannot pass silently."""
        import benchmark_panoptes_isolated as bench

        self.assertIn("spatial_filter", bench.ABLATABLE_COMPONENTS)
        self.assertEqual(
            bench.parse_ablate_components("spatial_filter"), ("spatial_filter",))

    def _call_count(self, ablate) -> int:
        """Run one real selection round and count calls into the veto.

        The invocation has to be realistic enough to REACH the filter. The first version
        of this test passed `evidence_records=[], strategy=None, wt_sequence=""`, which
        raised TypeError on an unknown keyword before the veto was ever reached - so both
        the ablated and unablated arms recorded zero calls and the test passed while
        proving nothing. Vacuous-pass is the failure mode this file warns about in its
        own docstring, and it occurred here on the first attempt. Hence
        `test_the_unablated_path_actually_reaches_the_veto` below, which fails if this
        invocation ever stops exercising the call site.
        """
        from controller.gatekeeper import Variant
        from controller.prompts import ResidueCluster, SearchStrategy
        from workers import spatial_filter as sf_module

        bench = _benchmark_module()
        library = bench.LIBRARIES["trpb"]
        wt = bench.resolve_library_metadata(library).sequence
        strategy = SearchStrategy(
            target_protein="TrpB", target_accession="P0A879",
            hypothesis="exercise the selection path",
            target_clusters=[ResidueCluster(
                cluster_id="sites", positions=list(library.positions),
                rationale="the ordered library")],
            max_mutation_depth=3,
        )
        pool = [Variant(notation=n) for n in
                ("V183L:F184I", "V183L:S228G", "F184I:V227K:S228G", "V183L")]
        evidence = [
            {"variant_sequence": Variant(notation="V183L").apply_to(wt),
             "computed_score": 0.4},
            {"variant_sequence": Variant(notation="F184I").apply_to(wt),
             "computed_score": 0.5},
        ]

        # Patch where the veto is DEFINED: the call site imports it inside the function
        # body, so patching the importing module would miss the binding and the spy would
        # record nothing regardless of the switch.
        with mock.patch.object(
            sf_module, "filter_spatially_viable_variants",
            side_effect=lambda notations, **kw: (list(notations), {}),
        ) as spy:
            math_worker.select_next_micro_batch(
                strategy=strategy, wt_sequence=wt, evidence=evidence,
                batch_size=2, candidate_pool=pool, ablate=ablate,
            )
            return spy.call_count

    def test_the_unablated_path_actually_reaches_the_veto(self) -> None:
        """Guards the test below against passing vacuously. If this fails, the selection
        path changed and the ablation test is no longer measuring anything."""
        self.assertGreater(
            self._call_count(()), 0,
            "the unablated path never reached the veto, so the ablation test below "
            "cannot distinguish a working switch from a broken one")

    def test_ablating_the_veto_skips_it(self) -> None:
        self.assertEqual(
            self._call_count(("spatial_filter",)), 0,
            "the veto ran despite being ablated: the switch does not reach the call site")


if __name__ == "__main__":
    unittest.main()
