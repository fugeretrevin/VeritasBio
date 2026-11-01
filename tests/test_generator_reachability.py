"""Reachability of the candidate generator: can it propose the answer at all?

Every test here was written against a defect measured in
`docs/day1-reachability.md`, on a 48-seed CR9114 run in which the production arm
proposed 757,829 candidates, had 25,712 rejected as `not_in_eligible_universe` with
no other reason appearing once, accepted 32, and found 0 of 651 top-1% hits.

The three defects, each with a test below:

  1. ALPHABET. `generate_candidate_variants` hardcodes a 19-residue alphabet
     (math_worker.py:1118). CR9114 is a BINARY library: each of its 16 positions
     admits exactly 2 residues, the germline one and the matured one. So ~94.7% of
     proposed singles do not exist in the library and the ineligibility compounds
     with depth. This is what starved the arm.

  2. DEPTH. `max_mutation_depth` is capped at 4 by a pydantic validator
     (controller/prompts.py:100, `le=4`). With `max_depth = max(..., 3)` at
     math_worker.py:1178 the effective depth is always in {3, 4}, for every library
     at every budget. CR9114's optimum is at depth 7 and 647 of its 651 top-1% hits
     are deeper than 4, so the ceiling makes them unreachable by construction.

  3. TASK SYMMETRY. `RandomSelector` and `RidgeSelector` both call
     `view.untested(revealed)` unconditionally, so they always rank the enumerated
     eligible pool. The production arm only receives that pool when the library
     declares `pool_restricted`, which exactly one library does - HER2, the only one
     that is NOT enumerable. The harness's own manifest text calls ranking an
     enumerated pool "strictly easier than generating under constraints", and that
     easier task was being handed to the baselines and withheld from the arm under
     test on eight of nine libraries.
"""
from __future__ import annotations

import unittest

from controller.prompts import ResidueCluster, SearchStrategy
from workers.math_worker import generate_candidate_variants


def _binary_strategy(positions, depth):
    """A strategy over `positions` in one cluster, as the benchmark planner builds it."""
    return SearchStrategy(
        target_protein="CR9114",
        target_accession="TEST",
        hypothesis="reachability fixture",
        target_clusters=[ResidueCluster(
            cluster_id="vh", positions=list(positions),
            rationale="all somatic positions are in scope",
        )],
        forbidden_positions=[],
        max_mutation_depth=depth,
    )


# A 16-site binary library in the shape of CR9114: parent residue at each site plus
# exactly one designated alternate allele. Positions and parent residues are the real
# ones from the LibrarySpec; the alternates are stand-ins, since the test only needs
# the ALPHABET SHAPE to be binary, not the biology to be right.
CR9114_POSITIONS = (29, 30, 31, 52, 57, 58, 59, 71, 74, 75, 76, 77, 84, 87, 95, 106)
CR9114_PARENT = "FSSITANTKSTSSRYY"
CR9114_ALTERNATE = "SNNVSGTSIFSNTHSS"

# Wild-type sequence long enough to index every position, carrying the parent residue
# at each variable site. Position numbering is 1-based.
_WT = ["A"] * 120
for _p, _aa in zip(CR9114_POSITIONS, CR9114_PARENT):
    _WT[_p - 1] = _aa
WT_SEQUENCE = "".join(_WT)

ALLOWED = {
    pos: (parent, alternate)
    for pos, parent, alternate in zip(CR9114_POSITIONS, CR9114_PARENT, CR9114_ALTERNATE)
}


class TestAlphabetIsRespected(unittest.TestCase):
    """Defect 1: a per-site allele restriction must be honoured by every construct."""

    def test_every_proposal_uses_only_declared_alleles(self):
        strategy = _binary_strategy(CR9114_POSITIONS, depth=4)
        candidates = generate_candidate_variants(
            strategy=strategy,
            wt_sequence=WT_SEQUENCE,
            max_candidates=2500,
            allowed_residues=ALLOWED,
        )
        self.assertGreater(len(candidates), 0, "generator emitted nothing")

        offenders = []
        for variant in candidates:
            for mutation in variant.mutations:
                permitted = ALLOWED.get(mutation.position)
                if permitted is None:
                    offenders.append((variant.notation, mutation.position, "position not in library"))
                elif mutation.mut_aa not in permitted:
                    offenders.append((variant.notation, mutation.position, mutation.mut_aa))
        self.assertEqual(
            offenders[:20], [],
            f"{len(offenders)} mutations used a residue the library does not contain. "
            f"On a binary library this is what produces `not_in_eligible_universe`.",
        )

    def test_multi_mutant_constructs_also_respect_the_alphabet(self):
        """The intra/inter-cluster constructs truncate the alphabet ALPHABETICALLY.

        `target_aas[:10]` and `target_aas[:8]` (math_worker.py:1247, 1275, 1278) keep
        `A C D E F G H I [K L]` and drop `M N Q R S T V W Y`. That cut is alphabetical,
        not chemical, and on a restricted library it discards most legal alleles. So
        this asserts specifically on candidates of depth >= 2.
        """
        strategy = _binary_strategy(CR9114_POSITIONS, depth=4)
        candidates = generate_candidate_variants(
            strategy=strategy,
            wt_sequence=WT_SEQUENCE,
            max_candidates=2500,
            allowed_residues=ALLOWED,
        )
        multi = [v for v in candidates if len(v.mutations) >= 2]
        self.assertGreater(len(multi), 0, "no multi-mutant candidates were generated at all")
        for variant in multi:
            for mutation in variant.mutations:
                self.assertIn(
                    mutation.mut_aa, ALLOWED[mutation.position],
                    f"{variant.notation} mutates {mutation.position} to "
                    f"{mutation.mut_aa}, which is not a declared allele",
                )

    def test_unrestricted_callers_are_unaffected(self):
        """Omitting the restriction must leave the existing behaviour alone.

        TrpB and GB1 are saturated 4-site libraries and their published numbers come
        from this path, so the default must not move.
        """
        strategy = _binary_strategy((183, 184, 227, 228), depth=4)
        wt = ["A"] * 300
        for p, aa in zip((183, 184, 227, 228), "VFVS"):
            wt[p - 1] = aa
        candidates = generate_candidate_variants(
            strategy=strategy, wt_sequence="".join(wt), max_candidates=2500,
        )
        self.assertGreater(len(candidates), 100)
        residues = {m.mut_aa for v in candidates for m in v.mutations}
        self.assertGreater(
            len(residues), 10,
            "with no restriction the generator should still use a wide alphabet",
        )


class TestDepthCeiling(unittest.TestCase):
    """Defect 2: the schema's `le=4` makes CR9114's optimum unreachable."""

    def test_strategy_accepts_a_depth_matching_a_sixteen_site_library(self):
        """CR9114 has 16 sites. A strategy must be able to say so.

        This is the single most binding constraint found: it is a validator, not an
        emergent property of the construct set, so no amount of budget or rounds
        lifts it.
        """
        strategy = _binary_strategy(CR9114_POSITIONS, depth=16)
        self.assertEqual(strategy.max_mutation_depth, 16)

    def test_generator_can_propose_beyond_depth_four(self):
        """With depth headroom and deep leads, proposals must exceed depth 4.

        647 of CR9114's 651 top-1% hits are deeper than 4, and its optimum is at
        depth 7, so a generator that cannot pass 4 cannot find them.
        """
        from controller.gatekeeper import Variant

        strategy = _binary_strategy(CR9114_POSITIONS, depth=16)
        # A depth-6 lead, as a campaign measuring deep variants would supply. The
        # initial random sample on CR9114 contains a depth>=9 variant with
        # probability 1.000, so deep leads are the normal case, not an edge case.
        lead_notation = ":".join(
            f"{CR9114_PARENT[i]}{CR9114_POSITIONS[i]}{CR9114_ALTERNATE[i]}"
            for i in range(6)
        )
        lead = Variant.from_sequence  # presence check; construction below uses notation
        self.assertTrue(callable(lead))
        candidates = generate_candidate_variants(
            strategy=strategy,
            wt_sequence=WT_SEQUENCE,
            max_candidates=2500,
            lead_variants=[Variant(notation=lead_notation)],
            allowed_residues=ALLOWED,
        )
        depths = {len(v.mutations) for v in candidates}
        self.assertTrue(
            any(d > 4 for d in depths),
            f"no candidate deeper than 4 was proposed; depths present: {sorted(depths)}",
        )


class TestPoolIsLargeEnoughToSpendTheBudget(unittest.TestCase):
    """Defect 3/5: the position slices were written for a 19-residue alphabet.

    `available_pos[:6]`, `combinations(available_pos[:4], 2)`, `cluster_allowed[:4]` and
    `c1_pos[:3]` bound the number of CANDIDATES a stage emits, not the number of
    positions: at 18 usable residues per site, 6 positions is ~108 candidates. On a
    binary library the same 6 positions yield 6, so after the alphabet fix the CR9114
    pool was 37-115 candidates against a 125-measurement budget and batches still ran
    short - measured at 57-60 of 125 spent.
    """

    #: A campaign measures 25 initial + 4 x 25 = 125 variants. The pool is rebuilt each
    #: round, but candidates already measured are rejected, so a pool barely larger than
    #: the budget still starves. Twice the budget is the floor this asserts.
    BUDGET = 125

    def test_binary_library_pool_can_cover_a_full_campaign(self):
        from controller.gatekeeper import Variant

        leads = [
            Variant(notation=f"{CR9114_PARENT[i]}{CR9114_POSITIONS[i]}{CR9114_ALTERNATE[i]}")
            for i in range(8)
        ]
        candidates = generate_candidate_variants(
            strategy=_binary_strategy(CR9114_POSITIONS, depth=16),
            wt_sequence=WT_SEQUENCE,
            max_candidates=2500,
            lead_variants=leads,
            allowed_residues=ALLOWED,
        )
        self.assertGreater(
            len(candidates), 2 * self.BUDGET,
            f"only {len(candidates)} distinct candidates for a {self.BUDGET}-measurement "
            f"campaign on a 16-site binary library; the batch will run short",
        )

    def test_saturated_library_pool_is_unchanged(self):
        """TrpB and GB1 are saturated, and their published numbers come from this path.

        Widening the position budget must be a no-op when every site already offers a
        full alphabet, or the fix silently moves the development-landscape results.
        """
        from controller.gatekeeper import Variant

        positions = (183, 184, 227, 228)
        wt = ["A"] * 300
        for p, aa in zip(positions, "VFVS"):
            wt[p - 1] = aa
        leads = [Variant(notation=f"{'VFVS'[i]}{positions[i]}W") for i in range(2)]
        candidates = generate_candidate_variants(
            strategy=_binary_strategy(positions, depth=4),
            wt_sequence="".join(wt),
            max_candidates=2500,
            lead_variants=leads,
        )
        # The cap binds on a saturated 4-site library, which is the pre-existing
        # behaviour: the assertion is that it still binds at exactly the same place.
        self.assertEqual(len(candidates), 2500)


class TestEncodingMatchesTheDeclaredAlphabets(unittest.TestCase):
    """The one-hot encoder must use the declared alphabets, like `one_hot_width` says.

    `_one_hot_matrix` built `n_sites * 20` columns unconditionally. On the 16-site binary
    CR9114 library that is 320 columns of which 288 are identically zero, and it made
    `LibrarySpec.one_hot_width` - which the manifest reports as `one_hot_features` -
    disagree with the matrix actually handed to the baselines.
    """

    def test_matrix_width_equals_declared_width_for_every_library(self):
        import benchmark_panoptes_isolated as bench

        for key, spec in bench.LIBRARIES.items():
            with self.subTest(library=key):
                parent_combo = spec.parent_residues
                matrix = bench._one_hot_matrix([parent_combo], spec)
                self.assertEqual(
                    matrix.shape[1], spec.one_hot_width,
                    f"{key}: encoder produced {matrix.shape[1]} columns but "
                    f"one_hot_width declares {spec.one_hot_width}; the manifest reports "
                    f"the latter",
                )
                # Exactly one active column per site, and the row sums to n_sites.
                self.assertEqual(matrix.sum(), float(spec.n_sites))

    def test_saturated_libraries_keep_the_twenty_wide_encoding(self):
        import benchmark_panoptes_isolated as bench

        for key in ("trpb", "gb1"):
            with self.subTest(library=key):
                spec = bench.LIBRARIES[key]
                matrix = bench._one_hot_matrix([spec.parent_residues], spec)
                self.assertEqual(matrix.shape[1], spec.n_sites * 20)

    def test_a_residue_outside_the_declared_alphabet_is_refused(self):
        """Silently encoding it would put a library/spec disagreement into the features."""
        import benchmark_panoptes_isolated as bench

        spec = bench.LIBRARIES["cr9114_h1"]
        # Replace the first site's residue with one its declared alphabet excludes.
        alphabet = spec.alphabet_for(0)
        intruder = next(a for a in "ACDEFGHIKLMNPQRSTVWY" if a not in alphabet)
        combo = intruder + spec.parent_residues[1:]
        with self.assertRaises(bench.BenchmarkError):
            bench._one_hot_matrix([combo], spec)


class TestTaskSymmetry(unittest.TestCase):
    """Defect 3: baselines rank an enumerated pool; the arm under test did not."""

    #: Above this many genotypes, enumerating the design space is impractical and a
    #: generative proposer is the only option. Below it, the published norm is to score
    #: the whole library every round - MLDE ranks "all possible sequences" (Wu et al.
    #: 2019) and ALDE ranks "all sequences in the design space" (Yang et al. 2025).
    #: GB1 at 160,000 is routinely enumerated in that literature.
    ENUMERABLE_CEILING = 1_000_000

    def test_the_asymmetry_is_real_and_is_recorded_here(self):
        """Documents which libraries are affected, rather than demanding a config change.

        Flipping `pool_restricted` in the specs would silently redefine what the
        benchmark measures and invalidate comparison with every prior run, so the fix is
        a run-time switch (next test) and this test exists to keep the list honest. If a
        library leaves this list, the switch is no longer needed for it.
        """
        from benchmark_panoptes_isolated import LIBRARIES

        affected = []
        for key, spec in LIBRARIES.items():
            space = 1
            for slot in range(spec.n_sites):
                space *= len(spec.alphabet_for(slot))
            if space <= self.ENUMERABLE_CEILING and not spec.pool_restricted:
                affected.append(key)

        self.assertEqual(
            sorted(affected),
            ["cr6261_h1", "cr6261_h9", "cr9114_h1", "creilov",
             "eqfp611_blue", "eqfp611_red", "gb1", "trpb"],
            "the set of libraries whose production arm generates while its baselines "
            "rank an enumerated pool has changed; confirm that is intended and update "
            "docs/day1-reachability.md",
        )
        # HER2 is the one library that is NOT enumerable (1.8e11) and is the only one
        # already declared pool-restricted - exactly backwards from what its size implies.
        self.assertTrue(LIBRARIES["her2"].pool_restricted)

    def test_production_arm_can_be_put_on_the_same_task_as_its_baselines(self):
        """The switch must exist, validate its input, and declare itself in the run."""
        import inspect

        from benchmark_panoptes_isolated import PanoptesSelector, build_parser

        parameters = inspect.signature(PanoptesSelector.__init__).parameters
        self.assertIn("candidate_source", parameters)
        self.assertEqual(
            parameters["candidate_source"].default, "generate",
            "the default must stay 'generate' or every published number silently "
            "changes meaning",
        )

        action = next(
            a for a in build_parser()._actions
            if a.dest == "production_candidate_source"
        )
        self.assertEqual(tuple(action.choices), ("generate", "enumerate"))
        self.assertEqual(action.default, "generate")


if __name__ == "__main__":
    unittest.main()
