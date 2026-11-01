"""The GP surrogate's feature space is a selectable ARM, and selecting it must move
exactly one thing - while NOT selecting it must move nothing at all.

What this pins. `build_gp_model` uses an ISOTROPIC Matern, so Euclidean distance is the
kernel's only input, while `featurize_variant_on_positions` appends a 6-dimensional global
block (`n_mut / 3.0` plus the five column sums of the per-site deltas) that is a
deterministic linear function of the per-site block plus the mutation count and was
measured to carry ~50% of total feature variance on every one of the nine benchmark
landscapes. It adds no information and only re-weights distance, so half of the GP's notion
of "nearby sequence" became "similar mutation burden". Over a FIXED design matrix one-hot
beats it on 8 of 9 landscapes.

But a better surrogate over a fixed matrix is not a better campaign - a sequential engine
chooses its own training distribution - so this is an ARM, not a fix:
`--ablate-components biophysical_features`, or the registered method
`panoptes_onehot_surrogate`. The engine default is unchanged, and the first class below is
the acceptance criterion for that: with the component absent, the GP is handed the same
array object it always was.

Two traps this file is written against, both of which have already fired here:

  * VACUOUS PASS. `tests/test_spatial_filter_ablation.py`'s first version raised TypeError
    before the stage under test was reached, so the ablated and unablated arms both
    recorded zero calls and the test passed proving nothing. Every ablation assertion below
    is therefore paired with a guard that fails if the unablated path stops reaching the
    stage.
  * PINNING A LITERAL. `one_hot_width == n_sites * 20` and `disabled == ["dpp", "prior"]`
    are both recorded in docs/HANDOFF.md as tests that turned a defect into an invariant.
    No width, count or component list is written as a number here: widths are derived from
    the declared alphabets, and the manifest's recorded feature space is checked against
    the matrix the GP was ACTUALLY handed rather than against the flag that asked for it
    (rule 2.12).
"""

from __future__ import annotations

import itertools
import unittest
from typing import Any, Dict, List, Optional, Sequence, Tuple
from unittest import mock

import numpy as np

from controller.gatekeeper import Variant
from controller.prompts import ResidueCluster, SearchStrategy
import workers.math_worker as mw


_BENCH = None


def _benchmark_module():
    """Import `benchmark_panoptes_isolated` by path, once.

    It must be registered in `sys.modules` BEFORE `exec_module`, or dataclass field
    resolution fails with `'NoneType' object has no attribute '__dict__'`.
    """
    global _BENCH
    if _BENCH is None:
        import importlib.util
        import sys

        if "benchmark_panoptes_isolated" in sys.modules:
            _BENCH = sys.modules["benchmark_panoptes_isolated"]
            return _BENCH
        spec = importlib.util.spec_from_file_location(
            "benchmark_panoptes_isolated", "benchmark_panoptes_isolated.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules["benchmark_panoptes_isolated"] = module
        spec.loader.exec_module(module)
        _BENCH = module
    return _BENCH


# --------------------------------------------------------------------------------------
# A small library with MIXED per-site alphabets, which is the case that distinguishes a
# declared-alphabet encoding from an unconditional 20-wide one, and the case where a
# truncating zip would silently succeed.
# --------------------------------------------------------------------------------------
POSITIONS: Tuple[int, ...] = (5, 9, 14)
ALLOWED_RESIDUES: Dict[int, str] = {5: "AV", 9: "VKR", 14: "DEFG"}
WT_SEQUENCE = "MKAVAHGYVSTGCDEWQLNP"  # position 5 -> 'A', 9 -> 'V', 14 -> 'D'
N_EVIDENCE = 6

#: The component whose ablation defines the arm. Imported from the engine rather than
#: written twice, so a rename cannot leave this file testing a name nothing reads.
COMPONENT = "biophysical_features"


def _strategy() -> SearchStrategy:
    return SearchStrategy(
        target_protein="synthetic",
        target_accession="SYN001",
        hypothesis="three variable sites with unequal alphabets",
        target_clusters=[
            ResidueCluster(
                cluster_id="sites",
                positions=list(POSITIONS),
                rationale="the declared library sites",
            ),
        ],
        forbidden_positions=[],
        max_mutation_depth=3,
    )


def _variant_for(genotype: str) -> Variant:
    """A Variant whose notation OMITS wild-type positions, as real notation does."""
    muts = [
        f"{WT_SEQUENCE[pos - 1]}{pos}{residue}"
        for pos, residue in zip(POSITIONS, genotype)
        if residue != WT_SEQUENCE[pos - 1]
    ]
    return Variant(notation=":".join(muts) if muts else "WT")


def _all_genotypes() -> List[str]:
    return ["".join(p) for p in itertools.product(*(ALLOWED_RESIDUES[q] for q in POSITIONS))]


def _evidence(genotypes: Sequence[str]) -> List[Dict[str, Any]]:
    out = []
    for index, genotype in enumerate(genotypes):
        sequence = list(WT_SEQUENCE)
        for pos, residue in zip(POSITIONS, genotype):
            sequence[pos - 1] = residue
        out.append({"variant_sequence": "".join(sequence),
                    "computed_score": 0.1 * (index + 1)})
    return out


class _SpyGP:
    """Records the matrices the GP is fitted and predicted on, then delegates.

    The fitted matrix is kept BY REFERENCE as well as by value: the acceptance criterion
    for the default path is that the GP is handed the array the biophysical path built,
    not a numerically equal copy of it.
    """

    def __init__(self, inner: Any) -> None:
        self._inner = inner
        self.fitted_ref: Optional[np.ndarray] = None
        self.fitted_on: Optional[np.ndarray] = None
        self.predicted_on: List[np.ndarray] = []

    def fit(self, X, y, **kwargs):
        self.fitted_ref = X
        self.fitted_on = np.array(X, dtype=np.float64, copy=True)
        self._inner.fit(X, y, **kwargs)
        return self

    def predict(self, X, return_std: bool = True):
        self.predicted_on.append(np.array(X, dtype=np.float64, copy=True))
        return self._inner.predict(X, return_std=return_std)


def _run_engine(
    ablate: Sequence[str] = (),
    **overrides: Any,
) -> Tuple[Any, List[_SpyGP], List[str]]:
    """One real selection round on a supplied pool, with the GP spied on.

    Returns the batch, the GP spies, and the notations handed to
    `featurize_variant_on_positions` in call order - which is how a test learns WHICH
    candidates survived the engine's own filters without reaching inside it.

    `dkl` is always ablated so the plain-GP branch is taken and the spy is reached
    whatever `TWO_STAGE_GP_ENABLED` happens to default to. The `assert spies` below then
    fails loudly rather than passing vacuously if that stops being true.
    """
    genotypes = _all_genotypes()
    spies: List[_SpyGP] = []
    real_build = mw.build_gp_model
    real_featurize = mw.featurize_variant_on_positions
    featurized: List[str] = []

    def _spy_build(*args, **kwargs):
        spy = _SpyGP(real_build(*args, **kwargs))
        spies.append(spy)
        return spy

    def _spy_featurize(variant, active_positions, wt_sequence, *args, **kwargs):
        featurized.append(variant.notation)
        return real_featurize(variant, active_positions, wt_sequence, *args, **kwargs)

    kwargs: Dict[str, Any] = dict(
        strategy=_strategy(),
        wt_sequence=WT_SEQUENCE,
        evidence=_evidence(genotypes[:N_EVIDENCE]),
        batch_size=3,
        random_seed=7,
        pdb_path=None,
        candidate_pool=[_variant_for(g) for g in genotypes[N_EVIDENCE:]],
        allowed_residues=dict(ALLOWED_RESIDUES),
        ablate={"dkl"} | set(ablate),
    )
    kwargs.update(overrides)

    with mock.patch.object(mw, "build_gp_model", _spy_build), \
            mock.patch.object(mw, "featurize_variant_on_positions", _spy_featurize):
        np.random.seed(0)
        batch = mw.select_next_micro_batch(**kwargs)
    assert spies, "the plain-GP branch was not taken; this test would prove nothing"
    return batch, spies, featurized


def _genotype_of(notation: str) -> str:
    """The full three-site genotype of a notation, wild type where it is silent."""
    return mw.genotype_on_positions(Variant(notation=notation), POSITIONS, WT_SEQUENCE)


def _expected_widths() -> Tuple[int, int]:
    alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
    one_hot_width = mw.one_hot_site_offsets(alphabets)[1]
    engine_width = len(POSITIONS) * mw.AA_FEATURE_DIM + mw.GLOBAL_EPISTATIC_DIM
    return one_hot_width, engine_width


# ======================================================================================
# 1. THE ACCEPTANCE CRITERION: with the component absent, nothing moves
# ======================================================================================
class TestTheDefaultPathIsUntouched(unittest.TestCase):
    """The whole point of making this an arm rather than a default.

    Every figure in docs/HANDOFF.md was measured with the GP on
    `featurize_variant_on_positions`. If an absent flag changed the GP's inputs even by a
    floating-point copy, every one of those figures would silently change meaning - which
    is what §5.2 means by "a genuinely new arm on every landscape including TrpB".
    """

    def test_the_engine_default_is_the_biophysical_encoding(self) -> None:
        self.assertEqual(mw.effective_surrogate_features(()), "engine_biophysical")
        self.assertEqual(mw.effective_surrogate_features(None), "engine_biophysical")
        # And it is the FIRST declared mode, which is what the manifest reports as the
        # engine default.
        self.assertEqual(mw.SURROGATE_FEATURE_MODES[0], "engine_biophysical")

    def test_with_the_component_absent_the_gp_is_handed_the_biophysical_array_itself(self) -> None:
        """Not an equal copy: the same object the biophysical path built.

        An equal copy would be correct today and is still refused, because it is the
        shape of change that makes a bitwise-identical run impossible to assert later.
        """
        _, spies, featurized = _run_engine()
        fitted = spies[0].fitted_ref
        self.assertIsInstance(fitted, np.ndarray)
        _, engine_width = _expected_widths()
        self.assertEqual(fitted.shape[1], engine_width)

        # And it IS the biophysical encoding, re-derived from the notations the engine
        # itself featurized, in its own call order.
        train_notations = featurized[:fitted.shape[0]]
        expected = np.asarray(
            [mw.featurize_variant_on_positions(
                Variant(notation=n), list(POSITIONS), WT_SEQUENCE)
             for n in train_notations],
            dtype=np.float64,
        )
        np.testing.assert_array_equal(spies[0].fitted_on, expected)

    def test_the_candidate_matrix_is_the_biophysical_one_in_the_default_mode(self) -> None:
        _, spies, featurized = _run_engine()
        candidate_matrix = max(spies[0].predicted_on, key=len)
        _, engine_width = _expected_widths()
        self.assertEqual(candidate_matrix.shape[1], engine_width)
        notations = featurized[-candidate_matrix.shape[0]:]
        expected = np.asarray(
            [mw.featurize_variant_on_positions(
                Variant(notation=n), list(POSITIONS), WT_SEQUENCE)
             for n in notations],
            dtype=np.float64,
        )
        np.testing.assert_array_equal(candidate_matrix, expected)

    def test_the_default_run_is_reproducible_to_the_bit(self) -> None:
        """Two default runs at one seed select the same batch. Weak on its own - it is
        the control for `test_the_two_arms_are_not_the_same_run` below, which is what
        makes the ablated/unablated difference attributable rather than noise."""
        first, _, _ = _run_engine()
        second, _, _ = _run_engine()
        self.assertEqual(
            [v.notation for v in first.variants], [v.notation for v in second.variants])

    def test_nothing_one_hot_is_even_constructed_in_the_default_mode(self) -> None:
        """A width computation that can raise must not run on the default path at all.

        The one-hot encoder is reached only in the ablated mode, so a library whose
        declared alphabets the encoder would reject cannot break a default run.
        """
        calls: List[int] = []
        real = mw.one_hot_encode_genotypes

        def _counted(*args, **kwargs):
            calls.append(1)
            return real(*args, **kwargs)

        with mock.patch.object(mw, "one_hot_encode_genotypes", _counted):
            _run_engine()
        self.assertEqual(calls, [], "the default path built one-hot columns it never uses")


# ======================================================================================
# 2. Selecting the arm moves the GP's inputs, and only those
# ======================================================================================
class TestTheArmMovesTheGPFeatureSpace(unittest.TestCase):

    def test_the_component_is_in_the_engine_vocabulary(self) -> None:
        self.assertIn(COMPONENT, mw.ABLATABLE_COMPONENTS)
        self.assertEqual(mw.effective_surrogate_features((COMPONENT,)), "one_hot")

    def test_it_is_distinct_from_the_components_it_could_be_confused_with(self) -> None:
        """`spatial` removes ONE cross-term from the biophysical encoding and leaves the
        rest in place; `dkl` changes the model class rather than its columns. Collapsing
        any of the three would make it impossible to ablate the feature space while
        holding the model fixed, which is the comparison this arm exists for."""
        for other in ("spatial", "spatial_filter", "dkl"):
            self.assertIn(other, mw.ABLATABLE_COMPONENTS)
            self.assertNotEqual(other, COMPONENT)
            self.assertEqual(
                mw.effective_surrogate_features((other,)), "engine_biophysical",
                f"ablating {other!r} moved the surrogate's feature space, which is a "
                f"different stage")

    def test_the_benchmark_accepts_the_name(self) -> None:
        bench = _benchmark_module()
        self.assertIn(COMPONENT, bench.ABLATABLE_COMPONENTS)
        self.assertEqual(bench.parse_ablate_components(COMPONENT), (COMPONENT,))

    def test_the_ablated_gp_sees_indicator_columns_over_the_declared_alphabets(self) -> None:
        _, spies, _ = _run_engine((COMPONENT,))
        one_hot_width, engine_width = _expected_widths()
        self.assertNotEqual(one_hot_width, engine_width,
                            "the two feature spaces must be distinguishable here")
        for index, matrix in enumerate([spies[0].fitted_on] + spies[0].predicted_on):
            with self.subTest(matrix=index):
                self.assertEqual(matrix.shape[1], one_hot_width)
                # Indicator columns: every entry 0 or 1 ...
                self.assertTrue(np.all((matrix == 0.0) | (matrix == 1.0)))
                # ... and exactly one allele active per site, on every row.
                self.assertTrue(np.all(matrix.sum(axis=1) == float(len(POSITIONS))))

    def test_the_ablated_gp_sees_exactly_the_baselines_one_hot_encoding(self) -> None:
        """The plumbed encoding, not a second implementation of the same idea.

        The claim is that this arm's GP sees the columns the benchmark's one-hot
        baselines see, so the expected matrix is built with the BENCHMARK's encoder from
        a benchmark library spec, and compared column for column.
        """
        import dataclasses

        bench = _benchmark_module()
        spec = dataclasses.replace(
            bench.LIBRARIES["cr9114_h1"],
            key="synthetic_mixed",
            label="synthetic mixed-alphabet library",
            positions=POSITIONS,
            parent_residues="".join(WT_SEQUENCE[p - 1] for p in POSITIONS),
            allowed_aa_per_position=tuple(ALLOWED_RESIDUES[p] for p in POSITIONS),
        )
        _, spies, featurized = _run_engine((COMPONENT,))

        train_notations = featurized[:spies[0].fitted_on.shape[0]]
        expected_train = bench._one_hot_matrix(
            [_genotype_of(n) for n in train_notations], spec).astype(np.float64)
        np.testing.assert_array_equal(spies[0].fitted_on, expected_train)

        candidate_matrix = max(spies[0].predicted_on, key=len)
        candidate_notations = featurized[-candidate_matrix.shape[0]:]
        expected_candidates = bench._one_hot_matrix(
            [_genotype_of(n) for n in candidate_notations], spec).astype(np.float64)
        np.testing.assert_array_equal(candidate_matrix, expected_candidates)

    def test_the_two_arms_are_not_the_same_run(self) -> None:
        """A switch that changed nothing would make every assertion above vacuous."""
        default_batch, default_spies, _ = _run_engine()
        ablated_batch, ablated_spies, _ = _run_engine((COMPONENT,))
        self.assertNotEqual(
            default_spies[0].fitted_on.shape[1], ablated_spies[0].fitted_on.shape[1])
        # Not asserted: that the SELECTED BATCH differs. On a 24-genotype toy library two
        # surrogates can agree on the top three, and pinning a difference there would be
        # pinning this fixture rather than the mechanism.
        self.assertEqual(len(default_batch.variants), len(ablated_batch.variants))

    def test_a_library_that_declares_no_alphabets_still_gets_one_hot_sites(self) -> None:
        """A saturated library declares its 20 alleles implicitly, not not-at-all."""
        _, spies, _ = _run_engine((COMPONENT,), allowed_residues=None)
        matrix = spies[0].fitted_on
        expected_width = mw.one_hot_site_offsets(
            mw.resolve_site_alphabets(POSITIONS, None))[1]
        self.assertEqual(matrix.shape[1], expected_width)
        self.assertEqual(expected_width, len(POSITIONS) * len(mw.CANONICAL_AA_ORDER))
        self.assertTrue(np.all(matrix.sum(axis=1) == float(len(POSITIONS))))


# ======================================================================================
# 3. featurize_variant_on_positions is untouched, and its other consumers still get it
# ======================================================================================
def _independent_biophysical_encoding(
    variant: Variant,
    positions: Sequence[int],
    wt_sequence: str,
) -> np.ndarray:
    """An independent oracle for `featurize_variant_on_positions`, from its contract.

    Deliberately NOT a call into the production function: docs/HANDOFF.md §5.7 says a
    test that calls the canonical implementation and compares it against itself proves
    nothing. This re-derives the documented contract - per-site 5D property deltas, six
    global summary deltas, four distance-weighted pairwise cross-terms - from the
    property table alone.
    """
    mutated = {m.position: m for m in variant.mutations}
    blocks: List[np.ndarray] = []
    for pos in positions:
        if pos in mutated:
            m = mutated[pos]
            blocks.append(
                mw.AA_PROPERTIES_NORMALIZED[m.mut_aa] - mw.AA_PROPERTIES_NORMALIZED[m.wt_aa])
        else:
            blocks.append(np.zeros(mw.AA_FEATURE_DIM, dtype=np.float64))
    linear = np.concatenate(blocks)

    muts = list(variant.mutations)
    if not muts:
        return np.concatenate([linear, np.zeros(mw.GLOBAL_EPISTATIC_DIM, dtype=np.float64)])

    deltas = [
        mw.AA_PROPERTIES_NORMALIZED[m.mut_aa] - mw.AA_PROPERTIES_NORMALIZED[m.wt_aa]
        for m in muts
    ]
    global_block = np.concatenate([[len(muts) / 3.0], np.sum(deltas, axis=0)])

    if len(muts) >= 2:
        charge = volume = hydro = dist = 0.0
        pairs = list(itertools.combinations(range(len(muts)), 2))
        for i, j in pairs:
            weight = 1.0 / np.sqrt(max(1, abs(muts[i].position - muts[j].position)))
            dist += weight
            charge += deltas[i][1] * deltas[j][1] * weight
            volume += deltas[i][2] * deltas[j][2] * weight
            hydro += deltas[i][0] * deltas[j][0] * weight
        n_pairs = float(len(pairs))
        interaction = np.array(
            [charge / n_pairs, volume / n_pairs, hydro / n_pairs, dist / n_pairs],
            dtype=np.float64)
    else:
        interaction = np.zeros(4, dtype=np.float64)

    return np.concatenate([linear, global_block, interaction])


class TestTheBiophysicalEncodingIsUnchanged(unittest.TestCase):
    """Its other callers - the prior, `AssayAdaptedHead`, the k-DPP kernel,
    `ui/landscape_plot.py` and four `scripts/run_*_benchmark.py` - must see no difference.
    With a LINEAR readout these features are not worse than one-hot and on GB1 are better
    by +0.103 rho, so moving them would be a regression there."""

    def test_output_matches_an_independent_re_derivation_of_its_contract(self) -> None:
        for genotype in _all_genotypes():
            variant = _variant_for(genotype)
            with self.subTest(genotype=genotype, notation=variant.notation):
                np.testing.assert_allclose(
                    mw.featurize_variant_on_positions(variant, list(POSITIONS), WT_SEQUENCE),
                    _independent_biophysical_encoding(variant, POSITIONS, WT_SEQUENCE),
                    rtol=0, atol=0)

    def test_the_global_block_is_still_present_and_still_redundant(self) -> None:
        """The block is NOT removed. It is a deterministic function of the rest, which is
        the diagnosis - and the reason this is a surrogate arm rather than an edit to the
        encoding."""
        width = len(POSITIONS) * mw.AA_FEATURE_DIM
        for genotype in _all_genotypes():
            variant = _variant_for(genotype)
            vector = mw.featurize_variant_on_positions(variant, list(POSITIONS), WT_SEQUENCE)
            with self.subTest(genotype=genotype):
                self.assertEqual(len(vector), width + mw.GLOBAL_EPISTATIC_DIM)
                per_site = vector[:width].reshape(len(POSITIONS), mw.AA_FEATURE_DIM)
                np.testing.assert_allclose(vector[width + 1:width + 6], per_site.sum(axis=0))
                self.assertAlmostEqual(vector[width], len(variant.mutations) / 3.0)

    def test_the_arm_does_not_change_which_variants_are_featurized_or_how(self) -> None:
        """Every `featurize_variant_on_positions` call inside the selector is identical
        in both arms: same variants, same positions, same keyword arguments. Only the GP's
        inputs move."""
        recorded: Dict[str, List[Tuple[str, Tuple[int, ...], Tuple[Any, ...]]]] = {}
        real = mw.featurize_variant_on_positions

        for label, ablate in (("default", ()), ("one_hot", (COMPONENT,))):
            calls: List[Tuple[str, Tuple[int, ...], Tuple[Any, ...]]] = []

            def _record(variant, active_positions, wt_sequence, *args, _calls=calls, **kwargs):
                _calls.append((variant.notation, tuple(active_positions),
                               tuple(sorted((k, repr(v)) for k, v in kwargs.items()))))
                return real(variant, active_positions, wt_sequence, *args, **kwargs)

            with mock.patch.object(mw, "featurize_variant_on_positions", _record):
                _run_engine(ablate)
            recorded[label] = calls

        self.assertTrue(recorded["default"],
                        "the biophysical encoding was not computed at all, so this proves nothing")
        self.assertEqual(recorded["default"], recorded["one_hot"])

    def test_the_dpp_similarity_kernel_keeps_the_biophysical_features(self) -> None:
        captured: Dict[str, np.ndarray] = {}

        def _spy_dpp(quality_scores, features, batch_size, **kwargs):
            captured["features"] = np.array(features, dtype=np.float64, copy=True)
            return list(range(min(batch_size, len(quality_scores))))

        with mock.patch.object(mw, "DPP_DIVERSITY_ENABLED", True), \
                mock.patch.object(mw, "greedy_dpp_sample", _spy_dpp):
            _, _, featurized = _run_engine((COMPONENT,))

        self.assertIn("features", captured, "k-DPP was not reached")
        notations = featurized[-captured["features"].shape[0]:]
        expected = np.asarray(
            [mw.featurize_variant_on_positions(
                Variant(notation=n), list(POSITIONS), WT_SEQUENCE)
             for n in notations],
            dtype=np.float64)
        np.testing.assert_array_equal(captured["features"], expected)
        _, engine_width = _expected_widths()
        self.assertEqual(captured["features"].shape[1], engine_width)

    def test_the_assay_adapted_head_keeps_the_biophysical_features(self) -> None:
        """It is a linear readout, and under a linear readout these features win on GB1."""
        seen: Dict[str, Any] = {}
        real_head = mw.AssayAdaptedHead

        class _SpyHead(real_head):  # type: ignore[misc, valid-type]
            def fit(self, X, y, **kwargs):
                seen["fit"] = np.array(X, dtype=np.float64, copy=True)
                return super().fit(X, y, **kwargs)

        with mock.patch.object(mw, "AssayAdaptedHead", _SpyHead):
            _run_engine((COMPONENT,))

        _, engine_width = _expected_widths()
        self.assertIn("fit", seen, "the head was never fitted")
        self.assertEqual(seen["fit"].shape[1], engine_width)


# ======================================================================================
# 4. The shared encoder: one implementation, two callers
# ======================================================================================
class TestTheOneHotEncoderIsShared(unittest.TestCase):

    def test_the_benchmark_and_the_engine_use_the_same_encoder(self) -> None:
        """Not 'produce the same answer today' - the same function object."""
        bench = _benchmark_module()
        calls: List[int] = []
        real = mw.one_hot_encode_genotypes

        def _counted(*args, **kwargs):
            calls.append(1)
            return real(*args, **kwargs)

        spec = bench.LIBRARIES["cr9114_h1"]
        with mock.patch.object(mw, "one_hot_encode_genotypes", _counted):
            bench._one_hot_matrix([spec.parent_residues], spec)
        self.assertEqual(len(calls), 1)

    def test_width_follows_the_declared_alphabets_for_every_real_library(self) -> None:
        bench = _benchmark_module()
        for key, spec in bench.LIBRARIES.items():
            with self.subTest(library=key):
                alphabets = [spec.alphabet_for(slot) for slot in range(spec.n_sites)]
                width = mw.one_hot_site_offsets(alphabets)[1]
                self.assertEqual(width, spec.one_hot_width)
                self.assertEqual(width, sum(len(a) for a in alphabets))

    def test_a_short_genotype_is_refused_rather_than_zipped(self) -> None:
        """Rule 2.1: notation omits wild-type positions, and a truncating zip silently
        inverted a published conclusion about hit diversity."""
        alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
        full = _all_genotypes()[0]
        with self.assertRaises(ValueError) as caught:
            mw.one_hot_encode_genotypes([full[:-1]], alphabets)
        self.assertIn(str(len(POSITIONS)), str(caught.exception))

    def test_wild_type_positions_are_filled_in_not_dropped(self) -> None:
        wt_genotype = "".join(WT_SEQUENCE[p - 1] for p in POSITIONS)
        # Notation naming only ONE site; the other two are wild type and silent.
        single = _variant_for(ALLOWED_RESIDUES[POSITIONS[0]][-1] + wt_genotype[1:])
        self.assertNotIn(":", single.notation)
        expanded = mw.genotype_on_positions(single, POSITIONS, WT_SEQUENCE)
        self.assertEqual(len(expanded), len(POSITIONS))
        self.assertEqual(expanded[1:], wt_genotype[1:])

    def test_an_undeclared_residue_can_be_refused_by_the_caller(self) -> None:
        alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
        intruder = next(a for a in mw.CANONICAL_AA_ORDER if a not in alphabets[0])
        genotype = intruder + _all_genotypes()[0][1:]

        def _refuse(row, slot, residue, alphabet):
            raise LookupError(f"{residue} not in {alphabet}")

        with self.assertRaises(LookupError):
            mw.one_hot_encode_genotypes([genotype], alphabets, on_undeclared=_refuse)

        # Without a callback the row still encodes, with that site's block left empty, so
        # an off-library engine proposal cannot abort a campaign mid-round.
        matrix = mw.one_hot_encode_genotypes([genotype], alphabets)
        self.assertEqual(matrix.sum(), float(len(POSITIONS) - 1))

    def test_column_order_is_deterministic_across_processes(self) -> None:
        """`CANONICAL_AMINO_ACIDS` is a set; iterating it would make the encoding depend
        on interpreter state and two runs would not be comparable."""
        self.assertEqual(mw.CANONICAL_AA_ORDER, "".join(sorted(mw.CANONICAL_AA_ORDER)))
        self.assertEqual(set(mw.CANONICAL_AA_ORDER), set("ACDEFGHIKLMNPQRSTVWY"))


# ======================================================================================
# 5. The manifest records the feature space that actually ran (rule 2.12)
# ======================================================================================
class TestTheManifestRecordsWhatTheGPWasHanded(unittest.TestCase):
    """Rule 2.12 in full: `component_states` derived `dpp` and `dkl` from
    `--ablate-components` alone while the engine ALSO gated them on a module global, so
    every manifest claimed two components were enabled when they were off, the
    configuration fingerprint was blind to the change, and the HER2 null control's
    t=-8.29 movement could not be attributed from its own artifacts.

    The lesson taken here is not "add a key". It is that the record and the behaviour must
    read the SAME gate. So:
      * the engine exposes exactly one resolver, `effective_surrogate_features`;
      * `select_next_micro_batch` calls it to decide;
      * the manifest calls it to record;
      * and the test below checks the recorded value against the width of the matrix the
        GP was ACTUALLY handed, which is neither of those two call sites.
    """

    def _manifest(self, methods: Sequence[str], disabled: Sequence[str] = ()):
        import sys
        sys.path.insert(0, "tests")
        import test_benchmark_manifest_provenance as helper

        bench = _benchmark_module()
        universe = helper._tiny_universe()
        args = helper._args(["--methods", ",".join(methods)])
        seeds = [0, 1]
        return bench.build_manifest(
            run_id="run_test_0001",
            args=args,
            universe=universe,
            evaluator=bench.Evaluator(universe),
            methods=list(methods),
            method_status={m: {"status": "available", "detail": ""} for m in methods},
            method_configs={m: {} for m in methods},
            components_disabled=tuple(disabled),
            seeds=seeds,
            initial_sets={s: universe.variant_ids[s * 5:s * 5 + 5] for s in seeds},
            outcomes=(),
            started=helper.datetime(2026, 11, 1, tzinfo=helper.timezone.utc),
            finished=helper.datetime(2026, 11, 1, tzinfo=helper.timezone.utc),
            repository_state_before=helper._FAKE_REPO_STATE,
            repository_state_after=helper._FAKE_REPO_STATE,
        )

    def _observed_mode(self, ablate: Sequence[str]) -> str:
        """Which feature space the GP was handed, read off the matrix widths.

        Independent of both the resolver and the manifest: it is the width of the array
        the spy saw, compared against the two widths derived from the declared alphabets.
        """
        _, spies, _ = _run_engine(ablate)
        one_hot_width, engine_width = _expected_widths()
        width = spies[0].fitted_on.shape[1]
        if width == one_hot_width:
            return "one_hot"
        if width == engine_width:
            return "engine_biophysical"
        raise AssertionError(f"GP was handed {width} columns, neither feature space")

    def test_the_recorded_value_matches_the_matrix_the_gp_received(self) -> None:
        bench = _benchmark_module()
        for disabled in ((), (COMPONENT,)):
            with self.subTest(disabled=disabled):
                recorded = self._manifest(
                    ["random", "panoptes_default"], disabled,
                )["configuration"]["model"]["surrogate_features"]["effective"]
                self.assertEqual(recorded, self._observed_mode(disabled))
                self.assertEqual(recorded, bench._engine_surrogate_features(disabled))

    def test_the_recorded_value_follows_the_engines_resolver_not_the_flag(self) -> None:
        """The discriminating case. If the manifest read the flag it would be right by
        coincidence in both rows above; here the resolver is forced to disagree with the
        flag, and a manifest that reads the flag fails."""
        bench = _benchmark_module()
        with mock.patch.object(mw, "effective_surrogate_features", lambda ablate=None: "one_hot"):
            recorded = self._manifest(
                ["random", "panoptes_default"], (),
            )["configuration"]["model"]["surrogate_features"]
        self.assertEqual(recorded["effective"], "one_hot")
        self.assertIn("one_hot_encode_genotypes", recorded["encoder"])
        self.assertEqual(recorded["resolver"], "workers.math_worker.effective_surrogate_features")
        # Sanity: with the patch lifted the same call records the real default again, so
        # the assertion above is about the patch and not about a sticky module state.
        self.assertEqual(
            self._manifest(["random", "panoptes_default"], ()
                           )["configuration"]["model"]["surrogate_features"]["effective"],
            bench._engine_surrogate_features(()))

    def test_the_surrogate_sentence_names_the_encoder_that_actually_runs(self) -> None:
        """Not 'the string says one-hot'. The previous value of this key said one-hot in
        every manifest this project ever wrote while the GP ran on the biophysical
        encoding, and that sentence is believed to have misdirected an entire
        investigation. A hard-coded sentence cannot satisfy both rows below."""
        symbols = {"engine_biophysical": "featurize_variant_on_positions",
                   "one_hot": "one_hot_encode_genotypes"}
        for disabled, mode in (((), "engine_biophysical"), ((COMPONENT,), "one_hot")):
            with self.subTest(mode=mode):
                model = self._manifest(
                    ["random", "panoptes_default"], disabled)["configuration"]["model"]
                self.assertEqual(model["surrogate_features"]["effective"], mode)
                self.assertIn(symbols[mode], model["surrogate"])
                encoder = model["surrogate_features"]["encoder"]
                self.assertIn(symbols[mode], encoder)
                other = next(m for m in symbols if m != mode)
                self.assertNotIn(symbols[other], encoder)

    def test_two_panoptes_arms_in_one_run_are_recorded_separately(self) -> None:
        """A run can now hold two arms whose ONLY difference is the feature space, so a
        single run-level value would be a lie for one of them."""
        bench = _benchmark_module()
        manifest = self._manifest(
            ["random", "ridge", "panoptes_default", bench.PANOPTES_ONEHOT_SURROGATE])
        per_arm = manifest["configuration"]["model"]["surrogate_features"]["per_arm"]
        self.assertEqual(set(per_arm), {"panoptes_default", bench.PANOPTES_ONEHOT_SURROGATE})
        self.assertEqual(per_arm["panoptes_default"]["effective"], "engine_biophysical")
        self.assertEqual(
            per_arm[bench.PANOPTES_ONEHOT_SURROGATE]["effective"], "one_hot")
        self.assertIn(COMPONENT,
                      per_arm[bench.PANOPTES_ONEHOT_SURROGATE]["components_disabled"])
        self.assertNotIn(COMPONENT, per_arm["panoptes_default"]["components_disabled"])
        # And each arm's recorded mode is the one the engine really hands that arm's GP.
        for arm, entry in per_arm.items():
            with self.subTest(arm=arm):
                self.assertEqual(
                    entry["effective"], self._observed_mode(entry["components_disabled"]))

    def test_the_setting_is_inside_the_hashed_configuration_subtree(self) -> None:
        """Not merely present somewhere in the manifest: inside `configuration`, which is
        the whole of what `configuration_fingerprint` hashes."""
        import json

        bench = _benchmark_module()
        default_only = self._manifest(["random", "panoptes_default"])
        with_arm = self._manifest(
            ["random", "panoptes_default", bench.PANOPTES_ONEHOT_SURROGATE])
        self.assertNotEqual(default_only["configuration_fingerprint"],
                            with_arm["configuration_fingerprint"])
        ablated = self._manifest(["random", "panoptes_default"], (COMPONENT,))
        self.assertNotEqual(default_only["configuration_fingerprint"],
                            ablated["configuration_fingerprint"])
        for manifest in (default_only, with_arm, ablated):
            configuration = manifest["configuration"]
            self.assertEqual(
                configuration["acquisition_policy"]["parameters"]["surrogate_features"],
                configuration["model"]["surrogate_features"]["effective"])
            reloaded = json.loads(json.dumps(manifest))
            self.assertEqual(
                bench.configuration_fingerprint(reloaded["configuration"]),
                reloaded["configuration_fingerprint"])
            bench.validate_manifest(manifest)  # must not raise

    def test_the_required_paths_demand_the_block(self) -> None:
        """A manifest missing the feature space cannot identify the arm, so its absence
        must be a validation failure rather than a quiet omission."""
        import copy

        bench = _benchmark_module()
        manifest = self._manifest(["random", "panoptes_default"])
        for path in ("configuration.model.surrogate_features",
                     "configuration.model.surrogate_features.effective",
                     "configuration.model.surrogate_features.per_arm"):
            with self.subTest(path=path):
                self.assertIn(path, bench.MANIFEST_REQUIRED_PATHS)
        broken = copy.deepcopy(manifest)
        del broken["configuration"]["model"]["surrogate_features"]
        with self.assertRaises(Exception):
            bench.validate_manifest(broken)


# ======================================================================================
# 6. The arm is selectable as a method, and it does not disturb the others
# ======================================================================================
class TestTheArmIsSelectableAsAMethod(unittest.TestCase):

    def test_it_is_a_registered_method_and_not_a_default_one(self) -> None:
        bench = _benchmark_module()
        self.assertIn(bench.PANOPTES_ONEHOT_SURROGATE, bench.ALL_METHODS)
        self.assertNotIn(bench.PANOPTES_ONEHOT_SURROGATE, bench.DEFAULT_METHODS)
        self.assertEqual(bench.PANOPTES_ONEHOT_SURROGATE_COMPONENT, COMPONENT)

    def test_its_selector_carries_its_own_name_and_the_component(self) -> None:
        """Both arms must appear in one run folder under distinct method keys, or they
        cannot be paired within seed - and the ablation is deliberately NOT folded into
        `panoptes_default`'s name, so the registered method is the only way to get both."""
        bench = _benchmark_module()
        library = bench.LIBRARIES["trpb"]
        metadata = bench.resolve_library_metadata(library)
        selector = bench.PanoptesSelector(
            0, metadata.sequence, metadata, False, None, library,
            ablate_components={COMPONENT},
            arm_name=bench.PANOPTES_ONEHOT_SURROGATE)
        self.assertEqual(selector.name, bench.PANOPTES_ONEHOT_SURROGATE)
        self.assertEqual(selector.surrogate_features, "one_hot")
        self.assertTrue(any("SURROGATE_FEATURES=one_hot" in n for n in selector.notes))

        default = bench.PanoptesSelector(
            0, metadata.sequence, metadata, False, None, library)
        self.assertEqual(default.name, "panoptes_default")
        self.assertEqual(default.surrogate_features, "engine_biophysical")
        self.assertFalse(any("SURROGATE_FEATURES" in n for n in default.notes),
                         "the default arm advertised a feature-space override it does not have")


if __name__ == "__main__":
    unittest.main()
