"""
tests/test_surrogate_feature_switch.py - the GP surrogate's feature space is a switch,
and moving it must not move anything else.

The defect this pins: `build_gp_model` uses an ISOTROPIC Matern, so Euclidean distance is
the kernel's only input, while `featurize_variant_on_positions` appends a 6-dimensional
global block (`n_mut / 3.0` plus the five column sums of the per-site deltas) that is a
deterministic linear function of the per-site block plus the mutation count and carries
~50% of total feature variance on every one of the nine benchmark landscapes. It adds no
information and only re-weights distance, so half of the GP's notion of "nearby sequence"
became "similar mutation burden". Measured paired over 10 shared splits with the model
held fixed, the engine's features are worse than one-hot on 8 of 9 landscapes.

Every assertion here is about MECHANISM. None pins a number, a width or a component list:
`one_hot_width == n_sites * 20` and `disabled == ["dpp", "prior"]` are both recorded in
docs/HANDOFF.md as tests that made a defect an invariant, and both cost real time. Widths
are derived from the declared alphabets; the engine encoding is checked against an
independent re-derivation of its documented contract rather than against itself.
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


# ----------------------------------------------------------------------------------------------
# A small library with MIXED per-site alphabets, which is the case that distinguishes a
# declared-alphabet encoding from an unconditional 20-wide one.
# ----------------------------------------------------------------------------------------------
POSITIONS: Tuple[int, ...] = (5, 9, 14)
# Unequal alphabet sizes, and each one contains its own parent residue, as a real
# library's does. 2 + 3 + 4 = 9 one-hot columns against 3 * 5 + 10 = 25 biophysical ones.
ALLOWED_RESIDUES: Dict[int, str] = {5: "AV", 9: "VKR", 14: "DEFG"}
WT_SEQUENCE = "MKAVAHGYVSTGCDEWQLNP"  # position 5 -> 'A', 9 -> 'V', 14 -> 'D'


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
    """A Variant whose notation omits wild-type positions, as real notation does."""
    muts = [
        f"{WT_SEQUENCE[pos - 1]}{pos}{residue}"
        for pos, residue in zip(POSITIONS, genotype)
        if residue != WT_SEQUENCE[pos - 1]
    ]
    return Variant(notation=":".join(muts) if muts else "WT")


def _all_genotypes() -> List[str]:
    return ["".join(parts) for parts in itertools.product(*(ALLOWED_RESIDUES[p] for p in POSITIONS))]


def _evidence(genotypes: Sequence[str]) -> List[Dict[str, Any]]:
    out = []
    for index, genotype in enumerate(genotypes):
        sequence = list(WT_SEQUENCE)
        for pos, residue in zip(POSITIONS, genotype):
            sequence[pos - 1] = residue
        out.append({
            "variant_sequence": "".join(sequence),
            "computed_score": 0.1 * (index + 1),
        })
    return out


class _SpyGP:
    """Records the matrices the GP is fitted and predicted on, then delegates."""

    def __init__(self, inner: Any) -> None:
        self._inner = inner
        self.fitted_on: Optional[np.ndarray] = None
        self.predicted_on: List[np.ndarray] = []

    def fit(self, X, y, **kwargs):
        self.fitted_on = np.array(X, dtype=np.float64, copy=True)
        self._inner.fit(X, y, **kwargs)
        return self

    def predict(self, X, return_std: bool = True):
        self.predicted_on.append(np.array(X, dtype=np.float64, copy=True))
        return self._inner.predict(X, return_std=return_std)


N_EVIDENCE = 6


def _run_engine(
    surrogate_features: Optional[str] = None,
    ablate: Optional[set] = None,
    **overrides: Any,
) -> Tuple[Any, List[_SpyGP], List[str]]:
    """Runs one selection round on a supplied pool, with the GP spied on.

    Returns the batch, the GP spies, and the notations handed to
    `featurize_variant_on_positions` in call order - which is how a test learns WHICH
    candidates survived the engine's own filters without reaching inside it.

    `ablate={"dkl"}` forces the plain-GP branch, so the spy is reached whatever
    TWO_STAGE_GP_ENABLED happens to default to - and `assert spies` below fails loudly
    rather than passing vacuously if that ever stops being true.
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
        ablate={"dkl"} | set(ablate or ()),
    )
    kwargs.update(overrides)
    if surrogate_features is not None:
        kwargs["surrogate_features"] = surrogate_features

    with mock.patch.object(mw, "build_gp_model", _spy_build), \
            mock.patch.object(mw, "featurize_variant_on_positions", _spy_featurize):
        np.random.seed(0)
        batch = mw.select_next_micro_batch(**kwargs)
    assert spies, "the plain-GP branch was not taken; this test would prove nothing"
    return batch, spies, featurized


def _genotype_of(notation: str) -> str:
    """The full three-site genotype of a notation, wild type where it is silent."""
    return mw.genotype_on_positions(Variant(notation=notation), POSITIONS, WT_SEQUENCE)


# ----------------------------------------------------------------------------------------------
# 1. The switch is honoured: what reaches the GP
# ----------------------------------------------------------------------------------------------
class TestTheGPReceivesTheSelectedFeatureSpace(unittest.TestCase):

    def _expected_widths(self) -> Tuple[int, int]:
        alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
        one_hot_width = mw.one_hot_site_offsets(alphabets)[1]
        engine_width = len(POSITIONS) * mw.AA_FEATURE_DIM + mw.GLOBAL_EPISTATIC_DIM
        return one_hot_width, engine_width

    def test_the_default_is_one_hot(self):
        import inspect

        default = inspect.signature(mw.select_next_micro_batch).parameters[
            "surrogate_features"].default
        self.assertEqual(default, "one_hot")
        self.assertIn(default, mw.SURROGATE_FEATURE_MODES)

    def test_one_hot_mode_feeds_the_gp_indicator_columns_over_the_declared_alphabets(self):
        _, spies, _ = _run_engine("one_hot")
        one_hot_width, engine_width = self._expected_widths()
        self.assertNotEqual(one_hot_width, engine_width,
                            "the two feature spaces must be distinguishable for this test")

        matrices = [spies[0].fitted_on] + spies[0].predicted_on
        for index, matrix in enumerate(matrices):
            with self.subTest(matrix=index):
                self.assertEqual(matrix.shape[1], one_hot_width)
                # Indicator columns: every entry is 0 or 1 ...
                self.assertTrue(np.all((matrix == 0.0) | (matrix == 1.0)))
                # ... and exactly one allele is active per site, on every row.
                self.assertTrue(np.all(matrix.sum(axis=1) == float(len(POSITIONS))))

    def test_one_hot_mode_feeds_the_gp_exactly_the_baselines_one_hot_encoding(self):
        """The plumbed encoding, not a second implementation of the same idea.

        The claim being pinned is that the GP sees the columns the benchmark's one-hot
        baselines see. So the expected matrix is built with the BENCHMARK's encoder, from
        the benchmark's own library spec, and compared column for column.
        """
        import dataclasses

        import benchmark_panoptes_isolated as bench

        spec = dataclasses.replace(
            bench.LIBRARIES["cr9114_h1"],
            key="synthetic_mixed",
            label="synthetic mixed-alphabet library",
            positions=POSITIONS,
            parent_residues="".join(WT_SEQUENCE[p - 1] for p in POSITIONS),
            allowed_aa_per_position=tuple(ALLOWED_RESIDUES[p] for p in POSITIONS),
        )
        _, spies, featurized = _run_engine("one_hot")

        # Training rows, in order: the WT reference, then the evidence genotypes. The
        # notations are taken from the engine's own call order, so the comparison is
        # against whatever it actually encoded rather than against an assumed pool.
        train_notations = featurized[:spies[0].fitted_on.shape[0]]
        expected_train = bench._one_hot_matrix(
            [_genotype_of(n) for n in train_notations], spec,
        ).astype(np.float64)
        np.testing.assert_array_equal(spies[0].fitted_on, expected_train)

        candidate_matrix = max(spies[0].predicted_on, key=len)
        candidate_notations = featurized[-candidate_matrix.shape[0]:]
        expected_candidates = bench._one_hot_matrix(
            [_genotype_of(n) for n in candidate_notations], spec,
        ).astype(np.float64)
        np.testing.assert_array_equal(candidate_matrix, expected_candidates)

    def test_engine_mode_feeds_the_gp_the_biophysical_encoding_instead(self):
        _, spies, featurized = _run_engine("engine_biophysical")
        _, engine_width = self._expected_widths()

        matrices = [spies[0].fitted_on] + spies[0].predicted_on
        for index, matrix in enumerate(matrices):
            with self.subTest(matrix=index):
                self.assertEqual(matrix.shape[1], engine_width)

        candidate_matrix = max(spies[0].predicted_on, key=len)
        candidate_notations = featurized[-candidate_matrix.shape[0]:]
        expected_candidates = np.asarray(
            [mw.featurize_variant_on_positions(
                Variant(notation=n), list(POSITIONS), WT_SEQUENCE)
             for n in candidate_notations],
            dtype=np.float64,
        )
        np.testing.assert_array_equal(candidate_matrix, expected_candidates)

    def test_the_two_modes_are_not_the_same_run(self):
        """A switch that changed nothing would make every assertion above vacuous."""
        _, one_hot_spies, _ = _run_engine("one_hot")
        _, engine_spies, _ = _run_engine("engine_biophysical")
        self.assertNotEqual(
            one_hot_spies[0].fitted_on.shape[1], engine_spies[0].fitted_on.shape[1],
        )

    def test_an_unknown_mode_is_refused_rather_than_silently_running_the_default(self):
        with self.assertRaises(ValueError) as caught:
            _run_engine("onehot")
        self.assertIn("surrogate_features", str(caught.exception))
        for mode in mw.SURROGATE_FEATURE_MODES:
            self.assertIn(mode, str(caught.exception))

    def test_a_library_that_declares_no_alphabets_still_gets_one_hot_sites(self):
        """A saturated library declares its 20 alleles implicitly, not not-at-all."""
        _, spies, _ = _run_engine("one_hot", allowed_residues=None)
        matrix = spies[0].fitted_on
        expected_width = mw.one_hot_site_offsets(
            mw.resolve_site_alphabets(POSITIONS, None))[1]
        self.assertEqual(matrix.shape[1], expected_width)
        self.assertEqual(expected_width, len(POSITIONS) * len(mw.CANONICAL_AA_ORDER))
        self.assertTrue(np.all(matrix.sum(axis=1) == float(len(POSITIONS))))


# ----------------------------------------------------------------------------------------------
# 2. featurize_variant_on_positions is untouched, and its other consumers still get it
# ----------------------------------------------------------------------------------------------
def _independent_biophysical_encoding(
    variant: Variant,
    positions: Sequence[int],
    wt_sequence: str,
) -> np.ndarray:
    """An independent oracle for `featurize_variant_on_positions`, from its docstring.

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
                mw.AA_PROPERTIES_NORMALIZED[m.mut_aa] - mw.AA_PROPERTIES_NORMALIZED[m.wt_aa]
            )
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
    column_sums = np.sum(deltas, axis=0)
    global_block = np.concatenate([[len(muts) / 3.0], column_sums])

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
            dtype=np.float64,
        )
    else:
        interaction = np.zeros(4, dtype=np.float64)

    return np.concatenate([linear, global_block, interaction])


class TestTheBiophysicalEncodingIsUnchanged(unittest.TestCase):
    """Its other callers - the prior, AssayAdaptedHead, the DPP kernel,
    ui/landscape_plot.py and four scripts/run_*_benchmark.py - must see no difference.
    With a LINEAR readout these features are not worse than one-hot and on GB1 are
    better by +0.103 rho, so changing them would be a regression there."""

    def test_output_matches_an_independent_re_derivation_of_its_contract(self):
        for genotype in _all_genotypes():
            variant = _variant_for(genotype)
            with self.subTest(genotype=genotype, notation=variant.notation):
                np.testing.assert_allclose(
                    mw.featurize_variant_on_positions(variant, list(POSITIONS), WT_SEQUENCE),
                    _independent_biophysical_encoding(variant, POSITIONS, WT_SEQUENCE),
                    rtol=0, atol=0,
                )

    def test_the_global_block_is_still_present_and_still_redundant(self):
        """The block is NOT removed. It is a deterministic function of the rest, which is
        the diagnosis - and the reason the fix is a surrogate switch rather than an edit
        to this function."""
        width = len(POSITIONS) * mw.AA_FEATURE_DIM
        for genotype in _all_genotypes():
            variant = _variant_for(genotype)
            vector = mw.featurize_variant_on_positions(variant, list(POSITIONS), WT_SEQUENCE)
            with self.subTest(genotype=genotype):
                self.assertEqual(len(vector), width + mw.GLOBAL_EPISTATIC_DIM)
                per_site = vector[:width].reshape(len(POSITIONS), mw.AA_FEATURE_DIM)
                # The five global summary deltas are the column sums of the per-site block.
                np.testing.assert_allclose(vector[width + 1:width + 6], per_site.sum(axis=0))
                # And the sixth is the mutation count, recoverable the same way.
                self.assertAlmostEqual(vector[width], len(variant.mutations) / 3.0)

    def test_the_switch_does_not_change_which_variants_are_featurized_or_how(self):
        """Every `featurize_variant_on_positions` call inside the selector is identical
        in both modes: same variants, same positions, same keyword arguments. Only the
        GP's inputs move."""
        recorded: Dict[str, List[Tuple[str, Tuple[int, ...], Tuple[Any, ...]]]] = {}
        real = mw.featurize_variant_on_positions

        for mode in ("one_hot", "engine_biophysical"):
            calls: List[Tuple[str, Tuple[int, ...], Tuple[Any, ...]]] = []

            def _record(variant, active_positions, wt_sequence, *args, _calls=calls, **kwargs):
                _calls.append((
                    variant.notation,
                    tuple(active_positions),
                    tuple(sorted((k, repr(v)) for k, v in kwargs.items())),
                ))
                return real(variant, active_positions, wt_sequence, *args, **kwargs)

            with mock.patch.object(mw, "featurize_variant_on_positions", _record):
                _run_engine(mode)
            recorded[mode] = calls

        self.assertTrue(recorded["one_hot"], "the biophysical encoding was not computed at all")
        self.assertEqual(recorded["one_hot"], recorded["engine_biophysical"])

    def test_the_dpp_similarity_kernel_keeps_the_biophysical_features_in_one_hot_mode(self):
        captured: Dict[str, np.ndarray] = {}

        def _spy_dpp(quality_scores, features, batch_size, **kwargs):
            captured["features"] = np.array(features, dtype=np.float64, copy=True)
            return list(range(min(batch_size, len(quality_scores))))

        with mock.patch.object(mw, "DPP_DIVERSITY_ENABLED", True), \
                mock.patch.object(mw, "greedy_dpp_sample", _spy_dpp):
            _, _, featurized = _run_engine("one_hot")

        self.assertIn("features", captured, "k-DPP was not reached")
        notations = featurized[-captured["features"].shape[0]:]
        expected = np.asarray(
            [mw.featurize_variant_on_positions(
                Variant(notation=n), list(POSITIONS), WT_SEQUENCE)
             for n in notations],
            dtype=np.float64,
        )
        np.testing.assert_array_equal(captured["features"], expected)
        engine_width = len(POSITIONS) * mw.AA_FEATURE_DIM + mw.GLOBAL_EPISTATIC_DIM
        self.assertEqual(captured["features"].shape[1], engine_width)

    def test_the_assay_adapted_head_keeps_the_biophysical_features_in_one_hot_mode(self):
        """It is a linear readout, and under a linear readout these features win on GB1."""
        seen: Dict[str, Any] = {}
        real_head = mw.AssayAdaptedHead

        class _SpyHead(real_head):  # type: ignore[misc, valid-type]
            def fit(self, X, y, **kwargs):
                seen["fit"] = np.array(X, dtype=np.float64, copy=True)
                return super().fit(X, y, **kwargs)

            def predict(self, X, **kwargs):
                seen["predict"] = np.array(X, dtype=np.float64, copy=True)
                return super().predict(X, **kwargs)

        with mock.patch.object(mw, "AssayAdaptedHead", _SpyHead):
            _run_engine("one_hot")

        engine_width = len(POSITIONS) * mw.AA_FEATURE_DIM + mw.GLOBAL_EPISTATIC_DIM
        self.assertIn("fit", seen, "the head was never fitted")
        self.assertEqual(seen["fit"].shape[1], engine_width)
        if "predict" in seen:
            self.assertEqual(seen["predict"].shape[1], engine_width)


# ----------------------------------------------------------------------------------------------
# 3. The shared encoder: one implementation, two callers
# ----------------------------------------------------------------------------------------------
class TestTheOneHotEncoderIsShared(unittest.TestCase):

    def test_the_benchmark_and_the_engine_use_the_same_encoder(self):
        """Not 'produce the same answer today' - the same function object."""
        import benchmark_panoptes_isolated as bench

        calls: List[int] = []
        real = mw.one_hot_encode_genotypes

        def _counted(*args, **kwargs):
            calls.append(1)
            return real(*args, **kwargs)

        spec = bench.LIBRARIES["cr9114_h1"]
        with mock.patch.object(mw, "one_hot_encode_genotypes", _counted):
            bench._one_hot_matrix([spec.parent_residues], spec)
        self.assertEqual(len(calls), 1)

    def test_width_follows_the_declared_alphabets_for_every_real_library(self):
        import benchmark_panoptes_isolated as bench

        for key, spec in bench.LIBRARIES.items():
            with self.subTest(library=key):
                alphabets = [spec.alphabet_for(slot) for slot in range(spec.n_sites)]
                width = mw.one_hot_site_offsets(alphabets)[1]
                self.assertEqual(width, spec.one_hot_width)
                self.assertEqual(width, sum(len(a) for a in alphabets))

    def test_a_short_genotype_is_refused_rather_than_zipped_against_the_site_list(self):
        """HANDOFF rule 2.1: notation omits wild-type positions, and a truncating zip
        silently inverted a published conclusion about hit diversity."""
        alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
        full = _all_genotypes()[0]
        with self.assertRaises(ValueError) as caught:
            mw.one_hot_encode_genotypes([full[:-1]], alphabets)
        self.assertIn(str(len(POSITIONS)), str(caught.exception))

    def test_wild_type_positions_are_filled_in_not_dropped(self):
        wt_genotype = "".join(WT_SEQUENCE[p - 1] for p in POSITIONS)
        # Notation naming only ONE site; the other two are wild type and silent.
        single = _variant_for(ALLOWED_RESIDUES[POSITIONS[0]][-1] + wt_genotype[1:])
        self.assertNotIn(":", single.notation)
        expanded = mw.genotype_on_positions(single, POSITIONS, WT_SEQUENCE)
        self.assertEqual(len(expanded), len(POSITIONS))
        self.assertEqual(expanded[1:], wt_genotype[1:])

    def test_an_undeclared_residue_can_be_refused_by_the_caller(self):
        alphabets = mw.resolve_site_alphabets(POSITIONS, ALLOWED_RESIDUES)
        intruder = next(a for a in mw.CANONICAL_AA_ORDER if a not in alphabets[0])
        genotype = intruder + _all_genotypes()[0][1:]

        def _refuse(row, slot, residue, alphabet):
            raise LookupError(f"{residue} not in {alphabet}")

        with self.assertRaises(LookupError):
            mw.one_hot_encode_genotypes([genotype], alphabets, on_undeclared=_refuse)

        # Without a callback the row still encodes, with that site's block left empty,
        # so an off-library engine proposal cannot abort a campaign.
        matrix = mw.one_hot_encode_genotypes([genotype], alphabets)
        self.assertEqual(matrix.sum(), float(len(POSITIONS) - 1))

    def test_column_order_is_deterministic_across_processes(self):
        """`CANONICAL_AMINO_ACIDS` is a set; iterating it would make the encoding depend
        on interpreter state, and two runs would not be comparable."""
        self.assertEqual(mw.CANONICAL_AA_ORDER, "".join(sorted(mw.CANONICAL_AA_ORDER)))
        self.assertEqual(set(mw.CANONICAL_AA_ORDER), set("ACDEFGHIKLMNPQRSTVWY"))


if __name__ == "__main__":
    unittest.main()
