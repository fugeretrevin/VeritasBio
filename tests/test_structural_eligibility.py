from __future__ import annotations

import itertools
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import numpy as np
from pydantic import ValidationError

from controller.gatekeeper import Variant
from controller.orchestrator import OrchestratorConfig, _apply_structural_filter
from controller.prompts import (
    ResidueCluster,
    SearchStrategy,
    build_strategy_user_prompt,
    get_strategy_system_prompt,
)
from workers import math_worker
from workers.math_worker import SelectedBatch, compute_biological_prior
from workers import spatial_filter
from workers.spatial_filter import (
    SpatialContactMap,
    filter_spatially_viable_variants,
    get_spatial_contact_map,
)


POSITIONS = (183, 184, 227, 228)
MUTATIONS = ("V183A", "F184I", "V227K", "S228G")
def write_pdb(path: Path, coords: dict[int, tuple[str, float]], chains: dict[int, str] | None = None) -> None:
    with path.open("w", encoding="ascii") as pdb_file:
        for serial, (position, (residue, x)) in enumerate(coords.items(), 1):
            chain = chains.get(position, "A") if chains else "A"
            pdb_file.write(
                f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{position:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          C  \n"
            )


def trpb_reference() -> str:
    sequence = list("A" * 228)
    sequence[182] = "V"
    sequence[183] = "F"
    sequence[226] = "V"
    sequence[227] = "S"
    return "".join(sequence)


def four_site_variants(wt_sequence: str) -> list[Variant]:
    variants: list[Variant] = []
    for depth in range(1, len(MUTATIONS) + 1):
        for subset in itertools.combinations(MUTATIONS, depth):
            variant = Variant(notation=":".join(subset))
            variant.validate_against_wt(wt_sequence)
            variants.append(variant)
    return variants


class StubSurrogate:
    def __init__(self, *args, **kwargs) -> None:
        pass

    def fit(self, features: np.ndarray, values: np.ndarray) -> None:
        pass

    def predict(self, features: np.ndarray, return_std: bool = False):
        means = np.zeros(len(features), dtype=np.float64)
        stds = np.ones(len(features), dtype=np.float64)
        return (means, stds) if return_std else means


class TestStructuralEligibility(unittest.TestCase):
    def test_structure_is_opt_in_and_does_not_constrain_planning(self) -> None:
        metadata = type(
            "Metadata",
            (),
            {
                "catalytic_positions": [],
                "binding_positions": [183, 184, 227, 228],
                "chromophore_residues": [],
                "pdb_id": "8EH0",
            },
        )()
        prompt = build_strategy_user_prompt("TrpB metadata", protein_metadata=metadata)
        self.assertNotIn("8EH0", prompt)
        self.assertNotIn("13 Angstroms", prompt)
        self.assertIn("Do NOT require cluster residues to be close in 3D space", get_strategy_system_prompt())
        self.assertNotIn("BURIED HYDROPHOBIC CORE STABILITY", get_strategy_system_prompt())
        self.assertIsNone(OrchestratorConfig(target_accession="TRPB_TM9D8").structure_pdb_id)
        self.assertEqual(OrchestratorConfig(structure_pdb_id="1EMA").structure_pdb_id, "1EMA")

    def test_all_deep_variants_reach_acquisition(self) -> None:
        wt_sequence = trpb_reference()
        variants = four_site_variants(wt_sequence)
        notations = [variant.notation for variant in variants]
        self.assertEqual(len(notations), 15)
        self.assertEqual(len(set(notations)), 15)
        stage_counts = {
            "generation": len(variants),
            "validated_unique_unmeasured": len(notations),
            "measured_removed": 0,
        }
        depth_histogram = Counter(len(variant.mutations) for variant in variants)

        with tempfile.TemporaryDirectory() as temp_dir:
            pdb_path = Path(temp_dir) / "distant_sites.pdb"
            write_pdb(
                pdb_path,
                {
                    183: ("VAL", 0.0),
                    184: ("PHE", 3.8),
                    227: ("VAL", 100.0),
                    228: ("SER", 104.0),
                },
            )
            legacy_approved, legacy_pruned = filter_spatially_viable_variants(
                ["V183A:F184I:V227K:S228G"], pdb_path=pdb_path, enable_heuristics=True
            )
            self.assertEqual(legacy_approved, [])
            self.assertTrue(legacy_pruned)

            strategy = SearchStrategy(
                target_protein="TrpB",
                target_accession="TRPB_TM9D8",
                hypothesis="Test combinatorial exploration across four assay sites.",
                target_clusters=[
                    ResidueCluster(
                        cluster_id="four_sites",
                        positions=list(POSITIONS),
                        rationale="Deterministic regression fixture.",
                    )
                ],
                max_mutation_depth=4,
            )
            structural_stage: list[str] = []
            real_filter = filter_spatially_viable_variants

            def capture_structural_stage(candidate_notations: list[str], **kwargs):
                stage_counts["before_structural"] = len(candidate_notations)
                result = real_filter(candidate_notations, **kwargs)
                stage_counts["after_structural"] = len(result[0])
                structural_stage.extend(result[0])
                return result

            real_expected_improvement = math_worker.compute_expected_improvement

            def check_acquisition_entry(means, stds, f_best, xi=0.01):
                stage_counts["acquisition_entry"] = len(structural_stage)
                self.assertEqual(structural_stage, notations)
                self.assertEqual(Counter(len(item.split(":")) for item in structural_stage), depth_histogram)
                return real_expected_improvement(means, stds, f_best, xi)

            with (
                patch.object(math_worker, "DeepKernelLearner", StubSurrogate),
                patch.object(math_worker, "generate_candidate_variants", return_value=variants),
                patch.object(math_worker, "sample_gflownet_candidates", return_value=[]),
                patch.object(spatial_filter, "filter_spatially_viable_variants", side_effect=capture_structural_stage),
                patch.object(math_worker, "compute_expected_improvement", side_effect=check_acquisition_entry),
                patch.object(math_worker, "greedy_dpp_sample", side_effect=lambda **kwargs: list(range(5))),
            ):
                selected = math_worker.select_next_micro_batch(
                    strategy=strategy,
                    wt_sequence=wt_sequence,
                    batch_size=5,
                    prior_mean_fn=lambda variant: 0.0,
                    acquisition_mode="ei",
                    pdb_path=pdb_path,
                )

        self.assertEqual(structural_stage, notations)
        self.assertIn("V183A:F184I:V227K:S228G", structural_stage)
        stage_counts["final_selection"] = len(selected.variants)
        self.assertEqual(len(selected.variants), 5)
        self.assertEqual(
            stage_counts,
            {
                "generation": 15,
                "validated_unique_unmeasured": 15,
                "measured_removed": 0,
                "before_structural": 15,
                "after_structural": 15,
                "acquisition_entry": 15,
                "final_selection": 5,
            },
        )

    def test_post_selection_default_preserves_batch_and_metadata(self) -> None:
        variants = [Variant(notation="V183A:F184I:V227K:S228G"), Variant(notation="V183A")]
        batch = SelectedBatch(
            variants=variants,
            expected_improvements=[0.1, 0.2],
            predicted_means=[1.0, 2.0],
            predicted_uncertainties=[0.3, 0.4],
            pareto_ranks=[1, 0],
        )
        result = _apply_structural_filter(
            batch, None, set(), enabled=OrchestratorConfig().enable_structural_filter, reference_sequence=trpb_reference()
        )
        self.assertIs(result, batch)
        self.assertEqual(result.variants, variants)
        self.assertEqual(result.expected_improvements, [0.1, 0.2])
        self.assertEqual(result.predicted_means, [1.0, 2.0])
        self.assertEqual(result.predicted_uncertainties, [0.3, 0.4])
        self.assertEqual(result.pareto_ranks, [1, 0])

    def test_opt_in_filter_keeps_metadata_aligned(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb_path = Path(temp_dir) / "distant_sites.pdb"
            write_pdb(
                pdb_path,
                {
                    183: ("VAL", 0.0),
                    184: ("PHE", 3.8),
                    227: ("VAL", 100.0),
                    228: ("SER", 104.0),
                },
            )
            batch = SelectedBatch(
                variants=[Variant(notation="V183A:F184I:V227K:S228G"), Variant(notation="V183A")],
                expected_improvements=[0.1, 0.2],
                predicted_means=[1.0, 2.0],
                predicted_uncertainties=[0.3, 0.4],
                pareto_ranks=[1, 0],
            )
            result = _apply_structural_filter(batch, pdb_path, set(), enabled=True, reference_sequence=trpb_reference())
            self.assertEqual([variant.notation for variant in result.variants], ["V183A"])
            self.assertEqual(result.expected_improvements, [0.2])
            self.assertEqual(result.predicted_means, [2.0])
            self.assertEqual(result.predicted_uncertainties, [0.4])
            self.assertEqual(result.pareto_ranks, [0])

    def test_missing_and_sequential_structures_do_not_reuse_geometry(self) -> None:
        self.assertEqual(SpatialContactMap().ca_coordinates, {})
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            missing = SpatialContactMap(temp_path / "absent.pdb")
            self.assertEqual(missing.ca_coordinates, {})
            first_path = temp_path / "first.pdb"
            second_path = temp_path / "second.pdb"
            write_pdb(first_path, {183: ("VAL", 0.0), 184: ("PHE", 5.0)})
            write_pdb(second_path, {183: ("VAL", 0.0), 184: ("PHE", 25.0)})
            first = get_spatial_contact_map(first_path)
            second = get_spatial_contact_map(second_path)
            self.assertAlmostEqual(first.get_distance(183, 184), 5.0)
            self.assertAlmostEqual(second.get_distance(183, 184), 25.0)
            self.assertIsNot(first, second)
            mismatched_path = temp_path / "mismatched.pdb"
            write_pdb(mismatched_path, {183: ("ALA", 0.0), 184: ("PHE", 5.0)})
            self.assertEqual(SpatialContactMap(mismatched_path, reference_sequence=trpb_reference()).ca_coordinates, {})

            split_chains_path = temp_path / "split_chains.pdb"
            write_pdb(
                split_chains_path,
                {
                    183: ("VAL", 0.0),
                    184: ("PHE", 5.0),
                    227: ("VAL", 100.0),
                    228: ("SER", 104.0),
                },
                chains={183: "A", 184: "A", 227: "B", 228: "B"},
            )
            self.assertNotIn(227, SpatialContactMap(split_chains_path).ca_coordinates)
            self.assertEqual(SpatialContactMap(split_chains_path, chain_id="B").ca_coordinates.keys(), {227, 228})
            self.assertEqual(
                filter_spatially_viable_variants(
                    ["V183A:F184I:V227K:S228G"], pdb_path=missing.pdb_path, enable_heuristics=True
                )[0],
                ["V183A:F184I:V227K:S228G"],
            )

    def test_structural_prior_contribution_is_disabled_with_or_without_context(self) -> None:
        variant = Variant(notation="V183A:F184I:V227K:S228G")
        baseline = compute_biological_prior(variant)
        with patch("workers.spatial_filter.compute_structural_surrogate_score", return_value=-100.0) as surrogate:
            contextual = compute_biological_prior(variant, pdb_path="explicit-structure.pdb")
        self.assertEqual(contextual, baseline)
        surrogate.assert_not_called()

    def test_sequence_and_explicit_restriction_validation_remain(self) -> None:
        wt_sequence = trpb_reference()
        with self.assertRaises(ValidationError):
            Variant(notation="X183A")
        with self.assertRaises(ValidationError):
            Variant(notation="V183A:V183K")
        with self.assertRaises(ValueError):
            Variant(notation="A999V").validate_against_wt(wt_sequence)
        with self.assertRaises(ValueError):
            Variant(notation="A183K").validate_against_wt(wt_sequence)
        with self.assertRaises(ValidationError):
            SearchStrategy(
                target_protein="TrpB",
                target_accession="TRPB_TM9D8",
                hypothesis="Forbidden positions are not searchable.",
                target_clusters=[
                    ResidueCluster(cluster_id="restricted", positions=[183], rationale="Test constraint.")
                ],
                forbidden_positions=[183],
            )


if __name__ == "__main__":
    unittest.main()