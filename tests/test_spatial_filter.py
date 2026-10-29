"""
tests/test_spatial_filter.py — Unit tests for 3D Spatial Distance & Steric Clash Filter.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from workers.spatial_filter import (
    SpatialContactMap,
    compute_structural_surrogate_score,
    filter_spatially_viable_variants,
    get_spatial_contact_map,
)


class TestSpatialFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.cmap = get_spatial_contact_map()

    def test_ca_coordinates_loaded(self) -> None:
        """Verifies that 1EMA.pdb C-alpha coordinates were parsed."""
        self.assertGreater(len(self.cmap.ca_coordinates), 200)
        self.assertIn(64, self.cmap.ca_coordinates)
        self.assertIn(65, self.cmap.ca_coordinates)
        self.assertIn(203, self.cmap.ca_coordinates)

    def test_distance_calculation(self) -> None:
        """Verifies accurate Euclidean distance calculation."""
        d_64_68 = self.cmap.get_distance(64, 68)
        self.assertIsNotNone(d_64_68)
        # Residues 64 and 68 are neighboring strands on the beta barrel (~7.32 A)
        self.assertAlmostEqual(d_64_68, 7.32, delta=0.2)

        # Distance from residue to itself is 0.0
        self.assertEqual(self.cmap.get_distance(64, 64), 0.0)

    def test_superfolder_hit_passes_cohesion(self) -> None:
        """The authentic superfolder GFP quadruple hit must pass spatial cohesion."""
        cohesion = self.cmap.evaluate_mutation_cluster_cohesion([64, 65, 68, 203])
        self.assertTrue(cohesion["is_cohesive"])
        self.assertLess(cohesion["mean_distance_angstroms"], 11.0)

        approved, pruned = filter_spatially_viable_variants(["F64L:S65T:V68L:T203Y"])
        self.assertEqual(approved, ["F64L:S65T:V68L:T203Y"])
        self.assertEqual(len(pruned), 0)

    def test_uncoupled_dollar_run_pruned(self) -> None:
        """The hallucinated distal cluster Y39K:F46A:M153A must be pruned."""
        approved, pruned = filter_spatially_viable_variants(["Y39K:F46A:M153A"])
        self.assertEqual(len(approved), 0)
        self.assertEqual(len(pruned), 1)
        self.assertIn("Spatially uncoupled", pruned[0][1])

    def test_buried_core_collapse_detected(self) -> None:
        """Mutating buried Phe46 to Ala without compensating partner must be rejected."""
        passes, reason = self.cmap.check_steric_and_core_clash([(46, "F", "A")])
        self.assertFalse(passes)
        self.assertIn("Buried core void collapse", reason)

    def test_steric_clash_detected(self) -> None:
        """Mutually expanding adjacent spatial neighbors to bulky aromatics must clash."""
        # Residues 64 and 68 are 7.32 A apart; mutating both to Trp causes steric collision
        passes, reason = self.cmap.check_steric_and_core_clash([(64, "F", "W"), (68, "V", "W")])
        self.assertFalse(passes)
        self.assertIn("Steric clash", reason)

    def test_electrostatic_clash_detected(self) -> None:
        """Introducing like charges at close spatial contact must be flagged."""
        # Residues 140 and 172 are surface loops near each other
        passes, reason = self.cmap.check_steric_and_core_clash([(64, "F", "R"), (68, "V", "R")])
        self.assertFalse(passes)
        self.assertTrue("Electrostatic" in reason or "desolvation" in reason)

    def test_catalytic_invariant_blocks_dark_mutants(self) -> None:
        """Catalytic invariant mask must prune non-aromatic 66 and all 67 mutants."""
        # Y66A destroys the conjugated fluorophore ring
        approved, pruned = filter_spatially_viable_variants(["F64L:Y66A:T203L"])
        self.assertEqual(len(approved), 0)
        self.assertEqual(len(pruned), 1)
        self.assertIn("Catalytic Core Invariant", pruned[0][1])

        # G67A breaks backbone twist
        approved_g, pruned_g = filter_spatially_viable_variants(["G67A"])
        self.assertEqual(len(approved_g), 0)
        self.assertEqual(len(pruned_g), 1)
        self.assertIn("Gly67 provides the tight steric backbone hinge", pruned_g[0][1])

        # Y66H retains aromatic ring (viable blue fluorescent protein)
        approved_h, pruned_h = filter_spatially_viable_variants(["Y66H"])
        self.assertEqual(approved_h, ["Y66H"])
        self.assertEqual(len(pruned_h), 0)

        # Authentic synergistic winners pass
        approved_syn, pruned_syn = filter_spatially_viable_variants(["S65T:T203Y", "F64L:S65T"])
        self.assertEqual(approved_syn, ["S65T:T203Y", "F64L:S65T"])
        self.assertEqual(len(pruned_syn), 0)

    def test_trpb_multiloop_active_site_quadruple_passes(self) -> None:
        """Authentic multi-loop active site quadruple hit (TrpB 183:184:227:228 at 18.95 A) passes cohesion."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if trpb_pdb.is_file():
            cmap_trpb = get_spatial_contact_map(trpb_pdb)
            cohesion = cmap_trpb.evaluate_mutation_cluster_cohesion([183, 184, 227, 228])
            self.assertTrue(cohesion["is_cohesive"], f"Failed cohesion: max distance {cohesion.get('max_distance_angstroms')} A")
            self.assertLess(cohesion["max_distance_angstroms"], 24.0)

            approved, pruned = filter_spatially_viable_variants(
                ["V183A:F184I:V227K:S228G"],
                pdb_path=trpb_pdb,
                enable_heuristics=True,
                allowed_active_positions={183, 184, 227, 228},
            )
            self.assertEqual(approved, ["V183A:F184I:V227K:S228G"])
            self.assertEqual(len(pruned), 0)

    def test_enable_heuristics_false_bypasses_filtering(self) -> None:
        """When enable_heuristics=False, filter acts as an identity function without pruning."""
        variants = ["Y39K:F46A:M153A", "F64L:Y66A:T203L"]
        approved, pruned = filter_spatially_viable_variants(variants, enable_heuristics=False)
        self.assertEqual(approved, variants)
        self.assertEqual(pruned, [])

    def test_reference_sequence_mismatch_clears_defensively(self) -> None:
        """Passing a mismatched reference sequence clears coordinates to prevent corrupt spatial mapping."""
        bad_ref = "A" * 250  # All alanines, does not match 1EMA GFP sequence
        cmap = get_spatial_contact_map(reference_sequence=bad_ref)
        self.assertEqual(cmap.ca_coordinates, {})

    def test_hetatm_cofactor_and_ligand_parsing(self) -> None:
        """Verifies that HETATM cofactors (0JO) and ligands (QOX) in 8EH0 are parsed."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        cmap = get_spatial_contact_map(trpb_pdb)
        self.assertIn("0JO", cmap.cofactor_coordinates)
        self.assertGreater(len(cmap.cofactor_coordinates["0JO"]), 15)
        self.assertIn("QOX", cmap.ligand_coordinates)
        self.assertGreater(len(cmap.catalytic_center_coordinates), 30)

    def test_cofactor_distance_and_zones(self) -> None:
        """Verifies 3D distance and zone classification to catalytic cofactor."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        cmap = get_spatial_contact_map(trpb_pdb)
        d_227 = cmap.get_cofactor_distance(227)
        d_228 = cmap.get_cofactor_distance(228)
        d_184 = cmap.get_cofactor_distance(184)

        self.assertIsNotNone(d_227)
        self.assertIsNotNone(d_228)
        self.assertIsNotNone(d_184)

        # 227 and 228 are in catalytic core shell (< 4.5 A)
        self.assertLess(d_227, 4.5)
        self.assertLess(d_228, 4.5)
        self.assertEqual(cmap.get_cofactor_zone(227), "catalytic_core")
        self.assertEqual(cmap.get_cofactor_zone(228), "catalytic_core")

        # 184 is in substrate entrance channel (4.5 A to 12.0 A)
        self.assertGreaterEqual(d_184, 4.5)
        self.assertLessEqual(d_184, 12.0)
        self.assertEqual(cmap.get_cofactor_zone(184), "entrance_channel")

    def test_cofactor_catalytic_clash_rejection(self) -> None:
        """Expanding residue 227 (< 4.2 A to 0JO) to Trp must trigger cofactor catalytic clash."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        cmap = get_spatial_contact_map(trpb_pdb)
        passes, reason = cmap.check_steric_and_core_clash([(227, "V", "W")])
        self.assertFalse(passes)
        self.assertIn("Cofactor catalytic clash", reason)

    def test_cofactor_geometry_surrogate_scoring(self) -> None:
        """Steric decongestion at entrance channel gets bonus; clash at core receives penalty."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        # F184I relieves entrance channel steric crowding -> +0.45 bonus
        score_decongest = compute_structural_surrogate_score("F184I", pdb_path=trpb_pdb)
        self.assertGreater(score_decongest, 0.0)

        # V227W expands into immediate cofactor pocket -> -3.5 penalty
        score_clash = compute_structural_surrogate_score("V227W", pdb_path=trpb_pdb)
        self.assertLess(score_clash, -3.0)


    def test_cofactor_charge_and_proline_disruption(self) -> None:
        """Proline and uncompensated acidic charge in immediate cofactor shell (< 4.2 A) must be rejected."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        cmap = get_spatial_contact_map(trpb_pdb)
        # Residue 227 is 3.8 A from 0JO cofactor; mutating to Pro kinks catalytic backbone
        passes_pro, reason_pro = cmap.check_steric_and_core_clash([(227, "V", "P")])
        self.assertFalse(passes_pro)
        self.assertIn("Cofactor catalytic loop disruption", reason_pro)

        # Mutating residue 227 to Asp introduces negative charge into phosphate binding sphere (even when exempt from buried core)
        passes_asp, reason_asp = cmap.check_steric_and_core_clash([(227, "V", "D")], allowed_active_positions={227})
        self.assertFalse(passes_asp)
        self.assertIn("Cofactor catalytic charge disruption", reason_asp)

    def test_extract_functional_pocket_with_cofactors(self) -> None:
        """extract_functional_pocket with include_cofactors=True pulls all cofactor-lining residues."""
        trpb_pdb = Path(__file__).resolve().parents[1] / "workers" / "data" / "8EH0.pdb"
        if not trpb_pdb.is_file():
            self.skipTest("8EH0.pdb not available")

        cmap = get_spatial_contact_map(trpb_pdb)
        pocket = cmap.extract_functional_pocket([83], cutoff_angstroms=8.5, include_cofactors=True)
        # Must include catalytic lysine 83 and active loop residues 227, 228
        self.assertIn(83, pocket)
        self.assertIn(227, pocket)
        self.assertIn(228, pocket)
        self.assertGreater(len(pocket), 25)


if __name__ == "__main__":
    unittest.main()

