"""
tests/test_spatial_filter.py — Unit tests for 3D Spatial Distance & Steric Clash Filter.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from workers.spatial_filter import (
    SpatialContactMap,
    filter_spatially_viable_variants,
    get_spatial_contact_map,
)


class TestSpatialFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.pdb_path = Path(__file__).resolve().parents[1] / "workers" / "data" / "1EMA.pdb"
        self.cmap = get_spatial_contact_map(self.pdb_path)

    def test_default_filter_is_identity_without_structure(self) -> None:
        variants = ["Y39K:F46A:M153A", "V183A:F184I:V227K:S228G"]
        approved, pruned = filter_spatially_viable_variants(variants)
        self.assertEqual(approved, variants)
        self.assertEqual(pruned, [])
        self.assertEqual(get_spatial_contact_map().ca_coordinates, {})

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

        approved, pruned = filter_spatially_viable_variants(
            ["F64L:S65T:V68L:T203Y"], pdb_path=self.pdb_path, enable_heuristics=True
        )
        self.assertEqual(approved, ["F64L:S65T:V68L:T203Y"])
        self.assertEqual(len(pruned), 0)

    def test_uncoupled_dollar_run_pruned(self) -> None:
        """The hallucinated distal cluster Y39K:F46A:M153A must be pruned."""
        approved, pruned = filter_spatially_viable_variants(
            ["Y39K:F46A:M153A"], pdb_path=self.pdb_path, enable_heuristics=True
        )
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
        approved, pruned = filter_spatially_viable_variants(
            ["F64L:Y66A:T203L"], pdb_path=self.pdb_path, enable_heuristics=True
        )
        self.assertEqual(len(approved), 0)
        self.assertEqual(len(pruned), 1)
        self.assertIn("Catalytic Core Invariant", pruned[0][1])

        # G67A breaks backbone twist
        approved_g, pruned_g = filter_spatially_viable_variants(
            ["G67A"], pdb_path=self.pdb_path, enable_heuristics=True
        )
        self.assertEqual(len(approved_g), 0)
        self.assertEqual(len(pruned_g), 1)
        self.assertIn("Gly67 provides the tight steric backbone hinge", pruned_g[0][1])

        # Y66H retains aromatic ring (viable blue fluorescent protein)
        approved_h, pruned_h = filter_spatially_viable_variants(
            ["Y66H"], pdb_path=self.pdb_path, enable_heuristics=True
        )
        self.assertEqual(approved_h, ["Y66H"])
        self.assertEqual(len(pruned_h), 0)

        # Authentic synergistic winners pass
        approved_syn, pruned_syn = filter_spatially_viable_variants(
            ["S65T:T203Y", "F64L:S65T"], pdb_path=self.pdb_path, enable_heuristics=True
        )
        self.assertEqual(approved_syn, ["S65T:T203Y", "F64L:S65T"])
        self.assertEqual(len(pruned_syn), 0)


if __name__ == "__main__":
    unittest.main()
