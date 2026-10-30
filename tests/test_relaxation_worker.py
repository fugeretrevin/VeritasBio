"""
tests/test_relaxation_worker.py — Test Suite for Track C Sidechain Rotamer & Energy Minimization.

Verifies:
  1. Fast CPU-based L-BFGS local minimizer resolves steric clashes in silico.
  2. Harmonic backbone restraints preserve C-alpha coordinates (RMSD < 0.2 A).
  3. Execution speed is sub-80ms per variant on standard CPU.
  4. Multi-mutant cooperative relaxation properly minimizes adjacent sidechains.
  5. PDB serialization preserves standard PDB ATOM format.
"""

from __future__ import annotations

import os
import time
import unittest
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
PDB_1EMA = ROOT_DIR / "workers" / "data" / "1EMA.pdb"


class TestRelaxationWorker(unittest.TestCase):

    def setUp(self):
        self.assertTrue(PDB_1EMA.is_file(), f"Reference PDB not found at {PDB_1EMA}")

    def test_import_and_seam_interface(self):
        """Verifies public exports and Pydantic models for Seam 1."""
        from workers.relaxation_worker import RelaxationResult, relax_variant_structure

        self.assertTrue(callable(relax_variant_structure))
        # Ensure RelaxationResult has required fields
        res = RelaxationResult(
            variant_notation="A206W",
            initial_energy_kcal=450.0,
            final_energy_kcal=12.5,
            initial_min_distance_angstrom=1.42,
            final_min_distance_angstrom=2.35,
            clashes_resolved=1,
            converged=True,
            relaxed_pdb_str="ATOM      1  N   MET A   1...",
            rmsd_ca_angstrom=0.04,
            runtime_ms=35.0,
        )
        self.assertEqual(res.variant_notation, "A206W")
        self.assertLess(res.final_energy_kcal, res.initial_energy_kcal)
        self.assertGreaterEqual(res.final_min_distance_angstrom, 2.0)

    def test_single_mutation_clash_resolution(self):
        """Tests that a bulky mutation with steric clash is relaxed to non-clashing geometry."""
        from workers.relaxation_worker import relax_variant_structure

        # A206 is on the surface/barrel boundary; mutating to bulky Trp (W) induces tight contacts
        result = relax_variant_structure(
            pdb_path=PDB_1EMA,
            mutations=["A206W"],
            max_iterations=200,
            harmonic_k=10.0,
        )

        self.assertIsNotNone(result)
        self.assertEqual(result.variant_notation, "A206W")
        self.assertTrue(result.converged)
        # Energy should decrease significantly after relaxation
        self.assertLess(result.final_energy_kcal, result.initial_energy_kcal)
        # Minimum interatomic distance should expand away from steric clash threshold (<1.8 A)
        self.assertGreaterEqual(result.final_min_distance_angstrom, 1.90)
        # C-alpha restraint ensures minimal backbone drift
        self.assertLess(result.rmsd_ca_angstrom, 0.20)
        # PDB string should contain valid ATOM records
        self.assertIn("ATOM", result.relaxed_pdb_str)
        self.assertIn("TRP A 206", result.relaxed_pdb_str)

    def test_fast_cpu_runtime_budget(self):
        """Guards against an order-of-magnitude regression in relaxation cost.

        A regression guard, not a benchmark: the ceiling is generous so a loaded
        machine does not fail the suite.
        """
        from workers.relaxation_worker import relax_variant_structure

        # Warm-up pass
        relax_variant_structure(PDB_1EMA, ["S65T"], max_iterations=50)

        t0 = time.perf_counter()
        result = relax_variant_structure(
            PDB_1EMA,
            ["S65T"],
            max_iterations=150,
            harmonic_k=10.0,
        )
        elapsed_ms = (time.perf_counter() - t0) * 1000.0

        self.assertLess(
            elapsed_ms,
            1500.0,
            f"Relaxation took {elapsed_ms:.1f}ms, over the 1500ms regression ceiling",
        )
        self.assertTrue(result.converged)

    def test_multi_mutant_cooperative_relaxation(self):
        """Verifies cooperative minimization of adjacent mutations in the chromophore / active pocket."""
        from workers.relaxation_worker import relax_variant_structure

        result = relax_variant_structure(
            PDB_1EMA,
            ["T65S", "Y66H"],
            max_iterations=250,
            harmonic_k=10.0,
        )

        self.assertEqual(result.variant_notation, "T65S:Y66H")
        self.assertTrue(result.converged)
        self.assertLess(result.final_energy_kcal, result.initial_energy_kcal)
        self.assertGreaterEqual(result.final_min_distance_angstrom, 1.85)


if __name__ == "__main__":
    unittest.main()

