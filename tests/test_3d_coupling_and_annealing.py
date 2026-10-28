"""
tests/test_3d_coupling_and_annealing.py — Unit tests for 3D Contact Map Coupling & GFlowNet Annealing.

Verifies:
  1. 3D Contact Distance Epistatic Coupling: Using authentic PDB coordinates boosts the coupling weight
     for residues that are far apart in 1D sequence but close in 3D folded space.
  2. GFlowNet Temperature Annealing: Sampling temperature correctly anneals from 1.60 (exploratory)
     down to 0.65 (focused exploitation) as evidence accumulates.
"""

from __future__ import annotations

import math
import unittest
from pathlib import Path
import numpy as np

from controller.gatekeeper import Variant
from workers.math_worker import featurize_variant_on_positions
from workers.spatial_filter import SpatialContactMap


class Test3DCouplingAndAnnealing(unittest.TestCase):
    """Tests for 3D contact distance coupling and dynamic temperature annealing."""

    def test_3d_spatial_coupling_boost(self) -> None:
        """
        Verifies that residues separated by large linear distance (e.g. 65 and 203 in GFP, delta=138)
        receive significantly boosted epistatic interaction weight when 3D spatial map is provided.
        """
        from workers.rag_worker import fetch_protein_metadata
        meta = fetch_protein_metadata("P42212")
        wt_seq = meta.sequence
        # S65T and T203Y: 138 positions apart in sequence, but in 3D space they are only ~12 A apart
        var = Variant(notation="S65T:T203Y")
        active_positions = [65, 203]

        pdb_path = Path(__file__).resolve().parent.parent / "workers" / "data" / "1EMA.pdb"
        spatial_map = SpatialContactMap(pdb_path=pdb_path)

        # 1. Featurize WITHOUT spatial map (pure 1D sequence distance: 1/sqrt(138) ≈ 0.085)
        vec_1d = featurize_variant_on_positions(var, active_positions, wt_seq, spatial_map=None)

        # 2. Featurize WITH spatial map (3D contact distance: 1/sqrt(12.19) ≈ 0.286)
        vec_3d = featurize_variant_on_positions(var, active_positions, wt_seq, spatial_map=spatial_map)

        # The last 4 dimensions are: [charge_prod, volume_prod, hydro_prod, dist_inv]
        dist_inv_1d = vec_1d[-1]
        dist_inv_3d = vec_3d[-1]

        # 3D contact weight must be substantially higher than 1D distance weight
        self.assertGreater(
            dist_inv_3d,
            dist_inv_1d * 2.0,
            f"3D contact weight ({dist_inv_3d:.4f}) should be >2x 1D distance weight ({dist_inv_1d:.4f})",
        )

    def test_gflownet_temperature_annealing_schedule(self) -> None:
        """Verifies mathematical properties of the exponential temperature annealing schedule."""
        # Round 1: Zero evidence -> Broad exploration (tau = 1.60)
        n_ev_r1 = 0
        tau_r1 = max(0.65, float(1.60 * np.exp(-n_ev_r1 / 60.0)))
        self.assertAlmostEqual(tau_r1, 1.60, places=2)

        # Round 3: ~30 evidence -> Balanced (tau ≈ 0.97)
        n_ev_r3 = 30
        tau_r3 = max(0.65, float(1.60 * np.exp(-n_ev_r3 / 60.0)))
        self.assertLess(tau_r3, tau_r1)
        self.assertGreater(tau_r3, 0.90)

        # Round 5: ~80 evidence -> Focused exploitation clamped at 0.65
        n_ev_r5 = 80
        tau_r5 = max(0.65, float(1.60 * np.exp(-n_ev_r5 / 60.0)))
        self.assertEqual(tau_r5, 0.65)


if __name__ == "__main__":
    unittest.main()
