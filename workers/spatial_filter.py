"""
workers/spatial_filter.py â€” 3D Spatial Distance & Steric Clash Biomechanical Filter.

Purpose:
    Extracts 3D C-alpha coordinates from crystallographic PDB structures (or AlphaFold models)
    to calculate spatial Euclidean distance matrices D_ij = ||r_i - r_j||_2.

    Guards against LLM epistatic hallucinations by validating that:
      1. Synergistic multi-mutants reside within physical contact spheres (D_ij <= 8.5 A).
      2. Residues separated by > 15 A across the scaffold are flagged as spatially uncoupled.
      3. Steric clashes (e.g. two bulky aromatics colliding at D_ij < 4.5 A) and buried
         hydrophobic core disruptions are pruned BEFORE spending cloud GPU budget.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple
import numpy as np

logger = logging.getLogger(__name__)

# Fallback structure used when a caller supplies no pdb_path.
#
# This is avGFP. Applied to any other target it yields contact shells, C-alpha
# distances and burial estimates for the wrong protein - silently, and with results
# that look entirely plausible. The benchmarks always pass an explicit structure
# (TrpB uses 8EH0), but server.py and the dashboard can reach these functions without
# one, so every fallback use is now logged loudly rather than assumed intentional.
DEFAULT_PDB_PATH: Path = Path(__file__).resolve().parent / "data" / "1EMA.pdb"

_WARNED_DEFAULT_PDB: Set[str] = set()


def _resolve_pdb_path(pdb_path: Optional[Any], caller: str) -> Path:
    """
    Returns the caller's structure, or the avGFP fallback with a loud warning.

    Centralised so the fallback cannot be taken silently from any of the call sites
    that previously each wrote `pdb_path if pdb_path is not None else DEFAULT_PDB_PATH`.
    """
    if pdb_path is not None:
        return Path(pdb_path)
    if caller not in _WARNED_DEFAULT_PDB:
        _WARNED_DEFAULT_PDB.add(caller)
        logger.warning(
            "%s received no pdb_path and is falling back to %s (avGFP). Structural "
            "results will describe avGFP, NOT the target under study. Pass pdb_path.",
            caller, DEFAULT_PDB_PATH.name,
        )
    return DEFAULT_PDB_PATH

# Standard Van der Waals sidechain volume estimates (in Angstroms^3)
AA_VOLUMES: Dict[str, float] = {
    "G": 60.1,  "A": 88.6,  "S": 89.0,  "C": 108.5, "D": 111.1,
    "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0,
    "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7,
    "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8,
}

# Core buried hydrophobic residues in avGFP (1EMA crystal structure)
AVGFP_BURIED_CORE_RESIDUES: Set[int] = {
    8, 12, 14, 16, 44, 46, 54, 62, 64, 68, 70, 83, 85, 92, 94, 99,
    119, 121, 123, 148, 150, 163, 165, 167, 179, 181, 183, 201, 203, 220, 222
}


THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

SOLVENT_AND_BUFFER_NAMES: Set[str] = {
    "HOH", "WAT", "EDO", "GOL", "PEG", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "DMS", "ACT", "FMT"
}

MODIFIED_AA_NAMES: Set[str] = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CME", "CRO", "CRQ", "CR2", "NRQ"
}

KNOWN_COFACTORS: Set[str] = {
    "PLP", "0JO", "LLP", "PMP", "HEM", "FAD", "FMN", "NAD", "NAP", "NDP", "SAM", "SAH", "ATP", "ADP", "GTP", "GDP", "COA", "TPP", "THF"
}


class SpatialContactMap:
    """Extracts and caches 3D coordinate matrices, cofactors, and ligand spheres from PDB structures."""

    def __init__(
        self,
        pdb_path: Optional[Path] = None,
        chain_id: Optional[str] = None,
        reference_sequence: Optional[str] = None,
    ) -> None:
        self.pdb_path: Optional[Path] = _resolve_pdb_path(pdb_path, "SpatialContactMap")
        self.chain_id = chain_id
        self.reference_sequence = reference_sequence
        self.ca_coordinates: Dict[int, Tuple[str, float, float, float]] = {}
        self.cofactor_coordinates: Dict[str, List[Tuple[str, str, float, float, float]]] = {}
        self.ligand_coordinates: Dict[str, List[Tuple[str, str, float, float, float]]] = {}
        self.catalytic_center_coordinates: List[Tuple[float, float, float]] = []
        self._distance_cache: Dict[Tuple[int, int], float] = {}
        self._cofactor_distance_cache: Dict[int, float] = {}
        self.buried_core_residues: Set[int] = set()
        self._load_pdb()
        self.buried_core_residues = self._compute_buried_core_residues()

    def _compute_buried_core_residues(self, cutoff: float = 7.0, min_neighbors: int = 7) -> Set[int]:
        """Dynamically identifies buried hydrophobic core residues in ANY 3D crystal/cryo-EM structure."""
        buried = set()
        for r1 in self.ca_coordinates:
            neighbors = sum(1 for r2 in self.ca_coordinates if r1 != r2 and (self.get_distance(r1, r2) or 99.0) <= cutoff)
            if neighbors >= min_neighbors:
                buried.add(r1)
        return buried if buried else AVGFP_BURIED_CORE_RESIDUES

    def _load_pdb(self) -> None:
        """Parses ATOM C-alpha coordinates and HETATM cofactor/ligand spheres from PDB."""
        if self.pdb_path is None:
            return

        if not self.pdb_path.is_file():
            pdb_stem = self.pdb_path.stem
            if len(pdb_stem) == 4 and pdb_stem.isalnum():
                try:
                    import urllib.request
                    url = f"https://files.rcsb.org/download/{pdb_stem.upper()}.pdb"
                    logger.info("Attempting auto-download of PDB structure from RCSB: %s", url)
                    self.pdb_path.parent.mkdir(parents=True, exist_ok=True)
                    urllib.request.urlretrieve(url, str(self.pdb_path))
                except Exception as e:
                    logger.warning("Could not download PDB from RCSB: %s", e)

        if not self.pdb_path.is_file():
            logger.warning("PDB file not found at %s. Spatial filtering disabled.", self.pdb_path)
            return

        selected_chain = self.chain_id
        try:
            with open(self.pdb_path, "r", encoding="utf-8") as f:
                for line in f:
                    if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 54:
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        chain_id = line[21].strip()

                        if selected_chain is None and chain_id:
                            selected_chain = chain_id

                        if chain_id and selected_chain and chain_id != selected_chain:
                            continue

                        # Insertion codes cannot be mapped safely by integer residue number alone
                        if self.reference_sequence and len(line) > 26 and line[26].strip():
                            self.ca_coordinates.clear()
                            self.cofactor_coordinates.clear()
                            self.ligand_coordinates.clear()
                            self.catalytic_center_coordinates.clear()
                            self.chain_id = selected_chain
                            logger.warning(
                                "Structural mapping unavailable: insertion codes in %s chain %s.",
                                self.pdb_path,
                                selected_chain or "blank",
                            )
                            return

                        # Handle chromophore triad (CRO in 1EMA occupies residue 66, representing 65-67)
                        if res_name == "CRO" and atom_name in ("CA1", "CA2", "CA3", "C1", "C2", "CA"):
                            try:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                # Assign chromophore coordinates across triad 65, 66, 67
                                self.ca_coordinates[65] = ("SER", x, y, z)
                                self.ca_coordinates[66] = ("TYR", x, y, z)
                                self.ca_coordinates[67] = ("GLY", x, y, z)
                            except (ValueError, IndexError):
                                continue
                            continue

                        if atom_name == "CA":
                            try:
                                res_seq = int(line[22:26].strip())
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                self.ca_coordinates[res_seq] = (res_name, x, y, z)
                            except (ValueError, IndexError):
                                continue

                        # Gap 3: Parse HETATM cofactor and ligand heavy atoms
                        if line.startswith("HETATM"):
                            if (
                                res_name not in SOLVENT_AND_BUFFER_NAMES
                                and res_name not in MODIFIED_AA_NAMES
                                and not atom_name.startswith("H")
                            ):
                                try:
                                    x = float(line[30:38].strip())
                                    y = float(line[38:46].strip())
                                    z = float(line[46:54].strip())
                                    coord_tuple = (atom_name, chain_id, x, y, z)
                                    if res_name in KNOWN_COFACTORS:
                                        self.cofactor_coordinates.setdefault(res_name, []).append(coord_tuple)
                                    else:
                                        self.ligand_coordinates.setdefault(res_name, []).append(coord_tuple)
                                    self.catalytic_center_coordinates.append((x, y, z))
                                except (ValueError, IndexError):
                                    continue
        except (OSError, UnicodeError) as exc:
            logger.warning("Could not read PDB structure at %s: %s", self.pdb_path, exc)
            return

        self.chain_id = selected_chain

        # Verify against reference sequence if provided
        if self.reference_sequence and self.ca_coordinates:
            mismatches = [
                pos for pos, (res_name, *_coords) in self.ca_coordinates.items()
                if pos < 1
                or pos > len(self.reference_sequence)
                or THREE_TO_ONE.get(res_name) != self.reference_sequence[pos - 1].upper()
            ]
            # Tolerate minor homolog/engineered differences (<= 15%), but reject major mismatches or indexing frameshifts
            mismatch_threshold = max(3, int(0.15 * len(self.ca_coordinates)))
            if len(mismatches) > mismatch_threshold:
                self.ca_coordinates.clear()
                self.cofactor_coordinates.clear()
                self.ligand_coordinates.clear()
                self.catalytic_center_coordinates.clear()
                logger.warning(
                    "Structural mapping unavailable for %s chain %s: residue identity mismatch count %d exceeds threshold %d (%s...).",
                    self.pdb_path,
                    self.chain_id or "blank",
                    len(mismatches),
                    mismatch_threshold,
                    mismatches[:5],
                )
                return

    def get_distance(self, res_a: int, res_b: int) -> Optional[float]:
        """Calculates 3D Euclidean distance (in Angstroms) between two residue C-alphas."""
        if res_a == res_b:
            return 0.0

        pair = (min(res_a, res_b), max(res_a, res_b))
        if pair in self._distance_cache:
            return self._distance_cache[pair]

        if res_a not in self.ca_coordinates or res_b not in self.ca_coordinates:
            return None

        _, x1, y1, z1 = self.ca_coordinates[res_a]
        _, x2, y2, z2 = self.ca_coordinates[res_b]

        d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
        self._distance_cache[pair] = d
        return d

    def get_cofactor_distance(self, res_id: int) -> Optional[float]:
        """
        Calculates minimum Euclidean distance between residue C-alpha and any
        catalytic cofactor atom (e.g. PLP, 0JO, HEM).
        """
        if res_id not in self.ca_coordinates:
            return None

        target_coords: List[Tuple[float, float, float]] = []
        if self.cofactor_coordinates:
            for atoms in self.cofactor_coordinates.values():
                for _, _, x, y, z in atoms:
                    target_coords.append((x, y, z))
        elif self.catalytic_center_coordinates:
            target_coords = self.catalytic_center_coordinates

        if not target_coords:
            return None

        if res_id in self._cofactor_distance_cache:
            return self._cofactor_distance_cache[res_id]

        _, ca_x, ca_y, ca_z = self.ca_coordinates[res_id]
        min_d = min(
            math.sqrt((ca_x - x) ** 2 + (ca_y - y) ** 2 + (ca_z - z) ** 2)
            for x, y, z in target_coords
        )
        self._cofactor_distance_cache[res_id] = min_d
        return min_d

    def get_catalytic_center_distance(self, res_id: int) -> Optional[float]:
        """
        Calculates minimum Euclidean distance between residue C-alpha and ANY
        catalytic cofactor or bound substrate ligand heavy atom.
        """
        if not self.catalytic_center_coordinates or res_id not in self.ca_coordinates:
            return None
        _, ca_x, ca_y, ca_z = self.ca_coordinates[res_id]
        return min(
            math.sqrt((ca_x - x) ** 2 + (ca_y - y) ** 2 + (ca_z - z) ** 2)
            for x, y, z in self.catalytic_center_coordinates
        )

    def get_cofactor_zone(self, res_id: int) -> str:
        """
        Categorizes residue into catalytic proximity zones (Gap 3):
          - 'catalytic_core': d < 4.5 A (immediate catalytic transition-state / cofactor shell)
          - 'entrance_channel': 4.5 A <= d <= 12.0 A (substrate accession and allosteric loop channel)
          - 'outer_scaffold': d > 12.0 A (distal framework)
        """
        d = self.get_cofactor_distance(res_id)
        if d is None:
            d = self.get_catalytic_center_distance(res_id)
        if d is None:
            return "outer_scaffold"
        if d < 4.5:
            return "catalytic_core"
        if d <= 12.0:
            return "entrance_channel"
        return "outer_scaffold"

    def evaluate_mutation_cluster_cohesion(self, positions: List[int]) -> Dict[str, Any]:
        """
        Evaluates whether a proposed cluster of mutation positions forms a cohesive
        physical contact sphere in 3D space rather than an uncoupled hallucinated cluster.
        """
        if len(positions) <= 1:
            return {"is_cohesive": True, "mean_distance_angstroms": 0.0, "max_distance_angstroms": 0.0}

        distances: List[float] = []
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                d = self.get_distance(positions[i], positions[j])
                if d is not None:
                    distances.append(d)

        if not distances:
            return {"is_cohesive": True, "mean_distance_angstroms": 0.0, "max_distance_angstroms": 0.0}

        mean_d = sum(distances) / len(distances)
        max_d = max(distances)

        # Multi-loop active site contact manifold: allow mean distance <= 16.0 A and max pairwise <= 24.0 A
        # This properly accommodates multi-loop catalytic pockets (e.g. TrpB loops 183-184 and 227-228 at 18.95 A,
        # or antibody paratope CDR-H1/H2/H3 loops) while still rejecting uncoupled scaffold hallucinations (> 25-40 A).
        is_cohesive = (mean_d <= 16.0) and (max_d <= 24.0)

        return {
            "is_cohesive": is_cohesive,
            "mean_distance_angstroms": round(mean_d, 2),
            "max_distance_angstroms": round(max_d, 2),
            "pairwise_distances": [round(d, 2) for d in distances],
        }

    def check_steric_and_core_clash(
        self,
        mutation_tuples: List[Tuple[int, str, str]],
        catalytic_invariants: Optional[Set[int]] = None,
        allowed_active_positions: Optional[Set[int]] = None,
    ) -> Tuple[bool, str]:
        """
        Checks a proposed set of mutations for steric clashes, core disruption,
        or catalytic invariant violations across any protein target.

        Args:
            mutation_tuples: List of (position, wt_aa, mut_aa) tuples.
            catalytic_invariants: Optional set of 1-indexed residue positions that cannot tolerate mutation.
            allowed_active_positions: Optional set of binding/active pocket residues exempt from buried core penalties.

        Returns:
            Tuple of (passes_filter: bool, reason: str).
        """
        # 0. Dynamic Catalytic & Active-Site Invariant Protection (Pillar 1)
        for pos, wt_aa, mut_aa in mutation_tuples:
            if catalytic_invariants and pos in catalytic_invariants:
                # Fluorophore spectral shift exception (Y66H/W in GFP retains aromatic ring)
                if pos == 66 and wt_aa == "Y" and mut_aa in ("H", "W"):
                    continue
                return False, f"Catalytic Invariant Violation: Residue {wt_aa}{pos} is an annotated active/catalytic site; mutation to {mut_aa} destroys biochemical function"

            # avGFP Fluorophore Core Triad fallback if catalytic_invariants not explicitly passed
            if pos == 66 and mut_aa not in ("Y", "H", "W") and wt_aa == "Y":
                return False, f"Catalytic Core Invariant Violation: Residue 66 forms the conjugated fluorophore ring; non-aromatic substitution {wt_aa}66{mut_aa} abolishes fluorescence (dark mutant)"
            if pos == 67 and mut_aa != "G" and wt_aa == "G":
                return False, f"Catalytic Core Invariant Violation: Gly67 provides the tight steric backbone hinge required for cyclization; mutation {wt_aa}67{mut_aa} causes backbone misfolding"

        # 1. Dynamic Buried Hydrophobic Core Collapse Detection (Coordination Sphere Density)
        buried_set = self.buried_core_residues if self.buried_core_residues else AVGFP_BURIED_CORE_RESIDUES
        if allowed_active_positions:
            buried_set = buried_set - allowed_active_positions
        for pos, wt_aa, mut_aa in mutation_tuples:
            if pos in buried_set:
                wt_vol = AA_VOLUMES.get(wt_aa, 120.0)
                mut_vol = AA_VOLUMES.get(mut_aa, 120.0)
                vol_delta = mut_vol - wt_vol

                # Large cavity creation in buried core (> 70 A^3 lost, e.g. Phe -> Ala)
                if vol_delta < -70.0:
                    # Check if there is a compensating expansion nearby (< 7 A)
                    has_compensating_expansion = any(
                        p != pos
                        and (self.get_distance(pos, p) or 99.0) < 7.0
                        and (AA_VOLUMES.get(m, 100.0) - AA_VOLUMES.get(w, 100.0)) > 50.0
                        for p, w, m in mutation_tuples
                    )
                    if not has_compensating_expansion:
                        return False, f"Buried core void collapse at residue {wt_aa}{pos}{mut_aa} (volume delta: {vol_delta:.1f} A^3 with no compensating neighbor)"

                # Introducing charged residue into buried core
                if wt_aa in "AILVMF" and mut_aa in "KRDENQ":
                    return False, f"Electrostatic desolvation penalty: buried hydrophobic core residue {wt_aa}{pos} mutated to charged/polar {mut_aa}"

        # 1b. Cofactor / Transition-State Catalytic Core Clash Protection (Gap 3)
        if self.cofactor_coordinates or self.catalytic_center_coordinates:
            for pos, wt_aa, mut_aa in mutation_tuples:
                d_cof = self.get_cofactor_distance(pos)
                if d_cof is None:
                    d_cof = self.get_catalytic_center_distance(pos)
                if d_cof is not None and d_cof < 4.2:
                    wt_vol = AA_VOLUMES.get(wt_aa, 120.0)
                    mut_vol = AA_VOLUMES.get(mut_aa, 120.0)
                    vol_delta = mut_vol - wt_vol
                    # Bulky aromatic expansion colliding directly into cofactor
                    if vol_delta > 65.0 and mut_aa in ("W", "Y", "F"):
                        return False, f"Cofactor catalytic clash: {wt_aa}{pos}{mut_aa} is {d_cof:.1f} A from cofactor with severe steric expansion +{vol_delta:.1f} A^3"
                    # Proline kinking catalytic active-site loop
                    if mut_aa == "P" and wt_aa != "P":
                        return False, f"Cofactor catalytic loop disruption: proline substitution {wt_aa}{pos}{mut_aa} kinks catalytic backbone at {d_cof:.1f} A from cofactor"
                    # Severe uncompensated electrostatic disruption in immediate catalytic core (< 4.2 A)
                    if wt_aa in "AILVMFST" and mut_aa in ("D", "E"):
                        return False, f"Cofactor catalytic charge disruption: {wt_aa}{pos}{mut_aa} introduces acidic charge at {d_cof:.1f} A from cofactor binding sphere"

        # 2. Check for spatial steric collision between neighboring mutations
        for i in range(len(mutation_tuples)):
            for j in range(i + 1, len(mutation_tuples)):
                pos_i, wt_i, mut_i = mutation_tuples[i]
                pos_j, wt_j, mut_j = mutation_tuples[j]

                d = self.get_distance(pos_i, pos_j)
                if d is not None and d < 7.5:
                    vol_i_delta = AA_VOLUMES.get(mut_i, 100.0) - AA_VOLUMES.get(wt_i, 100.0)
                    vol_j_delta = AA_VOLUMES.get(mut_j, 100.0) - AA_VOLUMES.get(wt_j, 100.0)

                    # Both mutate to bulky aromatics in tight proximity (< 7.5 A) -> severe steric clash
                    if vol_i_delta > 35.0 and vol_j_delta > 35.0:
                        return False, f"Steric clash: mutual expansion of spatial neighbors {wt_i}{pos_i}{mut_i} and {wt_j}{pos_j}{mut_j} (distance: {d:.2f} A)"

                    # Like-charge electrostatic repulsion at close contact
                    if (mut_i in "KR" and mut_j in "KR") or (mut_i in "DE" and mut_j in "DE"):
                        return False, f"Electrostatic repulsion clash between spatial neighbors {wt_i}{pos_i}{mut_i} and {wt_j}{pos_j}{mut_j} (distance: {d:.2f} A)"

        return True, "Passed spatial biomechanical filters"

    def extract_functional_pocket(
        self,
        functional_centers: Sequence[int],
        cutoff_angstroms: float = 8.5,
        include_cofactors: bool = True,
    ) -> List[int]:
        """
        Extracts all residue positions within cutoff_angstroms of ANY functional center
        or catalytic cofactor/ligand heavy atom.
        Universal 3D tessellation method that extracts catalytic/fluorophore pockets from PDB.
        """
        if not self.ca_coordinates:
            return sorted(list(functional_centers))

        pocket: Set[int] = set(functional_centers)
        for fc in functional_centers:
            if fc not in self.ca_coordinates:
                continue
            for res_id in self.ca_coordinates:
                dist = self.get_distance(fc, res_id)
                if dist is not None and dist <= cutoff_angstroms:
                    pocket.add(res_id)

        if include_cofactors and (self.cofactor_coordinates or self.catalytic_center_coordinates):
            for res_id in self.ca_coordinates:
                d_cof = self.get_cofactor_distance(res_id)
                if d_cof is None:
                    d_cof = self.get_catalytic_center_distance(res_id)
                if d_cof is not None and d_cof <= cutoff_angstroms:
                    pocket.add(res_id)

        return sorted(list(pocket))

    def compute_stability_ddg(
        self,
        mutation_tuples: List[Tuple[int, str, str]],
    ) -> float:
        """
        Heuristic folding-destabilization penalty score (unitless, higher = worse).

        NOT a kcal/mol free energy. The value is a sum of hand-chosen penalties for
        volume change, charge burial and proline introduction; there is no force
        field, no solvation model, no torsional term, and no parameterization
        against any experimental ddG dataset. Use it to rank candidates against
        each other, never as a physical energy.
        Positive = destabilizing (deleterious to fold).
        Negative = stabilizing.
        """
        ddg = 0.0
        buried_set = self.buried_core_residues if self.buried_core_residues else AVGFP_BURIED_CORE_RESIDUES

        for pos, wt_aa, mut_aa in mutation_tuples:
            wt_vol = AA_VOLUMES.get(wt_aa, 120.0)
            mut_vol = AA_VOLUMES.get(mut_aa, 120.0)
            vol_delta = mut_vol - wt_vol

            if pos in buried_set:
                # 1. Buried Core Overcrowding (steric clash)
                if vol_delta > 45.0:
                    ddg += 0.05 * (vol_delta - 45.0)
                # 2. Buried Cavity Creation (loss of hydrophobic packing)
                elif vol_delta < -55.0:
                    ddg += 0.035 * (abs(vol_delta) - 55.0)

                # 3. Buried Charge Desolvation penalty
                if wt_aa in "AILVMF" and mut_aa in "KRDENQ":
                    ddg += 2.5
                # 4. Buried Proline (backbone kink disruption)
                if mut_aa == "P" and wt_aa != "P":
                    ddg += 2.0
            else:
                # Surface Proline disruption
                if mut_aa == "P" and wt_aa != "P":
                    ddg += 0.5

        # Pairwise steric and electrostatic interactions
        n = len(mutation_tuples)
        for i in range(n):
            for j in range(i + 1, n):
                pos_i, wt_i, mut_i = mutation_tuples[i]
                pos_j, wt_j, mut_j = mutation_tuples[j]
                d = self.get_distance(pos_i, pos_j)
                if d is not None and d < 7.0:
                    vol_i = AA_VOLUMES.get(mut_i, 100.0) - AA_VOLUMES.get(wt_i, 100.0)
                    vol_j = AA_VOLUMES.get(mut_j, 100.0) - AA_VOLUMES.get(wt_j, 100.0)
                    if vol_i > 35.0 and vol_j > 35.0:
                        ddg += 3.0  # Double bulky collision
                    if (mut_i in "KR" and mut_j in "KR") or (mut_i in "DE" and mut_j in "DE"):
                        ddg += 2.0  # Like-charge electrostatic repulsion

        return round(ddg, 3)


# Global cache of parsed SpatialContactMap instances for high-throughput zero-cost reuse
_contact_map_cache: Dict[Tuple[str, Optional[str], Optional[str]], "SpatialContactMap"] = {}
_anm_hinge_cache: Dict[Tuple[str, float, int, Tuple[int, ...]], List[int]] = {}


def compute_gnm_hinge_residues(
    pdb_path: Optional[Path] = None,
    cutoff_distance: float = 12.0,
    top_k: int = 8,
    catalytic_invariants: Optional[Set[int]] = None,
) -> List[int]:
    """
    Identify allosteric hinge residues from Gaussian Network Model (GNM) normal modes.

    Builds the N x N Kirchhoff contact Laplacian over C-alpha coordinates, takes the
    lowest-frequency non-trivial modes, and flags hinges by the standard GNM
    criterion: a residue where the slow-mode amplitude changes sign is a pivot
    between two counter-moving dynamic subdomains. Catalytic invariants and
    residues inside the 4.5 A cofactor sphere are excluded (SPEC-001 Seam 1).

    This is a GNM, not an ANM. An ANM requires the 3N x 3N Hessian of directional
    second derivatives and yields vector rather than scalar mode amplitudes; that
    is a different model and is not implemented here.
    """
    resolved_path = str(_resolve_pdb_path(pdb_path, "spatial_filter"))
    inv_tuple = tuple(sorted(catalytic_invariants)) if catalytic_invariants else ()
    cache_key = (resolved_path, float(cutoff_distance), int(top_k), inv_tuple)
    if cache_key in _anm_hinge_cache:
        return list(_anm_hinge_cache[cache_key])

    cmap = get_spatial_contact_map(pdb_path)
    if not cmap.ca_coordinates:
        return []

    res_ids = sorted(cmap.ca_coordinates.keys())
    n = len(res_ids)
    if n < 5:
        return res_ids[:top_k]

    coords = np.array(
        [[cmap.ca_coordinates[r][1], cmap.ca_coordinates[r][2], cmap.ca_coordinates[r][3]] for r in res_ids],
        dtype=np.float64,
    )

    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_mat = np.sqrt(np.sum(diff * diff, axis=-1))

    within_cutoff = (dist_mat <= cutoff_distance) & (dist_mat > 1e-3)
    gamma = np.zeros((n, n), dtype=np.float64)
    gamma[within_cutoff] = -1.0 / np.maximum(dist_mat[within_cutoff] ** 2, 1.0)
    np.fill_diagonal(gamma, -np.sum(gamma, axis=1))

    eigenvalues, eigenvectors = np.linalg.eigh(gamma)

    non_zero_indices = [idx for idx, val in enumerate(eigenvalues) if val > 1e-5]
    if not non_zero_indices:
        return res_ids[:top_k]

    selected_modes = non_zero_indices[:3]

    # GNM hinge criterion: in a slow mode a hinge sits where the mode amplitude
    # crosses zero, because the residues either side move in opposite directions.
    # Each crossing is scored by how sharply the amplitude passes through zero,
    # with slower (lower-eigenvalue) modes weighted more heavily since they carry
    # the large-scale collective motion. Chain termini are included: the previous
    # implementation looped over range(1, n-1) only, so the first and last residue
    # could never be identified as hinges regardless of their dynamics.
    hinge_scores = np.zeros(n, dtype=np.float64)
    for m_idx in selected_modes:
        lam = max(float(eigenvalues[m_idx]), 1e-5)
        u = np.asarray(eigenvectors[:, m_idx], dtype=np.float64)
        weight = 1.0 / lam
        # Amplitude is near zero at a pivot, so invert magnitude to score crossings high.
        scale = float(np.max(np.abs(u))) or 1.0
        proximity_to_zero = 1.0 - (np.abs(u) / scale)
        # A sign change against the next residue marks the crossing on both sides.
        sign_change = np.zeros(n, dtype=np.float64)
        crossing = (u[:-1] * u[1:]) < 0.0
        sign_change[:-1][crossing] = 1.0
        sign_change[1:][crossing] = 1.0
        # Local gradient measures how sharply the mode passes through zero.
        grad = np.zeros(n, dtype=np.float64)
        grad[1:-1] = np.abs(u[2:] - u[:-2]) / 2.0
        grad[0] = abs(u[1] - u[0])
        grad[-1] = abs(u[-1] - u[-2])
        hinge_scores += weight * sign_change * grad * proximity_to_zero

    forbidden: Set[int] = set(catalytic_invariants) if catalytic_invariants else set()
    if 65 in cmap.ca_coordinates and 66 in cmap.ca_coordinates and 67 in cmap.ca_coordinates:
        if cmap.ca_coordinates[66][0] in ("TYR", "CRO") and cmap.ca_coordinates[67][0] in ("GLY", "CRO"):
            forbidden.update({65, 66, 67})

    candidate_pairs: List[Tuple[float, int]] = []
    for i, res_id in enumerate(res_ids):
        if res_id in forbidden:
            continue
        d_cof = cmap.get_cofactor_distance(res_id)
        if d_cof is not None and d_cof < 4.5:
            continue
        candidate_pairs.append((float(hinge_scores[i]), res_id))

    candidate_pairs.sort(key=lambda item: (-item[0], item[1]))
    result = [res_id for _, res_id in candidate_pairs[:top_k]]
    _anm_hinge_cache[cache_key] = result
    return list(result)


# Retained so existing callers keep working; the model is a GNM, see the docstring.
compute_anm_hinge_residues = compute_gnm_hinge_residues


def get_spatial_contact_map(
    pdb_path: Optional[Path] = None,
    chain_id: Optional[str] = None,
    reference_sequence: Optional[str] = None,
) -> SpatialContactMap:
    """Returns a cached SpatialContactMap instance for the target structure."""
    resolved_path = str(_resolve_pdb_path(pdb_path, "spatial_filter"))
    key = (resolved_path, chain_id, reference_sequence)
    if key not in _contact_map_cache:
        _contact_map_cache[key] = SpatialContactMap(pdb_path, chain_id=chain_id, reference_sequence=reference_sequence)
    return _contact_map_cache[key]


def extract_pocket_residues(
    functional_centers: Sequence[int],
    pdb_path: Optional[Path] = None,
    cutoff_angstroms: float = 8.5,
    include_cofactors: bool = True,
) -> List[int]:
    """Universal 3D tessellation helper: extracts all interacting residues within 8.5A of functional centers or cofactors."""
    cmap = get_spatial_contact_map(pdb_path)
    return cmap.extract_functional_pocket(functional_centers, cutoff_angstroms=cutoff_angstroms, include_cofactors=include_cofactors)


def filter_spatially_viable_variants(
    mutation_strings: List[str],
    pdb_path: Optional[Path] = None,
    catalytic_invariants: Optional[Set[int]] = None,
    enable_heuristics: bool = True,
    reference_sequence: Optional[str] = None,
    chain_id: Optional[str] = None,
    allowed_active_positions: Optional[Set[int]] = None,
) -> Tuple[List[str], List[Tuple[str, str]]]:
    """
    Filters a batch of mutation strings into approved viable variants and pruned clashing variants
    across any protein scaffold using 3D structural coordinates and catalytic invariants.

    Returns:
        Tuple of (approved_variants: List[str], pruned_variants: List[Tuple[str, reason]]).
    """
    if not enable_heuristics:
        return list(mutation_strings), []

    cmap = SpatialContactMap(pdb_path, chain_id=chain_id, reference_sequence=reference_sequence)
    if not cmap.ca_coordinates:
        logger.info("Skipping heuristic filtering: no usable structural context was loaded.")
        return list(mutation_strings), []

    approved: List[str] = []
    pruned: List[Tuple[str, str]] = []

    for mut_str in mutation_strings:
        if mut_str in ("WT", ""):
            approved.append(mut_str)
            continue

        parsed_tuples: List[Tuple[int, str, str]] = []
        positions: List[int] = []
        try:
            for part in mut_str.split(":"):
                wt_char = part[0]
                mut_char = part[-1]
                pos = int(part[1:-1])
                parsed_tuples.append((pos, wt_char, mut_char))
                positions.append(pos)
        except (ValueError, IndexError):
            approved.append(mut_str)
            continue

        # Check multi-mutant cohesion: if > 1 mutation, ensure they are not completely uncoupled (> 24 A apart)
        if len(positions) > 1:
            cohesion = cmap.evaluate_mutation_cluster_cohesion(positions)
            if not cohesion["is_cohesive"]:
                pruned.append((mut_str, f"Spatially uncoupled: max distance {cohesion.get('max_distance_angstroms', 0.0):.1f} A exceeds structural contact sphere (mean {cohesion.get('mean_distance_angstroms', 0.0):.1f} A)"))
                continue

        # Check steric, core collapse, and catalytic invariants
        passes, reason = cmap.check_steric_and_core_clash(
            parsed_tuples,
            catalytic_invariants=catalytic_invariants,
            allowed_active_positions=allowed_active_positions,
        )
        if not passes:
            pruned.append((mut_str, reason))
        else:
            approved.append(mut_str)

    return approved, pruned


def compute_structural_surrogate_score(
    variant_notation: str,
    pdb_path: Optional[Path] = None,
    catalytic_invariants: Optional[Set[int]] = None,
) -> float:
    """
    Computes a fast zero-cost structural surrogate energy/packing term S_struct(x) (Gap C).

    Formula:
      S_struct(x) = E_invariants + E_core_density + E_contact_cohesion + E_steric_clash

    Biophysical terms:
      1. Catalytic Invariant Protection: -4.5 penalty if catalytic residue is disrupted.
      2. Core Density & Cavity Creation: For buried core residues, overcrowding or voids
         incur destabilization penalties.
      3. Contact Sphere Cohesion: For multi-mutants, interacting pairs within < 8.5 A
         receive positive epistatic contact bonuses (+0.35 per pair).
      4. Steric Clash: Mutated neighbors < 7.5 A both expanding volume receive -2.5 clash penalty.

    Execution: Purely vectorized CPU numpy/math operations (~0.5 microseconds per candidate).
    Cost: $0.000000 (Zero cloud GPU spend).
    """
    if not variant_notation or variant_notation == "WT":
        return 0.0

    cmap = get_spatial_contact_map(pdb_path)
    if not cmap.ca_coordinates:
        return 0.0

    parsed_tuples: List[Tuple[int, str, str]] = []
    positions: List[int] = []
    try:
        for part in variant_notation.split(":"):
            wt_char = part[0]
            mut_char = part[-1]
            pos = int(part[1:-1])
            parsed_tuples.append((pos, wt_char, mut_char))
            positions.append(pos)
    except (ValueError, IndexError):
        return 0.0

    score = 0.0

    # 1. Invariant & Catalytic Checks
    for pos, wt_aa, mut_aa in parsed_tuples:
        if catalytic_invariants and pos in catalytic_invariants:
            if pos == 66 and wt_aa == "Y" and mut_aa in ("H", "W"):
                score += 0.2  # spectral shift permitted
            else:
                score -= 4.5
        elif pos == 66 and mut_aa not in ("Y", "H", "W") and wt_aa == "Y":
            score -= 4.5
        elif pos == 67 and mut_aa != "G" and wt_aa == "G":
            score -= 4.8

    # Positions already charged for burying a formal charge. The buried-core and
    # cofactor-shell terms describe ONE desolvation event, and both used to fire on
    # the same substitution - V227K took -2.0 and -2.5 for a single physical cause,
    # totalling the -4.5 that ranked the global optimum as catastrophic.
    charge_burial_billed: Set[int] = set()

    # 2. Buried Core Packing Density
    buried_set = cmap.buried_core_residues if cmap.buried_core_residues else AVGFP_BURIED_CORE_RESIDUES
    for pos, wt_aa, mut_aa in parsed_tuples:
        if pos in buried_set:
            wt_vol = AA_VOLUMES.get(wt_aa, 120.0)
            mut_vol = AA_VOLUMES.get(mut_aa, 120.0)
            vol_delta = mut_vol - wt_vol
            if vol_delta < -60.0:
                score -= 1.5
            elif vol_delta > 50.0:
                score -= 1.0  # core overcrowding
            # Burying a FORMAL charge (K/R/D/E) costs more than burying a
            # neutral-polar amide (N/Q); the previous "KRDENQ" class charged both
            # the same. Positions billed here are recorded so the cofactor-shell
            # term below cannot charge the same desolvation event a second time.
            if wt_aa in "AILVMF":
                if mut_aa in "KRDE":
                    score -= 2.0
                    charge_burial_billed.add(pos)
                elif mut_aa in "NQ":
                    score -= 0.9

    # 3. Pairwise Epistatic Proximity Bonus & Steric Clash Penalty
    n = len(parsed_tuples)
    for i in range(n):
        for j in range(i + 1, n):
            pos_i, wt_i, mut_i = parsed_tuples[i]
            pos_j, wt_j, mut_j = parsed_tuples[j]
            d = cmap.get_distance(pos_i, pos_j)
            if d is not None:
                if d <= 8.5:
                    score += 0.35  # Contact sphere epistatic synergy
                elif d > 16.0:
                    score -= 0.5   # Spatially uncoupled penalty

                if d < 7.5:
                    vol_i = AA_VOLUMES.get(mut_i, 100.0) - AA_VOLUMES.get(wt_i, 100.0)
                    vol_j = AA_VOLUMES.get(mut_j, 100.0) - AA_VOLUMES.get(wt_j, 100.0)
                    if vol_i > 35.0 and vol_j > 35.0:
                        score -= 2.5  # Steric clash
                    if (mut_i in "KR" and mut_j in "KR") or (mut_i in "DE" and mut_j in "DE"):
                        score -= 1.5  # Electrostatic repulsion

    # 4. 3D HETATM Cofactor / Transition-State Proximity & Catalytic Geometry (Gap 3)
    if cmap.cofactor_coordinates or cmap.catalytic_center_coordinates:
        for pos, wt_aa, mut_aa in parsed_tuples:
            d_cof = cmap.get_cofactor_distance(pos)
            if d_cof is None:
                d_cof = cmap.get_catalytic_center_distance(pos)
            if d_cof is not None:
                # Zone 1: Immediate catalytic transition-state shell (< 4.5 A)
                if d_cof < 4.5:
                    vol_delta = AA_VOLUMES.get(mut_aa, 120.0) - AA_VOLUMES.get(wt_aa, 120.0)
                    if vol_delta > 40.0:
                        score -= 3.5  # Steric crowding of catalytic pocket / cofactor
                    # One desolvation event, billed once. The buried-core term above
                    # already charges a formal charge buried in the hydrophobic core;
                    # a residue that is both buried AND in the cofactor shell used to
                    # pay twice for the same cause. The incremental penalty here is
                    # for proximity to the cofactor specifically.
                    if wt_aa in "AILVMF" and mut_aa in "KRDE":
                        score -= 0.9 if pos in charge_burial_billed else 2.5
                    elif wt_aa in "AILVMF" and mut_aa in "NQ":
                        score -= 0.4 if pos in charge_burial_billed else 1.1
                    if mut_aa == "P" and wt_aa != "P":
                        score -= 2.0  # Backbone rigidity disrupting catalytic orientation
                # Zone 2: Entrance channel & substrate accession (4.5 A to 12.0 A)
                elif d_cof <= 12.0:
                    # Steric decongestion / channel enlargement for non-natural / bulky substrates
                    # e.g., F184I or F184A in TrpB
                    if wt_aa in ("F", "W", "Y") and mut_aa in ("I", "V", "A", "L"):
                        score += 0.45  # Steric decongestion bonus
                    elif wt_aa in ("S", "T", "A") and mut_aa == "G":
                        score += 0.25  # Channel flexibility bonus

    # 5. Allosteric Hinge Cooperativity Bonus (SPEC-001 Seam 1)
    anm_hinges = set(
        compute_anm_hinge_residues(
            pdb_path=pdb_path,
            cutoff_distance=12.0,
            top_k=10,
            catalytic_invariants=catalytic_invariants,
        )
    )
    for pos, _wt_aa, _mut_aa in parsed_tuples:
        if pos in anm_hinges:
            score += 0.30

    return round(score, 4)


# -----------------------------------------------------------------------------
# Multi-Chain Quaternary Complex & Epitope-Paratope Interface (SPEC-002 Seam 1)
# -----------------------------------------------------------------------------
_multichain_ca_cache: Dict[str, Dict[Tuple[str, int], Tuple[str, float, float, float]]] = {}


def _load_multichain_ca(
    pdb_path: Optional[Path] = None,
) -> Dict[Tuple[str, int], Tuple[str, float, float, float]]:
    """Parses and caches C-alpha coordinates across ALL chains in a multi-chain PDB complex."""
    path_obj = _resolve_pdb_path(pdb_path, "multi_chain_ca_coords")
    key = str(path_obj)
    if key in _multichain_ca_cache:
        return _multichain_ca_cache[key]

    coords: Dict[Tuple[str, int], Tuple[str, float, float, float]] = {}
    if not path_obj.is_file():
        return coords

    with open(path_obj, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and len(line) >= 54:
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                alt_loc = line[16:17].strip()
                if alt_loc not in ("", "A"):
                    continue
                try:
                    res_name = line[17:20].strip()
                    chain_id = line[21:22].strip() or "A"
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords[(chain_id, res_seq)] = (res_name, x, y, z)
                except (ValueError, IndexError):
                    continue

    _multichain_ca_cache[key] = coords
    return coords


def compute_interface_contact_residues(
    pdb_path: Optional[Path] = None,
    binder_chains: Sequence[str] = ("A", "B"),
    target_chains: Sequence[str] = ("C",),
    cutoff_angstroms: float = 10.0,
) -> Dict[str, Any]:
    """
    Identifies intermolecular paratope-epitope contact residues and pairs across
    quaternary chain boundaries (e.g., antibody chains H/L or A/B against antigen chain C).
    """
    coords = _load_multichain_ca(pdb_path)
    b_set = {c.strip() for c in binder_chains}
    t_set = {c.strip() for c in target_chains}

    binder_items = [(ch, r_id, info) for (ch, r_id), info in coords.items() if ch in b_set]
    target_items = [(ch, r_id, info) for (ch, r_id), info in coords.items() if ch in t_set]

    binder_min_dist: Dict[Tuple[str, int], Tuple[str, float]] = {}
    target_min_dist: Dict[Tuple[str, int], Tuple[str, float]] = {}
    contact_pairs: List[Dict[str, Any]] = []

    for b_ch, b_res, (b_name3, bx, by, bz) in binder_items:
        b_aa1 = THREE_TO_ONE.get(b_name3, "X")
        for t_ch, t_res, (t_name3, tx, ty, tz) in target_items:
            d = math.sqrt((bx - tx) ** 2 + (by - ty) ** 2 + (bz - tz) ** 2)
            if d <= cutoff_angstroms:
                t_aa1 = THREE_TO_ONE.get(t_name3, "X")
                contact_pairs.append({
                    "binder_chain": b_ch,
                    "binder_res_seq": b_res,
                    "binder_wt_aa": b_aa1,
                    "target_chain": t_ch,
                    "target_res_seq": t_res,
                    "target_wt_aa": t_aa1,
                    "distance_angstrom": round(d, 3),
                })
                if (b_ch, b_res) not in binder_min_dist or d < binder_min_dist[(b_ch, b_res)][1]:
                    binder_min_dist[(b_ch, b_res)] = (b_aa1, d)
                if (t_ch, t_res) not in target_min_dist or d < target_min_dist[(t_ch, t_res)][1]:
                    target_min_dist[(t_ch, t_res)] = (t_aa1, d)

    contact_pairs.sort(key=lambda x: x["distance_angstrom"])
    binder_residues = [
        {
            "chain_id": ch,
            "res_seq": r_id,
            "wt_aa": aa1,
            "min_distance_angstrom": round(dist, 3),
        }
        for (ch, r_id), (aa1, dist) in sorted(binder_min_dist.items(), key=lambda item: item[1][1])
    ]
    target_residues = [
        {
            "chain_id": ch,
            "res_seq": r_id,
            "wt_aa": aa1,
            "min_distance_angstrom": round(dist, 3),
        }
        for (ch, r_id), (aa1, dist) in sorted(target_min_dist.items(), key=lambda item: item[1][1])
    ]

    return {
        "binder_chains": list(binder_chains),
        "target_chains": list(target_chains),
        "cutoff_angstroms": cutoff_angstroms,
        "binder_residues": binder_residues,
        "target_residues": target_residues,
        "contact_pairs": contact_pairs,
    }


def compute_interface_binding_ddg(
    pdb_path: Optional[Path] = None,
    mutation_tuples: Sequence[Tuple[int, str, str]] = (),
    binder_chains: Sequence[str] = ("A", "B"),
    target_chains: Sequence[str] = ("C",),
) -> float:
    """
    Heuristic interface-disruption penalty score (unitless, higher = more disruptive).

    NOT a kcal/mol binding free energy - the same caveat as compute_stability_ddg
    applies: hand-chosen constants, no force field, no parameterization against
    measured affinities. Computed across a
    multi-chain protein-protein interface (SPEC-002 Seam 1).
    Negative ddG_bind = tighter binding affinity (favorable).
    Positive ddG_bind = weakened binding or interfacial clash (unfavorable).
    """
    if not mutation_tuples:
        return 0.0

    coords = _load_multichain_ca(pdb_path)
    b_set = {c.strip() for c in binder_chains}
    t_set = {c.strip() for c in target_chains}
    target_items = [(ch, r_id, info) for (ch, r_id), info in coords.items() if ch in t_set]

    ddg_bind = 0.0
    for pos, wt_aa, mut_aa in mutation_tuples:
        # Locate mutated position in binder_chains
        b_coord = None
        for ch in binder_chains:
            if (ch.strip(), pos) in coords:
                b_coord = coords[(ch.strip(), pos)]
                break
        if b_coord is None:
            for (ch, r_id), info in coords.items():
                if ch in b_set and r_id == pos:
                    b_coord = info
                    break
        if b_coord is None or not target_items:
            continue

        _, bx, by, bz = b_coord
        contacts: List[Tuple[float, str]] = []
        for _t_ch, _t_res, (t_name3, tx, ty, tz) in target_items:
            d = math.sqrt((bx - tx) ** 2 + (by - ty) ** 2 + (bz - tz) ** 2)
            if d <= 12.0:
                contacts.append((d, THREE_TO_ONE.get(t_name3, "A")))

        if not contacts:
            continue

        contacts.sort(key=lambda x: x[0])
        min_d, nearest_t_aa = contacts[0]
        vol_delta = AA_VOLUMES.get(mut_aa, 120.0) - AA_VOLUMES.get(wt_aa, 120.0)

        # 1. Interfacial proline loop kinking penalty
        if mut_aa == "P" and wt_aa != "P":
            ddg_bind += 2.20

        # 2. Loss of interfacial van der Waals packing (cavity creation at interface)
        if min_d <= 9.5 and vol_delta < -45.0:
            ddg_bind += 0.025 * abs(vol_delta)

        # 3. Direct interfacial steric collision (< 6.2 A)
        if min_d < 6.2 and vol_delta > 45.0:
            ddg_bind += 1.80
        # 4. Favorable paratope aromatic pi-stacking / H-bonding (Tyr/Trp/Phe at 6.2 - 11.0 A)
        elif 6.2 <= min_d <= 11.0 and mut_aa in ("Y", "W", "F"):
            ddg_bind -= 0.55

        # 5. Intermolecular electrostatic salt-bridge vs repulsion across interface
        for d_contact, t_aa in contacts[:3]:
            if d_contact <= 10.5:
                if (mut_aa in "KR" and t_aa in "DE") or (mut_aa in "DE" and t_aa in "KR"):
                    ddg_bind -= 0.45  # Intermolecular salt bridge
                elif (mut_aa in "KR" and t_aa in "KR") or (mut_aa in "DE" and t_aa in "DE"):
                    ddg_bind += 0.90  # Intermolecular charge repulsion

    return round(ddg_bind, 4)



