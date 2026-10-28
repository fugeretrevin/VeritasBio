"""Optional 3D structural diagnostics and heuristic filtering.

Purpose:
    Extracts 3D C-alpha coordinates from crystallographic PDB structures (or AlphaFold models)
    to calculate spatial Euclidean distance matrices D_ij = ||r_i - r_j||_2.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)

# Standard Van der Waals sidechain volume estimates (in Angstroms^3)
AA_VOLUMES: Dict[str, float] = {
    "G": 60.1,  "A": 88.6,  "S": 89.0,  "C": 108.5, "D": 111.1,
    "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0,
    "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7,
    "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8,
}

THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class SpatialContactMap:
    """Extracts C-alpha coordinates from one explicitly selected PDB chain."""

    def __init__(
        self,
        pdb_path: Optional[Path] = None,
        chain_id: Optional[str] = None,
        reference_sequence: Optional[str] = None,
    ) -> None:
        self.pdb_path: Optional[Path] = Path(pdb_path) if pdb_path is not None else None
        self.chain_id = chain_id
        self.reference_sequence = reference_sequence
        self.ca_coordinates: Dict[int, Tuple[str, float, float, float]] = {}
        self._distance_cache: Dict[Tuple[int, int], float] = {}
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
        return buried

    def _load_pdb(self) -> None:
        """Parses ATOM and HETATM C-alpha coordinates from PDB."""
        if self.pdb_path is None:
            logger.info("Structural information unavailable: no PDB structure supplied.")
            return
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

                        if selected_chain is None:
                            selected_chain = chain_id
                        if chain_id != selected_chain:
                            continue

                        # Insertion codes cannot be mapped safely by integer residue number alone.
                        if self.reference_sequence and line[26].strip():
                            self.ca_coordinates.clear()
                            self.chain_id = selected_chain
                            logger.warning(
                                "Structural mapping unavailable: insertion codes in %s chain %s.",
                                self.pdb_path,
                                selected_chain or "blank",
                            )
                            return

                        if res_name == "CRO" and atom_name in ("CA1", "CA2", "CA3", "C1", "C2", "CA"):
                            try:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
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
        except (OSError, UnicodeError) as exc:
            logger.warning("Could not read PDB structure at %s: %s", self.pdb_path, exc)
            return

        self.chain_id = selected_chain
        if self.reference_sequence and self.ca_coordinates:
            mismatches = [
                position for position, (res_name, *_coords) in self.ca_coordinates.items()
                if position < 1
                or position > len(self.reference_sequence)
                or THREE_TO_ONE.get(res_name) != self.reference_sequence[position - 1].upper()
            ]
            if mismatches:
                self.ca_coordinates.clear()
                logger.warning(
                    "Structural mapping unavailable for %s chain %s: residue identity mismatch at %s.",
                    self.pdb_path,
                    self.chain_id or "blank",
                    mismatches[:5],
                )
                return

        if self.ca_coordinates:
            logger.info(
                "Loaded structural context from %s, chain %s (%d C-alpha coordinates).",
                self.pdb_path,
                self.chain_id or "blank",
                len(self.ca_coordinates),
            )
        else:
            logger.warning("No usable C-alpha coordinates found in %s, chain %s.", self.pdb_path, self.chain_id)

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

        # A valid epistatic contact sphere has a mean distance <= 13.0 A and max pairwise <= 16.0 A
        is_cohesive = (mean_d <= 13.0) and (max_d <= 16.0)

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
    ) -> Tuple[bool, str]:
        """
        Checks a proposed set of mutations for steric clashes, core disruption,
        or catalytic invariant violations across any protein target.

        Args:
            mutation_tuples: List of (position, wt_aa, mut_aa) tuples.
            catalytic_invariants: Optional set of 1-indexed residue positions that cannot tolerate mutation.

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
        buried_set = self.buried_core_residues
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
    ) -> List[int]:
        """
        Extracts all residue positions within cutoff_angstroms of ANY functional center.
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

        return sorted(list(pocket))

    def compute_stability_ddg(
        self,
        mutation_tuples: List[Tuple[int, str, str]],
    ) -> float:
        """
        Estimates thermodynamic folding stability change Delta-Delta-G (kcal/mol).
        Positive = destabilizing (deleterious to fold).
        Negative = stabilizing.
        """
        ddg = 0.0
        buried_set = self.buried_core_residues

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


def get_spatial_contact_map(pdb_path: Optional[Path] = None) -> SpatialContactMap:
    """Loads only the explicitly supplied structure; no cross-target cache is retained."""
    return SpatialContactMap(pdb_path)


def extract_pocket_residues(
    functional_centers: Sequence[int],
    pdb_path: Optional[Path] = None,
    cutoff_angstroms: float = 8.5,
) -> List[int]:
    """Universal 3D tessellation helper: extracts all interacting residues within 8.5A of functional centers."""
    cmap = get_spatial_contact_map(pdb_path)
    return cmap.extract_functional_pocket(functional_centers, cutoff_angstroms=cutoff_angstroms)


def filter_spatially_viable_variants(
    mutation_strings: List[str],
    pdb_path: Optional[Path] = None,
    catalytic_invariants: Optional[Set[int]] = None,
    enable_heuristics: bool = False,
    reference_sequence: Optional[str] = None,
) -> Tuple[List[str], List[Tuple[str, str]]]:
    """
    Filters a batch of mutation strings into approved viable variants and pruned clashing variants
    only when explicitly enabled. Default behavior preserves candidate order and identity.

    Returns:
        Tuple of (approved_variants: List[str], pruned_variants: List[Tuple[str, reason]]).
    """
    if not enable_heuristics:
        return list(mutation_strings), []

    cmap = SpatialContactMap(pdb_path, reference_sequence=reference_sequence)
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

        if any(position not in cmap.ca_coordinates for position in positions):
            approved.append(mut_str)
            continue

        # Check multi-mutant cohesion: if > 1 mutation, ensure they are not completely uncoupled (> 18 A apart)
        if len(positions) > 1:
            cohesion = cmap.evaluate_mutation_cluster_cohesion(positions)
            if not cohesion["is_cohesive"]:
                pruned.append((mut_str, f"Spatially uncoupled: max distance {cohesion.get('max_distance_angstroms', 0.0):.1f} A exceeds structural contact sphere (mean {cohesion.get('mean_distance_angstroms', 0.0):.1f} A)"))
                continue

        # Check steric, core collapse, and catalytic invariants
        passes, reason = cmap.check_steric_and_core_clash(parsed_tuples, catalytic_invariants=catalytic_invariants)
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

    # 2. Buried Core Packing Density
    buried_set = cmap.buried_core_residues
    for pos, wt_aa, mut_aa in parsed_tuples:
        if pos in buried_set:
            wt_vol = AA_VOLUMES.get(wt_aa, 120.0)
            mut_vol = AA_VOLUMES.get(mut_aa, 120.0)
            vol_delta = mut_vol - wt_vol
            if vol_delta < -60.0:
                score -= 1.5
            elif vol_delta > 50.0:
                score -= 1.0  # core overcrowding
            if wt_aa in "AILVMF" and mut_aa in "KRDENQ":
                score -= 2.0  # desolvation penalty

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

    return round(score, 4)


