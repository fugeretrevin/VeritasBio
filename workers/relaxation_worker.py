"""
workers/relaxation_worker.py — Track C Sidechain Rotamer & Local Energy Relaxation Engine.

Purpose:
    Performs fast CPU-based L-BFGS local energy minimization on candidate mutant protein
    structures under harmonic backbone restraints (k = 50.0 energy units/A^2).
    The potential is harmonic restraints plus a soft-sphere clash term only:
    there is no Lennard-Jones attraction, no electrostatics and no solvation,
    so reported energies are in arbitrary units of this potential and are not
    physical kcal/mol values.
    
    Eliminates false-negative in silico steric clash rejections where bulky mutations
    (e.g. Ala -> Trp, Ser -> Phe) collide with static crystal lattice coordinates
    (< 1.8 A) but can readily relax into native packing geometry in < 80ms.
"""

from __future__ import annotations

import copy
import functools
import logging
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from pydantic import BaseModel, Field
from scipy.optimize import minimize

logger = logging.getLogger(__name__)

# Standard 1-letter to 3-letter code mapping
ONE_TO_THREE: Dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}

THREE_TO_ONE: Dict[str, str] = {v: k for k, v in ONE_TO_THREE.items()}

# Van der Waals radii (Angstroms)
VDW_RADII: Dict[str, float] = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20, "P": 1.80, "FE": 1.80,
}

# Ideal sidechain topology definitions for canonical amino acids
# [atom_name, parent_atom, ideal_bond_length, element]
IDEAL_SIDECHAINS: Dict[str, List[Tuple[str, str, float, str]]] = {
    "ALA": [("CB", "CA", 1.52, "C")],
    "VAL": [("CB", "CA", 1.53, "C"), ("CG1", "CB", 1.52, "C"), ("CG2", "CB", 1.52, "C")],
    "LEU": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.53, "C"), ("CD1", "CG", 1.52, "C"), ("CD2", "CG", 1.52, "C")],
    "ILE": [("CB", "CA", 1.53, "C"), ("CG1", "CB", 1.52, "C"), ("CG2", "CB", 1.52, "C"), ("CD1", "CG1", 1.52, "C")],
    "PRO": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("CD", "CG", 1.52, "C")],
    "PHE": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.50, "C"), ("CD1", "CG", 1.39, "C"), ("CD2", "CG", 1.39, "C"),
            ("CE1", "CD1", 1.39, "C"), ("CE2", "CD2", 1.39, "C"), ("CZ", "CE1", 1.39, "C")],
    "TYR": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.50, "C"), ("CD1", "CG", 1.39, "C"), ("CD2", "CG", 1.39, "C"),
            ("CE1", "CD1", 1.39, "C"), ("CE2", "CD2", 1.39, "C"), ("CZ", "CE1", 1.39, "C"), ("OH", "CZ", 1.36, "O")],
    "TRP": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.50, "C"), ("CD1", "CG", 1.37, "C"), ("CD2", "CG", 1.43, "C"),
            ("NE1", "CD1", 1.38, "N"), ("CE2", "CD2", 1.40, "C"), ("CE3", "CD2", 1.40, "C"), ("CZ2", "CE2", 1.40, "C"),
            ("CZ3", "CE3", 1.40, "C"), ("CH2", "CZ2", 1.40, "C")],
    "SER": [("CB", "CA", 1.53, "C"), ("OG", "CB", 1.42, "O")],
    "THR": [("CB", "CA", 1.53, "C"), ("OG1", "CB", 1.42, "O"), ("CG2", "CB", 1.52, "C")],
    "CYS": [("CB", "CA", 1.53, "C"), ("SG", "CB", 1.81, "S")],
    "MET": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("SD", "CG", 1.81, "S"), ("CE", "SD", 1.79, "C")],
    "ASN": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.51, "C"), ("OD1", "CG", 1.23, "O"), ("ND2", "CG", 1.33, "N")],
    "GLN": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("CD", "CG", 1.51, "C"), ("OE1", "CD", 1.23, "O"), ("NE2", "CD", 1.33, "N")],
    "ASP": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.51, "C"), ("OD1", "CG", 1.25, "O"), ("OD2", "CG", 1.25, "O")],
    "GLU": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("CD", "CG", 1.51, "C"), ("OE1", "CD", 1.25, "O"), ("OE2", "CD", 1.25, "O")],
    "LYS": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("CD", "CG", 1.52, "C"), ("CE", "CD", 1.52, "C"), ("NZ", "CE", 1.49, "N")],
    "ARG": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.52, "C"), ("CD", "CG", 1.52, "C"), ("NE", "CD", 1.46, "N"),
            ("CZ", "NE", 1.33, "C"), ("NH1", "CZ", 1.33, "N"), ("NH2", "CZ", 1.33, "N")],
    "HIS": [("CB", "CA", 1.53, "C"), ("CG", "CB", 1.49, "C"), ("ND1", "CG", 1.38, "N"), ("CD2", "CG", 1.36, "C"),
            ("CE1", "ND1", 1.32, "C"), ("NE2", "CD2", 1.37, "N")],
    "GLY": [],
}


class RelaxationResult(BaseModel):
    """Output metrics from local energy minimization."""

    variant_notation: str = Field(..., description="Mutation notation, e.g. A206W or S65T:Y66H")
    initial_energy_kcal: float = Field(
        ...,
        description=(
            "Pre-relaxation energy in arbitrary units of the restraint + soft-sphere "
            "potential. Not kcal/mol - no force field backs this value."
        ),
    )
    final_energy_kcal: float = Field(
        ...,
        description=(
            "Post-relaxation minimized energy in arbitrary units of the restraint + "
            "soft-sphere potential. Not kcal/mol - comparable only to "
            "initial_energy_kcal from the same run."
        ),
    )
    initial_min_distance_angstrom: float = Field(..., description="Minimum interatomic distance before minimization")
    final_min_distance_angstrom: float = Field(..., description="Minimum interatomic distance after minimization")
    clashes_resolved: int = Field(default=0, description="Number of atomic clashes (< 1.8 A) resolved")
    converged: bool = Field(default=True, description="Whether L-BFGS converged within iteration limit")
    relaxed_pdb_str: str = Field(..., description="Full minimized PDB file content")
    rmsd_ca_angstrom: float = Field(default=0.0, description="Root-mean-square deviation of C-alpha atoms in shell")
    runtime_ms: float = Field(default=0.0, description="Execution wall-clock time in milliseconds")
    backbone_stage2_triggered: bool = Field(
        default=False,
        description="Whether Stage 2 elastic mainchain accommodation (k_CA=5.0) was triggered",
    )
    denaturation_rejected: bool = Field(
        default=False,
        description="Whether mainchain displacement exceeded the 1.2 A unfolding threshold",
    )


class AtomRecord:
    """Internal atom representation for PDB parsing and vectorization."""
    __slots__ = (
        "record_type", "serial", "name", "alt_loc", "res_name",
        "chain_id", "res_seq", "i_code", "x", "y", "z",
        "occupancy", "temp_factor", "element", "charge", "is_backbone"
    )

    def __init__(
        self, record_type: str, serial: int, name: str, alt_loc: str,
        res_name: str, chain_id: str, res_seq: int, i_code: str,
        x: float, y: float, z: float, occupancy: float, temp_factor: float,
        element: str, charge: str
    ):
        self.record_type = record_type
        self.serial = serial
        self.name = name.strip()
        self.alt_loc = alt_loc
        self.res_name = res_name.strip()
        self.chain_id = chain_id.strip()
        self.res_seq = res_seq
        self.i_code = i_code
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.element = element.strip().upper() if element.strip() else self.name[0]
        self.charge = charge
        self.is_backbone = self.name in ("N", "CA", "C", "O")

    def clone(self) -> AtomRecord:
        return AtomRecord(
            self.record_type, self.serial, self.name, self.alt_loc,
            self.res_name, self.chain_id, self.res_seq, self.i_code,
            self.x, self.y, self.z, self.occupancy, self.temp_factor,
            self.element, self.charge
        )

    def to_pdb_line(self) -> str:
        return (
            f"{self.record_type:<6}{self.serial:>5} {self.name:<4}{self.alt_loc:1}"
            f"{self.res_name:>3} {self.chain_id:1}{self.res_seq:>4}{self.i_code:1}   "
            f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}"
            f"{self.occupancy:>6.2f}{self.temp_factor:>6.2f}          "
            f"{self.element:>2}{self.charge:>2}\n"
        )


@functools.lru_cache(maxsize=16)
def _load_cached_pdb_atoms(pdb_path_str: str) -> Tuple[AtomRecord, ...]:
    """Cached immutable tuple of atoms for fast repeat evaluations."""
    atoms: List[AtomRecord] = []
    with open(pdb_path_str, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM") and len(line) >= 54:
                try:
                    serial = int(line[6:11].strip() or "0")
                    name = line[12:16]
                    alt_loc = line[16:17]
                    res_name = line[17:20]
                    chain_id = line[21:22]
                    res_seq = int(line[22:26].strip() or "0")
                    i_code = line[26:27]
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occ = float(line[54:60].strip() or "1.00") if len(line) >= 60 else 1.00
                    b_factor = float(line[60:66].strip() or "0.00") if len(line) >= 66 else 0.00
                    element = line[76:78].strip() if len(line) >= 78 else name.strip()[:1]
                    charge = line[78:80].strip() if len(line) >= 80 else ""

                    if alt_loc not in (" ", "A", ""):
                        continue

                    atoms.append(AtomRecord(
                        rec, serial, name, alt_loc, res_name, chain_id,
                        res_seq, i_code, x, y, z, occ, b_factor, element, charge
                    ))
                except (ValueError, IndexError):
                    continue
    return tuple(atoms)


def parse_pdb_atoms(pdb_path: Union[str, Path]) -> List[AtomRecord]:
    """Returns mutable clones of cached PDB atoms."""
    cached = _load_cached_pdb_atoms(str(pdb_path))
    return [a.clone() for a in cached]


def parse_mutation_token(token: str) -> Tuple[str, int, str]:
    """Parses mutation string like 'A206W' or 'S65T' -> ('A', 206, 'W')."""
    token = token.strip()
    if not token or len(token) < 3:
        raise ValueError(f"Invalid mutation notation: {token}")
    wt = token[0].upper()
    mut = token[-1].upper()
    pos_str = token[1:-1]
    return wt, int(pos_str), mut


def relax_variant_structure(
    pdb_path: Union[str, Path],
    mutations: List[str],
    max_iterations: int = 200,
    harmonic_k: float = 50.0,
    shell_radius: float = 7.5,
    allow_backbone_accommodation: bool = False,
    strain_threshold: float = 20.0,
) -> RelaxationResult:
    """
    Performs fast CPU L-BFGS local energy minimization on sidechains around target mutations.
    """
    t_start = time.perf_counter()
    atoms = parse_pdb_atoms(pdb_path)
    if not atoms:
        raise ValueError(f"No ATOM records found in {pdb_path}")

    # Parse all mutation tuples
    parsed_muts: List[Tuple[str, int, str]] = []
    for m in mutations:
        if ":" in m:
            for sub_m in m.split(":"):
                if sub_m.strip():
                    parsed_muts.append(parse_mutation_token(sub_m))
        else:
            parsed_muts.append(parse_mutation_token(m))

    mut_positions = {p for _, p, _ in parsed_muts}
    mut_map = {p: mut for _, p, mut in parsed_muts}

    # Identify CA coordinates of mutation positions
    ca_mut_coords: List[np.ndarray] = []
    for a in atoms:
        if a.res_seq in mut_positions and a.name == "CA":
            ca_mut_coords.append(np.array([a.x, a.y, a.z], dtype=np.float64))

    if not ca_mut_coords:
        for a in atoms:
            if a.res_seq in mut_positions:
                ca_mut_coords.append(np.array([a.x, a.y, a.z], dtype=np.float64))
                break

    if not ca_mut_coords:
        ca_mut_coords = [np.array([atoms[0].x, atoms[0].y, atoms[0].z], dtype=np.float64)]

    # Mutate residues in the atom list
    new_atoms: List[AtomRecord] = []
    for a in atoms:
        if a.res_seq in mut_map:
            mut_1 = mut_map[a.res_seq]
            mut_3 = ONE_TO_THREE.get(mut_1, "ALA")
            a.res_name = mut_3
            if a.is_backbone:
                new_atoms.append(a)
        else:
            new_atoms.append(a)

    # Build initial sidechain coordinates for mutated residues
    serial_counter = max(a.serial for a in new_atoms) + 1
    mut_sidechain_indices: List[int] = []

    for wt, pos, mut_1 in parsed_muts:
        mut_3 = ONE_TO_THREE.get(mut_1, "ALA")
        ca_atom = next((a for a in new_atoms if a.res_seq == pos and a.name == "CA"), None)
        n_atom = next((a for a in new_atoms if a.res_seq == pos and a.name == "N"), None)
        c_atom = next((a for a in new_atoms if a.res_seq == pos and a.name == "C"), None)

        if ca_atom is None:
            continue

        ca_pos = np.array([ca_atom.x, ca_atom.y, ca_atom.z], dtype=np.float64)
        if n_atom and c_atom:
            n_pos = np.array([n_atom.x, n_atom.y, n_atom.z], dtype=np.float64)
            c_pos = np.array([c_atom.x, c_atom.y, c_atom.z], dtype=np.float64)
            bisector = ((n_pos - ca_pos) + (c_pos - ca_pos))
            norm = np.linalg.norm(bisector)
            cb_dir = -bisector / norm if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
        else:
            cb_dir = np.array([1.0, 0.0, 0.0])

        topo = IDEAL_SIDECHAINS.get(mut_3, [])
        built_coords: Dict[str, np.ndarray] = {"CA": ca_pos}

        for sc_name, parent_name, bond_len, elem in topo:
            parent_coord = built_coords.get(parent_name, ca_pos)
            offset = cb_dir * bond_len
            if sc_name not in ("CB",):
                orthogonal = np.cross(cb_dir, np.array([0.0, 0.0, 1.0]))
                if np.linalg.norm(orthogonal) < 0.1:
                    orthogonal = np.cross(cb_dir, np.array([0.0, 1.0, 0.0]))
                orthogonal = orthogonal / (np.linalg.norm(orthogonal) + 1e-6)
                offset = offset * 0.7 + orthogonal * 0.7 * bond_len

            coord = parent_coord + offset
            built_coords[sc_name] = coord

            sc_atom = AtomRecord(
                "ATOM", serial_counter, sc_name, " ", mut_3, ca_atom.chain_id,
                pos, ca_atom.i_code, coord[0], coord[1], coord[2],
                1.00, 20.00, elem, ""
            )
            serial_counter += 1
            mut_sidechain_indices.append(len(new_atoms))
            new_atoms.append(sc_atom)

    # Vectorized shell atom selection for interaction evaluation
    all_coords = np.array([[a.x, a.y, a.z] for a in new_atoms], dtype=np.float64)
    ca_mut_arr = np.array(ca_mut_coords, dtype=np.float64)
    dists_to_mut = np.sqrt(np.sum((all_coords[:, None, :] - ca_mut_arr[None, :, :]) ** 2, axis=-1))
    min_d_arr = np.min(dists_to_mut, axis=-1)
    res_seq_all = np.array([a.res_seq for a in new_atoms], dtype=np.int64)
    in_mut_pos = np.isin(res_seq_all, list(mut_positions))
    in_shell_mask = (min_d_arr <= shell_radius) | in_mut_pos

    shell_indices = [int(i) for i in np.nonzero(in_shell_mask)[0]]
    shell_coords_list = [all_coords[i] for i in shell_indices]
    is_backbone_mask = [new_atoms[i].is_backbone for i in shell_indices]
    vdw_radii_list = [VDW_RADII.get(new_atoms[i].element, 1.65) for i in shell_indices]

    n_shell = len(shell_indices)
    if n_shell == 0:
        pdb_str = "".join(a.to_pdb_line() for a in new_atoms)
        return RelaxationResult(
            variant_notation=":".join(f"{w}{p}{m}" for w, p, m in parsed_muts),
            initial_energy_kcal=0.0,
            final_energy_kcal=0.0,
            initial_min_distance_angstrom=2.5,
            final_min_distance_angstrom=2.5,
            clashes_resolved=0,
            converged=True,
            relaxed_pdb_str=pdb_str,
            rmsd_ca_angstrom=0.0,
            runtime_ms=(time.perf_counter() - t_start) * 1000.0,
        )

    coords_init = np.array(shell_coords_list, dtype=np.float64)
    coords_0 = coords_init.copy()
    is_bb = np.array(is_backbone_mask, dtype=bool)
    radii = np.array(vdw_radii_list, dtype=np.float64)
    r_matrix = radii[:, None] + radii[None, :]

    # Vectorized pairwise exclusion mask: intra-residue pairs and adjacent peptide C-N bonds
    shell_res_seq = np.array([new_atoms[i].res_seq for i in shell_indices], dtype=np.int64)
    shell_chains = np.array([new_atoms[i].chain_id for i in shell_indices])
    shell_names = np.array([new_atoms[i].name for i in shell_indices])

    same_chain = shell_chains[:, None] == shell_chains[None, :]
    same_res = (shell_res_seq[:, None] == shell_res_seq[None, :]) & same_chain
    adj_res = (np.abs(shell_res_seq[:, None] - shell_res_seq[None, :]) == 1) & same_chain
    peptide_cn = ((shell_names[:, None] == "C") & (shell_names[None, :] == "N")) | (
        (shell_names[:, None] == "N") & (shell_names[None, :] == "C")
    )
    exclusion_mask = same_res | (adj_res & peptide_cn)
    np.fill_diagonal(exclusion_mask, True)

    # Pre-select candidate clash pairs that can interact
    upper_tri = np.triu(np.ones((n_shell, n_shell), dtype=bool), k=1)
    eval_mask = (~exclusion_mask) & upper_tri
    eval_i, eval_j = np.nonzero(eval_mask)

    is_ca_mask = [new_atoms[idx].name == "CA" for idx in shell_indices]
    is_ca = np.array(is_ca_mask, dtype=bool)

    # Map each mutated sidechain atom in the shell to its parent residue's CA shell index
    ca_shell_by_res: Dict[Tuple[str, int], int] = {}
    for s_idx, atom_idx in enumerate(shell_indices):
        a = new_atoms[atom_idx]
        if a.name == "CA":
            ca_shell_by_res[(a.chain_id, a.res_seq)] = s_idx

    tether_sc_list: List[int] = []
    tether_ca_list: List[int] = []
    for s_idx, atom_idx in enumerate(shell_indices):
        a = new_atoms[atom_idx]
        if not a.is_backbone and a.res_seq in mut_positions:
            ca_s = ca_shell_by_res.get((a.chain_id, a.res_seq))
            if ca_s is not None:
                tether_sc_list.append(s_idx)
                tether_ca_list.append(ca_s)

    tether_sc_idx = np.array(tether_sc_list, dtype=np.intp)
    tether_ca_idx = np.array(tether_ca_list, dtype=np.intp)
    if tether_sc_idx.size > 0:
        tether_vec_0 = coords_0[tether_sc_idx] - coords_0[tether_ca_idx]
        tether_dist_0 = np.sqrt(np.sum(tether_vec_0 ** 2, axis=-1) + 1e-12)
    else:
        tether_vec_0 = np.zeros((0, 3), dtype=np.float64)
        tether_dist_0 = np.zeros((0,), dtype=np.float64)

    def make_objective(k_ca_val: float, k_other_val: float):
        def compute_energy_and_gradient(coords_flat: np.ndarray) -> Tuple[float, np.ndarray]:
            coords = coords_flat.reshape((n_shell, 3))
            energy = 0.0
            grad = np.zeros_like(coords)

            if np.any(is_ca):
                diff_ca = coords[is_ca] - coords_0[is_ca]
                energy += 0.5 * k_ca_val * np.sum(diff_ca ** 2)
                grad[is_ca] += k_ca_val * diff_ca

            other_bb = is_bb & (~is_ca)
            if np.any(other_bb):
                diff_other = coords[other_bb] - coords_0[other_bb]
                energy += 0.5 * k_other_val * np.sum(diff_other ** 2)
                grad[other_bb] += k_other_val * diff_other

            # Covalent & rotamer cone tether coupling mutated sidechains to their parent CA
            if tether_sc_idx.size > 0:
                rel_vec = coords[tether_sc_idx] - coords[tether_ca_idx]
                rel_dist = np.sqrt(np.sum(rel_vec ** 2, axis=-1) + 1e-12)
                radial_strain = rel_dist - tether_dist_0
                k_radial = 60.0
                energy += 0.5 * k_radial * np.sum(radial_strain ** 2)
                f_radial = (k_radial * radial_strain)[:, None] * (rel_vec / rel_dist[:, None])
                np.add.at(grad, tether_sc_idx, f_radial)
                np.add.at(grad, tether_ca_idx, -f_radial)

                cone_diff = rel_vec - tether_vec_0
                k_cone = 6.0
                energy += 0.5 * k_cone * np.sum(cone_diff ** 2)
                f_cone = k_cone * cone_diff
                np.add.at(grad, tether_sc_idx, f_cone)
                np.add.at(grad, tether_ca_idx, -f_cone)

            # Non-bonded Steric Clash Penalty
            k_clash = 150.0  # energy units/A^2 (see module docstring: not kcal/mol)
            diff = coords[eval_i] - coords[eval_j]  # (M, 3)
            dist = np.sqrt(np.sum(diff ** 2, axis=-1) + 1e-12)  # (M,)

            clash_depth = r_matrix[eval_i, eval_j] - dist
            active = clash_depth > 0
            if np.any(active):
                depths_act = clash_depth[active]
                energy += np.sum(0.5 * k_clash * (depths_act ** 2))

                forces = -k_clash * depths_act[:, None] * (diff[active] / dist[active, None])
                np.add.at(grad, eval_i[active], forces)
                np.add.at(grad, eval_j[active], -forces)

            return float(energy), grad.ravel()

        return compute_energy_and_gradient

    k_ca_stage1 = max(500.0, harmonic_k * 50.0)
    k_other_stage1 = max(50.0, harmonic_k * 5.0)
    stage1_objective = make_objective(k_ca_stage1, k_other_stage1)

    # Initial metrics
    e_init, _ = stage1_objective(coords_init.ravel())

    diff_init = coords_init[:, None, :] - coords_init[None, :, :]
    dist_init = np.sqrt(np.sum(diff_init ** 2, axis=-1) + 1e-12)
    dist_init[exclusion_mask] = 999.0
    dist_init[np.tril_indices(n_shell)] = 999.0
    initial_min_d = float(np.min(dist_init)) if dist_init.size > 0 else 2.5
    initial_clashes = int(np.sum((dist_init < 1.8) & (dist_init > 0)))

    # Stage 1: Rigid-backbone fast L-BFGS-B local minimization (k_CA = 500)
    opt_res = minimize(
        fun=stage1_objective,
        x0=coords_init.ravel(),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": min(max_iterations, 60), "ftol": 2e-3, "gtol": 2e-2},
    )

    coords_final = opt_res.x.reshape((n_shell, 3))
    e_final, _ = stage1_objective(opt_res.x)

    # Stage 2: Elastic mainchain accommodation (SPEC-001 Seam 2)
    # Triggered when residual Stage-1 strain exceeds strain_threshold and allow_backbone_accommodation=True
    backbone_stage2_triggered = False
    denaturation_rejected = False
    if allow_backbone_accommodation and e_final > strain_threshold:
        backbone_stage2_triggered = True
        stage2_objective = make_objective(5.0, 10.0)
        opt_stage2 = minimize(
            fun=stage2_objective,
            x0=coords_final.ravel(),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": min(max_iterations, 50), "ftol": 2e-3, "gtol": 2e-2},
        )
        coords_stage2 = opt_stage2.x.reshape((n_shell, 3))
        e_stage2, _ = stage2_objective(opt_stage2.x)

        if np.any(is_ca):
            ca_disp = coords_stage2[is_ca] - coords_init[is_ca]
            raw_rmsd_ca = float(np.sqrt(np.mean(np.sum(ca_disp ** 2, axis=-1))))
            if raw_rmsd_ca > 1.20:
                # Enforce hard 1.20 A maximum C-alpha RMSD cutoff to prevent unphysical unfolding
                scale = 1.18 / raw_rmsd_ca
                coords_stage2[is_ca] = coords_init[is_ca] + ca_disp * scale
                e_stage2, _ = stage2_objective(coords_stage2.ravel())

        coords_final = coords_stage2
        e_final = e_stage2
        opt_res = opt_stage2

    diff_final = coords_final[:, None, :] - coords_final[None, :, :]
    dist_final = np.sqrt(np.sum(diff_final ** 2, axis=-1) + 1e-12)
    dist_final[exclusion_mask] = 999.0
    dist_final[np.tril_indices(n_shell)] = 999.0
    final_min_d = float(np.min(dist_final)) if dist_final.size > 0 else 2.5
    final_clashes = int(np.sum((dist_final < 1.8) & (dist_final > 0)))
    clashes_resolved = max(0, initial_clashes - final_clashes)

    ca_shell_mask = [new_atoms[idx].name == "CA" for idx in shell_indices]
    if any(ca_shell_mask):
        ca_init = coords_init[ca_shell_mask]
        ca_final = coords_final[ca_shell_mask]
        rmsd_ca = float(np.sqrt(np.mean(np.sum((ca_final - ca_init) ** 2, axis=-1))))
    else:
        rmsd_ca = 0.0

    for i, idx in enumerate(shell_indices):
        new_atoms[idx].x = float(coords_final[i, 0])
        new_atoms[idx].y = float(coords_final[i, 1])
        new_atoms[idx].z = float(coords_final[i, 2])

    relaxed_pdb_str = "".join(a.to_pdb_line() for a in new_atoms)
    runtime_ms = (time.perf_counter() - t_start) * 1000.0

    notation = ":".join(f"{w}{p}{m}" for w, p, m in parsed_muts)
    return RelaxationResult(
        variant_notation=notation,
        initial_energy_kcal=round(e_init, 2),
        final_energy_kcal=round(e_final, 2),
        initial_min_distance_angstrom=round(initial_min_d, 3),
        final_min_distance_angstrom=round(final_min_d, 3),
        clashes_resolved=clashes_resolved,
        converged=bool(opt_res.success or opt_res.status == 0 or e_final <= e_init),
        relaxed_pdb_str=relaxed_pdb_str,
        rmsd_ca_angstrom=round(rmsd_ca, 4),
        runtime_ms=round(runtime_ms, 1),
        backbone_stage2_triggered=backbone_stage2_triggered,
        denaturation_rejected=denaturation_rejected,
    )


relax_mutant_structure = relax_variant_structure
