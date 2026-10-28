"""
workers/mock_workers.py — DOMAIN: Science & Math Tools.

Purpose:
    Deterministic, high-fidelity biological mock stand-ins for RAG, LLM, and ESM workers.
    Allows the entire 6-step active learning cycle, UI dashboard, and testing suites
    to execute locally with 100% fidelity, zero network calls, and $0.00 cloud cost.

Behavior & Fidelity (INNOVATION_ROADMAP Phase 4):
    1. Realistic GFP Fitness Landscape: Simulates known avGFP deep mutational scanning (DMS)
       phenotypes (e.g. S65T maturation boost, T203Y yellow-shift stacking, Gly67 disruption).
    2. Epistatic Synergies (G2): Implements non-linear synergy bonuses for beneficial combinations
       (e.g. S65T:T203Y), demonstrating that combinatorial active learning beats single-mutant scans.
    3. Drop-in Replacement: Shares the exact interface with live Modal ESM-2 and Meta Llama 3.1 workers.
"""

from __future__ import annotations

import hashlib
import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple

from controller.gatekeeper import SingleMutation, Variant
from controller.prompts import BiophysicalHypothesis, ResidueCluster, SearchStrategy
from workers.rag_worker import ProteinMetadata, fetch_protein_metadata

logger = logging.getLogger(__name__)

# Known single-mutation fitness effects in avGFP (normalized relative to WT = 0.0)
_BENEFICIAL_EFFECTS: Dict[str, float] = {
    # Chromophore environment
    "S65T": 0.85,   # EGFP mutation (accelerated maturation, higher quantum yield)
    "S65A": 0.40,   # Enhanced folding kinetics
    "S65G": 0.35,   # Maturation variant
    "V68L": 0.30,   # Improved packing
    "Q69M": 0.25,   # Citrine-stabilizing mutation
    # Pocket & contact residues
    "T203Y": 0.95,  # EYFP mutation (pi-stacking with chromophore)
    "T203F": 0.60,  # Hydrophobic aromatic contact
    "T203H": 0.45,  # Partial ring interaction
    "A206K": 0.55,  # Monomeric GFP (mGFP), eliminates dimerization artifact
    "E222G": 0.70,  # Deprotonation enhancer, suppresses 395nm absorption
    "F64L": 0.40,   # Thermal stability at 37°C
    "S147P": 0.20,  # Improved loop flexibility
    # Spectral shift alternatives
    "Y66H": 0.20,   # BFP (blue fluorescent protein)
    "Y66W": 0.15,   # CFP (cyan fluorescent protein)
    # Superfolder GFP (sfGFP) folding & stability mutations
    "V163A": 0.65,  # Superfolder core packing enhancement
    "T97A": 0.45,   # Beta-strand 5 folding kinetics
    "M153T": 0.50,  # Beta-strand 7 hydrophobic core packing
    "S30R": 0.35,   # N-terminal loop stabilization
    "N212K": 0.35,  # C-terminal barrel electrostatic stabilization
}

# Synergistic epistatic pairs (non-linear positive epistasis)
_EPISTATIC_SYNERGIES: Dict[frozenset, float] = {
    frozenset(["S65T", "T203Y"]): 0.45,  # Classic dual-enhanced fluorophore
    frozenset(["S65T", "A206K"]): 0.30,  # Bright monomeric variant
    frozenset(["S65T", "E222G"]): 0.35,  # Single-peak green emitter
    frozenset(["T203Y", "E222G"]): 0.25, # High-quantum yellow emitter
    frozenset(["F64L", "S65T"]): 0.30,   # 37°C enhanced maturation
    frozenset(["S65T", "V163A"]): 0.35,  # Core fluorophore + superfolder packing
    frozenset(["T203Y", "V163A"]): 0.30, # Pi-stacking + stability
    frozenset(["T97A", "M153T"]): 0.35,  # Superfolder dual-strand cooperativity
    frozenset(["V163A", "M153T"]): 0.30, # Beta-sheet core stabilization
}


def compute_mock_gfp_fitness(
    mutations: Sequence[SingleMutation],
    deterministic_noise: bool = True,
) -> float:
    """
    Computes a realistic, continuous fitness score for an avGFP mutant.
    
    Returns:
        float score where 0.0 = WT baseline, >0 is brighter/improved, <0 is impaired.
    """
    if not mutations:
        return 0.0  # Wild-type reference

    base_score = 0.0
    mut_strs = [m.to_string() for m in mutations]
    mut_set = set(mut_strs)

    # 1. Add single-mutation effects
    for m in mutations:
        notation = m.to_string()
        if notation in _BENEFICIAL_EFFECTS:
            base_score += _BENEFICIAL_EFFECTS[notation]
        elif m.position == 66:
            # Tyr66 forms the conjugated ring; mutating to non-aromatic destroys fluorescence
            base_score -= 1.8
        elif m.position == 67:
            # Gly67 provides the tight steric backbone hinge; mutating it causes misfolding
            base_score -= 2.2
        elif m.mut_aa == "P":
            # Proline introduction inside beta-sheets disrupts hydrogen bonding
            base_score -= 0.65
        else:
            # Neutral / mild evolutionary penalty for random substitutions
            base_score -= 0.08

    # 2. Add non-linear epistatic interactions
    for pair, synergy in _EPISTATIC_SYNERGIES.items():
        if pair.issubset(mut_set):
            base_score += synergy

    # 3. Deterministic pseudo-noise for realism
    if deterministic_noise:
        seed_str = ":".join(sorted(mut_strs))
        hash_val = int(hashlib.md5(seed_str.encode()).hexdigest()[:8], 16)
        noise = ((hash_val % 1000) / 1000.0 - 0.5) * 0.04  # [-0.02, +0.02]
        base_score += noise

    return round(base_score, 4)


class MockESMWorker:
    """Mock scoring worker implementing the exact interface of workers.esm_worker."""

    def __init__(self, wt_sequence: Optional[str] = None) -> None:
        if wt_sequence is None:
            metadata = fetch_protein_metadata("P42212", use_cache=True)
            self.wt_sequence = metadata.sequence
        else:
            self.wt_sequence = wt_sequence

    def score_variants(
        self,
        wt_sequence: str,
        variant_sequences: Sequence[str],
    ) -> List[float]:
        """Scores a batch of full protein sequences against wild-type."""
        scores: List[float] = []
        for v_seq in variant_sequences:
            if v_seq == wt_sequence:
                scores.append(0.0)
                continue
            # Extract single mutations by comparing sequence differences
            muts: List[SingleMutation] = []
            for i, (wt_c, mut_c) in enumerate(zip(wt_sequence, v_seq)):
                if wt_c != mut_c:
                    muts.append(SingleMutation(wt_aa=wt_c, position=i + 1, mut_aa=mut_c))
            scores.append(compute_mock_gfp_fitness(muts))
        return scores

    def score_variants_cascade(
        self,
        wt_sequence: str,
        variant_sequences: Sequence[str],
        top_k_ratio: float = 0.20,
    ) -> Dict[str, Any]:
        """
        Simulates the Multi-Fidelity Cloud Cascade:
          - Tier 1: Fast 150M screening on all N variants with realistic mild approximation noise.
          - Promotes top quantile = ceil(N * top_k_ratio).
          - Tier 2: Precision 650M verification on promoted variants.
        """
        import math
        import numpy as np

        full_scores = self.score_variants(wt_sequence, variant_sequences)
        n = len(variant_sequences)
        if n == 0:
            return {
                "tier1_scores": [],
                "promoted_indices": [],
                "promoted_sequences": [],
                "tier2_scores": [],
            }

        # Simulated Tier-1 scores: true score + small approximation residual (correlation ~0.88 with 650M)
        tier1_scores: List[float] = []
        for i, (v_seq, s) in enumerate(zip(variant_sequences, full_scores)):
            seed = int(hashlib.sha256(f"t1_{v_seq}".encode()).hexdigest()[:8], 16)
            noise = ((seed % 1000) / 1000.0 - 0.5) * 0.15
            tier1_scores.append(round(s + noise, 4))

        k = max(1, math.ceil(n * top_k_ratio))
        promoted_indices = np.argsort(tier1_scores)[::-1][:k].tolist()
        promoted_seqs = [variant_sequences[i] for i in promoted_indices]
        tier2_scores = [full_scores[i] for i in promoted_indices]

        return {
            "tier1_scores": tier1_scores,
            "promoted_indices": promoted_indices,
            "promoted_sequences": promoted_seqs,
            "tier2_scores": tier2_scores,
        }

    def embed_sequences(self, sequences: Sequence[str], dim: int = 1280) -> List[List[float]]:
        """Generates deterministic latent sequence embeddings for surrogate models."""
        import numpy as np
        embeddings: List[List[float]] = []
        for seq in sequences:
            seed = int(hashlib.sha256(seq.encode("utf-8")).hexdigest()[:8], 16)
            rng = np.random.RandomState(seed)
            vec = rng.randn(dim).astype(np.float32)
            vec /= float(np.linalg.norm(vec) + 1e-8)
            embeddings.append(vec.tolist())
        return embeddings


class MockLLMWorker:
    """Mock strategic planner that produces valid SearchStrategy objects."""

    def __init__(self) -> None:
        pass

    def generate_strategy(
        self,
        metadata: ProteinMetadata,
        round_number: int = 1,
        top_evidence: Optional[Sequence[Dict[str, Any]]] = None,
        current_hop: int = 0,
    ) -> SearchStrategy:
        """Returns a biologically valid, round-progressive SearchStrategy dynamically for any protein target."""
        prot_name = metadata.protein_name
        acc = metadata.accession
        cat_pos = metadata.catalytic_positions
        bind_pos = metadata.binding_positions

        # Target: TrpB / Tryptophan synthase family
        if "trpb" in prot_name.lower() or "tryptophan synthase" in prot_name.lower() or acc in ("P50909", "TRPB_TM9D8"):
            if current_hop >= 1:
                hypothesis = f"Pocket Hop {current_hop}: Sweeping opposing catalytic communication loop (positions 227, 228) coupled to the locked lead anchor."
                clusters = [
                    ResidueCluster(
                        cluster_id="comm_loop_2",
                        positions=[227, 228],
                        rationale="Opposing loop residues directly contacting loop 1 to stabilize catalytic closure.",
                    )
                ]
                mechanistic_hyps = [
                    BiophysicalHypothesis(
                        pocket_name="communication_loop_2",
                        target_residues=[227, 228],
                        mechanism="Conformational lid closure stabilization coupled to locked Loop 1 lead.",
                        structural_rationale="Loop 2 residues (227, 228) form a dynamic cap over the PLP pocket, interacting directly with Loop 1 at 7.8 Angstroms to accelerate substrate turnover.",
                        suggested_chemistries={"227": ["M", "L", "Q"], "228": ["G", "A", "S"]},
                        epistatic_coupling="Positive synergistic coupling: Loop 2 mutations stabilize the closed catalytic conformation around the locked Loop 1 lead.",
                    )
                ]
                return SearchStrategy(
                    target_protein=metadata.protein_name,
                    target_accession=metadata.accession,
                    hypothesis=hypothesis,
                    target_clusters=clusters,
                    forbidden_positions=sorted(list(set([83] + cat_pos))),
                    max_mutation_depth=4,
                    mechanistic_hypotheses=mechanistic_hyps,
                )
            elif round_number == 1:
                hypothesis = "Round 1: Initial sampling of substrate entrance loop 1 (positions 183, 184) to explore active-site accessibility."
                clusters = [
                    ResidueCluster(
                        cluster_id="comm_loop_1",
                        positions=[183, 184],
                        rationale="Active site entrance loop gating substrate access to the PLP cofactor pocket.",
                    )
                ]
                mechanistic_hyps = [
                    BiophysicalHypothesis(
                        pocket_name="substrate_entrance_loop_1",
                        target_residues=[183, 184],
                        mechanism="Steric volume expansion gating the PLP cofactor entrance channel.",
                        structural_rationale="Phe184 is 6.2 Angstroms from the catalytic Lys87 Schiff base. Modulating bulky Phe to branched hydrophobics (Ile/Val) eliminates steric clash with substituted indoles.",
                        suggested_chemistries={"183": ["A", "V", "C"], "184": ["I", "V", "L"]},
                        epistatic_coupling="Gating step: establishes the substrate channel baseline before coupling to Loop 2.",
                    )
                ]
            else:
                hypothesis = f"Round {round_number}: Combinatorial coupling of substrate entrance loops (183, 184, 227, 228) to modulate catalytic turnover."
                clusters = [
                    ResidueCluster(
                        cluster_id="comm_loop_1",
                        positions=[183, 184],
                        rationale="Active site entrance loop gating substrate access.",
                    ),
                    ResidueCluster(
                        cluster_id="comm_loop_2",
                        positions=[227, 228],
                        rationale="Opposing loop residues directly contacting loop 1 to stabilize catalytic closure.",
                    ),
                ]
                mechanistic_hyps = [
                    BiophysicalHypothesis(
                        pocket_name="inter_loop_lid_complex",
                        target_residues=[183, 184, 227, 228],
                        mechanism="Dual-loop cooperative epistatic lid modulation gating catalytic turnover.",
                        structural_rationale="Loops 1 and 2 face each other across the PLP cofactor entrance (<8.5 Angstroms). Synergistic mutations balance open/closed transition rates.",
                        suggested_chemistries={
                            "183": ["A", "V", "C"],
                            "184": ["I", "V", "L"],
                            "227": ["M", "L", "Q"],
                            "228": ["G", "A", "S"],
                        },
                        epistatic_coupling="Strong positive non-linear epistasis between 184 and 227/228 yields superior turnover.",
                    )
                ]
            return SearchStrategy(
                target_protein=metadata.protein_name,
                target_accession=metadata.accession,
                hypothesis=hypothesis,
                target_clusters=clusters,
                forbidden_positions=sorted(list(set([83] + cat_pos))),
                max_mutation_depth=4 if round_number > 1 else 2,
                mechanistic_hypotheses=mechanistic_hyps,
            )

        # Target: General protein target with UniProt binding pocket annotations
        if bind_pos and len(bind_pos) >= 2 and acc != "P42212":
            c1 = sorted(bind_pos[:min(3, len(bind_pos))])
            c2 = sorted(bind_pos[min(3, len(bind_pos)):min(6, len(bind_pos))])
            clusters = [
                ResidueCluster(
                    cluster_id="active_pocket_shell",
                    positions=c1,
                    rationale="Substrate/cofactor binding pocket residues directly modulating active-site geometry.",
                )
            ]
            if c2:
                clusters.append(
                    ResidueCluster(
                        cluster_id="pocket_interaction_network",
                        positions=c2,
                        rationale="Secondary shell interacting residues modulating conformational flexibility.",
                    )
                )
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="binding_pocket_shell",
                    target_residues=c1 + (c2 if c2 else []),
                    mechanism="Binding pocket cavity and electrostatic tuning.",
                    structural_rationale=f"Residues {c1} directly coordinate substrate within 8 Angstroms of catalytic machinery {cat_pos}.",
                    suggested_chemistries={str(p): ["A", "L", "S", "T"] for p in c1},
                    epistatic_coupling="First-shell cavity packing directly modulates ligand binding affinity.",
                )
            ]
            return SearchStrategy(
                target_protein=metadata.protein_name,
                target_accession=metadata.accession,
                hypothesis=f"Round {round_number}: Targeting binding pocket contact network while strictly preserving catalytic machinery {cat_pos}.",
                target_clusters=clusters,
                forbidden_positions=sorted(list(set(cat_pos + [1, metadata.length]))),
                max_mutation_depth=2,
                mechanistic_hypotheses=mechanistic_hyps,
            )

        # Default: avGFP round-progressive strategy
        if round_number == 1:
            hypothesis = (
                "Round 1: Screening catalytic chromophore core (positions 64, 65, 68) "
                "to discover baseline fluorophore maturation and quantum yield enhancers."
            )
            clusters = [
                ResidueCluster(
                    cluster_id="chromophore_core",
                    positions=[64, 65, 68],
                    rationale="Direct chromophore environment; modulates excitation peaks and maturation rate.",
                ),
            ]
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="chromophore_core",
                    target_residues=[64, 65, 68],
                    mechanism="Fluorophore maturation kinetics and excited-state rigidification.",
                    structural_rationale="Phe64 and Ser65 flank the internal chromophore triad (65-67). S65T stabilizes the deprotonated phenolate anion.",
                    suggested_chemistries={"64": ["L", "M"], "65": ["T", "A", "G"], "68": ["V", "L"]},
                    epistatic_coupling="Ser65 controls chromophore ionization state and red-shifts excitation.",
                )
            ]
        elif round_number == 2:
            hypothesis = (
                "Round 2: Expanding to first-shell pi-stacking contacts (positions 203, 206, 222) "
                "coupled with chromophore core to unlock dual-enhanced yellow/green fluorophores."
            )
            clusters = [
                ResidueCluster(
                    cluster_id="chromophore_lead",
                    positions=[65, 68],
                    rationale="Retaining primary fluorophore maturation locus.",
                ),
                ResidueCluster(
                    cluster_id="pi_stacking_pocket",
                    positions=[203, 206, 222],
                    rationale="Spatial interaction partners in contact with internal chromophore (T203Y, A206K).",
                ),
            ]
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="pi_stacking_shell",
                    target_residues=[65, 203, 222],
                    mechanism="Aromatic pi-pi stacking over chromophore phenol ring.",
                    structural_rationale="Thr203 is positioned 3.4 Angstroms directly above the Tyr66 phenol ring. Substitution to Tyr (T203Y) establishes continuous pi-orbital overlap.",
                    suggested_chemistries={"65": ["T"], "203": ["Y", "F", "H"], "222": ["E", "U", "V"]},
                    epistatic_coupling="Synergy between 65 ionization and 203 pi-stacking produces enhanced quantum yield.",
                )
            ]
        elif round_number == 3:
            hypothesis = (
                "Round 3: Targeting beta-can barrel structural rigidity and loop dynamics (positions 147, 163, 167) "
                "to eliminate 37°C thermal quenching and enhance conformational stability."
            )
            clusters = [
                ResidueCluster(
                    cluster_id="active_pocket_lead",
                    positions=[65, 203],
                    rationale="Coupling discovered chromophore hits with barrel scaffold.",
                ),
                ResidueCluster(
                    cluster_id="beta_barrel_rigidity",
                    positions=[147, 163, 167],
                    rationale="Internal barrel packing and loop dynamics (V163A, S147P).",
                ),
            ]
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="beta_barrel_scaffold",
                    target_residues=[147, 163, 167],
                    mechanism="Beta-barrel tight packing and protection against solvent quenching.",
                    structural_rationale="Residues 163 and 167 pack against the internal central helix, shielding the chromophore from non-radiative vibrational relaxation.",
                    suggested_chemistries={"147": ["P", "S"], "163": ["A", "V"], "167": ["T", "I"]},
                    epistatic_coupling="Scaffold rigidity prevents thermal conformational breathing.",
                )
            ]
        elif round_number == 4:
            hypothesis = (
                "Round 4: Introducing canonical Superfolder GFP folding chaperones (positions 30, 97, 153, 212) "
                "to accelerate beta-sheet assembly and prevent off-pathway inclusion bodies."
            )
            clusters = [
                ResidueCluster(
                    cluster_id="scaffold_core_lead",
                    positions=[65, 163, 203],
                    rationale="Consolidating high-fitness core alleles.",
                ),
                ResidueCluster(
                    cluster_id="superfolder_chaperones",
                    positions=[30, 97, 153, 212],
                    rationale="Canonical superfolder GFP folding enhancers (T97A, M153T, S30R).",
                ),
            ]
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="superfolder_folding_nuclei",
                    target_residues=[30, 97, 153, 212],
                    mechanism="Accelerated beta-barrel folding kinetics and aggregation suppression.",
                    structural_rationale="M153T and T97A stabilize turns between strands 7-8 and 4-5, driving productive folding before chromophore cyclization.",
                    suggested_chemistries={"30": ["R"], "97": ["A"], "153": ["T"], "212": ["N"]},
                    epistatic_coupling="Fast folding kinetics synergizes with active-site quantum yield improvements.",
                )
            ]
        else:
            hypothesis = (
                f"Round {round_number}: Global multi-module combinatorial consolidation across all "
                "discovered beneficial hotspots (65, 97, 153, 163, 203) to assemble peak super-mutants."
            )
            clusters = [
                ResidueCluster(
                    cluster_id="global_hotspots",
                    positions=[65, 97, 153, 163, 203],
                    rationale="Global multi-module epistatic consolidation across all validated loci.",
                ),
            ]
            mechanistic_hyps = [
                BiophysicalHypothesis(
                    pocket_name="global_epistatic_network",
                    target_residues=[65, 97, 153, 163, 203],
                    mechanism="Multi-locus cooperative epistasis uniting chromophore, pi-stacking, and folding.",
                    structural_rationale="Consolidates validated alleles across all 3 structural domains into optimal combinatorial fitness maxima.",
                    suggested_chemistries={"65": ["T"], "97": ["A"], "153": ["T"], "163": ["A"], "203": ["Y"]},
                    epistatic_coupling="Peak multi-point epistatic synergy across the entire protein fold.",
                )
            ]

        return SearchStrategy(
            target_protein=metadata.protein_name,
            target_accession=metadata.accession,
            hypothesis=hypothesis,
            target_clusters=clusters,
            forbidden_positions=[1, 2, 67, 238],
            max_mutation_depth=2,
            mechanistic_hypotheses=mechanistic_hyps,
        )
