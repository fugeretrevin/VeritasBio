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
from controller.prompts import ResidueCluster, SearchStrategy
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
}

# Synergistic epistatic pairs (non-linear positive epistasis)
_EPISTATIC_SYNERGIES: Dict[frozenset, float] = {
    frozenset(["S65T", "T203Y"]): 0.45,  # Classic dual-enhanced fluorophore
    frozenset(["S65T", "A206K"]): 0.30,  # Bright monomeric variant
    frozenset(["S65T", "E222G"]): 0.35,  # Single-peak green emitter
    frozenset(["T203Y", "E222G"]): 0.25, # High-quantum yellow emitter
    frozenset(["F64L", "S65T"]): 0.30,   # 37°C enhanced maturation
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


class MockLLMWorker:
    """Mock strategic planner that produces valid SearchStrategy objects."""

    def __init__(self) -> None:
        pass

    def generate_strategy(
        self,
        metadata: ProteinMetadata,
        round_number: int = 1,
        top_evidence: Optional[Sequence[Dict[str, Any]]] = None,
    ) -> SearchStrategy:
        """Returns a biologically valid SearchStrategy for avGFP."""
        return SearchStrategy(
            target_protein=metadata.protein_name,
            target_accession=metadata.accession,
            hypothesis=(
                "Targeting the central chromophore triad (positions 65-67) and interacting "
                "steric pocket contacts (positions 203, 206, 222) to discover synergistic "
                "spectral and stability improvements."
            ),
            target_clusters=[
                ResidueCluster(
                    cluster_id="chromophore_triad",
                    positions=[64, 65, 66, 67, 68],
                    rationale="Direct chromophore environment; modulates excitation peaks and maturation rate.",
                ),
                ResidueCluster(
                    cluster_id="barrel_pocket",
                    positions=[203, 206, 222],
                    rationale="Spatial interaction partners in contact with the internal chromophore.",
                ),
            ],
            forbidden_positions=[1, 2, 238],
            max_mutation_depth=2,
        )
