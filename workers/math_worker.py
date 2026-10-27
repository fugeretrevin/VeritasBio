"""
workers/math_worker.py — DOMAIN: Science & Math Tools (Pillar 3: Quant Engine).

Purpose:
    Continuous physicochemical featurization + Gaussian Process Bayesian Optimization.
    Consumes the Strategic Planner's `SearchStrategy`, featurizes candidate sequence
    variants into continuous biophysical delta coordinates, fits a GP to ledger Evidence,
    and calculates vectorized Expected Improvement (EI) to select the next 50-variant micro-batch.

Invariants (ARCHITECTURE.md §3.3; INNOVATION_ROADMAP Phase 3):
    1. Rigorous Math Only: No LLM involvement in coordinate featurization, GP fitting,
       or EI calculation.
    2. Overcoming Discrete Bottleneck (G3): Uses normalized 5-factor amino acid properties
       (hydrophobicity, volume/mass, charge, polarity, pI) to create dense, continuous
       vectors that scale efficiently under Gaussian Process kernels.
    3. Combinatorial Exploration (G2): Explores multi-point combinatorial variants
       across proposed residue clusters up to `max_mutation_depth`.
"""

from __future__ import annotations

import itertools
import logging
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
from pydantic import BaseModel, Field
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from controller.gatekeeper import CANONICAL_AMINO_ACIDS, SingleMutation, Variant
from controller.prompts import SearchStrategy

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Biophysical Property Scales (AAindex Reference Values)
# 5 factors: [Hydrophobicity (Kyte-Doolittle), Molecular Mass (Da),
#             Net Charge (pH 7.4), Polarity (Grantham), Isoelectric Point (pI)]
# -----------------------------------------------------------------------------
_RAW_AA_PROPERTIES: Dict[str, np.ndarray] = {
    "A": np.array([1.8, 89.09, 0.0, 8.1, 6.00], dtype=np.float64),
    "C": np.array([2.5, 121.16, 0.0, 5.5, 5.07], dtype=np.float64),
    "D": np.array([-3.5, 133.10, -1.0, 13.0, 2.77], dtype=np.float64),
    "E": np.array([-3.5, 147.13, -1.0, 12.3, 3.22], dtype=np.float64),
    "F": np.array([2.8, 165.19, 0.0, 5.2, 5.48], dtype=np.float64),
    "G": np.array([-0.4, 75.07, 0.0, 9.0, 5.97], dtype=np.float64),
    "H": np.array([-3.2, 155.16, 0.1, 10.4, 7.59], dtype=np.float64),
    "I": np.array([4.5, 131.17, 0.0, 5.2, 6.02], dtype=np.float64),
    "K": np.array([-3.9, 146.19, 1.0, 11.3, 9.74], dtype=np.float64),
    "L": np.array([3.8, 131.17, 0.0, 4.9, 5.98], dtype=np.float64),
    "M": np.array([1.9, 149.21, 0.0, 5.7, 5.74], dtype=np.float64),
    "N": np.array([-3.5, 132.12, 0.0, 11.6, 5.41], dtype=np.float64),
    "P": np.array([-1.6, 115.13, 0.0, 8.0, 6.30], dtype=np.float64),
    "Q": np.array([-3.5, 146.15, 0.0, 10.5, 5.65], dtype=np.float64),
    "R": np.array([-4.5, 174.20, 1.0, 10.5, 10.76], dtype=np.float64),
    "S": np.array([-0.8, 105.09, 0.0, 9.2, 5.68], dtype=np.float64),
    "T": np.array([-0.7, 119.12, 0.0, 8.6, 5.60], dtype=np.float64),
    "V": np.array([4.2, 117.15, 0.0, 5.9, 5.96], dtype=np.float64),
    "W": np.array([-0.9, 204.23, 0.0, 5.4, 5.89], dtype=np.float64),
    "Y": np.array([-1.3, 181.19, 0.0, 6.2, 5.66], dtype=np.float64),
}

# Standardize properties to zero mean and unit variance across the 20 amino acids
_PROP_MATRIX = np.array(list(_RAW_AA_PROPERTIES.values()), dtype=np.float64)
_PROP_MEAN = np.mean(_PROP_MATRIX, axis=0)
_PROP_STD = np.std(_PROP_MATRIX, axis=0)
_PROP_STD[_PROP_STD == 0] = 1.0

AA_PROPERTIES_NORMALIZED: Dict[str, np.ndarray] = {
    aa: (vec - _PROP_MEAN) / _PROP_STD for aa, vec in _RAW_AA_PROPERTIES.items()
}
AA_FEATURE_DIM: int = 5


# -----------------------------------------------------------------------------
# Featurization
# -----------------------------------------------------------------------------
def featurize_single_delta(wt_aa: str, mut_aa: str) -> np.ndarray:
    """Computes normalized delta vector: mut_properties - wt_properties."""
    wt_vec = AA_PROPERTIES_NORMALIZED.get(wt_aa.upper())
    mut_vec = AA_PROPERTIES_NORMALIZED.get(mut_aa.upper())
    if wt_vec is None or mut_vec is None:
        raise ValueError(f"Invalid amino acid comparison: '{wt_aa}' -> '{mut_aa}'")
    return mut_vec - wt_vec


def featurize_variant_on_positions(
    variant: Variant,
    active_positions: Sequence[int],
    wt_sequence: str,
) -> np.ndarray:
    """
    Constructs a fixed-length continuous vector of length len(active_positions) * 5.
    For each position p in active_positions:
      - if mutated: inserts 5D property delta vector.
      - if unmutated: inserts zeros(5).
    """
    variant.validate_against_wt(wt_sequence)
    mut_dict: Dict[int, SingleMutation] = {m.position: m for m in variant.mutations}
    chunks: List[np.ndarray] = []

    for pos in active_positions:
        if pos in mut_dict:
            m = mut_dict[pos]
            chunks.append(featurize_single_delta(m.wt_aa, m.mut_aa))
        else:
            chunks.append(np.zeros(AA_FEATURE_DIM, dtype=np.float64))

    return np.concatenate(chunks)


# -----------------------------------------------------------------------------
# Bayesian Optimization & Expected Improvement
# -----------------------------------------------------------------------------
class SelectedBatch(BaseModel):
    """Output batch selected by the Quant Engine for Gatekeeper authorization."""

    variants: List[Variant] = Field(..., description="Selected candidate variants")
    expected_improvements: List[float] = Field(..., description="Computed EI value per variant")
    predicted_means: List[float] = Field(..., description="GP predicted mean fitness")
    predicted_uncertainties: List[float] = Field(..., description="GP predicted standard deviation")


def compute_expected_improvement(
    mean: np.ndarray,
    std: np.ndarray,
    f_best: float,
    xi: float = 0.01,
) -> np.ndarray:
    """
    Vectorized Expected Improvement (EI) acquisition function:
      EI(x) = (mu(x) - f_best - xi) * Phi(Z) + sigma(x) * phi(Z)
    """
    std = np.maximum(std, 1e-9)
    improvement = mean - f_best - xi
    z = improvement / std
    ei = improvement * norm.cdf(z) + std * norm.pdf(z)
    return np.maximum(ei, 0.0)


def build_gp_model(random_state: int = 42) -> GaussianProcessRegressor:
    """
    Constructs a well-conditioned Gaussian Process with Matern 5/2 kernel
    and robust alpha diagonal regularization for stable small-sample convergence.
    """
    kernel = ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
        length_scale=1.0, length_scale_bounds=(1e-1, 1e2), nu=2.5
    )
    return GaussianProcessRegressor(
        kernel=kernel,
        alpha=1e-2,
        n_restarts_optimizer=1,
        normalize_y=True,
        random_state=random_state,
    )


def extract_mutations_from_sequence(mutated_seq: str, wt_seq: str) -> List[SingleMutation]:
    """Infers single mutations by comparing mutated sequence to WT sequence."""
    if len(mutated_seq) != len(wt_seq):
        raise ValueError("Sequence length mismatch between mutated and WT")
    muts: List[SingleMutation] = []
    for i, (wt_char, mut_char) in enumerate(zip(wt_seq, mutated_seq)):
        if wt_char != mut_char:
            muts.append(SingleMutation(wt_aa=wt_char, position=i + 1, mut_aa=mut_char))
    return muts


def generate_candidate_variants(
    strategy: SearchStrategy,
    wt_sequence: str,
    max_candidates: int = 2500,
    diverse_amino_acids: Optional[Sequence[str]] = None,
    lead_variants: Optional[List[Variant]] = None,
) -> List[Variant]:
    """
    Generates combinatorial candidate variants across and between target clusters.
    Prioritizes high-impact physicochemical substitutions (T, Y, L, K, A, R)
    to discover synergistic epistatic multi-point improvements.
    Also executes iterative lead maturation on top of prior high-performing leads.
    """
    # Prioritize functionally impactful amino acid substitutions across size, polarity, and charge
    target_aas = list(diverse_amino_acids) if diverse_amino_acids else [
        "T", "Y", "L", "K", "A", "R", "D", "E", "G", "V", "F", "S", "H", "N", "Q"
    ]

    all_cluster_positions: Set[int] = set()
    for cluster in strategy.target_clusters:
        all_cluster_positions.update(cluster.positions)

    # Exclude any forbidden positions
    allowed_positions = sorted(list(all_cluster_positions - set(strategy.forbidden_positions)))

    if not allowed_positions:
        raise ValueError("No allowed positions available after filtering forbidden positions")

    candidates: List[Variant] = []
    seen_notations: Set[str] = set()

    # 1. Single mutants across all allowed cluster positions
    for pos in allowed_positions:
        wt_aa = wt_sequence[pos - 1]
        for mut_aa in target_aas:
            if mut_aa != wt_aa:
                notation = f"{wt_aa}{pos}{mut_aa}"
                if notation not in seen_notations:
                    seen_notations.add(notation)
                    candidates.append(Variant(notation=notation))

    # 2. Iterative Lead Maturation: Stacking new cluster mutations on top of discovered high-fitness leads
    if lead_variants:
        for lead in lead_variants[:5]:
            lead_positions = {m.position for m in lead.mutations}
            # Find candidate positions in current strategy not yet mutated in lead
            available_pos = [p for p in allowed_positions if p not in lead_positions]
            for pos in available_pos[:6]:
                wt_aa = wt_sequence[pos - 1]
                for mut_aa in target_aas[:8]:
                    if mut_aa != wt_aa:
                        new_muts = sorted(
                            lead.mutations + [SingleMutation(wt_aa=wt_aa, position=pos, mut_aa=mut_aa)],
                            key=lambda m: m.position,
                        )
                        notation = ":".join(m.to_string() for m in new_muts)
                        if notation not in seen_notations:
                            seen_notations.add(notation)
                            candidates.append(Variant(notation=notation, mutations=new_muts))
                            if len(candidates) >= max_candidates:
                                return candidates

    # 3. Combinatorial double mutants within clusters (intra-cluster)
    for cluster in strategy.target_clusters:
        cluster_allowed = [p for p in cluster.positions if p in allowed_positions]
        if len(cluster_allowed) >= 2 and strategy.max_mutation_depth >= 2:
            for pos1, pos2 in itertools.combinations(cluster_allowed[:4], 2):
                wt1 = wt_sequence[pos1 - 1]
                wt2 = wt_sequence[pos2 - 1]
                for mut1 in target_aas[:6]:
                    if mut1 == wt1:
                        continue
                    for mut2 in target_aas[:6]:
                        if mut2 == wt2:
                            continue
                        notation = f"{wt1}{pos1}{mut1}:{wt2}{pos2}{mut2}"
                        if notation not in seen_notations:
                            seen_notations.add(notation)
                            candidates.append(Variant(notation=notation))
                            if len(candidates) >= max_candidates:
                                return candidates

    # 4. Inter-cluster cross-talk double mutants (e.g. chromophore triad + pocket contact)
    if len(strategy.target_clusters) >= 2 and strategy.max_mutation_depth >= 2:
        for c1, c2 in itertools.combinations(strategy.target_clusters, 2):
            c1_pos = [p for p in c1.positions if p in allowed_positions][:3]
            c2_pos = [p for p in c2.positions if p in allowed_positions][:3]
            for pos1 in c1_pos:
                for pos2 in c2_pos:
                    if pos1 == pos2:
                        continue
                    wt1 = wt_sequence[pos1 - 1]
                    wt2 = wt_sequence[pos2 - 1]
                    for mut1 in target_aas[:5]:
                        if mut1 == wt1:
                            continue
                        for mut2 in target_aas[:5]:
                            if mut2 == wt2:
                                continue
                            notation = f"{wt1}{pos1}{mut1}:{wt2}{pos2}{mut2}"
                            if notation not in seen_notations:
                                seen_notations.add(notation)
                                candidates.append(Variant(notation=notation))
                                if len(candidates) >= max_candidates:
                                    return candidates

    return candidates


# -----------------------------------------------------------------------------
# Main Quant Engine Selection Interface
# -----------------------------------------------------------------------------
def select_next_micro_batch(
    strategy: SearchStrategy,
    wt_sequence: str,
    evidence: Optional[Sequence[Dict[str, Any]]] = None,
    batch_size: int = 50,
    exploration_xi: float = 0.01,
) -> SelectedBatch:
    """
    Executes the Gaussian Process active learning selection:
      1. Aggregates all active exploration positions.
      2. Featurizes existing evidence (if any) and fits the GP model.
      3. Generates candidate combinatorial variants.
      4. Predicts mean and uncertainty, computes Expected Improvement.
      5. Returns top-N variants with their EI values.
    """
    # Collect all positions from the strategy clusters AND prior evidence
    active_positions_set: Set[int] = set(p for c in strategy.target_clusters for p in c.positions)
    lead_variants: List[Variant] = []

    if evidence:
        sorted_ev = sorted(evidence, key=lambda e: float(e.get("computed_score", 0.0)), reverse=True)
        for ev in sorted_ev:
            score_val = float(ev.get("computed_score", 0.0))
            var_seq = ev.get("variant_sequence")
            if var_seq and len(var_seq) == len(wt_sequence):
                mut_indices = [i for i, (wc, vc) in enumerate(zip(wt_sequence, var_seq), 1) if wc != vc]
                active_positions_set.update(mut_indices)
                if score_val > 0.0 and len(lead_variants) < 5:
                    try:
                        lead_variants.append(Variant.from_sequence(var_seq, wt_sequence))
                    except Exception:
                        pass

    active_positions: List[int] = sorted(
        [p for p in active_positions_set if p not in strategy.forbidden_positions]
    )
    if not active_positions:
        raise ValueError("Active exploration positions cannot be empty")

    # Step 1: Featurize training evidence if present
    X_train_list: List[np.ndarray] = []
    y_train_list: List[float] = []

    # Include WT sequence baseline as score 0.0 (reference)
    wt_features = np.zeros(len(active_positions) * AA_FEATURE_DIM, dtype=np.float64)
    X_train_list.append(wt_features)
    y_train_list.append(0.0)

    f_best = 0.0

    if evidence:
        for ev in evidence:
            var_seq = ev.get("variant_sequence")
            score = ev.get("computed_score")
            if var_seq and score is not None:
                try:
                    muts = extract_mutations_from_sequence(var_seq, wt_sequence)
                    var_obj = Variant(
                        notation=":".join(m.to_string() for m in muts) if muts else "WT",
                        mutations=muts,
                    )
                    feat = featurize_variant_on_positions(var_obj, active_positions, wt_sequence)
                    X_train_list.append(feat)
                    y_train_list.append(float(score))
                    if float(score) > f_best:
                        f_best = float(score)
                except Exception as exc:
                    logger.debug("Skipping unparseable evidence sequence: %s", exc)

    X_train = np.array(X_train_list, dtype=np.float64)
    y_train = np.array(y_train_list, dtype=np.float64)

    # Step 2: Fit Gaussian Process model
    gp = build_gp_model()
    try:
        gp.fit(X_train, y_train)
    except Exception as exc:
        logger.warning("GP fit failed (%s). Retrying with high regularization...", exc)
        gp = GaussianProcessRegressor(alpha=0.1, random_state=42)
        gp.fit(X_train, y_train)

    # Step 3: Generate candidate variants (both exploration and iterative lead maturation)
    candidates = generate_candidate_variants(strategy, wt_sequence, lead_variants=lead_variants)

    # Filter out any candidates already evaluated in evidence
    evaluated_notations = set()
    if evidence:
        for ev in evidence:
            seq = ev.get("variant_sequence")
            if seq:
                evaluated_notations.add(seq)

    novel_candidates = [
        c for c in candidates if c.apply_to(wt_sequence) not in evaluated_notations
    ]
    if not novel_candidates:
        novel_candidates = candidates  # Fallback if all covered

    # Step 4: Featurize candidate pool and predict GP mean + uncertainty
    X_cand = np.array(
        [featurize_variant_on_positions(c, active_positions, wt_sequence) for c in novel_candidates],
        dtype=np.float64,
    )

    means, stds = gp.predict(X_cand, return_std=True)

    # Step 5: Compute Expected Improvement
    ei_scores = compute_expected_improvement(means, stds, f_best=f_best, xi=exploration_xi)

    # Rank all candidates by Expected Improvement (descending)
    sorted_candidate_indices = np.argsort(ei_scores)[::-1]

    # Enforce Batch Acquisition Diversity:
    # Caps the maximum representation of any single residue in a single batch to avoid cluster collapse
    max_per_position = max(2, batch_size // 4)
    position_counts: Dict[int, int] = {}
    selected_indices: List[int] = []

    for idx in sorted_candidate_indices:
        cand = novel_candidates[idx]
        cand_positions = [m.position for m in cand.mutations]
        # Check if adding this variant satisfies the positional diversity cap
        if all(position_counts.get(p, 0) < max_per_position for p in cand_positions):
            selected_indices.append(idx)
            for p in cand_positions:
                position_counts[p] = position_counts.get(p, 0) + 1
            if len(selected_indices) >= batch_size:
                break

    # If diversity cap leaves slots unfilled, backfill with remaining highest EI candidates
    if len(selected_indices) < batch_size:
        already_picked = set(selected_indices)
        for idx in sorted_candidate_indices:
            if idx not in already_picked:
                selected_indices.append(idx)
                if len(selected_indices) >= batch_size:
                    break

    selected_variants = [novel_candidates[i] for i in selected_indices]
    selected_eis = [round(float(ei_scores[i]), 6) for i in selected_indices]
    selected_means = [round(float(means[i]), 4) for i in selected_indices]
    selected_stds = [round(float(stds[i]), 4) for i in selected_indices]

    logger.info(
        "Quant Engine selected %d variants (Max EI: %.4f, Mean predicted: %.4f, f_best: %.4f, Diversity Spread: %d positions)",
        len(selected_variants),
        selected_eis[0] if selected_eis else 0.0,
        float(np.mean(selected_means)) if selected_means else 0.0,
        f_best,
        len(position_counts),
    )

    return SelectedBatch(
        variants=selected_variants,
        expected_improvements=selected_eis,
        predicted_means=selected_means,
        predicted_uncertainties=selected_stds,
    )
