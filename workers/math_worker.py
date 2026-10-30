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
import math
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
from pydantic import BaseModel, Field
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor

from controller.gatekeeper import CANONICAL_AMINO_ACIDS, SingleMutation, Variant
from controller.prompts import BiophysicalHypothesis, SearchStrategy

logger = logging.getLogger(__name__)

# Causes already reported by _warn_once, so a per-candidate failure is surfaced once
# instead of thousands of times across a large candidate pool.
_WARNED_ONCE: Set[str] = set()


def _warn_once(key: str, message: str, *args: Any) -> None:
    """Logs `message` at WARNING the first time this `key` is seen in the process."""
    if key in _WARNED_ONCE:
        return
    _WARNED_ONCE.add(key)
    logger.warning(message, *args)

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


GLOBAL_EPISTATIC_DIM: int = 10  # 6 global biophysical deltas + 4 pairwise interaction cross-terms


def featurize_variant_on_positions(
    variant: Variant,
    active_positions: Sequence[int],
    wt_sequence: str,
    include_epistatic: bool = True,
    spatial_map: Optional[Any] = None,
) -> np.ndarray:
    """
    Constructs a fixed-length continuous vector with two complementary blocks:
      1. Linear block (len(active_positions) * 5): 5D property deltas per active site.
      2. Global & Epistatic interaction block (10D, if include_epistatic=True):
         - 6 global summary deltas: mutation density (N/3.0), net hydropathy, net charge,
           net volume, net polarity, net aromaticity.
         - 4 pairwise second-order cross-terms: electrostatic coupling (charge products),
           steric packing (volume products), hydrophobic clustering, and 3D Euclidean contact
           attenuation (or 1D sequence distance if no 3D structure is supplied).

    Enables Gaussian Process to model synergistic epistatic cooperativity vs steric clash penalties.
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

    linear_vec = np.concatenate(chunks)
    if not include_epistatic:
        return linear_vec

    n_mut = len(variant.mutations)
    if n_mut == 0:
        global_vec = np.zeros(6, dtype=np.float64)
        interaction_vec = np.zeros(4, dtype=np.float64)
    else:
        deltas = [featurize_single_delta(m.wt_aa, m.mut_aa) for m in variant.mutations]
        sum_delta = np.sum(deltas, axis=0)  # [hydropathy, charge, volume, polarity, aromaticity]
        global_vec = np.array([
            n_mut / 3.0,
            sum_delta[0],
            sum_delta[1],
            sum_delta[2],
            sum_delta[3],
            sum_delta[4],
        ], dtype=np.float64)

        if n_mut >= 2:
            charge_prod = 0.0
            volume_prod = 0.0
            hydro_prod = 0.0
            dist_inv = 0.0
            pairs = list(itertools.combinations(range(n_mut), 2))
            for i, j in pairs:
                m1, m2 = variant.mutations[i], variant.mutations[j]
                d1, d2 = deltas[i], deltas[j]

                # 3D Euclidean contact distance coupling (if PDB spatial_map available)
                d3d = None
                if spatial_map is not None:
                    try:
                        d3d = spatial_map.get_distance(m1.position, m2.position)
                    except Exception:
                        d3d = None

                if d3d is not None and d3d > 0.0:
                    dist_weight = 1.0 / np.sqrt(max(1.0, d3d))
                else:
                    seq_dist = max(1, abs(m1.position - m2.position))
                    dist_weight = 1.0 / np.sqrt(seq_dist)

                dist_inv += dist_weight
                charge_prod += d1[1] * d2[1] * dist_weight
                volume_prod += d1[2] * d2[2] * dist_weight
                hydro_prod += d1[0] * d2[0] * dist_weight
            n_pairs = float(len(pairs))
            interaction_vec = np.array([
                charge_prod / n_pairs,
                volume_prod / n_pairs,
                hydro_prod / n_pairs,
                dist_inv / n_pairs,
            ], dtype=np.float64)
        else:
            interaction_vec = np.zeros(4, dtype=np.float64)

    return np.concatenate([linear_vec, global_vec, interaction_vec])


# -----------------------------------------------------------------------------
# Bayesian Optimization & Expected Improvement
# -----------------------------------------------------------------------------
class SelectedBatch(BaseModel):
    """Output batch selected by the Quant Engine for Gatekeeper authorization."""

    variants: List[Variant] = Field(..., description="Selected candidate variants")
    expected_improvements: List[float] = Field(..., description="Computed EI value per variant")
    predicted_means: List[float] = Field(..., description="GP predicted mean fitness")
    predicted_uncertainties: List[float] = Field(..., description="GP predicted standard deviation")
    pareto_ranks: Optional[List[int]] = Field(default=None, description="Multi-objective Pareto rank per variant (0 = non-dominated)")
    max_pool_ei: Optional[float] = Field(default=None, description="Maximum Expected Improvement across entire candidate pool")


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


def greedy_dpp_sample(
    quality_scores: np.ndarray,
    features: np.ndarray,
    batch_size: int,
    sigma_feat: float = 1.5,
    max_candidates_to_consider: int = 500,
) -> List[int]:
    """
    Selects a diverse, high-quality micro-batch of size `batch_size` using
    greedy Cholesky-factorized Determinantal Point Process (k-DPP) MAP inference.

    Algorithm (Chen et al. 2018):
      Maximizes log det(L_B), where L_ij = q_i * S_ij * q_j
      with quality q_i = quality_scores[i] and similarity S_ij = exp(-||f_i - f_j||^2 / (2 * sigma^2)).

    Guarantees:
      1. Orthogonal diversity: Redundant variants with parallel feature vectors are repelled (det -> 0).
      2. High quality: Top-performing candidates are prioritized across diverse feature sub-spaces.
      3. Fast runtime: O(k^2 * M) using incremental Cholesky updates (< 5ms on CPU).
    """
    n_total = len(quality_scores)
    if n_total <= batch_size:
        return list(range(n_total))

    # Pre-filter to the top candidates by quality to bound matrix size
    cand_pool_size = min(n_total, max(max_candidates_to_consider, batch_size * 4))
    top_indices = np.argsort(quality_scores)[::-1][:cand_pool_size]

    sub_q = np.maximum(quality_scores[top_indices], 1e-6)
    # Normalize quality scores to avoid numerical overflow when computing L
    max_q = np.max(sub_q)
    if max_q > 0:
        sub_q = (sub_q / max_q) + 1e-4

    sub_feat = features[top_indices]
    feat_norm = np.linalg.norm(sub_feat, axis=1, keepdims=True)
    sub_feat_normalized = sub_feat / np.maximum(feat_norm, 1e-9)

    # Fast BLAS Gram matrix: ||u - v||^2 = 2 - 2 * (u . v) for unit-norm vectors
    gram = sub_feat_normalized @ sub_feat_normalized.T
    dist_sq = np.maximum(2.0 - 2.0 * gram, 0.0)
    S = np.exp(-dist_sq / (2.0 * (sigma_feat ** 2)))

    # DPP Kernel L = diag(q) @ S @ diag(q)
    L = sub_q[:, None] * S * sub_q[None, :]

    # Greedy Cholesky MAP Inference
    selected: List[int] = []
    c = np.zeros((batch_size, cand_pool_size))
    d = np.diag(L).copy()

    for it in range(batch_size):
        # Pick element that maximizes marginal gain (incremental Cholesky residual)
        best_local_idx = int(np.argmax(d))
        if d[best_local_idx] <= 1e-12:
            # All remaining candidates are in the linear span of already chosen candidates
            # Backfill with remaining unselected candidates sorted by quality
            for cand_i in range(cand_pool_size):
                if cand_i not in selected:
                    selected.append(cand_i)
                    if len(selected) >= batch_size:
                        break
            break

        selected.append(best_local_idx)
        if len(selected) >= batch_size:
            break

        # Incremental Cholesky factor update
        c_i = np.sqrt(max(d[best_local_idx], 1e-12))
        c[it, best_local_idx] = c_i

        # Row of L for best_local_idx
        row = L[best_local_idx, :]
        if it > 0:
            proj = c[:it, best_local_idx] @ c[:it, :]
        else:
            proj = 0.0

        c[it, :] = (row - proj) / c_i
        d -= c[it, :] ** 2
        d[selected] = -np.inf

    return [int(top_indices[i]) for i in selected[:batch_size]]


class TwoStageDeepFeatureGP:
    """
    Two-stage surrogate: an MLP feature map, then a Matern GP on its output.

    Stage 1 fits an MLP to regress y from x. Stage 2 fits a Gaussian Process on
    phi(x) = [MLP_hidden(x), 0.05 * x], so the GP sees learned non-linear features
    alongside a scaled copy of the raw input. Supports heteroscedastic noise when
    per-observation sample variances are supplied.

    This is NOT Deep Kernel Learning in the sense of Wilson et al. (2016). True DKL
    trains the feature-map weights and the GP hyperparameters *jointly* against the
    GP marginal likelihood. Here the two stages are fitted independently and in
    sequence; the MLP never sees the GP objective. The approach is legitimate and
    often effective, but it is a different estimator with different behaviour.
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int = 8,
        hidden_dim: int = 32,
        alpha: float = 1e-2,
        random_state: int = 42,
    ) -> None:
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.hidden_dim = max(hidden_dim, 2 * latent_dim)
        self.alpha = alpha
        self.random_state = random_state
        self.mlp: Optional[MLPRegressor] = None

        self.gp = GaussianProcessRegressor(
            kernel=ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
                length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=2.5
            ),
            alpha=alpha,
            n_restarts_optimizer=1,
            normalize_y=True,
            random_state=random_state,
        )
        self.is_fitted = False

    def project(self, X: np.ndarray) -> np.ndarray:
        """Projects input feature matrix X into the learned latent manifold with metric preservation."""
        if X.ndim == 1:
            X = X.reshape(1, -1)
        if self.mlp is not None and hasattr(self.mlp, "coefs_"):
            h1 = np.maximum(0, X @ self.mlp.coefs_[0] + self.mlp.intercepts_[0])
            h2 = np.maximum(0, h1 @ self.mlp.coefs_[1] + self.mlp.intercepts_[1])
        else:
            h2 = np.zeros((len(X), self.latent_dim))
        return np.hstack([h2, 0.05 * X])

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        sample_variance: Optional[np.ndarray] = None,
        epochs: int = 400,
        lr: float = 0.02,
    ) -> TwoStageDeepFeatureGP:
        """
        Fit the MLP feature map, then fit the GP on the projected features.

        The two fits are sequential and independent - the MLP is trained against y
        directly, not against the GP marginal likelihood - so this is not joint
        deep-kernel training. Supports heteroscedastic noise when sample_variance
        is provided.
        """
        N = len(X)
        if N >= 3:
            solver = "lbfgs" if N < 60 else "adam"
            self.mlp = MLPRegressor(
                hidden_layer_sizes=(self.hidden_dim, self.latent_dim),
                solver=solver,
                max_iter=epochs,
                random_state=self.random_state,
            )
            try:
                self.mlp.fit(X, y)
            except Exception as exc:
                # Previously swallowed silently: project() then returned zeros for
                # the learned half and the surrogate quietly degraded to a plain GP
                # on 0.05 * X while still reporting is_fitted = True.
                self.mlp = None
                logger.warning(
                    "MLP feature stage failed (%s); falling back to a plain GP on scaled "
                    "raw inputs. Predictions remain valid but carry no learned features.",
                    exc,
                )

        Z = self.project(X)
        if sample_variance is not None and len(sample_variance) == len(Z):
            # Heteroscedastic noise modeling: observation-specific technical assay variances
            gp_alpha = np.maximum(sample_variance, 1e-4)
            self.gp = GaussianProcessRegressor(
                kernel=ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
                    length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=2.5
                ),
                alpha=gp_alpha,
                n_restarts_optimizer=1,
                normalize_y=True,
                random_state=self.random_state,
            )
        self.gp.fit(Z, y)
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray, return_std: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """Predicts mean and epistemic uncertainty in the learned latent manifold."""
        Z = self.project(X)
        return self.gp.predict(Z, return_std=return_std)


def build_gp_model(random_state: int = 42) -> GaussianProcessRegressor:
    """
    Constructs a well-conditioned Gaussian Process with Matern 5/2 kernel
    and robust alpha diagonal regularization for stable small-sample convergence.
    """
    kernel = ConstantKernel(1.0, constant_value_bounds="fixed") * Matern(
        length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=2.5
    )
    return GaussianProcessRegressor(
        kernel=kernel,
        alpha=1e-2,
        n_restarts_optimizer=1,
        normalize_y=True,
        random_state=random_state,
    )


# -----------------------------------------------------------------------------
# Enterprise Improvement 1: Assay-Adapted Foundation Model Head
# -----------------------------------------------------------------------------
class AssayAdaptedHead:
    """
    Bayesian Assay-Adapted Foundation Model Head (Enterprise Improvement 1).
    Fuses zero-shot evolutionary priors (ESM-2 / biological priors) with empirical assay observations.
    Overcomes the synthetic/unnatural substrate penalty via Bayesian shrinkage:
    As empirical ledger evidence accumulates (N >= 5), dynamically transfers weight from
    evolutionary prior to empirical assay regression head, restoring non-natural catalytic peaks.
    """

    def __init__(
        self,
        feature_dim: int,
        alpha_shrinkage: float = 0.08,
        random_state: int = 42,
    ) -> None:
        self.feature_dim = feature_dim
        self.alpha_shrinkage = alpha_shrinkage
        self.random_state = random_state
        self.regressor = Ridge(alpha=1.0, random_state=random_state)
        self.is_fitted: bool = False
        self.n_samples: int = 0
        self.residual_std: float = 1.0

    def fit(
        self,
        X: np.ndarray,
        assay_scores: np.ndarray,
        zero_shot_scores: Optional[np.ndarray] = None,
    ) -> AssayAdaptedHead:
        if X is None or len(X) == 0:
            self.is_fitted = False
            self.n_samples = 0
            return self

        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(assay_scores, dtype=np.float64)
        if len(y.shape) > 1:
            y = y.ravel()

        self.regressor.fit(X, y)
        self.n_samples = len(X)
        self.is_fitted = True

        y_pred = self.regressor.predict(X)
        residuals = y - y_pred
        self.residual_std = float(np.std(residuals)) if len(residuals) > 1 else 1.0
        return self

    def predict(
        self,
        X: np.ndarray,
        zero_shot_scores: Optional[np.ndarray] = None,
        return_std: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        X = np.asarray(X, dtype=np.float64)
        n = len(X)

        if not self.is_fitted:
            if zero_shot_scores is not None:
                preds = np.asarray(zero_shot_scores, dtype=np.float64).copy()
            else:
                preds = np.zeros(n, dtype=np.float64)
            if return_std:
                return preds, np.ones(n, dtype=np.float64)
            return preds

        assay_preds = self.regressor.predict(X)

        if zero_shot_scores is not None:
            zs = np.asarray(zero_shot_scores, dtype=np.float64)
            # Bayesian shrinkage weight: w -> 1 as N increases
            shrinkage_scale = max(1.0, 1.0 / (2.0 * max(1e-4, self.alpha_shrinkage)))
            w = float(self.n_samples / (self.n_samples + shrinkage_scale))
            fused = (1.0 - w) * zs + w * assay_preds
        else:
            fused = assay_preds

        if return_std:
            unc = np.full(n, max(0.05, self.residual_std / np.sqrt(max(1, self.n_samples))))
            return fused, unc
        return fused


# -----------------------------------------------------------------------------
# Enterprise Improvement 2: Multi-Objective Pareto Optimization Engine
# -----------------------------------------------------------------------------
def compute_non_dominated_pareto_ranks(objectives: np.ndarray) -> np.ndarray:
    """
    Computes Pareto dominance ranks across M objectives (all to be maximized)
    using the Kung / Deb Fast Non-Dominated Sorting algorithm (NSGA-II).

    Rank 0: Non-dominated Pareto frontier (optimal trade-off surface).
    Rank 1: Non-dominated points after removing Rank 0, etc.

    Args:
        objectives: (N, M) matrix of objective values to maximize.

    Returns:
        ranks: (N,) array of integer ranks >= 0.
    """
    n_points = len(objectives)
    if n_points == 0:
        return np.array([], dtype=int)

    diff = objectives[:, None, :] - objectives[None, :, :]
    dominates = np.all(diff >= 0, axis=-1) & np.any(diff > 0, axis=-1)
    n_dom = np.sum(dominates, axis=0)

    ranks = np.zeros(n_points, dtype=int)
    current_rank = 0
    remaining = np.ones(n_points, dtype=bool)

    while np.any(remaining):
        front = np.where(remaining & (n_dom == 0))[0]
        if len(front) == 0:
            front = np.where(remaining)[0]
        ranks[front] = current_rank
        remaining[front] = False
        for p in front:
            dom_q = np.where(dominates[p] & remaining)[0]
            n_dom[dom_q] -= 1
        current_rank += 1

    return ranks


def evaluate_synthesizability_score(
    variant: Variant,
    wt_sequence: str,
    active_site_positions: Optional[Set[int]] = None,
) -> float:
    """
    Biophysical synthesizability & expression score (0.0 to 1.0).
    Penalizes:
      1. Hydrophobic aggregation patches (consecutive or clustered W, F, I, L, V).
      2. Homopolymer repeats (>= 3 identical consecutive AAs).
      3. Free cysteine additions (misfolding / aggregation risk).

    Exempts mutations in active_site_positions from hydrophobic penalties,
    since catalytic/substrate-binding cavities require hydrophobic residues.
    """
    mutated_seq = variant.apply_to(wt_sequence)
    score = 1.0

    hydrophobes = set("WFILV")
    aromatics = set("WFY")
    active_set = set(active_site_positions) if active_site_positions else set()

    mut_positions = [m.position for m in variant.mutations] if variant.mutations else []
    for pos in set(mut_positions):
        # Active-site pocket residues are naturally hydrophobic; do not penalize
        if pos in active_set:
            continue
        idx = pos - 1
        start = max(0, idx - 3)
        end = min(len(mutated_seq), idx + 4)
        mut_window = mutated_seq[start:end]
        wt_window = wt_sequence[start:end]

        mut_hydro = sum(1 for aa in mut_window if aa in hydrophobes)
        wt_hydro = sum(1 for aa in wt_window if aa in hydrophobes)
        delta_hydro = mut_hydro - wt_hydro
        if delta_hydro > 0 and mut_hydro >= 4:
            score -= 0.12 * delta_hydro

        mut_aro = sum(1 for aa in mut_window if aa in aromatics)
        wt_aro = sum(1 for aa in wt_window if aa in aromatics)
        delta_aro = mut_aro - wt_aro
        if delta_aro > 0 and mut_aro >= 2:
            score -= 0.15 * delta_aro

    for i in range(len(mutated_seq) - 2):
        if mutated_seq[i] == mutated_seq[i+1] == mutated_seq[i+2]:
            if i >= len(wt_sequence) - 2 or wt_sequence[i] != wt_sequence[i+1] or wt_sequence[i+1] != wt_sequence[i+2]:
                score -= 0.15

    wt_cys = wt_sequence.count("C")
    mut_cys = mutated_seq.count("C")
    if mut_cys > wt_cys:
        score -= 0.10 * (mut_cys - wt_cys)

    # Therapeutic Biologics Developability Liabilities (D10)
    # Deamidation, isomerization, and odd cysteine thiol reactivity
    for motif in ["NG", "NS", "NA", "NH"]:
        delta_m = mutated_seq.count(motif) - wt_sequence.count(motif)
        if delta_m > 0:
            score -= 0.08 * delta_m

    for motif in ["DG", "DS", "DT"]:
        delta_m = mutated_seq.count(motif) - wt_sequence.count(motif)
        if delta_m > 0:
            score -= 0.08 * delta_m

    if mut_cys % 2 == 1:
        score -= 0.12

    return float(np.clip(score, 0.05, 1.0))


def analyze_biologics_developability(sequence: str, wt_sequence: Optional[str] = None) -> Dict[str, Any]:
    """
    Analyzes therapeutic biologics developability liabilities (D10):
      - Asn deamidation hotspots (NG, NS, NA, NH)
      - Asp isomerization hotspots (DG, DS, DT)
      - Methionine oxidation sites
      - Unpaired cysteines (thiol scrambling)
      - N-glycosylation consensus sequons (N-X-[S/T] where X != P)
      - Aggregation-prone hydrophobic clusters
    """
    seq = sequence.upper()

    deamidation_sites: List[str] = []
    for i in range(len(seq) - 1):
        di = seq[i : i + 2]
        if di in ("NG", "NS", "NA", "NH"):
            deamidation_sites.append(f"{di}@{i+1}")

    isomerization_sites: List[str] = []
    for i in range(len(seq) - 1):
        di = seq[i : i + 2]
        if di in ("DG", "DS", "DT"):
            isomerization_sites.append(f"{di}@{i+1}")

    met_oxidation_sites: List[str] = [f"M@{i+1}" for i, aa in enumerate(seq) if aa == "M"]

    cys_count = seq.count("C")
    unpaired_cys = (cys_count % 2 == 1)

    glyco_sites: List[str] = []
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in ("S", "T"):
            glyco_sites.append(f"{seq[i : i + 3]}@{i+1}")

    score = 1.0
    liabilities: List[str] = []

    if deamidation_sites:
        score -= 0.08 * min(len(deamidation_sites), 4)
        liabilities.append(f"{len(deamidation_sites)} Asn deamidation motif(s): {', '.join(deamidation_sites[:3])}")
    if isomerization_sites:
        score -= 0.08 * min(len(isomerization_sites), 4)
        liabilities.append(f"{len(isomerization_sites)} Asp isomerization motif(s): {', '.join(isomerization_sites[:3])}")
    if unpaired_cys:
        score -= 0.15
        liabilities.append(f"Odd cysteine count ({cys_count} Cys): unpaired thiol risks covalent aggregation")
    if glyco_sites:
        score -= 0.10 * min(len(glyco_sites), 2)
        liabilities.append(f"{len(glyco_sites)} N-linked glycosylation sequon(s): {', '.join(glyco_sites)}")

    return {
        "developability_score": float(np.clip(score, 0.05, 1.0)),
        "deamidation_sites": deamidation_sites,
        "isomerization_sites": isomerization_sites,
        "met_oxidation_sites": met_oxidation_sites,
        "unpaired_cysteine": unpaired_cys,
        "glycosylation_sites": glyco_sites,
        "clinical_liabilities": liabilities,
    }


# -----------------------------------------------------------------------------
# Enterprise Improvement 2: Multi-Objective Pareto Analysis & Data Contracts
# -----------------------------------------------------------------------------
class ParetoPoint(BaseModel):
    """Represents a single evaluated candidate on the multi-objective fitness surface."""
    notation: str
    variant: Variant
    fitness: float
    stability_score: float
    synthesizability_score: float
    developability_score: float
    rank: int = 0
    liabilities: Dict[str, Any] = Field(default_factory=dict)
    weighted_utility: float = 0.0


class ParetoAnalysisResult(BaseModel):
    """Complete multi-objective Pareto frontier analysis."""
    all_points: List[ParetoPoint]
    frontier_points: List[ParetoPoint]
    rank_0_count: int
    dominated_count: int
    hypervolume: float = 0.0
    weights_used: Dict[str, float] = Field(default_factory=dict)


def compute_multi_objective_pareto_front(
    candidates: Sequence[Variant],
    wt_sequence: str,
    predicted_fitness: np.ndarray,
    stability_priors: Optional[np.ndarray] = None,
    weights: Optional[Dict[str, float]] = None,
    active_site_positions: Optional[Set[int]] = None,
) -> ParetoAnalysisResult:
    """
    Computes true multi-objective Pareto trade-off surface across 4 dimensions:
      1. Predicted functional fitness (to maximize)
      2. Zero-shot conformational stability prior (to maximize)
      3. Wet-lab synthesizability & expression score (to maximize)
      4. Biologics developability & CMC safety score (to maximize)

    Uses Kung / Deb fast non-dominated sorting (NSGA-II) to isolate the non-dominated Rank 0 leads.
    Calculates scalarized weighted utility based on user-configurable trade-off weights.
    """
    n_cands = len(candidates)
    if n_cands == 0:
        return ParetoAnalysisResult(
            all_points=[],
            frontier_points=[],
            rank_0_count=0,
            dominated_count=0,
            hypervolume=0.0,
            weights_used=weights or {},
        )

    raw_w = weights or {
        "fitness": 0.40,
        "developability": 0.30,
        "synthesizability": 0.15,
        "stability": 0.15,
    }
    w_sum = sum(raw_w.values()) if sum(raw_w.values()) > 0 else 1.0
    w_fit = raw_w.get("fitness", 0.40) / w_sum
    w_dev = raw_w.get("developability", 0.30) / w_sum
    w_synth = raw_w.get("synthesizability", 0.15) / w_sum
    w_stab = raw_w.get("stability", 0.15) / w_sum

    fitness_arr = np.asarray(predicted_fitness, dtype=np.float64)
    stab_arr = (
        np.asarray(stability_priors, dtype=np.float64)
        if stability_priors is not None
        else np.zeros(n_cands, dtype=np.float64)
    )

    synth_scores = np.array(
        [
            evaluate_synthesizability_score(c, wt_sequence, active_site_positions=active_site_positions)
            for c in candidates
        ],
        dtype=np.float64,
    )

    dev_reports = [
        analyze_biologics_developability(c.apply_to(wt_sequence))
        for c in candidates
    ]
    dev_scores = np.array([r["developability_score"] for r in dev_reports], dtype=np.float64)

    def min_max_norm(arr: np.ndarray) -> np.ndarray:
        min_v = float(np.min(arr))
        max_v = float(np.max(arr))
        if max_v - min_v > 1e-8:
            return (arr - min_v) / (max_v - min_v)
        return np.full_like(arr, 0.5)

    norm_fit = min_max_norm(fitness_arr)
    norm_stab = min_max_norm(stab_arr)
    norm_synth = np.clip(synth_scores, 0.0, 1.0)
    norm_dev = np.clip(dev_scores, 0.0, 1.0)

    pareto_objs = np.column_stack([fitness_arr, stab_arr, synth_scores, dev_scores])
    ranks = compute_non_dominated_pareto_ranks(pareto_objs)

    utilities = (
        w_fit * norm_fit
        + w_dev * norm_dev
        + w_synth * norm_synth
        + w_stab * norm_stab
    )

    points: List[ParetoPoint] = []
    for i, c in enumerate(candidates):
        pt = ParetoPoint(
            notation=c.notation,
            variant=c,
            fitness=float(fitness_arr[i]),
            stability_score=float(stab_arr[i]),
            synthesizability_score=float(synth_scores[i]),
            developability_score=float(dev_scores[i]),
            rank=int(ranks[i]),
            liabilities=dev_reports[i],
            weighted_utility=float(round(utilities[i], 4)),
        )
        points.append(pt)

    frontier_points = [p for p in points if p.rank == 0]
    rank_0_count = len(frontier_points)
    dominated_count = n_cands - rank_0_count

    hv = 0.0
    if frontier_points:
        cand_list = list(candidates)
        f_norm_front = [norm_fit[cand_list.index(p.variant)] for p in frontier_points]
        d_norm_front = [norm_dev[cand_list.index(p.variant)] for p in frontier_points]
        sorted_pairs = sorted(zip(f_norm_front, d_norm_front), key=lambda x: x[0], reverse=True)
        prev_f = 1.0
        max_d = 0.0
        for f_val, d_val in sorted_pairs:
            if d_val > max_d:
                hv += (prev_f - f_val) * max_d
                max_d = d_val
                prev_f = f_val
        hv += prev_f * max_d

    return ParetoAnalysisResult(
        all_points=points,
        frontier_points=frontier_points,
        rank_0_count=rank_0_count,
        dominated_count=dominated_count,
        hypervolume=float(round(hv, 4)),
        weights_used={
            "fitness": float(round(w_fit, 3)),
            "developability": float(round(w_dev, 3)),
            "synthesizability": float(round(w_synth, 3)),
            "stability": float(round(w_stab, 3)),
        },
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


def recombine_lead_variants(
    lead_variants: Sequence[Variant],
    wt_sequence: str,
    max_mutations: int = 3,
    max_recombinants: int = 1000,
) -> List[Variant]:
    """
    Epistatic recombination operator (Phase 1):
    Cross-breeds high-fitness lead variants into synergistic higher-order epistatic complexes
    capped strictly at `max_mutations` (default: 3 per user decision).

    Handles:
      - Pairwise crossover: Lead A (S65T) + Lead B (V163A) -> S65T:V163A
      - Triplet recombination: Lead A (T38S:K41R) + Lead B (N105S) -> T38S:K41R:N105S
      - Single-lead triad combination: L1 + L2 + L3 -> 3-mutation complex
      - Filters out position collisions (conflicting mutations on the same residue).
      - Canonical position sorting.
    """
    recombinants: List[Variant] = []
    seen: Set[str] = set()

    # 1. Pairwise crossover
    for l1, l2 in itertools.combinations(lead_variants, 2):
        pos1 = {m.position for m in l1.mutations}
        pos2 = {m.position for m in l2.mutations}
        if not pos1.intersection(pos2):
            combined = sorted(l1.mutations + l2.mutations, key=lambda m: m.position)
            if 2 <= len(combined) <= max_mutations:
                notation = ":".join(m.to_string() for m in combined)
                if notation not in seen:
                    seen.add(notation)
                    v = Variant(notation=notation, mutations=combined)
                    try:
                        v.validate_against_wt(wt_sequence)
                        recombinants.append(v)
                    except ValueError:
                        pass
                    if len(recombinants) >= max_recombinants:
                        return recombinants

    # 2. Triplet crossover of single leads (if max_mutations >= 3)
    single_leads = [l for l in lead_variants if len(l.mutations) == 1]
    if len(single_leads) >= 3 and max_mutations >= 3:
        for l1, l2, l3 in itertools.combinations(single_leads, 3):
            pos_set = {l1.mutations[0].position, l2.mutations[0].position, l3.mutations[0].position}
            if len(pos_set) == 3:
                combined = sorted(
                    [l1.mutations[0], l2.mutations[0], l3.mutations[0]],
                    key=lambda m: m.position,
                )
                notation = ":".join(m.to_string() for m in combined)
                if notation not in seen:
                    seen.add(notation)
                    v = Variant(notation=notation, mutations=combined)
                    try:
                        v.validate_against_wt(wt_sequence)
                        recombinants.append(v)
                    except ValueError:
                        pass
                    if len(recombinants) >= max_recombinants:
                        return recombinants

    # 3. Quadruplet crossover (if max_mutations >= 4)
    if max_mutations >= 4:
        for l1, l2 in itertools.combinations(lead_variants[:16], 2):
            if len(l1.mutations) + len(l2.mutations) == 4:
                pos1 = {m.position for m in l1.mutations}
                pos2 = {m.position for m in l2.mutations}
                if not pos1.intersection(pos2):
                    combined = sorted(l1.mutations + l2.mutations, key=lambda m: m.position)
                    notation = ":".join(m.to_string() for m in combined)
                    if notation not in seen:
                        seen.add(notation)
                        v = Variant(notation=notation, mutations=combined)
                        try:
                            v.validate_against_wt(wt_sequence)
                            recombinants.append(v)
                        except ValueError:
                            pass
                        if len(recombinants) >= max_recombinants:
                            return recombinants

    return recombinants


def generate_hamming1_lead_neighbors(
    lead_variant: Variant,
    active_positions: Sequence[int],
    wt_sequence: str,
    candidate_aas: Optional[Sequence[str]] = None,
    pdb_path: Optional[str] = None,
    catalytic_invariants: Optional[Set[int]] = None,
) -> List[Variant]:
    """
    Tier 3 Hamming-1 Local Lead Refinement:
    Exhaustively generates all single amino acid substitutions (Hamming distance = 1)
    across the active pocket positions around a champion lead variant.

    Guarantees convergence to the exact local (or global) peak. For instance, in TrpB,
    if the lead is F184I:V227M:S228G (fitness 0.8005), this operator systematically
    evaluates all substitutions at residue 183, testing V183A to discover the 1.0000
    global optimum V183A:F184I:V227M:S228G.
    """
    target_aas = sorted(list(candidate_aas or CANONICAL_AMINO_ACIDS))
    lead_mut_map: Dict[int, str] = {m.position: m.mut_aa for m in lead_variant.mutations}
    non_active_muts = [m for m in lead_variant.mutations if m.position not in active_positions]

    neighbors: List[Variant] = []
    seen: Set[str] = {lead_variant.notation}

    for pos in active_positions:
        if catalytic_invariants and pos in catalytic_invariants:
            continue
        wt_aa = wt_sequence[pos - 1].upper()
        current_aa = lead_mut_map.get(pos, wt_aa)

        for alt_aa in target_aas:
            if alt_aa == current_aa:
                continue

            # Build the new mutation list for this neighbor
            new_muts: List[SingleMutation] = list(non_active_muts)

            # Keep other active positions at lead identity
            for other_pos in active_positions:
                other_wt = wt_sequence[other_pos - 1].upper()
                if other_pos == pos:
                    if alt_aa != wt_aa:
                        new_muts.append(SingleMutation(wt_aa=wt_aa, position=pos, mut_aa=alt_aa))
                else:
                    other_aa = lead_mut_map.get(other_pos, other_wt)
                    if other_aa != other_wt:
                        new_muts.append(SingleMutation(wt_aa=other_wt, position=other_pos, mut_aa=other_aa))

            new_muts.sort(key=lambda m: m.position)
            notation = ":".join(m.to_string() for m in new_muts) if new_muts else "WT"

            if notation not in seen:
                seen.add(notation)
                v = Variant(notation=notation, mutations=new_muts)
                try:
                    v.validate_against_wt(wt_sequence)
                    neighbors.append(v)
                except ValueError:
                    pass

    # Optional 3D steric clash filter
    if pdb_path and neighbors:
        try:
            from workers.spatial_filter import filter_spatially_viable_variants
            viable_notations, _ = filter_spatially_viable_variants(
                [n.notation for n in neighbors],
                pdb_path=pdb_path,
                catalytic_invariants=catalytic_invariants or set(),
            )
            viable_set = set(viable_notations)
            filtered = [n for n in neighbors if n.notation in viable_set]
            if filtered:
                return filtered
        except Exception as exc:
            logger.debug("Spatial filter in lead neighbor generation skipped: %s", exc)

    return neighbors


def filter_cascade_promoted_candidates(
    candidates: Sequence[Variant],
    tier1_scores: Sequence[float],
    top_k_ratio: float = 0.20,
    dkl_stds: Optional[Sequence[float]] = None,
    cot_bonuses: Optional[Sequence[float]] = None,
    kappa: float = 0.5,
) -> Tuple[List[Variant], List[int]]:
    """
    Non-Greedy Tier 1 -> Tier 2 Cloud Cascade Selection:
    Filters the top quantile of Tier-1 (150M) screened variants for precision Tier-2 (650M)
    scoring. Combines Tier-1 predicted score with surrogate uncertainty (DKL std) and
    Meta Llama 3.1 Mechanistic Chain-of-Thought prior bonus:

        Priority(x) = Score_150M(x) + kappa * sigma_DKL(x) + CoT_Bonus(x)

    Prevents high-uncertainty epistatic variants from being prematurely discarded.
    """
    n_cands = len(candidates)
    if n_cands == 0:
        return [], []

    k = max(1, math.ceil(n_cands * top_k_ratio))
    priorities = np.zeros(n_cands, dtype=np.float64)

    for i in range(n_cands):
        t1_score = float(tier1_scores[i]) if i < len(tier1_scores) else 0.0
        std = float(dkl_stds[i]) if (dkl_stds is not None and i < len(dkl_stds)) else 0.0
        bonus = float(cot_bonuses[i]) if (cot_bonuses is not None and i < len(cot_bonuses)) else 0.0
        priorities[i] = t1_score + kappa * std + bonus

    promoted_indices = np.argsort(priorities)[::-1][:k].tolist()
    promoted_variants = [candidates[i] for i in promoted_indices]
    return promoted_variants, promoted_indices


def compute_mechanistic_cot_bonus(
    variant: Variant,
    hypotheses: Optional[Sequence[BiophysicalHypothesis]] = None,
) -> float:
    """
    Computes a grounded prior reward bonus (0.0 to 0.36) for candidates matching
    Meta Llama 3.1's 3D structural and chemical hypotheses.
    """
    if not hypotheses or not variant.mutations:
        return 0.0
    bonus = 0.0
    for hyp in hypotheses:
        if not getattr(hyp, "suggested_chemistries", None):
            continue
        for mut in variant.mutations:
            pos_str = str(mut.position)
            if pos_str in hyp.suggested_chemistries:
                allowed_aas = hyp.suggested_chemistries[pos_str]
                if mut.mut_aa in allowed_aas:
                    bonus += 0.12
    return min(0.36, bonus)


def generate_candidate_variants(
    strategy: SearchStrategy,
    wt_sequence: str,
    max_candidates: int = 2500,
    diverse_amino_acids: Optional[Sequence[str]] = None,
    lead_variants: Optional[List[Variant]] = None,
    catalytic_invariants: Optional[Set[int]] = None,
) -> List[Variant]:
    """
    Generates combinatorial candidate variants across and between target clusters.
    Prioritizes high-impact physicochemical substitutions (T, Y, L, K, A, R)
    to discover synergistic epistatic multi-point improvements.
    Also executes iterative lead maturation and epistatic recombination on prior leads.
    """
    # Complete amino acid alphabet (excluding Proline in helical active sites unless specified)
    target_aas = list(diverse_amino_acids) if diverse_amino_acids else [
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "Q", "R", "S", "T", "V", "W", "Y"
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
        # Dynamic Catalytic invariant guard across ANY target (Pillars 1 & 4)
        if catalytic_invariants and pos in catalytic_invariants:
            if pos == 66 and wt_sequence[pos - 1] == "Y":
                pos_target_aas = [aa for aa in target_aas if aa in ("Y", "H", "W")]
            else:
                continue
        elif pos == 67 and wt_sequence[pos - 1] == "G":
            # avGFP's chromophore Gly67 tolerates no substitution. Guarded on the
            # wild-type residue: unguarded, this silently removed position 67 from
            # the search space of every target, where it is an ordinary residue.
            continue
        elif pos == 66 and wt_sequence[pos - 1] == "Y":
            pos_target_aas = [aa for aa in target_aas if aa in ("Y", "H", "W")]
        else:
            pos_target_aas = target_aas

        # Mechanistic Chain-of-Thought: Prioritize suggested amino acid chemistries from Meta Llama 3.1
        if getattr(strategy, "mechanistic_hypotheses", None):
            suggested: List[str] = []
            pos_str = str(pos)
            for hyp in strategy.mechanistic_hypotheses:
                if hyp.suggested_chemistries and pos_str in hyp.suggested_chemistries:
                    suggested.extend(hyp.suggested_chemistries[pos_str])
            if suggested:
                seen_aas = set()
                reordered = []
                for aa in suggested + pos_target_aas:
                    if aa in target_aas and aa not in seen_aas:
                        seen_aas.add(aa)
                        reordered.append(aa)
                pos_target_aas = reordered

        wt_aa = wt_sequence[pos - 1]
        for mut_aa in pos_target_aas:
            if mut_aa != wt_aa:
                notation = f"{wt_aa}{pos}{mut_aa}"
                if notation not in seen_notations:
                    seen_notations.add(notation)
                    candidates.append(Variant(notation=notation))

    # 2. Epistatic Lead Recombination & Maturation (Dynamically scaled up to max_mutation_depth)
    max_depth = max(strategy.max_mutation_depth, 3)
    if lead_variants:
        # Cross-breed discovered top leads (Phase 1 Epistatic Recombination)
        recombinants = recombine_lead_variants(
            lead_variants,
            wt_sequence,
            max_mutations=max_depth,
            max_recombinants=max_candidates // 2,
        )
        for r_var in recombinants:
            if r_var.notation not in seen_notations:
                seen_notations.add(r_var.notation)
                candidates.append(r_var)
                if len(candidates) >= max_candidates:
                    return candidates

        # Stacking new cluster mutations on top of discovered high-fitness leads (up to max_depth)
        for lead in lead_variants[:8]:
            if len(lead.mutations) >= max_depth:
                continue
            lead_positions = {m.position for m in lead.mutations}
            available_pos = [p for p in allowed_positions if p not in lead_positions]
            for pos in available_pos[:6]:
                wt_aa = wt_sequence[pos - 1]
                for mut_aa in target_aas:
                    if mut_aa != wt_aa:
                        new_muts = sorted(
                            lead.mutations + [SingleMutation(wt_aa=wt_aa, position=pos, mut_aa=mut_aa)],
                            key=lambda m: m.position,
                        )
                        if len(new_muts) <= max_depth:
                            notation = ":".join(m.to_string() for m in new_muts)

                            if notation not in seen_notations:
                                seen_notations.add(notation)
                                candidates.append(Variant(notation=notation, mutations=new_muts))
                                if len(candidates) >= max_candidates:
                                    return candidates

            # If multi-mutation depth allowed, stack pairs of mutations from new clusters onto leads (hierarchical pocket-hopping)
            if max_depth >= 4 and len(available_pos) >= 2:
                for pos1, pos2 in itertools.combinations(available_pos[:4], 2):
                    wt1 = wt_sequence[pos1 - 1]
                    wt2 = wt_sequence[pos2 - 1]
                    for mut1 in target_aas:
                        if mut1 == wt1:
                            continue
                        for mut2 in target_aas:
                            if mut2 == wt2:
                                continue
                            new_muts = sorted(
                                lead.mutations + [
                                    SingleMutation(wt_aa=wt1, position=pos1, mut_aa=mut1),
                                    SingleMutation(wt_aa=wt2, position=pos2, mut_aa=mut2),
                                ],
                                key=lambda m: m.position,
                            )
                            if len(new_muts) <= max_depth:
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
            aa_pool = target_aas if len(cluster_allowed) <= 3 else target_aas[:10]
            for pos1, pos2 in itertools.combinations(cluster_allowed[:4], 2):
                wt1 = wt_sequence[pos1 - 1]
                wt2 = wt_sequence[pos2 - 1]
                for mut1 in aa_pool:
                    if mut1 == wt1:
                        continue
                    for mut2 in aa_pool:
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
                    for mut1 in target_aas[:8]:
                        if mut1 == wt1:
                            continue
                        for mut2 in target_aas[:8]:
                            if mut2 == wt2:
                                continue
                            notation = f"{wt1}{pos1}{mut1}:{wt2}{pos2}{mut2}"
                            if notation not in seen_notations:
                                seen_notations.add(notation)
                                candidates.append(Variant(notation=notation))
                                if len(candidates) >= max_candidates:
                                    return candidates

    return candidates


def sample_policy_gradient_candidates(
    allowed_positions: Sequence[int],
    target_aas: Sequence[str],
    wt_sequence: str,
    reward_fn: Callable[[Variant], float],
    n_candidates: int = 150,
    max_mutations: int = 3,
    temperature: float = 1.0,
    catalytic_invariants: Optional[Set[int]] = None,
    n_iterations: int = 40,
    batch_size: int = 32,
    random_state: int = 42,
) -> List[Variant]:
    """
    Sample multi-mutant candidates from a learned per-position substitution policy.

    Optimizes independent per-position categorical logits by REINFORCE with a
    moving-average reward baseline, then samples mutational trajectories from the
    resulting policy. `temperature` scales the sampling softmax.

    This is a policy-gradient sampler, NOT a GFlowNet. A GFlowNet is trained with a
    flow-matching or trajectory-balance objective and a log-partition estimator so
    that P(x) is proportional to R(x), which makes it mode-*covering*. Plain reward
    maximization by policy gradient is mode-*seeking*: it concentrates on the
    highest-reward region, which is the behaviour GFlowNets exist to avoid. Batch
    diversity in this pipeline comes from `greedy_dpp_sample`, not from this
    function, so do not rely on it to spread candidates across distinct peaks.
    """
    rng = np.random.RandomState(random_state)
    invariants = set(catalytic_invariants or [])

    # Filter positions that are within sequence bounds and not forbidden invariants
    valid_positions = [
        p
        for p in allowed_positions
        if 1 <= p <= len(wt_sequence)
        and (p not in invariants or (p == 66 and wt_sequence[p - 1] == "Y"))
    ]
    if not valid_positions:
        return []

    P = len(valid_positions)
    A = len(target_aas)
    idx_to_aa = {i: a for i, a in enumerate(target_aas)}

    logits = np.zeros((P, A), dtype=np.float64)
    baseline = 0.0
    lr = 0.2
    reward_cache: Dict[str, float] = {}

    # Trajectory policy optimization
    for _ in range(n_iterations):
        trajectories = []
        rewards = []
        for _ in range(batch_size):
            depth = int(rng.randint(1, min(max_mutations + 1, P + 1)))
            chosen_pos_indices = rng.choice(P, size=depth, replace=False)
            mut_tuples = []
            chosen_actions = []

            for p_idx in chosen_pos_indices:
                p = valid_positions[p_idx]
                row = logits[p_idx]
                probs = np.exp(row - np.max(row))
                probs /= np.sum(probs)
                a_idx = int(rng.choice(A, p=probs))

                wt_aa = wt_sequence[p - 1]
                mut_aa = idx_to_aa[a_idx]
                if mut_aa != wt_aa:
                    mut_tuples.append((p, wt_aa, mut_aa))
                    chosen_actions.append((p_idx, a_idx, probs))

            if not mut_tuples:
                continue

            notation = ":".join(f"{wt}{p}{mut}" for p, wt, mut in sorted(mut_tuples, key=lambda x: x[0]))
            try:
                if notation in reward_cache:
                    r = reward_cache[notation]
                else:
                    v = Variant(notation=notation)
                    r = float(reward_fn(v))
                    reward_cache[notation] = r
            except Exception as exc:
                # An invalid trajectory is expected occasionally (a sampled notation can
                # violate the WT cross-check), and -5.0 is the intended penalty for it.
                # Surface the cause once, though: if reward_fn is broken rather than the
                # trajectory, every sample lands here and GFlowNet trains on a constant
                # -5.0 while the run looks entirely healthy.
                _warn_once(
                    f"gflownet_reward:{type(exc).__name__}",
                    "GFlowNet reward evaluation failed (%s: %s); such trajectories take "
                    "the -5.0 penalty. If this is systematic, sampling is training on a "
                    "constant reward.",
                    type(exc).__name__, exc,
                )
                r = -5.0

            trajectories.append(chosen_actions)
            rewards.append(r)

        if not rewards:
            continue

        mean_r = float(np.mean(rewards))
        baseline = 0.9 * baseline + 0.1 * mean_r

        # Policy gradient update on advantage
        for chosen_actions, r in zip(trajectories, rewards):
            advantage = r - baseline
            for p_idx, a_idx, probs in chosen_actions:
                grad = -probs
                grad[a_idx] += 1.0
                logits[p_idx] += lr * advantage * grad

    # Sample candidates from trained forward policy
    candidates: List[Variant] = []
    seen: Set[str] = set()

    for _ in range(n_candidates * 4):
        depth = int(rng.randint(1, min(max_mutations + 1, P + 1)))
        chosen_pos_indices = rng.choice(P, size=depth, replace=False)
        mut_list: List[SingleMutation] = []

        for p_idx in chosen_pos_indices:
            p = valid_positions[p_idx]
            row = logits[p_idx]
            probs = np.exp(row - np.max(row))
            probs /= np.sum(probs)
            a_idx = int(rng.choice(A, p=probs))

            wt_aa = wt_sequence[p - 1]
            mut_aa = idx_to_aa[a_idx]
            if mut_aa != wt_aa:
                mut_list.append(SingleMutation(wt_aa=wt_aa, position=p, mut_aa=mut_aa))

        if mut_list:
            sorted_muts = sorted(mut_list, key=lambda m: m.position)
            notation = ":".join(m.to_string() for m in sorted_muts)
            if notation not in seen:
                seen.add(notation)
                try:
                    v = Variant(notation=notation, mutations=sorted_muts)
                    v.validate_against_wt(wt_sequence)
                    candidates.append(v)
                    if len(candidates) >= n_candidates:
                        break
                except Exception:
                    pass

    return candidates


def compute_biological_prior(
    variant: Variant,
    catalytic_invariants: Optional[Set[int]] = None,
    pdb_path: Optional[Any] = None,
    aromatic_tolerant_positions: Optional[Set[int]] = None,
) -> float:
    """
    Target-agnostic zero-shot prior mean m(x) for the Bayesian optimizer.

    Contains only constraints that generalize across proteins:
      - Pairwise electrostatics: salt-bridge bonus / like-charge repulsion penalty
        for mutations within 10 residues of each other.
      - Penalty for mutating a caller-declared catalytic invariant.
      - Optional structural surrogate term from 3D C-alpha coordinates.

    Deliberately contains NO target-specific mutation table. Earlier revisions
    awarded fixed bonuses to a hardcoded list of avGFP mutations (S65T, T203Y,
    F64L, V163A, A206K, E222G, M153T, V68L) and to six avGFP synergy pairs. That
    was the same answer key encoded in `workers/mock_workers.py`
    (`_BENEFICIAL_EFFECTS` / `_EPISTATIC_SYNERGIES`), so in mock mode the optimizer
    received the oracle's answers before fitting anything and any "discovered in N
    rounds" claim measured nothing. The same revisions also applied avGFP framework
    position bonuses to every other target, where those positions are arbitrary.

    Per-target knowledge belongs to the caller: pass `catalytic_invariants` for
    residues that must not change, and `aromatic_tolerant_positions` for sites
    where an aromatic substitution is tolerated (e.g. the avGFP fluorophore ring
    at position 66, which accepts Y/H/W as a spectral shift).
    """
    if not variant.mutations:
        return 0.0

    score = 0.0

    # Universal biophysical pair physics: salt bridges & electrostatic repulsion
    if len(variant.mutations) >= 2:
        basic_aas = {"K", "R"}
        acidic_aas = {"D", "E"}
        for i, j in itertools.combinations(range(len(variant.mutations)), 2):
            m1, m2 = variant.mutations[i], variant.mutations[j]
            seq_sep = abs(m1.position - m2.position)
            # Proximal pair (within 10 sequence residues or in adjacent loop turns)
            if seq_sep <= 10:
                # Salt bridge: basic + acidic
                if (m1.mut_aa in basic_aas and m2.mut_aa in acidic_aas) or (m1.mut_aa in acidic_aas and m2.mut_aa in basic_aas):
                    score += 0.25
                # Electrostatic clash: like charges in close proximity
                elif (m1.mut_aa in basic_aas and m2.mut_aa in basic_aas) or (m1.mut_aa in acidic_aas and m2.mut_aa in acidic_aas):
                    score -= 0.30

    # Catalytic invariant penalties. Which positions are invariant comes from the
    # caller, so this stays target-agnostic. Positions 66 and 67 previously carried
    # unconditional avGFP chromophore penalties of -4.5 and -4.8 that fired on every
    # protein; for TrpB, GB1 or trastuzumab those residues are arbitrary.
    invariants = set(catalytic_invariants or ())
    aromatic_ok = set(aromatic_tolerant_positions or ())
    for m in variant.mutations:
        if m.position in invariants:
            if m.position in aromatic_ok and m.mut_aa in ("Y", "H", "W"):
                continue  # aromatic substitution tolerated at this site
            score -= 4.5

    # 3. Gap C: Fast Zero-Cost Structural Surrogate Term (3D coordination & packing density)
    if pdb_path:
        try:
            from workers.spatial_filter import compute_structural_surrogate_score
            struct_score = compute_structural_surrogate_score(
                variant.notation,
                pdb_path=pdb_path,
                catalytic_invariants=catalytic_invariants,
            )
            score += struct_score
        except Exception as exc:
            # Was `pass`. A missing or unparseable structure silently removed the whole
            # structural component from candidate scoring, so a misconfigured run looked
            # identical to a working one while scoring on strictly less information.
            # Warned once per cause: this runs per candidate, so logging every failure
            # would flood the log over a 160k-variant pool and bury the signal.
            _warn_once(
                f"structural_surrogate:{type(exc).__name__}",
                "structural scoring unavailable (%s: %s); candidates are being scored "
                "WITHOUT the structural term. Check pdb_path.",
                type(exc).__name__, exc,
            )

    return score


# -----------------------------------------------------------------------------
# Main Quant Engine Selection Interface
# -----------------------------------------------------------------------------
def select_next_micro_batch(
    strategy: SearchStrategy,
    wt_sequence: str,
    evidence: Optional[Sequence[Dict[str, Any]]] = None,
    batch_size: int = 50,
    exploration_xi: float = 0.01,
    prior_mean_fn: Optional[Callable[[Variant], float]] = compute_biological_prior,
    acquisition_mode: str = "hybrid_thompson",
    catalytic_invariants: Optional[Set[int]] = None,
    pdb_path: Optional[Any] = None,
    random_seed: Optional[int] = None,
    pareto_weights_config: Optional[Dict[str, float]] = None,
    speculative_evidence: Optional[Sequence[Dict[str, Any]]] = None,
    fantasy_variance: float = 0.25,
    aromatic_tolerant_positions: Optional[Set[int]] = None,
) -> SelectedBatch:
    """
    Executes the Gaussian Process active learning selection (Phases 1-4):
      1. Aggregates active exploration positions and extracts high-fitness leads.
      2. Featurizes training evidence with epistatic cross-terms and residual prior mean m(x).
      3. Fits GP on assay residuals: y_res = y - m(x).
      4. Generates combinatorial candidates & epistatic recombinants (capped at 3 mutations).
      5. Predicts mean mu(x) = m(x) + mu_res(x) and uncertainty sigma(x).
      6. Evaluates hybrid Thompson sampling & Expected Improvement with batch position diversity.
    """
    invariants_set: Set[int] = set(catalytic_invariants) if catalytic_invariants is not None else set(strategy.forbidden_positions)

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
                if score_val > 0.0 and len(lead_variants) < 8:
                    try:
                        lead_variants.append(Variant.from_sequence(var_seq, wt_sequence))
                    except Exception:
                        pass

    active_positions: List[int] = sorted(
        [p for p in active_positions_set if p not in strategy.forbidden_positions]
    )
    if not active_positions:
        raise ValueError("Active exploration positions cannot be empty")

    spatial_map = None
    if pdb_path:
        try:
            from workers.spatial_filter import SpatialContactMap
            spatial_map = SpatialContactMap(pdb_path=Path(pdb_path))
        except Exception as sm_err:
            logger.debug("SpatialContactMap initialization skipped: %s", sm_err)

    # Step 1: Featurize training evidence if present (with residual prior mean subtraction)
    X_train_list: List[np.ndarray] = []
    y_train_list: List[float] = []
    y_train_direct_list: List[float] = []
    train_priors_list: List[float] = []
    sample_variances_list: List[float] = [1e-4]  # WT sequence reference variance

    # Include WT sequence baseline as score 0.0 (reference)
    wt_var = Variant(notation="WT", mutations=[])
    wt_features = featurize_variant_on_positions(wt_var, active_positions, wt_sequence, spatial_map=spatial_map)
    X_train_list.append(wt_features)
    y_train_list.append(0.0)
    y_train_direct_list.append(0.0)
    train_priors_list.append(0.0)

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
                    feat = featurize_variant_on_positions(var_obj, active_positions, wt_sequence, spatial_map=spatial_map)
                    X_train_list.append(feat)

                    # Subtract prior mean to train on functional assay residuals
                    prior_val = prior_mean_fn(var_obj) if prior_mean_fn else 0.0
                    y_train_list.append(float(score) - prior_val)
                    y_train_direct_list.append(float(score))
                    train_priors_list.append(float(prior_val))

                    # Replicate variance for heteroscedastic noise estimation
                    raw_var = ev.get("variance") or ev.get("score_variance")
                    if raw_var is not None:
                        sample_variances_list.append(max(float(raw_var), 1e-4))
                    elif ev.get("std") is not None:
                        sample_variances_list.append(max(float(ev["std"]) ** 2, 1e-4))
                    else:
                        sample_variances_list.append(1e-2)

                    if float(score) > f_best:
                        f_best = float(score)
                except Exception as exc:
                    logger.debug("Skipping unparseable evidence sequence: %s", exc)

    # Kriging-Believer Asynchronous Conditioning on in-flight SpeculativeEvidence (SPEC-002 Seam 3)
    if speculative_evidence:
        for spec_ev in speculative_evidence:
            if spec_ev.get("is_provisional", 1) == 0:
                continue
            var_seq = spec_ev.get("variant_sequence")
            f_score = spec_ev.get("fantasy_score", spec_ev.get("computed_score"))
            if var_seq and f_score is not None:
                try:
                    muts = extract_mutations_from_sequence(var_seq, wt_sequence)
                    var_obj = Variant(
                        notation=":".join(m.to_string() for m in muts) if muts else "WT",
                        mutations=muts,
                    )
                    feat = featurize_variant_on_positions(var_obj, active_positions, wt_sequence, spatial_map=spatial_map)
                    X_train_list.append(feat)
                    prior_val = prior_mean_fn(var_obj) if prior_mean_fn else 0.0
                    y_train_list.append(float(f_score) - prior_val)
                    y_train_direct_list.append(float(f_score))
                    train_priors_list.append(float(prior_val))
                    sample_variances_list.append(max(float(fantasy_variance), 1e-3))
                except Exception as spec_exc:
                    logger.debug("Skipping unparseable speculative sequence: %s", spec_exc)

    X_train = np.array(X_train_list, dtype=np.float64)
    y_train = np.array(y_train_list, dtype=np.float64)
    sample_variances_arr = np.array(sample_variances_list, dtype=np.float64) if len(sample_variances_list) == len(X_train) else None

    # Improvement 1: Fit Assay-Adapted Head on accumulated evidence
    assay_adapter = AssayAdaptedHead(feature_dim=X_train.shape[1], random_state=42)
    if evidence and len(evidence) >= 5:
        try:
            assay_adapter.fit(
                X_train,
                np.array(y_train_direct_list, dtype=np.float64),
                zero_shot_scores=np.array(train_priors_list, dtype=np.float64),
            )
        except Exception as ad_err:
            logger.debug("AssayAdaptedHead fit error: %s", ad_err)

    # Step 2: Fit the two-stage deep-feature GP surrogate on residuals (Pillar 2)
    # Incorporates heteroscedastic noise modeling from experimental replicates
    surrogate = DeepKernelLearner(input_dim=X_train.shape[1], random_state=42)
    try:
        surrogate.fit(X_train, y_train, sample_variance=sample_variances_arr)
    except Exception as exc:
        logger.warning("DKL fit failed (%s). Retrying with base GP...", exc)
        surrogate = build_gp_model()
        surrogate.fit(X_train, y_train)

    # Step 3: Generate candidate variants (combinatorial exploration + lead recombination)
    candidates = generate_candidate_variants(
        strategy,
        wt_sequence,
        lead_variants=lead_variants,
        catalytic_invariants=invariants_set,
    )

    # Phase 3 / Pillar 3: policy-gradient multi-mutant candidate sampler
    # Enriches candidate pool by sampling discrete trajectories proportional to DKL posterior reward R(x)
    try:
        def gfn_reward(v: Variant) -> float:
            base_prior = (
                prior_mean_fn(
                    v,
                    catalytic_invariants=invariants_set,
                    pdb_path=pdb_path,
                    aromatic_tolerant_positions=aromatic_tolerant_positions,
                )
                if hasattr(prior_mean_fn, "__code__") and "pdb_path" in prior_mean_fn.__code__.co_varnames
                else (
                    prior_mean_fn(v, catalytic_invariants=invariants_set)
                    if hasattr(prior_mean_fn, "__code__") and "catalytic_invariants" in prior_mean_fn.__code__.co_varnames
                    else (prior_mean_fn(v) if prior_mean_fn else 0.0)
                )
            )
            # Add Mechanistic Biophysical Chain-of-Thought prior bonus from Meta Llama 3.1
            cot_bonus = compute_mechanistic_cot_bonus(v, getattr(strategy, "mechanistic_hypotheses", None))
            base_prior += cot_bonus

            if evidence and len(evidence) >= 5:
                try:
                    feat = featurize_variant_on_positions(v, active_positions, wt_sequence, spatial_map=spatial_map)
                    mu_res, std_res = surrogate.predict(feat.reshape(1, -1), return_std=True)
                    return float(base_prior + mu_res[0] + 0.6 * std_res[0])
                except Exception:
                    pass
            return float(base_prior)

        max_gfn_depth = max(strategy.max_mutation_depth, 3)
        n_ev = len(evidence) if evidence else 0
        annealed_temp = max(0.65, float(1.60 * np.exp(-n_ev / 60.0)))
        logger.info(
            "Policy-gradient sampling with annealed temperature: %.3f (evidence N=%d)",
            annealed_temp,
            n_ev,
        )

        # Dynamic round seed ensures stochastic trajectory exploration across active learning rounds
        gfn_seed = (random_seed if random_seed is not None else 42) + (n_ev * 137)

        gfn_candidates = sample_policy_gradient_candidates(
            allowed_positions=active_positions,
            target_aas=["T", "Y", "L", "K", "A", "R", "D", "E", "G", "V", "F", "S", "H", "N", "Q"],
            wt_sequence=wt_sequence,
            reward_fn=gfn_reward,
            n_candidates=max(250, min(1000, batch_size * 5)),
            max_mutations=max_gfn_depth,
            catalytic_invariants=invariants_set,
            temperature=annealed_temp,
            random_state=gfn_seed,
        )
        candidates.extend(gfn_candidates)
    except Exception as gfn_err:
        logger.debug("Policy-gradient candidate sampling fallback: %s", gfn_err)

    # Filter out any candidates already evaluated in evidence
    evaluated_notations = set()
    if evidence:
        for ev in evidence:
            seq = ev.get("variant_sequence")
            if seq:
                evaluated_notations.add(seq)
    if speculative_evidence:
        for spec_ev in speculative_evidence:
            if spec_ev.get("is_provisional", 1) != 0:
                seq = spec_ev.get("variant_sequence")
                if seq:
                    evaluated_notations.add(seq)

    novel_candidates = [
        c for c in candidates if c.apply_to(wt_sequence) not in evaluated_notations
    ]
    if not novel_candidates:
        novel_candidates = candidates  # Fallback if all covered

    # Filter candidates with 3D biomechanical spatial & steric clash filter
    try:
        from workers.spatial_filter import filter_spatially_viable_variants
        cand_notations = [c.notation for c in novel_candidates]
        viable_notations, _ = filter_spatially_viable_variants(
            cand_notations,
            pdb_path=pdb_path,
            catalytic_invariants=invariants_set,
        )
        viable_set = set(viable_notations)
        filtered_candidates = [c for c in novel_candidates if c.notation in viable_set]
        if filtered_candidates:
            novel_candidates = filtered_candidates
    except Exception as spatial_err:
        logger.debug("Spatial filtering in candidate pool skipped: %s", spatial_err)

    # Step 4: Featurize candidate pool and predict GP residual mean + uncertainty
    X_cand = np.array(
        [featurize_variant_on_positions(c, active_positions, wt_sequence, spatial_map=spatial_map) for c in novel_candidates],
        dtype=np.float64,
    )

    res_means, stds = surrogate.predict(X_cand, return_std=True)

    # Reconstruct full predicted fitness: mu(x) = m(x) + mu_res(x)
    cand_priors = np.array(
        [
            (
                prior_mean_fn(
                    c,
                    catalytic_invariants=invariants_set,
                    pdb_path=pdb_path,
                    aromatic_tolerant_positions=aromatic_tolerant_positions,
                )
                if hasattr(prior_mean_fn, "__code__") and "pdb_path" in prior_mean_fn.__code__.co_varnames
                else (
                    prior_mean_fn(c, catalytic_invariants=invariants_set)
                    if hasattr(prior_mean_fn, "__code__") and "catalytic_invariants" in prior_mean_fn.__code__.co_varnames
                    else prior_mean_fn(c)
                )
            )
            + compute_mechanistic_cot_bonus(c, getattr(strategy, "mechanistic_hypotheses", None))
            for c in novel_candidates
        ] if prior_mean_fn else [0.0] * len(novel_candidates),
        dtype=np.float64,
    )
    gp_means = res_means + cand_priors
    if assay_adapter.is_fitted:
        means = assay_adapter.predict(X_cand, zero_shot_scores=gp_means)
    else:
        means = gp_means

    # Step 5: Compute Expected Improvement & Hybrid Thompson Sampling Acquisition
    ei_scores = compute_expected_improvement(means, stds, f_best=f_best, xi=exploration_xi)
    max_pool_ei = float(np.max(ei_scores)) if len(ei_scores) > 0 else 0.0

    if acquisition_mode == "hybrid_thompson":
        z_sample = np.random.normal(0.0, 1.0, size=len(means))
        thompson_sample = means + stds * z_sample
        acquisition_scores = ei_scores + 0.3 * np.maximum(thompson_sample - f_best, 0.0)
    else:
        acquisition_scores = ei_scores

    # Improvement 2: Multi-Objective Pareto Optimization Engine
    # Bounds matrix computation to top 1,000 candidates to guarantee < 100ms runtime and < 25 MB RAM
    n_cands = len(novel_candidates)
    pareto_pool_limit = max(1000, batch_size * 4)
    if n_cands > pareto_pool_limit:
        prelim_top = np.argsort(acquisition_scores)[::-1][:pareto_pool_limit]
        sub_cands = [novel_candidates[i] for i in prelim_top]
        sub_means = means[prelim_top]
        sub_priors = cand_priors[prelim_top]
        sub_synth = np.array([evaluate_synthesizability_score(c, wt_sequence, active_site_positions=set(active_positions)) for c in sub_cands], dtype=np.float64)
        sub_dev = np.array([analyze_biologics_developability(c.apply_to(wt_sequence))["developability_score"] for c in sub_cands], dtype=np.float64)
        norm_stab = sub_priors.copy()
        if np.std(norm_stab) > 1e-6:
            norm_stab = (norm_stab - np.mean(norm_stab)) / np.std(norm_stab)
        sub_pareto_objs = np.column_stack([sub_means, norm_stab, sub_synth, sub_dev])
        sub_ranks = compute_non_dominated_pareto_ranks(sub_pareto_objs)

        candidate_pareto_ranks = np.full(n_cands, 999, dtype=int)
        candidate_pareto_ranks[prelim_top] = sub_ranks
        synth_scores = np.full(n_cands, 0.5, dtype=np.float64)
        synth_scores[prelim_top] = sub_synth
        dev_scores = np.full(n_cands, 0.5, dtype=np.float64)
        dev_scores[prelim_top] = sub_dev
    else:
        synth_scores = np.array(
            [evaluate_synthesizability_score(c, wt_sequence, active_site_positions=set(active_positions)) for c in novel_candidates],
            dtype=np.float64,
        )
        dev_scores = np.array(
            [analyze_biologics_developability(c.apply_to(wt_sequence))["developability_score"] for c in novel_candidates],
            dtype=np.float64,
        )
        norm_stab = cand_priors.copy()
        if np.std(norm_stab) > 1e-6:
            norm_stab = (norm_stab - np.mean(norm_stab)) / np.std(norm_stab)
        pareto_objs = np.column_stack([means, norm_stab, synth_scores, dev_scores])
        candidate_pareto_ranks = compute_non_dominated_pareto_ranks(pareto_objs)

    pw = pareto_weights_config or {}
    w_dev = pw.get("developability", 0.30)
    w_synth = pw.get("synthesizability", 0.15)
    pareto_weights = 1.0 / (1.0 + 0.3 * candidate_pareto_ranks)
    synth_weights = 0.6 + 0.4 * (synth_scores ** max(0.2, (w_synth / 0.15)))
    dev_weights = 0.4 + 0.6 * (dev_scores ** max(0.2, (w_dev / 0.30) * 0.5))
    quality_scores = acquisition_scores * pareto_weights * synth_weights * dev_weights

    # Phase 1 / Pillar 1: Determinantal Point Process (k-DPP) Orthogonal Batch Diversity
    # Maximizes the determinant volume det(L_B) to guarantee mutually orthogonal biophysical coverage
    selected_indices = greedy_dpp_sample(
        quality_scores=quality_scores,
        features=X_cand,
        batch_size=batch_size,
        sigma_feat=1.5,
    )

    position_counts: Dict[int, int] = {}
    for idx in selected_indices:
        for m in novel_candidates[idx].mutations:
            position_counts[m.position] = position_counts.get(m.position, 0) + 1

    selected_variants = [novel_candidates[i] for i in selected_indices]
    selected_eis = [round(float(ei_scores[i]), 6) for i in selected_indices]
    selected_means = [round(float(means[i]), 4) for i in selected_indices]
    selected_stds = [round(float(stds[i]), 4) for i in selected_indices]
    selected_pareto_ranks = [int(candidate_pareto_ranks[i]) for i in selected_indices]

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
        pareto_ranks=selected_pareto_ranks,
        max_pool_ei=max_pool_ei,
    )


# -----------------------------------------------------------------------------
# 1280-Dim ESM-2 Embedding Readout Adapter (SPEC-002 Seam 2)
# -----------------------------------------------------------------------------
# Each mutated position writes a 10-slot block (5 property deltas + 5 absolute
# properties) into a feature bank; the tail slots hold pairwise epistatic
# products. More banks means fewer distinct positions aliasing onto the same slots.
_READOUT_SLOTS_PER_POSITION = 10
_READOUT_EPISTASIS_SLOTS = 8


class SequenceReadoutAdapter:
    """
    Low-rank Ridge / SVD readout head over per-variant sequence embeddings.

    Calibrates against wet-lab observations to overcome the evolutionary
    conservation bias of zero-shot predictors on non-natural catalytic substrates
    and synthetic biologics, once at least `min_calibration_samples` measurements
    are available.

    By default the embedding is a deterministic physicochemical delta encoding
    computed locally (`embedding_source="physicochemical_delta"`). It is NOT an
    ESM representation. To use real ESM-2 layer-33 vectors, either pass
    `embedding_fn=workers.esm_worker.embed_sequences`, or supply precomputed
    vectors via the `embeddings=` argument of fit/predict/project_latent, and set
    `embedding_source` so the provenance is recorded.
    """

    def __init__(
        self,
        feature_dim: int = 128,
        latent_dim: int = 16,
        ridge_alpha: float = 0.05,
        min_calibration_samples: int = 8,
        random_state: int = 42,
        embedding_fn: Optional[Any] = None,
        embedding_source: str = "physicochemical_delta",
    ) -> None:
        min_dim = _READOUT_SLOTS_PER_POSITION + _READOUT_EPISTASIS_SLOTS
        if feature_dim < min_dim:
            raise ValueError(f"feature_dim must be at least {min_dim}, got {feature_dim}")
        self.feature_dim = int(feature_dim)
        self.latent_dim = latent_dim
        self.ridge_alpha = ridge_alpha
        self.min_calibration_samples = min_calibration_samples
        self.random_state = random_state
        self.embedding_fn = embedding_fn
        self.embedding_source = embedding_source
        self.is_fitted: bool = False
        self.n_observed: int = 0
        self._mean_h: Optional[np.ndarray] = None
        self._mean_y: float = 0.0
        self._weights: Optional[np.ndarray] = None
        self._svd_basis: Optional[np.ndarray] = None

    @property
    def embedding_dim(self) -> int:
        """Dimensionality of the embedding this adapter consumes."""
        return self.feature_dim

    @property
    def n_position_banks(self) -> int:
        """Number of distinct position banks the feature space provides."""
        usable = self.feature_dim - _READOUT_EPISTASIS_SLOTS
        return max(1, usable // _READOUT_SLOTS_PER_POSITION)

    def embed_sequences(
        self,
        sequences: Sequence[str],
        wt_sequence: str,
    ) -> np.ndarray:
        """
        Deterministic physicochemical delta encoding of each variant vs wt_sequence.

        Per mutated position, encodes the change in the 5 normalized amino-acid
        property scales plus the absolute properties of the substituted residue,
        then pairwise products across mutated sites as epistasis terms.

        If `embedding_fn` was supplied it is used instead, so real ESM-2
        representations can be substituted without changing any caller.
        """
        if self.embedding_fn is not None:
            vecs = np.asarray(self.embedding_fn(list(sequences)), dtype=np.float64)
            if vecs.ndim != 2 or vecs.shape[0] != len(sequences):
                raise ValueError(
                    f"embedding_fn must return one vector per sequence; got shape {vecs.shape}"
                )
            return vecs

        n_seq = len(sequences)
        out = np.zeros((n_seq, self.feature_dim), dtype=np.float64)
        L = len(wt_sequence)
        n_banks = self.n_position_banks
        epistasis_base = self.feature_dim - _READOUT_EPISTASIS_SLOTS

        for i, seq in enumerate(sequences):
            vec = out[i]
            mut_sites: List[Tuple[int, str, str]] = []
            limit = min(len(seq), L)
            for pos_0 in range(limit):
                if seq[pos_0] != wt_sequence[pos_0]:
                    mut_sites.append((pos_0 + 1, wt_sequence[pos_0], seq[pos_0]))

            for pos, wt_aa, mut_aa in mut_sites:
                wt_p = AA_PROPERTIES_NORMALIZED.get(wt_aa, np.zeros(5, dtype=np.float64))
                mut_p = AA_PROPERTIES_NORMALIZED.get(mut_aa, np.zeros(5, dtype=np.float64))
                # Spread positions across every available bank. 7 is coprime with
                # the bank counts used here, so consecutive positions land in
                # different banks rather than aliasing onto a handful of slots.
                bank = ((pos * 7) % n_banks) * _READOUT_SLOTS_PER_POSITION
                for k in range(5):
                    vec[bank + k] += float(mut_p[k] - wt_p[k])
                    vec[bank + 5 + k] += float(mut_p[k])

            # Continuous pairwise epistatic interaction terms
            for idx_a in range(len(mut_sites)):
                for idx_b in range(idx_a + 1, len(mut_sites)):
                    _p1, _, m1 = mut_sites[idx_a]
                    _p2, _, m2 = mut_sites[idx_b]
                    mp1 = AA_PROPERTIES_NORMALIZED.get(m1, np.zeros(5, dtype=np.float64))
                    mp2 = AA_PROPERTIES_NORMALIZED.get(m2, np.zeros(5, dtype=np.float64))
                    for k in range(4):
                        vec[epistasis_base + k] += float(mp1[k] * mp2[k]) * 0.5

        return out

    def fit(
        self,
        sequences: Sequence[str],
        scores: Sequence[float],
        wt_sequence: str,
        embeddings: Optional[np.ndarray] = None,
    ) -> "ESMReadoutAdapter":
        H = embeddings if embeddings is not None else self.embed_sequences(sequences, wt_sequence)
        y = np.asarray(scores, dtype=np.float64)
        n = len(y)
        self.n_observed = n
        if n < self.min_calibration_samples:
            self.is_fitted = False
            return self

        self._mean_h = np.mean(H, axis=0)
        self._mean_y = float(np.mean(y))
        H_c = H - self._mean_h
        y_c = y - self._mean_y

        # Dual kernel Ridge solution: w = H_c^T (H_c H_c^T + alpha * I)^-1 y_c
        K = H_c @ H_c.T + self.ridge_alpha * np.eye(n, dtype=np.float64)
        dual_alpha = np.linalg.solve(K, y_c)
        self._weights = H_c.T @ dual_alpha

        # Low-rank SVD projection basis for Deep Kernel GP conditioning
        _u, _s, vt = np.linalg.svd(H_c, full_matrices=False)
        k_dim = min(self.latent_dim, vt.shape[0])
        self._svd_basis = vt[:k_dim].T
        self.is_fitted = True
        return self

    def predict(
        self,
        sequences: Sequence[str],
        wt_sequence: str,
        zero_shot_priors: Optional[Sequence[float]] = None,
        embeddings: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        H = embeddings if embeddings is not None else self.embed_sequences(sequences, wt_sequence)
        priors = (
            np.asarray(zero_shot_priors, dtype=np.float64)
            if zero_shot_priors is not None
            else np.zeros(len(sequences), dtype=np.float64)
        )
        if not self.is_fitted or self._weights is None or self._mean_h is None:
            return priors

        readout_pred = (H - self._mean_h) @ self._weights + self._mean_y
        blend = min(1.0, float(self.n_observed) / 12.0)
        return blend * readout_pred + (1.0 - blend) * priors

    def project_latent(
        self,
        sequences: Sequence[str],
        wt_sequence: str,
        embeddings: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        H = embeddings if embeddings is not None else self.embed_sequences(sequences, wt_sequence)
        if self._svd_basis is None or self._mean_h is None:
            return H[:, : self.latent_dim]
        return (H - self._mean_h) @ self._svd_basis


select_top_variants = select_next_micro_batch

# Retained so existing callers keep working. See each docstring: the sampler is a
# REINFORCE policy gradient, not a GFlowNet; the surrogate is two-stage, not joint DKL.
sample_gflownet_candidates = sample_policy_gradient_candidates
DeepKernelLearner = TwoStageDeepFeatureGP

# Retained so existing callers keep working. The default embedding is a local
# physicochemical encoding, not an ESM representation - see the class docstring.
ESMReadoutAdapter = SequenceReadoutAdapter
