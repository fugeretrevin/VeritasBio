"""
ui/landscape_plot.py — Interactive Active-Learning Fitness Landscape Explorer (Plotly).

Purpose:
    Renders an interactive 2D latent physicochemical projection (PCA) of the variant space.
    Visualizes the Bayesian Optimization exploration/exploitation trade-off:
    - Explored variants colored by measured fitness.
    - Candidate pool colored by Expected Improvement (EI).
    - Marker sizes scaled by model uncertainty sigma(x).
    - Top hits highlighted with golden badges.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Set

import numpy as np
import plotly.graph_objects as go
from sklearn.decomposition import PCA

from controller.gatekeeper import Variant
from workers.math_worker import featurize_variant_on_positions

logger = logging.getLogger(__name__)


def generate_landscape_figure(
    evidence: List[Dict[str, Any]],
    candidate_variants: Optional[List[Variant]] = None,
    candidate_eis: Optional[List[float]] = None,
    candidate_means: Optional[List[float]] = None,
    candidate_stds: Optional[List[float]] = None,
    wt_sequence: Optional[str] = None,
) -> go.Figure:
    """
    Builds an interactive Plotly scatter plot showing the active learning search space.
    """
    fig = go.Figure()

    # Default WT if not provided
    if not wt_sequence:
        from workers.rag_worker import fetch_protein_metadata
        ctx = fetch_protein_metadata("P42212")
        wt_sequence = ctx.sequence

    # Collect all active mutation positions
    all_positions: Set[int] = {65, 66, 67, 163, 203}
    ev_variants: List[Tuple[Variant, float]] = []

    for ev in evidence:
        seq = ev.get("variant_sequence", "")
        if not seq or seq == wt_sequence:
            continue
        try:
            var = Variant.from_sequence(seq, wt_sequence)
            score = float(ev.get("computed_score", 0.0))
            ev_variants.append((var, score))
            for m in var.mutations:
                all_positions.add(m.position)
        except Exception:
            continue

    if candidate_variants:
        for var in candidate_variants[:250]:
            for m in var.mutations:
                all_positions.add(m.position)

    sorted_positions = sorted(list(all_positions))

    # 1. Featurize explored evidence variants
    ev_feats: List[np.ndarray] = []
    ev_scores: List[float] = []
    ev_labels: List[str] = []

    for var, score in ev_variants:
        try:
            feats = featurize_variant_on_positions(var, sorted_positions, wt_sequence)
            ev_feats.append(feats)
            ev_scores.append(score)
            mut_str = var.to_mutation_str() or "WT"
            ev_labels.append(f"<b>{mut_str}</b><br>Fitness: {score:.4f}<br>Status: Evaluated")
        except Exception:
            continue

    # 2. Featurize candidate pool if provided
    cand_feats: List[np.ndarray] = []
    cand_eis_clean: List[float] = []
    cand_labels: List[str] = []
    cand_sizes: List[float] = []

    if candidate_variants:
        for idx, var in enumerate(candidate_variants[:250]):
            try:
                feats = featurize_variant_on_positions(var, sorted_positions, wt_sequence)
                ei = candidate_eis[idx] if candidate_eis and idx < len(candidate_eis) else 0.05
                std = candidate_stds[idx] if candidate_stds and idx < len(candidate_stds) else 0.5
                mean = candidate_means[idx] if candidate_means and idx < len(candidate_means) else 0.0
                cand_feats.append(feats)
                cand_eis_clean.append(ei)
                cand_sizes.append(max(6, min(20, float(std) * 16)))
                cand_labels.append(
                    f"<b>{var.to_mutation_str()}</b><br>"
                    f"Expected Improvement: {ei:.6f}<br>"
                    f"Predicted Mean: {mean:.4f}<br>"
                    f"Uncertainty (σ): {std:.4f}<br>"
                    f"Status: Candidate"
                )
            except Exception:
                continue

    all_vectors = ev_feats + cand_feats
    if len(all_vectors) < 3:
        # Fallback informative graphic if cold start
        fig.add_annotation(
            text="Accumulating evidence... Run Active Learning round to project latent landscape.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="#94a3b8", family="Inter")
        )
        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor="#0b1329",
            plot_bgcolor="#0b1329",
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, visible=False),
            height=460,
            margin=dict(l=20, r=20, t=30, b=20),
        )
        return fig

    # Fit 2D PCA over 5D physicochemical space
    pca = PCA(n_components=2, random_state=42)
    coords_2d = pca.fit_transform(np.array(all_vectors))

    n_ev = len(ev_feats)
    ev_coords = coords_2d[:n_ev]
    cand_coords = coords_2d[n_ev:]

    # Trace 1: Candidate Pool (Unexplored Space)
    if len(cand_coords) > 0:
        fig.add_trace(go.Scatter(
            x=cand_coords[:, 0],
            y=cand_coords[:, 1],
            mode="markers",
            marker=dict(
                size=cand_sizes,
                color=cand_eis_clean,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(
                    title=dict(text="Expected<br>Improvement", side="top"),
                    thickness=12,
                    len=0.7,
                    x=1.02,
                ),
                opacity=0.65,
                line=dict(width=0.5, color="#38bdf8")
            ),
            hoverinfo="text",
            text=cand_labels,
            name="Candidate Pool (EI)",
        ))

    # Trace 2: Explored Evidence (Physical Scores)
    if len(ev_coords) > 0:
        fig.add_trace(go.Scatter(
            x=ev_coords[:, 0],
            y=ev_coords[:, 1],
            mode="markers",
            marker=dict(
                size=12,
                color=ev_scores,
                colorscale="Plasma",
                showscale=True,
                colorbar=dict(
                    title=dict(text="Measured<br>Fitness", side="top"),
                    thickness=12,
                    len=0.7,
                    x=1.15,
                ),
                symbol="circle",
                line=dict(width=2, color="#ffffff"),
            ),
            hoverinfo="text",
            text=ev_labels,
            name="Evaluated Evidence",
        ))

    # Trace 3: Highlight Top Hit
    if len(ev_scores) > 0:
        top_idx = int(np.argmax(ev_scores))
        fig.add_trace(go.Scatter(
            x=[ev_coords[top_idx, 0]],
            y=[ev_coords[top_idx, 1]],
            mode="markers+text",
            marker=dict(
                size=22,
                color="#fbbf24",
                symbol="star",
                line=dict(width=2, color="#000000"),
            ),
            text=["★ TOP HIT"],
            textposition="top center",
            hoverinfo="text",
            textfont=dict(color="#fbbf24", size=12, family="Inter"),
            name="Current Global Optimum",
        ))

    var_pc1 = round(pca.explained_variance_ratio_[0] * 100, 1)
    var_pc2 = round(pca.explained_variance_ratio_[1] * 100, 1)

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0b1329",
        plot_bgcolor="#0b1329",
        title=dict(
            text=f"Latent Physicochemical Landscape (PCA · {var_pc1}% + {var_pc2}% Variance)",
            font=dict(size=14, color="#f8fafc", family="Inter"),
            x=0.02,
            y=0.96,
        ),
        xaxis=dict(
            title="PC1 (Hydrophobicity & Electronic Descriptors)",
            gridcolor="#1e293b",
            zerolinecolor="#334155",
            color="#94a3b8",
        ),
        yaxis=dict(
            title="PC2 (Steric Bulk & Polarizability)",
            gridcolor="#1e293b",
            zerolinecolor="#334155",
            color="#94a3b8",
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0.0,
            font=dict(color="#cbd5e1", size=11),
        ),
        height=480,
        margin=dict(l=50, r=100, t=60, b=50),
    )

    return fig
