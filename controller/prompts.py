"""
controller/prompts.py — DOMAIN: Brains & Security (Pillar 2: Strategic Planner).

Purpose:
    System instructions and few-shot JSON templates for Meta Llama 3.1.
    The LLM consumes biological RAG context from UniProt and formulates a search
    boundary strategy targeting synergistic combinatorial residue clusters.

Constraints & Invariants (ARCHITECTURE.md §3.2, §7.4; INNOVATION_ROADMAP Phase 2):
    1. Strict JSON Mode: The LLM MUST output valid JSON conforming to SearchStrategy.
       No markdown conversational prose, commentary, or pleasantries.
    2. No Math / No Code: The LLM NEVER computes fitness scores, never runs optimization,
       and never decides the next system step.
    3. Combinatorial Hotspots (G2): Directs the model to focus on multi-residue epistatic
       clusters (e.g. residues 64-68 around the chromophore triad) rather than isolated single points.
    4. Fail-Closed Parsing: `parse_strategy_response()` strictly validates candidate positions
       against the target protein length.
"""

from __future__ import annotations

import json
import logging
import re
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Structured Strategy Schemas
# -----------------------------------------------------------------------------
class ResidueCluster(BaseModel):
    """A synergistic group of residues proposed for combinatorial exploration."""

    cluster_id: str = Field(..., description="Short identifier, e.g. 'chromophore_environment'")
    positions: List[int] = Field(..., min_length=1, description="1-indexed residue positions in protein")
    rationale: str = Field(..., description="Biological explanation of expected epistasis/synergy")

    @field_validator("positions")
    @classmethod
    def validate_positions(cls, v: List[int]) -> List[int]:
        if not v:
            raise ValueError("Residue cluster positions cannot be empty")
        for pos in v:
            if not isinstance(pos, int) or pos < 1:
                raise ValueError(f"Position {pos} must be a positive 1-indexed integer")
        # Ensure unique positions in cluster
        return sorted(list(set(v)))


class SearchStrategy(BaseModel):
    """
    Search boundary output proposed by Meta Llama 3.1.
    Passed to the Quant Engine (math_worker.py) to guide Bayesian Optimization.
    """

    target_protein: str = Field(..., description="Common protein name (e.g. avGFP)")
    target_accession: str = Field(..., description="UniProt accession (e.g. P42212)")
    hypothesis: str = Field(..., description="Scientific hypothesis explaining the chosen search boundaries")
    target_clusters: List[ResidueCluster] = Field(..., min_length=1, description="Residue clusters to sample combinations from")
    forbidden_positions: List[int] = Field(default_factory=list, description="Residues that should NOT be mutated")
    max_mutation_depth: int = Field(default=2, ge=1, le=4, description="Max simultaneous mutations (e.g. 2 for double mutants)")

    @model_validator(mode="after")
    def validate_no_overlap_forbidden(self) -> SearchStrategy:
        forbidden_set = set(self.forbidden_positions)
        for cluster in self.target_clusters:
            overlap = forbidden_set.intersection(set(cluster.positions))
            if overlap:
                raise ValueError(f"Cluster '{cluster.cluster_id}' contains forbidden positions: {sorted(overlap)}")
        return self


# -----------------------------------------------------------------------------
# Prompts
# -----------------------------------------------------------------------------
STRATEGY_SYSTEM_PROMPT = """You are the Strategic Biological Planner for VeritasBio, an autonomous cloud economics platform for computational biology.
Your sole role is to read biological metadata (UniProt functional domains, chromophore coordinates, mutagenesis literature) and define search boundaries for protein engineering.

CRITICAL RULES:
1. You MUST respond ONLY with a raw JSON object matching the requested schema.
2. Do NOT include markdown code fences (like ```json), commentary, greetings, or conversational text.
3. You NEVER calculate math, evaluate fitness equations, or pick individual amino acid mutations.
4. Focus on COMBINATORIAL EPISTASIS: identify clusters of 2 to 5 interacting residues whose simultaneous mutation can produce synergistic functional improvements (e.g. fluorophore core, packing contacts, dimer interface).
5. Specify forbidden positions if certain residues are essential for structural folding or cannot tolerate mutation.
"""

FEW_SHOT_EXAMPLE_JSON = {
    "target_protein": "avGFP",
    "target_accession": "P42212",
    "hypothesis": "Focusing on the chromophore triad (Ser65, Tyr66, Gly67) and adjacent steric contact residues (Val68, Thr203) to discover synergistic variants that enhance spectral quantum yield and thermal stability while avoiding disruption of the outer beta-barrel.",
    "target_clusters": [
        {
            "cluster_id": "chromophore_core",
            "positions": [64, 65, 66, 67, 68],
            "rationale": "Direct chromophore pocket residues that dictate fluorescence excitation/emission wavelengths and chromophore maturation kinetics."
        },
        {
            "cluster_id": "hydrophobic_pocket",
            "positions": [203, 206, 222],
            "rationale": "Spatial interaction partners in contact with the internal chromophore; known to modulate green-to-yellow shifts (T203Y) and monomeric stability."
        }
    ],
    "forbidden_positions": [1, 2, 238],
    "max_mutation_depth": 3
}


def get_strategy_system_prompt() -> str:
    """Returns the immutable system prompt for the Strategic Planner."""
    return STRATEGY_SYSTEM_PROMPT


def build_strategy_user_prompt(
    protein_summary: str,
    round_number: int = 1,
    top_evidence: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """
    Constructs the operational user prompt containing UniProt metadata
    and recent evidence from previous active-learning rounds.
    """
    lines = [
        f"Active Learning Iteration: Round {round_number}",
        "",
        "## Verified Biological Context (UniProt RAG):",
        protein_summary,
        "",
    ]

    if top_evidence:
        lines.append("## Top High-Performing Mutants Discovered in Previous Rounds (Evidence):")
        for ev in top_evidence[:5]:
            lines.append(f"- Variant: {ev.get('variant_sequence', 'N/A')} | Score: {ev.get('computed_score', 0.0):.4f}")
        lines.append("")

    lines.extend([
        "## Required JSON Output Schema:",
        json.dumps(FEW_SHOT_EXAMPLE_JSON, indent=2),
        "",
        "Now, analyze the biological context above and formulate the optimal SearchStrategy JSON. Return raw JSON only.",
    ])

    return "\n".join(lines)


def parse_strategy_response(raw_response: str, sequence_length: int) -> SearchStrategy:
    """
    Extracts and validates a SearchStrategy from LLM output.
    Fail-closed: raises ValueError if JSON is malformed or positions exceed sequence length.
    """
    text = raw_response.strip()

    # Strip markdown code blocks if the LLM hallucinated them despite instructions
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
        text = text.strip()

    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        logger.error("Failed to parse JSON from LLM strategy response: %s", exc)
        raise ValueError(f"LLM output is not valid JSON: {exc}") from exc

    strategy = SearchStrategy.model_validate(data)

    # Validate all positions against real protein length
    for cluster in strategy.target_clusters:
        for pos in cluster.positions:
            if pos > sequence_length:
                raise ValueError(
                    f"Position {pos} in cluster '{cluster.cluster_id}' exceeds sequence length ({sequence_length})"
                )

    for pos in strategy.forbidden_positions:
        if pos > sequence_length:
            raise ValueError(
                f"Forbidden position {pos} exceeds sequence length ({sequence_length})"
            )

    return strategy
