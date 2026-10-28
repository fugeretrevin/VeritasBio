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
Your sole role is to read target-specific biological metadata (UniProt functional domains, catalytic active sites, binding pockets, 3D structural coordinates, and mutagenesis literature) and define search boundaries for protein engineering across any protein family.

CRITICAL RULES:
1. You MUST respond ONLY with a raw JSON object matching the requested schema.
2. Do NOT include markdown code fences (like ```json), commentary, greetings, or conversational text.
3. You NEVER calculate math, evaluate fitness equations, or pick individual amino acid mutations.
4. Focus on COMBINATORIAL EPISTASIS: identify clusters of 2 to 5 interacting residues whose simultaneous mutation can produce synergistic functional improvements (e.g. active-site entrance loops, substrate binding pockets, cofactor coordination shells, allosteric interfaces, or packing contact shells).
5. CATALYTIC & BIOCHEMICAL INVARIANT GUARD:
   - Any residues annotated in the input context as "INVARIANT Catalytic Sites", covalent cofactor attachment sites (e.g. PLP-binding lysine), metal-coordinating residues, or strict conformational hinges MUST NEVER BE MUTATED.
   - You MUST place all annotated invariant positions into "forbidden_positions".
   - Destructive non-conservative mutations to catalytic residues abolish enzyme turnover, ligand binding, or chromophore cyclization.
6. 3D SPATIAL PROXIMITY CONSTRAINT:
   - Residues grouped in a cluster MUST be physically adjacent in 3D tertiary space (mean C-alpha distance <= 11 to 13 Angstroms).
   - Do NOT group residues separated across opposite sides of the protein scaffold that do not contact in 3D space.
   - Prioritize clusters located within the same structural pocket, adjacent secondary structure elements (e.g. adjacent beta-strands or alpha-helices), or interacting loop manifolds.
7. BURIED HYDROPHOBIC CORE STABILITY:
   - Deeply buried hydrophobic core residues maintain the tertiary folding energy of the scaffold.
   - Never mutate buried hydrophobic core residues to small cavity-creating residues (e.g. Trp/Phe/Leu/Ile to Ala/Gly) or charged residues without a compensatory packing partner, as this induces void collapse or severe desolvation penalties.
   - Concentrate mutagenesis on surface-accessible binding pockets, flexible loop regions, and pocket-lining secondary shells.
"""

FEW_SHOT_EXAMPLE_JSON = {
    "target_protein": "avGFP",
    "target_accession": "P42212",
    "hypothesis": "Focusing on the chromophore contact shell (Phe64, Ser65, Val68, Thr203) to discover synergistic variants that enhance spectral quantum yield and thermal stability while locking invariant catalytic backbone residue Gly67.",
    "target_clusters": [
        {
            "cluster_id": "chromophore_core",
            "positions": [64, 65, 68],
            "rationale": "Direct chromophore cavity residues that dictate fluorescence excitation/emission wavelengths and chromophore maturation kinetics."
        },
        {
            "cluster_id": "hydrophobic_pocket",
            "positions": [203, 206, 222],
            "rationale": "Spatial interaction partners in contact with the internal chromophore; known to modulate green-to-yellow shifts (T203Y) and monomeric stability."
        }
    ],
    "forbidden_positions": [1, 2, 67, 238],
    "max_mutation_depth": 3
}


def generate_few_shot_example(protein_metadata: Optional[Any] = None) -> Dict[str, Any]:
    """Generates a target-tailored few-shot SearchStrategy JSON conforming strictly to schema."""
    if protein_metadata is None:
        return FEW_SHOT_EXAMPLE_JSON

    prot_name = getattr(protein_metadata, "protein_name", "TargetProtein")
    acc = getattr(protein_metadata, "accession", "P00000")
    cat_pos = getattr(protein_metadata, "catalytic_positions", [])
    bind_pos = getattr(protein_metadata, "binding_positions", [])
    length = getattr(protein_metadata, "length", 200)

    # TrpB target family
    if "trpb" in prot_name.lower() or "tryptophan synthase" in prot_name.lower() or acc in ("P50909", "TRPB_TM9D8"):
        return {
            "target_protein": prot_name,
            "target_accession": acc,
            "hypothesis": "Targeting substrate entrance loop residues (183, 184, 227, 228) to modulate catalytic turnover while preserving invariant PLP-binding lysine 83.",
            "target_clusters": [
                {
                    "cluster_id": "substrate_entrance_loop_1",
                    "positions": [183, 184],
                    "rationale": "Active-site entrance loop gating substrate access to the PLP cofactor pocket."
                },
                {
                    "cluster_id": "substrate_entrance_loop_2",
                    "positions": [227, 228],
                    "rationale": "Opposing loop residues directly contacting loop 1 to stabilize the closed catalytic conformation."
                }
            ],
            "forbidden_positions": sorted(list(set([83] + cat_pos))),
            "max_mutation_depth": 4
        }

    # General target with annotated binding pocket
    if bind_pos and len(bind_pos) >= 2:
        c1 = sorted(bind_pos[:min(3, len(bind_pos))])
        c2 = sorted(bind_pos[min(3, len(bind_pos)):min(6, len(bind_pos))])
        clusters = [
            {
                "cluster_id": "active_pocket_shell",
                "positions": c1,
                "rationale": "Substrate/cofactor binding pocket residues directly modulating active-site geometry."
            }
        ]
        if c2:
            clusters.append({
                "cluster_id": "pocket_interaction_network",
                "positions": c2,
                "rationale": "Secondary shell interacting residues modulating conformational flexibility."
            })
        forbidden = sorted(list(set(cat_pos + [1, length])))
        return {
            "target_protein": prot_name,
            "target_accession": acc,
            "hypothesis": f"Targeting binding pocket shell residues while strictly preserving catalytic machinery at {cat_pos}.",
            "target_clusters": clusters,
            "forbidden_positions": forbidden,
            "max_mutation_depth": 3
        }

    return FEW_SHOT_EXAMPLE_JSON


def get_strategy_system_prompt() -> str:
    """Returns the immutable system prompt for the Strategic Planner."""
    return STRATEGY_SYSTEM_PROMPT


def build_strategy_user_prompt(
    protein_summary: str,
    round_number: int = 1,
    top_evidence: Optional[List[Dict[str, Any]]] = None,
    protein_metadata: Optional[Any] = None,
) -> str:
    """
    Constructs the operational user prompt containing UniProt metadata,
    dynamically grounded invariants, structural constraints, and recent evidence
    from previous active-learning rounds.
    """
    lines = [
        f"Active Learning Iteration: Round {round_number}",
        "",
        "## Verified Biological Context (UniProt RAG):",
        protein_summary,
        "",
    ]

    # Dynamically inject target-specific directives if protein_metadata is available
    if protein_metadata is not None:
        cat_pos = getattr(protein_metadata, "catalytic_positions", [])
        bind_pos = getattr(protein_metadata, "binding_positions", [])
        pdb_id = getattr(protein_metadata, "pdb_id", None)
        chrom_pos = getattr(protein_metadata, "chromophore_residues", [])

        directives = []
        if cat_pos:
            directives.append(
                f"- MANDATORY CATALYTIC INVARIANTS: You MUST include positions {cat_pos} in 'forbidden_positions'. Mutating these abolishes activity."
            )
        if chrom_pos:
            directives.append(
                f"- FLUOROPHORE INVARIANT: Invariant catalytic residues {chrom_pos[-1:]} must remain in 'forbidden_positions'."
            )
        if bind_pos:
            directives.append(
                f"- RECOMMENDED TARGET CLUSTERS (Binding Pocket / Entrance Loops): Prime candidate positions for combinatorial synergy: {bind_pos}."
            )
        if pdb_id:
            directives.append(
                f"- 3D SPATIAL GROUNDING: Proposed cluster residue pairs must be cohesive in structure PDB {pdb_id} (distance <= 13 Angstroms)."
            )

        if directives:
            lines.append("## Target-Specific Invariant & Structural Constraints:")
            lines.extend(directives)
            lines.append("")

    if top_evidence:
        lines.append("## Top High-Performing Mutants Discovered in Previous Rounds (Evidence):")
        for ev in top_evidence[:5]:
            lines.append(f"- Variant: {ev.get('variant_sequence', 'N/A')} | Score: {ev.get('computed_score', 0.0):.4f}")
        lines.append("")

    example_json = generate_few_shot_example(protein_metadata)
    lines.extend([
        "## Required JSON Output Schema:",
        json.dumps(example_json, indent=2),
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
