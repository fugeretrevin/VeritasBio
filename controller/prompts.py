"""
controller/prompts.py â€” DOMAIN: Brains & Security (Pillar 2: Strategic Planner).

Purpose:
    System instructions and few-shot JSON templates for Meta Llama 3.1.
    The LLM consumes biological RAG context from UniProt and formulates a search
    boundary strategy targeting synergistic combinatorial residue clusters.

Constraints & Invariants (ARCHITECTURE.md Â§3.2, Â§7.4; INNOVATION_ROADMAP Phase 2):
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
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set

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


class BiophysicalHypothesis(BaseModel):
    """
    Mechanistic structural hypothesis formulated by Meta Llama 3.1.
    Directly relates 3D pocket geometry, cofactor proximity, and residue chemistry
    to the active-learning exploration boundaries.
    """

    pocket_name: str = Field(..., description="Short name of pocket/loop, e.g. 'substrate_entrance_loop_1'")
    target_residues: List[int] = Field(..., min_length=1, description="1-indexed residue positions targeted")
    mechanism: str = Field(..., description="Biophysical mechanism (e.g. 'Steric volume relaxation', 'pi-pi stacking', 'Charge balance')")
    structural_rationale: str = Field(..., description="3D contact rationale referencing distance to catalytic core or cofactor")
    suggested_chemistries: Dict[str, List[str]] = Field(
        default_factory=dict,
        description="Suggested amino acid classes per position (e.g. {'184': ['I', 'V'], '227': ['M', 'Q']})",
    )
    epistatic_coupling: Optional[str] = Field(
        default=None,
        description="Predicted cooperative or antagonistic interaction with adjacent loops or cofactor",
    )

    @field_validator("suggested_chemistries", mode="before")
    @classmethod
    def clean_suggested_chemistries(cls, v: Any) -> Dict[str, List[str]]:
        if not isinstance(v, dict):
            return {}
        cleaned: Dict[str, List[str]] = {}
        for k, val in v.items():
            str_key = str(k).strip()
            if isinstance(val, list):
                cleaned[str_key] = [str(x).strip().upper() for x in val if isinstance(x, str) and len(str(x).strip()) == 1]
            elif isinstance(val, str):
                cleaned[str_key] = [c.upper() for c in val if c.isalpha() and len(c) == 1]
        return cleaned


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
    mechanistic_hypotheses: List[BiophysicalHypothesis] = Field(
        default_factory=list,
        description="Detailed biophysical hypotheses connecting 3D structure and chemistry",
    )

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
   - Residues grouped in a cluster MUST be physically or allosterically coupled in 3D tertiary space (mean C-alpha distance <= 16 Angstroms, max distance <= 24 Angstroms).
   - Interacting residues located within the same active-site pocket, coupled loop manifolds such as two entrance loops that face each other across a substrate channel, substrate channels, or interacting paratope CDR loops are prime targets for combinatorial epistasis.
   - Do NOT group residues separated across opposite poles of the scaffold that share no physical coupling (> 25 to 30 Angstroms).
7. BURIED HYDROPHOBIC CORE STABILITY:
   - Deeply buried hydrophobic core residues maintain the tertiary folding energy of the scaffold.
   - Never mutate buried hydrophobic core residues to small cavity-creating residues (e.g. Trp/Phe/Leu/Ile to Ala/Gly) or charged residues without a compensatory packing partner, as this induces void collapse or severe desolvation penalties.
   - Concentrate mutagenesis on surface-accessible binding pockets, flexible loop regions, and pocket-lining secondary shells.
8. MECHANISTIC BIOPHYSICAL CHAIN-OF-THOUGHT (CoT):
   - You MUST supply a structured "mechanistic_hypotheses" array containing entries for the proposed target clusters.
   - Formulate an explicit 3D biophysical mechanism (e.g. steric volume modification, electrostatic salt-bridge stabilization, aromatic pi-stacking, or loop conformational flexibility).
   - Reference distance and spatial orientation relative to the catalytic core or bound cofactor (e.g. Lys87 PLP cofactor in TrpB or Tyr66-Gly67 chromophore in GFP).
   - Propose specific recommended amino acid chemistry classes (e.g. branched hydrophobics [I, V, L], small neutrals [A, G, S], polar amides [Q, N], or aromatics [Y, F, W]) in "suggested_chemistries".
   - Predict positive or negative epistatic coupling between interacting positions in "epistatic_coupling".
9. SURFACE SHELL & ALLOSTERIC COOPERATIVITY EXPANSION:
   - When core active-site pockets have converged across iterations or when optimizing folding cooperativity, propose a secondary cluster of high-flexibility, solvent-exposed surface residues (e.g. outer loop, beta-sheet turn, or outer scaffold residues) that modulate folding kinetics and global stability without mutating core catalytic residues.
   - Surface positions tolerate diverse substitutions without disrupting the primary hydrophobic core, serving as key stability-rescuing permissive mutations for multi-mutant designs.
10. COFACTOR & TRANSITION-STATE CLASH REASONING:
   - When designing near catalytic cofactors (PLP, heme, flavin, NAD) or bound ligand spheres:
     - Immediate Transition-State Shell (< 4.5 Angstroms from cofactor/intermediate): preserve critical catalytic geometry; avoid bulky aromatic substitutions (W, Y, F), rigid prolines, or inappropriate charges that clash with or displace the cofactor transition state.
     - Substrate Entrance Channel (4.5 to 12.0 Angstroms): prime target for steric decongestion mutations (e.g. F/W/Y -> I/V/A/L or small neutrals) that relieve substrate bottlenecks and accommodate non-natural bulky substrates while maintaining loop closure.
"""

# A deliberately SYNTHETIC schema illustration.
#
# This example teaches the response format and nothing else. It must never describe a
# protein the engine is benchmarked on, because a few-shot example is read as guidance:
# the previous version used avGFP and named Ser65Thr and Thr203Tyr - the two best-known
# enhancing mutations of that landscape - and a target-tailored branch did the same for
# TrpB, naming three of the four residues of its published global optimum. Any campaign
# on those targets was therefore handed most of its answer before it began.
#
# Rules for editing this block:
#   - Use the placeholder accession and a generic protein description.
#   - Use positions that correspond to no real benchmark library.
#   - Keep suggested_chemistries to a neutral physicochemical spread, never a set
#     chosen because those substitutions are known to work.
# The target's real annotated sites already reach the model through the "Verified
# Biological Context" section of the user prompt, which is derived from UniProt/PDB
# rather than hand-written here.
FEW_SHOT_EXAMPLE_JSON = {
    "target_protein": "ExampleEnzyme (illustrative only - not a real target)",
    "target_accession": "X00000",
    "hypothesis": "Schema illustration: group annotated active-site contact residues into spatially coupled clusters, withhold annotated catalytic machinery, and propose a neutral physicochemical spread at each targeted position.",
    "target_clusters": [
        {
            "cluster_id": "primary_contact_shell",
            "positions": [41, 42],
            "rationale": "Annotated first-shell contact residues lining the substrate channel.",
        },
        {
            "cluster_id": "opposing_loop",
            "positions": [118, 119],
            "rationale": "Residues on the opposing loop within van der Waals contact of the first shell.",
        },
    ],
    "forbidden_positions": [1, 55, 150],
    "max_mutation_depth": 3,
    "mechanistic_hypotheses": [
        {
            "pocket_name": "substrate_channel",
            "target_residues": [41, 42, 118, 119],
            "mechanism": "Steric and electrostatic tuning of the substrate channel under conserved fold constraints.",
            "structural_rationale": "The two clusters face each other across the annotated channel, so substitutions in one are expected to couple to the other. Catalytic positions are withheld rather than mutated.",
            "suggested_chemistries": {
                "41": ["A", "S", "D", "K", "L", "F"],
                "42": ["A", "S", "D", "K", "L", "F"],
                "118": ["A", "S", "D", "K", "L", "F"],
                "119": ["A", "S", "D", "K", "L", "F"],
            },
            "epistatic_coupling": "Treat spatially coupled clusters as potentially epistatic and reason about combinations rather than single points.",
        }
    ],
}


def generate_few_shot_example(protein_metadata: Optional[Any] = None) -> Dict[str, Any]:
    """
    Returns the synthetic schema illustration, independent of the target.

    `protein_metadata` is accepted and ignored, for call-site compatibility. It used to
    tailor the example to the target, which is the wrong shape for a few-shot example:
    the model reads it as guidance, so tailoring it to the protein under study
    pre-answers the question. The target's own annotated sites already reach the model
    through the prompt's "Verified Biological Context" section.
    """
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
            try:
                from pathlib import Path
                from workers.spatial_filter import compute_anm_hinge_residues
                pdb_candidate = Path(__file__).resolve().parent.parent / "workers" / "data" / f"{pdb_id}.pdb"
                if pdb_candidate.is_file():
                    hinges = compute_anm_hinge_residues(
                        pdb_path=pdb_candidate,
                        cutoff_distance=12.0,
                        top_k=8,
                        catalytic_invariants=set(cat_pos) if cat_pos else None,
                    )
                    if hinges:
                        directives.append(
                            f"- ALLOSTERIC HINGE POSITIONS (GNM/ANM Normal Modes): Distal mechanical pivot residues {hinges} modulate global loop opening/closing dynamics."
                        )
            except Exception:
                pass

        if directives:
            lines.append("## Target-Specific Invariant & Structural Constraints:")
            lines.extend(directives)
            lines.append("")

    if top_evidence:
        lines.append("## Empirical Feedback from Previous Rounds (Iterative Hypothesis Refinement):")
        for ev in top_evidence[:5]:
            lines.append(f"- Variant: {ev.get('variant_sequence', 'N/A')} | Score: {ev.get('computed_score', 0.0):.4f}")
        lines.append("- Refine your mechanistic hypotheses to build on validated beneficial substitutions.")
        if round_number >= 3:
            lines.append("- STRATEGIC DIRECTIVE: Core active-site mutations have begun to converge. Consider proposing a secondary surface cluster or outer-loop flexible shell to unlock allosteric folding cooperativity.")
        lines.append("")

    example_json = generate_few_shot_example(protein_metadata)
    lines.extend([
        "## Required JSON Output Schema:",
        json.dumps(example_json, indent=2),
        "",
        "Now, analyze the biological context above and formulate the optimal SearchStrategy JSON. Return raw JSON only.",
    ])

    return "\n".join(lines)


def infer_assayable_positions(
    wt_sequence: str,
    evidence: Optional[Sequence[Dict[str, Any]]],
) -> List[int]:
    """
    Derives which positions the campaign's library actually varies, by diffing the
    measured variant sequences against wild type.

    This is observation, not configuration. A combinatorial library saturates a fixed
    set of sites, so every measured variant differs from WT only within that set: the
    scope is recoverable from the evidence the planner already holds, without anyone
    hardcoding position numbers. Returns [] when there is no evidence yet, in which
    case the caller has nothing to constrain on and should leave the search open.

    Accepts rows keyed `variant_sequence` (the ledger's evidence shape).
    """
    if not evidence or not wt_sequence:
        return []
    varying: Set[int] = set()
    for row in evidence:
        seq = (row or {}).get("variant_sequence") or ""
        if len(seq) != len(wt_sequence):
            # A length-changing construct tells us nothing about substitution scope.
            continue
        for index, (wt_aa, mut_aa) in enumerate(zip(wt_sequence, seq), start=1):
            if wt_aa != mut_aa:
                varying.add(index)
    return sorted(varying)


def restrict_strategy_to_assayable(
    strategy: SearchStrategy,
    assayable_positions: Optional[Iterable[int]],
) -> SearchStrategy:
    """
    Deterministically clips a SearchStrategy to positions the assay can measure.

    A planner reasoning about biophysics will legitimately propose distal allosteric
    residues, which is desirable science and unmeasurable on a library that saturates
    only a handful of sites: those proposals are discarded downstream as
    not_in_eligible_universe and the round's batch silently runs short. Prompt
    guidance alone cannot guarantee compliance, so the scope is also enforced here.

    Clips cluster positions, hypothesis target residues and suggested chemistries,
    drops anything left empty, and caps mutation depth at the number of measurable
    positions. The strategy is returned unchanged when no scope is supplied, which
    keeps the planner unconstrained by default. Raises ValueError when the scope and
    the strategy do not intersect at all, because silently substituting a different
    search space would misreport what was explored.
    """
    if not assayable_positions:
        return strategy

    allowed = {int(p) for p in assayable_positions}
    if not allowed:
        return strategy

    clipped_clusters = []
    for cluster in strategy.target_clusters:
        kept = [p for p in cluster.positions if p in allowed]
        if kept:
            clipped_clusters.append(cluster.model_copy(update={"positions": kept}))

    if not clipped_clusters:
        raise ValueError(
            "strategy proposes no assayable positions: "
            f"clusters targeted {sorted({p for c in strategy.target_clusters for p in c.positions})} "
            f"but only {sorted(allowed)} can be measured"
        )

    clipped_hypotheses = []
    for hyp in strategy.mechanistic_hypotheses or []:
        kept_residues = [p for p in hyp.target_residues if p in allowed]
        if not kept_residues:
            continue
        kept_chem = {
            pos: aas
            for pos, aas in (hyp.suggested_chemistries or {}).items()
            if str(pos).isdigit() and int(pos) in allowed
        }
        clipped_hypotheses.append(
            hyp.model_copy(update={"target_residues": kept_residues, "suggested_chemistries": kept_chem})
        )

    measurable = len({p for c in clipped_clusters for p in c.positions})
    return strategy.model_copy(update={
        "target_clusters": clipped_clusters,
        "mechanistic_hypotheses": clipped_hypotheses,
        "max_mutation_depth": max(1, min(strategy.max_mutation_depth, measurable)),
        # Keep only forbidden positions that are actually in scope, so the selector is
        # not handed a forbidden list referring to sites it cannot touch anyway.
        "forbidden_positions": [p for p in strategy.forbidden_positions if p in allowed],
    })


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


# -----------------------------------------------------------------------------
# Closed-Loop Mechanistic Residual Reflection (Gap 1 / Route C)
# -----------------------------------------------------------------------------
class ResidualOutlier(BaseModel):
    """An empirical outlier variant where GP surrogate prediction deviated from reality."""

    variant_notation: str = Field(..., description="Mutation notation (e.g. 'F184I:S228G')")
    predicted_score: float = Field(..., description="Prior surrogate predicted mean mu(x)")
    observed_score: float = Field(..., description="Empirical ground-truth score y(x)")
    residual: float = Field(..., description="Prediction residual e = y(x) - mu(x)")
    distance_angstroms: Optional[float] = Field(default=None, description="Max pairwise 3D distance between mutated residues")
    mutation_details: List[str] = Field(default_factory=list, description="Positional chemistries, e.g. ['F184I (aromatic->branched)', 'S228G (polar->flexible)']")


class ReflectionDiagnosis(BaseModel):
    """
    Structured biophysical diagnosis emitted by Meta Llama 3.1 upon analyzing
    active-learning prediction residuals between rounds.
    """

    target_protein: str = Field(..., description="Target protein name (e.g. TrpB TM9D8)")
    round_number: int = Field(..., description="Active learning round just completed")
    mechanistic_diagnosis: str = Field(
        ...,
        description="Biophysical and structural explanation for observed synergistic wins or epistatic crashes",
    )
    updated_suggested_chemistries: Dict[str, List[str]] = Field(
        default_factory=dict,
        description="Position-specific amino acid recommendations based on reflection (e.g. {'184': ['I', 'V']})",
    )
    forbidden_chemistries: Dict[str, List[str]] = Field(
        default_factory=dict,
        description="Position-specific amino acids to avoid based on observed crashes (e.g. {'228': ['K', 'R']})",
    )
    structural_cooperativity_insight: Optional[str] = Field(
        default=None,
        description="Specific 3D loop or pocket cooperativity insight guiding subsequent rounds",
    )

    @field_validator("updated_suggested_chemistries", "forbidden_chemistries", mode="before")
    @classmethod
    def clean_chemistries_map(cls, v: Any) -> Dict[str, List[str]]:
        if not isinstance(v, dict):
            return {}
        cleaned: Dict[str, List[str]] = {}
        for k, val in v.items():
            str_key = str(k).strip()
            if isinstance(val, list):
                cleaned[str_key] = [str(x).strip().upper() for x in val if isinstance(x, str) and len(str(x).strip()) == 1]
            elif isinstance(val, str):
                cleaned[str_key] = [c.upper() for c in val if c.isalpha() and len(c) == 1]
        return cleaned


REFLECTION_SYSTEM_PROMPT = """You are the Senior Biophysical Reasoning Engine for VeritasBio.
Your mission is Closed-Loop Residual Reflection: analyze active-learning surrogate prediction residuals (unpredicted synergistic wins and epistatic crashes) to guide subsequent rounds of protein engineering.

CRITICAL RULES:
1. You MUST respond ONLY with a raw JSON object matching the ReflectionDiagnosis schema.
2. Do NOT include markdown code fences (like ```json), commentary, or pleasantries.
3. You NEVER calculate math or evaluate equations directly.
4. Focus on 3D STRUCTURAL & CHEMICAL EXPLANATIONS:
   - Why did the positive outlier variants outperform the Gaussian Process model? (e.g. relieved steric hindrance, established new hydrogen bond, stabilized loop backbone).
   - Why did the negative outlier variants crash? (e.g. introduced buried charge desolvation, steric clash with cofactor, disrupted active site hinge).
5. Output actionable updated_suggested_chemistries and forbidden_chemistries per residue to steer candidate sampling in the next round.
"""


def get_reflection_system_prompt() -> str:
    """Returns the system prompt for closed-loop residual reflection."""
    return REFLECTION_SYSTEM_PROMPT.strip()


def build_reflection_prompt(
    target_protein: str,
    round_number: int,
    positive_outliers: List[ResidualOutlier],
    negative_outliers: List[ResidualOutlier],
    protein_summary: str,
    pdb_id: Optional[str] = None,
) -> str:
    """Constructs the prompt for Meta Llama 3.1 to reflect upon empirical prediction residuals."""
    lines = [
        f"Active Learning Iteration: Completed Round {round_number}",
        f"Target Protein: {target_protein}",
        f"Structural Coordinate PDB: {pdb_id or 'Local Structure'}",
        "",
        "## Verified Biological Context:",
        protein_summary,
        "",
        "## Empirical Gaussian Process Residual Outliers:",
        "The active-learning surrogate evaluated candidate variants. Below are the most significant prediction errors:",
        "",
        "### Top Positive Residuals (Unpredicted Synergistic Wins / Epistatic Goldmines):",
    ]

    if positive_outliers:
        for out in positive_outliers:
            dist_str = f" | 3D Max Distance: {out.distance_angstroms:.1f}A" if out.distance_angstroms else ""
            lines.append(
                f"- Variant: {out.variant_notation} | Predicted: {out.predicted_score:.4f} -> Observed: {out.observed_score:.4f} (Residual: +{out.residual:.4f}){dist_str}"
            )
            if out.mutation_details:
                lines.append(f"  Chemistries: {', '.join(out.mutation_details)}")
    else:
        lines.append("- (No strong positive outliers this round; predictions aligned well with observations)")

    lines.extend([
        "",
        "### Top Negative Residuals (Unpredicted Epistatic Crashes / Destabilizing Traps):",
    ])

    if negative_outliers:
        for out in negative_outliers:
            dist_str = f" | 3D Max Distance: {out.distance_angstroms:.1f}A" if out.distance_angstroms else ""
            lines.append(
                f"- Variant: {out.variant_notation} | Predicted: {out.predicted_score:.4f} -> Observed: {out.observed_score:.4f} (Residual: {out.residual:.4f}){dist_str}"
            )
            if out.mutation_details:
                lines.append(f"  Chemistries: {', '.join(out.mutation_details)}")
    else:
        lines.append("- (No severe negative crashes this round)")

    # Synthetic schema illustration, as in FEW_SHOT_EXAMPLE_JSON and for the same
    # reason. The previous version used TrpB's own positions and recommended "184": I
    # and "228": G - two residues of that landscape's published global optimum - in a
    # block the model is shown on every reflection call, whatever the target.
    # Placeholder positions here correspond to no real benchmark library.
    example_reflection = {
        "target_protein": "ExampleEnzyme (illustrative only - not a real target)",
        "round_number": round_number,
        "mechanistic_diagnosis": "Schema illustration: name the specific physical mechanism behind each residual. For a positive outlier, say which strain or clash the substitution relieved and which neighbouring position that enables. For a negative outlier, say which interaction it disrupted - charge placed near a cofactor, a cavity left unfilled, a buried polar atom left unsatisfied.",
        "updated_suggested_chemistries": {
            "41": ["A", "S", "L"],
            "118": ["G", "A", "V"],
        },
        "forbidden_chemistries": {
            "42": ["K", "R", "D", "E"],
        },
        "structural_cooperativity_insight": "Schema illustration: state the geometric relationship between the coupled clusters and what it implies for the next round's combinations.",
    }

    lines.extend([
        "",
        "## Required ReflectionDiagnosis JSON Output Schema:",
        json.dumps(example_reflection, indent=2),
        "",
        "Analyze the biophysical reasons behind these empirical residuals and output the ReflectionDiagnosis JSON only.",
    ])

    return "\n".join(lines)


def parse_reflection_response(raw_response: str) -> ReflectionDiagnosis:
    """Parses and validates ReflectionDiagnosis JSON from Meta Llama 3.1 output."""
    text = raw_response.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
        text = text.strip()

    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        logger.error("Failed to parse JSON from LLM reflection response: %s", exc)
        raise ValueError(f"LLM reflection output is not valid JSON: {exc}") from exc

    return ReflectionDiagnosis.model_validate(data)

