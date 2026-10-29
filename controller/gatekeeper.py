"""
controller/gatekeeper.py — DOMAIN: Brains & Security (Pillar 4: Deterministic Gatekeeper).

Purpose:
    Intercepts variant candidate batches proposed by the Quant Engine / Orchestrator,
    validates them against rigorous biophysical Pydantic schemas, calculates deterministic
    USD execution costs under a two-state (cold vs. warm container) model, and delegates
    atomic fund authorization to `ledger.database.authorize_spend()`.

Constraints & Invariants (ARCHITECTURE.md §3, §7; INNOVATION_ROADMAP.md Phase 1):
    1. Fail-Closed Only: Any parsing error, biophysical invalidity (e.g. mutation mismatch
       against WT), unknown tool/model, or non-positive cost results in immediate rejection.
    2. Never Bypass the Ledger: Only `authorize_spend()` in `ledger/database.py` moves money
       or creates APPROVED records.
    3. Deterministic: No LLMs or heuristics in validation or cost calculation.
    4. Two-State Cloud Economics: Amortizes Modal container boot time when warm within
       the `scaledown_window`.
"""

from __future__ import annotations

import logging
import math
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

from pydantic import BaseModel, Field, field_validator, model_validator

from ledger import database as ledger_db

logger = logging.getLogger(__name__)

# Standard 20 canonical amino acids
CANONICAL_AMINO_ACIDS: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# Permitted foundation models for scoring (fail-closed against unknown endpoints)
ALLOWED_TARGET_MODELS: Set[str] = {
    "esm2_t33_650M_UR50D",
    "esm2_650m",
    "esm2_t30_150M_UR50D",
    "esm2_150m",
    "esm_cascade",
    "mock_esm",
    "mock_cascade",
    "rfdiffusion_proteinmpnn_a10g",
    "mock_diffusion",
}

# Mutation regex: e.g. S65T or A123V
_SINGLE_MUTATION_REGEX = re.compile(r"^([A-Z])(\d+)([A-Z])$")


# -----------------------------------------------------------------------------
# Biophysical Domain Models & Parsers
# -----------------------------------------------------------------------------
class SingleMutation(BaseModel):
    """Represents a single-point amino acid substitution, e.g. S65T."""

    wt_aa: str = Field(..., min_length=1, max_length=1, description="Wild-type 1-letter amino acid")
    position: int = Field(..., ge=1, description="1-indexed position in protein sequence")
    mut_aa: str = Field(..., min_length=1, max_length=1, description="Mutated 1-letter amino acid")

    @field_validator("wt_aa", "mut_aa")
    @classmethod
    def validate_canonical_aa(cls, v: str) -> str:
        v_upper = v.upper()
        if v_upper not in CANONICAL_AMINO_ACIDS:
            raise ValueError(f"'{v}' is not one of the 20 canonical amino acids")
        return v_upper

    @model_validator(mode="after")
    def validate_substitution(self) -> SingleMutation:
        if self.wt_aa == self.mut_aa:
            raise ValueError(
                f"Wild-type and mutant amino acids must differ at position {self.position} ({self.wt_aa} == {self.mut_aa})"
            )
        return self

    @classmethod
    def from_string(cls, mut_str: str) -> SingleMutation:
        """Parse notation like 'S65T'."""
        match = _SINGLE_MUTATION_REGEX.match(mut_str.strip())
        if not match:
            raise ValueError(f"Malformed mutation string '{mut_str}'. Expected format like 'S65T'.")
        wt, pos_str, mut = match.groups()
        return cls(wt_aa=wt, position=int(pos_str), mut_aa=mut)

    def to_string(self) -> str:
        return f"{self.wt_aa}{self.position}{self.mut_aa}"


class Variant(BaseModel):
    """
    Represents a protein sequence variant, potentially with multiple mutations
    (e.g. 'S65T' or 'S65T:Y66H:V68A').
    """

    notation: str = Field(..., description="Mutation notation string, e.g. 'S65T:Y66H'")
    mutations: List[SingleMutation] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def parse_notation(cls, data: Any) -> Any:
        if isinstance(data, str):
            data = {"notation": data}
        if isinstance(data, dict) and "notation" in data and not data.get("mutations"):
            raw = str(data["notation"]).strip()
            if not raw:
                raise ValueError("Variant notation cannot be empty")
            if raw == "WT":
                data["mutations"] = []
                return data
            parts = raw.split(":")
            parsed_muts = [SingleMutation.from_string(p) for p in parts]
            # Verify no colliding positions in the same multi-point variant
            positions = [m.position for m in parsed_muts]
            if len(positions) != len(set(positions)):
                raise ValueError(f"Colliding mutation positions in variant '{raw}': {positions}")
            data["mutations"] = parsed_muts
        return data

    @classmethod
    def from_sequence(cls, variant_sequence: str, wt_sequence: str) -> Variant:
        """Derive a Variant from a mutated sequence and reference wild-type sequence."""
        if len(variant_sequence) != len(wt_sequence):
            raise ValueError(f"Sequence length mismatch: {len(variant_sequence)} vs {len(wt_sequence)}")
        muts: List[SingleMutation] = []
        for i, (wc, vc) in enumerate(zip(wt_sequence, variant_sequence), 1):
            if wc != vc:
                muts.append(SingleMutation(wt_aa=wc, position=i, mut_aa=vc))
        if not muts:
            return cls(notation="WT", mutations=[])
        notation = ":".join(m.to_string() for m in muts)
        return cls(notation=notation, mutations=muts)

    def validate_against_wt(self, wt_sequence: str) -> None:
        """
        Verify that:
        1. All positions exist in the wild-type sequence.
        2. The stated WT amino acid exactly matches the actual sequence residue.
        """
        if self.notation == "WT":
            return
        seq_len = len(wt_sequence)
        for m in self.mutations:
            if m.position > seq_len:
                raise ValueError(
                    f"Mutation {m.to_string()} specifies position {m.position}, "
                    f"which exceeds WT sequence length of {seq_len}."
                )
            actual_wt = wt_sequence[m.position - 1].upper()
            if actual_wt != m.wt_aa:
                raise ValueError(
                    f"Biophysical mismatch in mutation {m.to_string()}: "
                    f"WT sequence has '{actual_wt}' at position {m.position}, but mutation expected '{m.wt_aa}'."
                )

    def apply_to(self, wt_sequence: str) -> str:
        """Apply the variant mutations to the WT sequence and return the full mutated sequence."""
        if self.notation == "WT":
            return wt_sequence
        self.validate_against_wt(wt_sequence)
        chars = list(wt_sequence.upper())
        for m in self.mutations:
            chars[m.position - 1] = m.mut_aa
        return "".join(chars)


# -----------------------------------------------------------------------------
# Gatekeeper Configuration & Request / Decision Schemas
# -----------------------------------------------------------------------------
class GatekeeperConfig(BaseModel):
    """
    Economics parameters for compute cost estimation.
    Default rates match Modal A10G public pricing ($1.10/hour).
    """

    a10g_usd_per_sec: float = Field(
        default=0.000306, ge=0.000001, description="Modal A10G hourly rate ($1.10/3600s)"
    )
    gpu_sec_per_variant: float = Field(
        default=0.02, ge=0.0001, description="Average GPU seconds per variant for ESM forward pass"
    )
    gpu_sec_per_variant_150m: float = Field(
        default=0.005, ge=0.0001, description="Average GPU seconds per variant for ESM-2 150M Tier-1 forward pass"
    )
    cascade_top_k_ratio: float = Field(
        default=0.20, ge=0.01, le=1.0, description="Fraction of top candidates promoted from Tier 1 to Tier 2"
    )
    cold_start_sec: float = Field(
        default=20.0, ge=0.0, description="Estimated cold-boot & PyTorch model load time"
    )
    idle_timeout_sec: float = Field(
        default=45.0, ge=0.0, description="Modal scaledown window keep-alive duration"
    )
    mock_nominal_cost: float = Field(
        default=0.0001, ge=0.000001, description="Nominal non-zero USD cost for mock test jobs (D2)"
    )
    gpu_sec_per_diffusion_backbone: float = Field(
        default=8.0, ge=0.1, description="Average GPU seconds per de novo backbone trajectory"
    )
    mock_diffusion_nominal_cost: float = Field(
        default=0.0005, ge=0.000001, description="Nominal non-zero USD cost for mock diffusion jobs"
    )
    min_variants_per_batch: int = Field(default=1, ge=1)
    max_variants_per_batch: int = Field(default=500, ge=1)


class BatchEvaluationRequest(BaseModel):
    """Payload presented to the Gatekeeper before any execution occurs."""

    job_id: str = Field(..., min_length=1)
    target_model: str = Field(...)
    wildtype_sequence: str = Field(..., min_length=1)
    variants: List[Variant] = Field(..., min_length=1)
    expected_improvements: Optional[List[float]] = None
    budget_id: int = Field(default=ledger_db.DEFAULT_BUDGET_ID, ge=1)
    operator_uid: Optional[str] = Field(default=None, description="Authenticated Firebase UID")
    operator_email: Optional[str] = Field(default=None, description="Authenticated operator email")
    auth_token: Optional[str] = Field(default=None, description="Firebase ID token for FDA Part 11 authentication")

    @field_validator("wildtype_sequence")
    @classmethod
    def validate_wt_sequence(cls, v: str) -> str:
        clean = v.strip().upper()
        if not clean:
            raise ValueError("Wild-type sequence cannot be empty")
        for char in clean:
            if char not in CANONICAL_AMINO_ACIDS:
                raise ValueError(f"Invalid amino acid character '{char}' in WT sequence")
        return clean

    @field_validator("target_model")
    @classmethod
    def validate_target_model(cls, v: str) -> str:
        clean = v.strip()
        if clean not in ALLOWED_TARGET_MODELS:
            raise ValueError(f"Target model '{clean}' is not in allowed models: {sorted(ALLOWED_TARGET_MODELS)}")
        return clean

    @model_validator(mode="after")
    def validate_biophysical_consistency(self) -> BatchEvaluationRequest:
        # Cross-validate all variants against the WT sequence
        for var in self.variants:
            var.validate_against_wt(self.wildtype_sequence)
        # Check Expected Improvements list length matches variants if provided
        if self.expected_improvements is not None:
            if len(self.expected_improvements) != len(self.variants):
                raise ValueError(
                    f"Length mismatch: {len(self.expected_improvements)} expected improvements "
                    f"for {len(self.variants)} variants."
                )
            for ei in self.expected_improvements:
                if not (math.isfinite(ei) and ei >= 0):
                    raise ValueError(f"Expected Improvement must be non-negative and finite, got {ei}")
        return self


class DeNovoDiffusionRequest(BaseModel):
    """Payload presented to the Gatekeeper for de novo backbone diffusion & inverse folding."""

    job_id: str = Field(..., min_length=1)
    target_model: str = Field(default="mock_diffusion")
    num_designs: int = Field(default=5, ge=1, le=50)
    binder_length: int = Field(default=65, ge=20, le=250)
    target_pdb_id: Optional[str] = None
    target_sequence: Optional[str] = None
    hotspot_residues: Optional[List[int]] = None
    expected_improvements: Optional[List[float]] = None
    budget_id: int = Field(default=ledger_db.DEFAULT_BUDGET_ID, ge=1)
    operator_uid: Optional[str] = Field(default=None, description="Authenticated Firebase UID")
    operator_email: Optional[str] = Field(default=None, description="Authenticated operator email")
    auth_token: Optional[str] = Field(default=None, description="Firebase ID token for FDA Part 11 authentication")

    @field_validator("target_model")
    @classmethod
    def validate_target_model(cls, v: str) -> str:
        clean = v.strip()
        if clean not in ALLOWED_TARGET_MODELS:
            raise ValueError(f"Target model '{clean}' is not in allowed models: {sorted(ALLOWED_TARGET_MODELS)}")
        return clean


class GatekeeperDecision(BaseModel):
    """The immutable decision record returned by the Gatekeeper."""

    approved: bool
    job_id: str
    target_model: str
    variant_count: int
    estimated_cost: float
    is_warm: bool
    rejection_reason: Optional[str] = None
    roi_ei_per_dollar: Optional[float] = None
    mutated_sequences: Optional[List[str]] = None
    timestamp: float = Field(default_factory=time.time)


# -----------------------------------------------------------------------------
# Gatekeeper Implementation
# -----------------------------------------------------------------------------
class Gatekeeper:
    """
    Deterministic cloud economics gatekeeper.

    Responsibilities:
      - Validates candidate batches fail-closed.
      - Calculates exact USD cost under cold/warm container state.
      - Evaluates Expected Improvement ROI ($EI/$).
      - Executes atomic spend authorization against `ledger.database`.
    """

    def __init__(
        self,
        config: Optional[GatekeeperConfig] = None,
        db_path: Optional[Union[str, Path]] = None,
    ) -> None:
        self.config: GatekeeperConfig = config or GatekeeperConfig()
        self.db_path: Optional[Path] = Path(db_path) if db_path else None
        self._last_execution_time: Optional[float] = None

    @property
    def is_warm(self) -> bool:
        """Indicates whether the cloud GPU container is currently within its warm keep-alive window."""
        if self._last_execution_time is None:
            return False
        return (time.time() - self._last_execution_time) <= self.config.idle_timeout_sec

    def calculate_cost(
        self,
        variant_count: int,
        target_model: str,
        current_time: Optional[float] = None,
        force_warm: Optional[bool] = None,
    ) -> Tuple[float, bool]:
        """
        Determine execution cost in USD and whether container is warm.

        Formula (D3 & INNOVATION_ROADMAP Phase 1):
          Cost = (variant_count * gpu_sec_per_variant + boot_penalty) * a10g_usd_per_sec
        """
        if target_model == "mock_esm":
            # Nominal synthetic cost for local/mock mode (D2)
            return self.config.mock_nominal_cost, True
        if target_model == "mock_cascade":
            # Nominal synthetic cost for local/mock cascade (1.5x mock_nominal_cost)
            return round(self.config.mock_nominal_cost * 1.5, 6), True
        if target_model == "mock_diffusion":
            # Nominal synthetic cost for local/mock de novo diffusion mode
            return self.config.mock_diffusion_nominal_cost, True

        now = current_time if current_time is not None else time.time()

        if force_warm is not None:
            is_warm = force_warm
        elif self._last_execution_time is not None:
            is_warm = (now - self._last_execution_time) <= self.config.idle_timeout_sec
        else:
            is_warm = False

        boot_sec = 0.0 if is_warm else self.config.cold_start_sec

        if target_model == "esm_cascade":
            # Multi-fidelity cascade: Tier 1 (150M) runs on all N variants;
            # Tier 2 (650M) runs only on top_k = ceil(N * cascade_top_k_ratio).
            # Co-located in the same container, incurring at most one cold boot.
            tier1_variants = variant_count
            tier2_variants = max(1, math.ceil(variant_count * self.config.cascade_top_k_ratio))
            active_gpu_sec = (tier1_variants * self.config.gpu_sec_per_variant_150m) + (
                tier2_variants * self.config.gpu_sec_per_variant
            )
        elif target_model in ("esm2_150m", "esm2_t30_150m_ur50d"):
            active_gpu_sec = variant_count * self.config.gpu_sec_per_variant_150m
        elif target_model == "rfdiffusion_proteinmpnn_a10g":
            active_gpu_sec = variant_count * self.config.gpu_sec_per_diffusion_backbone
        else:
            active_gpu_sec = variant_count * self.config.gpu_sec_per_variant

        total_gpu_sec = active_gpu_sec + boot_sec
        cost_usd = total_gpu_sec * self.config.a10g_usd_per_sec

        # Ensure cost is strictly positive and finite
        if not (math.isfinite(cost_usd) and cost_usd > 0):
            cost_usd = self.config.mock_nominal_cost

        return round(cost_usd, 6), is_warm

    def calculate_roi(
        self,
        expected_improvements: Optional[Sequence[float]],
        cost_usd: float,
    ) -> Optional[float]:
        """Calculate Expected Improvement per Dollar: ROI = sum(EI) / Cost."""
        if expected_improvements and cost_usd > 0:
            total_ei = sum(expected_improvements)
            return round(total_ei / cost_usd, 4)
        return None

    def evaluate_and_authorize(
        self,
        payload: Union[BatchEvaluationRequest, DeNovoDiffusionRequest, Dict[str, Any], str],
        force_warm: Optional[bool] = None,
    ) -> GatekeeperDecision:
        """
        Main gatekeeper pipeline:
          1. Fail-closed parsing & biophysical validation.
          2. Batch size checks.
          3. Deterministic cost & ROI calculation.
          4. Atomic ledger spend authorization.
        """
        if isinstance(payload, DeNovoDiffusionRequest) or (
            isinstance(payload, dict) and "num_designs" in payload and "variants" not in payload
        ):
            return self.evaluate_and_authorize_diffusion(payload, force_warm=force_warm)
        job_id = "unknown"
        target_model = "unknown"
        variant_count = 0

        # Step 1: Fail-closed parse and schema validation
        try:
            if isinstance(payload, str):
                request = BatchEvaluationRequest.model_validate_json(payload)
            elif isinstance(payload, dict):
                request = BatchEvaluationRequest.model_validate(payload)
            elif isinstance(payload, BatchEvaluationRequest):
                request = payload
            else:
                raise ValueError(f"Unsupported payload type: {type(payload)}")

            job_id = request.job_id
            target_model = request.target_model
            variant_count = len(request.variants)

        except Exception as exc:
            logger.warning("Gatekeeper REJECTED payload: schema validation error: %s", exc)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=variant_count,
                estimated_cost=0.0,
                is_warm=False,
                rejection_reason=f"Validation error: {exc}",
            )

        # Step 2: Batch limits check
        if not (self.config.min_variants_per_batch <= variant_count <= self.config.max_variants_per_batch):
            reason = (
                f"Variant count {variant_count} outside allowed batch bounds "
                f"[{self.config.min_variants_per_batch}, {self.config.max_variants_per_batch}]"
            )
            logger.warning("Gatekeeper REJECTED job '%s': %s", job_id, reason)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=variant_count,
                estimated_cost=0.0,
                is_warm=False,
                rejection_reason=reason,
            )

        # Step 3: Cost & ROI calculation
        cost_usd, is_warm = self.calculate_cost(
            variant_count=variant_count,
            target_model=target_model,
            force_warm=force_warm,
        )
        roi = self.calculate_roi(request.expected_improvements, cost_usd)

        # Precompute mutated full sequences for approved execution
        mutated_seqs = [v.apply_to(request.wildtype_sequence) for v in request.variants]

        # Step 4: Atomic spend authorization against the SQLite ledger
        operator_uid = request.operator_uid
        operator_email = request.operator_email
        if request.auth_token:
            from controller.auth import verify_firebase_id_token
            try:
                user = verify_firebase_id_token(request.auth_token)
                operator_uid = user.uid
                operator_email = user.email
            except Exception as auth_err:
                logger.warning("Gatekeeper REJECTED job '%s': authentication failed: %s", job_id, auth_err)
                return GatekeeperDecision(
                    approved=False,
                    job_id=job_id,
                    target_model=target_model,
                    variant_count=variant_count,
                    estimated_cost=0.0,
                    is_warm=is_warm,
                    rejection_reason=f"Authentication failed: {auth_err}",
                )

        authorized = ledger_db.authorize_spend(
            job_id=job_id,
            target_model=target_model,
            variant_count=variant_count,
            estimated_cost=cost_usd,
            budget_id=request.budget_id,
            db_path=self.db_path,
            operator_uid=operator_uid,
            operator_email=operator_email,
        )

        if not authorized:
            reason = "Ledger authorization rejected (insufficient budget or ledger constraint)"
            logger.warning("Gatekeeper REJECTED job '%s': %s", job_id, reason)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=variant_count,
                estimated_cost=cost_usd,
                is_warm=is_warm,
                rejection_reason=reason,
                roi_ei_per_dollar=roi,
            )

        # Spend authorized: update warm container timestamp and return success
        self._last_execution_time = time.time()
        logger.info(
            "Gatekeeper APPROVED job '%s': %d variants on %s ($%.6f, warm=%s, ROI=%.2f)",
            job_id,
            variant_count,
            target_model,
            cost_usd,
            is_warm,
            roi if roi is not None else 0.0,
        )
        return GatekeeperDecision(
            approved=True,
            job_id=job_id,
            target_model=target_model,
            variant_count=variant_count,
            estimated_cost=cost_usd,
            is_warm=is_warm,
            roi_ei_per_dollar=roi,
            mutated_sequences=mutated_seqs,
        )

    def evaluate_and_authorize_diffusion(
        self,
        payload: Union[DeNovoDiffusionRequest, Dict[str, Any], str],
        force_warm: Optional[bool] = None,
    ) -> GatekeeperDecision:
        """
        Gatekeeper pipeline for de novo backbone diffusion & inverse folding jobs:
          1. Fail-closed parsing & structural parameter validation.
          2. Trajectory count bounds check.
          3. Deterministic compute cost calculation.
          4. Atomic ledger spend authorization against SQLite ledger.
        """
        job_id = "unknown"
        target_model = "unknown"
        num_designs = 0

        # Step 1: Parse and validate schema fail-closed
        try:
            if isinstance(payload, str):
                request = DeNovoDiffusionRequest.model_validate_json(payload)
            elif isinstance(payload, dict):
                request = DeNovoDiffusionRequest.model_validate(payload)
            elif isinstance(payload, DeNovoDiffusionRequest):
                request = payload
            else:
                raise ValueError(f"Unsupported payload type: {type(payload)}")

            job_id = request.job_id
            target_model = request.target_model
            num_designs = request.num_designs
        except Exception as exc:
            logger.warning("Gatekeeper REJECTED diffusion payload: validation error: %s", exc)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=num_designs,
                estimated_cost=0.0,
                is_warm=False,
                rejection_reason=f"Validation error: {exc}",
            )

        # Step 2: Batch limits check
        if not (self.config.min_variants_per_batch <= num_designs <= self.config.max_variants_per_batch):
            reason = (
                f"Design count {num_designs} outside allowed batch bounds "
                f"[{self.config.min_variants_per_batch}, {self.config.max_variants_per_batch}]"
            )
            logger.warning("Gatekeeper REJECTED diffusion job '%s': %s", job_id, reason)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=num_designs,
                estimated_cost=0.0,
                is_warm=False,
                rejection_reason=reason,
            )

        # Step 3: Cost calculation
        cost_usd, is_warm = self.calculate_cost(
            variant_count=num_designs,
            target_model=target_model,
            force_warm=force_warm,
        )
        roi = self.calculate_roi(request.expected_improvements, cost_usd)

        # Step 4: Atomic spend authorization against the SQLite ledger
        operator_uid = request.operator_uid
        operator_email = request.operator_email
        if request.auth_token:
            from controller.auth import verify_firebase_id_token
            try:
                user = verify_firebase_id_token(request.auth_token)
                operator_uid = user.uid
                operator_email = user.email
            except Exception as auth_err:
                logger.warning("Gatekeeper REJECTED diffusion job '%s': authentication failed: %s", job_id, auth_err)
                return GatekeeperDecision(
                    approved=False,
                    job_id=job_id,
                    target_model=target_model,
                    variant_count=num_designs,
                    estimated_cost=0.0,
                    is_warm=is_warm,
                    rejection_reason=f"Authentication failed: {auth_err}",
                )

        authorized = ledger_db.authorize_spend(
            job_id=job_id,
            target_model=target_model,
            variant_count=num_designs,
            estimated_cost=cost_usd,
            budget_id=request.budget_id,
            db_path=self.db_path,
            operator_uid=operator_uid,
            operator_email=operator_email,
        )

        if not authorized:
            reason = "Ledger authorization rejected (insufficient budget or ledger constraint)"
            logger.warning("Gatekeeper REJECTED diffusion job '%s': %s", job_id, reason)
            return GatekeeperDecision(
                approved=False,
                job_id=job_id,
                target_model=target_model,
                variant_count=num_designs,
                estimated_cost=cost_usd,
                is_warm=is_warm,
                rejection_reason=reason,
                roi_ei_per_dollar=roi,
            )

        self._last_execution_time = time.time()
        logger.info(
            "Gatekeeper APPROVED diffusion job '%s': %d backbones on %s ($%.6f, warm=%s)",
            job_id,
            num_designs,
            target_model,
            cost_usd,
            is_warm,
        )
        return GatekeeperDecision(
            approved=True,
            job_id=job_id,
            target_model=target_model,
            variant_count=num_designs,
            estimated_cost=cost_usd,
            is_warm=is_warm,
            roi_ei_per_dollar=roi,
        )
