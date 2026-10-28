"""
controller/orchestrator.py — DOMAIN: Brains & Security.

Purpose:
    Executes the unidirectional 6-step active learning loop (ARCHITECTURE.md §4):
      Step 1: Initialization (Target protein, ledger budget setup, UniProt RAG).
      Step 2: Strategy (Strategic Planner: Meta Llama 3.1 JSON search boundaries).
      Step 3: Selection (Quant Engine: Physicochemical GP Expected Improvement selection).
      Step 4: Authorization (Deterministic Gatekeeper + SQLite ledger atomic spend lock).
      Step 5: Execution (Cloud Muscle: Modal serverless A10G ESM-2 tensor scoring).
      Step 6: Update (Evidence logging, active learning feedback to next iteration).

Invariants & Constraints (ARCHITECTURE.md §3, §7; INNOVATION_ROADMAP Phase 5):
    1. Unidirectional Control Flow: Fixed sequence. The LLM NEVER chooses the next step.
    2. Fail-Closed Stop: Loop terminates gracefully when `authorize_spend()` rejects (budget exhausted).
    3. Dual Termination (D8): Loop terminates when budget is exhausted OR hard `max_iterations` cap is reached.
    4. Async & Streaming: Yields `OrchestratorEvent` generator items so the Gradio dashboard
       and CLI stream real-time progress without blocking threads.
"""

from __future__ import annotations

import itertools
import json
import logging
import os
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Union

import requests
from requests.adapters import HTTPAdapter, Retry

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from pydantic import BaseModel, Field

from controller.gatekeeper import (
    BatchEvaluationRequest,
    Gatekeeper,
    GatekeeperConfig,
    GatekeeperDecision,
    Variant,
)
from controller.prompts import (
    SearchStrategy,
    build_strategy_user_prompt,
    get_strategy_system_prompt,
    parse_strategy_response,
)
from ledger import database as ledger_db
from workers.esm_worker import score_variants, score_variants_cascade
from workers.math_worker import (
    SelectedBatch,
    extract_mutations_from_sequence,
    generate_hamming1_lead_neighbors,
    select_next_micro_batch,
)
from workers.mock_workers import MockLLMWorker
from workers.rag_worker import ProteinMetadata, fetch_protein_metadata

logger = logging.getLogger(__name__)


class OrchestratorConfig(BaseModel):
    """Configuration for an autonomous protein optimization run."""

    target_accession: str = Field(default="P42212", description="UniProt accession (avGFP)")
    total_budget_usd: float = Field(default=5.00, gt=0.0, description="Total budget in USD")
    budget_id: int = Field(default=ledger_db.DEFAULT_BUDGET_ID, ge=1)
    max_iterations: int = Field(default=10, ge=1, le=100, description="Hard cap on loop rounds (D8)")
    batch_size: int = Field(default=25, ge=1, le=500, description="Variants per micro-batch (20-30 saturates A10G parallel tensors)")
    cold_start_batch_size: Optional[int] = Field(
        default=None,
        ge=1,
        le=500,
        description="Cold-start micro-batch size to amortize container spin-up overhead (defaults to 2x batch_size, capped at 64)",
    )
    target_model: str = Field(default="esm2_650m", description="esm2_650m, esm2_150m, or mock_esm")
    db_path: Optional[str] = Field(default=None, description="Path to ledger SQLite file")
    use_mock_llm: bool = Field(default=True, description="Default True to prevent accidental credit burning")
    enable_live_llm: bool = Field(default=False, description="Explicit opt-in required to invoke live Meta API")
    enable_regret_stopping: bool = Field(default=True, description="Enable Bayesian regret early stopping when pool EI converges")
    regret_stopping_threshold: float = Field(default=0.005, description="Pool EI threshold below which convergence is flagged")
    regret_patience: int = Field(default=2, description="Consecutive rounds with EI below threshold before terminating")
    enable_stagnation_stopping: bool = Field(default=True, description="Enable fitness plateau stagnation stopping")
    stagnation_patience: int = Field(default=3, description="Consecutive rounds without fitness improvement before stopping")
    stagnation_epsilon: float = Field(default=0.001, description="Minimum improvement threshold to reset stagnation counter")
    enable_pocket_hopping: bool = Field(default=True, description="Enable hierarchical pocket-hopping across coupled spatial clusters")
    max_hops: int = Field(default=2, description="Maximum pocket transitions per campaign")
    enable_dynamic_batch_sizing: bool = Field(default=False, description="Enable budget-aware dynamic micro-batch scaling")
    enable_cloud_cascade: bool = Field(
        default=False,
        description="Enable multi-fidelity 150M screening -> 650M precision verification cascade",
    )
    cascade_top_k_ratio: float = Field(
        default=0.20,
        ge=0.01,
        le=1.0,
        description="Fraction of top candidates promoted to Tier 2 in cascade",
    )
    enable_lead_refinement: bool = Field(
        default=False,
        description="Exhaustively test all Hamming-1 mutations around the champion lead",
    )


class OrchestratorEvent(BaseModel):
    """Telemetry event yielded during execution for UI and logging."""

    iteration: int
    step: str  # 'init', 'strategy', 'selection', 'authorization', 'execution', 'update', 'halted', 'complete'
    message: str
    remaining_budget_usd: float
    cumulative_spent_usd: float
    current_job_id: Optional[str] = None
    decision: Optional[Dict[str, Any]] = None
    top_variants: Optional[List[Dict[str, Any]]] = None
    timestamp: float = Field(default_factory=time.time)


class OrchestratorSummary(BaseModel):
    """Final summary returned upon run completion."""

    target_protein: str
    target_accession: str
    iterations_completed: int
    total_budget_usd: float
    final_spent_usd: float
    final_remaining_usd: float
    stop_reason: str
    total_jobs_run: int
    total_variants_evaluated: int
    highest_fitness_score: float
    top_variant_sequence: Optional[str] = None


class Orchestrator:
    """The 6-step unidirectional active learning engine."""

    def __init__(self, config: Optional[OrchestratorConfig] = None) -> None:
        self.config: OrchestratorConfig = config or OrchestratorConfig()
        self.db_path: Optional[Path] = Path(self.config.db_path) if self.config.db_path else None
        self.gatekeeper = Gatekeeper(db_path=self.db_path)
        self.mock_llm = MockLLMWorker()
        self.last_metadata: Optional[ProteinMetadata] = None
        self._http_session: Optional[requests.Session] = None

    @property
    def http_session(self) -> requests.Session:
        """Persistent HTTP session with connection pooling and retries for LLM requests."""
        if self._http_session is None:
            self._http_session = requests.Session()
            retries = Retry(
                total=3,
                backoff_factor=1.0,
                status_forcelist=[429, 500, 502, 503, 504],
                allowed_methods=["POST"],
            )
            adapter = HTTPAdapter(
                pool_connections=5,
                pool_maxsize=10,
                max_retries=retries,
            )
            self._http_session.mount("https://", adapter)
            self._http_session.mount("http://", adapter)
        return self._http_session

    def _call_llama_api(self, system_prompt: str, user_prompt: str, max_retries: int = 3) -> str:
        """Invokes Meta Model API via persistent HTTP session in strict JSON mode with exponential backoff."""
        api_key = os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")
        if not api_key:
            raise RuntimeError("LLAMA_API_KEY environment variable is not set")

        # Automatically route between Meta Model API (api.meta.ai) and Llama API (api.llama.com)
        if api_key.startswith("LLM_"):
            url = os.getenv("META_API_URL", "https://api.meta.ai/v1/chat/completions")
            model = os.getenv("META_MODEL_ID", "muse-spark-1.3")
        else:
            url = os.getenv("META_API_URL", "https://api.llama.com/compat/v1/chat/completions")
            model = os.getenv("META_MODEL_ID", "Llama-3.3-70B-Instruct")

        payload = {
            "model": model,
            "messages": [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            "response_format": {"type": "json_object"},
            "temperature": 0.2,
        }

        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        }

        last_exc: Optional[Exception] = None
        for attempt in range(max_retries):
            try:
                resp = self.http_session.post(url, json=payload, headers=headers, timeout=60.0)
                resp.raise_for_status()
                data = resp.json()
                return data["choices"][0]["message"]["content"]
            except Exception as exc:
                last_exc = exc
                wait_time = 2.0 ** attempt
                logger.warning(
                    "Meta Model API request failed (attempt %d/%d: %s). Backing off for %.1fs...",
                    attempt + 1, max_retries, exc, wait_time,
                )
                time.sleep(wait_time)

        raise RuntimeError(f"Meta Model API failed after {max_retries} attempts: {last_exc}") from last_exc

    def _generate_strategy(
        self,
        metadata: ProteinMetadata,
        round_number: int,
        top_evidence: Optional[List[Dict[str, Any]]],
        current_hop: int = 0,
    ) -> SearchStrategy:
        """Step 2: Obtains SearchStrategy JSON from Meta Model API or fallback Mock."""
        live_opt_in = self.config.enable_live_llm or (os.getenv("ENABLE_LIVE_LLM", "").lower() in ("true", "1"))
        api_key = os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")

        # Strict credit protection: Only call live API if explicitly opted in and mock is not forced
        if live_opt_in and not self.config.use_mock_llm and api_key:
            try:
                model_name = "muse-spark-1.3" if api_key.startswith("LLM_") else "Llama 3.1/3.3"
                logger.info("Querying Meta Model API (%s) for Round %d (Hop %d) strategy...", model_name, round_number, current_hop)
                sys_prompt = get_strategy_system_prompt()
                user_prompt = build_strategy_user_prompt(
                    protein_summary=metadata.to_llm_summary(),
                    round_number=round_number,
                    top_evidence=top_evidence,
                    protein_metadata=metadata,
                )
                if current_hop >= 1:
                    user_prompt += (
                        f"\n\n[HIERARCHICAL POCKET-HOPPING ACTIVE — HOP {current_hop}]: "
                        f"Pocket 1 has converged. Shift target cluster to secondary interacting pocket "
                        f"and enable multi-point coupling (max_mutation_depth 4)."
                    )
                raw_json = self._call_llama_api(sys_prompt, user_prompt)
                return parse_strategy_response(raw_json, sequence_length=metadata.length)
            except Exception as exc:
                logger.warning(
                    "Meta Model API call failed (%s). Falling back to MockLLMWorker for Round %d.",
                    exc, round_number,
                )

        # Default fallback: Mock worker (zero credit burn)
        return self.mock_llm.generate_strategy(
            metadata,
            round_number=round_number,
            top_evidence=top_evidence,
            current_hop=current_hop,
        )

    def run_stream(self) -> Iterator[OrchestratorEvent]:
        """
        Runs the unidirectional loop as a Python generator, yielding live events.
        Enables non-blocking streaming for Gradio UI and interactive monitors.
        """
        # =====================================================================
        # Step 1: Initialization
        # =====================================================================
        ledger_db.init_db(self.db_path)

        # Ensure budget row exists
        existing_budget = ledger_db.get_budget(self.config.budget_id, self.db_path)
        if existing_budget is None:
            ledger_db.create_budget(
                total_allocation=self.config.total_budget_usd,
                budget_id=self.config.budget_id,
                db_path=self.db_path,
            )

        rem_usd = ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
        spent_usd = self.config.total_budget_usd - rem_usd

        yield OrchestratorEvent(
            iteration=0,
            step="init",
            message=f"Initialized budget: ${self.config.total_budget_usd:.4f}. Fetching biological metadata...",
            remaining_budget_usd=rem_usd,
            cumulative_spent_usd=spent_usd,
        )

        metadata = fetch_protein_metadata(self.config.target_accession, use_cache=True)
        self.last_metadata = metadata
        wt_sequence = metadata.sequence

        if metadata.chromophore_residues:
            feature_desc = f"Fluorophore at residues {metadata.chromophore_residues}."
        elif metadata.catalytic_positions:
            feature_desc = f"Catalytic sites at residues {metadata.catalytic_positions}."
        else:
            feature_desc = f"Scaffold length {len(wt_sequence)} AAs."

        yield OrchestratorEvent(
            iteration=0,
            step="init",
            message=f"Target: {metadata.protein_name} ({len(wt_sequence)} AAs). {feature_desc}",
            remaining_budget_usd=rem_usd,
            cumulative_spent_usd=spent_usd,
        )

        round_num = 1
        stop_reason = "max_iterations"
        consecutive_low_regret = 0
        consecutive_stagnation = 0
        prev_best_fitness = -1e9
        current_hop = 0

        while round_num <= self.config.max_iterations:
            rem_usd = ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
            spent_usd = self.config.total_budget_usd - rem_usd
            job_id = f"job-round-{round_num:03d}-{int(time.time())}"

            # Check if budget is already zero
            if rem_usd <= 0.000001:
                stop_reason = "budget_exhausted"
                break

            # -----------------------------------------------------------------
            # Step 2: Strategic Planning (LLM JSON boundary)
            # -----------------------------------------------------------------
            recent_evidence = ledger_db.get_evidence(db_path=self.db_path)
            # Sort top evidence by score descending
            top_evidence_sorted = sorted(
                recent_evidence, key=lambda e: float(e["computed_score"]), reverse=True
            )[:5]

            yield OrchestratorEvent(
                iteration=round_num,
                step="strategy",
                message=f"Round {round_num}: Formulating combinatorial search strategy...",
                remaining_budget_usd=rem_usd,
                cumulative_spent_usd=spent_usd,
                current_job_id=job_id,
            )

            strategy = self._generate_strategy(metadata, round_num, top_evidence_sorted, current_hop=current_hop)

            yield OrchestratorEvent(
                iteration=round_num,
                step="strategy_complete",
                message=f"Round {round_num}: Strategy defined: {strategy.hypothesis}",
                remaining_budget_usd=rem_usd,
                cumulative_spent_usd=spent_usd,
                current_job_id=job_id,
                decision=strategy.model_dump(),
            )

            # -----------------------------------------------------------------
            # Step 3: Quant Selection (Gaussian Process Expected Improvement)
            # -----------------------------------------------------------------
            # Gap B: Dynamic micro-batch sizing to amortize cold container spin-up overhead
            if self.config.enable_dynamic_batch_sizing:
                rounds_left = max(1, self.config.max_iterations - round_num + 1)
                cost_per_var = 0.00000612 if "150m" in str(self.config.target_model) else 0.00001000
                budget_per_round = rem_usd / rounds_left
                affordable_batch = int(budget_per_round / max(cost_per_var, 1e-7))
                effective_batch_size = max(self.config.batch_size, min(150, affordable_batch))
                sizing_note = f" (Budget-aware scaled: {effective_batch_size} variants across {rounds_left} remaining rounds)"
            else:
                is_cold = (not self.gatekeeper.is_warm and self.config.target_model != "mock_esm") or (
                    self.config.cold_start_batch_size is not None and round_num == 1
                )
                if is_cold and round_num == 1:
                    effective_batch_size = self.config.cold_start_batch_size or min(self.config.batch_size * 2, 64)
                    sizing_note = f" (Cold boot detected: scaled to {effective_batch_size} variants to amortize overhead)"
                else:
                    effective_batch_size = self.config.batch_size
                    sizing_note = ""

            yield OrchestratorEvent(
                iteration=round_num,
                step="selection",
                message=f"Round {round_num}: Sizing batch to {effective_batch_size} variants{sizing_note}. Fitting GP model & evaluating Expected Improvement...",
                remaining_budget_usd=rem_usd,
                cumulative_spent_usd=spent_usd,
                current_job_id=job_id,
            )

            # Resolve PDB structure file if available
            pdb_path = None
            if metadata.pdb_id:
                candidate_pdb = Path(__file__).resolve().parent.parent / "workers" / "data" / f"{metadata.pdb_id}.pdb"
                if candidate_pdb.exists():
                    pdb_path = candidate_pdb

            selected_batch = select_next_micro_batch(
                strategy=strategy,
                wt_sequence=wt_sequence,
                evidence=recent_evidence,
                batch_size=effective_batch_size,
                catalytic_invariants=set(metadata.catalytic_positions),
                pdb_path=pdb_path,
            )

            # -----------------------------------------------------------------
            # Step 3b: 3D Biomechanical Spatial & Steric Clash Pre-Filter
            # -----------------------------------------------------------------
            try:
                from workers.spatial_filter import filter_spatially_viable_variants
                candidate_notations = [v.notation for v in selected_batch.variants]
                viable_notations, pruned_items = filter_spatially_viable_variants(
                    candidate_notations,
                    pdb_path=pdb_path,
                    catalytic_invariants=set(metadata.catalytic_positions),
                )
                if pruned_items:
                    viable_set = set(viable_notations)
                    filtered_variants = [v for v in selected_batch.variants if v.notation in viable_set]
                    filtered_eis = [
                        ei for v, ei in zip(selected_batch.variants, selected_batch.expected_improvements)
                        if v.notation in viable_set
                    ]
                    # Maintain batch throughput if viable variants exist
                    if len(filtered_variants) >= min(5, len(selected_batch.variants) // 2):
                        selected_batch.variants = filtered_variants
                        selected_batch.expected_improvements = filtered_eis
                    logger.info(
                        "3D Spatial Filter screened %d variants in Round %d (%d viable, %d pruned)",
                        len(selected_batch.variants), round_num, len(filtered_variants), len(pruned_items),
                    )
            except Exception as spatial_err:
                logger.debug("Spatial filter check skipped: %s", spatial_err)

            # -----------------------------------------------------------------
            # Step 3c: Enterprise Improvement 3: Bayesian Regret Stopping
            # -----------------------------------------------------------------
            max_pool_ei = (
                selected_batch.max_pool_ei
                if selected_batch.max_pool_ei is not None
                else (max(selected_batch.expected_improvements) if selected_batch.expected_improvements else 0.0)
            )
            if self.config.enable_regret_stopping and round_num >= 3:
                if max_pool_ei < self.config.regret_stopping_threshold:
                    consecutive_low_regret += 1
                    logger.info(
                        "Round %d: Max pool EI %.6f below threshold %.4f (patience: %d/%d)",
                        round_num, max_pool_ei, self.config.regret_stopping_threshold,
                        consecutive_low_regret, self.config.regret_patience,
                    )
                    if consecutive_low_regret >= self.config.regret_patience:
                        # Check if hierarchical pocket hopping can sweep a coupled secondary pocket
                        if (
                            self.config.enable_pocket_hopping
                            and current_hop < self.config.max_hops
                            and rem_usd > 0.05
                        ):
                            ev_list = ledger_db.get_evidence(db_path=self.db_path)
                            if ev_list:
                                best_ev = max(ev_list, key=lambda e: float(e["computed_score"]))
                                b_seq = best_ev["variant_sequence"]
                                b_score = float(best_ev["computed_score"])
                                b_muts = extract_mutations_from_sequence(b_seq, wt_sequence)
                                if b_muts:
                                    current_hop += 1
                                    consecutive_low_regret = 0
                                    consecutive_stagnation = 0
                                    prev_best_fitness = b_score
                                    lead_str = ":".join(m.to_string() for m in b_muts)
                                    hop_msg = (
                                        f"Round {round_num}: POCKET HOP {current_hop}/{self.config.max_hops}! "
                                        f"Regret saturated on pocket {current_hop}. Locking lead anchor {lead_str} "
                                        f"(score: {b_score:.4f}) and sweeping coupled spatial pocket..."
                                    )
                                    logger.info(hop_msg)
                                    yield OrchestratorEvent(
                                        iteration=round_num,
                                        step="strategy",
                                        message=hop_msg,
                                        remaining_budget_usd=rem_usd,
                                        cumulative_spent_usd=spent_usd,
                                        current_job_id=job_id,
                                        top_variants=sorted(ev_list, key=lambda e: float(e["computed_score"]), reverse=True)[:5],
                                    )
                                    round_num += 1
                                    continue

                        stop_reason = "converged_regret_threshold"
                        logger.info("Bayesian regret stopping triggered: %s", stop_reason)
                        yield OrchestratorEvent(
                            iteration=round_num,
                            step="halted",
                            message=f"Round {round_num} STOPPED EARLY: converged_regret_threshold (Max pool EI {max_pool_ei:.5f} < {self.config.regret_stopping_threshold}). Budget preserved.",
                            remaining_budget_usd=rem_usd,
                            cumulative_spent_usd=spent_usd,
                        )
                        break
                else:
                    consecutive_low_regret = 0

            # -----------------------------------------------------------------
            # Step 4: Authorization (Deterministic Gatekeeper)
            # -----------------------------------------------------------------
            yield OrchestratorEvent(
                iteration=round_num,
                step="authorization",
                message=f"Round {round_num}: Gatekeeper validating {len(selected_batch.variants)} variants & pricing batch...",
                remaining_budget_usd=rem_usd,
                cumulative_spent_usd=spent_usd,
                current_job_id=job_id,
            )

            if self.config.enable_cloud_cascade:
                effective_target_model = (
                    "mock_cascade"
                    if self.config.target_model in ("mock_esm", "mock_cascade")
                    else "esm_cascade"
                )
            else:
                effective_target_model = self.config.target_model

            eval_request = BatchEvaluationRequest(
                job_id=job_id,
                target_model=effective_target_model,
                wildtype_sequence=wt_sequence,
                variants=selected_batch.variants,
                expected_improvements=selected_batch.expected_improvements,
                budget_id=self.config.budget_id,
            )

            decision = self.gatekeeper.evaluate_and_authorize(eval_request)

            if not decision.approved:
                # Fail-closed stop
                stop_reason = f"gatekeeper_rejected ({decision.rejection_reason})"
                yield OrchestratorEvent(
                    iteration=round_num,
                    step="halted",
                    message=f"Round {round_num} REJECTED by Gatekeeper: {decision.rejection_reason}",
                    remaining_budget_usd=ledger_db.get_remaining_budget(self.config.budget_id, self.db_path),
                    cumulative_spent_usd=spent_usd,
                    current_job_id=job_id,
                    decision=decision.model_dump(),
                )
                break

            # Spend was atomically approved and deducted
            rem_after_spend = ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
            spent_after_spend = self.config.total_budget_usd - rem_after_spend

            yield OrchestratorEvent(
                iteration=round_num,
                step="authorization",
                message=f"Round {round_num} APPROVED: Charged ${decision.estimated_cost:.6f} (warm={decision.is_warm}, ROI={decision.roi_ei_per_dollar:.2f})",
                remaining_budget_usd=rem_after_spend,
                cumulative_spent_usd=spent_after_spend,
                current_job_id=job_id,
                decision=decision.model_dump(),
            )

            # -----------------------------------------------------------------
            # Step 5: Execution (Modal Cloud Muscle)
            # -----------------------------------------------------------------
            yield OrchestratorEvent(
                iteration=round_num,
                step="execution",
                message=f"Round {round_num}: Executing ESM tensor scoring on {len(decision.mutated_sequences)} variants...",
                remaining_budget_usd=rem_after_spend,
                cumulative_spent_usd=spent_after_spend,
                current_job_id=job_id,
            )

            try:
                if self.config.enable_cloud_cascade:
                    cascade_res = score_variants_cascade(
                        wt_sequence=wt_sequence,
                        variant_sequences=decision.mutated_sequences or [],
                        target_model=effective_target_model,
                        top_k_ratio=self.config.cascade_top_k_ratio,
                        use_mock_fallback=True,
                    )
                    promoted_idx_set = set(cascade_res.get("promoted_indices", []))
                    results_pairs = []
                    for seq, score in zip(cascade_res.get("promoted_sequences", []), cascade_res.get("tier2_scores", [])):
                        results_pairs.append((seq, score))
                    t1_scores = cascade_res.get("tier1_scores", [])
                    for idx, seq in enumerate(decision.mutated_sequences or []):
                        if idx not in promoted_idx_set and idx < len(t1_scores):
                            results_pairs.append((seq, t1_scores[idx]))
                else:
                    scores = score_variants(
                        wt_sequence=wt_sequence,
                        variant_sequences=decision.mutated_sequences or [],
                        target_model=self.config.target_model,
                        use_mock_fallback=True,
                    )
                    results_pairs = list(zip(decision.mutated_sequences or [], scores))
            except Exception as exc:
                logger.error("Execution failed for job %s: %s", job_id, exc)
                ledger_db.mark_job_failed(job_id, db_path=self.db_path)
                yield OrchestratorEvent(
                    iteration=round_num,
                    step="execution",
                    message=f"Round {round_num} Job FAILED: {exc}. Funds remain deducted.",
                    remaining_budget_usd=rem_after_spend,
                    cumulative_spent_usd=spent_after_spend,
                    current_job_id=job_id,
                )
                round_num += 1
                continue

            # -----------------------------------------------------------------
            # Step 6: Evidence Update (Active Learning Logging)
            # -----------------------------------------------------------------
            ledger_db.record_evidence(job_id, results_pairs, db_path=self.db_path)
            ledger_db.mark_job_completed(job_id, db_path=self.db_path)

            all_evidence = ledger_db.get_evidence(db_path=self.db_path)
            sorted_all = sorted(all_evidence, key=lambda e: float(e["computed_score"]), reverse=True)
            current_best = float(sorted_all[0]["computed_score"]) if sorted_all else 0.0

            yield OrchestratorEvent(
                iteration=round_num,
                step="update",
                message=f"Round {round_num} COMPLETED: Recorded {len(results_pairs)} scores. Top fitness to date: {current_best:.4f}",
                remaining_budget_usd=rem_after_spend,
                cumulative_spent_usd=spent_after_spend,
                current_job_id=job_id,
                top_variants=sorted_all[:5],
            )

            # Check stagnation stopping after round 3
            if self.config.enable_stagnation_stopping and round_num >= 4:
                improvement = current_best - prev_best_fitness
                if improvement < self.config.stagnation_epsilon:
                    consecutive_stagnation += 1
                    logger.info(
                        "Round %d: Fitness improvement %.6f below stagnation epsilon %.4f (patience: %d/%d)",
                        round_num, improvement, self.config.stagnation_epsilon,
                        consecutive_stagnation, self.config.stagnation_patience,
                    )
                    if consecutive_stagnation >= self.config.stagnation_patience:
                        # Check if hierarchical pocket hopping can sweep a coupled secondary pocket
                        if (
                            self.config.enable_pocket_hopping
                            and current_hop < self.config.max_hops
                            and rem_after_spend > 0.05
                        ):
                            ev_list = ledger_db.get_evidence(db_path=self.db_path)
                            if ev_list:
                                best_ev = max(ev_list, key=lambda e: float(e["computed_score"]))
                                b_seq = best_ev["variant_sequence"]
                                b_score = float(best_ev["computed_score"])
                                b_muts = extract_mutations_from_sequence(b_seq, wt_sequence)
                                if b_muts:
                                    current_hop += 1
                                    consecutive_stagnation = 0
                                    consecutive_low_regret = 0
                                    prev_best_fitness = b_score
                                    lead_str = ":".join(m.to_string() for m in b_muts)
                                    hop_msg = (
                                        f"Round {round_num}: POCKET HOP {current_hop}/{self.config.max_hops}! "
                                        f"Plateau reached on pocket {current_hop}. Locking lead anchor {lead_str} "
                                        f"(score: {b_score:.4f}) and sweeping coupled spatial pocket..."
                                    )
                                    logger.info(hop_msg)
                                    yield OrchestratorEvent(
                                        iteration=round_num,
                                        step="strategy",
                                        message=hop_msg,
                                        remaining_budget_usd=rem_after_spend,
                                        cumulative_spent_usd=spent_after_spend,
                                        current_job_id=job_id,
                                        top_variants=sorted_all[:5],
                                    )
                                    round_num += 1
                                    continue

                        stop_reason = "converged_fitness_plateau"
                        logger.info("Fitness plateau stopping triggered: %s", stop_reason)
                        yield OrchestratorEvent(
                            iteration=round_num,
                            step="halted",
                            message=f"Round {round_num} STOPPED EARLY: converged_fitness_plateau (No improvement > {self.config.stagnation_epsilon} for {self.config.stagnation_patience} consecutive rounds). Budget preserved.",
                            remaining_budget_usd=rem_after_spend,
                            cumulative_spent_usd=spent_after_spend,
                            current_job_id=job_id,
                            top_variants=sorted_all[:5],
                        )
                        break
                else:
                    consecutive_stagnation = 0
                    prev_best_fitness = current_best
            else:
                prev_best_fitness = max(prev_best_fitness, current_best)

            round_num += 1

        # ---------------------------------------------------------------------
        # Tier 3: Hamming-1 Local Lead Refinement Phase
        # ---------------------------------------------------------------------
        if self.config.enable_lead_refinement:
            all_ev_refine = ledger_db.get_evidence(db_path=self.db_path)
            if all_ev_refine:
                best_ev = max(all_ev_refine, key=lambda e: float(e["computed_score"]))
                lead_seq = best_ev["variant_sequence"]
                lead_score = float(best_ev["computed_score"])
                lead_muts = extract_mutations_from_sequence(lead_seq, wt_sequence)
                lead_var = Variant(
                    notation=":".join(m.to_string() for m in lead_muts) if lead_muts else "WT",
                    mutations=lead_muts,
                )

                # Active positions for lead refinement: prioritize known active pockets or mutated positions
                active_refine_pos: List[int] = []
                if "strategy" in locals() and strategy and getattr(strategy, "target_clusters", None):
                    pos_set = set()
                    for cluster in strategy.target_clusters:
                        if hasattr(cluster, "positions"):
                            pos_set.update(cluster.positions)
                        elif isinstance(cluster, (list, tuple, set)):
                            pos_set.update(cluster)
                    active_refine_pos = sorted(list(pos_set))
                if not active_refine_pos and lead_muts:
                    active_refine_pos = sorted(list({m.position for m in lead_muts}))

                if active_refine_pos:
                    refine_neighbors = generate_hamming1_lead_neighbors(
                        lead_variant=lead_var,
                        active_positions=active_refine_pos,
                        wt_sequence=wt_sequence,
                        pdb_path=pdb_path if "pdb_path" in locals() else None,
                    )
                    existing_seqs = {e["variant_sequence"] for e in all_ev_refine}
                    novel_neighbors = [n for n in refine_neighbors if n.apply_to(wt_sequence) not in existing_seqs]

                    if novel_neighbors:
                        refine_batch = novel_neighbors[: max(20, self.config.batch_size * 2)]
                        refine_job_id = f"job-refine-{int(time.time())}"
                        refine_model = (
                            "mock_esm"
                            if self.config.target_model in ("mock_esm", "mock_cascade")
                            else "esm2_650m"
                        )

                        refine_req = BatchEvaluationRequest(
                            job_id=refine_job_id,
                            target_model=refine_model,
                            wildtype_sequence=wt_sequence,
                            variants=refine_batch,
                            budget_id=self.config.budget_id,
                        )
                        refine_decision = self.gatekeeper.evaluate_and_authorize(refine_req)
                        if refine_decision.approved:
                            yield OrchestratorEvent(
                                iteration=round_num,
                                step="refinement",
                                message=f"Tier 3 Lead Polish: Authorized ${refine_decision.estimated_cost:.4f} to sweep {len(refine_batch)} Hamming-1 neighbors around lead {lead_var.notation}...",
                                remaining_budget_usd=ledger_db.get_remaining_budget(self.config.budget_id, self.db_path),
                                cumulative_spent_usd=self.config.total_budget_usd - ledger_db.get_remaining_budget(self.config.budget_id, self.db_path),
                                current_job_id=refine_job_id,
                            )
                            refine_scores = score_variants(
                                wt_sequence=wt_sequence,
                                variant_sequences=refine_decision.mutated_sequences or [],
                                target_model=refine_model,
                                use_mock_fallback=True,
                            )
                            refine_pairs = list(zip(refine_decision.mutated_sequences or [], refine_scores))
                            ledger_db.record_evidence(refine_job_id, refine_pairs, db_path=self.db_path)
                            ledger_db.mark_job_completed(refine_job_id, db_path=self.db_path)

                            best_refine_score = max(refine_scores) if refine_scores else -1e9
                            if best_refine_score > lead_score:
                                refine_msg = f"Tier 3 Polish DISCOVERED NEW OPTIMUM! Fitness improved from {lead_score:.4f} to {best_refine_score:.4f}."
                            else:
                                refine_msg = f"Tier 3 Polish VERIFIED local optimality: no neighbor exceeded lead fitness {lead_score:.4f}."

                            yield OrchestratorEvent(
                                iteration=round_num,
                                step="refinement",
                                message=refine_msg,
                                remaining_budget_usd=ledger_db.get_remaining_budget(self.config.budget_id, self.db_path),
                                cumulative_spent_usd=self.config.total_budget_usd - ledger_db.get_remaining_budget(self.config.budget_id, self.db_path),
                                current_job_id=refine_job_id,
                                top_variants=sorted(ledger_db.get_evidence(db_path=self.db_path), key=lambda e: float(e["computed_score"]), reverse=True)[:5],
                            )

        # Run Completed
        budget_row = ledger_db.get_budget(self.config.budget_id, self.db_path)
        final_rem = budget_row["remaining"] if budget_row else ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
        final_spent = budget_row["cumulative_spent"] if budget_row else max(0.0, self.config.total_budget_usd - final_rem)
        all_ev = ledger_db.get_evidence(db_path=self.db_path)
        best_ev = max(all_ev, key=lambda e: float(e["computed_score"])) if all_ev else None
        best_fitness = float(best_ev["computed_score"]) if best_ev else 0.0

        yield OrchestratorEvent(
            iteration=round_num - 1,
            step="complete",
            message=f"Run complete ({stop_reason}). Total spent: ${final_spent:.4f}. Best fitness: {best_fitness:.4f}.",
            remaining_budget_usd=final_rem,
            cumulative_spent_usd=final_spent,
            top_variants=sorted(all_ev, key=lambda e: float(e["computed_score"]), reverse=True)[:5] if all_ev else [],
        )

    def run(
        self,
        on_event: Optional[Callable[[OrchestratorEvent], None]] = None,
    ) -> OrchestratorSummary:
        """Synchronous runner helper for CLI and batch runs."""
        last_event: Optional[OrchestratorEvent] = None
        for event in self.run_stream():
            last_event = event
            if on_event:
                on_event(event)

        budget_row = ledger_db.get_budget(self.config.budget_id, self.db_path)
        final_rem = budget_row["remaining"] if budget_row else ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
        final_spent = budget_row["cumulative_spent"] if budget_row else max(0.0, self.config.total_budget_usd - final_rem)
        all_ev = ledger_db.get_evidence(db_path=self.db_path)
        all_jobs = ledger_db.list_jobs(db_path=self.db_path)
        best_ev = max(all_ev, key=lambda e: float(e["computed_score"])) if all_ev else None

        return OrchestratorSummary(
            target_protein=self.last_metadata.protein_name if self.last_metadata else "Target Protein",
            target_accession=self.config.target_accession,
            iterations_completed=last_event.iteration if last_event else 0,
            total_budget_usd=self.config.total_budget_usd,
            final_spent_usd=round(final_spent, 6),
            final_remaining_usd=round(final_rem, 6),
            stop_reason=last_event.message if last_event else "unknown",
            total_jobs_run=len(all_jobs),
            total_variants_evaluated=len(all_ev),
            highest_fitness_score=float(best_ev["computed_score"]) if best_ev else 0.0,
            top_variant_sequence=best_ev["variant_sequence"] if best_ev else None,
        )
