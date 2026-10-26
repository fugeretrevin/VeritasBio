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

import json
import logging
import os
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Union

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
from workers.esm_worker import score_variants
from workers.math_worker import SelectedBatch, select_next_micro_batch
from workers.mock_workers import MockLLMWorker
from workers.rag_worker import ProteinMetadata, fetch_protein_metadata

logger = logging.getLogger(__name__)


class OrchestratorConfig(BaseModel):
    """Configuration for an autonomous protein optimization run."""

    target_accession: str = Field(default="P42212", description="UniProt accession (avGFP)")
    total_budget_usd: float = Field(default=5.00, gt=0.0, description="Total budget in USD")
    budget_id: int = Field(default=ledger_db.DEFAULT_BUDGET_ID, ge=1)
    max_iterations: int = Field(default=10, ge=1, le=100, description="Hard cap on loop rounds (D8)")
    batch_size: int = Field(default=50, ge=1, le=500, description="Variants per micro-batch")
    target_model: str = Field(default="esm2_650m", description="esm2_650m, esm2_150m, or mock_esm")
    db_path: Optional[str] = Field(default=None, description="Path to ledger SQLite file")
    use_mock_llm: bool = Field(default=True, description="Default True to prevent accidental credit burning")
    enable_live_llm: bool = Field(default=False, description="Explicit opt-in required to invoke live Meta API")


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

    def _call_llama_api(self, system_prompt: str, user_prompt: str, max_retries: int = 3) -> str:
        """Invokes Meta Model API via HTTP POST in strict JSON mode with exponential backoff."""
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

        last_exc: Optional[Exception] = None
        for attempt in range(max_retries):
            try:
                req = urllib.request.Request(
                    url,
                    data=json.dumps(payload).encode("utf-8"),
                    headers={
                        "Content-Type": "application/json",
                        "Authorization": f"Bearer {api_key}",
                    },
                )
                with urllib.request.urlopen(req, timeout=60.0) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                    return data["choices"][0]["message"]["content"]
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
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
    ) -> SearchStrategy:
        """Step 2: Obtains SearchStrategy JSON from Meta Model API or fallback Mock."""
        live_opt_in = self.config.enable_live_llm or (os.getenv("ENABLE_LIVE_LLM", "").lower() in ("true", "1"))
        api_key = os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")

        # Strict credit protection: Only call live API if explicitly opted in and mock is not forced
        if live_opt_in and not self.config.use_mock_llm and api_key:
            try:
                model_name = "muse-spark-1.3" if api_key.startswith("LLM_") else "Llama 3.1/3.3"
                logger.info("Querying Meta Model API (%s) for Round %d strategy...", model_name, round_number)
                sys_prompt = get_strategy_system_prompt()
                user_prompt = build_strategy_user_prompt(
                    protein_summary=metadata.to_llm_summary(),
                    round_number=round_number,
                    top_evidence=top_evidence,
                )
                raw_json = self._call_llama_api(sys_prompt, user_prompt)
                return parse_strategy_response(raw_json, sequence_length=metadata.length)
            except Exception as exc:
                logger.warning(
                    "Meta Model API call failed (%s). Falling back to MockLLMWorker for Round %d.",
                    exc, round_number,
                )

        # Default fallback: Mock worker (zero credit burn)
        return self.mock_llm.generate_strategy(metadata, round_number=round_number, top_evidence=top_evidence)

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
        wt_sequence = metadata.sequence

        yield OrchestratorEvent(
            iteration=0,
            step="init",
            message=f"Target: {metadata.protein_name} ({len(wt_sequence)} AAs). Chromophore at residues {metadata.chromophore_residues}.",
            remaining_budget_usd=rem_usd,
            cumulative_spent_usd=spent_usd,
        )

        round_num = 1
        stop_reason = "max_iterations"

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

            strategy = self._generate_strategy(metadata, round_num, top_evidence_sorted)

            # -----------------------------------------------------------------
            # Step 3: Quant Selection (Gaussian Process Expected Improvement)
            # -----------------------------------------------------------------
            yield OrchestratorEvent(
                iteration=round_num,
                step="selection",
                message=f"Round {round_num}: Fitting GP model & evaluating Expected Improvement over clusters...",
                remaining_budget_usd=rem_usd,
                cumulative_spent_usd=spent_usd,
                current_job_id=job_id,
            )

            selected_batch = select_next_micro_batch(
                strategy=strategy,
                wt_sequence=wt_sequence,
                evidence=recent_evidence,
                batch_size=self.config.batch_size,
            )

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

            eval_request = BatchEvaluationRequest(
                job_id=job_id,
                target_model=self.config.target_model,
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
                scores = score_variants(
                    wt_sequence=wt_sequence,
                    variant_sequences=decision.mutated_sequences or [],
                    target_model=self.config.target_model,
                    use_mock_fallback=True,
                )
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
            results_pairs = list(zip(decision.mutated_sequences or [], scores))
            ledger_db.record_evidence(job_id, results_pairs, db_path=self.db_path)
            ledger_db.mark_job_completed(job_id, db_path=self.db_path)

            all_evidence = ledger_db.get_evidence(db_path=self.db_path)
            sorted_all = sorted(all_evidence, key=lambda e: float(e["computed_score"]), reverse=True)

            yield OrchestratorEvent(
                iteration=round_num,
                step="update",
                message=f"Round {round_num} COMPLETED: Recorded {len(results_pairs)} scores. Top fitness to date: {float(sorted_all[0]['computed_score']):.4f}",
                remaining_budget_usd=rem_after_spend,
                cumulative_spent_usd=spent_after_spend,
                current_job_id=job_id,
                top_variants=sorted_all[:5],
            )

            round_num += 1

        # Run Completed
        final_rem = ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
        final_spent = self.config.total_budget_usd - final_rem
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

        final_rem = ledger_db.get_remaining_budget(self.config.budget_id, self.db_path)
        final_spent = self.config.total_budget_usd - final_rem
        all_ev = ledger_db.get_evidence(db_path=self.db_path)
        all_jobs = ledger_db.list_jobs(db_path=self.db_path)
        best_ev = max(all_ev, key=lambda e: float(e["computed_score"])) if all_ev else None

        return OrchestratorSummary(
            target_protein="avGFP",
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
