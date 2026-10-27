---
trigger: always_on
description: VeritasBio project memory — confirmed scope decisions, deviations from ARCHITECTURE.md, invariants, and build status. Read before any change.
---

# VeritasBio — Project Context (agent memory)

> Source of truth: `ARCHITECTURE.md` (v3.1, synced with D1–D8 on 2026-10-26). This file
> records **decisions made with the user**. If the two ever conflict, the decision log
> below wins and `ARCHITECTURE.md` must be updated to match. Last updated: 2026-10-26.

## What we are building
A **deterministic cloud-economics engine** for computational biology (Meta Global AI
Developer Hackathon → future B2B SaaS). An LLM only *proposes* search boundaries; math
selects variants; a fail-closed Python Gatekeeper + SQLite ledger decides whether paid GPU
compute (Modal A10G + ESM) may run. **Not an open-ended AI agent.**

## Non-negotiable rules (ARCHITECTURE.md §3, §7)
1. **Fail-closed:** reject on parse error, unknown tool, missing/invalid cost, any exception.
2. **Never bypass the ledger:** no worker is invoked by the LLM or UI. Every paid action
   must follow `ledger.database.authorize_spend(...) -> True`.
3. **LLM:** Llama 3.1, strict JSON mode, no conversational text, **never does math or
   executes code**, never chooses which step/tool runs next.
4. **Orchestrator:** fixed, unidirectional 6-step loop (§4). No agentic tool-selection loops.
5. **UI:** Gradio must never block on Modal; poll SQLite via `yield`/background tasks.
6. **Separation of concerns:** RAG / Planner / Quant / Gatekeeper / Cloud Muscle stay in
   their own modules (`workers/`, `controller/`, `ledger/`).

## Decision log (confirmed with user)
| # | Topic | Decision |
|---|-------|----------|
| D1 | Money precision | SQLite stores **INTEGER micro-USD** (1 USD = 1,000,000). Column names unchanged from §6. Public Python API uses float USD. Costs round **up**, allocations round **down**. |
| D2 | $0.00 jobs | `authorize_spend()` **keeps rejecting cost ≤ 0**. Mock mode charges a small **nominal synthetic cost** so the full budget flow is exercised. |
| D3 | Cost model | Deterministic: `(est_gpu_seconds_per_variant × variants + cold_start_seconds) × A10G_usd_per_second`. All constants in config, none from the LLM. |
| D4 | Modal failure after approval | New status **`FAILED`** (APPROVED → FAILED). **No refund** (Modal may have billed). Loop continues. |
| D5 | Budget scope | **One run = one fresh ledger file**, `Budget.id = 1`. No run_id column on Jobs. |
| D6 | Demo target | **avGFP** + **ESM-2 650M**; hidden answer key = **ProteinGym GFP DMS (Sarkisyan 2016)**; metric = Spearman (`scripts/evaluate.py`). |
| D7 | LLM host | **Meta Llama API** (`LLAMA_API_KEY`). |
| D8 | Variants & stopping | Point/multi-point mutations as strings, e.g. `S65T:Y66H`, applied to the WT sequence. Loop stops on budget exhaustion **or** a hard **max-iteration cap**. |
| D9 | Second benchmark | **TrpB** active-site landscape (Johnston et al. 2024), combinatorially complete 20⁴ at residues 183/184/227/228 of the engineered parent Tm9D8\*. Reference sequence is the **parent, not UniProt P50909** (they differ by 10 substitutions + a His-tag). Reported metrics are **held-out** $\rho$ and **enrichment vs. a random-screening baseline**, never in-sample fit alone. |

## Changes vs. the original v3.0 §6 (now reflected in ARCHITECTURE.md v3.1)
- Money columns are `INTEGER` micro-USD instead of `float` (D1).
- `Jobs.status` adds `FAILED` (D4).
- Added CHECK constraints and immutability triggers (no deletes, no evidence edits,
  monotonic spend, legal status transitions only, evidence only for APPROVED jobs).

## Ledger API (implemented: `ledger/database.py`)
- `init_db`, `create_budget`, `get_budget`, `get_remaining_budget`
- `authorize_spend(job_id, target_model, variant_count, estimated_cost) -> bool` (single spend gate, atomic, never raises)
- `mark_job_completed`, `mark_job_failed`, `get_job`, `list_jobs`
- `record_evidence(job_id, [(sequence, score), ...])` (atomic batch), `get_evidence`
- Job lifecycle: `APPROVED|REJECTED` on insert → `APPROVED → COMPLETED|FAILED`. `PENDING` reserved.

## Build status
- [x] Scaffold per §5 (purpose docstrings only). `.env.example` instead of committing `.env`.
- [x] `ledger/schema.sql`, `ledger/database.py` (46/46 smoke checks passing).
- [x] `INNOVATION_ROADMAP.md` (5 strategic gap solutions & 6-phase task breakdown).
- [x] `controller/gatekeeper.py` (Pydantic schemas + two-state cold/warm D3 cost model + ROI evaluator + 16 unit tests passing).
- [x] `workers/rag_worker.py`, `controller/prompts.py` (UniProt RAG + Llama 3.1 combinatorial prompts + 9 unit tests passing).
- [x] `workers/math_worker.py` (Physicochemical featurization + Gaussian Process active learning + 7 unit tests passing).
- [x] `workers/esm_worker.py`, `workers/mock_workers.py` (Modal A10G ESM-2 650M + high-fidelity mock + 8 unit tests passing).
- [x] `controller/orchestrator.py`, `ui/dashboard.py`, `main.py` (6-step active learning loop, async Gradio UI, 4 unit tests passing).
- [x] `scripts/evaluate.py`, `scripts/data/sarkisyan_gfp_benchmark.json` (ProteinGym DMS benchmark & cloud economics metrics, 6 unit tests passing).
- [x] Wet-lab translation & enterprise innovations: `workers/wetlab_worker.py`, `ui/structure_viewer.py`, `ui/landscape_plot.py` (codon optimization, Twist/IDT CSV, Opentrons OT-2 script, FDA 21 CFR Part 11 audit certificate, 3Dmol WebGL, 6 unit tests passing).
- [x] Iterative lead maturation & epistatic synergy: `workers/math_worker.py`, `controller/gatekeeper.py`, `workers/mock_workers.py` (lead recombination, sequence-to-variant parsing, round-progressive directed evolution).
- [x] Full test discovery passing: 60/60 unit tests across all test suites.
- [x] TrpB second benchmark: `scripts/download_trpb_benchmark.py`, `scripts/run_trpb_benchmark.py`,
      `tests/test_trpb_benchmark.py` (159,999 empirical variants, dataset-integrity assertions,
      held-out evaluation + random-screening baseline, 8 unit tests passing).

## Open questions (ask the user when relevant)
- Numeric constants for D3: A10G $/sec, GPU-seconds/variant, cold-start seconds.
- Value of the max-iteration cap (D8) and the mock nominal cost (D2).
- Llama 3.1 variant (8B / 70B / 405B) and exact Meta Llama API model ID.
- Python version and whether to pin `requirements.txt`.

