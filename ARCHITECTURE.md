# PanoptesBio: System Architecture & Execution Blueprint
**Version:** 3.1 (Frontier Architecture — ledger decisions D1–D8 incorporated, see §8)
**Objective:** Meta Global AI Developer Hackathon & Future B2B SaaS

## 1. Vision & Value Proposition
**The Pitch:** PanoptesBio is an enterprise-grade cloud economics platform for computational biology. It solves the "brute-force GPU waste" problem by treating LLMs not as physics engines, but as heavily restricted orchestration layers.
**The Value:** By combining biological RAG, Bayesian Optimization, and a strict deterministic financial Gatekeeper, PanoptesBio achieves frontier-level scientific discovery (target binding/fitness) while reducing a lab's cloud compute bill by 80%.

## 2. Problem & Solution
*   **The Problem:** LLMs hallucinate spatial thermodynamics. On the other hand, brute-force foundation model screening (like Anthropic's scorpion venom test generating 300k variants) requires unrestricted GPU compute that is economically impossible for most labs. 
*   **The Solution:** An autonomous "Chief Scientist" loop. The AI reads literature (RAG) to set search boundaries. A Gaussian Process handles the rigorous math to select the best micro-batch of variants. A Python Gatekeeper intercepts the batch, calculates the exact USD cost, checks an immutable ledger, and only spins up the remote GPU (Modal) if the budget allows. The system learns from the results and iteratively improves.

---

## 3. Core System Pillars (Separation of Concerns)
To guarantee enterprise-grade security and prevent runaway AI loops, the system is strictly separated into five domains:

1.  **The Biological Priors (RAG):** UniProt REST API. Provides verified biological constraints (e.g., active binding sites).
2.  **The Strategic Planner (LLM):** Meta Llama 3.1. Operates in strict JSON-mode. Reads the RAG context and formulates a search boundary. **Constraint: The LLM is never allowed to do math or directly execute code.**
3.  **The Quant Engine (Math):** `scikit-learn` Gaussian Process (GP) & Bayesian Optimization (BO). Takes the LLM's boundary, maps the fitness landscape, calculates *Expected Improvement*, and mathematically selects the next 50 variants.
4.  **The Deterministic Gatekeeper (Security):** Python + Pydantic + SQLite. Intercepts the selected variants, validates the schema, calculates the USD cost, checks the remaining budget, and explicitly approves or hard-rejects the batch.
5.  **The Cloud Muscle (Execution):** Modal serverless A10G GPUs. Boots up, downloads the ESM foundation model, scores the variants via batched PyTorch tensors, and spins down instantly.

---

## 4. The 6-Step Active Learning Data Flow
This is the unidirectional `while` loop that the `orchestrator.py` manages. It stops when **either** the budget can no longer cover the next batch **or** a hard max-iteration cap is reached (D8). The step order is fixed; the LLM never chooses what runs next.

1.  **Initialization:** User sets target protein (demo: avGFP) and a budget (e.g., $50). A fresh ledger file is created with `Budget.id = 1` (D5). The RAG engine fetches UniProt metadata.
2.  **Strategy:** Llama 3.1 (Meta Llama API, D7) reads the UniProt metadata and outputs a JSON search strategy.
3.  **Selection:** The Gaussian Process (BO) evaluates the strategy and selects a highly targeted micro-batch of 50 variants, expressed as mutation strings on the wild-type (e.g., `S65T:Y66H`, D8).
4.  **Authorization:** The Gatekeeper validates the batch, calculates the cost with the deterministic formula (D3), and calls `authorize_spend()`, which atomically checks the `Budget` table, deducts the cost, and logs the job as `APPROVED` — or logs it as `REJECTED`. (Fails closed on any error.)
5.  **Execution:** The batch is sent to Modal. ESM-2 650M scores the sequences and returns the biological fitness values. If Modal fails, the job is marked `FAILED`; the funds stay spent and the loop continues (D4).
6.  **Update (Active Learning):** The scores are logged in the SQLite `Evidence` table, the job is marked `COMPLETED`, and the scores are fed *back* into the Gaussian Process to update its mathematical map. Loop returns to Step 2.

---

## 5. Directory & File Scaffolding
```text
panoptesbio/
├── .env                  # Secrets (git-ignored): LLAMA_API_KEY, HF_TOKEN, MODAL_TOKEN_ID
├── .env.example          # Committed template for .env (no real values).
├── .agents/rules/        # AI assistant memory: project_context.md (always-on rule).
├── requirements.txt      # modal, gradio, pydantic, torch, scikit-learn, requests
├── main.py               # Entry point. Inits DB, starts UI, wires Orchestrator.
│
├── controller/           # DOMAIN: Brains & Security
│   ├── gatekeeper.py     # Pydantic schemas, cost calculation, budget blocking.
│   ├── orchestrator.py   # The While-loop managing the 6-step data flow.
│   └── prompts.py        # System instructions and few-shot JSON examples for LLM.
│
├── ledger/               # DOMAIN: Immutable Memory
│   ├── schema.sql        # DDL for Budget, Jobs, Evidence + CHECK constraints & immutability triggers.
│   └── database.py       # SQLite CRUD wrapper (init_db, create_budget, authorize_spend,
│                         #   mark_job_completed/failed, record_evidence, get_*/list_*).
│
├── workers/              # DOMAIN: Science & Math Tools
│   ├── rag_worker.py     # Hits rest.uniprot.org API for biological metadata.
│   ├── math_worker.py    # Gaussian Process / Expected Improvement selection.
│   ├── esm_worker.py     # Modal @function for remote GPU scoring.
│   └── mock_workers.py   # Offline stand-ins for local UI/flow testing (nominal synthetic cost, D2).
│
├── ui/                   # DOMAIN: Frontend
│   └── dashboard.py      # Gradio interface with async DB polling (no freezing).
│
└── scripts/              # DOMAIN: Hackathon Evaluation
    └── evaluate.py       # Spearman rank vs hidden answer key (ProteinGym GFP DMS, D6).
```

> `deduct_funds` / `log_job` from v3.0 are intentionally **not** public: they are folded into `authorize_spend()` so that money can only move together with an approval record.

---

## 6. Database Schema Specifications (SQLite)
The application relies on a local SQLite database with three core tables. **All money columns are `INTEGER` micro-USD** (1 USD = 1,000,000) for exact arithmetic (D1); the Python API in `database.py` accepts and returns float USD, rounding costs **up** and allocations **down**.
*   **Budget Table:** `id` (int, pk), `total_allocation` (int micro-USD), `cumulative_spent` (int micro-USD). Invariant: `cumulative_spent <= total_allocation`.
*   **Jobs Table:** `job_id` (str, pk), `target_model` (str), `variant_count` (int > 0), `estimated_cost` (int micro-USD), `status` (str: PENDING, APPROVED, REJECTED, COMPLETED, FAILED).
*   **Evidence Table:** `evidence_id` (pk), `job_id` (fk), `variant_sequence` (str), `computed_score` (float), `timestamp` (ISO-8601 UTC).

**Job lifecycle:** jobs are inserted as `APPROVED` or `REJECTED` by `authorize_spend()`; legal transitions are `PENDING → APPROVED | REJECTED` and `APPROVED → COMPLETED | FAILED`. `PENDING` is reserved for future use.

**Immutability (enforced by triggers in `schema.sql`):**
*   No rows may be deleted from any table.
*   `cumulative_spent` can only increase (no refunds, D4).
*   Job identity and cost fields can never change; only `status`, along legal transitions.
*   Evidence is append-only and may only be recorded for `APPROVED` jobs.

---

## 7. IDE Agent Directives (RULES FOR AI CODING ASSISTANT)
*Read these rules carefully before generating code for this workspace.*
1.  **Fail-Closed Only:** When writing `gatekeeper.py`, default to rejecting the payload if there is a parsing error, an unknown tool, or a missing cost estimate.
2.  **Never Bypass the Ledger:** No worker function (Modal, RAG, etc.) can be invoked directly by the LLM or UI. Everything MUST pass through the `authorize_spend()` function in `database.py`.
3.  **UI Asynchrony:** When writing `dashboard.py` in Gradio, never block the main thread waiting for a Modal job. Use Gradio's `yield` or background tasks to poll the SQLite database for job completion.
4.  **LLM Constraints:** When writing API calls to Meta Llama 3.1, enforce JSON output mode. The LLM must not generate conversational text in its operational responses.

---

## 8. Confirmed Design Decisions (v3.1)
| # | Topic | Decision |
|---|-------|----------|
| D1 | Money precision | Ledger stores INTEGER micro-USD; API uses float USD. Costs round up, allocations round down. |
| D2 | $0.00 jobs | `authorize_spend()` rejects any cost ≤ 0. Mock mode charges a small nominal synthetic cost. |
| D3 | Cost model | `(gpu_seconds_per_variant × variants + cold_start_seconds) × A10G_usd_per_second`; all constants in config, never from the LLM. |
| D4 | Modal failure | Job → `FAILED`; no refund (compute may have been billed); loop continues. |
| D5 | Budget scope | One run = one fresh ledger file, `Budget.id = 1`. |
| D6 | Demo target | avGFP + ESM-2 650M; answer key = ProteinGym GFP DMS (Sarkisyan 2016); Spearman metric. |
| D7 | LLM host | Meta Llama API (`LLAMA_API_KEY`). |
| D8 | Variants & stopping | Mutation strings on WT (e.g., `S65T:Y66H`); stop on budget exhaustion or a hard max-iteration cap. |

**Open items:** A10G $/sec, GPU-seconds per variant, cold-start seconds, max-iteration cap, mock nominal cost, Llama 3.1 size/model ID, Python version & dependency pins.