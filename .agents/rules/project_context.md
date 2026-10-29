---
trigger: always_on
description: PanoptesBio project memory — confirmed scope decisions, deviations from ARCHITECTURE.md, invariants, and build status. Read before any change.
---

# PanoptesBio — Project Context (agent memory)

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
| D10 | Target modality & commercial roadmap | **Therapeutic biologics first** (antibodies, nanobodies, binders, folding stability) for immediate go-to-market where ESM-2 and DKL GP natively attain proven $\rho \approx 0.80$ predictive skill. Unnatural catalytic enzyme adapters, enterprise value-based pricing, and live CRO lab execution staged for post-hackathon expansion. |
| D11 | Scorpion Antivenom Neutralizing Engine | **Scorpion antivenom scaffold-based active learning** (`SCFV_9004G_CN2`, PDB: `2YBR`), targeting *Centruroides noxius* Cn2 neurotoxin. Solves antivenom development without de novo hallucination or $30k+ unguided brute force, discovering high-potency neutralizing leads under $0.01 compute with Twist CSV and Opentrons robot protocols. |
| D12 | PanoptesBio Firebase Authentication & FDA Part 11 Electronic Identity | **Firebase Authentication & Cloud Storage integration** (`panoptesbio`, `panoptesbio.firebaseapp.com`, `G-3CQM0RGVH3`). `controller/firebase_config.py` supplies credentials; `controller/auth.py` provides `AuthenticatedUser` model, token verification, and FDA Part 11 dual-identifier HMAC-SHA256 electronic signature stamping; `ledger/database.py` and `controller/gatekeeper.py` stamp operator identity directly into immutable SQLite Merkle audit blocks. |


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
- [x] Full test discovery passing: 117/117 unit tests across all test suites in 17.8s.
- [x] TrpB second benchmark: `scripts/download_trpb_benchmark.py`, `scripts/run_trpb_benchmark.py`,
      `tests/test_trpb_benchmark.py` (159,999 empirical variants, dataset-integrity assertions,
      held-out evaluation + random-screening baseline, 8 unit tests passing).
- [x] Enterprise hardening & biologics developability: `workers/math_worker.py`, `workers/wetlab_worker.py`,
      `tests/test_enterprise_hardening.py` (Heteroscedastic Gaussian Process with replicate variance,
      therapeutic antibody developability liability scoring, 96/384-well anti-edge evaporation barrier plate layout).
- [x] Gap 4 Enterprise Compliance & LIMS Integration: `ledger/schema.sql`, `ledger/database.py`,
      `workers/lims_worker.py`, `ledger/postgres_schema.sql`, `ui/dashboard.py`, `tests/test_enterprise_compliance.py`
      (FDA 21 CFR Part 11 cryptographic Merkle hash chain with HMAC-SHA256, schema immutability triggers,
      Benchling LIMS client with cDNA registration + virtual plate layout + manifest export,
      multi-tenant PostgreSQL schema with Row-Level Security, 8 unit tests passing, full suite 148/148 passing).
- [x] Scorpion Antivenom Neutralizing Antibody Campaign: `workers/data/SCFV_9004G_CN2.json`, `workers/data/2YBR.pdb`,
      `scripts/run_scorpion_antivenom_campaign.py`, `ui/structure_viewer.py`, `ui/dashboard.py`, `tests/test_scorpion_antivenom.py`
      (Authentic PDB 2YBR Cα contact graph, CDR-H3 paratope active learning, biologics developability liability screening,
      Twist cDNA CSV, Opentrons OT-2 neutralization protocol, FDA Part 11 certificate, 5 unit tests passing, full suite 153/153 passing).
- [x] Budgeted De Novo Diffusion & Inverse Folding Engine: `controller/gatekeeper.py`, `workers/diffusion_worker.py`,
      `tests/test_diffusion_worker.py` (Parametric 3-helix bundle geometric backbone sampling, ProteinMPNN inverse folding,
      zero unpaired cysteines & developability liability filtering, epitope hotspot docking, fail-closed spend gate,
      Twist cDNA CSV synthesis orders, Opentrons OT-2 scripts, FDA Part 11 audit certificates, 5 unit tests passing).
- [x] PanoptesBio Firebase Authentication & FDA Part 11 Electronic Signatures: `controller/firebase_config.py`,
      `controller/auth.py`, `ledger/database.py`, `controller/gatekeeper.py`, `ui/dashboard.py`, `tests/test_firebase_auth.py`
      (Multi-tenant Firebase project integration `panoptesbio.firebaseapp.com`, AuthenticatedUser model,
      FDA 21 CFR §11.200 HMAC-SHA256 electronic signature stamping into immutable SQLite Merkle audit blocks,
      Dashboard Tab 9 and connection status badge, 5 unit tests passing, full suite 169/169 passing in 114s).
- [x] Frontier Gap Resolution & Full-Stack Platform: `workers/spatial_filter.py`, `controller/prompts.py`, `server.py`, `scripts/ingest_plate_data.py`, `ui/structure_viewer.py`, `web/src/components/CampaignStudio.jsx`, `web/src/components/WetlabExportStudio.jsx`, `web/src/components/Protein3DViewer.jsx` (Closed-loop Llama 3.1 residual reflection, HETATM cofactor catalytic core clash checks with 4.2 Å proximity shell, real-time Server-Sent Events campaign streaming, drag-and-drop microplate reader CSV ingestion with hydration barrier validation, dual-viewport WT vs Lead 3Dmol WebGL comparison, and Twist/Opentrons/Benchling/FDA Part 11 export studio; full suite 180/180 passing in 40s).
- [x] Track B Multi-Objective Pareto Optimization & Developability Loss: `workers/math_worker.py`, `server.py`, `web/src/components/ParetoFrontierStudio.jsx`, `web/src/App.jsx`, `tests/test_pareto_optimization.py`, `tests/test_server_api.py` (NSGA-II Kung/Deb fast non-dominated sorting, 4-objective trade-off surface [Fitness, Developability, Synthesizability, Stability], 2D hypervolume indicator, `/api/campaign/pareto-frontier` REST endpoint, interactive trade-off SVG canvas with live weight sliders & presets, 183/183 unit tests passing in 39.2s).


## Open questions (ask the user when relevant)
- Numeric constants for D3: A10G $/sec, GPU-seconds/variant, cold-start seconds.
- Value of the max-iteration cap (D8) and the mock nominal cost (D2).
- Llama 3.1 variant (8B / 70B / 405B) and exact Meta Llama API model ID.
- Python version and whether to pin `requirements.txt`.

