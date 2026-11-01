> **SUPERSEDED — ARCHIVED 2026-11-01. DO NOT QUOTE THIS FILE.**
>
> This is a planning or positioning document written before the measurements were taken.
> It is kept for history, not for reference, and parts of it are now known to be wrong:
> it describes gaps that have since been closed, and it carries claims that
> [`README.md`](../../README.md) and [`docs/HANDOFF.md`](../HANDOFF.md) explicitly
> withdraw. Where any of it disagrees with those two files, those two are correct.
>
> In particular, nothing here should be read as a statement of what the system was
> measured to do. The measured results, with their limitations, are in `README.md`; the
> run artifacts behind them are in [`evidence/`](../../evidence/README.md).
# VeritasBio â€” Frontier Gap Resolution & Productionization Execution Plan

**Document Status:** Approved Engineering & Commercial Execution Plan (v1.0)  
**Target Horizon:** Meta Global AI Developer Hackathon Victory $\to$ B2B Enterprise Bio-SaaS Spinout  
**Author:** VeritasBio Core Engineering & DeepMind Antigravity  
**Date:** October 28, 2026  

---

## 1. Executive Summary & System Objectives

VeritasBio has validated its foundational premise: **coupling a fail-closed micro-USD ledger with Gaussian Process active learning slashes computational biology cloud spend by 73.9xâ€“134x while attaining 80.1%â€“91.2% of library global optima across 4 authentic experimental landscapes (avGFP, TrpB, GB1, Trastuzumab/HER2).**

To elevate VeritasBio from an elite hackathon project into a **generationally impactful, venture-backed enterprise platform**, this document establishes the definitive execution plan to resolve the four fundamental operational and biophysical gaps:

```
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                                   THE 4 FRONTIER GAPS & STRATEGIC SOLUTIONS                            â”‚
â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
â”‚ Gap                           â”‚ Root Limitation                  â”‚ Technical Solution                  â”‚
â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
â”‚ 1. The "LLM as Decoy" Risk    â”‚ Llama 3.1 only suggests spots    â”‚ Closed-Loop Residual Reflection     â”‚
â”‚                               â”‚ at round 0; GP does 90% of math. â”‚ (Llama analyzes unpredicted errors) â”‚
â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
â”‚ 2. Turnaround & Batching      â”‚ 10-variant rounds take 4 months  â”‚ Dual-Mode Batching: In-Silico Rapid â”‚
â”‚    Scale Asymmetry            â”‚ in wet labs; 96-well plates need â”‚ ($B=15$) vs. Plate-Scale CRO        â”‚
â”‚                               â”‚ 88 variants per 10-day cycle.    â”‚ ($B=88$ via $k$-DPP Diversity)      â”‚
â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
â”‚ 3. Catalytic vs. Evolutionary â”‚ PLMs (ESM-2) measure natural     â”‚ 3D HETATM Cofactor / Ligand Sphere  â”‚
â”‚    Unnatural Substrate Gap    â”‚ fitness, not unnatural substrate â”‚ Geometry & Catalytic Clash Penalty  â”‚
â”‚                               â”‚ transition-state catalysis.      â”‚ (<4.5Ã… to PLP / transition state)   â”‚
â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
â”‚ 4. Enterprise Compliance &    â”‚ Local SQLite & Gradio cannot be  â”‚ PostgreSQL + Supabase RLS,          â”‚
â”‚    LIMS Multi-Tenancy         â”‚ deployed inside enterprise pharmaâ”‚ Benchling REST API integration, and â”‚
â”‚                               â”‚ R&D without audit logging.       â”‚ HMAC-SHA256 21 CFR Part 11 chain    â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”´â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”´â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
```

---

## 2. Technical & Mathematical Architecture per Gap

### Gap 1: Closed-Loop Mechanistic Residual Reflection (Meta Llama 3.1)

#### Objective
Elevate Meta Llama 3.1 from a one-time boundary proposer into an **active biophysical reasoner** that diagnoses active-learning failures and steers subsequent search trajectories.

#### Algorithmic Flow
At the completion of each active learning round $t$:
1. **GP Residual Computation:**
   For all measured variants $x_i \in \mathcal{D}_t$, compute the surrogate model prediction residual:
   $$e_i = y_{\text{empirical}}(x_i) - \hat{\mu}_{\text{GP}}(x_i)$$
   Rank variants into:
   - **Synergistic Outliers (Epistatic Goldmines):** Top positive residuals $e_i > +1.5\sigma_n$ where observed fitness vastly exceeds additive/smooth GP predictions.
   - **Destabilizing Outliers (Epistatic Crashes):** Top negative residuals $e_i < -1.5\sigma_n$ where predicted benign substitutions destroyed function.
2. **3D Structural Reflection Prompt:**
   Construct a prompt injecting the 3D atomic coordinates, inter-residue distances $d_{ij}$, and side-chain chemistry deltas for these outlier variants:
   ```json
   {
     "positive_outliers": [
       {"variant": "F184I:S228G", "predicted": 0.420, "observed": 0.8005, "residual": 0.3805, "distance_angstroms": 18.95}
     ],
     "negative_outliers": [
       {"variant": "V183K:F184W", "predicted": 0.510, "observed": 0.0820, "residual": -0.4280, "distance_angstroms": 3.82}
     ]
   }
   ```
3. **Llama Mechanistic Reflection Response:**
   Llama 3.1 outputs a structured `ReflectionSummary` JSON:
   - `mechanistic_diagnosis`: Biophysical explanation (e.g. *"F184I alleviated steric congestion at the loop entry, permitting S228G to relieve backbone strain without disrupting PLP Schiff base"*).
   - `updated_suggested_chemistries`: Refined amino acid preferences for remaining rounds.
   - `forbidden_chemistries`: Forbidden residue-chemistry combinations.
4. **Surrogate Prior Reweighting:**
   The GFlowNet trajectory sampler and candidate generator apply dynamic multiplicative bonuses:
   $$R_{\text{guided}}(x) = R_{\text{surrogate}}(x) \times \exp\left(\sum_{k} \gamma_k \cdot \mathbb{I}(x \in \text{ReflectionHypothesis}_k)\right)$$
   where $\gamma_k \in [0.10, 0.35]$.

---

### Gap 2: Dual-Mode Plate-Scale Asynchronous Batching

#### Objective
Eliminate the mismatch between rapid in-silico simulation ($B=10\text{â€“}25$, seconds) and physical wet-lab economics ($B=88$ or $B=308$, 10-day cycles in 96/384-well plates).

#### Mathematical & Operational Modes
The system introduces a top-level execution configuration:
```python
class CampaignMode(str, Enum):
    IN_SILICO_RAPID = "in_silico_rapid"  # Micro-batches (B=10-25), 6-8 rounds, live streaming UI
    PLATE_SCALE_CRO = "plate_scale_cro"  # Plate batches (B=88/308), k-DPP diversity, 1-2 wet-lab cycles
```

#### Plate-Scale ($B=88$) Specifications:
1. **Determinantal Point Process ($k$-DPP) Sizing:**
   Instead of greedily picking the top 88 Expected Improvement candidates (which would yield redundant sequence clones), $k$-DPP samples a subset $Y \subset \mathcal{X}$ of size $k=88$ proportional to:
   $$\mathcal{P}(Y) \propto \det(L_Y), \quad L_{ij} = q_i \cdot k_{\text{physicochemical}}(x_i, x_j) \cdot q_j$$
   where $q_i = \text{EI}(x_i)$ and $k(x_i, x_j) = \exp(-\|\mathbf{f}(x_i) - \mathbf{f}(x_j)\|^2 / 2\sigma^2)$.
2. **Robotic Anti-Evaporation Plate Architecture:**
   - **Outer Perimeter (36 Wells: A1-A12, H1-H12, B1-G1, B12-G12):** Dedicated to sterile $200\mu\text{L}$ water/PBS hydration barriers to stop edge-effect evaporation.
   - **Controls (8 Wells):** 1x WT parent standard, 1x null mutant, 2x known stabilizing variants, 4x media blanks.
   - **Screening Candidates:** Exactly 52 internal wells in 96-well format (or 88 wells in border-compensated format) and 308 wells in 384-well format.
3. **Asynchronous Ledger State Handshake:**
   - Ledger Job State transitions: `APPROVED` $\to$ `SYNTHESIS_ORDERED` $\to$ `WAITING_FOR_ASSAY` $\to$ `COMPLETED`.
   - Generates Twist DNA order CSV and Opentrons OT-2 script, pauses the loop, and awaits empirical plate reader data via webhook/CLI ingestion (`scripts/ingest_plate_data.py`).

---

### Gap 3: 3D HETATM Cofactor / Ligand Coordinate Sphere Geometry

#### Objective
Resolve the low zero-shot accuracy of foundation PLMs on unnatural catalytic reactions (e.g. TrpB synthetic indoles) by incorporating explicit ligand/cofactor 3D geometry into pre-scoring filters and the DKL surrogate.

#### Implementation Architecture:
1. **HETATM Parsing in `SpatialContactMap`:**
   Extend `workers/spatial_filter.py` to parse both `ATOM` and `HETATM` records:
   - Pyridoxal 5'-phosphate (`PLP`) in TrpB (PDB 8EH0).
   - Heme groups (`HEM`), Flavin (`FAD`/`FMN`), and synthetic ligands (`LIG`).
2. **Transition-State Proximity Filter:**
   Calculate the minimal Euclidean distance between mutated side-chain atoms and cofactor active atoms:
   $$d_{\text{cofactor}}(pos) = \min_{a \in \text{Cofactor}} \|r_{pos, C\alpha} - r_a\|_2$$
   - **Zone 1: Catalytic Core ($d < 4.5\text{ Ã…}$):** Mutations to bulky aromatics (Trp/Tyr/Phe) or charged residues that disrupt electrostatic binding to the phosphate/ring trigger an immediate **Catalytic Geometry Clash Penalty** ($-3.5$ in surrogate score, or fail-closed pruning).
   - **Zone 2: Entrance Channel ($4.5\text{ Ã…} \le d \le 12.0\text{ Ã…}$):** Favorable region for steric relaxation (e.g. F184I/V) to accommodate non-cognate synthetic substrates.
3. **Substrate Physicochemical Conditioning in DKL Head:**
   Concatenate substrate descriptor vectors (molecular weight, logP, TPSA, H-bond donors/acceptors) with the 5-factor residue delta embeddings, conditioning the Gaussian Process on the specific chemical substrate.

---

### Gap 4: Enterprise Productionization, LIMS & Compliance

#### Objective
Transition VeritasBio from a single-user local tool to an enterprise-grade cloud B2B platform deployable in pharma VPCs with complete regulatory auditability.

#### Production Architecture:
1. **Dual-Backend Ledger Interface:**
   Refactor `ledger/database.py` with an abstract `LedgerBackend` interface:
   - `SQLiteLedger`: Zero-config local execution for unit tests, hackathon evaluations, and developer workstations.
   - `PostgresLedger`: High-concurrency enterprise cloud ledger with Row-Level Security (RLS) scoping jobs and evidence by `organization_id` and `project_id`.
2. **Benchling Registry & Inventory Integration:**
   - Service module `workers/lims_worker.py`:
     - Calls `POST /api/v2/dna-sequences` to register codon-optimized lead constructs with full gene annotations and mutation tags.
     - Calls `POST /api/v2/plates` to populate 96-well virtual plate layouts mirroring the physical Opentrons layout.
3. **FDA 21 CFR Part 11 Cryptographic Audit Trail:**
   - Every ledger transaction generates an HMAC-SHA256 digital signature:
     $$\text{Signature} = \text{HMAC}_{\text{KMS}}\left(\text{timestamp} \parallel \text{job\_id} \parallel \text{variant\_count} \parallel \text{cost\_micro\_usd} \parallel \text{status}\right)$$
   - Implemented as a tamper-evident hash-linked block chain in SQL (`audit_ledger_blocks`), satisfying FDA predicate rules for electronic records and patent inventorship.

---

## 3. Phased Implementation Roadmap

```
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                                           EXECUTION SCHEDULE                                           â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
  PHASE 1: IMMEDIATE GAP RESOLUTION (Today)
  â”œâ”€â”€ Task 1.1: Implement Llama 3.1 Closed-Loop Residual Reflection in controller/
  â”‚             - Prompts: Add ReflectionSummary, residual ranking, and reflection prompt builder.
  â”‚             - Orchestrator: Hook residual computation between Step 5 and Step 2 of next round.
  â”œâ”€â”€ Task 1.2: Upgrade workers/spatial_filter.py for HETATM Cofactor / Ligand Parsing
  â”‚             - Parse HETATM PLP coordinates in PDB 8EH0.
  â”‚             - Add cofactor proximity & catalytic distortion scoring.
  â”œâ”€â”€ Task 1.3: Add Dual-Mode Batching (In-Silico Rapid vs Plate-Scale CRO)
  â”‚             - Support campaign_mode in OrchestratorConfig.
  â”‚             - Implement 96-well plate-scale generator with k-DPP diversity.
  â””â”€â”€ Task 1.4: Verify Full Test Suite & Benchmark Stability (Zero Regressions).

  PHASE 2: ENTERPRISE PRODUCTIONIZATION (Tomorrow / Day 2)
  â”œâ”€â”€ Task 2.1: Abstract Ledger Interface & Implement PostgreSQL / Supabase Schema.
  â”œâ”€â”€ Task 2.2: Implement workers/lims_worker.py with Benchling API mock/live adapters.
  â””â”€â”€ Task 2.3: Cryptographic 21 CFR Part 11 Hash Chain Audit Implementation.

  PHASE 3: COMMERCIALIZATION & GO-TO-MARKET (Post-Hackathon)
  â”œâ”€â”€ Task 3.1: Launch Commercial Biologics Wedge (Antibody affinity maturation & developability).
  â”œâ”€â”€ Task 3.2: Formalize Twist Bioscience & Automated CRO partnerships (Compute-to-Plate).
  â””â”€â”€ Task 3.3: Deploy "FinOps for Bio" value-based pricing (Base SaaS + 20% compute savings).
```

---

## 4. Acceptance Criteria & Test Verification Plan

Every phase must satisfy deterministic automated verification before promotion:

| Component | Test File | Target Metric / Acceptance Criteria |
| :--- | :--- | :--- |
| **Llama Reflection** | `tests/test_llama_reflection.py` | Validates positive/negative residual extraction, reflection prompt structure, and prior bonus injection into GFlowNet. |
| **HETATM Cofactor** | `tests/test_spatial_filter.py` | Confirms PLP cofactor atom parsing from PDB 8EH0; rejects catalytic distortion mutations $<4.5\text{ Ã…}$ from PLP. |
| **Plate-Scale Batching** | `tests/test_plate_scale_batching.py` | Verifies $B=88$ batch generation, $k$-DPP determinant volume $>0$, anti-evaporation perimeter preservation, and Opentrons well assignment. |
| **Full Benchmark Suite** | `scripts/run_side_by_side_benchmarks.py` | Zero performance degradation on TrpB ($\rho \ge +0.64$) and avGFP ($\rho \ge +0.79$). |

---
*Document approved for immediate codebase integration.*


