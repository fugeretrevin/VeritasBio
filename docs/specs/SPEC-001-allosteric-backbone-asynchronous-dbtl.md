# SPEC-001: Next-Generation Biophysical Hardening & Asynchronous DBTL Engine

- **Status:** Draft / Ready for Review
- **Author:** VeritasBio Core Engineering Team
- **Date:** October 29, 2026
- **Domain:** Structural Biophysics, Bayesian Optimization, FinOps Ledger, Automated Physical Biology
- **Triage Label:** `ready-for-agent`

---

## Problem Statement

VeritasBio v3.2 demonstrated elite active-learning metrics across four benchmark datasets (avGFP, TrpB, GB1, Trastuzumab), achieving empirical global optimum discovery on TrpB (7.1597 vs 1.0000 WT) with 13.0x top-1% enrichment and sub-four-cent compute economics ($0.0376 total spend).

However, rigorous analysis of the production results revealed six operational and biophysical limitations:

1. **Local-Shell Myopia:** Search boundary proposals from the LLM are confined to primary and secondary contact shells (< 8 Å from active sites or fluorophores), remaining blind to long-range allosteric hinge residues (> 15 Å) known to modulate catalytic turnover.
2. **Backbone Rigidity Bottleneck:** Track C sidechain relaxation maintains an ultra-stiff harmonic C-alpha restraint (k = 500 kcal/mol·Å²), which successfully prevents false-positive clashes but incorrectly penalizes loop remodelings or proline/glycine mutations requiring slight backbone accommodation (> 0.2 Å displacement).
3. **Synchronous Run Blocking:** In silico iteration takes 2–3 minutes, whereas physical cDNA synthesis (Twist) and liquid-handling assay execution (Opentrons/Tecan) require 24 to 120 hours. The current orchestrator requires synchronous completion of each round before proposing the next batch, stalling active discovery.
4. **Wet-Lab Experimental Replicate Noise:** Ground-truth benchmark assays evaluate idealized floating-point fitness scores. Real plate-reader absorbance/fluorescence data exhibits edge-well evaporation, pipetting drift, and replicate variance (heteroscedastic noise) that can distort Expected Improvement acquisition if not filtered.
5. **Manual Remote Execution Boundary:** While the system exports Twist CSV, Opentrons OT-2 Python scripts, and Benchling manifests, execution on remote cloud laboratories (Emerald Cloud Lab, Strateos) requires human file-handling rather than automated webhook dispatch.

---

## Solution

Build the **Biophysical Hardening & Asynchronous DBTL Engine (v3.3)**:

1. **Elastic Network / Allosteric Hinge Detection (`workers/spatial_filter.py`):**
   Compute anisotropic network model (ANM) / Gaussian network model (GNM) normal mode fluctuations from PDB coordinates. Identify dynamic hinge residues and solvent-exposed allosteric communication channels, injecting them directly into Llama 3.1 combinatorial search prompts alongside active-site contact residues.
2. **Two-Stage Elastic Backbone Minimizer (`workers/relaxation_worker.py`):**
   Introduce a two-stage local relaxation pipeline:
   - *Stage 1 (Sidechain Fast L-BFGS):* Rigid backbone ($k_{CA} = 500$), optimizing sidechain chi angles in < 40ms.
   - *Stage 2 (Elastic Mainchain Accommodation):* If residual Lennard-Jones repulsive strain exceeds $20.0\text{ kcal/mol}$, relax backbone coordinates under soft harmonic restraints ($k_{CA} = 5.0\text{ kcal/mol}\cdot\text{\AA}^2$) with a maximum permissible $C_\alpha$ displacement threshold of $1.2\text{ \AA}$, preventing unphysical unfolding while rescuing loop remodeling mutations.
3. **Asynchronous Ledger State Machine (`ledger/database.py`, `controller/orchestrator.py`):**
   Add asynchronous job lifecycle states (`STAGED_FOR_WETLAB`, `IN_FLIGHT_SYNTHESIS`, `INGESTION_PENDING`). Allow the Bayesian quant engine to stage speculative exploratory batches for Round $N+1$ without blocking on physical validation of Round $N$, using "hallucinated fantasy observations" (Kriging believer heuristic) with atomic ledger rollbacks upon real plate data ingestion.
4. **Robust Replicate Ingestion & Outlier Rejection (`scripts/ingest_plate_data.py`, `server.py`):**
   Implement automated replicate outlier filtering (Huber robust M-estimation and Grubbs' test, $\alpha = 0.05$) across multi-well replicates before normalizing against the WT baseline (`B02`), writing variance-weighted observations into the heteroscedastic GP kernel.
5. **Programmatic Cloud-Lab Adapter (`workers/cloudlab_worker.py`):**
   Implement client adapters for Emerald Cloud Lab (Symbolic Lab Language generator) and Strateos JSON protocol dispatch, complete with webhook ingestion (`/api/cloudlab/webhook`) for automated microplate reader callback.

---

## User Stories

1. As a computational protein engineer, I want the system to calculate normal mode hinge dynamics from the target PDB, so that Llama 3.1 proposes allosteric mutations outside the immediate active site.
2. As a structural biologist, I want sidechain relaxation to accommodate subtle backbone loop movements (up to 1.2 Å) when bulky mutations are introduced, so that valid candidate leads are not falsely rejected due to rigid crystal structure constraints.
3. As a wet-lab director, I want the active learning orchestrator to generate speculative candidate batches while physical plates are in the incubator, so that our team does not lose 5 business days waiting between design cycles.
4. As a quality assurance officer, I want multi-well plate reader data to automatically identify and drop contaminated or evaporated edge wells before calculating fitness, so that noisy experimental outliers do not distort the Gaussian Process model.
5. As an enterprise CRO customer, I want to dispatch synthesis and assay protocols directly to Strateos or Emerald Cloud Lab via REST API/SLL, so that the entire design-build-test-learn cycle executes autonomously without manual CSV handling.
6. As a laboratory automation scientist, I want in-flight batches to use Kriging believer fantasy sampling, so that proposed exploratory batches do not redundantly re-sample already-synthesized sequence motifs.
7. As a finance officer, I want speculative in-flight batches to require cryptographic budget pre-authorization in the SQLite ledger, so that pipelined synthesis orders never exceed authorized capital allocation.
8. As a compliance auditor, I want every asynchronous state transition (`STAGED_FOR_WETLAB` -> `COMPLETED`) to append an immutable HMAC-SHA256 Merkle block with operator timestamp, so that FDA 21 CFR Part 11 integrity is maintained across weeks of wet-lab execution.
9. As a frontend user, I want the web dashboard to visually distinguish between verified experimental points, speculative in-flight candidates, and in silico models on the Pareto frontier plot, so that I can clearly track physical validation progress.
10. As a biocatalysis researcher, I want the allosteric hinge filter to flag residues that couple to catalytic cofactors (PLP, HEM, NADP), so that distal mutations that destabilize cofactor binding are automatically penalized.
11. As a machine learning engineer, I want the plate ingestion endpoint to supply per-variant standard error ($\sigma_i$) to the heteroscedastic Gaussian Process, so that uncertain experimental measurements are appropriately down-weighted during acquisition optimization.
12. As a platform administrator, I want cloud-lab webhook callbacks to authenticate via signed HMAC headers, so that unauthorized external services cannot inject unverified fitness scores into the ledger.

---

## Implementation Decisions

### 1. Module Boundaries & Interfaces

- **`workers/spatial_filter.py`**:
  - Expose `compute_anm_hinge_residues(pdb_path: Path, cutoff_distance: float = 12.0) -> List[int]`.
  - Computes Kirchohoff / Hessian matrix on $C_\alpha$ coordinates, solves low-frequency eigenvectors (modes 7–12), and identifies residues with high mean-square fluctuation and hinge transition points.
  - Return `allosteric_hinge_positions` alongside existing `contact_shell_positions`.
- **`workers/relaxation_worker.py`**:
  - Update `relax_mutant_structure` to accept `allow_backbone_accommodation: bool = True` and `strain_threshold: float = 20.0`.
  - When initial steric clash energy $E_{\text{clash}} > \text{strain\_threshold}$, trigger Stage 2 relaxation with $k_{CA} = 5.0\text{ kcal/mol}\cdot\text{\AA}^2$.
  - Enforce hard cutoff: if resulting $C_\alpha\text{ RMSD} > 1.2\text{ \AA}$, reject the structural candidate to prevent unphysical denaturation.
- **`ledger/schema.sql` & `ledger/database.py`**:
  - Add `Jobs.status` enum values: `STAGED_FOR_WETLAB`, `IN_FLIGHT_SYNTHESIS`.
  - Add table `SpeculativeEvidence` for Kriging fantasy observations with `is_provisional = 1`.
  - Add function `commit_physical_evidence_and_reconcile(job_id, observed_results)` which atomically converts provisional rows to verified `Evidence` and triggers Merkle audit chain update.
- **`workers/cloudlab_worker.py`**:
  - Implement `StrateosClient` and `EmeraldCloudLabClient` with mock fallback for testing.
  - Generate protocol specifications matching Opentrons OT-2 volume transfers and Tecan absorbance/fluorescence plate reader layouts.
- **`server.py`**:
  - Add endpoint `POST /api/cloudlab/dispatch` (dispatches protocol and updates ledger status to `IN_FLIGHT_SYNTHESIS`).
  - Add endpoint `POST /api/cloudlab/webhook` (authenticates payload signature and ingests finished plate reader data).
  - Add robust outlier rejection to `POST /api/campaign/ingest-plate`.

---

## Testing Decisions

- **Test Behavior, Not Internal Solvers:**
  - Test that allosteric hinge calculations return physically valid residue indices within sequence boundaries and correctly identify published hinge residues on standard benchmark structures (PDB: `1EMA`, `8EH0`).
  - Test that the two-stage relaxation successfully relaxes a deliberately jammed bulky mutation (e.g. `G67W` in GFP) without causing $C_\alpha\text{ RMSD}$ to exceed $1.2\text{ \AA}$.
  - Test that asynchronous batch staging prevents duplicate exploration of in-flight sequence space while adhering strictly to fail-closed ledger spend authorization.
  - Test that noisy plate reader data with deliberate outlier wells (e.g., $10\times$ baseline noise in an edge well) is cleansed by Huber/Grubbs filtering without corrupting the WT reference well (`B02`).
  - Test that FDA Part 11 Merkle hash chains remain unbroken across asynchronous provisional-to-permanent state transitions.
- **Prior Art in Codebase:**
  - `tests/test_relaxation_worker.py` (L-BFGS minimization and $C_\alpha$ RMSD assertions).
  - `tests/test_enterprise_compliance.py` (Merkle audit block validation and SQLite triggers).
  - `tests/test_server_api.py` (FastAPI route validation and error contracts).

---

## Out of Scope

- Full-atom microsecond Molecular Dynamics (MD) trajectories with explicit solvent (GROMACS/AMBER) — this engine prioritizes sub-second deterministic FinOps screening.
- Live physical API credential provisioning for proprietary commercial CROs during test suites (all external calls use high-fidelity mock adapters with contract parity).
- Complete de novo backbone generation from scratch (handled separately by `workers/diffusion_worker.py`).

---

## Further Notes

- Maintains 100% backward compatibility with existing SQLite ledger schema and Gradio/React UI.
- All monetary quantities continue to use integer micro-USD ($1\text{ USD} = 1,000,000\ \mu\text{USD}$) with monotonic spend invariant.
- All new tests must pass alongside the existing 188-test regression suite.

