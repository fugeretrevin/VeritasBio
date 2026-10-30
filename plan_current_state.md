Analyzed diffusion_worker.py:L1-180
Analyzed diffusion_worker.py:L181-450
Analyzed diffusion_worker.py:L451-597

Here is a complete overview of the **PanoptesBio / VeritasBio** architecture as it stands today, how our **low-cost binder discovery & synthesis engine** compares to high-cost brute-force generative approaches, and the **remaining technical gaps & next steps**.

---

## 1. Full Inventory of What We Have Built

### A. Step 0–6 Active Learning & Biophysical Core
| Layer | Module & Key Symbols | Current Capabilities |
| :--- | :--- | :--- |
| **Step 0: Spatial & Allosteric Pre-Filter** | [`workers/spatial_filter.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py) ([`SpatialContactMap`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py#L43-L250), [`compute_anm_hinge_residues`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py#L253-L320), [`compute_interface_binding_ddg`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py#L561-L664)) | Parses multi-chain `ATOM` & `HETATM` records, protects cofactor catalytic shells ($4.2\text{ \AA}$), computes Kirchhoff Laplacian normal-mode allosteric hinges (`SPEC-001`), and scores quaternary antibody–target interface $\Delta\Delta G_{\text{bind}}$ (`SPEC-002`). |
| **Step 0b: 2-Stage Elastic Structure Relaxation** | [`workers/relaxation_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/relaxation_worker.py) ([`relax_variant_structure`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/relaxation_worker.py#L140-L380)) | Fast CPU L-BFGS local minimizer (`<50ms`) with Stage 1 rigid-$C_\alpha$ sidechain rotamer cone relaxation ($k=500$) and Stage 2 soft mainchain accommodation ($k=5.0$, capped at $C_\alpha\text{ RMSD} \le 1.20\text{ \AA}$) to rescue false-negative packing mutations (`SPEC-001`). |
| **Step 1–2: RAG + Closed-Loop Llama 3.1 Planner** | [`workers/rag_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/rag_worker.py), [`controller/prompts.py`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/prompts.py) ([`build_closed_loop_reflection_block`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/prompts.py#L96-L185)) | Fetches UniProt/PDB domain constraints, injects prior round GP prediction residuals ($y_{\text{obs}} - \hat{\mu}_{\text{GP}}$) and ANM hinge indices, and outputs strict JSON residue search boundaries. |
| **Step 3: GFlowNet + DKL GP + Pareto Quant Engine** | [`workers/math_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py) ([`sample_gflownet_candidates`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L465-L620), [`ESMReadoutAdapter`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L1408-L1529), [`select_next_micro_batch`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L1057-L1200), [`compute_pareto_frontier`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L1205-L1350)) | Samples multi-mutant trajectories via GFlowNet, projects 1280-dim embeddings through a supervised Ridge + SVD adapter (`SPEC-002`), conditions the Heteroscedastic GP on Kriging-believer in-flight assays (`SPEC-002`), and ranks candidates across a 4-objective NSGA-II Pareto frontier (Fitness, Developability, Synthesizability, Stability). |
| **Step 4: Fail-Closed Gatekeeper & Ledger** | [`controller/gatekeeper.py`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/gatekeeper.py) ([`Gatekeeper`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/gatekeeper.py#L150-L390), [`estimate_wetlab_synthesis_cost`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/gatekeeper.py#L238-L275)), [`ledger/database.py`](file:///c:/Users/stoof/source/repos/VeritasBio/ledger/database.py) | Integer micro-USD accounting, two-state cold/warm GPU cost model, wet-lab synthesis + CRO cost model (`SPEC-002`), `SpeculativeEvidence` in-flight staging (`SPEC-001`), and FDA 21 CFR Part 11 HMAC-SHA256 Merkle hash chain. |
| **Step 5: Cloud GPU Scoring & De Novo Binder Engine** | [`workers/esm_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/esm_worker.py), [`workers/diffusion_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py) ([`generate_budgeted_de_novo_binders`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L417-L597)) | Modal A10G ESM-2 650M masked marginal scoring + parametric 3-helix bundle backbone generator (`GeometricBackboneSampler`) and `ProteinMPNNInverseFolder`. |
| **Step 6: Wet-Lab CRO, Cloud-Lab & Curve Ingestion** | [`workers/wetlab_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/wetlab_worker.py), [`workers/cloudlab_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/cloudlab_worker.py), [`scripts/ingest_plate_data.py`](file:///c:/Users/stoof/source/repos/VeritasBio/scripts/ingest_plate_data.py) | Twist cDNA CSV export, Opentrons OT-2 Python protocol, Benchling LIMS sync, Emerald Cloud Lab SLL & Strateos Autoprotocol dispatch (`SPEC-001`), Grubbs/Huber outlier filter (`SPEC-001`), and 4PL $\text{IC}_{50}/K_d$ + Michaelis-Menten kinetic curve fitting (`SPEC-002`). |

---

## 2. How We Approach Ultra-Low-Cost De Novo Binder Discovery & Synthesis

Standard generative protein binder campaigns typically spend **$\$15,000\text{--}\$40,000+$** because they generate $1,000\text{--}10,000$ unscreened backbone/sequence hallucinations and synthesize hundreds of individual clonal genes at $\$18\text{--}\$45$ each, only to see $95\%+$ fail in wet-lab binding assays due to:
1. **Unpaired cysteines, deamidation (`NG`/`NS`), or isomerization (`DG`/`DS`)** causing aggregation on the plate.
2. **Steric clashes or buried unsatisfied polar atoms** at the target–binder interface.
3. **Lack of closed-loop affinity maturation** (doing one massive one-shot synthesis batch instead of small, high-information micro-batches).

Our architecture solves this with a **3-tier cost-compression funnel**:

1. **Pre-Compute Biophysical Sieve (`$0.00` CPU Cost)**:
   - Instead of sending thousands of raw diffusion trajectories to GPU folding or wet-lab synthesis, we dock compact 65-residue 3-helix bundle mini-proteins ([`GeometricBackboneSampler`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L129-L318)) or existing humanized antibody frameworks (`2YBR.pdb`) against the target epitope hotspot centroid, inverse-fold sequences free of liability motifs ([`ProteinMPNNInverseFolder`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L323-L411)), relax interfacial sidechains ([`relax_variant_structure`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/relaxation_worker.py#L140-L380)), and rank by quaternary interface $\Delta\Delta G_{\text{bind}}$ ([`compute_interface_binding_ddg`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py#L561-L664)).
2. **Micro-Batch GPU Verification (`< $0.01` per Round)**:
   - Only the top Pareto-optimal candidates pass the [`Gatekeeper`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/gatekeeper.py#L150-L390) spend check to run on Modal A10G GPUs (`$0.000306/sec`).
3. **Active-Learning Micro-Plates (`12–24` Variants Instead of `1,000+`)**:
   - By pairing the [`ESMReadoutAdapter`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L1408-L1529) with Kriging-believer batch selection ([`select_next_micro_batch`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/math_worker.py#L1057-L1200)), we only need to synthesize **8 to 16 high-confidence variants per DBTL round** rather than screening thousands of random designs.

---

## 3. Remaining Architectural Gaps & What We Should Do Next

Inspecting [`workers/diffusion_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py), [`workers/wetlab_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/wetlab_worker.py), and [`web/src/`](file:///c:/Users/stoof/source/repos/VeritasBio/web/src/) reveals **4 concrete gaps** that will complete the de novo synthesis pipeline and unify the frontend with the `SPEC-001`/`SPEC-002` backend:

### Gap 1: De Novo Physics & Self-Consistency Scoring in [`workers/diffusion_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py)
- **Current State**:
  1. In [`GeometricBackboneSampler.sample_3helix_bundle`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L284-L310), the generated PDB string hardcodes `GLY` as the residue name for every atom (`ATOM ... CA GLY B...`) because sequence inverse folding happens *after* coordinate generation.
  2. In [`generate_budgeted_de_novo_binders`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L521-L525), `predicted_plddt` and `estimated_rmsd_to_target` use a random uniform heuristic (`sampler.rng.uniform(-2.0, 3.0)`) rather than computing actual **Quaternary Interface $\Delta\Delta G_{\text{bind}}$**, **2-Stage L-BFGS Relaxation Energy**, and **Self-Consistency Packing Score**.
- **Concrete Fix**:
  - Thread the inverse-folded 1-letter amino acid sequence back into the PDB `ATOM` records as authentic 3-letter residue codes (`ALA`, `LEU`, `TYR`, `TRP`, etc.) and construct a combined **Target (Chain A) + Designed Binder (Chain B) complex PDB**.
  - Score every de novo binder with [`compute_interface_binding_ddg`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py#L561-L664) and [`relax_variant_structure`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/relaxation_worker.py#L140-L380) so only designs with favorable interfacial free energy ($\Delta\Delta G_{\text{bind}} < 0$) and low self-consistency RMSD are accepted.

### Gap 2: De Novo $\to$ Active Learning Hand-Off (Two-Stage Discovery Funnel)
- **Current State**: De novo generation ([`workers/diffusion_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py)) stops after generating initial scaffolds; it cannot hand off the top designed mini-binder into the 6-step GFlowNet + DKL GP active learning loop ([`controller/orchestrator.py`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/orchestrator.py)) for iterative affinity maturation.
- **Concrete Fix**:
  - Add a `seed_active_learning_from_denovo(design: DeNovoBinderDesign, target_pdb: str)` bridge that registers the top de novo 3-helix bundle as a dynamic campaign target (with its hotspot interface residues automatically mapped as the search space) so the 6-step GP/GFlowNet engine can affinity-mature the de novo binder over 3 rounds.

### Gap 3: Multiplexed Oligo-Pool & Fragment Assembly Synthesis Mode (`10x` Cheaper Physical DNA)
- **Current State**: [`workers/wetlab_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/wetlab_worker.py) and [`workers/diffusion_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py#L527) always price synthesis as individual linear clonal genes at $\$0.09/\text{bp}$ ($\approx \$17.55$ per 65-aa mini-binder or $\$65+$ per full-length scFv).
- **Concrete Fix**:
  - Add a **Multiplexed Oligo Pool / Golden Gate Assembly mode** to [`workers/wetlab_worker.py`](file:///c:/Users/stoof/source/repos/VeritasBio/workers/wetlab_worker.py) and [`estimate_wetlab_synthesis_cost`](file:///c:/Users/stoof/source/repos/VeritasBio/controller/gatekeeper.py#L238-L275) for short mini-binders ($\le 80\text{ aa}$ / $240\text{ bp}$) and localized CDR cassettes, with automatic Golden Gate `BsaI` overhang flanking sequences. This drops physical DNA cost from **$\sim\$18\text{--}\$65/\text{variant}$ down to $\sim\$0.85\text{--}\$2.50/\text{variant}$**.

### Gap 4: Frontend Studio Integration for `SPEC-001`, `SPEC-002`, and De Novo Binder Design
- **Current State**: The backend has endpoints for `/api/cloudlab/dispatch`, `/api/campaign/ingest-curves`, and `generate_budgeted_de_novo_binders`, but [`web/src/components/WetlabExportStudio.jsx`](file:///c:/Users/stoof/source/repos/VeritasBio/web/src/components/WetlabExportStudio.jsx) and [`web/src/components/StudioView.jsx`](file:///c:/Users/stoof/source/repos/VeritasBio/web/src/components/StudioView.jsx) do not yet expose interactive panels for:
  1. Launching a **Budgeted De Novo Mini-Binder** run and viewing the docked Chain A + Chain B 3D complex in [`Protein3DViewer.jsx`](file:///c:/Users/stoof/source/repos/VeritasBio/web/src/components/Protein3DViewer.jsx).
  2. Fitting **4PL $\text{IC}_{50}/K_d$ & Michaelis-Menten curves** in the Wet-Lab Studio.
  3. Dispatching **Emerald Cloud Lab SLL / Strateos Autoprotocol** payloads and viewing live `SpeculativeEvidence` in-flight status.

---

Would you like me to draft **`SPEC-003`** and implement these 4 upgrades (De Novo Quaternary $\Delta\Delta G_{\text{bind}}$ & Complex PDB Threading, De Novo $\to$ Active Learning Maturation Bridge, Multiplexed Oligo-Pool Golden Gate Synthesis, and Studio UI Wiring)?