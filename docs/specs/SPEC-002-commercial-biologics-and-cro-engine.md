# SPEC-002: Commercial Biologics, Multi-Chain Interface ΔΔG, 1280-Dim ESM Readout & Kinetic Curve Engine

- **Status:** Approved for Implementation
- **Author:** VeritasBio Core Engineering Team
- **Date:** October 29, 2026
- **Domain:** Multi-Chain Structural Biophysics, Foundation Model Readout Adapters, Dose-Response Kinetics, Dual Compute/Synthesis FinOps
- **Triage Label:** `ready-for-agent`

---

## Overview

SPEC-002 upgrades the 6-step PanoptesBio workflow for commercial B2B SaaS and autonomous CRO deployment across four vertical slices:

1. **Seam 1 (Multi-Chain Complex & Epitope-Paratope Interface $\Delta\Delta G_{\text{bind}}$):**
   Extend `workers/spatial_filter.py` to parse multi-chain PDB complexes (e.g., antibody `H`/`L` chains bound to antigen `A` chain in `1N8Z` and `2YBR`), exposing `compute_interface_contact_residues()` and `compute_interface_binding_ddg()`.
2. **Seam 2 (1280-Dim ESM-2 Embedding Readout Adapter):**
   Implement `ESMReadoutAdapter` in `workers/math_worker.py` to project high-dimensional foundation model representations into low-rank biophysical features trained on empirical assay data, overcoming zero-shot natural evolutionary bias on non-natural substrates.
3. **Seam 3 (Kriging-Believer Asynchronous GP Conditioning):**
   Update `select_top_variants()` in `workers/math_worker.py` to condition the Heteroscedastic GP on provisional `SpeculativeEvidence` from `ledger/database.py`, guaranteeing zero redundant re-sampling while physical plates are in-flight.
4. **Seam 4 (4PL Dose-Response $\text{IC}_{50}/K_d$ & Kinetic Curve Ingestion + Synthesis Cost Gate):**
   Add `fit_dose_response_4pl()` and `fit_michaelis_menten_kinetics()` in `scripts/ingest_plate_data.py`, `estimate_wetlab_synthesis_cost()` in `controller/gatekeeper.py`, and `POST /api/campaign/ingest-curves` in `server.py`.

