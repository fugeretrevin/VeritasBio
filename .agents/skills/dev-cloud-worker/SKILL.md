---
name: dev-cloud-worker
description: >-
  Use this skill when developing, refactoring, or testing cloud GPU workers,
  foundation model scoring (ESM-2), backbone diffusion (RFdiffusion), sequence inverse folding
  (ProteinMPNN), or high-fidelity mock implementations in workers/.
---

# Cloud Worker Development Skill

This skill governs the development of high-throughput computational biology workers in `workers/`.

## Core Invariants
1. **Never Bypass the Spend Gate:**
   - Every worker execution MUST require a `job_id` and be called ONLY after `ledger_db.authorize_spend()` returns `True`.
2. **Dual Execution (Modal + High-Fidelity Mock):**
   - Modal serverless app decorator pattern (`_modal_function`) must be used for cloud deployment.
   - A deterministic, high-fidelity local physics/mock fallback MUST always be implemented so unit tests and local runs never hard-crash when Modal API keys are absent.
3. **Container Amortization:**
   - Configure Modal apps with `scaledown_window` (e.g. 120s) to reuse warm containers across active learning iterations.
4. **Biophysical Sanity:**
   - Always run developability filtering (via `analyze_biologics_developability()`) on designed or mutated sequences.
   - Enforce zero unpaired cysteines on de novo designs.

## File Manifest
- [workers/esm_worker.py](file:///c:/Users/stoof/source/repos/VeritasBio/workers/esm_worker.py): ESM-2 650M / 150M zero-shot masked marginal scoring.
- [workers/diffusion_worker.py](file:///c:/Users/stoof/source/repos/VeritasBio/workers/diffusion_worker.py): De novo 3-helix bundle geometric backbone generator & ProteinMPNN inverse folding.
- [workers/mock_workers.py](file:///c:/Users/stoof/source/repos/VeritasBio/workers/mock_workers.py): High-fidelity biological landscape simulators with non-linear epistasis.
- [workers/spatial_filter.py](file:///c:/Users/stoof/source/repos/VeritasBio/workers/spatial_filter.py): 3D PDB coordinate extractor and contact sphere calculator.

## Verification Checklist
```powershell
python -m unittest tests/test_diffusion_worker.py
python -m unittest tests/test_esm_worker.py
```
Ensure:
- Budget authorization rejection correctly halts worker execution.
- PDB coordinates and secondary structural metrics match expected ranges.
