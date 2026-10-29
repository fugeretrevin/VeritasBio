---
name: dev-orchestrator
description: >-
  Use this skill to orchestrate complex feature developments across PanoptesBio.
  It decomposes user requests into non-overlapping parallel tracks, coordinates subagents,
  runs automated verification gates, and synthesizes approval summaries for the user.
---

# Development Orchestrator Skill

The **Development Orchestrator** is the meta-level coordinator for PanoptesBio engineering tasks. It enforces modularity, grounding, parallel development, and human approval gates.

## Step-by-Step Orchestration Workflow

### Step 1: Decompose into Disjoint Tracks
Analyze the user request and divide the implementation across the 4 orthogonal domains:
- **Track 1: Cloud & Science Workers** (`workers/`)
- **Track 2: Brains & Economics Gatekeeper** (`controller/`, `ledger/`)
- **Track 3: UI & 3D Visualization** (`ui/`, `main.py`)
- **Track 4: Benchmarks & Wet-Lab Compliance** (`scripts/`, `tests/`)

### Step 2: Establish the Interface Contract First
Before modifying backend logic or writing UI:
1. Verify existing schemas in [gatekeeper.py](file:///c:/Users/stoof/source/repos/PanoptesBio/controller/gatekeeper.py).
2. If new endpoints or target models are needed, add them to `ALLOWED_TARGET_MODELS` and define strongly typed Pydantic models with field validators.

### Step 3: Execute Disjoint Tasks
Implement changes within the track boundaries. Subagents working on Track 1 must not modify Track 3, avoiding race conditions and merge conflicts.

### Step 4: Run Grounded Verification Gates
Execute the workspace verification script:
```powershell
python .agents/scripts/verify_workspace.py
```
Or execute the full unit test suite:
```powershell
python -m unittest discover tests
```
Ensure:
- 100% test pass rate.
- Cryptographic Merkle chain validation passes (`ledger_db.verify_audit_integrity()`).
- Zero hardcoded mock bypasses in production Gatekeeper paths.

### Step 5: Synthesize Human Approval Checkpoint
Present a concise report to the user:
- What was built and why.
- Verification test results (number of tests passed, execution time).
- Empirical biophysical proof (e.g. developability score, radius of gyration, predicted pLDDT, or benchmark correlation $\rho$).
- Proposed commit message, requesting explicit permission before staging/committing.
