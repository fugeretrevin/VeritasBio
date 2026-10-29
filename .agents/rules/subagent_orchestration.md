# Subagent Orchestration & Parallel Decomposition Rules

## 1. Principle of Disjoint Decomposition
To parallelize development across agents or subagents without merge conflicts or regressions, tasks must be strictly mapped to **orthogonal, non-overlapping workstreams**:

| Workstream | Scope / Directories | Primary Artifacts | Concurrency Rule |
| :--- | :--- | :--- | :--- |
| **Track 1: Cloud & Muscle** | `workers/` | ESM, Diffusion, Spatial contact | Must not modify `controller/` or `ui/` |
| **Track 2: Brains & Economics** | `controller/`, `ledger/` | Gatekeeper, Prompts, SQLite schema | Must not touch frontend or robotics protocols |
| **Track 3: UI & Visualization** | `ui/`, `main.py` | Gradio blocks, 3Dmol WebGL viewers | Consumes worker APIs; never calls Modal directly |
| **Track 4: Benchmarks & CRO** | `scripts/`, `tests/` | Benchmark scripts, unit tests, Twist CSV | Read-only access to engines; emits test reports |

## 2. Orchestration Protocol
When tackling multi-component features or broad architectural expansions:
1. **Step 1 — Partitioning:** Decompose the user request into disjoint tasks matching the 4 tracks above.
2. **Step 2 — Contract First:** Define shared schemas in `controller/gatekeeper.py` or interface types before coding implementations.
3. **Step 3 — Parallel Execution:** Work on Track 1, Track 2, Track 3, or Track 4 independently.
4. **Step 4 — Verification Gate:** Run `tests/` to verify that all components integrate cleanly.
5. **Step 5 — Approval Checkpoint:** Request user sign-off with clear verification output.

## 3. Grounding Verification
Subagents must never assume an API or function signature exists without checking the code using `view_file` or `grep_search`. Always inspect parameter names, return types, and default values.
