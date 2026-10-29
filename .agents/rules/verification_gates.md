# Verification Gates & Grounding Mandate for Development Agents

## 1. Grounding Mandate (Anti-Hallucination)
1. **Empirical Grounding Only:**
   - Agents must never invent biophysical constants, protein sequences, residue numbers, or benchmark scores.
   - Every residue position must be cross-referenced with local reference files in `workers/data/` (e.g., `1EMA.pdb`, `2YBR.pdb`, `8EH0.pdb`, `P42212.json`, `SCFV_9004G_CN2.json`, `sarkisyan_gfp_benchmark.json`).
2. **Deterministic Schemas:**
   - All payloads passing between brains, quant engine, and cloud muscle must be strictly validated with Pydantic schemas in `controller/gatekeeper.py`.
   - Fail-closed: invalid schemas or unexpected fields must trigger immediate rejection.

## 2. Mandatory Verification Gates
Before presenting any work as completed or asking the user for commit approval, the agent **MUST** execute and pass:
1. **Targeted Unit Tests:**
   - Run the specific test file covering the modified module (e.g. `python -m unittest tests/test_XYZ.py`).
2. **Full Repository Regression Suite:**
   - Run `python -m unittest discover tests`.
   - All tests must pass (100% green). Zero regressions allowed.
3. **FDA 21 CFR Part 11 Audit Chain Check:**
   - If any change touches `ledger/`, `workers/wetlab_worker.py`, or `workers/lims_worker.py`, run `ledger_db.verify_audit_integrity()` to guarantee that the SHA-256 Merkle hash chain and HMAC signatures are valid.

## 3. Human Approval Invariants
1. **No Silent Git Pushes:** Never execute `git push` without explicit user permission.
2. **Pre-Commit Review:** Always provide a clear, concise summary of modified files, new tests, and empirical metrics before asking the user to commit.
3. **Fail-Closed Spend Authorization:** Even in mock testing, cloud compute spend must follow `ledger_db.authorize_spend()`. Never bypass the ledger.
