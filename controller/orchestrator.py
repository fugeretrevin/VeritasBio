"""
controller/orchestrator.py — DOMAIN: Brains & Security.

Purpose:
    Runs the unidirectional 6-step active-learning loop (ARCHITECTURE.md §4):
    Initialization -> Strategy -> Selection -> Authorization -> Execution -> Update,
    repeating until the ledger budget is exhausted.

Constraints:
    - Deterministic, bounded control flow. NOT an open-ended AI agent loop:
      the LLM never chooses which step/tool runs next.
    - Execution (Step 5) is only reachable after authorize_spend() returns True.

Status: scaffold only — not implemented yet.
"""
