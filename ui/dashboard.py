"""
ui/dashboard.py — DOMAIN: Frontend.

Purpose:
    Gradio interface to set target protein + budget and observe progress
    (budget remaining, job statuses, evidence/fitness over time).

Constraints (ARCHITECTURE.md §7.2, §7.3):
    - Never block the main thread on a Modal job; use `yield`/background tasks
      to poll the SQLite ledger.
    - Never invokes workers directly; read-only view of the ledger plus a
      "start run" hand-off to the Orchestrator.

Status: scaffold only — not implemented yet.
"""
