"""
main.py — Application entry point.

Purpose:
    Initializes the SQLite ledger (ledger.database.init_db), starts the Gradio UI
    (ui.dashboard), and wires the UI to the Orchestrator (controller.orchestrator).

Constraints (ARCHITECTURE.md §3, §7):
    - Does not invoke any worker directly. All spend flows through
      ledger.database.authorize_spend() via the Gatekeeper/Orchestrator.

Status: scaffold only — not implemented yet.
"""
