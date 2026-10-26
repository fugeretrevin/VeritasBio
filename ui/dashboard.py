"""
ui/dashboard.py — DOMAIN: Frontend.

Purpose:
    Real-time Gradio dashboard for VeritasBio.
    Provides non-blocking telemetry monitoring, budget controls, live audit trails,
    and protein fitness progression tracking.

Constraints (ARCHITECTURE.md §7.2, §7.3; INNOVATION_ROADMAP Phase 5):
    1. UI Asynchrony: Uses generator `yield` streaming. Never blocks the main thread
       waiting on Modal or model forward passes.
    2. Read-Only Ledger View: Displays immutable data from the SQLite database;
       never bypasses `authorize_spend()`.
    3. Visual Excellence: Sleek dark theme with glowing metric cards, live status pulses,
       and real-time audit tables.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import gradio as gr

from controller.orchestrator import Orchestrator, OrchestratorConfig, OrchestratorEvent
from ledger import database as ledger_db

logger = logging.getLogger(__name__)


def format_jobs_table(jobs: List[Dict[str, Any]]) -> List[List[Any]]:
    """Formats jobs for the Gradio Dataframe."""
    rows = []
    for j in reversed(jobs):
        rows.append([
            j.get("job_id", ""),
            j.get("target_model", ""),
            j.get("variant_count", 0),
            f"${float(j.get('estimated_cost', 0.0)):.6f}",
            j.get("status", ""),
        ])
    return rows


def format_evidence_table(evidence: List[Dict[str, Any]]) -> List[List[Any]]:
    """Formats top evidence for the Gradio Dataframe."""
    sorted_ev = sorted(evidence, key=lambda e: float(e["computed_score"]), reverse=True)
    rows = []
    for rank, ev in enumerate(sorted_ev[:25], 1):
        seq = ev.get("variant_sequence", "")
        # Truncate sequence display if full string
        display_seq = f"{seq[:16]}...{seq[-8:]}" if len(seq) > 28 else seq
        rows.append([
            rank,
            display_seq,
            f"{float(ev.get('computed_score', 0.0)):.4f}",
            ev.get("job_id", ""),
            ev.get("timestamp", ""),
        ])
    return rows


def create_dashboard(db_path: Optional[str] = None) -> gr.Blocks:
    """Builds and returns the complete Gradio Blocks application."""

    custom_css = """
    .gradio-container {
        background-color: #07090e !important;
        font-family: 'Inter', sans-serif !important;
    }
    .metric-card {
        background: rgba(18, 23, 37, 0.75) !important;
        border: 1px solid rgba(255, 255, 255, 0.08) !important;
        border-radius: 12px !important;
        padding: 1.25rem !important;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.4) !important;
    }
    .brand-title {
        font-family: 'Outfit', sans-serif;
        font-size: 1.8rem;
        font-weight: 700;
        color: #ffffff;
        letter-spacing: -0.02em;
    }
    .brand-sub {
        color: #94a3b8;
        font-size: 0.95rem;
    }
    """

    with gr.Blocks(title="VeritasBio — Cloud Economics Engine") as demo:
        gr.HTML(f"<style>{custom_css}</style>")

        # Header Section
        with gr.Row():
            with gr.Column(scale=8):
                gr.HTML("""
                <div style="display:flex; align-items:center; gap:1rem; margin-bottom: 0.5rem;">
                    <div style="width:40px; height:40px; border-radius:10px; background:linear-gradient(135deg, #00f0ff, #8b5cf6); display:flex; align-items:center; justify-content:center; font-weight:800; font-size:1.3rem; color:#05070c;">V</div>
                    <div>
                        <div class="brand-title">VeritasBio</div>
                        <div class="brand-sub">Deterministic Cloud Economics Engine for Computational Biology (v3.1)</div>
                    </div>
                </div>
                """)
            with gr.Column(scale=4):
                status_box = gr.Textbox(
                    label="Engine State",
                    value="Ready to initialize",
                    interactive=False,
                )

        # High-Level Metric Tiles
        with gr.Row():
            with gr.Column(scale=4):
                metric_rem_budget = gr.Number(label="Remaining Budget ($)", value=5.00, precision=4)
            with gr.Column(scale=4):
                metric_spent = gr.Number(label="Cumulative Spent ($)", value=0.00, precision=4)
            with gr.Column(scale=4):
                metric_top_fitness = gr.Number(label="Top Discovered Fitness", value=0.00, precision=4)

        # Main Workspace
        with gr.Row():
            # Left Panel: Controls & Telemetry
            with gr.Column(scale=4):
                gr.Markdown("### Run Configuration")
                input_protein = gr.Textbox(label="Target Accession", value="P42212 (avGFP)", interactive=False)
                input_budget = gr.Slider(label="Budget (USD)", minimum=0.05, maximum=50.0, step=0.05, value=1.00)
                input_iterations = gr.Slider(label="Max Iterations (Cap)", minimum=1, maximum=20, step=1, value=5)
                input_batch_size = gr.Slider(label="Micro-Batch Size", minimum=10, maximum=100, step=10, value=50)
                input_model = gr.Dropdown(
                    label="Target Foundation Model",
                    choices=["mock_esm", "esm2_650m", "esm2_150m"],
                    value="mock_esm",
                )

                btn_start = gr.Button("Launch Optimization Loop", variant="primary")
                btn_refresh = gr.Button("Refresh Ledger Telemetry")

                log_output = gr.Textbox(
                    label="Real-Time Event Stream",
                    lines=8,
                    max_lines=12,
                    interactive=False,
                )

            # Right Panel: Tables & Telemetry
            with gr.Column(scale=8):
                with gr.Tabs():
                    with gr.TabItem("Top Variants (Evidence)"):
                        table_variants = gr.Dataframe(
                            headers=["Rank", "Variant Sequence", "Fitness Score", "Job ID", "Timestamp"],
                            datatype=["number", "str", "str", "str", "str"],
                            interactive=False,
                        )

                    with gr.TabItem("Ledger Audit Trail (Jobs)"):
                        table_jobs = gr.Dataframe(
                            headers=["Job ID", "Model", "Variant Count", "Estimated Cost", "Status"],
                            datatype=["str", "str", "number", "str", "str"],
                            interactive=False,
                        )

        # ---------------------------------------------------------------------
        # Event Handlers & Streaming Generator
        # ---------------------------------------------------------------------
        def on_refresh() -> Tuple[float, float, float, str, List[List[Any]], List[List[Any]]]:
            rem = ledger_db.get_remaining_budget(db_path=db_path)
            jobs = ledger_db.list_jobs(db_path=db_path)
            evidence = ledger_db.get_evidence(db_path=db_path)

            spent = sum(float(j["estimated_cost"]) for j in jobs if j["status"] in ("APPROVED", "COMPLETED", "FAILED"))
            top_fit = max([float(e["computed_score"]) for e in evidence], default=0.0)

            return (
                rem,
                spent,
                top_fit,
                f"Ledger synced: {len(jobs)} jobs, {len(evidence)} evidence entries.",
                format_evidence_table(evidence),
                format_jobs_table(jobs),
            )

        def stream_orchestrator(
            budget_val: float,
            iter_cap: int,
            batch_sz: int,
            target_mdl: str,
        ) -> Iterator[Tuple[float, float, float, str, str, List[List[Any]], List[List[Any]]]]:
            """Streams live events using Python generator yields."""
            config = OrchestratorConfig(
                target_accession="P42212",
                total_budget_usd=budget_val,
                max_iterations=iter_cap,
                batch_size=batch_sz,
                target_model=target_mdl,
                db_path=db_path,
                use_mock_llm=(target_mdl == "mock_esm"),
            )

            orchestrator = Orchestrator(config)
            accumulated_logs: List[str] = []

            for event in orchestrator.run_stream():
                accumulated_logs.append(f"[{event.step.upper()}] {event.message}")
                log_str = "\n".join(accumulated_logs[-12:])

                jobs = ledger_db.list_jobs(db_path=db_path)
                evidence = ledger_db.get_evidence(db_path=db_path)
                top_fit = max([float(e["computed_score"]) for e in evidence], default=0.0)

                yield (
                    round(event.remaining_budget_usd, 4),
                    round(event.cumulative_spent_usd, 4),
                    round(top_fit, 4),
                    f"Round {event.iteration} — {event.step.capitalize()}",
                    log_str,
                    format_evidence_table(evidence),
                    format_jobs_table(jobs),
                )

        btn_refresh.click(
            fn=on_refresh,
            outputs=[
                metric_rem_budget,
                metric_spent,
                metric_top_fitness,
                status_box,
                table_variants,
                table_jobs,
            ],
        )

        btn_start.click(
            fn=stream_orchestrator,
            inputs=[input_budget, input_iterations, input_batch_size, input_model],
            outputs=[
                metric_rem_budget,
                metric_spent,
                metric_top_fitness,
                status_box,
                log_output,
                table_variants,
                table_jobs,
            ],
        )

    return demo


if __name__ == "__main__":
    app_demo = create_dashboard()
    app_demo.launch(server_name="127.0.0.1", server_port=7860, inbrowser=False)
