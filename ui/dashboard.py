"""
ui/dashboard.py — Real-Time Telemetry & Economics Dashboard.

Purpose:
    Enterprise Gradio frontend for VeritasBio.
    Provides non-blocking telemetry monitoring, budget controls, live audit trails,
    empirical benchmark Pareto visualization, and protein fitness progression tracking.

Constraints (ARCHITECTURE.md §7.2, §7.3; INNOVATION_ROADMAP Phase 5):
    1. UI Asynchrony: Uses generator `yield` streaming. Never blocks the main thread
       waiting on Modal or model forward passes.
    2. Read-Only Ledger View: Displays immutable data from the SQLite database;
       never bypasses `authorize_spend()`.
    3. High-Contrast Enterprise Aesthetics: Clean modern canvas (#f8fafc), crisp white
       panels (#ffffff), high-contrast slate typography (#0f172a), high-tech dark terminal
       event stream, and interactive Pareto curve visualization.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

import gradio as gr

from controller.orchestrator import Orchestrator, OrchestratorConfig, OrchestratorEvent
from ledger import database as ledger_db

logger = logging.getLogger(__name__)

BENCHMARK_PLOT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "benchmark_results.png"


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
    for rank, ev in enumerate(sorted_ev[:35], 1):
        seq = ev.get("variant_sequence", "")
        # Truncate sequence display if full string for readability
        display_seq = f"{seq[:14]}...{seq[-8:]}" if len(seq) > 26 else seq
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
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&family=Outfit:wght@600;700;800;900&display=swap');

    :root {
        --canvas-bg: #f8fafc;
        --card-bg: #ffffff;
        --border-subtle: #e2e8f0;
        --border-accent: #cbd5e1;
        --text-primary: #0f172a;
        --text-secondary: #334155;
        --text-muted: #64748b;
        --brand-blue: #0284c7;
        --brand-purple: #7c3aed;
        --brand-emerald: #059669;
        --brand-amber: #d97706;
    }

    body, .gradio-container {
        background-color: var(--canvas-bg) !important;
        font-family: 'Inter', -apple-system, sans-serif !important;
        color: var(--text-primary) !important;
    }

    /* Brand Header */
    .brand-container {
        display: flex;
        align-items: center;
        gap: 1.25rem;
        padding: 0.75rem 0 1.25rem 0;
        border-bottom: 2px solid #e2e8f0;
        margin-bottom: 1.5rem;
    }
    .brand-logo-mark {
        width: 48px;
        height: 48px;
        border-radius: 12px;
        background: linear-gradient(135deg, #0284c7, #7c3aed);
        display: flex;
        align-items: center;
        justify-content: center;
        font-family: 'Outfit', sans-serif;
        font-weight: 800;
        font-size: 1.6rem;
        color: #ffffff;
        box-shadow: 0 4px 14px rgba(2, 132, 199, 0.35);
    }
    .brand-title {
        font-family: 'Outfit', sans-serif;
        font-size: 2.1rem;
        font-weight: 800;
        color: #0f172a !important;
        letter-spacing: -0.03em;
        line-height: 1.1;
    }
    .brand-title span {
        color: #0284c7 !important;
    }
    .brand-subtitle {
        color: #475569 !important;
        font-size: 0.95rem;
        font-weight: 500;
        letter-spacing: -0.01em;
        margin-top: 0.2rem;
    }

    /* High-Level Stat Cards */
    .stat-card-wrap {
        background: #ffffff !important;
        border: 1px solid #e2e8f0 !important;
        border-radius: 12px !important;
        box-shadow: 0 1px 4px rgba(0, 0, 0, 0.05) !important;
        padding: 1rem 1.25rem !important;
        transition: transform 0.2s ease, box-shadow 0.2s ease !important;
    }
    .stat-card-wrap:hover {
        box-shadow: 0 6px 16px rgba(0, 0, 0, 0.08) !important;
        transform: translateY(-1px) !important;
    }

    /* High-Impact Launch Button */
    .btn-launch {
        background: linear-gradient(135deg, #0284c7, #6366f1) !important;
        color: #ffffff !important;
        font-family: 'Outfit', sans-serif !important;
        font-weight: 700 !important;
        font-size: 1.05rem !important;
        letter-spacing: 0.01em !important;
        border: none !important;
        border-radius: 10px !important;
        padding: 0.85rem !important;
        box-shadow: 0 4px 14px rgba(2, 132, 199, 0.35) !important;
        transition: all 0.2s ease !important;
        cursor: pointer !important;
    }
    .btn-launch:hover {
        box-shadow: 0 6px 20px rgba(2, 132, 199, 0.5) !important;
        transform: translateY(-2px) !important;
    }

    /* Terminal Console for Events */
    .terminal-container {
        border-radius: 10px;
        overflow: hidden;
        border: 1px solid #1e293b;
        background: #090d16;
        margin-top: 1rem;
    }
    .terminal-header {
        background: #0f172a;
        padding: 8px 14px;
        display: flex;
        align-items: center;
        gap: 6px;
        border-bottom: 1px solid rgba(255, 255, 255, 0.08);
    }
    .terminal-dot {
        width: 10px;
        height: 10px;
        border-radius: 50%;
        display: inline-block;
    }
    .terminal-title {
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.78rem;
        color: #94a3b8;
        margin-left: 8px;
        font-weight: 600;
    }

    /* Event Stream Logbox */
    .monospace-log textarea {
        background: #090d16 !important;
        color: #38bdf8 !important;
        font-family: 'JetBrains Mono', monospace !important;
        font-size: 0.86rem !important;
        line-height: 1.55 !important;
        border: 1px solid #1e293b !important;
        border-radius: 0 0 10px 10px !important;
    }

    /* Form and Panel Inputs */
    label, .label-wrap {
        color: #1e293b !important;
        font-weight: 600 !important;
        font-size: 0.92rem !important;
    }
    input, select, .gr-input {
        background-color: #ffffff !important;
        border: 1px solid #cbd5e1 !important;
        color: #0f172a !important;
        border-radius: 8px !important;
        font-weight: 500 !important;
    }
    input:focus, select:focus {
        border-color: #0284c7 !important;
        box-shadow: 0 0 0 2px rgba(2, 132, 199, 0.2) !important;
    }

    /* Tables & Dataframes */
    table, .table-wrap {
        background: #ffffff !important;
        border: 1px solid #e2e8f0 !important;
        border-radius: 8px !important;
        color: #0f172a !important;
    }
    th {
        background: #f1f5f9 !important;
        color: #0f172a !important;
        font-weight: 700 !important;
        border-bottom: 2px solid #cbd5e1 !important;
    }
    td {
        border-bottom: 1px solid #f1f5f9 !important;
        color: #1e293b !important;
    }
    tr:nth-child(even) td {
        background: #f8fafc !important;
    }
    """

    with gr.Blocks(title="VeritasBio — Cloud Economics Engine") as demo:
        gr.HTML(f"<style>{custom_css}</style>")

        # Top Navigation Header
        with gr.Row():
            with gr.Column(scale=8):
                gr.HTML("""
                <div class="brand-container">
                    <div class="brand-logo-mark">V</div>
                    <div>
                        <div class="brand-title">Veritas<span>Bio</span></div>
                        <div class="brand-subtitle">Deterministic Cloud Economics Engine for Computational Biology (v3.1)</div>
                    </div>
                </div>
                """)
            with gr.Column(scale=4):
                status_box = gr.Textbox(
                    label="Active Engine Status",
                    value="● Ready to initialize",
                    interactive=False,
                )

        # High-Contrast KPI Stat Cards
        with gr.Row():
            with gr.Column(scale=4, elem_classes=["stat-card-wrap"]):
                metric_rem_budget = gr.Number(
                    label="Remaining Capital ($ USD)",
                    value=5.00,
                    precision=4,
                    interactive=False,
                )
            with gr.Column(scale=4, elem_classes=["stat-card-wrap"]):
                metric_spent = gr.Number(
                    label="Cumulative Cloud Spend ($ USD)",
                    value=0.00,
                    precision=4,
                    interactive=False,
                )
            with gr.Column(scale=4, elem_classes=["stat-card-wrap"]):
                metric_top_fitness = gr.Number(
                    label="Top Discovered Fitness Score",
                    value=0.00,
                    precision=4,
                    interactive=False,
                )

        # Main Workspace: Controls vs Tabs
        with gr.Row():
            # Left Panel: Run Configuration & Terminal
            with gr.Column(scale=4):
                with gr.Group():
                    gr.Markdown("### Autonomous Search Controls")
                    input_protein = gr.Textbox(
                        label="Target Protein Accession",
                        value="P42212 (Aequorea victoria avGFP)",
                        interactive=False,
                    )
                    input_budget = gr.Slider(
                        label="Compute Budget Allocation ($ USD)",
                        minimum=0.05,
                        maximum=50.0,
                        step=0.05,
                        value=1.00,
                    )
                    input_iterations = gr.Slider(
                        label="Maximum Active Learning Iterations",
                        minimum=1,
                        maximum=20,
                        step=1,
                        value=5,
                    )
                    input_batch_size = gr.Slider(
                        label="Micro-Batch Size (Variants per Loop)",
                        minimum=10,
                        maximum=100,
                        step=10,
                        value=50,
                    )
                    input_model = gr.Dropdown(
                        label="Foundation Inference Model",
                        choices=["mock_esm", "esm2_650m", "esm2_150m"],
                        value="mock_esm",
                    )
                    check_live_llm = gr.Checkbox(
                        label="Enable Live Meta AI Strategic Planner (Calls Meta API)",
                        value=False,
                        interactive=True,
                    )

                    btn_start = gr.Button("Launch Active Learning Loop", variant="primary", elem_classes=["btn-launch"])
                    btn_refresh = gr.Button("Sync Ledger State")

                # Terminal Event Stream Box
                with gr.Group():
                    gr.HTML("""
                    <div class="terminal-container">
                        <div class="terminal-header">
                            <span class="terminal-dot" style="background:#ef4444;"></span>
                            <span class="terminal-dot" style="background:#f59e0b;"></span>
                            <span class="terminal-dot" style="background:#10b981;"></span>
                            <span class="terminal-title">veritasbio-kernel :: live-event-stream.log</span>
                        </div>
                    </div>
                    """)
                    log_output = gr.Textbox(
                        label="Real-Time Event Stream",
                        lines=10,
                        max_lines=14,
                        interactive=False,
                        elem_classes=["monospace-log"],
                    )

            # Right Panel: Telemetry, Tabs, and Pareto Plots
            with gr.Column(scale=8):
                with gr.Tabs():
                    with gr.TabItem("Top Discovered Variants"):
                        table_variants = gr.Dataframe(
                            headers=["Rank", "Variant Sequence", "Fitness Score", "Job ID", "Timestamp"],
                            datatype=["number", "str", "str", "str", "str"],
                            interactive=False,
                            wrap=True,
                        )

                    with gr.TabItem("Ledger Audit Trail (Spend & Authorization)"):
                        table_jobs = gr.Dataframe(
                            headers=["Job ID", "Model", "Variant Count", "Estimated Cost", "Status"],
                            datatype=["str", "str", "number", "str", "str"],
                            interactive=False,
                            wrap=True,
                        )

                    with gr.TabItem("Pareto Frontier & Cloud Savings"):
                        plot_preview = gr.Image(
                            value=str(BENCHMARK_PLOT_PATH) if BENCHMARK_PLOT_PATH.is_file() else None,
                            label="Empirical Benchmark & Cloud Valuation (Sarkisyan 2016 Nature avGFP)",
                            type="filepath",
                            interactive=False,
                        )

        # ---------------------------------------------------------------------
        # Event Handlers & Streaming Generator
        # ---------------------------------------------------------------------
        def on_refresh() -> Tuple[float, float, float, str, List[List[Any]], List[List[Any]], Optional[str]]:
            rem = ledger_db.get_remaining_budget(db_path=db_path)
            jobs = ledger_db.list_jobs(db_path=db_path)
            evidence = ledger_db.get_evidence(db_path=db_path)

            spent = sum(float(j["estimated_cost"]) for j in jobs if j["status"] in ("APPROVED", "COMPLETED", "FAILED"))
            top_fit = max([float(e["computed_score"]) for e in evidence], default=0.0)
            plot_file = str(BENCHMARK_PLOT_PATH) if BENCHMARK_PLOT_PATH.is_file() else None

            return (
                round(rem, 4),
                round(spent, 4),
                round(top_fit, 4),
                f"Ledger synced: {len(jobs)} jobs, {len(evidence)} evidence entries.",
                format_evidence_table(evidence),
                format_jobs_table(jobs),
                plot_file,
            )

        def stream_orchestrator(
            budget_val: float,
            iter_cap: int,
            batch_sz: int,
            target_mdl: str,
            enable_live_llm_val: bool = False,
        ) -> Iterator[Tuple[float, float, float, str, str, List[List[Any]], List[List[Any]], Optional[str]]]:
            """Streams live events using Python generator yields."""
            config = OrchestratorConfig(
                target_accession="P42212",
                total_budget_usd=budget_val,
                max_iterations=iter_cap,
                batch_size=batch_sz,
                target_model=target_mdl,
                db_path=db_path,
                use_mock_llm=not enable_live_llm_val,
                enable_live_llm=enable_live_llm_val,
            )

            orchestrator = Orchestrator(config)
            accumulated_logs: List[str] = []
            plot_file = str(BENCHMARK_PLOT_PATH) if BENCHMARK_PLOT_PATH.is_file() else None

            for event in orchestrator.run_stream():
                accumulated_logs.append(f"[{event.step.upper():12s}] {event.message}")
                log_str = "\n".join(accumulated_logs[-14:])

                jobs = ledger_db.list_jobs(db_path=db_path)
                evidence = ledger_db.get_evidence(db_path=db_path)
                top_fit = max([float(e["computed_score"]) for e in evidence], default=0.0)

                yield (
                    round(event.remaining_budget_usd, 4),
                    round(event.cumulative_spent_usd, 4),
                    round(top_fit, 4),
                    f"● Round {event.iteration} — {event.step.upper()}",
                    log_str,
                    format_evidence_table(evidence),
                    format_jobs_table(jobs),
                    plot_file,
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
                plot_preview,
            ],
        )

        btn_start.click(
            fn=stream_orchestrator,
            inputs=[input_budget, input_iterations, input_batch_size, input_model, check_live_llm],
            outputs=[
                metric_rem_budget,
                metric_spent,
                metric_top_fitness,
                status_box,
                log_output,
                table_variants,
                table_jobs,
                plot_preview,
            ],
        )

    return demo


if __name__ == "__main__":
    app_demo = create_dashboard()
    app_demo.launch(server_name="127.0.0.1", server_port=7860, inbrowser=False)
