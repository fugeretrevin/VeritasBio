"""
main.py — Application entry point for PanoptesBio.

Purpose:
    Initializes the SQLite ledger, starts the Gradio telemetry dashboard,
    or runs headless active learning loops via CLI arguments.

Usage:
    # Launch interactive Gradio web dashboard:
    python main.py

    # Run headless autonomous CLI loop:
    python main.py --cli --budget 2.00 --iterations 5 --model mock_esm

Invariants (ARCHITECTURE.md §3, §7; INNOVATION_ROADMAP Phase 5):
    - Paid compute always flows through Gatekeeper and `authorize_spend()`.
    - Fail-closed execution.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from controller.orchestrator import Orchestrator, OrchestratorConfig, OrchestratorEvent
from ledger import database as ledger_db
from ui.dashboard import create_dashboard

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("panoptesbio")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="PanoptesBio: Deterministic Cloud Economics Engine for Computational Biology",
    )
    parser.add_argument(
        "--cli",
        action="store_true",
        help="Run in headless CLI mode instead of launching Gradio UI",
    )
    parser.add_argument(
        "--budget",
        type=float,
        default=5.00,
        help="Total run budget in USD (default: $5.00)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="Maximum loop iteration cap (default: 5)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Variants per micro-batch (default: 50)",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="mock_esm",
        choices=["mock_esm", "esm2_650m", "esm2_150m"],
        help="Target foundation scoring model (default: mock_esm)",
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Custom path to SQLite ledger file",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=7860,
        help="Port for Gradio web UI (default: 7860)",
    )
    parser.add_argument(
        "--live-llm",
        action="store_true",
        default=False,
        help="Explicitly enable live Meta Model API (default: False, prevents credit burn)",
    )

    args = parser.parse_args()

    # Step 1: Initialize Ledger Database
    db_file = ledger_db.init_db(args.db)
    logger.info("Ledger initialized at %s", db_file)

    if args.cli:
        logger.info(
            "Starting Headless Orchestrator (Budget: $%.2f, Max Iterations: %d, Model: %s, Live LLM: %s)...",
            args.budget, args.iterations, args.model, args.live_llm,
        )
        config = OrchestratorConfig(
            target_accession="P42212",
            total_budget_usd=args.budget,
            max_iterations=args.iterations,
            batch_size=args.batch_size,
            target_model=args.model,
            db_path=args.db,
            use_mock_llm=not args.live_llm,
            enable_live_llm=args.live_llm,
        )
        orchestrator = Orchestrator(config)

        def print_event(ev: OrchestratorEvent) -> None:
            print(f"[{ev.step.upper():12s}] Round {ev.iteration:02d} | Left: ${ev.remaining_budget_usd:.4f} | {ev.message}")

        summary = orchestrator.run(on_event=print_event)
        print("\n" + "=" * 60)
        print("PANOPTESBIO RUN SUMMARY")
        print("=" * 60)
        print(f"Target:               {summary.target_protein} ({summary.target_accession})")
        print(f"Rounds Completed:     {summary.iterations_completed}")
        print(f"Total Spent:          ${summary.final_spent_usd:.6f}")
        print(f"Remaining Budget:     ${summary.final_remaining_usd:.6f}")
        print(f"Total Variants Evaluated: {summary.total_variants_evaluated}")
        print(f"Highest Fitness Score:    {summary.highest_fitness_score:.4f}")
        print("=" * 60)
    else:
        logger.info("Launching PanoptesBio Gradio Dashboard on port %d...", args.port)
        import gradio as gr
        demo = create_dashboard(db_path=args.db)
        demo.launch(
            server_name="127.0.0.1",
            server_port=args.port,
            inbrowser=False,
            css=getattr(demo, "custom_css", None),
            theme=gr.themes.Default(),
        )


if __name__ == "__main__":
    main()
