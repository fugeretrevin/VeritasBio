"""
scripts/run_full_dollar_experiment.py — Full 7-Layer Production Run with Live LLM & Modal A10G.

Executes a complete, live $1.00 budget experiment through all 7 architectural layers:
  Layer 1: Biological Prior Ingestion (UniProt RAG)
  Layer 2: Strategic Planning (Live Meta Model API in strict JSON mode)
  Layer 3: Quant Engine Active Learning (Physicochemical featurization + GP Expected Improvement)
  Layer 4: Deterministic Financial Gatekeeper (Fail-closed schema validation + cold/warm pricing)
  Layer 5: Immutable Financial Ledger (SQLite micro-USD atomic lock)
  Layer 6: Cloud Muscle (Modal serverless A10G GPU inference)
  Layer 7: Evidence Ingestion & Feedback + Wet-Lab Translation (Twist CSV, Opentrons OT-2, FDA SHA-256 certificate)
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from dotenv import load_dotenv
load_dotenv(ROOT_DIR / ".env")

from controller.orchestrator import Orchestrator, OrchestratorConfig, OrchestratorEvent
from ledger import database as ledger_db
from scripts.evaluate import evaluate_run, print_report, sequence_to_mutation_str
from workers.rag_worker import fetch_protein_metadata
from workers.wetlab_worker import (
    generate_cryptographic_audit_certificate,
    generate_opentrons_protocol,
    generate_twist_order_csv,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("veritasbio.dollar_experiment")


import argparse

def run_dollar_experiment(
    target_accession: str = "TRPB_TM9D8",
    budget_usd: float = 1.00,
    rounds: int = 4,
    batch_size: int = 25,
    target_model: str = "esm2_150m",
    enable_live_llm: bool = True,
    db_path: Optional[str] = None,
) -> None:
    meta = fetch_protein_metadata(target_accession, use_cache=True)
    db_slug = "".join(c for c in target_accession.lower() if c.isalnum() or c == "_")
    db_file = Path(db_path) if db_path else ROOT_DIR / "ledger" / f"dollar_live_{db_slug}.db"
    if db_file.exists():
        try:
            db_file.unlink()
        except Exception:
            pass

    ledger_db.init_db(db_file)

    print("\n" + "=" * 70)
    print(" VERITASBIO FULL 7-LAYER LIVE PRODUCTION EXPERIMENT ($1.00 BUDGET)")
    print("=" * 70)
    print(f"Target Protein:           {meta.protein_name} ({target_accession})")
    print(f"Total Budget:             ${budget_usd:.2f} USD")
    print(f"Max Rounds:               {rounds}")
    print(f"Variants Per Batch:       {batch_size}")
    print(f"Cloud Scoring Engine:     Modal Serverless A10G GPU ({target_model})")
    print(f"Strategic Planner:        {'Meta Model API (LLAMA_API_KEY)' if enable_live_llm else 'Mock Strategic Planner'}")
    print(f"Ledger File:              {db_file}")
    print("=" * 70 + "\n")

    config = OrchestratorConfig(
        target_accession=target_accession,
        total_budget_usd=budget_usd,
        max_iterations=rounds,
        batch_size=batch_size,
        target_model=target_model,
        db_path=str(db_file),
        use_mock_llm=not enable_live_llm,
        enable_live_llm=enable_live_llm,
        enable_dynamic_batch_sizing=True,
    )

    orchestrator = Orchestrator(config)

    def on_event(ev: OrchestratorEvent) -> None:
        step_tag = f"LAYER: {ev.step.upper()}"
        print(f"[{step_tag:18s}] Round {ev.iteration:02d} | Remaining: ${ev.remaining_budget_usd:.6f} | {ev.message}")

    summary = orchestrator.run(on_event=on_event)

    print("\n" + "=" * 70)
    print(" ACTIVE LEARNING LOOP COMPLETE — RUN SUMMARY")
    print("=" * 70)
    print(f"Target Protein:           {meta.protein_name} ({target_accession})")
    print(f"Iterations Completed:     {summary.iterations_completed}")
    print(f"Termination Reason:       {summary.stop_reason}")
    print(f"Total Compute Spent:      ${summary.final_spent_usd:.6f}")
    print(f"Remaining Unspent Funds:  ${summary.final_remaining_usd:.6f}")
    print(f"Total Variants Evaluated: {summary.total_variants_evaluated}")
    print(f"Highest Model Score:      {summary.highest_fitness_score:.4f}")
    print(f"Top Sequence Discovered:  {summary.top_variant_sequence[:35]}...")
    print("=" * 70 + "\n")

    # -----------------------------------------------------------------
    # Benchmark Evaluation against Ground Truth (TrpB or GFP)
    # -----------------------------------------------------------------
    print(">>> Running Benchmark Evaluation Against Empirical Ground Truth...")
    bench_csv = None
    if "trpb" in target_accession.lower():
        bench_csv = str(ROOT_DIR / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv")
    elif "p42212" in target_accession.lower() or "gfp" in target_accession.lower():
        bench_csv = str(ROOT_DIR / "scripts" / "data" / "GFP_AEQVI_Sarkisyan_2016.csv")

    metrics = evaluate_run(
        db_path=str(db_file),
        benchmark_path=bench_csv,
        target_name=meta.protein_name,
        wt_sequence=meta.sequence,
    )
    print_report(metrics)

    # -----------------------------------------------------------------
    # Layer 7: Generate Enterprise Wet-Lab & Compliance Deliverables
    # -----------------------------------------------------------------
    print(">>> Generating Layer 7 Physical Wet-Lab Translation Deliverables...")
    all_evidence = ledger_db.get_evidence(db_path=str(db_file))
    wt_seq = meta.sequence

    # Format top variants with mutation descriptions
    top_candidates = []
    for ev in sorted(all_evidence, key=lambda e: float(e["computed_score"]), reverse=True)[:10]:
        seq = ev.get("variant_sequence", "")
        mut_str = sequence_to_mutation_str(wt_seq, seq)
        top_candidates.append({
            "variant_sequence": seq,
            "mutations": mut_str,
            "computed_score": float(ev["computed_score"]),
        })

    # 1. Twist Bioscience / IDT Order CSV
    twist_csv = generate_twist_order_csv(top_candidates, protein_name=meta.protein_name)
    twist_path = ROOT_DIR / "ledger" / f"twist_idt_{db_slug}_synthesis_order.csv"
    twist_path.write_text(twist_csv, encoding="utf-8")
    print(f" [+] Twist/IDT Synthesis Order CSV written to: {twist_path}")

    # 2. Opentrons OT-2 Pipetting Script
    var_names = [f"{meta.protein_name[:4].upper()}_{c['mutations']}" for c in top_candidates]
    opentrons_script = generate_opentrons_protocol(var_names, target_protein=meta.protein_name)
    opentrons_path = ROOT_DIR / "ledger" / f"opentrons_ot2_{db_slug}_protocol.py"
    opentrons_path.write_text(opentrons_script, encoding="utf-8")
    print(f" [+] Opentrons OT-2 Robotics Protocol written to: {opentrons_path}")

    # 3. FDA 21 CFR Part 11 SHA-256 Audit Certificate
    cert_dict = generate_cryptographic_audit_certificate(
        db_path=Path(db_file),
        target_protein=f"{meta.protein_name} ({target_accession})",
    )
    cert_path = ROOT_DIR / "ledger" / f"fda_21cfr11_{db_slug}_certificate.json"
    cert_path.write_text(json.dumps(cert_dict, indent=2), encoding="utf-8")
    print(f" [+] FDA 21 CFR Part 11 SHA-256 Certificate written to: {cert_path}")
    print("\n" + "=" * 70)
    print(" ALL 7 LAYERS EXECUTED SUCCESSFULLY -- EXPERIMENT FULLY COMPLETE!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run full live 7-layer VeritasBio production experiment")
    parser.add_argument("--target", type=str, default="TRPB_TM9D8", help="Target UniProt accession (e.g. TRPB_TM9D8 or P42212)")
    parser.add_argument("--budget", type=float, default=1.00, help="Total USD budget allocation (default: 1.00)")
    parser.add_argument("--rounds", type=int, default=4, help="Maximum active learning rounds (default: 4)")
    parser.add_argument("--batch-size", type=int, default=25, help="Batch size per round (default: 25)")
    parser.add_argument("--model", type=str, default="esm2_150m", help="Target scoring model (default: esm2_150m)")
    parser.add_argument("--mock-llm", action="store_true", help="Force mock LLM worker instead of Meta Model API")

    args = parser.parse_args()
    run_dollar_experiment(
        target_accession=args.target,
        budget_usd=args.budget,
        rounds=args.rounds,
        batch_size=args.batch_size,
        target_model=args.model,
        enable_live_llm=not args.mock_llm,
    )
