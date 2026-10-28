"""
scripts/run_production_end_to_end.py — Full Production End-to-End Evaluation.

Executes a live end-to-end active learning campaign across all 4 benchmark datasets:
  1. avGFP (P42212)
  2. TrpB (TRPB_TM9D8)
  3. GB1 (SPG1_STRSG)
  4. Trastuzumab (HER2_TRASTUZUMAB)

Configuration:
  - Live Meta Model API (Llama 3.1 combinatorial hypothesis & residual reflection)
  - Live Modal Cloud GPU (Modal A10G ESM-2 650M tensor inference)
  - Budget: $0.50 per target
  - Batch size: 10 variants per round (3 rounds)
  - Verification: FDA 21 CFR Part 11 SHA-256 Merkle chain, Benchling LIMS, Opentrons OT-2
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

try:
    from dotenv import load_dotenv
    load_dotenv(ROOT_DIR / ".env")
except ImportError:
    pass

from controller.orchestrator import Orchestrator, OrchestratorConfig, OrchestratorEvent
from ledger import database as ledger_db
from workers.lims_worker import BenchlingClient
from workers.math_worker import analyze_biologics_developability
from workers.wetlab_worker import (
    generate_opentrons_protocol,
    generate_plate_layout,
    generate_twist_order_csv,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("veritasbio.e2e_production")

TARGETS = [
    {
        "id": "P42212",
        "name": "avGFP (Green Fluorescent Protein)",
        "pdb_id": "1EMA",
        "category": "Fluorescent Biosensor / Folding Stability",
        "benchmark_source": "Sarkisyan et al. 2016 (51,714 variants)",
    },
    {
        "id": "TRPB_TM9D8",
        "name": "TrpB Tryptophan Synthase Beta Subunit",
        "pdb_id": "8EH0",
        "category": "Industrial Biocatalysis / Catalytic Active Site",
        "benchmark_source": "Johnston et al. 2024 (159,999 variants)",
    },
    {
        "id": "SPG1_STRSG",
        "name": "Protein G (GB1 Binding Domain)",
        "pdb_id": "2QMT",
        "category": "Binding Affinity / Epistatic Fitness Peak",
        "benchmark_source": "Wu et al. 2016 (149,361 variants)",
    },
    {
        "id": "HER2_TRASTUZUMAB",
        "name": "Trastuzumab (Herceptin) CDR-H3 Fab",
        "pdb_id": "1N8Z",
        "category": "Therapeutic Biologics / Clinical Developability",
        "benchmark_source": "Mason et al. 2021 (36,391 variants)",
    },
]


def run_target_e2e(
    target_info: Dict[str, str],
    budget_usd: float = 0.50,
    rounds: int = 3,
    batch_size: int = 10,
) -> Dict[str, Any]:
    acc = target_info["id"]
    name = target_info["name"]
    pdb = target_info["pdb_id"]

    logger.info("=" * 80)
    logger.info("STARTING LIVE END-TO-END RUN FOR: %s (%s)", name, acc)
    logger.info("Target Model: Modal A10G ESM-2 650M | LLM: Live Meta Llama API | Budget: $%.2f", budget_usd)
    logger.info("=" * 80)

    db_file = ROOT_DIR / "ledger" / f"e2e_{acc.lower()}.db"
    if db_file.exists():
        try:
            os.remove(db_file)
        except OSError:
            pass

    config = OrchestratorConfig(
        target_accession=acc,
        total_budget_usd=budget_usd,
        max_iterations=rounds,
        batch_size=batch_size,
        target_model="esm2_650m",
        db_path=str(db_file),
        use_mock_llm=False,
        enable_live_llm=True,
        enable_reflection=True,
        structure_pdb_id=pdb,
        enable_regret_stopping=True,
    )

    events_log: List[Dict[str, Any]] = []

    def on_event(ev: OrchestratorEvent) -> None:
        events_log.append({
            "iteration": ev.iteration,
            "step": ev.step,
            "message": ev.message,
            "remaining_budget_usd": ev.remaining_budget_usd,
            "cumulative_spent_usd": ev.cumulative_spent_usd,
        })
        if ev.step in ("init", "strategy", "strategy_complete", "authorization", "execution", "reflection_complete", "complete"):
            logger.info("[%s] Round %02d | Left: $%.4f | %s", ev.step.upper()[:12], ev.iteration, ev.remaining_budget_usd, ev.message[:120])

    t_start = time.time()
    orchestrator = Orchestrator(config)
    summary = orchestrator.run(on_event=on_event)
    duration_s = time.time() - t_start

    # Retrieve all jobs, evidence, and audit trail blocks
    jobs = ledger_db.list_jobs(db_path=str(db_file))
    evidence = ledger_db.get_evidence(db_path=str(db_file))
    blocks = ledger_db.list_audit_blocks(db_path=str(db_file))
    is_audit_valid, block_count, audit_msg = ledger_db.verify_audit_integrity(db_path=str(db_file))

    # Wet-lab protocol generation
    top_5_evidence = sorted(evidence, key=lambda e: float(e["computed_score"]), reverse=True)[:5]
    wt_seq = orchestrator.last_metadata.sequence if orchestrator.last_metadata else ""
    
    twist_csv = ""
    ot2_script = ""
    plate_layout: Dict[str, Any] = {}
    developability_score = 1.0
    liability_summary = "N/A"

    if top_5_evidence and wt_seq:
        var_records = [
            {"variant_name": f"{acc}_lead_{i+1}", "protein_sequence": ev["variant_sequence"], "computed_score": float(ev["computed_score"])}
            for i, ev in enumerate(top_5_evidence)
        ]
        try:
            twist_csv = generate_twist_order_csv(var_records, target_organism="e_coli")
        except Exception:
            pass

        try:
            plate_layout = generate_plate_layout(var_records, plate_format="96_well", anti_edge_barrier=True)
            ot2_script = generate_opentrons_protocol([r["variant_name"] for r in var_records], target_protein=name)
        except Exception:
            pass

        if "HER2" in acc or "TRASTUZUMAB" in acc:
            best_seq = top_5_evidence[0]["variant_sequence"]
            dev_res = analyze_biologics_developability(best_seq)
            developability_score = dev_res.get("developability_score", 1.0)
            liability_summary = f"{len(dev_res.get('liabilities_flagged', []))} liabilities flagged"

    # LIMS synchronization test
    benchling_synced = False
    try:
        b_client = BenchlingClient()
        if top_5_evidence:
            manifest = b_client.export_lims_manifest(f"VeritasBio_{acc}_Plate", plate_layout if "plate_layout" in locals() else {})
            benchling_synced = bool(manifest and len(manifest) > 20)
    except Exception:
        benchling_synced = False

    return {
        "accession": acc,
        "name": name,
        "category": target_info["category"],
        "benchmark_source": target_info["benchmark_source"],
        "duration_seconds": round(duration_s, 2),
        "iterations_completed": summary.iterations_completed,
        "variants_evaluated": summary.total_variants_evaluated,
        "total_budget_usd": budget_usd,
        "total_spent_usd": summary.final_spent_usd,
        "remaining_budget_usd": summary.final_remaining_usd,
        "highest_fitness_score": summary.highest_fitness_score,
        "top_variant_sequence": summary.top_variant_sequence,
        "stop_reason": summary.stop_reason,
        "jobs_count": len(jobs),
        "evidence_count": len(evidence),
        "audit_blocks_count": len(blocks),
        "audit_valid": is_audit_valid,
        "audit_message": audit_msg,
        "twist_csv_generated": bool(twist_csv),
        "ot2_script_generated": bool(ot2_script),
        "benchling_synced": benchling_synced,
        "developability_score": developability_score,
        "events_log": events_log,
        "last_hypothesis": orchestrator.last_hypothesis.hypothesis if hasattr(orchestrator, "last_hypothesis") and orchestrator.last_hypothesis else "N/A",
    }


def main() -> None:
    print("\n" + "=" * 100)
    print("VERITASBIO FULL END-TO-END LIVE DEMO & PRODUCTION BENCHMARK")
    print("Execution Modality: Live Meta Llama API (Step 2) + Live Modal Cloud A10G GPU (Step 5)")
    print("Ledger Budget: $0.50 per dataset | Targets: 4 Distinct Empirical Biology Datasets")
    print("=" * 100 + "\n")

    results: List[Dict[str, Any]] = []

    for t_info in TARGETS:
        try:
            res = run_target_e2e(t_info, budget_usd=0.50, rounds=3, batch_size=10)
            results.append(res)
        except Exception as e:
            logger.error("Failed end-to-end run for %s: %s", t_info["name"], e, exc_info=True)
            results.append({
                "accession": t_info["id"],
                "name": t_info["name"],
                "error": str(e),
            })

    # Print summary table
    print("\n" + "=" * 115)
    print("PRODUCTION END-TO-END RESULTS SUMMARY TABLE")
    print("=" * 115)
    print(f"{'Target Protein':<28} {'Category':<22} {'Screened':<10} {'Top Fitness':<14} {'Spent ($)':<12} {'Left ($)':<10} {'Audit Chain':<12}")
    print("-" * 115)
    for r in results:
        if "error" in r:
            print(f"{r['name']:<28} ERROR: {r['error']}")
            continue
        audit_str = "VERIFIED" if r["audit_valid"] else "FAILED"
        print(
            f"{r['name'][:26]:<28} "
            f"{r['category'][:20]:<22} "
            f"{r['variants_evaluated']:<10} "
            f"{r['highest_fitness_score']:>10.4f}    "
            f"${r['total_spent_usd']:<10.6f} "
            f"${r['remaining_budget_usd']:<8.4f} "
            f"{audit_str:<12}"
        )
    print("=" * 115)

    # Save results JSON
    out_json = ROOT_DIR / "scripts" / "data" / "production_e2e_results.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed execution data saved to {out_json}")


if __name__ == "__main__":
    main()
