"""
scripts/run_side_by_side_benchmarks.py — Runs both production benchmarks side-by-side.

Executes:
  1. TrpB (TRPB_TM9D8) — 5 rounds, $1.00 budget, 3D contact map & assay adapter.
  2. avGFP (P42212) — 5 rounds, $1.00 budget, 3D contact map & assay adapter.
Then prints a rigorous, comparative side-by-side performance table.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from scripts.run_full_dollar_experiment import run_dollar_experiment
from scripts.evaluate import evaluate_run, BenchmarkMetrics


def main() -> None:
    print("\n" + "=" * 80)
    print(" VERITASBIO DUAL END-TO-END PRODUCTION BENCHMARK RUN")
    print(" Active Learning Loop with 3D Contact Coupling, GFlowNet Annealing, & Pareto Front")
    print("=" * 80 + "\n")

    # -------------------------------------------------------------
    # 1. RUN 1: TrpB Active-Site Landscape (TRPB_TM9D8)
    # -------------------------------------------------------------
    print("\n>>> [1/2] LAUNCHING BENCHMARK 1: TRPB_TM9D8 (Tryptophan Synthase Beta)...")
    t0_trpb = time.perf_counter()
    run_dollar_experiment(
        target_accession="TRPB_TM9D8",
        budget_usd=1.00,
        rounds=12,
        batch_size=20,
        target_model="esm2_150m",
        enable_live_llm=True,
    )
    t1_trpb = time.perf_counter()

    # -------------------------------------------------------------
    # 2. RUN 2: avGFP Fluorescent Protein (P42212)
    # -------------------------------------------------------------
    print("\n>>> [2/2] LAUNCHING BENCHMARK 2: P42212 (avGFP Green Fluorescent Protein)...")
    t0_gfp = time.perf_counter()
    run_dollar_experiment(
        target_accession="P42212",
        budget_usd=1.00,
        rounds=12,
        batch_size=20,
        target_model="esm2_150m",
        enable_live_llm=True,
    )
    t1_gfp = time.perf_counter()

    # -------------------------------------------------------------
    # 3. SIDE-BY-SIDE EVALUATION
    # -------------------------------------------------------------
    trpb_db = ROOT_DIR / "ledger" / "dollar_live_trpb_tm9d8.db"
    gfp_db = ROOT_DIR / "ledger" / "dollar_live_p42212.db"
    trpb_csv = str(ROOT_DIR / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv")
    gfp_csv = str(ROOT_DIR / "scripts" / "data" / "GFP_AEQVI_Sarkisyan_2016.csv")

    m_trpb = evaluate_run(
        db_path=str(trpb_db),
        benchmark_path=trpb_csv,
        target_name="TRPB_TM9D8",
    )
    m_gfp = evaluate_run(
        db_path=str(gfp_db),
        benchmark_path=gfp_csv,
        target_name="P42212",
    )

    print("\n" + "=" * 80)
    print(" VERITASBIO FULL PRODUCTION BENCHMARK — SIDE-BY-SIDE PERFORMANCE")
    print("=" * 80)
    header = f"{'Metric':<36} | {'TrpB (TRPB_TM9D8)':<20} | {'avGFP (P42212)':<20}"
    print(header)
    print("-" * len(header))
    print(f"{'Combinatorial Ground-Truth Size':<36} | {'159,999 variants':<20} | {'51,714 variants':<20}")
    print(f"{'Execution Time':<36} | {f'{t1_trpb - t0_trpb:.1f}s':<20} | {f'{t1_gfp - t0_gfp:.1f}s':<20}")
    print(f"{'Variants Sampled':<36} | {f'{m_trpb.matched_benchmark_variants} variants':<20} | {f'{m_gfp.matched_benchmark_variants} variants':<20}")
    print(f"{'Percentage of Space Screened':<36} | {f'{m_trpb.matched_benchmark_variants / 159999 * 100:.3f}%':<20} | {f'{m_gfp.matched_benchmark_variants / 51714 * 100:.3f}%':<20}")
    print(f"{'Total Dollars Allocated':<36} | {'$1.000000':<20} | {'$1.000000':<20}")
    print(f"{'Actual Dollars Spent':<36} | {f'${m_trpb.total_spent_usd:.6f}':<20} | {f'${m_gfp.total_spent_usd:.6f}':<20}")
    print(f"{'Budget Preserved (Unspent)':<36} | {f'${1.00 - m_trpb.total_spent_usd:.6f}':<20} | {f'${1.00 - m_gfp.total_spent_usd:.6f}':<20}")
    print(f"{'Top Discovered Variant':<36} | {m_trpb.best_discovered_variant[:20]:<20} | {m_gfp.best_discovered_variant[:20]:<20}")
    print(f"{'Best Model Score':<36} | {f'{m_trpb.best_discovered_score:.4f}':<20} | {f'{m_gfp.best_discovered_score:.4f}':<20}")
    print(f"{'Top Lead Ground Truth':<36} | {f'{m_trpb.best_ground_truth_score:.4f}':<20} | {f'{m_gfp.best_ground_truth_score:.4f}':<20}")
    print(f"{'True Library Global Optimum':<36} | {f'{m_trpb.global_ground_truth_optimum:.4f} ({m_trpb.global_optimum_variant[:9]})':<20} | {f'{m_gfp.global_ground_truth_optimum:.4f} ({m_gfp.global_optimum_variant[:9]})':<20}")
    print(f"{'Pct of Global Optimum Attained':<36} | {f'{m_trpb.pct_of_optimum_attained:.1f}%':<20} | {f'{m_gfp.pct_of_optimum_attained:.1f}%':<20}")
    print(f"{'Spearman Rank Correlation (rho)':<36} | {f'{m_trpb.spearman_rho:+.4f} (p={m_trpb.spearman_pvalue:.2e})':<20} | {f'{m_gfp.spearman_rho:+.4f} (p={m_gfp.spearman_pvalue:.2e})':<20}")
    print(f"{'Pearson Linear Correlation (r)':<36} | {f'{m_trpb.pearson_r:+.4f} (p={m_trpb.pearson_pvalue:.2e})':<20} | {f'{m_gfp.pearson_r:+.4f} (p={m_gfp.pearson_pvalue:.2e})':<20}")
    print(f"{'Top-Decile Hit Rate':<36} | {f'{m_trpb.top_decile_hit_rate * 100:.1f}%':<20} | {f'{m_gfp.top_decile_hit_rate * 100:.1f}%':<20}")
    print(f"{'Fold-Enrichment over Random':<36} | {f'{m_trpb.enrichment_vs_random:.2f}x':<20} | {f'{m_gfp.enrichment_vs_random:.2f}x':<20}")
    print(f"{'Cloud Capital Cost Savings':<36} | {f'{m_trpb.cloud_cost_savings_pct:.1f}%':<20} | {f'{m_gfp.cloud_cost_savings_pct:.1f}%':<20}")
    print(f"{'Cloud Efficiency Multiplier':<36} | {f'{m_trpb.cloud_efficiency_multiple:.1f}x cheaper':<20} | {f'{m_gfp.cloud_efficiency_multiple:.1f}x cheaper':<20}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
