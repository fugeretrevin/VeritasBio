"""
scripts/plot_benchmark.py — Publication-Quality Benchmark Visualization Engine.

Generates high-resolution dark-mode charts comparing VeritasBio active learning
discovery metrics and cloud economics against brute-force screening.

Usage:
    python scripts/plot_benchmark.py --db ledger/veritasbio_ledger.db --out scripts/benchmark_results.png
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # Non-interactive headless backend
import matplotlib.pyplot as plt
import numpy as np

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from ledger import database as ledger_db
from scripts.evaluate import evaluate_run, sequence_to_mutation_str
from workers.rag_worker import fetch_protein_metadata


def generate_benchmark_figures(
    db_path: Optional[str] = None,
    output_path: str = "scripts/benchmark_results.png",
) -> str:
    """Renders a dual-panel dark-mode figure of discovery trajectory and cloud ROI."""
    out_file = Path(output_path).resolve()
    out_file.parent.mkdir(parents=True, exist_ok=True)

    metrics = evaluate_run(db_path=db_path)
    all_evidence = ledger_db.get_evidence(db_path=db_path)
    all_jobs = ledger_db.list_jobs(db_path=db_path)

    try:
        wt_seq = fetch_protein_metadata("P42212").sequence
    except Exception:
        wt_seq = ""

    # Style configuration: VeritasBio Dark Aesthetic
    plt.style.use("dark_background")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), facecolor="#090d16")

    # Panel 1: Pareto Discovery Curve (Cumulative Spend vs Max Fitness)
    ax1.set_facecolor("#0f172a")
    ax1.grid(True, linestyle="--", alpha=0.2, color="#94a3b8")

    # Reconstruct trajectory per round
    job_costs = {j["job_id"]: float(j["estimated_cost"]) for j in all_jobs}
    cum_spend = 0.0
    spends = [0.0]
    best_scores = [0.0]
    current_max = 0.0

    # Group evidence by job
    job_ids = [j["job_id"] for j in all_jobs]
    for jid in job_ids:
        job_ev = [e for e in all_evidence if e["job_id"] == jid]
        if not job_ev:
            continue
        cum_spend += job_costs.get(jid, 0.0001)
        for e in job_ev:
            score = float(e["computed_score"])
            if score > current_max:
                current_max = score
        spends.append(cum_spend)
        best_scores.append(current_max)

    # If only 1 job or small list, ensure at least 2 points
    if len(spends) == 1:
        spends.append(metrics.total_spent_usd)
        best_scores.append(metrics.best_discovered_score)

    ax1.step(spends, best_scores, where="post", color="#00f0ff", linewidth=3, label="VeritasBio Discovery Frontier")
    ax1.scatter(spends, best_scores, color="#38bdf8", s=80, zorder=5, edgecolor="#ffffff", linewidth=1.5)

    # Highlight top discovered variant
    ax1.annotate(
        f"Top Variant: {metrics.best_discovered_variant}\nFitness: +{metrics.best_discovered_score:.4f}",
        xy=(spends[-1], best_scores[-1]),
        xytext=(spends[-1] * 0.5, best_scores[-1] * 0.85),
        arrowprops=dict(facecolor="#00f0ff", shrink=0.08, width=2, headwidth=8),
        fontsize=10,
        fontweight="bold",
        color="#ffffff",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#1e293b", edgecolor="#00f0ff", alpha=0.9),
    )

    ax1.set_title("Pareto Discovery Frontier: Fitness vs. Spend", fontsize=13, fontweight="bold", pad=15, color="#f8fafc")
    ax1.set_xlabel("Cumulative Compute Spent (USD)", fontsize=11, color="#cbd5e1")
    ax1.set_ylabel("Peak Protein Fitness Score (Normalized)", fontsize=11, color="#cbd5e1")
    ax1.legend(loc="lower right", facecolor="#1e293b", edgecolor="#334155")

    # Panel 2: Cloud Economics Comparison (VeritasBio vs Brute Force)
    ax2.set_facecolor("#0f172a")
    ax2.grid(True, linestyle="--", alpha=0.2, color="#94a3b8", axis="y")

    categories = ["Brute-Force\nCombinatorial", "VeritasBio\nDeterministic Engine"]
    costs = [metrics.estimated_bruteforce_cost_usd, metrics.total_spent_usd]
    colors = ["#ef4444", "#10b981"]

    bars = ax2.bar(categories, costs, color=colors, width=0.45, edgecolor="#ffffff", linewidth=1.2, alpha=0.9)

    for bar, cost in zip(bars, costs):
        height = bar.get_height()
        ax2.annotate(
            f"${cost:.4f}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color="#f8fafc",
        )

    # Savings Callout Badge
    ax2.text(
        0.5,
        max(costs) * 0.70,
        f"{metrics.cloud_cost_savings_pct:.1f}% Compute Savings\n({metrics.cloud_efficiency_multiple:.1f}x Cloud Efficiency)",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        color="#10b981",
        bbox=dict(boxstyle="round,pad=0.7", facecolor="#064e3b", edgecolor="#10b981", alpha=0.9),
    )

    ax2.set_title("Cloud Infrastructure Cost Comparison (Modal A10G)", fontsize=13, fontweight="bold", pad=15, color="#f8fafc")
    ax2.set_ylabel("Total Compute Cost (USD)", fontsize=11, color="#cbd5e1")

    # Header / Meta
    plt.suptitle(
        f"VeritasBio Optimization Report — avGFP (Sarkisyan et al. DMS Benchmark) | Spearman ρ = {metrics.spearman_rho:+.4f}",
        fontsize=14,
        fontweight="bold",
        color="#ffffff",
        y=0.98,
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.94])
    plt.savefig(out_file, dpi=200, facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)

    return str(out_file)


def main() -> None:
    parser = argparse.ArgumentParser(description="VeritasBio Benchmark Visualization Engine")
    parser.add_argument("--db", type=str, default=None, help="Path to SQLite ledger database")
    parser.add_argument("--out", type=str, default="scripts/benchmark_results.png", help="Output PNG path")

    args = parser.parse_args()
    saved = generate_benchmark_figures(db_path=args.db, output_path=args.out)
    print(f"Benchmark figure generated successfully: {saved}")


if __name__ == "__main__":
    main()
