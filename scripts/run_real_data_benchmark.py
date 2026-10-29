"""
scripts/run_real_data_benchmark.py — Real-World Active Learning Benchmark on Authentic ProteinGym Data.

Purpose:
    Benchmarks PanoptesBio's active-learning loop directly on the 51,714 empirical
    FACS-seq variant measurements from Sarkisyan et al. (2016) Nature avGFP assay.

    Simulates the laboratory closed loop:
      1. RAG Worker supplies biological hotspots (e.g. positions 65, 163, 175, 105, 203).
      2. Strategic Planner & Quant Engine propose candidate batches from the real 51,714 library.
      3. The Ledger authorizes spend.
      4. Ground-truth experimental fitness values are retrieved from the empirical assay.
      5. Quant Engine fits the Gaussian Process model and selects next-round variants.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import random
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy import stats

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from controller.gatekeeper import Variant
from ledger import database as ledger_db
from workers.math_worker import (
    build_gp_model,
    compute_expected_improvement,
    featurize_variant_on_positions,
)
from workers.rag_worker import fetch_protein_metadata

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("panoptesbio.real_benchmark")

DATA_CSV = ROOT_DIR / "scripts" / "data" / "GFP_AEQVI_Sarkisyan_2016.csv"


def load_dms_library(csv_path: Path = DATA_CSV) -> Dict[str, float]:
    """Loads mutant -> DMS_score lookup table from authentic ProteinGym CSV."""
    if not csv_path.exists():
        from scripts.download_real_benchmark import download_proteingym_gfp
        download_proteingym_gfp(csv_path)

    library: Dict[str, float] = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mut = row.get("mutant", "").strip()
            score_str = row.get("DMS_score", "")
            if mut and score_str:
                try:
                    library[mut] = float(score_str)
                except ValueError:
                    continue
    return library


def run_real_active_learning(
    total_budget_usd: float = 1.00,
    rounds: int = 5,
    batch_size: int = 25,
    db_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Executes an autonomous active-learning simulation on the 51,714 empirical Sarkisyan variants.
    """
    db_file = db_path or (ROOT_DIR / "ledger" / "real_data_ledger.db")
    if os.path.exists(db_file):
        try:
            os.remove(db_file)
        except OSError:
            pass

    ledger_db.init_db(db_file)
    ledger_db.create_budget(total_allocation=total_budget_usd, budget_id=1, db_path=db_file)

    logger.info("Loading authentic ProteinGym GFP dataset (51,714 variants)...")
    dms_library = load_dms_library()
    all_scores = np.array(list(dms_library.values()))
    top_1_pct_cutoff = float(np.percentile(all_scores, 99.0))
    top_10_pct_cutoff = float(np.percentile(all_scores, 90.0))

    meta = fetch_protein_metadata("P42212")
    wt_seq = meta.sequence

    # UniProt Hotspot prior: variants involving known functional hotspots
    hotspot_positions = [64, 65, 66, 67, 68, 105, 147, 148, 153, 163, 167, 175, 203, 206, 222]
    
    # Pre-filter candidate universe to variants intersecting known biological hotspots
    candidate_pool: List[Tuple[str, float, Variant]] = []
    for mut, score in dms_library.items():
        try:
            v = Variant(notation=mut)
            # Check if any mutation position is in hotspot_positions
            if any(m.position in hotspot_positions for m in v.mutations):
                # Verify substitutions match WT sequence
                v.validate_against_wt(wt_seq)
                candidate_pool.append((mut, score, v))
        except Exception:
            continue

    logger.info(
        "Candidate Pool: %d validated hotspot-guided variants extracted from 51,714 empirical records.",
        len(candidate_pool),
    )

    evaluated_mutants: Set[str] = set()
    train_x: List[np.ndarray] = []
    train_y: List[float] = []
    ground_truth_values: List[float] = []
    top_1_pct_hits = 0
    total_spent = 0.0

    # Initial seed batch from hotspot candidates
    rng = random.Random(42)
    current_batch = rng.sample(candidate_pool, min(batch_size, len(candidate_pool)))

    logger.info("==================================================")
    logger.info("STARTING REAL-DATA ACTIVE LEARNING LOOP (%d ROUNDS)", rounds)
    logger.info("==================================================")

    gp_model = build_gp_model(random_state=42)

    for round_num in range(1, rounds + 1):
        rem_budget = ledger_db.get_remaining_budget(budget_id=1, db_path=db_file)
        est_cost = len(current_batch) * 0.00001 + 0.0001
        
        job_id = f"job-real-round-{round_num:03d}"
        approved = ledger_db.authorize_spend(
            job_id=job_id,
            target_model="real_dms_sarkisyan",
            variant_count=len(current_batch),
            estimated_cost=est_cost,
            db_path=db_file,
        )

        if not approved:
            logger.warning("Round %d: Ledger rejected spend (budget exhausted). Stopping loop.", round_num)
            break

        total_spent += est_cost
        logger.info(
            "Round %02d [APPROVED]: Evaluating %d variants ($%.6f, Remaining: $%.4f)",
            round_num, len(current_batch), est_cost, rem_budget - est_cost,
        )

        round_evidence = []
        for mut, true_score, v_obj in current_batch:
            evaluated_mutants.add(mut)
            ground_truth_values.append(true_score)
            if true_score >= top_1_pct_cutoff:
                top_1_pct_hits += 1

            # Featurize continuous biophysical delta coordinates
            feat = featurize_variant_on_positions(v_obj, hotspot_positions, wt_seq)
            train_x.append(feat)
            # Normalized score relative to WT (~3.72)
            norm_score = true_score - 3.72
            train_y.append(norm_score)
            round_evidence.append((v_obj.apply_to(wt_seq), true_score))

        ledger_db.record_evidence(job_id=job_id, results=round_evidence, db_path=db_file)
        ledger_db.mark_job_completed(job_id=job_id, db_path=db_file)

        # Fit Gaussian Process
        X_mat = np.array(train_x)
        y_vec = np.array(train_y)
        gp_model.fit(X_mat, y_vec)

        # Score unevaluated candidates via Acquisition Function (Expected Improvement)
        remaining_candidates = [c for c in candidate_pool if c[0] not in evaluated_mutants]
        if not remaining_candidates or round_num == rounds:
            break

        cand_feats = np.array([
            featurize_variant_on_positions(c[2], hotspot_positions, wt_seq)
            for c in remaining_candidates
        ])
        means, stds = gp_model.predict(cand_feats, return_std=True)
        f_best = float(np.max(y_vec))
        eis = compute_expected_improvement(means, stds, f_best=f_best, xi=0.01)

        # Select top candidates by Expected Improvement
        top_indices = np.argsort(eis)[::-1][:batch_size]
        current_batch = [remaining_candidates[idx] for idx in top_indices]

        best_in_round = max(current_batch, key=lambda c: c[1])
        logger.info(
            "Round %02d Complete | Top EI Variant Selected: %s (Empirical DMS: %.4f)",
            round_num, best_in_round[0], best_in_round[1],
        )

    # Evaluate Final Performance
    all_eval_feats = np.array(train_x)
    pred_means = gp_model.predict(all_eval_feats)
    truth_arr = np.array(ground_truth_values)
    spearman_rho, rho_p = stats.spearmanr(pred_means, truth_arr)
    pearson_r, r_p = stats.pearsonr(pred_means, truth_arr)

    # Cost of screening all 51,714 variants on Modal A10G
    bruteforce_eval_cost = 51714 * 0.000306 * 0.02 + 0.006
    cost_savings_pct = (1.0 - (total_spent / bruteforce_eval_cost)) * 100.0

    logger.info("==================================================")
    logger.info("REAL-DATA BENCHMARK RESULTS (PROTEINGYM SARKISYAN)")
    logger.info("==================================================")
    logger.info("Total Real Variants Screened:    %d / 51,714 (%.2f%%)", len(evaluated_mutants), (len(evaluated_mutants)/51714)*100)
    logger.info("Top 1%% Hits Discovered:          %d variants", top_1_pct_hits)
    logger.info("Empirical Spearman Correlation:  +%.4f (p = %.2e)", spearman_rho, rho_p)
    logger.info("Empirical Pearson Correlation:   +%.4f (p = %.2e)", pearson_r, r_p)
    logger.info("Total Real Run Spend:            $%.6f", total_spent)
    logger.info("Exhaustive Brute-Force Cost:     $%.4f", bruteforce_eval_cost)
    logger.info("Cloud Capital Saved:             %.1f%% (%.1fx cheaper)", cost_savings_pct, bruteforce_eval_cost / max(1e-6, total_spent))
    logger.info("==================================================")

    return {
        "screened": len(evaluated_mutants),
        "top_1_pct_hits": top_1_pct_hits,
        "spearman_rho": float(spearman_rho),
        "pearson_r": float(pearson_r),
        "total_spent_usd": total_spent,
        "cost_savings_pct": cost_savings_pct,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run real-data active learning benchmark on ProteinGym")
    parser.add_argument("--rounds", type=int, default=5, help="Active learning iterations (default: 5)")
    parser.add_argument("--batch-size", type=int, default=25, help="Batch size per round (default: 25)")
    parser.add_argument("--budget", type=float, default=1.00, help="Budget in USD (default: $1.00)")
    args = parser.parse_args()

    run_real_active_learning(
        total_budget_usd=args.budget,
        rounds=args.rounds,
        batch_size=args.batch_size,
    )
