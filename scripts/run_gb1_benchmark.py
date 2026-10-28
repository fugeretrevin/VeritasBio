"""
scripts/run_gb1_benchmark.py — Real-World Active Learning Benchmark on the Protein G (GB1) Landscape.

Purpose:
    Benchmarks VeritasBio's deterministic cloud-economics active-learning loop
    against the 149,361 empirical variant measurements from the four-site combinatorial
    landscape of Streptococcal Protein G B1 domain (Wu, Dai, Olson et al. 2016, eLife).

    Every one of the 20^4 = 160,000 combinations at residues 39, 40, 41, 54 was assayed for
    IgG Fc binding affinity. 149,361 variants yielded high-confidence sequencing counts.
    This provides an authentic, epistatically dense fitness landscape with zero survivorship bias.

Simulates the laboratory closed loop:
    1. RAG Worker supplies the wildtype GB1 sequence and active IgG Fc interface residues.
    2. Quant Engine featurizes candidates and fits Gaussian Process models.
    3. Ledger authorizes compute spend under fail-closed economics.
    4. Ground-truth experimental binding fitness values are retrieved.
    5. Expected Improvement (EI) selects the optimal next micro-batch.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import random
import sys
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
logger = logging.getLogger("veritasbio.gb1_benchmark")

DATA_CSV = ROOT_DIR / "scripts" / "data" / "GB1_STRSG_Wu_2016.csv"
PARENT_ACCESSION = "SPG1_STRSG"
PARENT_METADATA_JSON = ROOT_DIR / "workers" / "data" / f"{PARENT_ACCESSION}.json"

ACTIVE_SITE_POSITIONS: List[int] = [39, 40, 41, 54]
LIBRARY_SIZE = 149_361


def load_fitness_library(csv_path: Path = DATA_CSV) -> Dict[str, float]:
    """Loads mutant -> measured fitness lookup table from the prepared GB1 CSV."""
    if not csv_path.exists():
        from scripts.download_gb1_benchmark import download_gb1_dataset
        download_gb1_dataset()

    library: Dict[str, float] = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            mut = row.get("mutant", "").strip()
            score_str = row.get("DMS_score", "")
            if mut and score_str:
                try:
                    library[mut] = float(score_str)
                except ValueError:
                    continue
    return library


def run_gb1_active_learning(
    total_budget_usd: float = 1.00,
    rounds: int = 5,
    batch_size: int = 25,
    pool_size: Optional[int] = None,
    db_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Executes an autonomous active-learning simulation on the empirical GB1 landscape.

    Args:
        total_budget_usd: Compute budget allocation.
        rounds: Number of active learning iterations.
        batch_size: Variants evaluated per micro-batch.
        pool_size: Optional candidate subsample (seeded) for fast evaluation.
        db_path: Optional custom path for the SQLite ledger.
    """
    db_file = db_path or (ROOT_DIR / "ledger" / "gb1_ledger.db")
    if os.path.exists(db_file):
        try:
            os.remove(db_file)
        except OSError:
            pass

    ledger_db.init_db(db_file)
    ledger_db.create_budget(total_allocation=total_budget_usd, budget_id=1, db_path=db_file)

    logger.info("Loading empirical GB1 landscape (Wu et al. 2016, eLife)...")
    fitness_library = load_fitness_library()
    all_scores = np.array(list(fitness_library.values()))
    top_1_pct_cutoff = float(np.percentile(all_scores, 99.0))
    global_best = float(np.max(all_scores))
    wt_fitness = fitness_library.get("WT", 1.0)

    meta = fetch_protein_metadata(PARENT_ACCESSION)
    wt_seq = meta.sequence

    rng = random.Random(42)
    items = sorted(fitness_library.items())
    # Exclude WT from candidate mutation sampling if present
    items = [it for it in items if it[0] != "WT"]

    if pool_size is not None and pool_size < len(items):
        items = rng.sample(items, pool_size)
        logger.info("Subsampled candidate pool to %d variants (seeded).", pool_size)

    candidate_pool: List[Tuple[str, float, Variant]] = []
    for mut, score in items:
        try:
            v = Variant(notation=mut)
            v.validate_against_wt(wt_seq)
            candidate_pool.append((mut, score, v))
        except Exception:
            continue

    logger.info(
        "Candidate Pool: %d validated active-site variants from %d empirical records.",
        len(candidate_pool), len(fitness_library),
    )

    logger.info("Featurizing candidate pool on positions %s...", ACTIVE_SITE_POSITIONS)
    pool_feats = np.array([
        featurize_variant_on_positions(c[2], ACTIVE_SITE_POSITIONS, wt_seq)
        for c in candidate_pool
    ])

    evaluated_mutants: Set[str] = set()
    train_idx: List[int] = []
    train_y: List[float] = []
    ground_truth_values: List[float] = []
    top_1_pct_hits = 0
    total_spent = 0.0

    current_idx = rng.sample(range(len(candidate_pool)), min(batch_size, len(candidate_pool)))

    logger.info("==================================================")
    logger.info("STARTING GB1 ACTIVE LEARNING LOOP (%d ROUNDS)", rounds)
    logger.info("==================================================")

    gp_model = build_gp_model(random_state=42)

    consecutive_low_ei = 0
    consecutive_stagnation = 0
    prev_best_truth = -1e9

    for round_num in range(1, rounds + 1):
        rem_budget = ledger_db.get_remaining_budget(budget_id=1, db_path=db_file)
        est_cost = len(current_idx) * 0.00001 + 0.0001

        job_id = f"job-gb1-round-{round_num:03d}"
        approved = ledger_db.authorize_spend(
            job_id=job_id,
            target_model="real_landscape_gb1",
            variant_count=len(current_idx),
            estimated_cost=est_cost,
            db_path=db_file,
        )

        if not approved:
            logger.warning("Round %d: Ledger rejected spend (budget exhausted). Stopping loop.", round_num)
            break

        total_spent += est_cost
        logger.info(
            "Round %02d [APPROVED]: Evaluating %d variants ($%.6f, Remaining: $%.4f)",
            round_num, len(current_idx), est_cost, rem_budget - est_cost,
        )

        round_evidence = []
        for idx in current_idx:
            mut, true_score, v_obj = candidate_pool[idx]
            evaluated_mutants.add(mut)
            ground_truth_values.append(true_score)
            if true_score >= top_1_pct_cutoff:
                top_1_pct_hits += 1

            train_idx.append(idx)
            # Centre targets on WT fitness (1.0)
            train_y.append(true_score - wt_fitness)
            round_evidence.append((v_obj.apply_to(wt_seq), true_score))

        ledger_db.record_evidence(job_id=job_id, results=round_evidence, db_path=db_file)
        ledger_db.mark_job_completed(job_id=job_id, db_path=db_file)

        # Fit Gaussian Process
        X_mat = pool_feats[train_idx]
        y_vec = np.array(train_y)
        gp_model.fit(X_mat, y_vec)

        remaining_idx = np.array(
            [i for i in range(len(candidate_pool)) if candidate_pool[i][0] not in evaluated_mutants]
        )
        if remaining_idx.size == 0 or round_num == rounds:
            break

        # Score unevaluated candidates via Acquisition Function (Expected Improvement)
        means, stds = gp_model.predict(pool_feats[remaining_idx], return_std=True)
        f_best = float(np.max(y_vec))
        eis = compute_expected_improvement(means, stds, f_best=f_best, xi=0.01)

        # 1. Regret stopping check: If Expected Improvement drops to near zero
        max_ei = float(np.max(eis)) if len(eis) > 0 else 0.0
        if round_num >= 3 and max_ei < 0.005:
            consecutive_low_ei += 1
            if consecutive_low_ei >= 2:
                logger.info(
                    "Bayesian regret stopping triggered at round %d: max EI %.5f < 0.005 (patience: %d/2). Preserving unspent budget.",
                    round_num, max_ei, consecutive_low_ei,
                )
                break
        else:
            consecutive_low_ei = 0

        # 2. Fitness plateau check: If empirical top fitness stalls over consecutive rounds
        current_best_truth = float(np.max(ground_truth_values))
        if round_num >= 4 and (current_best_truth - prev_best_truth) < 0.001:
            consecutive_stagnation += 1
            if consecutive_stagnation >= 3:
                logger.info(
                    "Fitness plateau stopping triggered at round %d: no improvement > 0.001 over 3 consecutive rounds. Preserving unspent budget.",
                    round_num,
                )
                break
        else:
            consecutive_stagnation = 0
            prev_best_truth = current_best_truth

        top_local = np.argsort(eis)[::-1][:batch_size]
        current_idx = [int(remaining_idx[i]) for i in top_local]

        best_in_round = max(current_idx, key=lambda i: candidate_pool[i][1])
        logger.info(
            "Round %02d Complete | Top EI Variant Selected: %s (Empirical Fitness: %.4f)",
            round_num, candidate_pool[best_in_round][0], candidate_pool[best_in_round][1],
        )

    if not train_idx:
        logger.error("No variants were evaluated — budget of $%.4f is too small to run.", total_budget_usd)
        return {
            "screened": 0,
            "top_1_pct_hits": 0,
            "total_spent_usd": total_spent,
            "cost_savings_pct": 0.0,
        }

    # Evaluate Final Performance (held-out vs in-sample)
    truth_arr = np.array(ground_truth_values)
    insample_means = gp_model.predict(pool_feats[train_idx])
    insample_rho, insample_p = stats.spearmanr(insample_means, truth_arr)
    insample_r, _ = stats.pearsonr(insample_means, truth_arr)

    holdout_idx = np.array(
        [i for i in range(len(candidate_pool)) if candidate_pool[i][0] not in evaluated_mutants]
    )
    if holdout_idx.size >= 2:
        holdout_means = gp_model.predict(pool_feats[holdout_idx])
        holdout_truth = np.array([candidate_pool[int(i)][1] for i in holdout_idx])
        spearman_rho, rho_p = stats.spearmanr(holdout_means, holdout_truth)
        pearson_r, r_p = stats.pearsonr(holdout_means, holdout_truth)
    else:
        spearman_rho = rho_p = pearson_r = r_p = float("nan")

    best_found = float(np.max(truth_arr))

    # Random-screening baseline
    n_screened = len(evaluated_mutants)
    baseline_rng = random.Random(1234)
    baseline_hits: List[int] = []
    baseline_best: List[float] = []
    for _ in range(200):
        draw = baseline_rng.sample(range(len(candidate_pool)), n_screened)
        vals = [candidate_pool[i][1] for i in draw]
        baseline_hits.append(sum(1 for v in vals if v >= top_1_pct_cutoff))
        baseline_best.append(max(vals))
    rand_hits = float(np.mean(baseline_hits))
    rand_best = float(np.mean(baseline_best))
    enrichment = top_1_pct_hits / rand_hits if rand_hits > 0 else float("inf")

    bruteforce_eval_cost = LIBRARY_SIZE * 0.000306 * 0.02 + 0.006
    cost_savings_pct = (1.0 - (total_spent / bruteforce_eval_cost)) * 100.0

    logger.info("==================================================")
    logger.info("PROTEIN G (GB1) BENCHMARK RESULTS (WU ET AL. 2016)")
    logger.info("==================================================")
    logger.info(
        "Variants Screened:               %d / %d (%.3f%%)",
        len(evaluated_mutants), len(candidate_pool),
        (len(evaluated_mutants) / max(1, len(candidate_pool))) * 100,
    )
    logger.info("Top 1%% Hits Discovered:          %d variants", top_1_pct_hits)
    logger.info("  vs random screening (n=200):   %.2f hits  -> %.2fx enrichment", rand_hits, enrichment)
    logger.info("Best Fitness Found:              %.4f", best_found)
    logger.info("  vs random screening (n=200):   %.4f", rand_best)
    logger.info("Global Optimum in Library:       %.4f", global_best)
    logger.info("Wildtype GB1 Fitness:            %.4f", wt_fitness)
    logger.info("--- Predictive skill on %d UNSEEN library members ---", holdout_idx.size)
    logger.info("Held-Out Spearman Correlation:   %+.4f (p = %.2e)", spearman_rho, rho_p)
    logger.info("Held-Out Pearson Correlation:    %+.4f (p = %.2e)", pearson_r, r_p)
    logger.info("--- Fit on the %d training points (NOT predictive skill) ---", len(train_idx))
    logger.info("In-Sample Spearman Correlation:  %+.4f (p = %.2e)", insample_rho, insample_p)
    logger.info("In-Sample Pearson Correlation:   %+.4f", insample_r)
    logger.info("Total Real Run Spend:            $%.6f", total_spent)
    logger.info("Exhaustive Brute-Force Cost:     $%.4f", bruteforce_eval_cost)
    logger.info(
        "Cloud Capital Saved:             %.1f%% (%.1fx cheaper)",
        cost_savings_pct, bruteforce_eval_cost / max(1e-6, total_spent),
    )
    logger.info("==================================================")

    return {
        "screened": len(evaluated_mutants),
        "top_1_pct_hits": top_1_pct_hits,
        "best_fitness_found": best_found,
        "random_baseline_hits": rand_hits,
        "random_baseline_best": rand_best,
        "enrichment_vs_random": enrichment,
        "global_optimum": global_best,
        "parent_fitness": wt_fitness,
        "holdout_n": int(holdout_idx.size),
        "spearman_rho": float(spearman_rho),
        "pearson_r": float(pearson_r),
        "insample_spearman_rho": float(insample_rho),
        "insample_pearson_r": float(insample_r),
        "total_spent_usd": total_spent,
        "cost_savings_pct": cost_savings_pct,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run active learning benchmark on the GB1 landscape")
    parser.add_argument("--rounds", type=int, default=5, help="Active learning iterations (default: 5)")
    parser.add_argument("--batch-size", type=int, default=25, help="Batch size per round (default: 25)")
    parser.add_argument("--budget", type=float, default=1.00, help="Budget in USD (default: $1.00)")
    parser.add_argument(
        "--pool-size",
        type=int,
        default=None,
        help="Subsample the candidate pool to N variants for a faster demo (default: all 149,361)",
    )
    args = parser.parse_args()

    run_gb1_active_learning(
        total_budget_usd=args.budget,
        rounds=args.rounds,
        batch_size=args.batch_size,
        pool_size=args.pool_size,
    )
