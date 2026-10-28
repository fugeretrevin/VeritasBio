"""
scripts/run_trpb_benchmark.py — Real-World Active Learning Benchmark on the TrpB Fitness Landscape.

Purpose:
    Benchmarks VeritasBio's active-learning loop against the 159,999 empirical variant
    measurements from the combinatorially complete TrpB active-site landscape
    (Johnston, Iqbal et al. 2024).

    This is a harder and cleaner test than the avGFP assay: every one of the 20^4 combinations
    at residues 183/184/227/228 was measured, so the loop is scored against a landscape with
    no missing cells and no survivorship bias — the optimum is known to be present in the pool.

    Simulates the laboratory closed loop:
      1. RAG Worker supplies the Tm9D8* parent sequence and its active-site residues.
      2. Strategic Planner & Quant Engine propose candidate batches from the real library.
      3. The Ledger authorizes spend.
      4. Ground-truth experimental fitness values are retrieved from the empirical assay.
      5. Quant Engine fits the Gaussian Process model and selects next-round variants.

Prepare the data first:
    python scripts/download_trpb_benchmark.py
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
logger = logging.getLogger("veritasbio.trpb_benchmark")

DATA_CSV = ROOT_DIR / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv"
PARENT_ACCESSION = "TRPB_TM9D8"
PARENT_METADATA_JSON = ROOT_DIR / "workers" / "data" / f"{PARENT_ACCESSION}.json"

# The four saturated active-site residues of the published library.
ACTIVE_SITE_POSITIONS: List[int] = [183, 184, 227, 228]

LIBRARY_SIZE = 159_999


def load_fitness_library(csv_path: Path = DATA_CSV) -> Dict[str, float]:
    """Loads mutant -> measured fitness lookup table from the prepared TrpB CSV."""
    if not csv_path.exists():
        from scripts.download_trpb_benchmark import (
            build_benchmark_csv,
            download_trpb_raw,
            download_uniprot_trpb,
        )
        download_uniprot_trpb()
        download_trpb_raw()
        build_benchmark_csv()

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


def load_parent_fitness(path: Path = PARENT_METADATA_JSON) -> float:
    """
    Reads the Tm9D8* parent's measured fitness from the prepared metadata cache.
    Used only to centre the GP targets; it is never substituted for a measurement.
    """
    record = json.loads(path.read_text(encoding="utf-8"))
    return float(record["veritasbio_provenance"]["parent_measured_fitness"])


def run_trpb_active_learning(
    total_budget_usd: float = 1.00,
    rounds: int = 5,
    batch_size: int = 25,
    pool_size: Optional[int] = None,
    db_path: Optional[str] = None,
) -> Dict[str, float]:
    """
    Executes an autonomous active-learning simulation on the empirical TrpB landscape.

    `pool_size` optionally subsamples the candidate universe (seeded, reproducible) to keep
    acquisition scoring fast during a live demo. Default: the complete 159,999-variant library.
    """
    db_file = db_path or (ROOT_DIR / "ledger" / "trpb_ledger.db")
    if os.path.exists(db_file):
        try:
            os.remove(db_file)
        except OSError:
            pass

    ledger_db.init_db(db_file)
    ledger_db.create_budget(total_allocation=total_budget_usd, budget_id=1, db_path=db_file)

    logger.info("Loading empirical TrpB landscape (Johnston et al. 2024)...")
    fitness_library = load_fitness_library()
    parent_fitness = load_parent_fitness()
    all_scores = np.array(list(fitness_library.values()))
    top_1_pct_cutoff = float(np.percentile(all_scores, 99.0))
    global_best = float(np.max(all_scores))

    meta = fetch_protein_metadata(PARENT_ACCESSION)
    wt_seq = meta.sequence

    # Every library member is an active-site variant, so the pool is the whole landscape.
    rng = random.Random(42)
    items = sorted(fitness_library.items())
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

    # Featurize once; the acquisition step re-scores the same vectors every round.
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
    logger.info("STARTING TRPB ACTIVE LEARNING LOOP (%d ROUNDS)", rounds)
    logger.info("==================================================")

    gp_model = build_gp_model(random_state=42)

    consecutive_low_ei = 0
    consecutive_stagnation = 0
    prev_best_truth = -1e9

    for round_num in range(1, rounds + 1):
        rem_budget = ledger_db.get_remaining_budget(budget_id=1, db_path=db_file)
        est_cost = len(current_idx) * 0.00001 + 0.0001

        job_id = f"job-trpb-round-{round_num:03d}"
        approved = ledger_db.authorize_spend(
            job_id=job_id,
            target_model="real_landscape_trpb",
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
            # Centre targets on the Tm9D8* parent so the GP models improvement over the parent.
            train_y.append(true_score - parent_fitness)
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
        # The ledger refused the very first batch, so no model was ever fitted.
        logger.error("No variants were evaluated — budget of $%.4f is too small to run.", total_budget_usd)
        return {
            "screened": 0,
            "top_1_pct_hits": 0,
            "total_spent_usd": total_spent,
            "cost_savings_pct": 0.0,
        }

    # Evaluate Final Performance.
    #
    # Two numbers are reported, and the distinction matters when quoting them:
    #   in-sample  — the GP scored on the very points it was trained on. It measures fit,
    #                not predictive skill, and will sit near 1.0 almost regardless of the model.
    #   held-out   — the GP scored on library members it never saw. This is the honest number.
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

    # Random-screening baseline under an identical budget. Without this comparator, "13 top-1%
    # hits" is not interpretable — the loop has to beat blind sampling to be worth anything.
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

    # Cost of screening the entire library on Modal A10G, same unit rate as the avGFP benchmark.
    bruteforce_eval_cost = LIBRARY_SIZE * 0.000306 * 0.02 + 0.006
    cost_savings_pct = (1.0 - (total_spent / bruteforce_eval_cost)) * 100.0

    logger.info("==================================================")
    logger.info("TRPB LANDSCAPE BENCHMARK RESULTS (JOHNSTON 2024)")
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
    logger.info("Parent Tm9D8* Fitness:           %.4f", parent_fitness)
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
        "parent_fitness": parent_fitness,
        "holdout_n": int(holdout_idx.size),
        "spearman_rho": float(spearman_rho),
        "pearson_r": float(pearson_r),
        "insample_spearman_rho": float(insample_rho),
        "insample_pearson_r": float(insample_r),
        "total_spent_usd": total_spent,
        "cost_savings_pct": cost_savings_pct,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run active learning benchmark on the TrpB landscape")
    parser.add_argument("--rounds", type=int, default=5, help="Active learning iterations (default: 5)")
    parser.add_argument("--batch-size", type=int, default=25, help="Batch size per round (default: 25)")
    parser.add_argument("--budget", type=float, default=1.00, help="Budget in USD (default: $1.00)")
    parser.add_argument(
        "--pool-size",
        type=int,
        default=None,
        help="Subsample the candidate pool to N variants for a faster demo (default: all 159,999)",
    )
    args = parser.parse_args()

    run_trpb_active_learning(
        total_budget_usd=args.budget,
        rounds=args.rounds,
        batch_size=args.batch_size,
        pool_size=args.pool_size,
    )
