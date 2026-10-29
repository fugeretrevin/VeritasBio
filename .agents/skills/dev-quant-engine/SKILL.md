---
name: dev-quant-engine
description: >-
  Use this skill when modifying active learning algorithms, Gaussian Process regression,
  GFlowNet generative candidate sampling, Expected Improvement acquisition functions,
  or multi-objective Pareto optimization in workers/math_worker.py.
---

# Quant Engine Development Skill

This skill governs the mathematical and statistical engine of PanoptesBio (`workers/math_worker.py`).

## Core Invariants
1. **Separation of Concerns:**
   - The LLM only proposes search boundaries, biophysical hypotheses, or pocket seed coordinates.
   - The math engine (`workers/math_worker.py`) performs ALL numerical selections, acquisition function scoring, uncertainty estimation, and ranking.
2. **Held-Out Generalization over In-Sample Fit:**
   - Never report in-sample $R^2$ or $\rho$ as the primary metric. Always validate on unseen holdout test sets.
   - Compare enrichment factors against a random screening baseline (e.g. 200 random draws).
3. **Replicate Variance & Heteroscedastic Noise:**
   - Support heteroscedastic observation noise when wet-lab replicates are available.
4. **Multi-Objective Pareto Tradeoffs:**
   - Use Kung / Deb fast non-dominated sorting (`compute_non_dominated_pareto_ranks`) to balance fitness, synthesizability, and developability liabilities.

## Key APIs & Functions
- `build_gp_model(random_state=42)`: Gaussian Process Regressor with Matern 5/2 kernel.
- `featurize_variants(variants, wt_sequence)`: Physicochemical and 3D contact embedding.
- `compute_expected_improvement(mean, std, f_best, xi=0.01)`: Closed-form analytical EI.
- `sample_gflownet_candidates(...)`: Energy-guided Markov sampling with temperature annealing.
- `compute_non_dominated_pareto_ranks(objectives)`: Multi-objective NSGA-II frontier.
- `analyze_biologics_developability(sequence)`: Detection of deamidation, isomerization, and odd cysteines.

## Verification Checklist
```powershell
python -m unittest tests/test_math_worker.py
python -m unittest tests/test_enterprise_hardening.py
```
Ensure:
- Gaussian Process fit executes deterministically.
- GFlowNet annealing temperature lowers over iterations.
- Pareto ranking correctly isolates non-dominated Rank 0 leads.
