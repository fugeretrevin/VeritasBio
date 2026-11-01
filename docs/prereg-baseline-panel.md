# Pre-registration — the baseline panel

**Dated 2026-11-01.** Written and committed **before any arm in this panel was run and
before any of its outputs were read.** The commit that carries this file is the
pre-registration record; anything in it that is later contradicted is a falsified
prediction and is reported as one.

## 1. The gap this closes

Every comparative claim in this repository is against **exactly one baseline**: a one-hot
ridge regression (`RidgeSelector`, `sklearn.linear_model.Ridge(alpha=1.0)`, untuned). The
cross-protein sweep (`docs/HANDOFF.md` §1, commit `406ef58`) reports six landscapes as
`panoptes_default − ridge`. "Beats ridge" is being read as "beats the baselines", and
nobody has measured whether ridge is a strong baseline or a weak one.

Two facts make the question live rather than rhetorical:

1. On **CR9114-H1** — the one win, +3.17 at t=+6.10 — `ridge` was measured as **not
   separable from `random`**. A win over a comparator that is itself no better than
   chance is a much weaker result than it reads, and that limitation is currently stated
   nowhere.
2. On **TrpB**, `ridge` **found the global optimum** that the engine never did. So ridge
   is not uniformly weak either. It is weak somewhere and strong elsewhere, and nobody
   has characterised where.

## 2. What this study is NOT

Under `docs/HANDOFF.md` §5.4, promoting a second model to the status of **reference
baseline** amends a frozen hypothesis, which is not permitted without an explicit,
dated amendment record. Therefore:

- **`ridge` remains the registered reference baseline.** Nothing here changes a default,
  changes `DEFAULT_METHODS`, or re-points any existing comparison.
- Every model added below is an **additional arm**, measured alongside ridge.
- If the panel shows the registered reference is too weak to support the claims built on
  it, that is reported as a **finding requiring a documented amendment**, with the
  amendment text drafted but **not enacted**.

## 3. Arms

All arms select from the identical pool (`CandidateView.untested(revealed)` — the
`oracle_evaluable_pool`, applied identically to every method), receive the identical
per-seed initial sample, are budget-matched at `initial + 4 × batch`, and fall back to
`RandomSelector` in round 1 when nothing has been revealed yet — i.e. they match
`RidgeSelector`'s budget accounting and initial-sample handling exactly.

| arm | status | what it is | why a reviewer would demand it |
| --- | --- | --- | --- |
| `random` | pre-existing | uniform without replacement | **Without `random` the panel cannot calibrate anything.** It is the floor that says whether a landscape's top-1% count can discriminate methods at all. |
| `ridge` | pre-existing, **registered reference** | one-hot `Ridge(alpha=1.0)`, greedy on prediction | the incumbent; unchanged |
| `random_forest` | pre-existing, never swept | one-hot `RandomForestRegressor(n_estimators=300)` | explicitly required; the standard nonlinear comparator, and the one §7 Phase 4 already names |
| `gbm_onehot` | **new** | one-hot `GradientBoostingRegressor()`, sklearn defaults | gradient boosting is the default-choice regressor for small-n tabular problems. If the engine cannot beat an out-of-the-box GBM, "optimisation advantage" is not the right description of it. |
| `knn_hamming` | **new** | `KNeighborsRegressor(n_neighbors=5, weights="distance")` on one-hot (Euclidean on one-hot is a monotone function of Hamming distance, so the neighbour sets are Hamming neighbour sets) | asks whether anything beyond "predict a sequence's fitness from its nearest measured sequences" is needed. A memorisation baseline with no global model. |
| `gp_onehot` | **new** | `GaussianProcessRegressor`, `C()*RBF() + WhiteKernel()` with the noise term **fitted**, `normalize_y=True`, greedy on the posterior **mean** | the engine *is* a GP-based Bayesian-optimisation method. A plain GP on one-hot isolates what the engine adds *on top of a GP* from what the GP was doing anyway. Acquisition is held at greedy-on-mean deliberately: switching to UCB would vary the surrogate and the acquisition function at once. |
| `greedy_measured` | **new** | **no model at all.** Ranks untested candidates lexicographically by Hamming distance to the fitness-ordered list of measured variants: nearest to the best measured variant first, ties broken by distance to the second best, and so on, then by universe index. | the "is the model doing anything?" control. A local-exploit rule with zero learned parameters. If it matches the engine, the engine's surrogate and acquisition machinery are not what is producing the result. |
| `panoptes_default` | pre-existing | the engine, flags as of commit `7415846` (`DPP_DIVERSITY_ENABLED=False`, `TWO_STAGE_GP_ENABLED=False`) | the thing under test |

`greedy_measured` operates on full-length per-site genotype strings (`CandidateView.combos`),
which are already WT-filled by construction, and asserts equal length before differencing.
It does not `zip()` (rule 2.1).

## 4. Endpoint and test, fixed in advance

- **Primary endpoint:** `cum_top1_hits_excluding_initial` at each (method, seed)'s final
  reached round — the same endpoint as the §1 table, so the rows are comparable to it.
- **Unit of independence:** the seed. 48 paired seeds per landscape.
- **Test:** paired over seeds. Reported for every pair: mean difference, paired *t*,
  seeds better / worse / tied, Wilcoxon signed-rank *p*, sign-test *p*, and a 10,000-draw
  whole-seed bootstrap 95% CI. **A mean difference with no paired test is not a verdict**
  (rule 2.11). `|t| < 2` with roughly half the seeds better is recorded as **not
  separable**, never as a win.
- **Gate:** `veritas.preflight` is run per landscape at the budget actually about to run,
  and its verdict is reported next to that landscape's table. A failed gate means the
  landscape is not quoted (rule 2.7).

**Pre-declared selection-bias warning.** The "best baseline" row is a **maximum over
seven arms** and is therefore upward-biased by the winner's curse. It is reported *with
every individual pairwise row beside it*, and the bias is stated wherever the max is
quoted. The best baseline is chosen by its own mean endpoint — **not** by which choice
is most or least favourable to the engine. That rule is fixed here, in advance.

## 5. Predictions

Direction and rough magnitude of `arm − ridge` on the primary endpoint. Magnitudes are
deliberately coarse; the direction and the separability call are the pre-registered
content.

### 5.1 TrpB (held-out seeds 48–95, 4 saturated sites, 125 measurements)

Rationale: TrpB is a 4-site combinatorial library with strong epistasis, and ridge is an
additive model. Nonlinear surrogates should have room. Ridge nonetheless found the global
optimum across the campaign, so it is not weak here.

| arm | predicted Δ vs ridge | predicted call |
| --- | --- | --- |
| `random` | −10 to −30 | clear loss |
| `random_forest` | **+5 to +20** | **beats ridge** |
| `gbm_onehot` | **+5 to +20** | **beats ridge** |
| `gp_onehot` | 0 to +8 | ≥ ridge, possibly not separable |
| `knn_hamming` | −5 to −20 | loss |
| `greedy_measured` | **+5 to +25** | **beats ridge**; predicted the strongest trivial arm |

**Headline TrpB prediction:** the engine's reported **+7.93 over ridge does NOT survive**
against the best of {`random_forest`, `gbm_onehot`, `greedy_measured`}. The engine ties
or loses against the best baseline on TrpB.

### 5.2 CR9114-H1 (16 binary sites, development data under rule 2.4)

Rationale: a bnAb −log₁₀(K_D) assay floored at 7.0, so most of the library is censored at
the floor and the signal is a sparse, strongly epistatic step. That is the worst case for
an additive linear model, which is consistent with ridge having measured as not separable
from random here.

| arm | predicted Δ vs ridge | predicted call |
| --- | --- | --- |
| `random` | ≈ 0 | not separable (already measured) |
| `random_forest` | **+1 to +6** | **beats ridge** |
| `gbm_onehot` | **+1 to +6** | **beats ridge** |
| `gp_onehot` | 0 to +3 | ≥ ridge |
| `knn_hamming` | 0 to +4 | ≥ ridge |
| `greedy_measured` | **+1 to +6** | **beats ridge** |

**Headline CR9114 prediction:** the engine's **+3.17 over ridge does NOT survive** against
the best of the panel. The single win in this project is predicted to be a win over a
weak comparator on development data — a *mechanism* observation, not a generalization
result (rule 2.4), and one whose effect size is predicted to be an artefact of the
comparator's weakness.

### 5.3 CreiLOV (15 sites, mixed alphabet, the worst loss: −50.29 on 48/48 seeds)

Rationale: ridge is extremely strong here and the engine lost on every single seed. A
15-site DMS landscape of this kind is expected to be largely additive, which is exactly
the regime a linear model owns and a 125-point nonlinear fit does not.

| arm | predicted Δ vs ridge | predicted call |
| --- | --- | --- |
| `random` | ≪ ridge | clear loss |
| `random_forest` | −10 to +5 | **does not beat ridge** |
| `gbm_onehot` | −10 to +5 | **does not beat ridge** |
| `gp_onehot` | −5 to +8 | ≈ ridge |
| `knn_hamming` | ≪ ridge | clear loss |
| `greedy_measured` | < ridge | loss |

**Headline CreiLOV prediction:** **ridge is at or near the top of its own neighbourhood
here.** No panel arm beats it with |t| ≥ 2. The engine's −50.29 is therefore *not* an
artefact of an unfairly strong comparator; it is a real loss, and it gets worse, not
better, when the panel is widened.

### 5.4 Cross-landscape prediction

**Ridge is a fair baseline on CreiLOV and an unfair (too weak) one on CR9114-H1**, with
TrpB intermediate. The honest summary the panel is predicted to produce: of the engine's
results, **one win, and that win does not survive the panel.**

## 6. Falsification criteria, fixed in advance

"**Ridge is a fair reference baseline**" on a given landscape is **falsified** there if:

- **F1.** Any arm in this panel beats `ridge` by a margin **≥ the engine's reported
  margin over ridge on that landscape**, at |t| ≥ 2 and with ≥ 2/3 of seeds better. Then
  on that landscape "beats ridge" does not mean "beats the baselines", and every
  comparative claim scoped to it must be restated against the stronger arm.
- **F2.** `ridge` is not separable from `random` (|t| < 2) **and** some other arm beats
  ridge at |t| ≥ 2. Ridge is then demonstrably below the achievable level on that
  landscape and is a weak reference there.

"**Ridge is a fair reference baseline**" is **confirmed** on a landscape if no panel arm
beats it at |t| ≥ 2 — i.e. ridge sits at the top of its neighbourhood.

**Landscape informativeness (carried over from the HER2 lesson, §1).** If *no*
model-based arm is separable from `random` on a landscape, that landscape's top-1% count
cannot discriminate methods, and **every** verdict on it — including the engine's — is
weak regardless of sign. That is a property of the library and no budget change fixes it.

**My predictions above are falsified** if: nonlinear arms fail to beat ridge on CR9114
(§5.2); or any arm beats ridge at |t| ≥ 2 on CreiLOV (§5.3); or the engine's margin over
ridge survives intact against the whole panel on TrpB (§5.1).

## 7. Protocol

- Harness: `benchmark_panoptes_isolated.py`, invoked through
  `scripts/analysis/generality_sweep.ps1` (its `.sh` twin cannot run on this machine —
  no bash, no WSL). **No second harness is written.**
- Landscapes in priority order: **trpb** (`--base-seed 48`, held-out seeds 48–95 — the
  one landscape where nothing in the 2026-10-31 fix set changed the variant set), then
  **cr9114_h1** (`--base-seed 0`), then **creilov** (`--base-seed 0`). `--seeds 48`,
  `--rounds 4`, `--initial-size 25`, `--batch-size 25` throughout.
- Read-out: `scripts/analysis/paired_stats.py` with
  `--metric cum_top1_hits_excluding_initial --reference ridge`, plus the engine-vs-panel
  read-out with `--reference panoptes_default`.
- The tree is committed before each scored run; the harness hashes the source tree and a
  mid-run edit voids the run.
