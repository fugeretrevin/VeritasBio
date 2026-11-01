# ADR 0001 — Estimate the biophysical prior's coefficient instead of fixing it at 1.0

**Status:** proposed, not implemented. No default changes until the measurements in §7 exist.
**Date:** 2026-11-01
**Supersedes:** nothing. **Corrects:** an earlier verbal diagnosis recorded in §2.

---

## 1. The measurement that forces this

On CreiLOV, across all 48 campaigns of `run_20261101_024544_e82e6214`, the stage-wise
Spearman rho of each ranking stage against true fitness
(`scripts/analysis/ranking_decomposition.py`, branch `worktree-agent-a83565bcde59a81f3`
@ `bb925c9`):

| stage | rho |
| --- | --- |
| GP residual alone | **+0.148** |
| composite µ(x) | **−0.032** |
| after acquisition | −0.104 |
| after Pareto | −0.121 |

and the biophysical prior on its own at **rho = −0.249, negative in 48 of 48 campaigns**.
CreiLOV is the only landscape measured so far where `|rho_prior| > rho_resid`.

The committed batch scores **0.81× its own ballot's base rate** — worse than drawing
uniformly from the candidate pool the engine itself built. Ablating the prior is worth
**+4.00 top-1% hits at 48 paired seeds** (replicating +4.25 at 12 seeds). Deleting it
outright, however, costs **12.8 top-1% hits on TrpB** (`docs/545.md` §4). So the prior
helps on one landscape and actively harms on another, and the current code cannot express
that.

---

## 2. The diagnosis that was wrong, recorded because it was acted on

It was stated in conversation on 2026-11-01 that the fix was to *"make the prior a GP mean
function the residual can learn to correct, rather than adding them on one scale."*

**That is already what the code does.** `workers/math_worker.py` trains the GP on
residuals and adds the prior back at prediction:

```python
y_train_list.append(float(score) - prior_val)        # ~line 2007, training
...
gp_means = res_means + cand_priors                   # ~line 2308, prediction
```

This is textbook universal kriging with a fixed mean function. The architecture is not the
defect. Re-deriving it as "add a mean function" would have produced a no-op change and a
claim of improvement, which is the failure mode this repository exists to document.

---

## 3. The actual defect

**The mean function's coefficient is hardcoded to exactly 1.0 and is never estimated.**

Written out, the model is `score ≈ 1.0·prior(x) + residual(x)`. The GP is handed
`score − 1.0·prior` as its target and `1.0·prior` is added back. There is no path by
which the data can say the prior is wrong.

The consequences are asymmetric and both bad:

* **A useless prior injects noise into the training targets.** Subtracting an
  uninformative signal from `y` does not cancel at predict time; it raises the variance
  the GP must explain with a fixed noise model.
* **An anti-correlated prior is actively inverted.** At rho = −0.249 the subtraction
  *adds* signal of the wrong sign to the target, and the addition at predict re-applies
  it. This is the mechanism behind composite µ(x) scoring **worse than the GP residual
  alone** (+0.148 → −0.032): the prior does not merely dilute a weakly-correct surrogate,
  it cancels and then reverses it.

Two details that matter for any implementation:

* `PRIOR_CALIBRATION_ENABLED` is **`False`** by default (`math_worker.py:1613`), so the
  tanh squashing at `:1690` is **not** in the active path. The prior is used raw, on
  whatever scale `compute_biological_prior` happens to return.
* There is a **second** composition layer downstream: when fitted,
  `means = assay_adapter.predict(X_cand, zero_shot_scores=gp_means)` (`:2310`). So
  `AssayAdaptedHead` consumes the composed mean. Any change here has two consumers, and
  one independent measurement already suggests the head is doing more than assumed —
  see §6.

---

## 4. Decision

Replace the fixed unit coefficient with an **estimated, shrunk coefficient** `β̂`, fitted
from the measured evidence each round:

```
y_train   = score − β̂ · prior
µ(x)      = β̂ · prior(x) + µ_res(x)
```

`β̂` is the shrunk least-squares coefficient of measured score on prior, over the rows the
GP is about to be fitted on.

This is the smallest change that lets the data contradict the prior. It is not a new
architecture; it is the one free parameter the existing architecture should always have
had.

### Why this has a built-in negative control

The decision is attractive mainly because its failure mode is bounded by construction:

| landscape | prior's measured rho | expected β̂ | expected effect |
| --- | --- | --- | --- |
| TrpB | positive (prior ablation costs 12.8 hits) | ≈ 1 | **behaviour approximately unchanged** |
| CreiLOV | **−0.249**, negative in 48/48 | shrunk toward 0, possibly < 0 | the damage stops |
| a landscape where the prior is noise | ≈ 0 | ≈ 0 | the prior drops out |

TrpB is therefore a self-executing control: if β̂ does not land near 1 there, the estimator
is wrong and the change must not ship, independently of any yield figure. **This is the
property that makes the change worth making rather than merely plausible.**

---

## 5. Constraints on the implementation

1. **Shrinkage is mandatory, not optional.** Round 1 fits on ~25 points and the estimate
   will be unstable. An unshrunk β̂ would thrash between rounds and the campaign would
   inherit the variance. Specify the shrinkage target (1.0, the status quo) and its
   strength *before* measuring, and record β̂ per round in the selection trace.
2. **The status quo must be exactly recoverable.** `β = 1.0` must reproduce today's
   behaviour byte-for-byte. The standard set on branch
   `worktree-agent-a14c69a80e0afa44a` is the bar: 15,840 cells across all 22 behavioural
   columns of `round_metrics.csv` agreeing exactly against the parent commit over 3 arms ×
   48 seeds × 5 rounds. Reproduce that, do not assert it.
3. **Ship as an arm first, never as a default.** Changing the composition changes every
   landscape, and `docs/HANDOFF.md` §1 already carries one withdrawn claim and two
   superseded figures caused by arms changing meaning under one name.
4. **Rule 2.12.** The manifest must record the β̂ actually used, verified against the code
   path and not against a CLI flag. Two manifests in this repository wrongly claimed `dpp`
   and `dkl` were enabled because the manifest read the flag.
5. **`AssayAdaptedHead` must be held fixed** for the first measurement, so the effect is
   attributable to the coefficient and not to the head's response to a changed input.

---

## 6. Open questions, which block nothing but must be stated

* **Is `β̂` the right functional form?** A single scalar assumes the prior is uniformly
  mis-scaled. The prior is a sum of components (steric, buried-core, charge, cohesion) and
  it is plausible that some help while others harm — in which case a per-component
  coefficient vector is the correct object and a scalar will underperform it. The scalar is
  proposed first because it is one parameter, testable at current sample sizes, and
  strictly better-specified than the fixed 1.0. **Do not extend to a vector before the
  scalar is measured.**
* **Does `AssayAdaptedHead` already absorb some of this?** An independent 48-seed result
  (branch `worktree-agent-a14c69a80e0afa44a`) suggests the head may be responsible for an
  effect previously attributed to the GP: with the prior ablated, `zero_shot_scores`
  becomes zero and the head becomes the dominant non-GP term. A measured 12-seed
  interaction of +3.67 did **not** replicate at 48 seeds (−0.19, CI [−1.04, +0.65]). If the
  head is already performing an implicit re-weighting, a second one may be redundant.
  **This should be checked before implementation, not after.**
* **What does `β̂ < 0` mean, and should it be allowed?** A negative coefficient says the
  biophysical prior is predictive with the wrong sign, which is a claim about the prior's
  construction, not a tuning outcome. Clamping at 0 is the conservative choice and
  discards information; allowing negatives is more honest and risks fitting noise at 25
  points. **Decide before measuring, and record the decision here.**

---

## 7. Pre-registration — to be committed before any arm runs

The arm is `panoptes_learned_prior_beta` (or the ablation component
`fixed_prior_coefficient`, whichever composes better with the existing vocabulary).

Primary endpoint `cum_top1_hits_excluding_initial`, 48 paired seeds, budgets matched per
arm, seed as the unit of independence, one shared starting sample per seed, 10,000-sample
whole-seed bootstrap, paired t and Wilcoxon both reported.

| landscape | prediction | falsifies the decision if |
| --- | --- | --- |
| **TrpB** | **β̂ ≈ 1**, Δ within ±3 hits of default | β̂ far from 1, **or** Δ worse than −3 — the prior demonstrably helps here, so the estimator must recover that |
| **CreiLOV** | Δ ≥ **+3** vs default, i.e. recovering a useful share of the +4.00 that full ablation gives | Δ ≤ 0 — the coefficient is not the mechanism |
| **GB1** | Δ within ±3 (prior's role unmeasured here) | — |
| **CR9114** | no prediction; the engine beats the registered reference there while violating the label scale by 14.14 SDs, and the panel shows it loses to four other baselines | — |

**Report β̂ per round per landscape regardless of the yield outcome.** If β̂ on CreiLOV is
near 1 despite rho = −0.249, the estimator is broken and no yield figure from it means
anything.

**The decision is falsified if the TrpB row fails**, whatever CreiLOV does. A change that
fixes the worst landscape by breaking the only one whose result survives a seven-arm
baseline panel is not an improvement.

---

## 8. What is explicitly out of scope

* Removing the prior. Measured: costs 12.8 top-1% hits on TrpB.
* Re-deriving the prior's components or its biophysics.
* Changing `AssayAdaptedHead`, in this change.
* Touching `controller/orchestrator.py`. Note that `test_regret_early_stopping` is
  **independently known to be flaky on unchanged code** — 6 identical pytest invocations at
  one commit gave 5 passed / 1 failed — so a failure there is not evidence about this
  change, and the test should be de-flaked before it is used to judge anything.
