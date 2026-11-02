# Pre-registration — does the ranking stack beat a coin on its own ballot?

**Written and committed before any arm of this campaign ran.** Branch
`worktree-agent-aaefbfccd52f7bce2`, implementation commit `4685572`, base `014e3cb`.
Nothing below was written with a CreiLOV, TrpB or GB1 outcome in hand; the only outcomes
in hand are the ones cited as priors, all of them already in `docs/HANDOFF.md` or on
branch `worktree-agent-a83565bcde59a81f3` (`bb925c9`).

Rule 2.11 of the handoff is the reason this file exists: `eqFP611-blue: 9.50 vs 8.79` was
published as a win from a mean difference with no paired test and was a coin flip. Every
verdict here is a paired statistic over seeds or it is not a verdict.

---

## 1. The claim under test

`scripts/analysis/ranking_decomposition.py` (on `bb925c9`) measured, over all 48 CreiLOV
campaigns, the Spearman rho of each ranking stage against true fitness:

| stage | rho |
| --- | --- |
| GP residual alone | **+0.148** |
| composite mu(x) = m(x) + mu_res(x) | −0.032 |
| after acquisition | −0.104 |
| after Pareto | **−0.121** |

and that the committed batch scores **0.81×** its own ballot's base rate — worse than
drawing uniformly from the candidate pool the engine itself built. The prior's own rho is
−0.249, **negative in 48 of 48 campaigns**.

Three competing explanations for CreiLOV's −50.29 loss against ridge are already
falsified: the spatial veto (ablating it *costs* 2.04 hits), the candidate-pool depth
ceiling (costs 1 hit of 1,655; the pre-veto ballot spans depths 1–12 and reaches 99.9% of
the hit set), and the GP's feature basis (+0.79, p=0.33, and −0.19 with the prior ablated).
What remains is the **ranking stack** and the **fabricated `y = 0.0` wild-type row**.

No existing arm can test the ranking stack as a whole. `--ablate-components` removes
stages one at a time, but something still has to order the ballot.

## 2. The arm

`panoptes_random_from_ballot` (engine `selection_rule="random_from_ballot"`). The engine's
full candidate pipeline runs unchanged — annotation-driven planning, combinatorial
generation, lead recombination, the GFlowNet sampler, the novelty filter, the spatial
viability veto — and the batch is then drawn **uniformly without replacement from the
surviving ballot**. `mu(x)`, `sigma(x)`, the prior, the acquisition value, the Pareto
multiplier and the quality score are all still computed, still written to the candidate
trace, and **never consulted**. It is the engine with its ranking stack replaced by a coin.

**What it does not isolate, registered in advance so it cannot be quietly dropped from the
write-up:** the surrogate still shapes the **ballot**, because `gfn_reward` scores GFlowNet
trajectories through it, and the GFlowNet sampler is the majority candidate source. This
arm isolates the **ranking stack**, not the surrogate.

Budget accounting and initial-sample handling are `panoptes_default`'s, unchanged:
identity-deduplicated batch fill, the same backfill contract (a longer request returns a
strict prefix-extension of the shorter draw), the same shared starting sample per seed.
Budget matching is the entire basis of the paired test, so it is **gate-verified** with
`veritas.preflight --run-dir`, not asserted.

## 3. Protocol

| item | value |
| --- | --- |
| landscapes | `creilov` (the motivating landscape), then `trpb` (the control where ranking works), then `gb1` |
| arms, one run folder each | `random`, `ridge`, `panoptes_default`, `panoptes_random_from_ballot` |
| seeds | 0–47, **48 paired**, identical starting sample per seed across arms |
| budget | 25 initial + 4 × 25 = 125 unique measurements per arm per seed |
| endpoint | `cum_top1_excluding_initial` at the final round |
| statistic | paired mean difference over seeds, paired *t*, Wilcoxon signed-rank, sign test, seeds better/worse/tied, and a 10,000-sample whole-seed bootstrap interval |
| unit of independence | the seed |
| pre-flight gate | run **with** `--run-dir`, and `--corpus-dir <real repo root>\downloads\combingym`; failures reported as failures |

Δ always means **`panoptes_random_from_ballot` − `panoptes_default`** on that endpoint.

## 4. The asymmetry, and why a split result is the informative one

These two predictions point in opposite directions, and that is the point of running both.

**CreiLOV predicts the coin WINS.** 0.81× is a direct measurement of this contrast: if the
committed batch's precision is 0.81 of its ballot's base rate, a uniform draw from that
ballot scores `1/0.81 = 1.235×` the default. At `panoptes_default` = 6.94 that is **8.57
hits, Δ = +1.6**.

**TrpB predicts the coin LOSES badly.** TrpB is the one landscape where the engine beats
the strongest of a seven-arm baseline panel (+12.88 over `gp_onehot`, t=+3.57, 31/16/1), at
51.42 hits against a 1.0% chance base rate. Whatever the ballot's base rate is, it cannot
be 51%; the ranking stack is doing most of that work. **Δ ≈ −45.**

**A result where the coin wins on CreiLOV and loses on TrpB is the informative outcome.**
It would mean the ranking stack is **mis-conditioned, not broken**: the same algebra has
real skill where its inputs are sound (TrpB, where the `y = 0.0` row is 0.35 label SDs out
and the prior is not anti-correlated) and *negative* skill where they are not (CreiLOV,
15.68 label SDs, prior rho −0.249 in 48 of 48). The remedy would then be the inputs — the
prior's sign and the WT row's scale — and not the ranking algorithm.

Committed in advance, which outcome means what:

| outcome | what it means |
| --- | --- |
| coin **wins** CreiLOV, **loses** TrpB | **the registered prediction.** The ranking stack is mis-conditioned, not broken. Fix the inputs. |
| coin wins on **all three** | the ranking stack is broken in general. This would also have to explain TrpB's +12.88, which it cannot — so this outcome would put the TrpB figure itself in doubt. |
| coin **loses on all three**, CreiLOV included | the 0.81× figure **does not reach the campaign endpoint**. The rho decomposition would then be measuring something that does not transfer — the same lesson already on record for `biophysical_features`, where rho gain and campaign gain ordered in reverse at Spearman −0.5. CreiLOV's loss would be back to unexplained by ranking. |
| coin ≈ default on CreiLOV (not separable) | the ranking stack is a **no-op** there; the loss is upstream — pool, WT row, or prior — and the 0.81× is real but endpoint-irrelevant. |

### An unresolved disagreement in my own priors, registered rather than hidden

Two figures in the handoff imply different CreiLOV answers and **I cannot reconcile them in
advance**:

- the **0.81× ratio** gives coin ≈ 8.6 hits, Δ ≈ **+1.6**;
- the ballot being "**5.2× enriched over chance**" on a 42,242-variant eligible set with
  1,655 hits (3.92% chance) gives a ballot base rate near 20% and so coin ≈ 20 hits,
  Δ ≈ **+13**. But 0.81 × 20% would make the default's precision 16%, i.e. 16 hits, and the
  default measures 6.94.

They cannot both describe the same ballot. Most likely the 5.2× is measured over the
*union* of all rounds' ballots (which is larger and more enriched than any single round's),
or pre-veto, while 0.81× is per-round and post-veto. **Registered consequence:** the
*sign* is the falsifiable claim on CreiLOV; the magnitude interval is widened to span both
derivations, and whichever derivation the result selects is itself a finding about which
of those two published figures describes the engine's real ballot.

## 5. Predictions, with intervals

Δ = `panoptes_random_from_ballot` − `panoptes_default`, top-1% hits excluding the initial
sample, 48 paired seeds. Intervals are 80% credible ranges for the paired mean.

| landscape | `panoptes_default` (known) | predicted Δ | 80% interval | predicted sign | separable? |
| --- | --- | --- | --- | --- | --- |
| **CreiLOV** | 6.94 | **+1.6** | **[−1.0, +14.0]** | **positive** | yes, t ≥ 2 |
| **TrpB** | 51.42 | **−45** | **[−50, −28]** | **negative** | yes, t ≤ −5 |
| **GB1** | 20.40 | **−16** | **[−22, −6]** | **negative** | yes, t ≤ −3 |

Reasoning per landscape:

- **CreiLOV** — §4. Widened downward because the coin's trajectory diverges: its leads are
  random ballot members rather than high-fitness ones, and `recombine_lead_variants` builds
  the next round's deep candidates *from the leads*, so a coin can degrade its own ballot
  in a way the frozen 0.81× ratio cannot see. Widened upward to cover the enrichment
  derivation.
- **TrpB** — 1,600 hits of 159,999 is 1.0% chance; the post-veto ballot passes 1,002 of
  those 1,600, and the shallow shell the engine prefers runs at 7.49% (d ≤ 2) against 0.91%
  (d ≥ 3), so a uniform ballot draw plausibly lands at 3–8 hits per 100 measurements
  against the default's 51.42.
- **GB1** — 1,494 hits of 149,360, shallow shell 3.32%. The `y = 0.0` row is 0.20 label SDs
  out here, so the surrogate is **not** poisoned and the ranking should do real work. Less
  confident than TrpB because removing k-DPP alone is worth +6.25 on GB1 (t=3.92, 36/48),
  i.e. one batch-construction stage there is already known to be harmful.

### Falsification condition

**Primary (CreiLOV).** The prediction is falsified if Δ_CreiLOV ≤ 0, **or** if Δ_CreiLOV > 0
but not separable (paired t < 2, or Wilcoxon p ≥ 0.05, or better on ≤ 24 of 48 seeds). A
positive-but-unseparable Δ is recorded as **falsified**, not as a weak win — that is rule
2.11's whole content.

**Primary (TrpB).** Falsified if Δ_TrpB ≥ 0, or if Δ_TrpB < 0 but not separable.

**Joint.** The *mis-conditioning* interpretation in §4 is falsified unless both primaries
hold in the registered directions. If CreiLOV and TrpB come out the same sign, §4's table
says what is then concluded and no third story is admissible.

**Secondary (GB1).** Sign only. A GB1 null is reportable and does not falsify the primaries.

## 6. Stage ladder — where in the stack the damage is

On CreiLOV, 48 paired seeds, cumulative ablation against `panoptes_default`. Separate run
folders, each with the identical seed block (the starting sample is
`random.Random(seed).sample(universe.variant_ids, 25)`, a pure function of the seed and the
universe, independent of the arm list — verified by comparing
`configuration.seeding.seed_block_sha256` across folders).

The rungs are chosen by reading the decomposition backwards from the endpoint, and the
third is identified from the code rather than taken from the flag list: `means = m(x) +
mu_res(x)` is where the composite loses the GP residual's +0.148, and the prior's own rho
is −0.249, negative in 48 of 48 campaigns. So the prior is the largest single step in that
table and the rung that should matter most.

| rung | ablated | rho after this stage is removed | predicted Δ vs default | interval |
| --- | --- | --- | --- | --- |
| 1 | `pareto` | −0.104 | **+0.8** | [−1, +3] |
| 2 | `pareto,acquisition` | −0.032 | **+2.5** | [0, +6] |
| 3 | `pareto,acquisition,prior` | +0.148 | **+7** | [+2, +14] |

Registered predictions:

1. **Signs all non-negative**, and the **order monotone**: Δ₁ < Δ₂ < Δ₃.
2. Rung 3 should land near the already-measured prior-only ablation, **+4.00 at 48 seeds**
   (replicating +4.25 at 12), since the prior is the dominant term.
3. **The agreement test.** If the campaign effect sizes do not match the rho decomposition
   in sign and order, then the rho decomposition is measuring something that does not reach
   the endpoint. That is registered in advance as a **reportable finding, not a failure** —
   it is the same result already on record for `biophysical_features`, where surrogate rho
   gain and campaign gain ordered in reverse at Spearman −0.5.
4. **The residual test.** If Δ(coin) > Δ₃, there is harm in the ranking stack beyond these
   three stages — the GP residual's own sign on this landscape, or `AssayAdaptedHead`, which
   sits downstream of `mu(x)` and which no rung here removes.

## 7. The confound, and the arm that separates it

**Every CreiLOV figure in this campaign is measured *through* the fabricated `y = 0.0`
wild-type row.** `select_next_micro_batch` injects one synthetic training row per round
carrying the literal label 0.0 at variance 1e-4, while a real measurement defaults to 1e-2.
CreiLOV's real labels occupy 2.793–4.196, so that row sits 11.6 units below the minimum —
8.3× the entire label range, **15.68 label standard deviations** from centre — and every
round's surrogate is fitted through it. Ablating it is worth **+29.75 hits** (12 seeds,
t=12.17, 12/0/0); a 48-seed confirmation is in flight on branch `confirm-wt-anchor-48`.

**A coin is immune to a poisoned surrogate and the default arm is not**, so part of any
random-from-ballot win on CreiLOV is the WT row rather than the ranking stack. The 2 × 2
that separates them:

|  | anchor present | `wt_anchor` ablated |
| --- | --- | --- |
| `panoptes_default` | run A | run G |
| `panoptes_random_from_ballot` | run A | run G |

Registered predictions:

1. `wt_anchor` ablation moves the **default** a lot — Δ ≈ **+30**, interval [+15, +45] at 48
   seeds, i.e. the 12-seed figure replicates.
2. `wt_anchor` ablation moves the **coin** by far less, Δ in **[−2, +4]**. It is *not*
   predicted to be exactly zero, because the surrogate shapes the ballot through
   `gfn_reward` even when it does not shape the selection; that residual movement is the
   measurement of how much the surrogate reaches the ballot.
3. Therefore, **with `wt_anchor` ablated the coin should LOSE to `panoptes_default` on
   CreiLOV.** This is the decisive separation: a coin that wins only while the surrogate is
   poisoned indicts the WT row, and a coin that wins either way indicts the ranking stack.
4. If prediction 3 fails and the coin still wins with the anchor gone, the ranking stack is
   harmful on CreiLOV **independently** of the fabricated row, and the two defects are
   additive rather than one explaining the other.

## 8. Disclosures carried in advance

- **TrpB and GB1 are development data** under handoff rule 2.4. Results on them are
  mechanism evidence, not generalization evidence. CreiLOV is the only non-development
  landscape here.
- **`threshold.<lib>` is expected to FAIL the gate** on every run folder, from the known
  reader/schema mismatch in `veritas/preflight.py::QUALIFICATION_MANIFEST_PATHS` (it
  resolves at `configuration.components.qualify_threshold`, a scalar, while the manifest
  records `qualify_quantile` and `qualify_threshold_source` as siblings). The bar it
  concerns is consumed only by the `qualify_then_optimize` / `conditional_probe` arms,
  which this campaign does not run, and the failure is symmetric across arms — so it cannot
  bias any paired contrast. It is nonetheless an un-passed gate check these figures carry
  and it is stated wherever they are.
- **`--ablate-components dpp` and `dkl` are no-ops** on today's defaults: both globals are
  False, so `panoptes_default` is already `top_n_quality`, not k-DPP. The manifest's record
  of which batch branch ran is read from `SelectedBatch.selection_rule_applied`, stamped
  inside the branch, not from the flag (rule 2.12).
- **A comparative claim must carry its comparator**, per the reporting requirement in the
  baseline-panel amendment: `ridge` is 1st of seven on CreiLOV, 2nd on TrpB and 6th of
  seven on CR9114, so "beats ridge" is not "beats the baselines". `random` and `ridge` are
  both run in every folder here for that reason.
- **No default is changed by this campaign.** `DEFAULT_METHODS`, `DPP_DIVERSITY_ENABLED`,
  `TWO_STAGE_GP_ENABLED`, `PRIOR_CALIBRATION_ENABLED` and the WT anchor are all untouched;
  `wt_anchor` is an ablatable component, not a removal.
