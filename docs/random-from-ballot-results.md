# Results — does the engine's ranking stack beat a coin on its own ballot?

Pre-registered in [`prereg-random-from-ballot.md`](prereg-random-from-ballot.md), committed
at `3584b7b` **before any arm of this campaign ran**. Implementation `4685572`, base
`014e3cb`. Every run below carries `commit = 3584b7b` and
`repository_unchanged = true`.

Seven run folders, 48 paired seeds each, 25 initial + 4 × 25 = 125 measurements per arm per
seed. Budgets gate-verified: **every arm spent all 125 distinct measurements on all 48
seeds, in all seven folders**.

**Every contrast is reported on both endpoints**: raw top-1% hits excluding the initial
sample, and the diversity-adjusted count at d ≥ 2 (fitness-descending greedy, keep a hit
only if ≥ 2 substitutions from everything kept — `veritas.calibration.endpoint.greedy_diverse`,
not a sixth implementation). They disagree in three places and each disagreement is called
out. The reason this is not optional: the combined head-to-head campaign measured a one-hot
surrogate gain of **+5.00 raw hits on TrpB at t=+5.33 that becomes +0.10 at d ≥ 2 at
t=+0.29, not separable**, with redundancy rising 3.63× → 3.96×. A raw count cannot tell
"found more good variants" from "found the same basin again", and a coin drawing from a
dense ballot is exactly the intervention where that distinction can invert a conclusion.

---

## 1. The answer

**The ranking stack beats a coin on two landscapes of three, and loses to one — and where
it loses, the cause is not the ranking stack.**

| landscape | raw Δ (coin − default) | d ≥ 2 Δ | verdict |
| --- | --- | --- | --- |
| **CreiLOV** | **+5.438** | **+4.604** | **the coin WINS, on both endpoints** |
| **TrpB** | **−37.250** | **−4.938** | the coin loses; ranking does real work |
| **GB1** | **−10.979** | **−0.604, NOT separable** | raw loss only; **nothing distinct** is gained |

And the follow-up arm settles what the CreiLOV win means: with the fabricated `y = 0.0`
wild-type training row ablated, the **same coin loses to the same default by −25.583 raw /
−10.312 at d ≥ 2, on 48 of 48 seeds**. The ranking stack is not broken. On CreiLOV it is
**disabled by its own training data**, and a coin is immune to that where the stack is not.

---

## 2. Pre-registration outcome

| registered prediction | predicted | measured | outcome |
| --- | --- | --- | --- |
| CreiLOV Δ positive and separable; point +1.6, 80% [−1.0, +14.0] | positive | **+5.438**, t=+6.11, 35/10/3 | **CONFIRMED**, inside the interval |
| TrpB Δ negative and separable; point −45, 80% [−50, −28] | negative | **−37.250**, t=−38.38, 0/48/0 | **CONFIRMED**, inside the interval |
| GB1 Δ negative (sign only, secondary) | negative | **−10.979** raw, t=−7.55 | **CONFIRMED on raw**; see §4, the d ≥ 2 contrast is null |
| the §4 *joint* prediction: coin wins CreiLOV **and** loses TrpB → the stack is mis-conditioned, not broken | split result | split result | **CONFIRMED** — this was registered as the informative outcome |
| ladder: all rungs non-negative, order monotone Δ₁ < Δ₂ < Δ₃ | monotone | monotone on raw, **NOT on d ≥ 2** | **PARTLY FALSIFIED**, §5 |
| rung 3 lands near the known prior-only +4.00 | ≈ +4 to +7 | **+5.083** | CONFIRMED |
| `wt_anchor` moves the default a lot: +30, [+15, +45] | +30 | **+31.771** | **CONFIRMED**; replicates the 12-seed +29.75 at 48 seeds |
| `wt_anchor` moves the coin little: [−2, +4] | ≈ 0 | **+0.750**, t=+0.93, p=0.356, 24/22/2 | **CONFIRMED**, not separable |
| with `wt_anchor` ablated the coin LOSES to the default | loses | **−25.583**, 0/48/0 | **CONFIRMED** |
| the registered unresolved disagreement: does the ballot's base rate match 0.81× or "5.2× enriched"? | unresolved in advance | **neither; it is 12.4× chance**, §7 | RESOLVED empirically |

**Nothing in the primary registration was falsified.** The one partial falsification is the
ladder's monotonicity on the diversity-adjusted endpoint, and it is informative rather than
inconvenient — see §5.

---

## 3. Primary paired tables

Δ = `panoptes_random_from_ballot` − `panoptes_default`, paired on seed, n = 48.
`b/w/t` = seeds better / worse / tied. CI is a 10,000-sample whole-seed bootstrap on the
paired mean. Three tests are quoted because the project has twice called a coin flip a win
from one of them.

### Arm means

| landscape | arm | raw | d ≥ 2 | redundancy (raw / d≥2) |
| --- | --- | --- | --- | --- |
| CreiLOV | ridge | 57.229 ± 1.642 | 23.208 ± 0.636 | 2.47 |
| | **coin** | **12.375 ± 0.644** | **9.667 ± 0.472** | **1.28** |
| | panoptes_default | 6.938 ± 0.773 | 5.062 ± 0.540 | 1.37 |
| | random | 1.167 ± 0.138 | 1.167 ± 0.138 | 1.00 |
| TrpB | panoptes_default | 49.646 ± 0.819 | 13.688 ± 0.256 | **3.63** |
| | ridge | 31.104 ± 2.248 | 7.229 ± 0.529 | **4.30** |
| | **coin** | **12.396 ± 0.561** | **8.750 ± 0.383** | **1.42** |
| | random | 0.938 ± 0.117 | 0.917 ± 0.115 | 1.02 |
| GB1 | ridge | 33.292 ± 1.998 | 6.896 ± 0.372 | **4.83** |
| | panoptes_default | 20.396 ± 1.302 | 6.938 ± 0.304 | 2.94 |
| | **coin** | **9.417 ± 0.613** | **6.333 ± 0.335** | **1.49** |
| | random | 1.083 ± 0.136 | 1.083 ± 0.136 | 1.00 |

### The primary contrast

| landscape | endpoint | Δ | t | p(t) | Wilcoxon | sign | CI95 | b/w/t | separable |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CreiLOV | raw | **+5.438** | +6.11 | 1.9e−07 | 3.0e−06 | 2.5e−04 | [+3.71, +7.15] | 35/10/3 | **yes** |
| CreiLOV | d ≥ 2 | **+4.604** | +6.57 | 3.7e−08 | 9.9e−07 | 4.2e−05 | [+3.23, +5.96] | 35/8/5 | **yes** |
| TrpB | raw | **−37.250** | −38.38 | 3.9e−37 | 1.6e−09 | 7.1e−15 | [−39.06, −35.31] | 0/48/0 | **yes** |
| TrpB | d ≥ 2 | **−4.938** | −11.25 | 6.3e−15 | 6.9e−09 | 5.9e−11 | [−5.79, −4.08] | 2/43/3 | **yes** |
| GB1 | raw | **−10.979** | −7.55 | 1.2e−09 | 5.3e−08 | 2.5e−08 | [−13.79, −8.25] | 5/42/1 | **yes** |
| GB1 | d ≥ 2 | **−0.604** | −1.34 | 0.185 | 0.234 | 1.000 | [−1.48, +0.25] | 24/23/1 | **NO** |

### Where the two endpoints disagree

1. **GB1: the ranking stack's entire advantage over a coin is redundant hits.** −10.979
   raw collapses to −0.604 at d ≥ 2, t=−1.34, better on 24 of 48 seeds — a coin flip. The
   stack finds 2.2× as many hits as the coin and **not one additional distinct one**.
2. **TrpB: the advantage is real but 7.5× smaller than raw claims.** −37.250 → −4.938.
   `panoptes_default`'s TrpB redundancy is 3.63 against the coin's 1.42, so most of its
   raw margin is re-measurement inside one basin.
3. **TrpB, `coin − ridge`: the sign INVERTS.** −18.708 raw (t=−7.92, 7/41/0, separable)
   becomes **+1.521 at d ≥ 2** (t=+2.26, Wilcoxon 0.031, **sign test 0.174**, 27/17/4).
   The t and Wilcoxon tests cross 0.05 and the **sign test does not**: 27 of 44 non-tied
   seeds is a coin flip on a seed-majority basis. **Reported as NOT separable.** This is
   exactly the shape of result rule 2.11 exists for, and the baseline-panel reporting
   requirement is why both tests are quoted. The defensible statement is: at d ≥ 2 on
   TrpB, a coin drawing from the engine's ballot is **indistinguishable from one-hot
   ridge**, while `ridge` holds a 4.30 redundancy against the coin's 1.42.
4. **CreiLOV: the endpoints AGREE.** +5.438 → +4.604, and the d ≥ 2 contrast is *more*
   significant. The coin's CreiLOV win is not a pile of near-duplicates; it is distinct
   variants. This is the contrast the coordinator's warning was aimed at, and it survives.

**Comparators, named.** On CreiLOV `ridge` is the strongest of the seven-arm panel (57.23),
so "the coin still loses to ridge by −44.854 raw / −13.542 at d ≥ 2, 0/48/0" is a loss
against the strongest available baseline. On **GB1 `ridge` is NOT the strongest raw
baseline — `gp_onehot` is, at 46.02 against ridge's 33.29** (carried from today's
head-to-head, not re-derived here), so every GB1 contrast below is *against ridge*, with
comparator strength known to be understated. On TrpB `ridge` is 2nd of seven.

---

## 4. The direct answer

**Does the engine's ranking stack beat a coin on its own ballot, and on which landscapes?**

- **TrpB: yes, decisively, on both endpoints.** +37.250 raw and +4.938 at d ≥ 2, 48 of 48
  seeds raw. The stack is doing real work where its inputs are sound (the `y = 0.0` row is
  0.35 label SDs out, and `wt_anchor` ablation there is +0.542, t=+0.41 — carried).
- **GB1: on raw yes, on distinct hits NO.** +10.979 raw, but **−0.604 at d ≥ 2 is not
  separable from zero**. Measured at the endpoint a lab can use, the stack buys nothing a
  coin would not.
- **CreiLOV: no — the coin beats it, on both endpoints**, by +5.438 raw and +4.604 at d ≥ 2.
  §6 shows this is the fabricated wild-type row, not the ranking algebra: remove that row
  and the stack beats the coin by 25.583 raw hits on 48 of 48 seeds.

So the 0.81× figure that motivated this campaign **reaches the endpoint on CreiLOV and
reverses everywhere else**, and the conclusion the campaign supports is the one the
pre-registration named as informative: the ranking stack is **mis-conditioned, not broken**.

One number worth stating plainly. On CreiLOV the gap to ridge is −50.292 raw. The coin
closes **5.438 of it, 10.8%**; ablating the wild-type row closes **31.771, 63.2%**. The
ranking stack is the smaller of the two defects by a factor of six.

---

## 5. The stage ladder, and whether it agrees with the rho decomposition

CreiLOV, 48 paired seeds, cumulative ablation. Each rung is its own run folder;
**all five CreiLOV folders carry the identical seed block `29639d2992de02b2…`, and `ridge`
is bit-identical per seed across all of them** — which is the direct test that the starting
samples and the evaluator really are the same, not an assumption.

| rung | ablated | Δ raw vs default | Δ d ≥ 2 | separable |
| --- | --- | --- | --- | --- |
| 1 | `pareto` | +0.417 (t=+0.96, 23/19/6) | +0.625 (t=+1.85, 24/17/7) | **no, either endpoint** |
| 2 | `pareto,acquisition` | +0.521 (t=+1.15, 22/18/8) | −0.292 (t=−1.00, 17/25/6) | **no, either endpoint** |
| 3 | `pareto,acquisition,prior` | **+5.083** (t=+4.57, 34/11/3) | **+2.688** (t=+4.13, 32/13/3) | **yes, both** |

The third rung is the **prior**, identified from the code and not from the flag list:
`means = m(x) + mu_res(x)` is where the composite loses the GP residual's +0.148, and the
prior's own rho is −0.249, negative in 48 of 48 campaigns. It is the only rung that reaches
the endpoint.

### Per-stage increments, against the rho decomposition

| stage removed | Δ rho (from the decomposition) | Δ raw hits | Δ d ≥ 2 hits |
| --- | --- | --- | --- |
| Pareto (−0.121 → −0.104) | +0.017 | +0.417, t=+0.96, **ns** | +0.625, t=+1.85, **ns** |
| acquisition (−0.104 → −0.032) | +0.072 | +0.104, t=+0.19, **ns** | **−0.917, t=−2.59, w=0.0088, sign=0.0051, 10/28/10 — separable and NEGATIVE** |
| prior (−0.032 → +0.148) | +0.180 | **+4.562**, t=+3.97 | **+2.979**, t=+4.16 |

**Verdict: the decomposition agrees with the campaign only for its largest step, and has
the sign backwards for one stage on the endpoint that matters.**

- **Sign.** All three rho steps are positive (removing the stage improves rho). The
  campaign agrees for Pareto and the prior. For **acquisition it disagrees on the
  diversity-adjusted endpoint**: removing the acquisition function *costs* 0.917 distinct
  hits, separably (t=−2.59, Wilcoxon 0.0088, sign 0.0051, better on 10 of 38 non-tied
  seeds). Rho says removing acquisition should help second-most of the three; the endpoint
  says it hurts.
- **Order.** Rho orders the stages Pareto < acquisition < prior. The campaign orders them
  acquisition < Pareto < prior on raw, and acquisition < Pareto < prior on d ≥ 2 — the two
  small stages are **transposed** relative to rho on both endpoints, and on d ≥ 2 the
  transposition is separable because the acquisition rung is separably negative.
- **Magnitude.** Two of the three rho steps do not reach the endpoint at all. Pareto's
  +0.017 rho is worth 0.4 raw hits, which 48 seeds cannot separate from zero.

This was registered in advance (§6.3 of the pre-registration) as a **reportable finding,
not a failure**: it means the stage-wise rho decomposition is measuring something that does
not transfer to the campaign endpoint for its two smaller stages. It is the third
independent instance of that lesson in this project — the first being `biophysical_features`,
where surrogate rho gain and campaign gain ordered in *reverse* at Spearman −0.5.

### The residual test (registered §6.4)

Does the coin beat the deepest ladder rung? If so, there is harm beyond these three stages.

| contrast | raw | d ≥ 2 |
| --- | --- | --- |
| coin − rung 3 (`pareto,acquisition,prior`) | +0.354, t=+0.25, p=0.80, 27/21/0 — **ns** | +1.917, t=+2.19, w=0.037, **sign=0.461**, 26/20/2 |
| coin − rung 1 (`pareto`) | +5.021, t=+5.74, 34/10/4 | +3.979, t=+5.79, 37/9/2 |
| coin − rung 2 (`pareto,acquisition`) | +4.917, t=+4.42, 34/10/4 | +4.896, t=+6.58, 37/8/3 |

**No residual is established.** On raw, removing those three stages fully accounts for the
coin's advantage (+5.083 of +5.438). On d ≥ 2 a +1.917 gap remains with t=+2.26-grade
evidence but **a sign test of 0.461 — 26 of 46 non-tied seeds, a coin flip**, so it is
reported as not separable. Had it been separable it would have implicated the GP residual's
own sign on CreiLOV or `AssayAdaptedHead`, which sits downstream of `mu(x)` and which no
rung removes; that remains a live possibility, unmeasured.

### On GB1

This ladder was run on CreiLOV only, so it says **nothing** about GB1's gap. Noting for the
next session: GB1's loss now has four mechanisms eliminated — the surrogate's feature space
(one-hot closes 33.4% of the ridge gap and it remains a decisive loss on five of six
comparators), the spatial veto (ablating it *costs* 4.27 hits; the cohesion rule fires zero
GB1 rejections at any depth), the fabricated WT row (0.20 label SDs; 48-seed ablation
+0.542 at t=+0.41), and now **the ranking stack, which on the d ≥ 2 endpoint contributes
nothing separable at all (−0.604, t=−1.34)**. Against that, §3's observation that
`panoptes_default − ridge` on GB1 is **−12.896 raw but +0.042 at d ≥ 2 (t=+0.09, 18/19/11,
not separable)** suggests the gap being chased may be substantially an artifact of the raw
endpoint, with ridge's 4.83 redundancy against the engine's 2.94 doing the work. That is a
hypothesis this campaign did not set out to test and did not pre-register.

---

## 6. The wild-type row, separated

The 2 × 2, both arms in one folder per cell, one seed block, `ridge` bit-identical across
folders.

| | anchor present (run A) | `wt_anchor` ablated (run G) | Δ of ablating it |
| --- | --- | --- | --- |
| `panoptes_default` raw | 6.938 | **38.708** | **+31.771**, t=+35.88, 48/0/0, CI [+30.06, +33.50] |
| `panoptes_default` d ≥ 2 | 5.062 | **20.688** | **+15.625**, t=+24.04, 48/0/0, CI [+14.38, +16.88] |
| **coin** raw | 12.375 | 13.125 | **+0.750**, t=+0.93, p=0.356, sign 0.883, 24/22/2 — **not separable** |
| **coin** d ≥ 2 | 9.667 | 10.375 | **+0.708**, t=+1.15, p=0.257, sign 0.280, 25/17/6 — **not separable** |
| coin − default, **anchor present** | — | — | **+5.438** raw / **+4.604** d ≥ 2 (coin wins) |
| coin − default, **anchor ablated** | — | — | **−25.583** raw / **−10.312** d ≥ 2, **0/48/0** (coin loses) |

**The CreiLOV coin win is the wild-type row, not the ranking stack.** The interaction is a
swing of −31.021 raw / −14.916 at d ≥ 2 — the ranking stack is worth −5.4 hits when its
surrogate is fitted through a label 15.68 SDs out of range, and **+25.6 hits when it is
not**. Both endpoints give the same answer, so this is not a redundancy effect.

Two further things this measures:

- **It replicates the 12-seed figure at 48 seeds**: +31.771 against the reported +29.75,
  48 of 48 seeds, on an independent seed block from that measurement.
- **It bounds how much the surrogate reaches the ballot.** The arm was registered as
  isolating the ranking stack *and not the surrogate*, because `gfn_reward` scores GFlowNet
  trajectories through the same surrogate. The coin's response to removing the poisoned row
  is **+0.750 raw, not separable on either endpoint** — so the surrogate's influence on the
  ballot is below this campaign's resolution. The isolation is cleaner than the
  pre-registration was willing to claim.
- `wt_anchor` is **not** promoted to a default here. It remains an ablatable component.

---

## 7. The ballot's base rate — the registered disagreement, resolved

The pre-registration recorded two of this project's own figures that imply CreiLOV answers
an order of magnitude apart and cannot both describe the same ballot: the **0.81× ratio**
(coin ≈ 8.6 hits) and "**5.2× enriched over chance**" (coin ≈ 20 hits).

**Resolved: neither. The post-veto ballot's base rate on CreiLOV is 12.4× chance.**

Because the coin draws uniformly from that ballot, its committed precision *is* an unbiased
estimate of the ballot's base rate — no modelling needed: **12.375 hits per 100 selected
measurements = 12.4%, against a chance base rate of 1.000% (1,655 hits of 165,427).** The
"5.2×" figure therefore does not describe the post-veto ballot the ranking stack is handed.

The campaign-level analogue of 0.81× is **6.938 / 12.375 = 0.561**: the committed batch
scores 0.56× what a uniform draw from the engine's own ballot achieves. That is worse than
0.81, and the two are *different quantities* — 0.81 is per-round and frozen, 0.56 is an
endpoint ratio between two arms whose ballots diverge after round 1. Both say the same
thing about direction.

**`candidate_trace.csv` cannot be used for this.** `BoundedCandidateTrace` retains
`selected + top_50_by_quality + 25_sampled_rejections`, so its "ballot" over-represents
high-quality members. Computing a ballot base rate from it gave 7.90% for the default and
15.10% for the coin — both biased upward, and the coin's implied ratio of 0.864 is an
artifact of that retention policy, not a measurement. Quoted here only so the next session
does not repeat it.

---

## 8. Verification of the arm itself

- **The default path is provably untouched.** Against the published 48-seed CreiLOV sweep
  run `run_20261101_024544_e82e6214` (commit `406ef58`, same seed block), **15,120 cells
  over 720 (arm, seed, round) rows × 21 behavioural columns are IDENTICAL** for `random`,
  `ridge` and `panoptes_default`. A matching mean would have been weak evidence; this is
  cell-by-cell. `panoptes_default − ridge = −50.292, t=−32.01, 0/48/0` reproduces the
  published figure exactly.
- **The coin is uniform, measured on the real selection path.** Controlled probe on
  CreiLOV, 12 seeds, one fixed ballot per comparison, 300 committed picks per rule: mean
  selected depth / mean ballot depth = **1.0295 for the coin against 1.2293 for the ranking
  stack**. The ranking stack prefers deep candidates by 23%; the coin does not prefer
  anything.
- **`generation_coverage.csv`'s cross-arm depth ratio is NOT a uniformity test** and put
  the coin at 0.826. The two arms' ballots differ after round 1, and the file's two
  histograms accumulate per `propose()` call including backfill, so they are not the same
  population. Recorded as an instrument caveat so the figure is not mistaken for a defect
  in the arm.
- **Rule 2.12, non-vacuously.** `SelectedBatch.selection_rule_applied` is stamped inside
  the branch that chose the indices, and `PanoptesSelector.propose` refuses the round if it
  disagrees with the request. The check is proven non-vacuous by the fact that a caller
  asking for `"ranked"` is told **`"top_n_quality"`** — an echo of the argument could not
  produce that. Every manifest here records, under `results.notes`,
  `panoptes_default -> top_n_quality` and `panoptes_random_from_ballot -> random_from_ballot`.
  `configuration.components.disabled` correctly reads `['dkl','dpp']` from the module
  globals, and `['dkl','dpp','wt_anchor']` in run G.
- **25 unit tests**, mechanism-only, each behavioural one paired with a non-vacuity guard.
  The first fixture had a post-veto ballot of **four** against a batch of four, which makes
  every selection rule identical — both perturbation tests passed while measuring nothing,
  and the guards caught it.

---

## 9. Gate verdicts, including the failures

The gate was run **twice per run folder** — statically before the scored run (rule 2.7) and
again **with `--run-dir`**, without which `budget` and `threshold` report NOT-APPLICABLE.
`--corpus-dir` was pointed at the real repository root's `downloads\combingym`; a
`FAIL corpus.lock` inside a worktree is an environment artifact of
`scripts/freeze_cross_protein_corpus.py:55` resolving against its own `__file__`. Verified
independently: **CORPUS VERIFIED 29/29 at `c5dccd514a9a`** (rule 2.5), and the rule 2.2
genotype anchors reproduce — TrpB 168/2,242 = 7.4933%, GB1 72/2,167 = 3.3226%.

On this branch, all seven folders:

| check | verdict |
| --- | --- |
| `corpus.lock`, `corpus.derived_table`, `corpus.prepared` | **PASS** |
| `encoding` | NOT-APPLICABLE — no encoded-design-space arm was requested |
| `budget.<arm>` | **PASS for every arm**: all 48 seeds spent all 125 distinct measurements |
| `endpoint.<lib>` | **PASS** — creilov d≥2 ceiling 647/1655 (60.9% removed); trpb 166/1600 (89.6%); gb1 146/1494 (90.2%) |
| `threshold.<lib>` | **FAIL on all seven** |
| `row_order.<lib>` | NOT-APPLICABLE — no `--invariance-runner` supplied |

**The `threshold.<lib>` failure is the known reader/schema bug, and it is now fixed on
`main`.** `30c1e94` ("the threshold check read one key of four and failed correct runs") and
`017de36` landed after this branch was cut. Re-running the gate over these **same, unchanged
run artifacts** from a throwaway worktree at `017de36`:

> **All seven folders exit 0, and `threshold.<lib>` PASSES**: creilov "the bar is quantile
> 99.96, resolving to 4.165688884891826 on this landscape's own values, from the array that
> defines top1.threshold"; trpb 0.6998334526751873; gb1 5.362064688799873.

So: **six of six checks pass on every folder under the fixed checker, and the only failure
on this branch is a checker defect already repaired upstream.** The bar it concerns is
consumed solely by the `qualify_then_optimize` / `conditional_probe` arms, which this
campaign does not run, and the failure is symmetric across arms, so it could not bias any
paired contrast either way.

---

## 10. What this does not show

- **All three landscapes are development data.** TrpB and GB1 were already so under rule
  2.4; CreiLOV became so on `ec96417` ("CreiLOV is now development data"). Every figure
  here is therefore **mechanism evidence, not generalization evidence**, and this campaign
  adds no independent-protein validation.
- **The arm does not isolate the surrogate**, only the ranking stack. §6 bounds the
  surrogate's reach into the ballot at +0.750 hits, not separable — but that is a bound, not
  a zero.
- **The ladder is CreiLOV-only.** Nothing here attributes TrpB's or GB1's stage structure.
- **`random_forest`, `gp_onehot`, `knn_hamming`, `gbm_onehot` and `greedy_measured` were
  not run.** Comparator strength is therefore stated, not measured, in every contrast
  except against `ridge` and `random` — and on GB1 `ridge` is known to be the weaker
  reference.
- **d ≥ 2 is a sequence-level proxy.** Hamming separation is not evidence of functional or
  biophysical independence, and `endpoint.<lib>` reports that the adjustment removes 89.6%
  of TrpB's and 90.2% of GB1's hit set — so the d ≥ 2 counts are small and the two
  endpoints are answering genuinely different questions, not one correcting the other.
- **No default was changed.** `DEFAULT_METHODS`, `DPP_DIVERSITY_ENABLED`,
  `TWO_STAGE_GP_ENABLED`, `PRIOR_CALIBRATION_ENABLED` and the wild-type anchor row are all
  as they were.

## Run folders

All under the real repository root's gitignored `benchmark_comparison_results\`.

| tag | run folder | library | arms | ablated |
| --- | --- | --- | --- | --- |
| A | `run_20261101_162530_f42cda5d` | creilov | random, ridge, default, **coin** | — |
| B | `run_20261101_164652_e8ba9aef` | trpb | random, ridge, default, **coin** | — |
| C | `run_20261101_165622_1f358506` | gb1 | random, ridge, default, **coin** | — |
| G | `run_20261101_165937_6bb1cf9f` | creilov | ridge, default, **coin** | `wt_anchor` |
| D | `run_20261101_170655_71fd0bcd` | creilov | ridge, default | `pareto` |
| E | `run_20261101_171659_e7539f21` | creilov | ridge, default | `pareto,acquisition` |
| F | `run_20261101_172412_f60f4635` | creilov | ridge, default | `pareto,acquisition,prior` |
