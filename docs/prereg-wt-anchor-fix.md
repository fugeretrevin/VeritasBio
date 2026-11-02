# Pre-registration — repairing the wild-type anchor instead of deleting it

**Written and committed 2026-11-01, BEFORE any arm of this campaign ran and before any
result of it was read.** Branch `wt-anchor-measured-label`, based on `bb925c9` (where
`wt_anchor` became ablatable) merged cleanly with `main` at `30c1e94`.

Falsifiability is the point. Every number below is a commitment, and §7 says exactly what
would make each one wrong. A result that contradicts a prediction here is reported as a
falsification, not re-narrated.

> **OUTCOME POINTER, added 2026-11-01 after the campaign, the only edit to this file.**
> Nothing below is changed, removed or softened — no bar is moved and no prediction is
> restated. The scored result, prediction by prediction, is in
> [`docs/wt-anchor-fix-results.md`](wt-anchor-fix-results.md) §5. In short: §5a's CreiLOV
> prediction and §5b's and §5c's aggregates **held**; §5a's CR9114 and §5b's TrpB
> predictions are **falsified**; and the `err-SD` metric this document reasons from is
> **falsified as an explanatory metric** — CR9114 ranks 2nd of nine on it and is a tight
> null, GB1 ranks 7th and carries the only separable repair-specific gain. The
> label-scale theory in either form does not predict the effect.

---

## 1. The defect, and why ablation is not the fix

`workers/math_worker.py::select_next_micro_batch` injects one synthetic training row per
round for the wild type carrying the literal label **`y = 0.0`** whatever the assay's
scale is, at variance `1e-4` — which asserts that a fabricated row is **100× more certain
than a real unreplicated assay reading**, which the same function puts at `1e-2`. Every
round's surrogate is fitted through it.

Ablating the row (`--ablate-components wt_anchor`) is worth **+29.75 top-1% hits on
CreiLOV** (t = 12.17, better on 12 of 12 paired seeds, ridge held constant at 59.167;
`docs/HANDOFF.md` §1).

But **ablation deletes the anchor, and the anchor exists for a reason**: it pins the
surrogate at a known point so that round 1 is not fitting an unconstrained GP to 25
observations. Deleting it costs that constraint everywhere. The principled repair is to
anchor at the wild type's **measured** label.

## 2. Step 1 result — the wild type's measured label, per landscape

Established by `scripts/analysis/wt_label_availability.py` **before any fix was written**.
Three independent sources, each cross-checked against a fourth: the bracket the prepared
table's own `DMS_score_bin` column implies, since `docs/dataset-inventory.md` §7 records
the rule `DMS_score_bin = 1 iff DMS_score >= the wild-type value`.

| landscape | in the prepared table? | WT measured label | where it comes from | 0-SD | **err-SD** | WT pct | WT a hit? |
| --- | --- | --- | --- | --- | --- | --- | --- |
| creilov | no | **4.023229753** | upstream CombinGym `CreiLOV_clean.xlsx` `wt` row, `mean_log` | 15.68 | **16.74** | 87.99 | no |
| cr9114_h1 | no | **8.424622015** | in-repo `prepare_cr9114_h1.py::CR9114_GERMLINE_SCORE` | 14.14 | **13.18** | 19.18 | no |
| cr6261_h1 | no | **7.0** (assay FLOOR, censored) | upstream `bnAbs_CR6261_H1_clean.xlsx` `n_mut == 0` | 11.00 | **8.75** | 0.00 | no |
| cr6261_h9 | no | **7.0** (assay FLOOR, censored) | upstream `bnAbs_CR6261_H9_clean.xlsx` `n_mut == 0` | 8.35 | **6.72** | 0.00 | no |
| eqfp611_blue | no | **1.574628255** | upstream `eqFP611_clean.xlsx` `n_mut == 0`, `blue` | 1.19 | **4.25** | 99.84 | **HIT** |
| eqfp611_red | no | **0.085316563** | upstream `eqFP611_clean.xlsx` `n_mut == 0`, `red` | 0.80 | **0.29** | 23.62 | no |
| her2 | YES (literal `WT`) | **0.1** (detection FLOOR, censored) | the prepared CSV's own `WT` row | 0.56 | **0.19** | 0.00 | no |
| trpb | no | **0.408073925** | in-repo `workers/data/TRPB_TM9D8.json` provenance block | 0.35 | **6.81** | 99.32 | **HIT** |
| gb1 | YES (literal `WT`) | **1.0** (assay definition: parent = 1.0) | the prepared CSV's own `WT` row | 0.20 | **2.53** | 97.56 | no |

* `0-SD` = `|0 − label mean| / label sd` — how much of an **outlier** the fabricated row
  is among the labels. This is the ordering `docs/HANDOFF.md` records.
* `err-SD` = `|wt_label − 0| / label sd` — how **wrong** the fabricated label is about the
  point it claims to pin. **This is the quantity the repair corrects.**

**Three findings, and the second changes the experimental design.**

1. **A measured WT label exists on all nine landscapes** — but on only **two** is it in
   the prepared table. The repo convention drops the parent ("the parent is not a library
   member"), so seven of the nine need the pinned upstream workbook or an in-repo
   constant. On three of the nine (cr6261 ×2, her2) the value is a **censored floor**, not
   a resolved measurement, and is recorded as such.

2. **The two orderings disagree, and `TrpB IS NOT A NEGATIVE CONTROL` for this repair.**
   Under `0-SD` TrpB is 0.35 and reads as a null. Under `err-SD` it is **6.81**, the
   fifth-largest of nine, and the parent measures 0.408 against a top-1% threshold of
   0.3494 — **the TrpB parent is itself a hit, at the 99.32nd percentile.** Replacing
   `y = 0.0` with `y = 0.408` on TrpB does not nudge a near-true label; it tells the
   surrogate that the all-parent corner of the design space is excellent. The genuine
   negative controls for the measured-label repair are **eqFP611-red (0.29)** and
   **HER2 (0.19)**, and HER2 is un-gatable on informativeness.

3. **GB1 is not a null either** (2.53), though it is a smaller change than TrpB.

## 3. The three behaviours

Mutually exclusive, independently runnable, pairing within seed on one shared starting
sample:

| arm | how | anchor row | label | variance |
| --- | --- | --- | --- | --- |
| **zero** (today's default) | no flag | present | literal `0.0` | `1e-4` |
| **ablated** | `--ablate-components wt_anchor` | absent | — | — |
| **measured** (the repair) | `--wt-anchor-mode measured` | present | the library's recorded parent label | `1e-2` |

**The variance is `1e-2`, and it is a parity statement rather than a tuning knob.**
`1e-4` claims 100× the precision of a real unreplicated reading, which this module puts
at `1e-2`. Once the row carries a real measurement it is entitled to exactly a real
measurement's precision and no more. The replicate-derived variances for the landscapes
that publish them are all **below** `1e-2` (CreiLOV 1.76e-3 over three log-scale
replicates; CR9114 7.67e-3 from its reported SEM), so `1e-2` is the conservative choice:
it cannot over-trust the anchor, and an override is clamped at `1e-4` so the
over-confidence cannot be reintroduced by accident.

The measured anchor is also placed in **residual space** (`y − m(x)`) like every measured
row, with the prior **evaluated** at the wild type. The `zero` path hardcodes the
anchor's prior to `0.0`, which is only coincidentally right.

### 3a. Disclosure: the measured arm receives one real fitness value outside its budget

This must be stated wherever a result from it is quoted. The `measured` arm is handed the
parent's true label without spending a measurement on it. On four landscapes that label is
high — TrpB 99.32nd percentile **and a hit**, eqFP611-blue 99.84th **and a hit**, GB1
97.56th, CreiLOV 87.99th — so it is real information, not merely a corrected constant.
`ridge` and `random` get no equivalent.

Three things make this disclosable rather than disqualifying, and the disclosure is
required anyway:

* the `zero` arm already asserts a wild-type label without spending budget, so the arms
  are comparable to each other; what changes is only whether the asserted label is true;
* a prospective lab always measures its parent, so the information is legitimately
  available at campaign start — unlike the landscape quantile of §5.2 of the handoff;
* it is nonetheless **one measurement of asymmetric information against the baselines**,
  so a `measured`-arm win must be reported as *"the engine plus one free parent
  measurement, against ridge at the full budget"*, never as a bare win.

The clean fix is to charge the parent measurement to the budget (125 → 124 + 1). That is
**not** done here, because it would change the budget match and confound the contrast.
Recorded as the first item of §6's cost.

## 4. The campaign

* 48 paired seeds, `--base-seed 0`, seed as the unit of independence.
* One shared starting sample per seed, identical across arms; verified by
  `configuration.seeding.seed_block_sha256` being equal across the three run folders of a
  landscape **and** by `ridge` and `random` being numerically identical across them.
* Budgets matched to the registered protocol: creilov / gb1 / cr9114_h1 / trpb all at
  `--initial-size 25 --batch-size 25 --rounds 4` = **125**.
* Three run folders per landscape, one per anchor mode, each carrying
  `ridge,random,panoptes_default`.
* Endpoint: **top-1% hits excluding the initial sample**, the endpoint every figure in
  `docs/HANDOFF.md` §1 uses. Thresholds read from each run's own `manifest.json`
  (rule 2.3), never recomputed.
* Statistics: paired over seeds — Δ, paired t, seeds better / worse / tied, and a
  10,000-sample whole-seed bootstrap 95% CI. **A difference of means is not a verdict**
  (rule 2.11).
* Pre-flight gate run **with `--run-dir`** after each run, because without it `threshold`
  and `budget` report NOT-APPLICABLE. Failures reported as failures.
* Landscape priority: **creilov**, then **trpb** and **gb1**, then **cr9114_h1**.

## 5. Predictions

Δ is in **top-1% hits excluding the initial sample**, at 48 paired seeds. Intervals are
the range I expect the bootstrap 95% CI to lie within.

### 5a. Contrast 1 — `measured − zero` (the repair against today's default)

| landscape | point | interval | confidence it is separable from 0 | reasoning |
| --- | --- | --- | --- | --- |
| **creilov** | **+30** | **+20 to +38** | **95%** | err-SD 16.74. Deleting the row is already worth +29.75 at 12 seeds; replacing it removes the same pathology and adds a valid observation, so the repair should recover at least as much. |
| **cr9114_h1** | **+6** | **+1 to +13** | 75% | err-SD 13.18, and the correct anchor (19.18th pct) points the surrogate AWAY from the parent corner — which is right, since CR9114's mean hit depth is 8.91. |
| **trpb** | **+2**, sign uncertain | **−6 to +10** | **70% that `\|Δ\| > 2` and p < 0.05** | **Pre-registered as NOT a null.** err-SD 6.81 and the anchor flips from a near-mean pseudo-observation to a top-0.7% value at the shallow corner, on the one landscape where shallow bias has already been measured as harmful (round 1 selected 0 of 193 depth-4 candidates). The magnitude is predicted; the SIGN is not. |
| **gb1** | **+2** | **−3 to +7** | 50% | err-SD 2.53, WT at the 97.56th pct but below the 2.152 threshold. A smaller change than TrpB's. |

### 5b. Contrast 2 — `measured − ablated` (the contrast that decides what ships)

This is the one that matters, so the reasoning is stated before the numbers.

The `ablated` arm has 25 real observations in round 1. The `measured` arm has those 25
**plus one more real observation, at the unique depth-0 corner of the feature space** —
the point every mutation string in the library is defined relative to. One extra correct
point out of 26 is a small effect in general. But the wild type is not a general point:
it is the only variant at depth 0, it is maximally distant from the sampled region in
Hamming terms, and on four of these landscapes its label is in the top 12% of the
distribution.

So **I predict the measured anchor is BETTER than no anchor** — ablation *understated*
the repair, and the fix is worth more than +29.75 — with the effect concentrated where
the parent is informative.

| landscape | point | interval | confidence Δ > 0 | confidence separable |
| --- | --- | --- | --- | --- |
| **creilov** | **+2** | **−2 to +6** | 70% | 40% |
| **trpb** | **+4** | **+1 to +9** | **80%** | 60% |
| **gb1** | **+2** | **0 to +6** | 70% | 45% |
| **cr9114_h1** | **+1** | **−1 to +5** | 60% | 30% |

**Aggregate primary prediction:** `measured ≥ ablated` on **at least 3 of the 4**
landscapes, and strictly better (p < 0.05) on **at least 1**.

### 5c. The three possible answers, and which I predict

| if | then | my probability |
| --- | --- | --- |
| measured **≈** ablated everywhere | the anchor is useless and should simply go — the simpler change, and it needs no parent label, which is the fix's whole provenance burden | 30% |
| measured **>** ablated | ablation understated the fix; the repair is worth more than +29.75 and the anchor should be **repaired, not removed** | **55%** |
| measured **<** ablated | anchoring at all is harmful. The mechanism would be LEVERAGE: one point at the far corner of the feature space, at a variance 1e-2 against the GP's fitted noise, drags the length-scale and pulls predictions toward the parent corner — which on TrpB is exactly the shallow bias already measured as harmful. If this is what happens, expect it to show up on TrpB and eqFP611-blue (the two landscapes where the parent is a hit) before CreiLOV. | 15% |

### 5d. The identity prediction

The flag-absent default path is **byte-identical**. Against the parent commit, over
`random, ridge, random_forest, panoptes_default` × 48 seeds × 5 rounds, every behavioural
column of `round_metrics.csv` agrees exactly, differing only in wall clock. Predicted with
**99%** confidence, because the `zero` branch is the previous code verbatim including the
order of its appends, and the new arguments are omitted entirely from the engine call on
that arm rather than passed as `"zero"`.

## 6. What this costs if it ships — written in advance so it cannot be minimised later

Changing the default **re-bases every number in this repository**. The following would
have to be re-measured, and are listed now so the list cannot be trimmed after the result
is known:

1. Every TrpB figure in `docs/HANDOFF.md` §1 "What survives and can be quoted" —
   `panoptes_default` 51.42 raw, `+7.75 at d≥2 (t=+13.4)`, and `+12.88 (t=+3.57, 31/16/1)`
   against the seven-arm baseline panel. **This is the repository's only surviving win**,
   and §2 of this document predicts TrpB moves.
2. The whole cross-protein sweep table (CreiLOV −50.29, GB1 −12.90, eqFP611-red 0 hits,
   eqFP611-blue +0.73, CR9114 +3.17, HER2 gate-failed).
3. The seven-arm baseline panel on TrpB, CR9114-H1 and CreiLOV
   (`worktree-agent-a8325f2a02fbcdc90` @ `493826b`).
4. The `biophysical_features` one-hot arm (`worktree-agent-a14c69a80e0afa44a` @ `a779422`).
5. The spatial-veto ablation (`worktree-agent-a038e3b42e920227a` @ `3804237`).
6. The `wt_anchor` ablation itself, including the 48-seed confirmation on
   `confirm-wt-anchor-48`.
7. Every per-depth and per-stage diagnostic computed from a `panoptes_default` run:
   `generation_coverage.json`, `depth_returns.csv`, `candidate_trace.csv`,
   `ranking_decomposition.py` and `pool_hit_reachability.py` outputs.
8. The budget curves in `docs/preregistration-cross-protein.md` §63 and §1434, which are
   **frozen and must not be edited** — a withdrawal would have to be cited alongside them.

And two new obligations the default change would create:

9. **A provenance obligation on seven of the nine landscapes.** The parent label is not in
   those prepared tables. A default that reads it makes the gitignored upstream corpus a
   runtime dependency of the engine's default path, which it is not today. Shipping would
   require the nine values to be committed into the prepared CSVs or into a committed
   registry with digests — not read from `downloads/`.
10. **The §3a disclosure becomes permanent**, and the budget should then be charged
    (125 → 1 parent + 124), which is itself a protocol amendment.

## 7. Falsification conditions

Stated as tests, with the consequence named.

* **F1.** `measured − zero` on CreiLOV has a bootstrap 95% CI containing 0.
  → §5a is falsified and the label-scale theory of the CreiLOV loss does not survive;
  whatever the `wt_anchor` ablation measured, it was not the label.
* **F2.** `measured − ablated` is **negative and separable** (p < 0.05) on **two or more**
  landscapes. → §5c's middle row is falsified, anchoring at all is harmful, and the
  recommendation becomes **remove the anchor**, not repair it.
* **F3.** TrpB moves on **neither** contrast (`|Δ| < 1.0` and p > 0.1 on both).
  → §2's `err-SD` reframe is falsified, the `0-SD`/centre framing was the right one, and
  TrpB really is the negative control `docs/HANDOFF.md` treats it as.
* **F4.** `ridge` or `random` differ at all between the three run folders of one
  landscape, or `seed_block_sha256` differs. → the arms are not paired; that landscape's
  contrasts are **void**.
* **F5.** Any scored run folder whose `run.wt_anchor_verified.agrees_with_configuration`
  is false, or whose `vacuous` is true for an arm that called the engine.
  → that folder is **void** under rule 2.12; a manifest that misdescribes its own
  component state is worse than no manifest.
* **F6.** Any behavioural cell of `round_metrics.csv` differs on the flag-absent default
  arm against the parent commit. → **the default moved, and the entire campaign is void.**
* **F7.** The `measured` arm's observed anchor label differs from the library's recorded
  `parent_measured_fitness` by more than 1e-12 in any round. → the plumbing is wrong and
  the arm is not measuring its own label. (Enforced at the call site: it raises.)

## 8. What is NOT done here, deliberately

* **The default is not changed.** `DEFAULT_WT_ANCHOR_MODE` stays `"zero"`, pinned by
  `tests/test_wt_anchor_measured_mode.py` so the one-line edit that would flip it
  silently fails CI.
* **No measured WT label is invented.** `--wt-anchor-mode measured` refuses to start on a
  library with none, rather than falling back to `0.0` — which would be the same defect
  wearing the fix's name.
* **`f_best = 0.0` is left alone.** It is the same literal-zero defect in the EI
  reference point, on the assay's scale, and it is identical in all three modes. Changing
  it here would confound the anchor contrast. Reported as a finding, not fixed.
* **eqFP611-red is not run**, despite being the one clean negative control, because it was
  not in the brief's priority order and a 48-seed run there is the obvious next
  experiment. Named so the omission is visible.
