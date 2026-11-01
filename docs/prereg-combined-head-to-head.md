# Pre-registration — the combined head-to-head: engine + one-hot surrogate against a seven-arm baseline panel, one seed block

**Dated 2026-11-01. Written and committed BEFORE any arm of the combined campaign ran and
before any of its artifacts were read.** Branch `worktree-agent-aa53d65b7e6a10b8f`.

---

## 1. Why this exists: a splice that must not be made

Two 48-seed TrpB results exist and they are **on different seed blocks**, so they cannot be
combined:

| measurement | branch | seed block | figures |
| --- | --- | --- | --- |
| seven-arm baseline panel | `worktree-agent-a8325f2a02fbcdc90` (`493826b`) | TrpB **48–95** | `panoptes_default` 51.42; strongest baseline `gp_onehot` 38.54; Δ **+12.88**, t=+3.57, 31/16/1 |
| one-hot GP surrogate arm | `worktree-agent-a14c69a80e0afa44a` (`a779422`) | TrpB **0–47** (`--base-seed 0`) | `panoptes_default` 49.65; `panoptes_onehot_surrogate` **57.35**; Δ **+7.71**, p=2.6e−08, 41/3/4 |

Writing "engine + one-hot 57.35 beats best-baseline 38.54, so +18.8" would be a splice of two
seed blocks. **An exactly analogous splice produced the withdrawn "6.3× the best baseline"
claim** (`docs/HANDOFF.md` §1): it divided one arm's score by a comparator fixed at a
different budget, and the two ends of its quoted range used two *different engine arms*, so
"6.3× decaying to 1.7×" was never one curve. The correction is not a caveat — it is to
**measure the head-to-head directly**: all nine arms, one run folder per landscape, one
shared starting sample per seed, budgets matched, paired within seed.

Note also that the two premise numbers disagree about `panoptes_default` itself (51.42 on
seeds 48–95 against 49.65 on seeds 0–47, a 1.77-hit seed-block effect on the same arm at the
same commit-era defaults). That gap is the size of the error a splice silently absorbs, and
it is the reason no quantity below is formed by addition.

## 2. Design, fixed here

- **Arms (9), one run folder per landscape:** `random`, `ridge`, `random_forest`,
  `gbm_onehot`, `knn_hamming`, `gp_onehot`, `greedy_measured`, `panoptes_default`,
  `panoptes_onehot_surrogate`.
- **Landscapes, in priority order:** `trpb`, then `gb1`.
- **TrpB runs on the HELD-OUT block, seeds 48–95** (`--base-seed 48`, 48 seeds). TrpB is the
  one landscape with a designated held-out block — those seeds were never looked at during
  the 2^5 configuration scan — and it is the one landscape where the corrected qualification
  bar accepts the identical variant set, so nothing in the 2026-10-31 fix set changed what is
  searched there. **GB1 has no designated held-out block and runs `--base-seed 0`**, matching
  every GB1 figure already in `docs/HANDOFF.md`; its seed block is therefore *not* held out
  and that is stated wherever its numbers are.
- **Budget 125 = 25 initial + 4 × 25**, identical for every arm. Every arm takes its first
  batch from `RandomSelector(seed)` on the same pool, so round 0 is the same shared sample
  per seed for all nine arms. The read-out **checks** this (identical round-0 top-1% count
  per seed across arms) rather than assuming it, and checks per-arm
  `n_unique_measurements` for budget match.
- **Unit of independence: the seed.** 48 campaigns on one library are not 48 independent
  replicates of a biological claim (`HANDOFF` §5.6); they are 48 paired replicates of an
  algorithmic contrast, which is all that is claimed.
- **Both landscapes are DEVELOPMENT DATA under rule 2.4.** Nothing here is generalization
  evidence, whatever its sign or size.

### Endpoints — both reported, neither privileged

1. **Raw** `cum_top1_hits_excluding_initial` at the final round.
2. **Diversity-adjusted at d≥2** (and d≥3 for context): per seed, hits sorted by fitness
   descending, greedily keep one only when it is at least `d` Hamming away from everything
   already kept. Applied after measurement, so it cannot advantage a method that knew of it.
   Genotypes are expanded to full library length with wild type filled
   (`veritas.calibration.genotypes.canonical_combo` / `hamming`, which raises on unequal
   length). Rule 2.1: mutation strings omit wild-type positions; nothing is zipped.

**Both endpoints are pre-registered as co-primary because the panel already proved the answer
depends on which one is used:** on raw TrpB `gp_onehot` beats `ridge` by +9.48, but at d≥2 it
*loses* to ridge by −0.94 (t=−2.0), with redundancy 7.06× against ridge's 4.54×. "The
strongest baseline" is therefore endpoint-dependent, and a single-endpoint claim here would be
a choice of comparator disguised as a measurement.

### Statistics

For each endpoint: paired mean difference over the shared seeds, paired t, seeds
better/worse/tied, Wilcoxon signed-rank p, sign-test p, and a 10,000-sample whole-seed
bootstrap 95% CI. **A mean difference with no paired test is not a verdict** (rule 2.11 — the
most repeated error in this project; it published `eqFP611-blue` 9.50 vs 8.79 as a WIN at
t=0.77, 25/48 seeds). `|t| < 2.0` with no clear seed majority is reported as **NOT SEPARABLE**
and given no sign.

"Strongest baseline" means the highest-mean arm among the seven baselines **under the endpoint
being reported**, chosen within the same run. A maximum over seven arms is upward-biased for
the baseline by the winner's curse; that bias is declared here, it favours the baseline rather
than the engine, and every pairwise row is reported alongside so the maximum is never the only
number shown.

## 3. Predictions, with intervals

Point predictions with 80% intervals. These are predictions about the *measured* values on
**one** seed block, not restatements of the premise numbers.

### 3A. TrpB, raw endpoint

| quantity | prediction | 80% interval |
| --- | --- | --- |
| strongest raw baseline | `gp_onehot` | — (alternative: `random_forest` or `greedy_measured`) |
| its mean | 38.5 | [33, 44] |
| `panoptes_default` | 51.4 | [46, 57] |
| `panoptes_onehot_surrogate` | 59.1 | [52, 66] |
| **Δ = `panoptes_onehot_surrogate` − strongest baseline** | **+20.6** | **[+10, +31]** |

The point value is deliberately the naive spliced figure (51.42 + 7.71 − 38.54 = +20.59), so
that the measurement can contradict it. **Predicted verdict: WIN, t > 2, seed majority better.**

### 3B. TrpB, diversity-adjusted at d≥2

| quantity | prediction | 80% interval |
| --- | --- | --- |
| strongest baseline at d≥2 | **`ridge`**, *not* `gp_onehot` | — |
| **Δ = `panoptes_onehot_surrogate` − strongest baseline at d≥2** | **+6.0** | **[0, +14]** |

Rationale: `panoptes_default − ridge` at d≥2 was +7.75 (t=+13.4) on this block, and `gp_onehot`
falls below ridge once redundancy is removed. The interval is wider and its lower bound touches
zero because the one-hot surrogate's effect on *redundancy* has never been measured — it could
raise hit count by concentrating on one basin, in which case the raw gain shrinks or vanishes
here. **Predicted verdict: WIN, but less than the raw margin.**

### 3C. GB1, raw endpoint — the open question this settles

On GB1 the engine **loses** to ridge by −12.90 (t=−5.63, 6/41/1) and the one-hot surrogate
recovers +4.31 (p=0.0031, 34/12/2). The question is whether one-hot closes that gap.

**Prediction: it does NOT.** It closes roughly one third of it.

| quantity | prediction | 80% interval |
| --- | --- | --- |
| strongest raw baseline | `ridge` | — (alternative: `gbm_onehot` or `gp_onehot`) |
| its mean | 33.3 | [28, 38] |
| `panoptes_default` | 20.4 | [16, 25] |
| `panoptes_onehot_surrogate` | 24.7 | [20, 30] |
| **Δ = `panoptes_onehot_surrogate` − strongest baseline** | **−8.6** | **[−17, −1]** |
| fraction of the −12.90 ridge gap closed | 33% | [10%, 60%] |

**Predicted verdict: LOSS, t < −2.** At d≥2 the predicted Δ is **−5.0**, 80% interval
**[−13, +3]**, i.e. a loss that may become not-separable once ridge's own redundancy is removed;
no sign is predicted with confidence there.

### What it means if one-hot does NOT close the GB1 gap

It means **the GB1 loss is not a feature-space defect**, and the two candidate explanations on
the table are not the explanation either:

- it is not the spatial veto — ablating the veto *costs* GB1 4.27 hits (t=−3.22), so the veto
  was helping, and on GB1 the cohesion rule contributes **zero** rejections at any depth;
- it is not the fabricated `y = 0.0` wild-type training row — GB1 sits 0.20 label SDs from
  centre, which is one of the two landscapes where that defect is arithmetically invisible.

With the surrogate's feature space, the veto and the WT row all eliminated, the residual GB1
loss would be attributable only to **acquisition and selection** — the stage TrpB diagnostics
already name as the bottleneck (round 1 was offered 193 depth-4 candidates and selected zero;
greedy acquisition measures +11.27 against EI+Thompson). That is a *narrowing*, not an
explanation, and it must be reported as an open problem rather than a diagnosis.

If instead one-hot **does** close the gap (Δ ≥ 0 on GB1 at t > 2), the correct reading is not
"the engine is fine" but that the engine's biophysical featurisation was costing more than
12.9 hits on GB1 while costing nothing visible on TrpB — and then `biophysical_features`
stops being an ablation component and becomes a defect requiring a default change and a
re-basing of every published figure.

### 3D. The splice itself

**Prediction: the measured TrpB Δ lies within ±6 hits of the naive spliced +20.59.** If it
does, the splice's error was a seed-block mismatch of ordinary size and the direction of the
spliced claim happened to be right while its arithmetic was unlicensed. If the measured Δ
falls **outside** ±6, the splice was wrong in kind and not merely in precision — which is the
substantive finding this campaign exists to produce either way.

## 4. Falsification conditions

The pre-registration is **falsified** if any of the following is measured:

- **F1.** TrpB: `panoptes_onehot_surrogate` − strongest raw baseline is **≤ 0**, or is not
  separable (`|t| < 2`). (Predicted +20.6.)
- **F2.** TrpB at d≥2: that same contrast is **≤ 0** or not separable. (Predicted +6.0.) F1
  and F2 are separate: the panel already showed an arm can win on one endpoint and lose on the
  other, so a split outcome falsifies one and not the other and must be reported that way.
- **F3.** TrpB: the measured Δ falls outside **[+14.6, +26.6]** (±6 of the spliced +20.59).
- **F4.** GB1: `panoptes_onehot_surrogate` − strongest raw baseline is **≥ 0**, i.e. one-hot
  closes the gap. (Predicted −8.6.)
- **F5.** The strongest baseline is the **same arm under both endpoints on TrpB**, i.e. the
  endpoint does not change the comparator after all. (Predicted: `gp_onehot` raw, `ridge` at
  d≥2.)
- **F6.** Any arm's `panoptes_default` figure on TrpB seeds 48–95 departs from 51.42 by more
  than ±3 hits, or on GB1 seeds 0–47 from 20.40 by more than ±3 — which would mean the merge
  perturbed the default path and **the whole campaign is void**, not merely surprising.

A falsified prediction is reported as falsified. It is written here precisely so that it cannot
be re-narrated afterwards.

## 5. Validity conditions — the campaign is VOID if any fails

1. **The default path must be byte-identical after the merge.** The two branches both modify
   the harness and one modifies `workers/math_worker.py`; `a14c69a8` proved its own change
   harmless by comparing 15,840 cells across all 22 behavioural columns of
   `round_metrics.csv` against its parent commit over 3 arms × 48 seeds × 5 rounds. An
   equivalent identity check against `main` (`014e3cb`) is run on the merged tree **before**
   any campaign arm, and a perturbed default is a finding that voids the campaign.
2. **`DEFAULT_METHODS` unchanged**, pinned by `tests/test_baseline_panel_arms.py`. No arm here
   is promoted to a default and `ridge` remains the registered reference baseline; `HANDOFF`
   §5.4 makes promoting another model an amendment to a frozen hypothesis, which needs a dated
   amendment record and not a tuple edit.
3. **Rule 2.2 parser gate reproduced**: TrpB shallow shell 168/2,242 = 7.4933% and GB1
   72/2,167 = 3.3226%, depth 0 excluded from both shells and retained in the universe.
4. **Corpus lock verifies 29/29** at CombinGym `c5dccd514a9a`.
5. **Pre-flight gate run WITH `--run-dir`** after each run. Without it `threshold` and `budget`
   report NOT-APPLICABLE, which is how a `threshold.<lib>` failure stayed invisible on the
   published runs. A NOT-APPLICABLE is **not** a pass. `threshold.<lib>` is *expected* to fail
   on both landscapes — `veritas/preflight.py::QUALIFICATION_MANIFEST_PATHS` reads
   `configuration.components.qualify_threshold` as a mapping when it is a scalar, while the
   manifest records `qualify_quantile` and `qualify_threshold_source` as siblings. It is a
   known, unfixed reader/schema mismatch, symmetric across arms, so it cannot bias a paired
   contrast — and it is nonetheless an un-passed check that these figures carry, and it is
   stated wherever they are.
6. **No repository edit during a scored run.** The harness hashes the source tree before and
   after; a mid-run change voids the run. Everything in this campaign is committed first.
7. **Output filenames stamped.** Other sessions run scored campaigns in the same shared
   results root; a fixed filename there is a shared mutable name and has already truncated an
   artifact. A run directory without a `manifest.json` is a **partial launch, not data**.
8. **Budgets matched and pairing verified** per §2, reported as problems rather than assumed.

## 6. What this campaign cannot support

- Any generalization claim. TrpB and GB1 are both development data (rule 2.4).
- Any claim about `panoptes_onehot_surrogate` being a *corrected* `panoptes_default`. It is a
  **separate arm** — the same engine with `biophysical_features` ablated — and every figure in
  `docs/HANDOFF.md` was measured on the biophysical encoding.
- Any attribution of a GB1 difference to the surrogate **alone**. The GFlowNet's reward scores
  candidates *through* this same GP, so the proposal distribution is downstream of the feature
  space: on GB1 the one-hot arm proposed 7,833 distinct candidates against the default arm's
  8,832, at 0.801 deep-enrichment against 0.690. GB1's one-hot effect is confounded with
  candidate-pool composition. TrpB's was measured as not pool-explained.
- Any front-loading multiple. The withdrawn 6.3× must not be restated until it is re-derived
  with a per-budget-best comparator and a **single** arm, and nothing here does that.
