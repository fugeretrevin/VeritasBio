# Results — the combined head-to-head: nine arms, one seed block, both endpoints

**Measured 2026-11-01.** Branch `worktree-agent-aa53d65b7e6a10b8f`, commit `fd9361a` (clean
tree, `dirty_tree: false`). Pre-registered in
[`prereg-combined-head-to-head.md`](prereg-combined-head-to-head.md), committed at `fd9361a`
**before any arm ran**.

| landscape | run directory | seed block | state |
| --- | --- | --- | --- |
| **TrpB** | `run_20261101_163432_c6faa21c` | **48–95, the designated HELD-OUT block** | SUCCESS |
| **GB1** | `run_20261101_164651_45ecc6b2` | 0–47, **not** a held-out block | SUCCESS |

Nine arms, one run folder per landscape, 48 paired seeds, budget 125 = 25 initial + 4 × 25,
one shared starting sample per seed. Artifacts are in the **gitignored** results root at the
real repository root: `C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results`.

Both landscapes are **development data** under rule 2.4. Everything here is mechanism
evidence, never generalization evidence.

---

## 1. The splice, and what it was worth

Two 48-seed TrpB results existed on **different seed blocks** and so could not be combined:
the seven-arm panel on seeds 48–95 (`panoptes_default` 51.42, strongest baseline `gp_onehot`
38.54) and the one-hot surrogate arm on seeds 0–47 (`panoptes_default` 49.65,
`panoptes_onehot_surrogate` 57.35). The naive splice gives
51.42 + 7.71 − 38.54 = **+20.59**. Measured directly on one block it is **+17.88**.

So the spliced figure was **2.7 hits too high**, inside the ±6 band the pre-registration
declared in advance (F3). The splice got the sign and the order of magnitude right and the
number wrong — which is exactly what cannot be known without measuring, and is why the
withdrawn 6.3× claim is not rehabilitated by this result.

The two premise numbers also disagree about `panoptes_onehot_surrogate` itself: **57.35 on
seeds 0–47 against 56.42 on seeds 48–95**, and about the one-hot *effect*: **+7.71 on
seeds 0–47 against +5.00 on seeds 48–95**. The seed block moves the headline arm by ~0.9
hits and its effect size by ~2.7.

## 2. Reproduction — every premise figure reproduces exactly

Measured here, on one run, against what was previously published on two:

| figure | previously | here | agree |
| --- | --- | --- | --- |
| TrpB `panoptes_default` (seeds 48–95) | 51.42 | **51.42** | exact |
| TrpB `gp_onehot` | 38.54 | **38.54** | exact |
| TrpB `panoptes_default − gp_onehot` | +12.88, t=+3.57, 31/16/1 | **+12.875, t=+3.57, 31/16/1** | exact |
| TrpB `panoptes_default − ridge` at d≥2 | +7.75, t=+13.4 | **+7.750, t=+13.40, 47/1/0** | exact |
| TrpB `gp_onehot` redundancy / `ridge` redundancy | 7.06× / 4.54× | **7.06× / 4.54×** | exact |
| GB1 `panoptes_default` | 20.40 | **20.40** | exact |
| GB1 `panoptes_default − ridge` | −12.90, t=−5.63, 6/41/1 | **−12.896, t=−5.63, 6/41/1** | exact |
| GB1 `panoptes_onehot_surrogate` | 24.71 | **24.71** | exact |
| GB1 one-hot effect | +4.31, 34/12/2 | **+4.312, t=+3.05, 34/12/2** | exact |

Validity conditions (§5 of the pre-registration), all met:

- **The default path is byte-identical after the merge.** `main` (`014e3cb`) and the merged
  tree (`fd9361a`) were run on TrpB seeds 48–95 with `random, ridge, random_forest,
  panoptes_default`: **21,120 cells across all 22 behavioural columns** of
  `round_metrics.csv` — 4 arms × 48 seeds × 5 rounds — agree **exactly**, the only exempt
  columns being `round_seconds` and `cumulative_seconds`. `seed_block_sha256` identical
  (`5f9cfea438cc7dd1…`), `source_tree_sha256` **different** (`2c2b7c11…` → `0ed5a9bb…`), so
  the code genuinely changed and the behaviour did not.
- **Adding five arms perturbs nothing.** The same 21,120 cells agree between that 4-arm run
  and the 9-arm TrpB campaign.
- **Every arm implementation is byte-identical to the branch its published figures came
  from** (AST-level source hashes): the five panel selectors to `493826b`, `PanoptesSelector`
  and `_one_hot_matrix` to `a779422`, `RandomSelector`/`RidgeSelector`/
  `RandomForestSelector`/`_rank_select` to `main`.
- **`_one_hot_matrix`'s refactor is a bitwise no-op.** The merged version (which now delegates
  to `workers.math_worker.one_hot_encode_genotypes`) was compared against `main`'s inline loop
  as an independent oracle on the **full** universes: **12,799,920 TrpB cells and 11,948,800
  GB1 cells bitwise identical**, same width 80, same `float32`, every row carrying exactly
  `n_sites` ones. This matters because all five one-hot baselines and `greedy_measured`'s
  Hamming identity read that matrix and were measured on the inline version.
- **`DEFAULT_METHODS` unchanged**, pinned by `tests/test_baseline_panel_arms.py`.
- **Rule 2.2 parser gate**: TrpB shallow shell 168/2,242 = 7.4933%, GB1 72/2,167 = 3.3226%.
- **Corpus lock 29/29** at CombinGym `c5dccd514a9a`.

## 3. TrpB — held-out seeds 48–95, all nine arms

Raw is `cum_top1_hits_excluding_initial`; `redund` is raw hits per diversity-adjusted hit.
Every arm spent all 125 distinct measurements on all 48 seeds (gate-verified).

| arm | raw | d≥2 | redund (d≥2) | d≥3 | best fitness | role |
| --- | --- | --- | --- | --- | --- | --- |
| `random` | 0.98 | 0.98 | 1.00× | 0.94 | 0.4118 | baseline |
| `ridge` | 29.06 | **6.40** | 4.54× | **2.06** | 0.6993 | registered reference |
| `random_forest` | 26.46 | 3.92 | 6.76× | 1.52 | 0.6620 | baseline |
| `gbm_onehot` | 27.21 | 4.79 | 5.68× | 1.75 | 0.6865 | baseline |
| `knn_hamming` | 21.06 | 5.44 | 3.87× | 2.00 | 0.6459 | baseline |
| `gp_onehot` | **38.54** | 5.46 | 7.06× | 1.77 | 0.6807 | baseline |
| `greedy_measured` | 21.69 | 2.60 | 8.33× | 1.12 | 0.6082 | baseline |
| `panoptes_default` | 51.42 | 14.15 | 3.63× | 4.50 | 0.7341 | engine |
| **`panoptes_onehot_surrogate`** | **56.42** | **14.25** | 3.96× | **4.79** | **0.7455** | engine |

**The strongest baseline is a different arm under each endpoint** — `gp_onehot` on raw,
`ridge` at d≥2 and d≥3 — exactly as pre-registered (F5). `gp_onehot` buys its raw lead with
redundancy: 7.06 raw hits per distinct one against ridge's 4.54.

### Engine with the one-hot surrogate, against the strongest baseline under each endpoint

| endpoint | strongest baseline | it scores | engine+one-hot | Δ | t | b/w/t | 95% bootstrap CI | Wilcoxon p | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **raw** | `gp_onehot` | 38.54 | 56.42 | **+17.88** | +5.08 | 37/11/0 | [+11.37, +24.88] | 2.7e−05 | **WIN** |
| **d≥2** | `ridge` | 6.40 | 14.25 | **+7.85** | +15.41 | **48/0/0** | [+6.88, +8.85] | 1.6e−09 | **WIN** |
| d≥3 | `ridge` | 2.06 | 4.79 | **+2.73** | +10.75 | 43/2/3 | [+2.25, +3.21] | 9.1e−09 | **WIN** |

And `panoptes_default` for comparison: raw **+12.88** (t=+3.57, 31/16/1,
[+6.19, +20.02]); d≥2 **+7.75** (t=+13.40, 47/1/0, [+6.62, +8.88]); d≥3 **+2.44**
(t=+12.08, 45/1/2, [+2.04, +2.83]).

The engine beats **all seven** baselines on all three endpoints, every contrast separable.
Its weakest raw contrast is against `gp_onehot` (+17.88) and its weakest d≥2 contrast against
`ridge` (+7.85).

### The one-hot surrogate's gain on TrpB is REDUNDANT hits

| contrast | raw | d≥2 | d≥3 |
| --- | --- | --- | --- |
| `panoptes_onehot_surrogate − panoptes_default` | **+5.00**, t=+5.33, 40/6/2, [+3.19, +6.83] | **+0.10**, t=+0.29, 22/19/7, [−0.58, +0.79] — **NOT SEPARABLE** | +0.29, t=+1.55, 21/12/15 — not separable |

This is a finding the splice could never have produced, and it is the sharpest result here.
The one-hot surrogate adds five raw top-1% hits per campaign on TrpB at t=+5.33 and **adds
nothing at all once hits within Hamming 1 of an already-counted hit are removed** — the
engine's redundancy rises 3.63× → 3.96× and the diversity-adjusted count is flat at 14.15 →
14.25. On the endpoint the pre-registration treats as co-primary, **the two engine arms are
indistinguishable on TrpB**. Any claim that the one-hot surrogate improves the engine must be
scoped to the raw count.

## 4. GB1 — seeds 0–47, all nine arms

| arm | raw | d≥2 | redund (d≥2) | d≥3 | best fitness | role |
| --- | --- | --- | --- | --- | --- | --- |
| `random` | 1.08 | 1.08 | 1.00× | 1.00 | 3.1738 | baseline |
| `ridge` | 33.29 | **6.90** | 4.83× | **2.29** | 6.5194 | registered reference |
| `random_forest` | 35.98 | 4.92 | 7.32× | 2.00 | 6.0657 | baseline |
| `gbm_onehot` | 31.40 | 5.27 | 5.96× | 1.88 | 6.0115 | baseline |
| `knn_hamming` | 27.90 | 6.08 | 4.59× | 2.21 | 6.0392 | baseline |
| `gp_onehot` | **46.02** | 6.12 | 7.51× | 2.00 | 6.8323 | baseline |
| `greedy_measured` | 31.17 | 3.90 | 8.00× | 1.58 | 6.0382 | baseline |
| `panoptes_default` | 20.40 | 6.94 | 2.94× | 2.96 | 5.5197 | engine |
| **`panoptes_onehot_surrogate`** | **24.71** | **7.69** | 3.21× | **3.04** | 5.6259 | engine |

**`ridge` is NOT the strongest baseline on GB1 raw. `gp_onehot` is, at 46.02 — 12.7 hits
above ridge.** The panel never ran GB1, so every GB1 figure in `docs/HANDOFF.md` is a
contrast against the project's third-strongest raw arm. The loss is substantially worse than
recorded.

| endpoint | strongest baseline | it scores | engine+one-hot | Δ | t | b/w/t | 95% bootstrap CI | Wilcoxon p | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **raw** | `gp_onehot` | 46.02 | 24.71 | **−21.31** | −9.10 | 7/41/0 | [−25.79, −16.62] | 3.3e−08 | **LOSS** |
| **d≥2** | `ridge` | 6.90 | 7.69 | **+0.79** | +1.81 | 24/16/8 | [−0.06, +1.65] | 0.109 | **NOT SEPARABLE** |
| d≥3 | `ridge` | 2.29 | 3.04 | **+0.75** | +4.23 | 28/7/13 | [+0.42, +1.10] | 2.5e−04 | WIN |

`panoptes_default`: raw **−25.62** (t=−11.39, 2/45/1); d≥2 **+0.04** (t=+0.09, 18/19/11) —
not separable; d≥3 **+0.67** (t=+3.44, 26/8/14) — win.

### Does one-hot close GB1's gap? No — it closes a third of it, as pre-registered

Against `ridge`, the comparator the open question was framed in:

| arm | GB1 raw − ridge | t | b/w/t |
| --- | --- | --- | --- |
| `panoptes_default` | −12.896 | −5.63 | 6/41/1 |
| `panoptes_onehot_surrogate` | **−8.583** | −3.66 | 13/34/1 |

**33.4% of the gap closed** — the pre-registration predicted 33%, 80% interval [10%, 60%].
It remains a decisive loss on every comparator: against ridge −8.58, against `random_forest`
−11.27, against `gp_onehot` −21.31, against `gbm_onehot` −6.69, against `greedy_measured`
−6.46; only against `knn_hamming` is it not separable (−3.19, t=−1.80), and it beats only
`random`.

**What that means**, as the pre-registration committed in advance: the GB1 loss is **not a
feature-space defect**, and together with the two existing results it is not the spatial veto
(ablating the veto *costs* GB1 4.27 hits, and the cohesion rule contributes zero GB1
rejections at any depth) and not the fabricated `y = 0.0` wild-type training row (GB1 sits
0.20 label SDs from centre, one of only two landscapes where that defect is arithmetically
invisible). With the surrogate's feature space, the veto and the WT row all eliminated, the
residual GB1 loss is attributable only to **acquisition and selection** — the stage TrpB
diagnostics already name as the bottleneck. That is a **narrowing, not a diagnosis**, and the
GB1 loss remains unexplained.

**Caveat carried from the pre-registration:** GB1's one-hot effect is **confounded with
candidate-pool composition** — the GFlowNet's reward scores candidates through this same GP,
and the one-hot arm proposed 7,833 distinct candidates against the default arm's 8,832 at
0.801 deep-enrichment against 0.690. GB1's +4.31 must not be attributed to the surrogate
alone. TrpB's was measured as not pool-explained.

### GB1's loss is RAW-ONLY, and that is a new finding

At d≥2 the engine is **at parity with the strongest baseline** (+0.04 for the default arm,
+0.79 for the one-hot arm, both not separable) and at d≥3 it **wins** (+0.67 and +0.75). The
baselines' GB1 raw leads are heavily redundant: `gp_onehot` 7.51 raw hits per distinct one,
`random_forest` 7.32×, `greedy_measured` 8.00×, against the engine's 2.94–3.21×.

This does **not** rescue the engine. The raw loss is real, it is the pre-registered primary
endpoint alongside d≥2, and the comparator ranking moved *against* the engine when the panel
arms were added. What it establishes is that **"the engine loses on GB1" is endpoint-dependent
and must always be stated with its endpoint** — a −25.62 raw loss and a statistical tie at
d≥2 are the same 48 campaigns.

`ridge` is separable from `random` on both landscapes under all three endpoints (TrpB
+28.08 t=+12.38; GB1 +32.21 t=+16.11), so unlike CR9114 these are contrasts against real
comparators.

## 5. The pre-registration, scored

| # | condition | outcome |
| --- | --- | --- |
| **F1** | TrpB raw Δ ≤ 0 or not separable | **not triggered** — +17.88, t=+5.08. Predicted +20.6, 80% [+10, +31]: **inside** |
| **F2** | TrpB d≥2 Δ ≤ 0 or not separable | **not triggered** — +7.85, t=+15.41, 48/0/0. Predicted +6.0, 80% [0, +14]: **inside** |
| **F3** | TrpB Δ outside [+14.6, +26.6] | **not triggered** — +17.88 |
| **F4** | GB1 Δ ≥ 0 (one-hot closes the gap) | **not triggered** — −21.31. The directional prediction holds |
| **F5** | same strongest baseline under both TrpB endpoints | **not triggered** — `gp_onehot` raw, `ridge` at d≥2, as predicted; replicates on GB1 |
| **F6** | `panoptes_default` off by > ±3 from 51.42 (TrpB) or 20.40 (GB1) | **not triggered** — 51.42 and 20.40 exactly |

**One prediction was wrong and it is recorded as wrong.** §3C named `ridge` as GB1's strongest
raw baseline at 33.3 and predicted Δ = −8.6 with 80% interval [−17, −1]. The strongest raw
baseline is **`gp_onehot` at 46.02** and Δ is **−21.31**, outside that interval. The sign and
the "one-hot does not close the gap" claim hold; the magnitude does not, because the
comparator was mis-predicted. The −8.6 figure is exactly right *against ridge* (−8.58), which
is the trap: the prediction was correct about the quantity it had seen measured and wrong
about the quantity it registered.

Three predictions were not registered and are reported as post-hoc observations: that the
one-hot gain is redundancy-only on TrpB at d≥2; that GB1's loss inverts to parity at d≥2; and
that `gp_onehot` outranks `ridge` on GB1 raw.

## 6. Gate verdicts, including the failure

Pre-run gate (rule 2.7), **both landscapes**: 4 passed, **0 FAILED**, 4 not applicable.
Post-run gate **with `--run-dir`**, both landscapes: **13 passed, 1 FAILED, 10 not
applicable**.

- **PASS `budget.<arm>` for all nine arms on both landscapes**: every arm spent all 125
  distinct measurements on every one of 48 seeds. This is the check the paired tests depend
  on, and it is the check that reports NOT-APPLICABLE when the gate is run without
  `--run-dir`.
- **PASS** `corpus.lock`, `corpus.derived_table`, `corpus.prepared`, `endpoint.trpb`
  (d≥2 ceiling 166 of 1,600, 89.6% of the hit set removed), `endpoint.gb1` (146 of 1,494,
  90.2%).
- **FAIL `threshold.trpb` and `threshold.gb1`.** `qualification.quantile` and
  `qualification.source` are not recorded where the reader looks. This is the known, unfixed
  reader/schema mismatch: `veritas/preflight.py::QUALIFICATION_MANIFEST_PATHS` resolves at
  `configuration.components.qualify_threshold`, which is a **scalar**, while the manifest
  records `qualify_quantile` and `qualify_threshold_source` as **siblings**. It is not missing
  data. The bar it concerns is consumed only by the `qualify_then_optimize` /
  `conditional_probe` arms, neither of which ran here, and the failure is symmetric across
  arms, so it **cannot bias any paired contrast above**. It is nonetheless an **un-passed
  check that every figure in this document carries**, and it is not called a pass.
- NOT-APPLICABLE: `encoding.<arm>` × 9 (no arm normalises an encoded design space by column
  range) and `row_order.<lib>` (needs an invariance runner). A NOT-APPLICABLE is not a pass.

## 7. The sentence

> On **TrpB**, held-out seeds 48–95, the engine with the one-hot GP surrogate beats the
> strongest of seven baselines by **+17.9 top-1% hits** raw (56.4 vs `gp_onehot`'s 38.5,
> t=+5.08, better on 37 of 48 seeds, 95% CI [+11.4, +24.9]) and by **+7.9** diversity-adjusted
> at d≥2 (14.3 vs `ridge`'s 6.4, t=+15.41, better on **48 of 48** seeds, CI [+6.9, +8.9]); on
> **GB1** it **loses** by **−21.3** raw (24.7 vs `gp_onehot`'s 46.0, t=−9.10, worse on 41 of
> 48 seeds, CI [−25.8, −16.6]) while being **statistically indistinguishable** at d≥2
> (+0.8 vs `ridge`'s 6.9, t=+1.81, 24/16/8) — and both landscapes are development data, so
> this is mechanism evidence and not generalization evidence.

## 8. Not supported by this campaign

- Any generalization claim: TrpB and GB1 are both development data (rule 2.4).
- Any front-loading multiple. The withdrawn 6.3× is not re-derived here and must not be
  restated until it is, with a per-budget-best comparator and a single arm.
- Any claim that `panoptes_onehot_surrogate` is a *corrected* `panoptes_default`. It is a
  separate arm — the same engine with `biophysical_features` ablated — and on TrpB's
  co-primary d≥2 endpoint the two are not separable.
- Any attribution of GB1's one-hot effect to the surrogate alone (pool confound, §4).
- `48 campaigns on one library are not 48 independent replicates` of any biological claim;
  they are 48 paired replicates of an algorithmic contrast.
