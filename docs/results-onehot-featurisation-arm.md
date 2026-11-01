# Results — the one-hot surrogate arm on CreiLOV, GB1 and TrpB

Answers [`prereg-onehot-featurisation-arm.md`](prereg-onehot-featurisation-arm.md), which was
committed at `ce9b358` **before any of these runs was read**. Code under test: `e6c11c7`
(the arm) and `098982f` (the sweep script); the sweep log records
`head: 098982fea5c3b1a6de9cf55a3960db51d69fb1b5` and an empty `git status`.

Branch `worktree-agent-a14c69a80e0afa44a`. Not merged.

---

## 0. Read this before quoting anything here

**The pre-flight gate FAILS on all three landscapes.** Run with `--run-dir` and the
corpus directory resolved against the real repository root, every landscape reports
**8 passed, 1 FAILED, 5 not applicable**, and the failure is `threshold.<library>` on
each. Rule 2.7 says a non-zero exit means do not trust the run, so these figures are
reported as **provisional**, not as gated.

What the failure is, precisely, because it is not what its message says:

* the message is "`qualification.quantile` is not recorded … `qualification.source` is not
  recorded";
* the manifest **does** record both, as
  `configuration.components.qualify_quantile = 99.96` and
  `configuration.components.qualify_threshold_source = 'quantile'`;
* `veritas/preflight.py`'s `QUALIFICATION_MANIFEST_PATHS` resolves at
  `configuration.components.qualify_threshold`, which is a **scalar** `4.1657`, and
  `qualification_from_manifest` wraps a scalar as
  `{"absolute": …, "quantile": None, "source": None}` — discarding the two sibling keys it
  then reports as missing. It only recognises a *mapping* named `qualification`.

So this is a reader/schema mismatch in the gate, it fails identically on the three runs
HANDOFF §1's published figures come from, and it is **not** caused by this arm. It is also
about the **qualification bar**, which no arm in these runs consumes: `--methods` names
`random, ridge, panoptes_default, panoptes_onehot_surrogate`, and the bar is read only by
`panoptes_qualify_then_optimize` / `panoptes_conditional_probe`. The scoring threshold is a
separate, recorded quantity (`results.top1.threshold_value`, the 99.0th percentile over the
eligible set). **None of that makes the check a pass, and it is not reported as one.**
Deliberately not fixed here: changing a gate check requires validating it against the real
artifact rather than the fixture beside it (rule 2.10), which is its own task.

The checks that did pass are the ones that bear on these numbers:

| check | result |
| --- | --- |
| `corpus.lock` | **PASS** on all three, with `--corpus-dir <real root>\downloads\combingym` |
| `corpus.derived_table`, `corpus.prepared` | PASS on all three |
| `budget.<arm>` for **all four arms** | PASS — every one of 48 seeds spent all 125 distinct measurements |
| `endpoint.<library>` | PASS — creilov 647/1655 (60.9% removed), gb1 146/1494 (90.2%), trpb 166/1600 (89.6%) |
| `threshold.<library>` | **FAIL** (above) |
| `encoding.*`, `row_order.*` | NOT-APPLICABLE, with the missing input named |

`corpus.lock` FAILED in the sweep's own in-line gate. That is the worktree artifact:
`scripts/freeze_cross_protein_corpus.py` derives `DEFAULT_CORPUS_DIR` from its own
`__file__`, so from this worktree it looks in
`…\.claude\worktrees\agent-a14c69a80e0afa44a\downloads\combingym`, which is gitignored and
absent, and all 29 pinned files report `absent`. The real directory is
`C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\downloads\combingym`, and pointing the
same verifier there passes 29/29. **Not corpus drift, and rule 2.5 is intact.**

---

## 1. The acceptance criterion: with the component absent, nothing moved

Three GB1 runs at 48 seeds, `--base-seed 0`, 25 + 4×25:

| run | code | arms | artifacts |
| --- | --- | --- | --- |
| **A** | parent commit `7415846`, i.e. before this work | `random, ridge, panoptes_default` | `…\scratchpad\baseline-7415846\benchmark_comparison_results\run_20261101_142511_24587b09` |
| **B** | `e6c11c7`, flag absent | same three | `…\agent-a14c69a80e0afa44a\benchmark_comparison_results\run_20261101_142627_4db47fab` |
| **C** | `098982f`, the arm present | the three **plus** `panoptes_onehot_surrogate` | `…\agent-a14c69a80e0afa44a\benchmark_comparison_results\run_20261101_143628_d91f53ab` |

**A ≡ B.** 15,840 cells across the 22 behavioural columns of `round_metrics.csv` —
every method, every one of 48 seeds, every one of 5 rounds — agree **exactly**. The only
two columns that differ are `round_seconds` and `cumulative_seconds`, which are wall clock.
Both runs carry the same shared-starting-sample digest `e85ea769c7863906`.

**B ≡ C.** Identical on the same 15,840 cells for the three shared arms, so adding the
fourth arm to a run does not perturb the others. Same sample digest.

**And the published figures reproduce.** `panoptes_default` in these runs gives
CreiLOV 6.9375 and GB1 20.3958 against ridge's 57.2292 and 33.2917, i.e. **−50.2917 on
0/48 seeds** and **−12.8958 on 6/41/1** — HANDOFF §1's −50.29 and −12.90 to four decimals.

The mechanism behind the identity, pinned by
`tests/test_surrogate_feature_arm.py::TestTheDefaultPathIsUntouched`: with the component
absent the GP is handed the **same array object** the biophysical path built, not a
numerically equal copy, and **no one-hot column is constructed at all**.

---

## 2. The 48-seed paired table

Primary endpoint `cum_top1_hits_excluding_initial` at each (method, seed)'s final round.
Seed is the unit of independence; every arm in a row shares one starting sample per seed.
95% CI is a 10,000-draw whole-seed bootstrap; `p` is Wilcoxon signed-rank; W/L/T counts
seeds. Read out with `scripts/analysis/paired_stats.py`.

### Arm means

| landscape | run | random | ridge | `panoptes_default` | `panoptes_onehot_surrogate` |
| --- | --- | --- | --- | --- | --- |
| CreiLOV | `run_20261101_142751_db6cc014` | 1.1667 | **57.2292** | 6.9375 | 7.7292 |
| GB1 | `run_20261101_143628_d91f53ab` | 1.0833 | **33.2917** | 20.3958 | 24.7083 |
| TrpB | `run_20261101_143938_7776a616` | 0.9375 | 31.1042 | 49.6458 | **57.3542** |

### `panoptes_onehot_surrogate` − `panoptes_default`

| landscape | Δ | 95% CI | Wilcoxon p | sign p | W/L/T | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| **CreiLOV** | **+0.7917** | [−0.292, +2.000] | 0.33 | 0.211 | 25/16/7 | **not separable** |
| **GB1** | **+4.3125** | [+1.562, +7.042] | **0.0031** | 0.0016 | **34/12/2** | **improvement** |
| **TrpB** | **+7.7083** | [+6.104, +9.333] | **2.6e-08** | 1.6e-09 | **41/3/4** | **large improvement** |

### `panoptes_onehot_surrogate` − `ridge`, with the default arm for reference

| landscape | Δ default−ridge | W/L/T | Δ one-hot−ridge | 95% CI | W/L/T |
| --- | --- | --- | --- | --- | --- |
| CreiLOV | −50.2917 | 0/48/0 | **−49.5000** | [−52.750, −46.167] | **0/48/0** |
| GB1 | −12.8958 | 6/41/1 | **−8.5833** | [−13.104, −4.083] | 13/34/1 |
| TrpB | +18.5417 | 40/7/1 | **+26.2500** | [+21.896, +30.750] | **48/0/0** |

**The CreiLOV gap does not close.** 49.5 hits, losing on every one of 48 seeds, against
50.3 for the default arm. On GB1 the arm recovers a third of the deficit and still loses on
34 of 48. On TrpB it beats ridge on 48 of 48.

---

## 3. Pre-registration verdict

**No falsification condition fired.** F1 (CreiLOV Δ vs ridge > −35, or ≥24/48 seeds
better): no, −49.50 at 0/48. F2 (CreiLOV Δ ≥ +10, p<0.01): no, +0.79 at p=0.33. F3
(GB1 Δ ≥ +5, p<0.01): **no, narrowly** — +4.31 at p=0.0031, so the p-value cleared and the
effect size missed by 0.69 hits. F4 (Δ positive and p<0.05 on **all three**): no — two of
three fired, CreiLOV did not. F5 (Δ ≤ −5, p<0.01 anywhere): no.

**But two of the three point predictions were wrong, and that is the real result.**

| landscape | predicted Δ vs default (90%) | observed | inside the interval? |
| --- | --- | --- | --- |
| CreiLOV | +1.0, [−3, +8] | +0.79 | yes — and the point estimate was accurate |
| GB1 | −1.0, [−6, +3] | **+4.31** | **no, and the sign was wrong** |
| TrpB | 0.0, [−5, +5] | **+7.71** | **no** |

The CreiLOV prediction, including "the gap to ridge stays near −49 and is not closed", was
right for the reason given. The GB1 and TrpB predictions were wrong, and the F3/F4
thresholds were set loosely enough relative to those point predictions that being wrong did
not trip them. **The honest reading is not "the pre-registration held"**; it is that the
§1.2 argument — a fixed-matrix rho advantage need not survive the engine's self-selected
training distribution — is **correct on CreiLOV and wrong on the two saturated libraries**,
and that F3 and F4 were mis-calibrated as tests of it.

---

## 4. Does better surrogate rho translate into a better campaign?

**Partly, and not in the order rho predicts.**

| landscape | Δ rho, fixed matrix (one-hot − engine) | Δ top-1% hits, campaign | p |
| --- | --- | --- | --- |
| GB1 | **+0.107** (0.306 vs 0.199) | +4.31 | 0.0031 |
| CreiLOV | +0.077 (0.914 vs 0.837) | +0.79 | 0.33 |
| TrpB | **+0.025** (0.174 vs 0.149) | **+7.71** | 2.6e-08 |

The rho ordering is GB1 > CreiLOV > TrpB. The campaign ordering is TrpB > GB1 > CreiLOV.
Over these three landscapes the two are, descriptively, **anti-ordered** (Spearman −0.5 on
n=3, which is a description and not an inference): the landscape with the *smallest*
fixed-matrix gain produced the *largest* campaign gain, and the landscape with the second
largest rho gain produced nothing measurable. **A fixed-matrix rho comparison did not
predict either the size or the existence of the campaign effect.**

Three measured reasons, in order of how much they explain:

**(a) Where the deficit is large enough, ranking cannot be the binding constraint.** On
CreiLOV the engine converts about 7 of ~100 selections into top-1% hits where ridge
converts ~57. No re-ranking of the ballot it holds can find 50 more hits. The arm is the
clean test of that: it hands the GP the baselines' own encoding, keeps the model class, the
prior, the acquisition, the Pareto stage and the generator fixed, and buys +0.79 at p=0.33.

**(b) The variance-allocation pathology is a function of library geometry, and it is
strongest exactly where the arm helped.** The redundant 10-dimensional global block is a
deterministic linear function of the rest, so its damage scales with the share of columns it
occupies under an isotropic kernel:

| landscape | sites | one-hot width | biophysical width | redundant block's share |
| --- | --- | --- | --- | --- |
| TrpB | 4 (saturated) | 80 | 30 | **33%** |
| GB1 | 4 (saturated) | 80 | 30 | **33%** |
| CreiLOV | 15 (narrow alphabets) | 35 | 85 | **12%** |

On the two saturated libraries one-hot is the *wider* encoding and replaces one in which a
third of the columns are redundant. On CreiLOV one-hot is the *narrower* encoding and the
redundant share is a third as large, so there is much less for the switch to repair. This
predicts the observed pattern and the rho protocol did not see it, because rho was measured
on a uniform random sample rather than on the engine's own selections.

**(c) Regime.** Where the engine is already competitive — TrpB ahead of ridge, GB1 within
13 hits — better ranking converts into hits. Where it is 50 behind, it does not.

### Confound, disclosed: the proposal distribution is downstream of the GP

The component is documented as moving "the GP surrogate's fit and predict inputs only",
which is literally true, but the **GFlowNet's reward scores candidates through that same
GP**, so the candidate pool is downstream of it. Measured from `generation_coverage.csv`:

| landscape | arm | distinct proposed | share of universe | mean proposed depth | deep enrichment |
| --- | --- | --- | --- | --- | --- |
| CreiLOV | default | 5,154 | 3.12% | 5.573 | 0 |
| CreiLOV | one-hot | 5,091 | 3.08% | 5.588 | 0 |
| GB1 | default | 8,832 | 5.91% | 3.464 | 0.690 |
| GB1 | one-hot | **7,833** | 5.24% | 3.553 | **0.801** |
| TrpB | default | 6,802 | 4.25% | 3.624 | 0.874 |
| TrpB | one-hot | **6,281** | 3.93% | 3.608 | **0.863** |

* **CreiLOV: no confound.** The pool is unchanged to within 1%, so the null is a clean
  ranking null and not a pool artifact.
* **GB1: confounded.** The arm proposed 11% fewer distinct candidates at 16% higher
  deep-enrichment. Pool composition moved in a direction that is independently known to
  matter — the spatial-veto study found it *helped* GB1 by raising pool hit density 1.8× —
  so **+4.31 must not be attributed to the surrogate alone.**
* **TrpB: not explained by pool density.** The pool got 8% narrower and slightly *less*
  deep-enriched while the result improved by 7.71 hits, so the gain runs against what pool
  density would predict and is attributable to ranking.

---

## 5. How this composes with the candidate-pool audit

On **CreiLOV** the −49.5 gap to one-hot ridge survives an arm that (i) is strictly better
on a fixed design matrix, (ii) uses the baselines' own encoder, (iii) holds every other
stage fixed, and (iv) leaves the candidate pool unchanged to within 1%. Together with the
spatial-veto result — ablating the veto *costs* CreiLOV 2.04 hits (t=−2.33), so the veto was
helping by raising pool hit density 2.8× — two of the live hypotheses are now eliminated:

* **the surrogate's feature basis** — ruled out here, with the pool held constant;
* **the spatial viability veto** — ruled out by the veto study, in the opposite direction
  to the one predicted.

What remains on CreiLOV, and what this result points at:

1. **The content of the candidate pool** — not its density, which both studies now
   implicate as *helping*, but what the generator can reach at all. CreiLOV is
   reachability-limited at 13.3%, the engine covers 3.1% of the universe at mean proposed
   depth 5.57, and its deep-enrichment is **0** on both arms. A candidate-pool depth audit
   speaks to exactly this, and this result strengthens its hypothesis by removing a
   competing one rather than merely failing to contradict it.
2. **The acquisition and the non-GP ranking terms** — the calibrated biophysical prior
   `m(x)` the posterior is built on, the Pareto filter, and greedy EI+Thompson. HANDOFF
   already records greedy beating EI+Thompson by +11.27 on TrpB and round 1 selecting zero
   of 193 offered depth-4 candidates.

On **GB1** and **TrpB** the surrogate basis is *not* ruled out — it is a real, significant
campaign effect there, confounded with pool composition on GB1 and not on TrpB.

**TrpB is a new arm, not a correction.** `panoptes_onehot_surrogate`'s TrpB figures stand
beside `panoptes_default`'s and must never be substituted for, averaged with, or described
as a revision of any published TrpB number. HANDOFF's 45.21 and +7.93 were measured on
`engine_biophysical`.

---

## 6. `test_regret_stopping` — measured, with a recommendation and no decision taken

### What it pins

`tests/test_regret_stopping.py::test_regret_early_stopping` asserts
`summary.iterations_completed < 10` and `final_remaining_usd > 0.50` for a mock-ESM
orchestrator run on P42212 with `regret_stopping_threshold=0.01`. The rule it exercises is
`max_pool_ei < regret_stopping_threshold` for `regret_patience` consecutive rounds from
round 3 on, where `max_pool_ei = max(EI(means, stds))` and `stds` comes straight from the GP
posterior. So the threshold is a comparison against an **absolute** EI, on a scale the GP's
kernel distances set.

### The measurement

Eight repeats per arm per case, each in a **fresh process** (numpy's global RNG state carries
across runs in one interpreter and the Thompson draw uses it). `select_next_micro_batch` is
wrapped in `controller.orchestrator` — where the name is bound — and the one-hot arm is
obtained by injecting `ablate={"biophysical_features"}`; nothing else differs.

| case | arm | iterations (min..max) | `iterations < max` | stopped on regret | budget assertion |
| --- | --- | --- | --- | --- | --- |
| early_stopping | `engine_biophysical` (today's default) | 8..10 | **5/8** | 2/8 | 8/8 |
| early_stopping | `one_hot` | 8..10 | **5/8** | 5/8 | 8/8 |
| pocket_hopping | `engine_biophysical` | 6..6 | 8/8 | 8/8 | 8/8 (hops 1/1, 8/8) |
| pocket_hopping | `one_hot` | 5..6 | 8/8 | 8/8 | 8/8 (hops 1/1, 8/8) |

Per-repeat iteration counts: biophysical `[8,8,9,10,9,10,9,10]`, one-hot `[8,10,10,8,10,9,8,9]`.

### Three findings

1. **`test_regret_early_stopping` is flaky on the completely unchanged default.** It hit the
   10-iteration cap in **3 of 8** repeats with no code change at all, so its assertion fails
   roughly a third of the time on today's `engine_biophysical` path. Both arms pass at
   **exactly the same rate, 5/8**. An observation that flipping the featurisation "broke
   `test_regret_stopping`" is therefore fully consistent with one unlucky draw and is not
   evidence about the featurisation.
2. **The featurisation does change a user-visible behaviour, not only a threshold the test
   happened to encode.** Round-1 `max_pool_ei` is the same to four decimals (median 1.43659
   vs 1.43632 — at round 1 the GP has almost no data and the prior dominates), but the
   final-round median collapses from **3.45e-03 to 4.00e-06**, an ~860× difference. The loop
   consequently reaches `stop_reason = converged_regret_threshold` in **5 of 8** one-hot runs
   against **2 of 8** biophysical runs: the campaign stops earlier and spends less budget.
   The config's absolute default of 0.005 means a different thing in each feature space.
3. **The orchestrator declares no alphabets, so one-hot there is very wide.** It calls
   `select_next_micro_batch` without `allowed_residues`, so every active position falls back
   to the canonical 20. Measured widths actually used across rounds: **60, 120, 160, 200 and
   240 columns**, against the biophysical **25–70** on the same rounds, with tens of training
   points. Nothing has measured whether that helps on an ESM-pseudo-likelihood objective.

### Recommended decision — not taken here

* **Leave the orchestrator on `engine_biophysical`.** It already is, and by taking the
  engine default rather than by a pin: this arm defaults off, so `controller/orchestrator.py`
  is untouched by this work and `test_regret_stopping` is unaffected by it. No pin is needed,
  and adding one would be a line asserting a default that is already the default.
* **De-flake `test_regret_early_stopping` before it is used as evidence about anything.**
  Either seed the loop so the EI trajectory is reproducible, or assert on the trajectory —
  that `max_pool_ei` falls below the threshold and stays there — rather than on a round count
  that varies by two rounds run to run. Until then a pass or a fail observed while changing a
  default is uninterpretable.
* **If the orchestrator's feature space is ever changed, make `regret_stopping_threshold`
  relative first** — a fraction of round-1 `max_pool_ei`, say — because an absolute EI
  threshold cannot transfer across feature spaces. This is the same objection rule 2.7's
  `threshold` check makes about the qualification bar: an absolute does not transfer, and
  0.005 is no more portable across feature spaces than 0.70 was across landscapes.

---

## 7. Reproduction

```powershell
# the arm, beside the default, in one run folder so the two pair within seed
powershell -File scripts\analysis\surrogate_arm_sweep.ps1 -Seeds 48 -Tag prereg

# the gate, WITH --run-dir and the corpus resolved against the REAL repository root
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library creilov `
  --budget 125 --initial-size 25 --run-dir <run_dir> `
  --corpus-dir <real repo root>\downloads\combingym `
  --arms random,ridge,panoptes_default,panoptes_onehot_surrogate <--acknowledge-derived ...>

# the readout
.venv\Scripts\python.exe -X utf8 scripts\analysis\paired_stats.py <run_dir> `
  --reference panoptes_default --metric cum_top1_hits_excluding_initial
```

Run folders, absolute, as `veritas.paths` reports them — note the benchmark's default
`--output-dir` is **checkout**-relative, so these landed under the worktree and not under
the real repository root:

* CreiLOV `C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\.claude\worktrees\agent-a14c69a80e0afa44a\benchmark_comparison_results\run_20261101_142751_db6cc014`
* GB1 `C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\.claude\worktrees\agent-a14c69a80e0afa44a\benchmark_comparison_results\run_20261101_143628_d91f53ab`
* TrpB `C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\.claude\worktrees\agent-a14c69a80e0afa44a\benchmark_comparison_results\run_20261101_143938_7776a616`
* identity B `…\benchmark_comparison_results\run_20261101_142627_4db47fab`
* identity A `C:\Users\stoof\AppData\Local\Temp\claude\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad\baseline-7415846\benchmark_comparison_results\run_20261101_142511_24587b09`
  — a detached `git worktree` at `7415846`, i.e. a temporary checkout; if it has been pruned
  the run folder is gone with it and run A must be reproduced from that commit.
