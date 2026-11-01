# Handoff — 2026-10-31 — candidate reachability

Read [`docs/HANDOFF.md`](HANDOFF.md) first for the standing rules. This file covers one
session and supersedes nothing in it. The supporting analysis is
[`docs/day1-reachability.md`](day1-reachability.md), which carries the full depth tables,
the provenance review and the void-run record.

---

## 1. One-line state

The candidate generator was hardcoded to one library geometry — four saturated sites — and
was therefore blind on seven of nine landscapes. That is fixed and measured. The engine now
beats one-hot ridge on two landscapes, loses badly on a third, and the surrogate itself is
shown to be **at parity with ridge**, which redirects the next step away from selection
math and onto feature representation.

## 2. Quotable results

48 paired seeds, `--base-seed 0`, 4 rounds. `top1` is `mean_top1_hits_excluding_initial`.

> **VERDICTS CORRECTED 2026-11-01.** Every number below reproduced exactly in the full
> generality sweep (commit `406ef58`, 48 paired seeds). Two of the verdicts did not
> survive a paired test, which none of them had been given:
>
> | landscape | Δ vs ridge | t | seeds b/w/t | corrected verdict |
> |---|---|---|---|---|
> | CR9114-H1 | +3.17 | +6.10 | 41/5/2 | WIN stands — but development data (rule 2.4) |
> | **eqFP611-blue** | **+0.73** | **+0.77** | **25/18/5** | **NOT SEPARABLE — not a win** |
> | eqFP611-red | 0 vs 15.35 | — | — | loss stands |
>
> `eqFP611-blue` was the only landscape here offering generalization evidence, and it is
> a coin flip. The sweep also added CreiLOV (**−50.29, t=−32.01, 0 of 48 seeds better**)
> and GB1 (**−12.90, t=−5.63**), so the full picture is one development-data win, one
> tie and three losses. See `HANDOFF.md` §1 and rule 2.11.

| landscape | geometry | random | ridge | **panoptes** | verdict as first published |
|---|---|---|---|---|---|
| **CR9114-H1** | 16 sites, binary | 0.875 | 1.833 | **5.000** | **WIN, 2.7x** |
| **eqFP611-blue** | 13 sites, binary | 0.729 | 8.792 | **9.500** | ~~WIN~~ → not separable |
| **eqFP611-red** | 13 sites, binary | 0.688 | **15.354** | **0** | **LOSS, total** |
| TrpB | 4 sites, saturated | — | — | 49.6458 | unchanged (control) |

- CR9114: best fitness 9.70007 vs ridge 9.62781; **100 % of seeds found >= 1 hit** vs
  ridge's 66.7 %; one seed reached 9.83451, the **global optimum**
  (`S31N:T57S:T71S:K74I:S75F:T76S:S77N`, depth 7).
- eqFP611-blue: best fitness 1.58473 vs ridge 1.54926; 100 % of seeds vs 89.6 %.
- eqFP611-red: spent 80.27 of 81 measurements and found nothing. **Same library as blue**,
  same 13 binary sites, same budget — only the phenotype differs. See §5.
- TrpB is **bit-identical** to two prior default runs (49.6458 top-1 %, 0.7356 best). It is
  the negative control and it passed: these changes are provably a no-op on a saturated
  library.

**CR9114 is development data** under rule 2.4 — the intervention was designed against its
measured failure. eqFP611-blue and -red are the only generalization evidence here, and they
disagree.

## 3. What changed in the code

| # | change | location |
|---|---|---|
| 1 | `allowed_residues` — per-site allele restriction honoured by every construct | `math_worker.generate_candidate_variants` |
| 2 | depth bound `le=4` -> `le=64`; real limit applied downstream per library | `controller/prompts.py:100` |
| 3 | `_position_budget` — position slices scale with alphabet size, so a stage's *candidate* budget stays fixed | `math_worker` |
| 4 | **per-site restriction in the policy-gradient sampler** — the majority candidate source | `math_worker.sample_policy_gradient_candidates` |
| 5 | planner no longer pre-clips depth to the literal 4 | `benchmark_panoptes_isolated._AnnotationDrivenPlanner` |
| 6 | `_one_hot_matrix` uses the declared alphabets; width now equals `one_hot_width` | `benchmark_panoptes_isolated` |
| 7 | `allowed_aa_per_position` declared for the six non-saturated libraries | `LIBRARIES` registry |
| 8 | `--production-candidate-source {generate,enumerate}` | `benchmark_panoptes_isolated` |

Tests: **674 pass, 2 skipped**. New: `tests/test_generator_reachability.py` (12 tests).

Why these mattered, measured on CR9114: the arm went from **26 of 125** measurements and
**0** hits to **125 of 125** and 5.00. `le=4` alone put 647 of 651 top-1 % hits out of
reach at any budget; the hardcoded 19-residue alphabet made ~99.9 % of proposals sequences
the library does not contain.

**Fix 4 was the decisive one and the code audit missed it**, because the audit was scoped to
`generate_candidate_variants`. The run's own `candidate_trace.csv` showed the real split:
**6,506 of 7,415** candidates came from the `gflownet` sampler and 909 from the
combinatorial generator. After fixes 1–3 and 5, 84.85 % of the ballot still used impossible
residues. *Scope an audit to the measured dominant path, not the apparent one.*

## 4. Deep analysis: what will and will not work next

Run offline, no campaigns needed. Script:
`scripts/analysis/surrogate_headroom.py` (train on a realistic budget, predict on 3,000 random
held-out variants, 5 repeats, identical one-hot features for every model).

| landscape | budget | ridge | GP | RandomForest | GP − ridge |
|---|---|---|---|---|---|
| trpb | 125 | 0.182 | 0.185 | 0.174 | +0.002 |
| gb1 | 125 | 0.238 | **0.289** | 0.314 | +0.051 |
| her2 | 113 | 0.062 | 0.045 | −0.001 | −0.017 |
| creilov | 125 | 0.908 | 0.912 | 0.869 | +0.004 |
| eqfp611_blue | 81 | 0.701 | 0.694 | **0.780** | −0.007 |
| eqfp611_red | 81 | 0.457 | 0.458 | **0.587** | +0.001 |
| cr9114_h1 | 125 | 0.836 | 0.846 | **0.906** | +0.011 |
| cr6261_h1 | 19 | **0.790** | 0.731 | 0.768 | **−0.058** |
| cr6261_h9 | 19 | 0.824 | 0.824 | 0.777 | −0.000 |

**Conclusions that change the plan:**

1. **The GP is not a differentiator.** Parity with ridge on 7 of 9 landscapes, better on 1,
   worse on 1. So `--production-candidate-source enumerate` would buy *parity*, not a win.
   Do not build a product claim on the surrogate being better.
2. **RandomForest beats both on 3 of 9** — and the benchmark already ships `random_forest`
   as a baseline. Any future claim must be measured against it, not just ridge.
3. **The engine's FEATURES are the suspect, not its model class.** A GP on plain one-hot
   reaches rho = 0.458 held-out on eqFP611-red. The engine's GP, on its own
   physicochemical/spatial featurisation, reached rho = −0.0069 *in the actual run* on the
   same landscape. Same model family, same landscape, opposite result. This is consistent
   with `SequenceReadoutAdapter`'s default embedding being a hand-rolled 64-dim hash that
   the code itself documents as **not** an ESM representation.

**Next experiment, highest value:** add one-hot to the surrogate's feature set (or replace
the physicochemical encoding with it) and re-measure. That is a bounded change with a
pre-registered prediction: eqFP611-red should move off zero. If it does not, the problem is
the training distribution the engine selects for itself, not the features.

**Caveat on the restricted-range argument:** the proposed-candidate fitness spread on
eqFP611-red (sd 0.2581) is close to the library's (sd 0.2906), so attenuation is mild and
the in-run rho of −0.0069 cannot be explained away by range restriction. It is a real
signal.

## 5. The open puzzle

eqFP611 **blue and red are the same library** — identical 13 binary sites, identical
positions, identical budget. The engine wins blue (9.50 vs 8.79) and scores zero on red
(0 vs 15.35). Geometry cannot explain it. Candidate diagnosis is in
`day1-reachability.md` §: on red the generator proposed **301 of 8,191** variants whose
mean fitness (0.1940) is *below* the library mean (0.2316), containing **0 hits where
chance alone predicts 3.01**. Local expansion from random low-fitness leads does not reach
red's hit set, which sits at 6.1x the library mean.

Two pre-registered landscape predictors were computed before these results existed and
**both are falsified** (`scripts/analysis/landscape_predictors.py`):

- rho_1 (neighbour fitness autocorrelation) predicted local search would fail where it is
  low. eqFP611-red's rho_1 is **0.687**, among the highest; HER2 (0.031), TrpB (0.400) and
  GB1 (0.426) are far lower.
- additivity predicted we beat ridge where it is low. We beat ridge on CR9114, whose
  additivity is **0.650** — higher than eqFP611-red (0.409), TrpB (0.257), GB1 (0.194).

Do not substitute a third post-hoc predictor and present it as a finding. With 9
landscapes and CR6261-H1/H9 non-independent, effective n = 8 and a correlation needs
|r| > 0.71 for p < 0.05.

## 6. Not done

- **Generality sweep incomplete.** CreiLOV, GB1, CR9114 (fair-encoding re-run) and HER2
  were not run. Script: `scripts/analysis/generality_sweep.sh`, which now refuses to start if a
  benchmark is already running.
- **Enumerate mode never run.** Implemented, tested for wiring, zero campaign data.
- **The CR9114 win is confounded** across fixes 1–6, applied together. Attribution needs
  them switchable one at a time; only 8 is a flag today.
- Dead code at `math_worker.py:941-961` (quadruplet crossover, provably subsumed by lines
  901–918). Budget starvation: stages demand ~17,954 slots against `max_candidates=2500`,
  so the cluster stages get zero whenever a lead exists. No WT-reversion move, though
  `generate_hamming1_lead_neighbors` implements it correctly and is never called.
- `max(strategy.max_mutation_depth, 3)` at `math_worker.py:1178` silently raises the schema
  default of 2 to 3. Left alone deliberately: changing it moves TrpB and GB1.

## 7. Rules this session added

1. **A task-stop is not a process kill.** Stopping a background sweep left its shell alive;
   it ran for three hours alongside a replacement and voided ten run directories. Use
   `taskkill /F /T` on the shell PID, then verify zero `benchmark_panoptes` processes.
2. **No repository edits while a scored run executes** — including markdown. The harness
   hashes the source tree and voids a run whose hash moves mid-flight. It is right to.
3. **Read progress from the run you started**, by manifest, never from "the newest
   directory". During contention that was the other sweep's, and TrpB was reported as
   CreiLOV for an hour.
4. **Run IDs are UTC; the shell is local.** A directory named `_211552` was created at
   17:15 local.
5. **Scope a code audit to the measured dominant path.** `candidate_trace.csv` records
   which generator actually produced the ballot. Check it first.
6. **A test can pin a bug.** `test_library_registry` asserted
   `one_hot_width == n_sites * 20`, requiring every library to be saturated — the exact
   assumption that caused the defect.
7. **Pre-register landscape predictors before the deltas exist.** Both of ours were
   falsified; post-hoc, something would certainly have fit nine points.
