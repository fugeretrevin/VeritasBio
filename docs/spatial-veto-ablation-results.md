# Result — the spatial viability veto is NOT the cause of the engine's losses

Read [`prereg-spatial-veto-ablation.md`](prereg-spatial-veto-ablation.md) first. It was
committed in `a62b782` before any ablated arm was run, and it fixed every threshold below.

**H-VETO — "the veto is the cause of the three measured losses against one-hot ridge" — is
FALSIFIED on all three landscapes.** On two of the three the veto is *net beneficial*:
removing it makes the engine measurably worse.

Campaign: 2026-11-01, 48 paired seeds per landscape, `--base-seed 0`, budgets matched per
arm, seed as the unit of independence, 20,000-sample whole-seed bootstrap, read out with
`scripts/analysis/veto_ablation_readout.py` (which imports `paired_stats`'s own functions
rather than deriving a second recipe). Metric `cum_top1_hits_excluding_initial` at the
final round. Harness `scripts/analysis/spatial_veto_ablation.ps1`.

## 1. Every published baseline reproduced bit-for-bit

Before anything else, because it is what makes the rest readable. The unablated arm and
ridge reproduce `docs/HANDOFF.md` §1 to four decimals, including the paired statistics:

| landscape | panoptes_default | ridge | Δ | t | W/L/T | HANDOFF §1 |
| --- | --- | --- | --- | --- | --- | --- |
| eqfp611_red | 0.0000 | 15.3542 | −15.3542 | −7.86 | 0/36/12 | 0 vs 15.35 ✓ |
| creilov | 6.9375 | 57.2292 | **−50.2917** | **−32.01** | **0/48/0** | −50.29, −32.01, 0/48/0 ✓ |
| gb1 | 20.3958 | 33.2917 | **−12.8958** | **−5.63** | **6/41/1** | −12.90, −5.63, 6/41/1 ✓ |

## 2. The paired table

| landscape | arm | mean | sd |
| --- | --- | --- | --- |
| **eqfp611_red** | ridge | 15.3542 | 13.5340 |
| | panoptes_default, veto ON | 0.0000 | 0.0000 |
| | panoptes_default, **veto ABLATED** | **4.2083** | 4.8332 |
| **creilov** | ridge | 57.2292 | 11.3742 |
| | panoptes_default, veto ON | 6.9375 | 5.3571 |
| | panoptes_default, **veto ABLATED** | **4.8958** | 5.1498 |
| **gb1** | ridge | 33.2917 | 13.8425 |
| | panoptes_default, veto ON | 20.3958 | 9.0218 |
| | panoptes_default, **veto ABLATED** | **16.1250** | 8.3834 |

| landscape | paired comparison | n | Δ | t | 95% boot CI | Wilcoxon p | sign p | W/L/T |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eqfp611_red | ablated − unablated | 48 | **+4.2083** | **+6.03** | [+2.896, +5.604] | 3.7e−06 | 7.45e−09 | 28/0/20 |
| | unablated − ridge | 48 | −15.3542 | −7.86 | [−19.146, −11.646] | 1.67e−07 | 2.91e−11 | 0/36/12 |
| | ablated − ridge | 48 | −11.1458 | −5.76 | [−14.875, −7.375] | 5.96e−06 | 1.51e−05 | 7/35/6 |
| creilov | ablated − unablated | 48 | **−2.0417** | **−2.33** | [−3.750, −0.333] | 0.0126 | 0.00661 | 13/32/3 |
| | unablated − ridge | 48 | −50.2917 | −32.01 | [−53.333, −47.188] | 1.61e−09 | 7.11e−15 | 0/48/0 |
| | ablated − ridge | 48 | −52.3333 | −34.49 | [−55.229, −49.354] | 1.61e−09 | 7.11e−15 | 0/48/0 |
| gb1 | ablated − unablated | 48 | **−4.2708** | **−3.22** | [−6.792, −1.688] | 0.00323 | 0.00794 | 14/33/1 |
| | unablated − ridge | 48 | −12.8958 | −5.63 | [−17.396, −8.500] | 3.12e−06 | 1.77e−07 | 6/41/1 |
| | ablated − ridge | 48 | −17.1667 | −7.91 | [−21.312, −12.938] | 9.57e−08 | 5.54e−06 | 8/39/1 |

## 3. Verdict per landscape, against the pre-registered bars

| landscape | gap, veto ON | gap, veto ABLATED | closed | bar was | **verdict** |
| --- | --- | --- | --- | --- | --- |
| eqfp611_red | −15.3542 | −11.1458 | **+4.2083 = +27.4%** | ≥ 50% | **CONTRIBUTING CAUSE** — not the cause |
| creilov | −50.2917 | −52.3333 | **−2.0417 = −4.1%** | ≥ 50% | **NOT A CAUSE** — veto net beneficial |
| gb1 | −12.8958 | −17.1667 | **−4.2708 = −33.1%** | ≥ 50% | **NOT A CAUSE** — veto net beneficial |

- **eqFP611-red: the veto is the complete cause of the ZERO, and a contributing cause of
  about a quarter of the loss.** The ceiling was exactly 0 of 82 hits and the arm scored
  exactly 0.0000 on all 48 seeds — arithmetic, not a measurement. Ablating it produces hits
  on 28 seeds and never costs a hit on any seed (28 better, **0 worse**, 20 tied at zero;
  sign p = 7.45e−09). But 4.21 against ridge's 15.35 leaves −11.15 at t = −5.76, so **73%
  of that loss is not the veto.**
- **CreiLOV and GB1: the veto is not the cause, and it was helping.** Both movements are
  negative with bootstrap intervals excluding zero. The gap to ridge *widens*.

The pre-registered claim that decided this was written before the campaign: **a ceiling
that is not binding cannot be the cause.** The veto let 250 of CreiLOV's 1,655 hits through
and the engine found 6.94 — 2.8% of its own ceiling; on GB1, 3.7%. Removing a ceiling the
engine was nowhere near could not help, and it hurt, because on those two landscapes the
veto *raises* the candidate pool's hit density (2.8× on CreiLOV, 1.8× on GB1). Ablating it
dilutes the pool. **The veto was doing the engine's ranking a favour it could not do for
itself.**

### Was the prediction falsified?

**No. Every pre-registered interval was hit**, and the predicted shape — red a contributing
cause, CreiLOV and GB1 not causes — is what happened.

| prediction | predicted | observed | |
| --- | --- | --- | --- |
| red, ablated mean | point 2.5, [0.5, 8.0] | 4.2083 | ✓ |
| red, Δ ablated−unablated | > 0, 0 seeds worse | +4.21, 0 worse | ✓ |
| red, Δ ablated−ridge | [−15.0, −7.0], t < −5 | −11.15, t −5.76 | ✓ |
| creilov, Δ ablated−unablated | [−5.0, +8.0] pt +1.0, \|t\| < 4 | −2.04, t −2.33 | ✓ |
| creilov, Δ ablated−ridge | [−55, −38], \|t\| > 15, ≥46 worse | −52.33, t −34.49, 48 worse | ✓ |
| gb1, Δ ablated−unablated | [−6.0, +6.0] pt −1.0, \|t\| < 4 | −4.27, t −3.22 | ✓ |
| gb1, Δ ablated−ridge | [−18, −7], ≥38 worse | −17.17, 39 worse | ✓ |
| red, seeds better ≥ 30 | ≥ 30 | **28** | ✗ |
| red, seeds worse vs ridge ≥ 40 | ≥ 40 | **35** | ✗ |

**Two sub-predictions missed, both seed counts on eqFP611-red, and both for the same
reason: I did not anticipate structural ties.** The unablated arm is pinned at exactly 0,
so 20 of 48 seeds tie at zero and 6 tie against ridge; ties are neither wins nor losses, so
a threshold stated on the count of wins was mis-specified for a floored arm. **The bars are
not being moved.** The primary criterion was `fraction_closed ≥ 0.50`, it was stated first,
and red fails it at 0.274 either way. For the record, the three-way vocabulary's
`seeds-better ≥ 30` clause also fails at 28 — so read literally, the pre-registration puts
red in neither "confirmed" nor "contributing" but in the residual bucket. I am reporting it
as a contributing cause on the strength of the fraction (27.4%) and a sign test at
p = 7.45e−09 with zero losing seeds, and flagging that the clause was badly drafted rather
than quietly satisfying it.

## 4. The veto also starves the arm of its budget, and ablation fixes that

A finding from the gate rather than from the endpoint. On eqFP611-red:

| run | state | `budget.panoptes_default` | mean unique measurements of 81 |
| --- | --- | --- | --- |
| veto ON | **PARTIAL**, 7 short batches, 2 DEGENERATE seeds | **FAIL** — 71 of 81 on the worst seed | 80.27 |
| veto ABLATED | SUCCESS | **PASS** — all 81 on every seed | 81.00 |

The unablated arm reproduces the published run's shortfall exactly (also 71 of 81). Per
HANDOFF §5.8 the unablated red arm is `VALID_BUDGET_LIMITED`, so strictly it should not be
mean-compared against a fully budgeted arm. **It does not change the direction.** The mean
shortfall is 0.73 measurements; even ridge, the better ranker, earns 0.24 hits per selected
measurement here, so 0.73 extra measurements buys at most ~0.18 hits against an observed
+4.21. CreiLOV and GB1 spent all 125 on every seed in both arms.

## 5. Pre-flight gate (rule 2.7) — run on all six runs, and what it says

| check | eqfp611_red ON / ABL | creilov ON / ABL | gb1 ON / ABL |
| --- | --- | --- | --- |
| corpus.lock | PASS / PASS | PASS / PASS | PASS / PASS |
| corpus.derived_table | PASS / PASS | PASS / PASS | PASS / PASS |
| corpus.prepared | PASS / PASS | PASS / PASS | PASS / PASS |
| budget.ridge | PASS / PASS | PASS / PASS | PASS / PASS |
| budget.panoptes_default | **FAIL** / PASS | PASS / PASS | PASS / PASS |
| endpoint.`<lib>` | PASS / PASS | PASS / PASS | PASS / PASS |
| threshold.`<lib>` | **FAIL** / **FAIL** | **FAIL** / **FAIL** | **FAIL** / **FAIL** |
| encoding, row_order | NOT-APPLICABLE | NOT-APPLICABLE | NOT-APPLICABLE |

**`threshold` fails on every run, both arms, all three landscapes — and it fails
identically on the three runs that produced the published HANDOFF §1 table.** Verified by
re-running the gate against `run_20261101_024225_8fa1d86d`,
`run_20261101_024544_e82e6214` and `run_20261101_025014_5e79e3ae`. The reason is
`qualification.quantile is not recorded: the bar is an absolute constant … source is not
recorded`. It applies to a qualification bar that only `panoptes_qualify_then_optimize` and
`panoptes_conditional_probe` consume — neither arm here uses one — and it is symmetric
across the two arms, so it cannot bias the ablated-vs-unablated contrast. **But it is
un-passed, and HANDOFF §1's three loss figures carry it too.** Nobody had seen it because
`generality_sweep.ps1` runs the gate *without* `--run-dir`, which makes `threshold` and
`budget` NOT-APPLICABLE. That is a gap in the sweep harness, reported here, not fixed here.

### corpus.lock was failing for an environment reason, now fixed

On the first pass the gate reported `FAIL corpus.lock` on all three landscapes — "a corpus
that is not the one the protocol was frozen against", which reads as rule-2.5 drift.
It was not. `freeze_cross_protein_corpus.py` derived its default corpus directory from its
own `__file__`, i.e. the **checkout**; `downloads/` is gitignored, so from a worktree all 29
pinned files reported `absent`. Running the same verifier against the real root's mirror
gives **`CORPUS VERIFIED: 29/29 files match CombinGym @ c5dccd514a9a`**. Fixed by resolving
it through the module that owns root resolution (`veritas.paths.corpus_dir`, `search=
"repository"`), with four regression tests in `tests/test_path_resolution.py` — including
one on the script's own `DEFAULT_CORPUS_DIR`, because the gate passes that value and a
correct `paths.corpus_dir` with a stale script default would still fail. Paths checked:
`…\.claude\worktrees\agent-a038e3b42e920227a\downloads\combingym` (absent) and
`…\VeritasBio\downloads\combingym` (present, 29/29).

## 6. "The veto is wrong here" vs "the veto is harmful in principle"

Kept apart as the pre-registration required, because they have different evidence.

**The veto is WRONG on CreiLOV, independently of any campaign.** `8QI9` is
phototropin/CrPhotLOV1: 11 of 15 assayed positions map, and **all 11 disagree with the
library's own parent residue** (29:L≠A, 34:L≠G, 47:Y≠Q, 60:L≠R, 61:Q≠D, 92:K≠I, …).
Positions 3, 4, 5, 7 are absent, so pairs touching them drop out of the cohesion mean,
weakening the veto further. Production passes no `reference_sequence`, so
`SpatialContactMap`'s own 15% identity guard (`spatial_filter.py:221-243`) never fires. The
veto removes 94.6% of CreiLOV's universe **using a different protein's distances**. So
CreiLOV's result — that ablation costs 2.04 hits — says a filter reading unrelated geometry
was still, net, useful. That is a statement about pool density, not about biophysics, and it
must not be read as "the biophysical constraint helps".

**Whether the veto is harmful in principle is NOT answered by this campaign.** The only
landscape here with every assayed position mapped and every residue identity matching is
**GB1** — and on GB1 the cohesion rule contributes **zero** rejections at any depth (the
whole library fits inside ~8.2 Å). All 118,766 GB1 rejections are single-site buried-core
and steric calls firing at depth 1, e.g. `V39D` → "Electrostatic desolvation penalty". So
GB1 tests the chemistry rules on correct geometry; it does not test the depth-monotone
cohesion veto at all. eqFP611-red's `3M24` is the closest to clean (11 of 13 mapped, **one**
identity mismatch, 174:I≠A) and 98.4% of its rejections are cohesion — which is why its
0-of-82 ceiling is the strongest of the three results.

**TrpB (`8EH0`, 62.6% ceiling, correctly mapped) is the only landscape that could test the
cohesion veto on correct geometry, and it was not in this campaign.** That is the gap.

## 7. What the honest default should be

**Leave the veto ON by default. Do not disable it, and do not present it as validated.**

- Disabling it by default would be a **scientific claim** that distal co-mutations are not
  less likely to be jointly viable. This campaign supplies no evidence for that claim and
  two measurements against it: ablation cost 2.04 hits on CreiLOV and 4.27 on GB1, both with
  intervals excluding zero.
- Those two measurements do **not** vindicate the constraint either. On CreiLOV the veto is
  computing distances on the wrong protein and still helps, which means its benefit there is
  **pool-density enrichment, not biophysics**. A random filter that happened to enrich hit
  density 2.8× would have scored the same. That is not a defence of a biophysical prior.
- What *is* established and should change: **the veto must not be the only path to a
  candidate, and it must not be able to empty a landscape's reachable hit set.** On
  eqFP611-red its ceiling is literally 0 of 82 — the arm was structurally incapable of
  scoring, reported for a day as a performance loss. A veto whose reachable-hit ceiling is
  zero should be a **hard refusal to run that arm**, not a silent zero. `spatial_veto_profile.py`
  already computes the number; nothing consumes it.
- The two structure-mapping defects should be fixed before the veto is argued about again.
  Passing `reference_sequence` and an explicit `chain_id` from the benchmark call site would
  turn both into loud failures instead of plausible-looking wrong distances; the guard
  already exists and is simply never given its input.

Follow-on that would actually settle the principle, in priority order:

1. **Ablate the veto on TrpB** — correct geometry, 62.6% ceiling, cohesion live, and the one
   landscape the engine wins on. This is the missing cell.
2. **Split the ablation**: cohesion separately from buried-core/steric. They are one flag and
   three unrelated rules with opposite depth profiles, and GB1 proves a landscape can be
   dominated by the rules nobody is arguing about.
3. **Fix `1N8Z` chain selection and `8QI9`'s identity**, then re-measure. Until then every
   HER2 and CreiLOV structural number describes the wrong molecule.

## 8. Artifacts

| landscape | veto ON | veto ABLATED |
| --- | --- | --- |
| eqfp611_red | `run_20261101_141026_9fa07d9d` | `run_20261101_141145_be885c2f` |
| creilov | `run_20261101_141315_70988c06` | `run_20261101_141741_c089fb20` |
| gb1 | `run_20261101_142214_850df815` | `run_20261101_142330_2f8d06a9` |

In the **checkout's** gitignored `benchmark_comparison_results/` — the benchmark writes
beside the sources it ran, which for a worktree is the worktree. Campaign log
`vetoablation_20261101_091023.log`.

Pairing and labelling were verified before any number was read, and the read-out refuses
rather than warns if they fail:

- `ridge` is **bit-identical across each landscape's two runs** (48/48 seeds), so the two
  engine arms saw the same initial samples and the cross-run pairing is sound.
- `components.states.spatial_filter` is `true` in every unablated manifest and `false` in
  every ablated one; `requested_spec` is `"spatial_filter"`; configuration fingerprints
  differ between the arms (e.g. red `c092fc3ef7f2756e` vs `128473252ccf8aaf`).
- Rule 2.12, checked against the module and not the flag: `workers/math_worker.py` has no
  module-level boolean global mentioning SPATIAL (only `DPP_DIVERSITY_ENABLED`,
  `TWO_STAGE_GP_ENABLED`, `PRIOR_CALIBRATION_ENABLED`), so the single gate is
  `if "spatial_filter" in ablated:` at line 2233 and the manifest's CLI-derived state is the
  whole story. Verified empirically on eqFP611-red's own `3M24` before the campaign:
  unablated reaches the veto (1 call), ablated does not (0 calls) — not a vacuous pass.

**One un-gated second call site**, found while tracing the switch and reported rather than
changed: `math_worker.py:1029`, inside `generate_hamming1_lead_neighbors`, calls the veto
and is **not** covered by `--ablate-components spatial_filter`. It does not affect these
results — that function is called only from `controller/orchestrator.py:1353`, never from
the benchmark path (`docs/HANDOFF-2026-10-31-reachability.md` records it as "implemented
correctly and never called"). It would matter the moment the production pipeline is
ablated, and `ABLATABLE_COMPONENTS` currently promises more coverage than it delivers.
