# Pre-registration — confirming the fabricated WT training row at 48 seeds

**Written and committed before any arm of this campaign ran.** Branch
`confirm-wt-anchor-48`, based on `bb925c9a3c740dc935728aae3072fc7f85a385d5`
(the commit carrying `wt_anchor`, `physchem_features` and
`featurize_variant_one_hot`), merged cleanly against `main` at `b5f090b`.

This document exists so the result cannot be re-narrated afterwards. HANDOFF §1
records that the project's previous flagship verdicts — eqFP611-blue as a "win",
GB1 as "parity" — were both produced by reading a mean difference without a
paired test, and that the one prediction which *was* written down in advance
(§8) is the only one that survived contact with the data.

---

## 1. The claim under test

`select_next_micro_batch` injects one synthetic training row per round for the
wild type, carrying the literal label `y = 0.0` whatever the assay's scale, at
variance `1e-4` while a real measurement defaults to `1e-2`. CreiLOV's real
labels occupy 2.793–4.196, so that row sits 11.6 units below the minimum: 8.3×
the entire label range, 15.68 label standard deviations from centre. Every
round's surrogate is fitted through it.

**The claim:** that row is the *dominant* cause of the engine's measured
−50.29-hit loss against one-hot ridge on CreiLOV.

Measured once, at **12** paired seeds: `--ablate-components wt_anchor` moves the
engine from 9.083 to 38.833 (+29.750, t = 12.17, better on 12 of 12 seeds),
which is 59% of the deficit. This campaign settles it at **48**.

"Dominant" is defined numerically, once, here: **Δ ≥ 25.04 hits**, i.e. at
least half of the 50.08-hit deficit the 12-seed campaign measured against its
own ridge constant of 59.167.

### Two distinct violations, deliberately not conflated

The row is wrong in two independent ways, and this campaign cannot separate
them because `wt_anchor` removes both at once:

| violation | CreiLOV | TrpB | GB1 |
| --- | --- | --- | --- |
| **label scale** — distance from 0.0 to the label centre, in label SDs | 15.68 | 0.35 | 0.20 |
| **variance** — asserted `1e-4` against a real measurement's `1e-2` | 100× over-confident | 100× | 100× |

The label violation is landscape-specific. **The variance violation is
identical on every landscape.** So a small non-zero effect on TrpB or GB1 is
*not* automatically evidence against the scale theory — it may be the variance
term alone. That is written down here, in advance, precisely so it cannot be
deployed afterwards as an unfalsifiable escape hatch: the falsification
thresholds in §4 and §5 are binding whatever mechanism is proposed.

---

## 2. Design

- **Metric:** `cum_top1_hits_excluding_initial` at each seed's final round.
- **Unit of independence:** the seed. 48 paired seeds, `--base-seed 0`.
- **Arms per landscape**, two runs because `--ablate-components` is global:
  1. `--methods ridge,panoptes_default` (unablated)
  2. `--methods ridge,panoptes_default --ablate-components wt_anchor`
- **Primary contrast:** Δ = `panoptes_default`(ablated) − `panoptes_default`(default), paired by seed.
- **Secondary contrast:** each engine arm − `ridge`, paired by seed.
- **Statistics (rule 2.11):** paired t over seeds, seeds better/worse/tied, and
  a 10,000-sample whole-seed bootstrap 95% interval. A mean difference alone is
  not a verdict.
- **Budgets:** 25 initial + 4 × 25 = 125 on creilov, trpb, gb1, cr9114_h1.
  Matched across arms within a landscape by construction.

### Pairing proof — a refusal, not a warning

`ridge` is independent of `--ablate-components`, so within each landscape it
**must be bit-identical across both arms on all 48 seeds**, and
`configuration.seeding.seed_block_sha256` must match. If either fails the
landscape is reported as **refused**, not as a result with a caveat.

### Gate

Rule 2.7, run **with `--run-dir`** (HANDOFF §1: `generality_sweep.ps1` omits it,
which makes `threshold` and `budget` report NOT-APPLICABLE and is why
`threshold.<lib>` was never seen to fail on the three published runs). From a
worktree the gate also needs `--corpus-dir <real root>\downloads\combingym`,
because `freeze_cross_protein_corpus.py:55` resolves the corpus against its own
`__file__` and `downloads/` is gitignored; without it all 29 pinned files report
`absent` and `corpus.lock` FAILS as an environment artifact rather than as
rule-2.5 drift. Verified before this document was committed: with the flag,
`CORPUS VERIFIED: 29/29 files match CombinGym @ c5dccd514a9a`.

**`threshold.<lib>` is expected to FAIL** on every run here, for the recorded
reason — `qualification.quantile` and `qualification.source` are absent from the
manifest. It is symmetric across arms and the bar it concerns is consumed only
by the `qualify_then_optimize`/`conditional_probe` arms, neither of which runs
here, so it cannot bias a paired contrast. It is nonetheless an un-passed check
and is reported as a failure, never as a pass.

---

## 3. Instrument audit — completed before this document was committed

Recorded here because an instrument audited after the fact is not an audit.

1. **The ablation removes the synthetic row *and* its variance entry.** Both
   `sample_variances_list.append(1e-4)` and the four row appends sit inside the
   same `if "wt_anchor" not in ablated:` block. This matters because
   `math_worker.py:2170` sets `sample_variances_arr = None` when the lengths
   disagree — a misalignment would silently disable heteroscedastic noise for
   the whole arm rather than raise.
2. **That test is not vacuous.** Verified by mutation: moving the variance
   append outside the guard makes
   `test_the_sample_variance_vector_stays_aligned` fail with
   `AssertionError: None != 2`, which is exactly the silent-`None` path. The
   mutation was reverted and the tree is clean. 24 of 24 tests pass across
   `test_wt_anchor_ablation.py`, `test_onehot_feature_ablation.py` and
   `test_spatial_filter_ablation.py`.
3. **Rule 2.12 — the manifest's component state is honest for this component.**
   `wt_anchor` has **no module-level global**; `component_states()` applies
   `globally_off` only to `dpp` (`DPP_DIVERSITY_ENABLED`) and `dkl`
   (`TWO_STAGE_GP_ENABLED`). So the `--ablate-components` set is the sole and
   authoritative gate, and `configuration.components.states.wt_anchor`
   correctly reads `True` on a default run. Verified in a real manifest.
4. **The flag-absent path is unchanged.** A 4-seed default CreiLOV run on this
   branch is metric-identical to the same run on `main`'s engine (`b5f090b`):
   Δ = 0.000, tied on 4 of 4 seeds, `ridge` bit-identical, `seed_block_sha256`
   identical. The `configuration_fingerprint` does differ, as it must — the
   ablation vocabulary gained two entries — which is a recorded configuration
   change, not a behaviour change. Ablating nothing changes nothing.

---

## 4. Predictions — CreiLOV, the confirmation

The 12-seed campaign gives mean +29.750 at t = 12.17 on n = 12, hence a
between-seed SD of 8.47 and a 12-seed standard error of 2.45. Carrying both the
48-seed sampling error (SE 1.22) and the uncertainty in the 12-seed estimate
itself gives a combined SE of ≈ 2.7. Seeds 0–11 are a subset of 0–47, which
shrinks this somewhat; the interval below is deliberately not shrunk for it.

| quantity | prediction |
| --- | --- |
| Δ (`wt_anchor` − default), 48 seeds | **+29.8, 95% interval [+24, +36]** |
| paired t | > 8 |
| seeds better / worse / tied | ≥ 44 / ≤ 2 / ≤ 4 |
| ablated arm vs `ridge` | ≈ **−20.3**, still a decisive loss |

**The engine is predicted to still lose CreiLOV after this fix.** 59% of a
50-hit deficit is not a win, and this document predicts the headline stays a
loss. Anything in the eventual write-up that reads as "the CreiLOV loss is
fixed" contradicts this line.

### Falsification of "the WT row is the dominant cause"

Falsified if **any** of:

- Δ < +25.04 hits (below half the 50.08-hit deficit), **or** the bootstrap 95%
  interval's lower bound falls below +25.04;
- paired t < 2;
- seeds-better < 32 of 48 (i.e. not better on two thirds).

Falsified in the *opposite* direction — the claim was understated — if
Δ > +45, which would make the row very nearly the whole deficit and would
mean the 12-seed estimate was badly biased low.

---

## 5. Predictions — TrpB and GB1, the built-in negative controls

These two landscapes are the reason this campaign matters beyond CreiLOV.
Their assays really are centred near zero (0.35 and 0.20 label SDs), so the
label violation is absent, and **they are exactly the two landscapes every
quotable engine figure in this repository comes from** (HANDOFF §1, "What
survives and can be quoted"). The scale-violation theory therefore predicts
`wt_anchor` is approximately a no-op on both — a negative control the theory
built for itself rather than one chosen for convenience.

| landscape | predicted Δ | predicted t | verdict predicted |
| --- | --- | --- | --- |
| **trpb** | **0.0, 95% interval [−3, +3]** | \|t\| < 2 | not separable |
| **gb1** | **0.0, 95% interval [−3, +3]** | \|t\| < 2 | not separable |

### Falsification, and what it would cost

- **Falsified if \|Δ\| > 5 with \|t\| > 3 on either landscape.** The
  scale-violation theory is then wrong: the row does damage where the label is
  nearly true, so something other than label distance is doing the work —
  most plausibly the 100× variance violation, which is identical everywhere.
- **3 < \|Δ\| ≤ 5, or \|t\| > 3 at \|Δ\| ≤ 3:** equivocal. Reported as
  equivocal, resolved by neither reading.

**If TrpB moves materially, every TrpB figure in this repository changes
meaning** — the +7.93 / t = +17.1 diversity-adjusted result, the 6.3× at a
50-measurement budget, the 45.21 top-1% hits, and the claim that TrpB is the
one landscape untouched by the 2026-10-31 fixes. That is the single most
consequential outcome this campaign could produce, and it is to be reported
**first and loudly**, ahead of the CreiLOV confirmation, if it happens.

TrpB and GB1 are development data under rule 2.4.

---

## 6. Predictions — CR9114-H1, the counter-example

CR9114 violates the label scale by **14.14 SDs**, nearly as badly as CreiLOV's
15.68 — and the engine **wins** there (+3.17 vs ridge, t = +6.10, better on 41
of 48 seeds). A simple "the WT row causes the losses" story does not survive
that, and this campaign must say what it predicts rather than discover a reading
afterwards.

**Primary prediction: Δ positive but far smaller than CreiLOV's, in [0, +10].**

The three outcomes and what each would mean:

**(a) Materially positive (Δ > +3, t > 3) — predicted most likely.** The row
costs hits wherever the scale is violated, including where the engine wins.
The theory is then not "the row causes losses" but the weaker and more
defensible "the row costs hits in proportion to its scale violation; whether
that produces a *net* loss against ridge depends on the landscape's other
properties." CR9114 ceases to be a counter-example and becomes a second
confirmation — the engine won there *despite* the defect.

**(b) Approximately a no-op (|Δ| ≤ 2, |t| < 2).** Then label-scale violation
alone does **not** predict harm, because CR9114 violates it by 14.14 SDs and
is unaffected. CreiLOV's +29.75 would need a CreiLOV-specific explanation —
its 15 mixed-alphabet sites, or its 1,655-variant hit tier — and the headline
finding would have to be restated as a CreiLOV result, not a scale result.
This is the outcome that would most damage the finding while leaving the
CreiLOV number intact.

**(c) Materially negative (Δ < −2, t < −2).** The synthetic row is
*load-bearing* on CR9114. On a binary 16-site bnAb library whose assay is
−log₁₀(K_D) floored at 7.0, a zero-labelled reference may be supplying a
genuine low anchor the measured set never provides. This would mean the
correct fix is **landscape-dependent**, and that a global default change is
unsafe — the strongest possible argument against enacting the fix.

CR9114 is **development data under rule 2.4**: the reachability intervention
was designed against its measured failure. Whatever this arm shows is mechanism
evidence, never generalization evidence.

---

## 7. Secondary, only if time allows

Pre-registered now so that running them later is not a new, unregistered
decision:

- `--ablate-components wt_anchor,prior,physchem_features` on creilov at 48
  seeds. 12-seed point: +41.083. Predicted Δ +41, interval [+34, +48];
  predicted still at or just below ridge (12 seeds gave −9.00).
- eqfp611_red (17 + 4 × 16) and eqfp611_blue. Red's veto ceiling is 0 of 82
  hits, so `wt_anchor` alone is predicted to be **exactly 0.0000 on all 48
  seeds** — arithmetic, not a measurement. Blue, at 1.19 label SDs, is
  predicted a near-no-op, |Δ| ≤ 3.

## 8. What this campaign will not do

It will **not** change the default. Correcting or dropping the row re-bases
every number in the repository, TrpB and GB1 included, and that decision
belongs to whoever owns the freeze. The recommendation and the amendment text
it would require are to be written and left unexecuted.
