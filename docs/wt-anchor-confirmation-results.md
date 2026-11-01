# The fabricated WT training row, at 48 paired seeds

Pre-registered in [`prereg-wt-anchor-confirmation.md`](prereg-wt-anchor-confirmation.md),
committed at `5474ea5` **before any arm ran**, amended at `e3812c9` (CR9114's comparator)
before any CR9114 outcome existed on disk.

Branch `confirm-wt-anchor-48`, based on `bb925c9a3c740dc935728aae3072fc7f85a385d5` and
merged cleanly against `main` at `b5f090b`. 48 paired seeds, `--base-seed 0`, budgets
matched per arm, seed as the unit of independence, 10,000-sample whole-seed bootstrap.
Metric: `cum_top1_hits_excluding_initial` at each seed's final round.

---

## 1. The three answers, in one page

**The negative controls held. TrpB did not move.** Δ = +1.458 at t = +1.60, better on 27
of 48 seeds, bootstrap [−0.292, +3.250] — straddling zero, inside the pre-registered
[−3, +3] band. GB1 likewise: +0.542 at t = +0.41, bootstrap [−2.021, +3.105]. **Every
quotable TrpB and GB1 figure in this repository survives**, with one precision below
(§5): "not separable" is not "identical".

**On CreiLOV the finding is CONFIRMED and was not falsified on any criterion.**
Δ = **+31.771**, t = **+35.88**, better on **48 of 48** seeds, bootstrap
**[+30.000, +33.438]**. The pre-registration predicted +29.8 with interval [+24, +36];
the measurement lands inside it. The dominance threshold was set in advance at +25.04
(half the 50.08-hit deficit), and the bootstrap interval's *lower bound* of +30.000
clears it. The row accounts for **63.2%** of CreiLOV's −50.292-hit deficit against ridge.

**CR9114 is still the counter-example, and it is now a much sharper one.** Δ = **+0.042**,
t = **+0.09**, 18/21/9, bootstrap **[−0.854, +0.938]**. This is pre-registered outcome
**(b)**, which §6 named as "the outcome that would most damage the finding while leaving
the CreiLOV number intact". **My primary CR9114 prediction — outcome (a), a materially
positive effect — was wrong**, and the pre-registration said in advance what being wrong
there would mean.

So the headline must be stated carefully, because the obvious version of it is false:

> The fabricated row is the dominant cause of the **CreiLOV** loss. It is **not** a
> general explanation of the engine's losses, and **label-scale violation does not
> predict its effect.**

---

## 2. The primary table

Δ = `panoptes_default --ablate-components wt_anchor` − `panoptes_default`, paired by seed.
`MDE@t2` is the smallest |Δ| this 48-seed design would have called significant, printed
beside every null so a null cannot be mistaken for a shortage of power.

| landscape | label SDs | ridge | default | ablated | **Δ** | **t** | **B/W/T** | **boot95** | MDE@t2 | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **creilov** | 15.68 | 57.229 | 6.938 | 38.708 | **+31.771** | **+35.88** | **48/0/0** | **[+30.000, +33.438]** | 1.771 | **CONFIRMED** |
| **trpb** | 0.35 | 31.104 | 49.646 | 51.104 | +1.458 | +1.60 | 27/18/3 | [−0.292, +3.250] | 1.826 | not separable |
| **gb1** | 0.20 | 33.292 | 20.396 | 20.938 | +0.542 | +0.41 | 26/20/2 | [−2.021, +3.105] | 2.673 | not separable |
| **cr9114_h1** | 14.14 | 1.833 | 5.000 | 5.042 | **+0.042** | **+0.09** | 18/21/9 | **[−0.854, +0.938]** | 0.924 | **null** |
| eqfp611_blue | 1.19 | 8.792 | 9.521 | 11.146 | +1.625 | +7.83 | 38/1/9 | [+1.229, +2.021] | 0.415 | small, real |
| eqfp611_red | 0.80 | — | — | — | — | — | — | — | — | **REFUSED** |

**eqFP611-red is refused, not reported.** Both arms ran 7 and 9 seeds short of the
81-measurement budget with status `short_batch`, so it is not budget-matched and not
quotable. This independently reproduces HANDOFF §1's "7 short batches" note. The refusal
is the analyser's, not a judgement call: it checks every seed's final
`n_unique_measurements` and `status` before computing anything.

### Secondary, also pre-registered (§7)

| arm, creilov, 48 seeds | engine | Δ vs default | t | B/W/T | boot95 | vs ridge |
| --- | --- | --- | --- | --- | --- | --- |
| `wt_anchor` | 38.708 | +31.771 | +35.88 | 48/0/0 | [+30.000, +33.438] | −18.521 (t −10.87) |
| `wt_anchor,prior,physchem_features` | 49.521 | **+42.583** | **+46.24** | 48/0/0 | [+40.771, +44.333] | **−7.708** (t −5.59) |

Predicted +41 with interval [+34, +48]; measured +42.583. Inside.

### The engine still loses CreiLOV, as pre-registered

§4 predicted the ablated arm would reach ≈ −20.3 against ridge and that "the engine is
predicted to still lose CreiLOV after this fix". Measured: **−18.521**, t = −10.87,
better on 3 of 48 seeds. Even with all three components ablated it is −7.708 on 10 of 48.
**Nothing here is a CreiLOV win**, and any later summary saying "the CreiLOV loss is
fixed" contradicts both the prediction and the measurement.

---

## 3. Why the scale theory is dead, and what killed it

The theory was that the row does damage in proportion to how far `y = 0.0` sits from the
assay's real labels. That ordering is refuted by its own two best cases pointing opposite
ways. Relative effect is Δ as a percentage of the default arm's own score, which removes
the objection that CR9114's endpoint is simply too compressed to show a large absolute
gain:

| landscape | label SDs | 0.0 below min, in label ranges | Δ absolute | Δ relative |
| --- | --- | --- | --- | --- |
| creilov | 15.68 | 1.99× | +31.771 | **+457.9%** |
| cr9114_h1 | **14.14** | **2.47×** | **+0.042** | **+0.8%** |
| eqfp611_blue | 1.19 | — | +1.625 | +17.1% |
| trpb | 0.35 | 0.00× (0.0 is inside the range) | +1.458 | +2.9% |
| gb1 | 0.20 | 0.00× (0.0 is the minimum) | +0.542 | +2.7% |

- The label-SD covariate predicts **creilov ≈ cr9114 ≫ blue > trpb > gb1**.
- The measured ordering is **creilov ≫ blue > trpb ≈ gb1 > cr9114**.
- CR9114 is at the **wrong end** of the measured ordering.

Range-normalising does not rescue it: CR9114's 0.0 sits **2.47** label-ranges below its
minimum against CreiLOV's **1.99**, so the *more* violated landscape on that measure is
the one with no effect. And this is not a power failure — CR9114's MDE is **0.924 hits**,
i.e. 18% of the engine's own 5.000, and the measurement is +0.042 with a bootstrap
interval of ±0.9.

**Both candidate scale covariates order CR9114 at or above CreiLOV, and CR9114 shows
nothing.** Whatever makes CreiLOV sensitive, it is not distance from zero.

This also means the two mechanisms §1 of the pre-registration deliberately refused to
conflate can now be partly separated: the **variance** violation (`1e-4` against `1e-2`)
is identical on all nine landscapes, and four of the six landscapes here show effects at
or under 3% relative. So the variance violation alone is small, and the large CreiLOV
effect needs a CreiLOV-specific interaction. **What that interaction is remains
unidentified, and this campaign does not claim it.** A plausible next covariate — offered
as a hypothesis, not a result — is that CreiLOV is the only landscape here whose labels
are both uncensored at the floor *and* tightly packed (sd 0.240 within a 1.403 range),
so the surrogate has no genuine low-fitness evidence to compete with the fabricated row.
CR9114 floors at exactly 7.000 and GB1 at exactly 0.000.

---

## 4. What the campaign proves about its own instrument

Reported because an unaudited instrument makes every number above a guess.

**Three of HANDOFF §1's published engine-vs-ridge figures reproduced to three decimals**,
from run directories produced today on a different branch:

| contrast | HANDOFF §1 | measured here |
| --- | --- | --- |
| creilov default − ridge | −50.29, t −32.01, 0/48/0 | **−50.292, t −32.01, 0/48/0** |
| gb1 default − ridge | −12.90, t −5.63, 6/41/1 | **−12.896, t −5.63, 6/41/1** |
| cr9114_h1 default − ridge | +3.17, t +6.10, 41/5/2 | **+3.167, t +6.10, 41/5/2** |
| eqfp611_blue default − ridge | +0.73, t +0.77, 25/18/5 | **+0.729, t +0.77, 25/18/5** |

CreiLOV's −50.292 reproduced **three** times from three independent run directories.

**Pairing proved, not assumed.** Within every landscape, `ridge`'s full per-round
sub-table — all 240 rows, every column except wall-clock timings — is **bit-identical**
across the two arms, `configuration.seeding.seed_block_sha256` matches, and
`results.top1.threshold_value` matches. The analyser refuses a landscape rather than
warning. It also verifies that `configuration.components.states.wt_anchor` flipped
`True → False` between the arms, so the manifest is confirmed to describe the arm that
actually ran.

**Instrument audit, completed before the pre-registration was committed** (§3 there):
the ablation removes the synthetic row *and* its `1e-4` variance entry, which matters
because `math_worker.py:2170` silently sets `sample_variances_arr = None` on a length
mismatch rather than raising; the test pinning that is **non-vacuous**, verified by
mutation (moving the variance append outside the guard produces
`AssertionError: None != 2`, then reverted); `wt_anchor` has **no module-level global**,
so under rule 2.12 the `--ablate-components` set is the authoritative gate and the
manifest reads it correctly; and the flag-absent path is **metric-identical** to `main`'s
engine on 4 CreiLOV seeds (Δ = 0.000, tied 4/4). Ablating nothing changes nothing.

Rules 2.1 and 2.2 do not bind this analysis: nothing here parses a mutation string. All
figures come from the harness's own `round_metrics.csv` and `manifest.json`.

---

## 5. An arithmetic error in the published characterisation

**The direction of the finding is unaffected. The magnitude as written is wrong.**

HANDOFF §1, `docs/demo-script.md`, the comment at `workers/math_worker.py:2087`, the
docstring of `tests/test_wt_anchor_ablation.py`, and (inherited) §1 of this campaign's own
pre-registration all say the row

> "sits **11.6 units below the minimum — 8.3× the entire label range**, 15.68 label SDs
> from centre".

Measured directly from `CREILOV_CHLRE_Chen_2023.csv` (n = 165,427; labels 2.793–4.196,
mean 3.769, sd 0.240, range 1.403):

| as published | actual |
| --- | --- |
| 11.6 **units** below the minimum | **2.793** raw units below the minimum |
| — | which is **11.62 label SDs** below the minimum |
| 8.3× the entire label range | **1.99×** the label range |
| 15.68 label SDs from centre | **15.68** — correct |

"11.6" is the distance to the minimum **in label SDs**, mislabelled as raw units; "8.3×"
then divides that SD figure by a **raw** range (11.62 / 1.403 = 8.30), mixing units. The
range-based overstatement is ≈ 4.2×. The 15.68-SD figure, which is the one the ablation
argument actually rests on, is correct.

**Not corrected here**, because four of the five sites are inside the freeze's published
prose and one is a test docstring; the correction is listed with the amendment text in §7.

---

## 6. Gate verdicts, including the failures

Run **with `--run-dir`**, which is the point: `generality_sweep.ps1` omits it, and that is
why nobody had seen `threshold.<lib>` fail. Confirmed here — the same gate without
`--run-dir` reports 4 passed / 4 NOT-APPLICABLE and **exit 0**; with it, 6 passed /
1 FAILED and **exit 1**.

All **14** post-run gate invocations (one per run, across six landscapes in three
campaigns) returned the same pattern, with eqFP611-red the only exception (below):

| check | verdict |
| --- | --- |
| `corpus.lock` | **PASS** — every pinned file present, both hashes match |
| `corpus.derived_table` | PASS |
| `corpus.prepared` | PASS |
| `budget.panoptes_default` | **PASS** — all 48 seeds spent all 125 distinct measurements |
| `budget.ridge` | **PASS** — all 48 seeds spent all 125 |
| `endpoint.<lib>` | PASS |
| **`threshold.<lib>`** | **FAIL** |
| `encoding.<arm>` | NOT-APPLICABLE — neither `ridge` nor `panoptes_default` normalises an encoded design space by column range (only `alde` does) |
| `row_order.<lib>` | NOT-APPLICABLE — needs an `--invariance-runner`; not supplied |

**eqFP611-red additionally FAILED `budget.panoptes_default`** — "71 of 81" and "70 of 81"
distinct measurements on the two arms. The gate caught the short budget independently of
the analyser's refusal in §2, which is the two instruments agreeing: red is not
budget-matched and is not reportable. Its `threshold` check failed as well, so red is
5 passed / 2 FAILED.

**`threshold.<lib>` FAILED on all 14 runs**, for the reason HANDOFF §1 already records:
`qualification.quantile` and `qualification.source` are not written to the manifest, so
the bar is an unattributed absolute. This is **not** reported as a pass anywhere in this
document. It cannot bias any paired contrast here — the bar is consumed only by the
`qualify_then_optimize` / `conditional_probe` arms, neither of which ran, and the failure
is symmetric across arms within each landscape (identical `threshold_value`, verified) —
but it is an un-passed check that these figures carry, and it is stated wherever they are.

**`corpus.lock` requires `--corpus-dir <real root>\downloads\combingym` from a worktree.**
`scripts/freeze_cross_protein_corpus.py:55` resolves the corpus against its own
`__file__` and `downloads/` is gitignored, so without the flag all 29 pinned files report
`absent` and the check FAILS as an **environment artifact**, not as rule-2.5 drift.
Verified directly: `CORPUS VERIFIED: 29/29 files match CombinGym @ c5dccd514a9a`. Note
the flag takes the **corpus directory**, not the repository root — pointing it at the real
root still fails.

---

## 7. Recommendation — written, deliberately NOT executed

### The recommendation

**Do not change the default yet, and do not change it as a silent default flip.** The
evidence supports a change; the evidence does **not** yet support it being safe to
re-base the repository on.

What the change should be, when someone with authority over the freeze takes it:

1. **Use the measured wild-type label where the library provides one; omit the row
   entirely where it does not.** A synthetic `y = 0.0` is a fabricated measurement and
   has no defence on any landscape. Note that omitting it is never measurably harmful
   here: the worst point estimate across six landscapes is **+0.042** and no interval
   reaches materially below zero.
2. **Give the row the same variance a real measurement gets** (`1e-2`, not `1e-4`). A
   true label asserted 100× over-confidently is still a defect, and it is the one part of
   this defect that is identical on all nine landscapes.
3. **Ship it as a named policy version, not an edit to the current one.** HANDOFF §5.1 is
   explicit: keep the old ID alongside the new and never silently modify the old one and
   reuse its results. `wt_anchor` should become e.g. `wt_anchor_measured_v2` with the
   legacy behaviour retained as `wt_anchor_synthetic_v1`.
4. **Re-measure all nine landscapes at 48 paired seeds before adopting it.** "Not
   separable" is not "identical": TrpB's point estimate is **+1.458** against an MDE of
   **1.826**, so a real effect of one to three hits there is not excluded. Every quotable
   TrpB figure would need re-stating even though the control held.
5. **Fix the gate and the arithmetic first**, since both are cheap and both currently
   travel with the figures: record `qualification.quantile` and `qualification.source` in
   the manifest so `threshold.<lib>` can pass, and correct the five sites carrying the
   unit-mixing error in §5.

### Rule 2.4 has fired, and this is the most important consequence

**CreiLOV is now development data.** The `wt_anchor` component exists *because* of
CreiLOV's measured failure, and the algorithm was changed because of what was seen there.
Rule 2.4 is unambiguous: "If you inspect a held-out target and change the algorithm
because of what you saw, that target is reclassified as development data and a new one
must be reserved." So **+31.771 is mechanism evidence, not generalization evidence**, and
a new independent landscape must be reserved before this is presented as a general fix.

The only generalization evidence this campaign produced is **eqFP611-blue: +1.625 at
t = +7.83 on 38 of 48 seeds** — small, real, and on a landscape that played no part in
designing the fix.

### What is now unexplained, and was not before

HANDOFF §1 said "The CreiLOV and GB1 losses are therefore unexplained" after the veto
hypothesis fell. CreiLOV is now 63.2% explained. **GB1's −12.896 remains unexplained**:
the row is worth +0.542 there, i.e. 4.2% of its deficit. GB1's loss has now survived two
pre-registered root-cause campaigns (the spatial veto, which *helped*, and this one).

### Amendment text, if the change is adopted

For HANDOFF §1, replacing the "Next:" line of the fabricated-row blockquote:

> **Confirmed at 48 paired seeds**, branch `confirm-wt-anchor-48`, pre-registered in
> `docs/prereg-wt-anchor-confirmation.md` before any arm ran. On CreiLOV
> `--ablate-components wt_anchor` is worth **+31.771** top-1% hits (t = +35.88, 48 of 48
> seeds, bootstrap [+30.000, +33.438]), which is **63.2%** of the −50.292 deficit — the
> pre-registered dominance threshold of +25.04 is cleared by the interval's lower bound.
> The engine **still loses** CreiLOV at −18.521, and −7.708 with
> `wt_anchor,prior,physchem_features`.
>
> **The negative controls held.** TrpB +1.458 (t = +1.60, 27/18/3) and GB1 +0.542
> (t = +0.41, 26/20/2) are not separable, so every quotable TrpB and GB1 figure stands —
> but neither is *identical*, and TrpB's point estimate sits just under its 1.826 minimum
> detectable effect.
>
> **The scale-violation ordering is FALSIFIED.** CR9114 violates the label scale by 14.14
> SDs — and by 2.47 label-ranges, *more* than CreiLOV's 1.99 — and
> `--ablate-components wt_anchor` is worth **+0.042** there (t = +0.09, 18/21/9, bootstrap
> [−0.854, +0.938]) against a minimum detectable effect of 0.924. The measured ordering is
> creilov ≫ eqfp611_blue > trpb ≈ gb1 > cr9114, which puts CR9114 at the wrong end. The
> row is the dominant cause of the **CreiLOV** loss specifically; it is **not** a general
> explanation of the engine's losses, and distance from zero does not predict its effect.
> **GB1's −12.896 remains unexplained.**
>
> **Rule 2.4 has fired: CreiLOV is now development data**, because `wt_anchor` was built
> against its measured failure. The +31.771 is mechanism evidence. The only
> generalization evidence is eqFP611-blue, +1.625 at t = +7.83 on 38 of 48 seeds.
> eqFP611-red is **refused** — both arms ran short of budget on 7 and 9 seeds.

And a correction to be made wherever the characterisation appears (HANDOFF §1,
`docs/demo-script.md`, `workers/math_worker.py:2087`,
`tests/test_wt_anchor_ablation.py`):

> ~~sits 11.6 units below the minimum — 8.3× the entire label range~~ → **sits 2.793
> units below the minimum, which is 11.62 label SDs or 1.99× the entire label range**,
> 15.68 label SDs from centre.

---

## 8. Reproducing this

```powershell
# From the worktree. -Command, not -File: -File binds a bare array positionally and
# `trpb` lands in -Seeds.
powershell -NoProfile -ExecutionPolicy Bypass -Command `
  "& 'scripts\analysis\wt_anchor_confirmation.ps1' -Only creilov,trpb,gb1,cr9114_h1 -Seeds 48"
```

Run directories, resolved by name from each campaign's own log (never "the newest
directory"):

| landscape | default | ablated |
| --- | --- | --- |
| creilov | `run_20261101_153459_0ffe66ff` | `run_20261101_153933_1371483b` |
| trpb | `run_20261101_154322_214bb5c6` | `run_20261101_154533_0e42b04f` |
| gb1 | `run_20261101_154808_51d7d012` | `run_20261101_154927_7df31826` |
| cr9114_h1 | `run_20261101_155046_465ea4b2` | `run_20261101_155309_7c564936` |
| creilov (triple) | `run_20261101_155738_e84f4ef2` | `run_20261101_160322_549dce7c` |
| eqfp611_blue | `run_20261101_161002_efadfb52` | `run_20261101_161139_d3f2ffcb` |
| eqfp611_red | `run_20261101_160653_ed6a76f6` | `run_20261101_160827_5a760dd0` (both REFUSED) |

They live in the **gitignored** `benchmark_comparison_results/` of this worktree. The
results root is shared across worktrees, so every file this campaign wrote is stamped
`wtanchor*`, and `run_20261101_144437_27fc1a1f` — an orphan with one file and no manifest
left by another session — is excluded by construction, since directories are resolved by
name from this campaign's log and the analyser refuses any run without a `manifest.json`.
