# The wild-type anchor, repaired rather than deleted — measured

Branch `wt-anchor-measured-label`. Pre-registered in
[`docs/prereg-wt-anchor-fix.md`](prereg-wt-anchor-fix.md), committed at `a032170`
**before any arm ran**. 48 paired seeds per landscape, `--base-seed 0`, budget 125
(`--initial-size 25 --batch-size 25 --rounds 4`), three run folders per landscape sharing
one seed block.

Reproduce the availability table with
`.venv\Scripts\python.exe -X utf8 scripts/analysis/wt_label_availability.py`.

---

## 0. The verdicts, first

1. **Anchor at the measured label. Do not remove the anchor.** The repair is **equal or
   better than removal on every landscape measured, strictly better on one, and never
   worse anywhere** — on either endpoint. On GB1 removal is worth nothing
   (+0.54 raw, p=0.69) while the repair is worth **+7.04 raw (t=+4.31) and +0.90 at d≥2
   (t=+2.89)**, so the two are not interchangeable and ablation is not the fix.
2. **But the answer is not the same on all four landscapes, and it cannot be known for
   five of the nine.** Four landscapes gave four qualitatively different answers. On three
   of the five unmeasured ones (CR6261 ×2, HER2) the parent's label is a **censored assay
   floor**, so the repair as implemented would anchor at a known-wrong value there and
   must not be used until that is handled.
3. **The label-scale theory is FALSIFIED, and so is my own reframe of it.** Neither
   distance-from-centre nor anchor-error orders the effects. CR9114 is second-largest on
   both metrics and is a tight null; GB1 is second-smallest on both and carries the only
   repair-specific win.
4. **The two contrasts are independent mechanisms and never co-occur.**
   `ablated − zero` fires on CreiLOV alone; `measured − ablated` fires on GB1 alone.
   Removing a poisonous row and adding a true observation are different things, and no
   landscape measured shows both.
5. **TrpB moved, and the movement is redundancy.** `measured − zero` is +1.94 raw at
   p=0.0497 — separable — but **+0.42 at d≥2 at p=0.347, not separable.** The raw gain is
   the same variants found again. Reported first and loudly below.
6. **The flag-absent default is byte-identical**: 18,240 cells, 19 behavioural columns,
   4 arms × 48 seeds × 5 rounds, **zero differences**.
7. **The default is NOT changed.** §7 holds the amendment and it is unexecuted.

---

## 1. TrpB FIRST, because it is the designated negative control and it MOVED

| contrast | raw | d≥2 |
| --- | --- | --- |
| `measured − zero` | **+1.938, t=+2.01, p=0.0497, 31/16/1, CI95 [+0.062, +3.812] — SEPARABLE** | +0.417, t=+0.95, p=0.347, 26/18/4, CI95 [−0.438, +1.250] — not separable |
| `measured − ablated` | +0.479, t=+0.48, p=0.631, 22/22/4, CI95 [−1.376, +2.438] | +0.292, t=+0.60, p=0.555, 26/18/4, CI95 [−0.667, +1.229] |
| `ablated − zero` | +1.458, t=+1.60, p=0.117, 27/18/3, CI95 [−0.292, +3.208] | +0.125, t=+0.31, p=0.761, 20/19/9, CI95 [−0.688, +0.917] |

**What this means, stated carefully.** The raw contrast is separable, with a 95% CI whose
lower bound is **+0.062** — it clears zero by six hundredths of a hit. The
diversity-adjusted contrast is not separable and its CI straddles zero. So:

* **the repair does change TrpB**, which contradicts the brief's premise that TrpB is a
  negative control for this work, and which
  [`docs/prereg-wt-anchor-fix.md`](prereg-wt-anchor-fix.md) §5a predicted in advance;
* **but the change is entirely redundant hits.** Raw and d≥2 disagree, and under the
  rule the coordinator set — a raw-only contrast cannot tell "found more good variants"
  from "found the same variants again" — the honest reading is **no demonstrated gain in
  distinct solutions on TrpB.**
* The engine's TrpB advantage over ridge is *unharmed* in all three modes:
  +18.54 / +20.00 / +20.48 raw and +6.46 / +6.58 / +6.88 at d≥2 for zero / ablated /
  measured. **Nothing here threatens the repository's one surviving win** — but every
  TrpB figure would still have to be re-measured if the default changed, because the arm
  mean moves 49.65 → 51.58 raw.

**On whether TrpB is a null.** Pre-registration F3 said the `err-SD` reframe is falsified
if TrpB moves on *neither* contrast (|Δ| < 1.0 and p > 0.1 on both). It moved on one
(+1.94, p=0.0497), so **F3 did not fire** — but the margin is one endpoint and 0.0497,
and the d≥2 result goes the other way. The reframe is **not** vindicated; see §5.

---

## 2. Step 1 — is the wild type's measured label available, per landscape?

Established **before a line of the fix was written**. Three independent sources, each
cross-checked against a fourth: the bracket the prepared table's own `DMS_score_bin`
column implies, since `docs/dataset-inventory.md` §7 records the rule
`DMS_score_bin = 1 iff DMS_score >= the wild-type value`.

| landscape | in the prepared table? | WT measured label | source | censored? | 0-SD | err-SD | WT pct | WT a hit? |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| creilov | **no** | 4.023229753 | upstream `CreiLOV_clean.xlsx` `wt` row, `mean_log` | no | 15.68 | 16.74 | 87.99 | no |
| cr9114_h1 | **no** | 8.424622015 | in-repo `prepare_cr9114_h1.py::CR9114_GERMLINE_SCORE` | no | 14.14 | 13.18 | 19.18 | no |
| cr6261_h1 | **no** | 7.0 | upstream `bnAbs_CR6261_H1_clean.xlsx` `n_mut == 0` | **YES, assay floor** | 11.00 | 8.75 | 0.00 | no |
| cr6261_h9 | **no** | 7.0 | upstream `bnAbs_CR6261_H9_clean.xlsx` `n_mut == 0` | **YES, assay floor** | 8.35 | 6.72 | 0.00 | no |
| eqfp611_blue | **no** | 1.574628255 | upstream `eqFP611_clean.xlsx` `n_mut == 0`, `blue` | no | 1.19 | 4.25 | 99.84 | **HIT** |
| eqfp611_red | **no** | 0.085316563 | upstream `eqFP611_clean.xlsx` `n_mut == 0`, `red` | no | 0.80 | 0.29 | 23.62 | no |
| her2 | **YES** (literal `WT`) | 0.1 | the prepared CSV's own `WT` row | **YES, detection floor** | 0.56 | 0.19 | 0.00 | no |
| trpb | **no** | 0.408073925 | in-repo `workers/data/TRPB_TM9D8.json` provenance block | no | 0.35 | 6.81 | 99.32 | **HIT** |
| gb1 | **YES** (literal `WT`) | 1.0 | the prepared CSV's own `WT` row (parent ≡ 1.0 by assay definition) | no | 0.20 | 2.53 | 97.56 | no |

* `0-SD` = `|0 − label mean| / label sd` — how much of an **outlier** the fabricated row
  is. This is the ordering `docs/HANDOFF.md` records.
* `err-SD` = `|wt_label − 0| / label sd` — how **wrong** the fabricated label is about the
  point it claims to pin.
* Units are named on every figure, because mixing them is what produced the withdrawn
  "8.3× the entire label range" (see §8).

**Findings.**

1. **A measured WT label exists on all nine landscapes — but on only two is it in the
   prepared table.** The repo convention drops the parent ("the parent is not a library
   member"), so seven need the pinned upstream workbook or an in-repo constant. **Five of
   the nine therefore depend on the gitignored `downloads/combingym/`**, which is a
   provenance obligation, not a detail (§7).
2. **On three of the nine the value is a CENSORED FLOOR, not a measurement.** CR6261-H1
   and H9 sit at the 7.0 `−log₁₀K_D` floor with all three replicates identical; HER2's
   parent sits at the 0.1 detection floor that 5,833 of its 11,277 labelled rows share,
   at the **0.00th percentile** — not credible for the parent of an affinity-matured
   therapeutic antibody. **Anchoring at a censored bound is anchoring at a known-wrong
   value**, so the repair as implemented must not be applied to those three until
   censoring is handled. This is where the fix needs more than a label.
3. **Every value is independently bracketed.** The GB1 upstream join was *refused* at
   5e-7 (our prepared CSV stores six decimal places) and the value still resolves from
   the prepared table's own row — the refusal is the check working, not a failure.

### 2a. A documentation inaccuracy found on the way

`docs/dataset-inventory.md` §7 names the TrpB `DMS_score_bin` cutoff as **0.408085**,
in a list presented as the per-landscape wild-type values. 0.408085 lies **outside** the
bracket the column itself implies, `(0.40785355, 0.408084865]`, by 1.4e-7: it is
`min{score : bin == 1}` rounded — the smallest *mutant* at or above the cutoff
(`V183E:F184A:S228N`) — not the parent. The parent is **0.408073925**, which does lie
inside. Not fixed here; the inventory is another agent's artifact.

---

## 3. The three behaviours, and what actually ran

| arm | how | anchor row | label | variance |
| --- | --- | --- | --- | --- |
| **zero** (today's default) | no flag | present | literal `0.0` | `1e-4` |
| **ablated** | `--ablate-components wt_anchor` | absent | — | — |
| **measured** (the repair) | `--wt-anchor-mode measured` | present | the library's recorded parent label | `1e-2` |

**The variance is `1e-2` and that is a parity statement, not a tuning knob.** `1e-4`
claims 100× the precision of a real unreplicated reading, which this module itself puts at
`1e-2`. Replicate-derived variances where published are all **below** `1e-2` — CreiLOV
1.76e-3 over three log-scale replicates, CR9114 7.67e-3 from its reported SEM — so `1e-2`
cannot over-trust the anchor, and an override is clamped at `1e-4` so the over-confidence
cannot return by accident. The measured anchor also sits in **residual space** with the
prior **evaluated** at the wild type; the `zero` path hardcodes that prior to `0.0`, which
is only coincidentally right.

### 3a. Rule 2.12 — what ran, verified against the code path, and not vacuous

Two manifests in this project wrongly claimed `dpp` and `dkl` were enabled because the
recorder read the CLI flag rather than the module-level global. So here:

* `select_next_micro_batch` returns what it actually did on
  `SelectedBatch.wt_anchor_state`, **written by the branch that executed**;
* the harness compares that to the arm's configuration **every round** and raises
  `BenchmarkError` on a mismatch, so a misdescribing run cannot finish;
* `run.wt_anchor_verified` in the manifest is built from that return value;
* it declares itself **`vacuous`** when no arm in the run called the engine, because a
  check that cannot fail must say so (rule 2.10, where a threshold check passed the exact
  defect it was written for).

All 12 campaign folders: `agrees_with_configuration = true`, `vacuous = false`,
`rounds_observed` = 192 each (48 seeds × 4 rounds). **The check is demonstrably not
vacuous**: a unit test drives it with an engine report that disagrees with the
configuration and asserts it reports the disagreement, and another drives it with a
baseline-only run and asserts it reports `vacuous`.

The three modes also produce **three distinct configuration fingerprints** — `0c991820…`,
`aae05372…`, `debc7956…` — so a zero-anchor and a measured-anchor run can no longer hash
identically, which was the concrete harm rule 2.12 records.

### 3b. Variance/row alignment, which used to fail silently

The four per-row training lists are built across three append sites that each append a
feature vector **before** its label and variance, inside a broad `except Exception`. A
raise between the two left `sample_variances_arr` as `None` — which
`DeepKernelLearner.fit` *also* drops silently — so the round quietly became homoscedastic
while `y_train`, bent by the same event, stayed misaligned, and nothing in the artifacts
said so. It now raises and voids the round, and a test forces the desynchronisation.

---

## 4. Results, 48 paired seeds, both endpoints

`d≥2` is the pre-registered diversity adjustment: hits fitness-descending, keep one only
when it is Hamming ≥ 2 from everything already kept, computed with
`veritas.calibration.endpoint.greedy_diverse` on genotypes canonicalised through
`veritas.calibration.genotypes.canonical_combo` (rule 2.1 — wild type filled at every
omitted position, nothing zipped, no sixth implementation). The raw endpoint is taken from
`round_metrics.csv` **and** independently recomputed from `selection_trace.csv`: they
agree on **192 of 192** (landscape × mode × seed) cells.

### Validity, by construction rather than by judgement

All 12 folders: `manifest.json` present, `run.state == SUCCESS`, **48 seeds completed,
every seed at exactly 125 unique measurements, for all three arms**, anchor state verified.
Nothing was discarded, because nothing failed. Pairing premises, checked not assumed:
`seed_block_sha256` identical across each landscape's three folders; `initial_sets_sha256`
identical for all 48 seeds; the top-1% threshold **read from each manifest** (rule 2.3) and
identical; and `ridge` and `random` **numerically identical cell-for-cell across all three
folders** of every landscape — those arms never read the anchor, so this is the direct test
that the arms are paired.

### 4a. CreiLOV — threshold 4.122311456637142, ridge 57.229 raw / 23.208 d≥2

| mode | raw mean (sd) | d≥2 mean (sd) |
| --- | --- | --- |
| zero | 6.938 (5.357) | 5.062 (3.744) |
| ablated | 38.708 (4.172) | 20.688 (3.440) |
| measured | 38.646 (4.275) | 21.104 (3.075) |

| contrast | raw | d≥2 |
| --- | --- | --- |
| `measured − zero` | **+31.708, t=+32.20, p=1.1e−33, 48/0/0, CI95 [+29.792, +33.583]** | **+16.042, t=+24.29, p=3.0e−28, 48/0/0, CI95 [+14.729, +17.292]** |
| `measured − ablated` | −0.062, t=−0.10, p=0.922, 20/23/5, CI95 [−1.312, +1.167], MDE 1.785 | +0.417, t=+0.82, p=0.414, 28/17/3, CI95 [−0.542, +1.354], MDE 1.414 |
| `ablated − zero` | **+31.771, t=+35.88, p=8.3e−36, 48/0/0, CI95 [+30.062, +33.500]** | **+15.625, t=+24.04, p=4.7e−28, 48/0/0** |
| engine − ridge, zero | −50.292, t=−32.01 | −18.146, t=−24.58 |
| engine − ridge, ablated | −18.521, t=−10.87 | −2.521, t=−3.31 |
| engine − ridge, measured | −18.583, t=−11.36 | −2.104, t=−2.91 |

**`ablated − zero` reproduces the independent `confirm-wt-anchor-48` figure to three
decimals: +31.771 against +31.771.** Different branch, different driver, same seed block.

On CreiLOV the repair and the removal are **indistinguishable** on both endpoints, and the
engine still loses to ridge in both (−18.6 raw, −2.1 at d≥2). The defect accounts for
**63.2%** of the −50.29 deficit on raw and **86.1%** of the −18.15 deficit at d≥2; the
remainder is not the anchor.

### 4b. TrpB — threshold 0.3493993794399995, ridge 31.104 raw / 7.229 d≥2

| mode | raw mean (sd) | d≥2 mean (sd) |
| --- | --- | --- |
| zero | 49.646 (5.674) | 13.688 (1.776) |
| ablated | 51.104 (6.846) | 13.812 (2.598) |
| measured | 51.583 (6.424) | 14.104 (2.684) |

Contrasts in §1. **Raw and d≥2 disagree**: the only separable contrast is raw
`measured − zero` at p=0.0497, and its d≥2 counterpart is p=0.347.

### 4c. GB1 — threshold 2.152377110000002, ridge 33.292 raw / 6.896 d≥2

| mode | raw mean (sd) | d≥2 mean (sd) |
| --- | --- | --- |
| zero | 20.396 (9.022) | 6.938 (2.108) |
| ablated | 20.938 (9.416) | 7.104 (2.195) |
| measured | 27.438 (9.088) | 7.833 (1.993) |

| contrast | raw | d≥2 |
| --- | --- | --- |
| `measured − zero` | **+7.042, t=+4.31, p=8.3e−05, 34/12/2, CI95 [+3.938, +10.313]** | **+0.896, t=+2.89, p=0.0058, 27/12/9, CI95 [+0.292, +1.500]** |
| **`measured − ablated`** | **+6.500, t=+3.70, p=0.00056, 33/15/0, CI95 [+3.125, +9.958]** | **+0.729, t=+1.99, p=0.0526, 23/13/12, CI95 [+0.021, +1.458]** |
| `ablated − zero` | +0.542, t=+0.41, p=0.687, 26/20/2, CI95 [−2.042, +3.083], MDE 3.742 | +0.167, t=+0.53, p=0.597, 23/20/5, CI95 [−0.438, +0.771], MDE 0.876 |
| engine − ridge, zero | −12.896, t=−5.63 | +0.042, t=+0.09, **not separable** |
| engine − ridge, ablated | −12.354, t=−5.02 | +0.208, t=+0.48, not separable |
| engine − ridge, measured | −5.854, t=−2.45 | **+0.938, t=+2.12, p=0.0392, SEPARABLE** |

**This is the landscape that decides what ships.** Removing the anchor is worth nothing
(+0.54 raw against an MDE of 3.74; +0.17 at d≥2 against an MDE of 0.88). Repairing it is
worth **+7.04 raw and +0.90 at d≥2, separable on both**, and the two endpoints **agree**,
so it is not redundancy. The gain is therefore **information in the true label, not
removal of a false one** — and ablation could never have found it.

Note also that `measured − ablated` is separable at d≥2 by bootstrap CI
([+0.021, +1.458]) while its t-test p is 0.0526, just over the conventional line. Both are
reported; the verdict rests on the whole-seed bootstrap, which is the protocol's test.

### 4d. CR9114-H1 — threshold 9.608730734985462, ridge 1.833 raw / 1.479 d≥2

| mode | raw mean (sd) | d≥2 mean (sd) |
| --- | --- | --- |
| zero | 5.000 (3.101) | 4.354 (2.401) |
| ablated | 5.042 (2.917) | 4.292 (2.153) |
| measured | 5.083 (3.079) | 4.312 (2.299) |

| contrast | raw | d≥2 |
| --- | --- | --- |
| `measured − zero` | +0.083, t=+0.18, p=0.857, 24/18/6, CI95 [−0.812, +0.958], MDE 1.291 | −0.042, t=−0.12, p=0.905, CI95 [−0.729, +0.625], MDE 0.972 |
| `measured − ablated` | +0.042, t=+0.16, p=0.870, 19/16/13, CI95 [−0.458, +0.542], MDE 0.707 | +0.021, t=+0.10, p=0.919, CI95 [−0.375, +0.417], MDE 0.568 |
| `ablated − zero` | +0.042, t=+0.09, p=0.929, CI95 [−0.833, +0.958], MDE 1.294 | −0.062, t=−0.18, p=0.862, CI95 [−0.750, +0.646] |
| engine − ridge, all three modes | +3.167 / +3.208 / +3.250 | +2.875 / +2.812 / +2.833 |

**A tight null on all three contrasts and both endpoints**, reproducing
`confirm-wt-anchor-48`'s `ablated − zero` of +0.042 exactly. CR9114 violates the label
scale by 14.14 SDs — second only to CreiLOV — and the anchor is worth nothing there at all.
The engine's win over ridge is unchanged in every mode.

### 4e. The pattern, and what it is not

| landscape | err-SD rank (of 9) | `ablated − zero` | `measured − ablated` |
| --- | --- | --- | --- |
| creilov | 1 (16.74) | **+31.77 raw / +15.63 d≥2** | null |
| cr9114_h1 | 2 (13.18) | null | null |
| trpb | 5 (6.81) | null (p=0.117) | null |
| gb1 | 7 (2.53) | null | **+6.50 raw / +0.73 d≥2** |

The two effects **never co-occur**, and the scalar that was supposed to order them
inverts: the landscape ranked 2nd is a null and the landscape ranked 7th carries the only
repair-specific win. See §5.

---

## 5. The pre-registration: what held, what was falsified

| prediction | predicted | measured (raw) | verdict |
| --- | --- | --- | --- |
| §5a creilov `measured − zero` | +30, [+20, +38], 95% separable | **+31.708, CI [+29.79, +33.58]** | **HELD**, point estimate within 1.7 |
| §5a cr9114 `measured − zero` | +6, [+1, +13], 75% separable | +0.083, CI [−0.81, +0.96] | **FALSIFIED** — outside the interval, null |
| §5a trpb `measured − zero` | not a null; 70% that \|Δ\|>2 and p<0.05 | +1.938, p=0.0497 | **narrowly held on raw, NOT on d≥2** |
| §5a gb1 `measured − zero` | +2, [−3, +7], 50% separable | **+7.042, CI [+3.94, +10.31]** | interval's upper edge; separability call **right** |
| §5b creilov `measured − ablated` | +2, [−2, +6], 40% separable | −0.062 | within interval; **"not separable" was right** |
| §5b trpb `measured − ablated` | +4, [+1, +9], 60% separable | +0.479 | **FALSIFIED** — below the interval, null |
| §5b gb1 `measured − ablated` | +2, [0, +6], 45% separable | **+6.500, CI [+3.13, +9.96]** | above the interval; **separability call right** |
| §5b cr9114 `measured − ablated` | +1, [−1, +5], 30% separable | +0.042 | within interval; null as the modal call |
| §5b aggregate | `measured ≥ ablated` on ≥3 of 4; strictly better on ≥1 | 3 of 4 (+0.479, +6.500, +0.042; CreiLOV −0.062); **strictly better on GB1** | **HELD** |
| §5c which of three outcomes | **measured > ablated**, 55% | measured > ablated, separable on GB1, never worse | **HELD** |
| §5d identity | 99% byte-identical | 18,240 / 18,240 cells identical | **HELD** |

### Falsification conditions

* **F1** (CreiLOV CI contains 0) — **did not fire.** CI [+29.79, +33.58].
* **F2** (`measured − ablated` negative and separable on ≥2) — **did not fire.** No
  negative separable contrast anywhere, on either endpoint. So the recommendation does
  **not** flip to removal.
* **F3** (TrpB moves on neither contrast) — **did not fire**, by one endpoint and
  p=0.0497. The reframe survives on raw and fails at d≥2; it is not vindicated.
* **F4** (pairing broken) — **did not fire.** All four landscapes: seed block, starting
  samples and `ridge`/`random` identical across folders.
* **F5** (manifest misdescribes its own anchor state) — **did not fire** in any of 12.
* **F6** (default moved) — **did not fire.** §6.
* **F7** (anchor label mismatch) — **did not fire**; it is enforced at the call site.

### What is falsified, beyond the per-landscape bars

**The label-scale theory, and my own reframe of it.** The pre-registration's reasoning
leans on `err-SD` as the quantity that predicts the effect. It does not:

* CR9114 is **2nd of nine** on `err-SD` (13.18) and **2nd of nine** on `0-SD` (14.14), and
  it is a tight null on every contrast and both endpoints — Δ ≤ 0.083 against an MDE of
  0.707–1.294. `confirm-wt-anchor-48` found the same thing independently.
* GB1 is **7th of nine** on `err-SD` (2.53) and **9th** on `0-SD` (0.20), and it carries
  the only separable repair-specific gain in the campaign.

So **neither distance from zero nor the anchor's error predicts the damage**, and the
`err-SD` reframe in `docs/prereg-wt-anchor-fix.md` §2 is falsified as an explanatory
metric. What it *did* correctly predict is narrower and still worth keeping: that **TrpB
is not a null for the repair** (it moved, p=0.0497 raw) and that **GB1 is not a null
either** (it carries the largest repair effect measured). The two orderings do disagree,
and the `0-SD` framing would have predicted both to be nulls. The bars are not moved.

**A post-hoc mechanism, offered as a hypothesis and NOT tested here.** The two contrasts
measure different things. `ablated − zero` removes a high-leverage false point and can
only help where that point dominates the fit; `measured − ablated` adds a true observation
at the unique depth-0 corner and can only help where that observation is informative. They
fire on disjoint landscapes because those are different conditions. What would
discriminate it: a frozen-state intervention (`docs/HANDOFF.md` §5.5) measuring the
anchor's leverage on the fitted length-scale, and eqFP611-red — the only clean negative
control for the repair at 0.29 `err-SD`, which has **not** been run.

---

## 6. The identity check — the flag-absent default is untouched

`trpb`, arms `random, ridge, random_forest, panoptes_default`, 48 seeds, 5 rounds
(0–4), parent commit `fe19713` against `a032170`, same seed block
`5567fef93e8a220d892ddefca9229dd3f0280d0bb83d5f9dc50c6ab0e1c52a47`.

```
columns 24 = 3 key (method, seed, round) + 19 behavioural + 2 wall clock
             (round_seconds, cumulative_seconds, excluded and named)
rows      960 in each, 960 shared keys, 0 only-A, 0 only-B
CELLS COMPARED: 18,240   DIFFERING: 0
VERDICT: IDENTICAL on every behavioural cell. The default path is untouched.
```

Equivalent to branch `worktree-agent-a14c69a80e0afa44a`'s standard (15,840 cells over
3 arms × 48 × 5 × 22 columns, where its 22 counts the three key columns); this is 4 arms
and the key columns are the join, so they match by construction.

**One caveat, stated because it would otherwise be trusted.** The *configuration
fingerprint* DOES change — `ac75da42…` → `e6b651a1…` — because
`configuration.components.wt_anchor` is a new subtree. That is deliberate and is the
rule-2.12 fix: a zero-anchor and a measured-anchor run must not hash identically. But it
means **a fingerprint comparison across `a032170` is meaningless**, exactly as for
`e083333`. Quote the commit, not the fingerprint.

---

## 7. The amendment — what should ship, what it costs. **NOT ENACTED.**

### Status

**POST-HOC.** `docs/HANDOFF.md` §5.4 requires an amendment to be *dated before analysis*,
and that condition cannot be met: the `wt_anchor` ablation result was known before this
work began and is why it exists. A false blindness claim is worse than no amendment, so
none is made. What *is* pre-registered, before any arm ran, is
`docs/prereg-wt-anchor-fix.md`.

### It promotes nothing

`DEFAULT_WT_ANCHOR_MODE` stays `"zero"`, pinned by
`tests/test_wt_anchor_measured_mode.py::test_the_default_mode_is_still_zero` so the
one-line edit that would flip it silently fails CI. `ridge` stays the registered
reference. `DEFAULT_METHODS` is unchanged. No frozen document is edited.

### The recommendation

**Ship `measured`, not `ablated` — but not yet, and not as a default.** The evidence for
preferring repair over removal is one landscape (GB1, on both endpoints), and the evidence
that repair is *free* is one landscape (CreiLOV, where the two tie). Two of the four
measured landscapes are development data (TrpB, GB1) under rule 2.4, and **GB1 — which
carries the whole case — is one of them.** So this is mechanism evidence, not
generalization evidence.

Three things gate a default change, and two are not met:

1. **eqFP611-red must be run.** It is the only clean negative control for the repair
   (0.29 `err-SD`) and the scale theory's last testable prediction. `confirm-wt-anchor-48`
   reports it was **REFUSED by the budget gate at 71 of 81**, so the budget has to be
   fixed first. **Not met.**
2. **The three censored landscapes need a censoring decision.** CR6261-H1, CR6261-H9 and
   HER2 have a parent label at the assay floor. Anchoring at a censored bound is anchoring
   at a known-wrong value, and the current implementation cannot tell the difference. The
   honest options are to refuse `measured` on a censored parent, or to treat it as a
   left-censored observation with an appropriately inflated variance. **Not met, and not
   decided here.**
3. **Provenance.** On seven of nine landscapes the parent label is not in the prepared
   table; five need the **gitignored** `downloads/combingym/`. A default that reads it
   makes that corpus a runtime dependency of the engine's default path, which it is not
   today. The nine values must be committed — into the prepared CSVs, or a registry with
   digests — not read from `downloads/`. This compounds `docs/HANDOFF.md` §7's open item
   that the lock hashes only upstream files and only CR9114 has a digest. **Not met.**

### What it adds — a reporting requirement, not a change to the test

1. **Every production-arm figure must carry its anchor mode.** A bare `panoptes_default`
   no longer describes a configuration: one name now covers three training sets. Same
   hazard as `DPP_DIVERSITY_ENABLED` / `TWO_STAGE_GP_ENABLED` changing default on
   2026-10-31 and making "45.21" and "51.42" two different arms under one name. The mode
   is in `configuration.components.wt_anchor.mode_requested` and verified in
   `run.wt_anchor_verified`.
2. **A `measured`-mode figure must carry this disclosure.** That arm receives the parent's
   true label without spending a measurement on it, and on TrpB and eqFP611-blue that
   label is **itself a top-1% hit** (99.32nd and 99.84th percentile; GB1 97.56th, CreiLOV
   87.99th). `ridge` and `random` get no equivalent. Required phrasing: *"the engine plus
   one free parent measurement, against &lt;comparator&gt; at the full budget"* — never a
   bare win. The clean fix is to charge it (125 → 1 parent + 124), which is itself a
   protocol amendment and is not done here.
3. **Every anchor contrast must be reported on raw AND d≥2.** TrpB is the case: +1.94 raw
   at p=0.0497 and +0.42 at d≥2 at p=0.347. A raw-only contrast there would have been
   reported as a win and it is redundancy.
4. **`docs/HANDOFF.md`'s "TrpB and GB1 are the only two landscapes where this defect is
   invisible" is WITHDRAWN.** Their anchor errors are 6.81 and 2.53 label SD, the TrpB
   parent is a top-1% hit, and GB1 carries the largest repair effect measured. The
   sentence was derived from `0-SD` and the repair responds to neither metric.
5. **The label-scale theory must not be restated.** §5.

### What would have to be re-measured if the default changed

1. Every TrpB figure in `docs/HANDOFF.md` §1 "What survives and can be quoted" — 51.42
   raw, +7.75 at d≥2 (t=+13.4), +12.88 (t=+3.57, 31/16/1) against the seven-arm panel.
   **The repository's only surviving win**, and the arm mean moves 49.65 → 51.58 raw here.
2. The whole cross-protein sweep table: CreiLOV −50.29, GB1 −12.90, eqFP611-red 0 hits,
   eqFP611-blue +0.73, CR9114 +3.17, HER2 (gate-failed).
3. The seven-arm baseline panel on TrpB, CR9114-H1 and CreiLOV
   (`worktree-agent-a8325f2a02fbcdc90` @ `493826b`).
4. The `biophysical_features` one-hot arm (`worktree-agent-a14c69a80e0afa44a` @ `a779422`)
   and the combined one-hot head-to-head.
5. The spatial-veto ablation (`worktree-agent-a038e3b42e920227a` @ `3804237`).
6. The `wt_anchor` ablation itself, including `confirm-wt-anchor-48` @ `90472c2`.
7. Every `panoptes_default`-derived diagnostic: `generation_coverage.json`,
   `depth_returns.csv`, `candidate_trace.csv`, and the outputs of
   `ranking_decomposition.py` and `pool_hit_reachability.py`.
8. `docs/preregistration-cross-protein.md` §63 and §1434 — **frozen; cite a withdrawal,
   do not edit.**

And one figure that would move *favourably* and must not be quoted until re-measured
under the amendment: **GB1 goes from a non-result against ridge at d≥2 (+0.042, p=0.93) to
a separable win (+0.938, p=0.039)** in `measured` mode. On development data.

---

## 8. Corrections and findings carried out of this work

* **"11.6 units below the minimum, 8.3× the entire label range" is WRONG and is
  withdrawn.** It mixed units: 11.6 is the figure in *standard deviations*, printed as raw
  units, and 8.3× followed from it. Correct: 0.0 sits **2.793 raw units** below CreiLOV's
  minimum, which is **11.62 label SD** and **1.99×** the range. The figure the argument
  rests on is the **15.68 SD from centre**. Corrected on this branch in
  `workers/math_worker.py` and `tests/test_wt_anchor_ablation.py`, where it had
  propagated from `bb925c9`.
* **`f_best = 0.0` is the same literal-zero defect in the EI reference point**, on the
  assay's scale, and it is identical in all three modes so it cannot confound this
  contrast. On CreiLOV it starts 2.79 units below the worst observation; on a landscape
  with negative labels it would start above the best. **Not fixed** — a separate,
  independently measurable change.
* **The `zero` anchor's prior is hardcoded to `0.0`** rather than evaluated at the wild
  type. Only coincidentally equivalent; the `measured` mode evaluates it.
* **`docs/dataset-inventory.md` §7's TrpB wild-type figure 0.408085 is off by one row.**
  §2a.
* **Pre-flight gate, all 12 folders, run WITH `--run-dir`**: `8 passed, 0 FAILED, 2 not
  applicable`. The two NOT-APPLICABLE are `encoding` ("no arms were supplied") and
  `row_order.<lib>` ("needs `--invariance-runner module:function`") — both structural, and
  neither silently passing. `threshold.<lib>` now **PASSES** (main's `30c1e94`), and
  `budget.*` passes at 125 on all 48 seeds for all three arms. The campaign driver's log
  shows `corpus.lock` FAILING: that is the **environment artifact**, not rule-2.5 drift —
  `scripts/freeze_cross_protein_corpus.py` resolves `downloads/combingym` against its own
  `__file__` and `downloads/` is gitignored, so a worktree cannot see it. Re-run with
  `--corpus-dir <real repo root>\downloads\combingym` it verifies **29/29 at
  `c5dccd514a9a`** and the gate passes.
