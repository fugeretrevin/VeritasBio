# The baseline panel — results

**Pre-registration:** [`prereg-baseline-panel.md`](prereg-baseline-panel.md), committed as
`83abd97` before any arm of this panel was run or read. Read it first; this file is
scored against it.

**Runs.** Three landscapes, 48 paired seeds each, budgets matched per arm, seed as the
unit of independence, source tree `9c4cb9a`. Harness
`scripts/analysis/baseline_panel.ps1` (generality_sweep.ps1's invocation with the panel's
arms), read out by `scripts/analysis/baseline_panel_readout.py`.

| landscape | run | seeds | state | gate |
| --- | --- | --- | --- | --- |
| TrpB | `run_20261101_150147_603cd6ad` | 48–95 (held out) | SUCCESS | 12 pass, **1 FAIL: `threshold.trpb`**, 9 n/a |
| CR9114-H1 | `run_20261101_150944_7df437bb` | 0–47 | SUCCESS | 12 pass, **1 FAIL: `threshold.cr9114_h1`**, 9 n/a |
| CreiLOV | `run_20261101_151449_50e33efc` | 0–47 | SUCCESS | 12 pass, **1 FAIL: `threshold.creilov`**, 9 n/a |

Every arm on every seed of every landscape spent **all 125 distinct measurements** —
`budget.<arm>` passes for all eight arms on all three runs, and the round-0 top-1% hit
count is identical across arms per seed. The pairing premises are therefore checked, not
assumed.

> **The gate does not fully pass, and that is reported rather than smoothed.**
> `threshold.<lib>` FAILS on all three runs: `qualification.quantile` and
> `qualification.source` are not recorded in the manifest, so the bar reads as an absolute
> constant and the array it was resolved from is unknown. It is **symmetric across arms**,
> so it cannot bias a paired contrast — but an un-passed check must not be reported as a
> pass. The same check fails on the three runs `HANDOFF.md` §1's published loss figures
> come from; it was invisible there because `generality_sweep.ps1` invokes the gate
> **without `--run-dir`**, which makes `threshold` and `budget` report NOT-APPLICABLE.
> This panel runs the gate a second time against the run it produced, which is why the
> failure is visible at all.

**Reproduction check — the panel's own validity evidence.** The two landscapes with
published `panoptes_default − ridge` figures reproduce to three decimals on independent
runs of a 48-seed block:

| landscape | published (§1, commit `406ef58`) | measured here | seeds b/w/t |
| --- | --- | --- | --- |
| CreiLOV | −50.29, t=−32.01, 0/48/0 | **−50.292, t=−32.01** | 0/48/0 |
| CR9114-H1 | +3.17, t=+6.10, 41/5/2 | **+3.167, t=+6.10** | 41/5/2 |

So the panel's `ridge` and engine arms are the same arms that produced the published
numbers, and the new arms sit beside them in the same runs.

---

## 1. The panel

Primary endpoint `cum_top1_hits_excluding_initial` at the final round, as pre-registered.
`Δ` is `arm − ridge`; `ridge` is the **registered reference baseline and stays it**.

### TrpB — held-out seeds 48–95 (DEVELOPMENT DATA, rule 2.4)

| arm | mean | sd | Δ vs ridge | t | b/w/t | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `random` | 0.979 | 0.86 | −28.08 | −12.38 | 1/44/3 | loss |
| **`ridge`** | **29.062** | 16.03 | — | — | — | reference |
| `random_forest` | 26.458 | 17.29 | −2.60 | −1.07 | 24/21/3 | **not separable** |
| `gbm_onehot` | 27.208 | 16.80 | −1.85 | −0.82 | 21/23/4 | **not separable** |
| `knn_hamming` | 21.062 | 13.44 | −8.00 | −3.50 | 13/33/2 | loss |
| **`gp_onehot`** | **38.542** | 21.97 | **+9.48** | **+3.41** | 34/13/1 | **beats ridge** |
| `greedy_measured` | 21.688 | 16.12 | −7.38 | −2.97 | 16/30/2 | loss |
| `panoptes_default` | 51.417 | 5.83 | +22.35 | +9.10 | 44/4/0 | win |

`ridge − random` = **+28.08, t=+12.38, 44/1/3** — ridge is decisively better than chance here.

### CR9114-H1 (DEVELOPMENT DATA, rule 2.4 — the project's one win)

| arm | mean | sd | Δ vs ridge | t | b/w/t | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `random` | 0.875 | 0.79 | −0.96 | −2.46 | 19/21/8 | see below |
| **`ridge`** | **1.833** | 2.60 | — | — | — | reference |
| `random_forest` | 11.396 | 8.39 | **+9.56** | **+8.31** | 40/6/2 | **beats ridge** |
| `gbm_onehot` | 9.000 | 6.84 | **+7.17** | **+6.59** | 39/5/4 | **beats ridge** |
| `knn_hamming` | 10.979 | 4.46 | **+9.15** | **+12.81** | 46/1/1 | **beats ridge** |
| `gp_onehot` | 3.167 | 2.93 | +1.33 | +2.86 | 29/12/7 | beats ridge |
| **`greedy_measured`** | **11.604** | 7.32 | **+9.77** | **+9.25** | 42/2/4 | **beats ridge** |
| `panoptes_default` | 5.000 | 3.10 | +3.17 | +6.10 | 41/5/2 | win over ridge only |

**`ridge` versus `random` is the sentence this landscape turns on.** Δ = **+0.958**,
paired t = +2.46, Wilcoxon *p* = 0.044 — but **better on 21 of 48 seeds, worse on 19,
tied on 8, sign test *p* = 0.875.** The t and the Wilcoxon are carried by a few seeds
where ridge's heavy tail pays off (sd 2.60 against a mean of 1.83); on a seed-majority
basis it is a coin flip. Every other arm is separable from `random` at sign-test
*p* ≤ 1.8 × 10⁻⁶. **Ridge is not robustly better than guessing on CR9114**, which is what
`README.md` already says from 12 seeds and what 48 seeds does not overturn.

### CreiLOV (the only landscape here that is NOT development data)

| arm | mean | sd | Δ vs ridge | t | b/w/t | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| `random` | 1.167 | 0.95 | −56.06 | −34.55 | 0/48/0 | loss |
| **`ridge`** | **57.229** | 11.37 | — | — | — | reference, **and the panel's best** |
| `random_forest` | 54.958 | 17.79 | −2.27 | −0.83 | 25/22/1 | not separable |
| `gbm_onehot` | 51.479 | 13.79 | −5.75 | −2.47 | 19/29/0 | loss |
| `knn_hamming` | 32.479 | 9.91 | −24.75 | −11.96 | 2/46/0 | loss |
| `gp_onehot` | 53.708 | 10.87 | −3.52 | −1.83 | 19/26/3 | not separable |
| `greedy_measured` | 35.188 | 13.54 | −22.04 | −9.71 | 2/45/1 | loss |
| `panoptes_default` | 6.938 | 5.36 | −50.29 | −32.01 | 0/48/0 | loss |

`ridge − random` = +56.06, t=+34.55, 48/0/0. **No panel arm beats ridge here**, so ridge
is confirmed fair on CreiLOV, and the engine's −50.29 is a real loss rather than an
artefact of an unfairly strong comparator.

**Landscape informativeness.** On all three landscapes the model-based arms are separable
from `random`, so all three can discriminate methods — unlike HER2, whose endpoint cannot
(§1, §3C).

---

## 2. Survival: the engine against the strongest baseline in the panel

| landscape | best baseline | base mean | engine | Δ | t | b/w/t | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TrpB | `gp_onehot` | 38.54 | 51.42 | **+12.88** | **+3.57** | 31/16/1 | **win** |
| CR9114-H1 | `greedy_measured` | 11.60 | 5.00 | **−6.60** | **−6.25** | 10/36/2 | **loss** |
| CreiLOV | `ridge` | 57.23 | 6.94 | **−50.29** | **−32.01** | 0/48/0 | **loss** |

The best baseline is chosen by its **own mean**, as fixed in advance, not by which choice
flatters the engine. A maximum over seven arms is upward-biased by the winner's curse;
every pairwise row is in §1 and in the read-out's section B.

**CR9114 in full, because this is the project's only win:**

| comparison | Δ | t | b/w/t |
| --- | --- | --- | --- |
| engine − `ridge` | **+3.17** | +6.10 | 41/5/2 |
| engine − `random_forest` | **−6.40** | −4.84 | 16/29/3 |
| engine − `gbm_onehot` | **−4.00** | −3.90 | 14/28/6 |
| engine − `knn_hamming` | **−5.98** | −7.18 | 7/40/1 |
| engine − `greedy_measured` | **−6.60** | −6.25 | 10/36/2 |

A rule with **zero learned parameters** — rank untested candidates by Hamming distance to
the fitness-ordered measured variants — beats the engine by 6.60 top-1% hits on the one
landscape the engine won, on 36 of 48 seeds.

### 2.1 The same question with hits required to be distinct

A raw top-1% count treats N variants in one mutational neighbourhood as N wins, and three
of the arms that beat the engine on CR9114 are local-exploit rules, so this had to be
checked rather than argued. Greedy diversity-adjusted yield at Hamming separation d,
applied after measurement, thresholds from each run's own manifest, genotypes expanded to
full library length with `veritas.calibration.genotypes` (rule 2.1):

**CR9114, arm − engine:**

| arm | raw | d≥2 | d≥3 | redundancy |
| --- | --- | --- | --- | --- |
| `random_forest` | +6.40 (t+4.8) | +1.90 (t+2.7) | +0.40 (t+1.1, n.s.) | 1.82× |
| `knn_hamming` | +5.98 (t+7.2) | **+3.31 (t+5.5)** | **+1.79 (t+5.6, 32/11/5)** | 1.43× |
| `greedy_measured` | +6.60 (t+6.2) | +3.52 (t+4.9) | −0.44 (t−2.0) | 4.84× |
| `panoptes_default` | — | — | — | 1.15× |

So `greedy_measured`'s and `random_forest`'s advantages are **partly basin-concentrated**
and dissolve at d≥3 — a fair correction, and the engine's hits are the most distinct of
any arm (1.15× redundancy). But **`knn_hamming` beats the engine at every d**, so the
conclusion does not change: the engine is not the best arm on CR9114 under raw, d≥2 or
d≥3.

**TrpB, the correction runs the other way and it matters.** `gp_onehot` beats `ridge` on
the raw count (+9.48, t=+3.41) but **not once hits must be distinct** (d≥2: −0.94,
t=−2.0; redundancy 7.06× against ridge's 4.54×). The engine's advantage *strengthens*
under adjustment: engine − `gp_onehot` = **+8.69 at d≥2 (t=+14.7, 48/48 seeds)** and
+2.73 at d≥3 (45/48). On TrpB, therefore:

- on the **pre-registered primary endpoint**, ridge is *not* the strongest baseline —
  `gp_onehot` is;
- on the **diversity-adjusted endpoint**, ridge *is* the strongest baseline, and
  `gp_onehot`'s raw lead was redundancy.

Both are reported. Neither is chosen after the fact.

**CreiLOV:** the engine loses to every baseline except `random` at every d (ridge +18.15
at d≥2, t=+24.6, 48/48). Ridge stays top or tied-top; `gp_onehot` at d≥2 is +0.29, not
separable.

---

## 3. Was the pre-registration falsified?

**Partly, and the parts that failed are the parts about the baselines, not the parts about
ridge.**

| prediction | measured | outcome |
| --- | --- | --- |
| **TrpB** `random_forest` beats ridge by +5…+20 | −2.60, not separable | **falsified** |
| **TrpB** `gbm_onehot` beats ridge by +5…+20 | −1.85, not separable | **falsified** |
| **TrpB** `greedy_measured` beats ridge by +5…+25, "the strongest trivial arm" | −7.38, loss | **falsified** |
| **TrpB** `gp_onehot` ≥ ridge, 0…+8 | +9.48, beats ridge | direction right, magnitude under-predicted |
| **TrpB** `knn_hamming` loses −5…−20 | −8.00 | confirmed |
| **TrpB** `random` −10…−30 | −28.08 | confirmed |
| **TrpB headline:** engine's margin does NOT survive the panel | survives, +12.88 t+3.57 | **falsified** |
| **CR9114** RF / gbm / knn / greedy all beat ridge | +9.56 / +7.17 / +9.15 / +9.77 | confirmed (magnitudes under-predicted ~2×) |
| **CR9114** `gp_onehot` 0…+3 | +1.33 | confirmed |
| **CR9114** ridge ≈ random, not separable | +0.96, t=+2.46 but sign *p*=0.875 | split: falsified on the t, confirmed on the sign test |
| **CR9114 headline:** the +3.17 does NOT survive the panel | −6.60 vs the best arm | **confirmed** |
| **CreiLOV** no arm beats ridge at \|t\|≥2 | none does | confirmed |
| **CreiLOV** RF / gbm / gp do not beat ridge | −2.27 / −5.75 / −3.52 | confirmed |
| **CreiLOV headline:** the −50.29 is real, not a comparator artefact | confirmed, and it worsens against nothing | **confirmed** |
| **Cross-landscape:** ridge fair on CreiLOV, too weak on CR9114, TrpB intermediate | exactly this | **confirmed** |

I was **wrong about which baselines are strong**: I predicted tree ensembles and local
exploit would be the strong arms and the GP would be marginal, and on TrpB the opposite
happened — the GP was the only arm to beat ridge and the two tree models were
indistinguishable from it. I was **right about where ridge is weak**, which is the
question the panel was built to answer.

### Falsification criteria from §6 of the pre-registration

- **CR9114 — F1 TRIGGERED.** `greedy_measured` beats ridge by **+9.77**, which exceeds the
  engine's own reported margin over ridge there (+3.17), at t=+9.25 with **42 of 48 seeds
  better** (87.5%, above the 2/3 bar). Three further arms clear the same bar. **"Ridge is
  a fair reference baseline" is FALSIFIED on CR9114.**
- **CR9114 — F2 not met on its letter.** F2 required ridge to be *not* separable from
  random at |t| < 2; measured t = +2.46. On the sign test it is a coin flip
  (*p* = 0.875), so F2 is met in substance and not in form. F1 is sufficient on its own.
- **TrpB — F1 NOT triggered.** `gp_onehot` beats ridge by +9.48 at t=+3.41 with 34/48
  (70.8%) better, but the engine's margin over ridge on the same metric is +22.35, so the
  stronger arm does not close the engine's gap. Ridge is nonetheless **not** the strongest
  baseline on the primary endpoint here.
- **CreiLOV — CONFIRMED FAIR.** No arm beats ridge at |t| ≥ 2.

---

## 4. Is ridge a fair reference baseline, and where is it weak?

**It depends on the landscape, and that is the finding.** A single-reference comparison
reports a method's advantage against a comparator whose own strength varies by a factor of
six across these three landscapes relative to the achievable baseline level.

| landscape | ridge's standing | verdict on ridge as reference |
| --- | --- | --- |
| **CreiLOV** | best of seven arms; nothing beats it | **fair, and strong.** Comparing against ridge here is comparing against the best thing in the panel. |
| **TrpB** | second of seven; `gp_onehot` beats it by +9.48 raw, but ridge is top again once hits must be distinct | **adequate, with a named exception.** Any raw-count claim must carry `gp_onehot`. |
| **CR9114-H1** | sixth of seven; four arms beat it by +7 to +10, and it is a coin flip against `random` on the sign test | **NOT fair. Too weak to support a comparative claim.** |

**Where ridge is weak, mechanistically.** Ridge is additive. It is strong where the
landscape is close to additive over a wide dynamic range (CreiLOV: 15 sites, labels spanning
2.79–4.20) and weak where the assay is censored at a floor so that most of the library
carries no gradient (CR9114: −log₁₀(K_D) floored, ridge reaches 1.83 hits of a 651-hit tier
while four nonlinear or memorising arms reach 9–11.6). **The predictor of ridge's weakness
is assay censoring, not library size or site count** — CR9114 has 16 sites and CreiLOV 15,
and they sit at opposite ends of this table.

---

## 5. How many of the engine's results survive, and the split that matters

The coordinator's finding that `select_next_micro_batch` injects a fabricated wild-type
training row at the literal label `y = 0.0` regardless of assay scale forces a second cut
through the survival table. **TrpB (0.35 label SDs from centre) and GB1 (0.20) are the only
two landscapes whose assays are genuinely centred near zero** — the only two where that
defect is invisible — and they are exactly the two landscapes every quotable engine figure
in this repository comes from.

| landscape | survives the panel? | development data (2.4)? | is the label-scale defect visible here? |
| --- | --- | --- | --- |
| **TrpB** | **yes**, +12.88 raw / +8.69 at d≥2 | **yes** | **no** — 0.35 label SDs |
| **CR9114-H1** | **no**, −6.60 | **yes** | yes, 14.14 SDs |
| **CreiLOV** | **no**, −50.29 | no | yes, 15.68 SDs |

Read down the three columns together:

1. **Zero engine results survive the panel on data that is not development data.** CreiLOV
   is the only non-development landscape in this panel, and there the engine loses to every
   baseline except `random`.
2. **The one surviving result is on development data AND on the one landscape where the
   label-scale defect cannot show itself.** That is not a reason to discard it — the TrpB
   margin is large, diversity-robust, and reproduces — but it is a reason to stop
   describing it as evidence that the method generalises. It is evidence that the method
   works on a 4-site saturated library with a zero-centred assay, which is the single
   configuration in which nothing known to be broken is exercised.
3. **CR9114 is a genuine counter-example to "the WT row explains the losses".** Its scale
   violation is 14.14 SDs, among the worst, and the engine still *beats ridge* there
   (+3.17, t=+6.10). What the panel adds is that beating ridge there is a low bar: four
   standard baselines clear it by three times as much. So the WT-anchor defect and the
   weak-comparator problem are **two independent deficiencies**, and fixing either one
   leaves the other.

---

## 6. Claims that must be narrowed or withdrawn

Nothing below is enacted by this document. Each is a claim and the narrowest honest
restatement of it.

### `docs/HANDOFF.md`

| § | claim | required change |
| --- | --- | --- |
| §1 table | `CR9114-H1 \| +3.17 \| +6.10 \| 41/5/2 \| win — but development data` | **Narrow.** Add: *win over `ridge` only. Against `random_forest`, `gbm_onehot`, `knn_hamming` and `greedy_measured` the engine LOSES by 4.00 to 6.60 hits (t −3.90 to −7.18); `ridge` itself is a coin flip against `random` on the sign test (p=0.875).* |
| §1 | "One win on development data, one tie, three losses." | **Narrow.** The win is against one weak comparator. Restate as *"no win against the strongest available baseline on any landscape, and none at all on data that is not development data."* |
| §1 quotable | `panoptes_default ("reliability"), raw top-1% hits = 45.21` | **Withdraw as current.** On the same held-out seed block 48–95 at today's defaults, `panoptes_default` measures **51.417**. 45.21 is a different arm under one name — exactly what §1's own warning describes. Re-quote with the commit. |
| §1 quotable | `vs ridge, diversity-adjusted at d≥2: +7.93, t=+17.1` | **Restate.** Today's arm measures **+7.75, t=+13.4** (47/1/0) on the same block. Also note ridge is *not* the strongest raw-count baseline on TrpB. |
| §1 quotable | `at a 50-measurement budget, 6.3× the best baseline, decaying to 1.7× at 125` | **Withdraw pending re-derivation.** 6.3× is `frontier − dkl,dpp,gflownet` (14.12) divided by **ALDE** (2.25), and at budget 50 ALDE is **not** the best baseline — `ridge` is, at 2.83, which gives 5.0×, and `panoptes_default` itself gives 10.58/2.83 = **3.7×**. `budget_report.txt` fixes its comparator at budget 125 and then quotes ratios at every budget, which is not what "the best baseline" says. The panel adds `gp_onehot`, so the comparator must be re-chosen per budget. |
| §2.11 | the rule itself | **Reinforce, no change.** CR9114's `ridge − random` is a new instance: t=+2.46 and Wilcoxon p=0.044 against a sign test of p=0.875 on 21/19/8. Two tests on one contrast disagree; the verdict must carry both. |
| §7 Phase 4 | "Validated ridge, random forest, random, other compatible baselines…" | **No change needed — it already asks for this.** Record that `random_forest`, `gbm_onehot`, `knn_hamming`, `gp_onehot` and `greedy_measured` now exist and that the cross-protein sweep ran without them. |

### `README.md`

| where | claim | required change |
| --- | --- | --- |
| "CR9114-H1 antibody — 12 paired seeds" table | `engine vs ridge \| top-1% \| +4.67 \| +7.39 \| 12/0/0 \| engine better` | **Narrow to "engine vs ridge".** At 48 seeds the engine loses to four standard baselines on this landscape. The row is not wrong; the impression it leaves is. |
| same section | "on this landscape ridge, the standard baseline, is **no better than guessing**" | **Keep and strengthen the test.** At 48 seeds: Δ=+0.96, t=+2.46, Wilcoxon p=0.044, **sign p=0.875, 21/19/8**. State the two tests' disagreement rather than only the one that supports the sentence. |
| "The live agent on TrpB — 12 paired seeds" | `live agent vs ridge \| top-1% \| +20.25` | **Add the comparator caveat.** On TrpB `ridge` is not the strongest baseline on raw count — `gp_onehot` beats it by +9.48 (t=+3.41). The agent's lead survives, reduced. |
| "Across six landscapes" table | `Δ vs ridge` | **Add a column or a footnote** naming the panel's best baseline per landscape. The table's own header is honest; the surrounding prose reads it as "the baselines". |
| Baselines paragraph (~line 426) | "Baselines are the published implementations or standard ones: scikit-learn `Ridge(alpha=1.0)` and `RandomForestRegressor`…" | **Extend** to list the five arms added here, and say plainly that the six-landscape sweep used only `random`, `ridge` and the engine. |

---

## 7. The amendment this finding requires — DRAFT, NOT ENACTED

Writing the text, not applying it. Enacting it needs a human decision, because it touches
a frozen protocol.

> ### Amendment 1 to `docs/preregistration-cross-protein.md`, dated 2026-11-01
>
> **Status: POST-HOC. Outcomes have been inspected.**
>
> §5.4 of `docs/HANDOFF.md` permits an amendment "provided the original is archived, the
> amendment is dated before analysis, and results are reported whichever covariate is more
> favourable". **The second condition cannot be met here**: the panel of
> `docs/prereg-baseline-panel.md` has been run and read. This amendment is therefore
> explicitly post-hoc and must be labelled so wherever it is cited. It is written this way
> rather than back-dated because a blindness claim that is not true is worse than no
> amendment.
>
> **What is NOT amended.** The registered reference baseline remains the one-hot ridge
> regression `sklearn.linear_model.Ridge(alpha=1.0)`. No baseline is promoted. `DEFAULT_METHODS`
> is unchanged. Every figure already reported as `arm − ridge` remains valid as a statement
> about ridge, and none is recomputed against a different comparator.
>
> **What is amended — a reporting requirement, not a change of test.** From this date, a
> comparative claim about an engine arm on landscape *L* may be stated only if it carries,
> in the same table or the adjacent sentence:
>
> 1. the paired contrast against `ridge`, as now; **and**
> 2. the paired contrast against the highest-mean arm of `BASELINE_PANEL_METHODS` measured
>    on *L* at matched budget, chosen by that arm's own mean and not by favourability; **and**
> 3. `ridge − random` on *L* with **both** the paired t and the sign test, because the two
>    disagree on CR9114-H1 and a reader cannot tell which was used; **and**
> 4. where the landscape's panel has been measured, the diversity-adjusted contrast at
>    d≥2, because the raw and adjusted orderings disagree on TrpB.
>
> A claim for which (2) has not been measured is reported as **"advantage over ridge,
> comparator strength unmeasured"** — never as an advantage over "the baselines".
>
> **Rationale.** `ridge`'s strength relative to the achievable baseline level is not
> stable across landscapes: it is the best of seven arms on CreiLOV, second of seven on
> TrpB, and sixth of seven on CR9114-H1, where four arms beat it by +7.17 to +9.77 at
> t ≥ 6.59 and it is itself a coin flip against `random` on the sign test. A single-reference
> design therefore produces advantages that are not comparable across landscapes, and
> aggregating or narrating them as a single result is not supported.
>
> **Why a reporting requirement and not a new reference.** Changing the reference would
> amend the registered test after its outcome is known, which is the thing §5.4 guards
> against. A reporting requirement narrows what may be claimed from the existing test
> without altering the test, so no published figure becomes irreproducible and nothing is
> re-scored against a comparator chosen with knowledge of the result.
>
> **Scope.** Applies to `README.md`, `docs/HANDOFF.md`, and any new analysis. It does not
> require re-running the six-landscape sweep; it requires that the sweep's rows be labelled
> for what they measured, which is a comparison against one baseline.

---

## 8. Housekeeping and known defects in this work

- **`threshold.<lib>` fails the gate on all three runs** (§ above). Symmetric across arms;
  cannot bias a paired contrast; still an un-passed check.
- **An orphaned run directory, `run_20261101_144437_27fc1a1f`, holds one file and no
  manifest.** It is a partial launch of this panel that a task-stop ended — *a task-stop is
  not a process kill* (`day1-reachability.md` §7). Left in place rather than deleted, and
  recorded here so nobody reads it as a result.
- **Two other sessions were running scored sweeps in the shared results root during this
  work**, at one point two concurrently on different landscapes. This panel's three runs
  were executed with nothing else running — the harness's own guard verified it at launch,
  and all three report `SUCCESS` with full budgets — but the concurrency is noted because
  it is the condition that voided ten run directories before.
- **`panoptes_default` is the engine as shipped**, which on CreiLOV and CR9114 includes the
  wild-type-anchor training row described in §5. The baseline-versus-baseline comparisons —
  the substance of this panel — do not touch the engine's planner and are unaffected by it.
- From a worktree, `scripts/freeze_cross_protein_corpus.py verify` and
  `veritas.preflight` both default their corpus directory to the checkout they stand in and
  report all 29 pinned files `absent`. Against the real repository root the same verifier
  reports **CORPUS VERIFIED: 29/29 @ `c5dccd514a9a`** (rule 2.5 holds). Rule 2.2's parser
  gate also passes: TrpB 168/2242 = 7.4933%, GB1 72/2167 = 3.3226%.
