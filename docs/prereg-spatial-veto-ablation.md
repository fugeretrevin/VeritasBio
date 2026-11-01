# Pre-registration — is the spatial viability veto the CAUSE of the engine's losses?

**Written 2026-11-01, before any ablated arm had been run or read.** Committed before the
campaign started. Branch `worktree-agent-a038e3b42e920227a`.

The hypothesis under test, stated so it can fail:

> **H-VETO:** `workers/spatial_filter.py::filter_spatially_viable_variants`, applied to the
> candidate pool inside `select_next_micro_batch`, is the cause of the three measured
> losses of `panoptes_default` against one-hot ridge on eqFP611-red, CreiLOV and GB1.

---

## 1. The arms, the metric, and the protocol

Endpoint: **`cum_top1_hits_excluding_initial`** at the final round, which is the metric
`docs/HANDOFF.md` §1 reports Δ on. 48 paired seeds, `--base-seed 0`, budgets matched per
arm, **seed as the unit of independence**. Read out with `scripts/analysis/paired_stats.py`
(HANDOFF §4's recipe: paired difference over seeds, 20,000-sample whole-seed bootstrap CI,
Wilcoxon signed-rank, sign test, W/L/T) — no second recipe is invented.

| arm | invocation |
| --- | --- |
| `ridge` | one-hot ridge baseline, present in **both** runs per landscape |
| `panoptes_default` **unablated** | no `--ablate-components` |
| `panoptes_default` **ablated** | `--ablate-components spatial_filter` |

`--ablate-components` is a single global flag applied to every panoptes arm in a run, so
the ablated arm cannot share a run directory with the unablated one. Two runs per
landscape, same seed block. This is legitimate pairing because the manifest states the
round-0 rule: `random.Random(seed).sample over the eligible set, identical for every
method within a seed`. **Pairing check, pre-registered:** `ridge` must come out
**bit-identical across the two runs of a landscape**. If it does not, the two runs are not
paired and nothing below may be quoted.

Budgets, matching `scripts/analysis/generality_sweep.ps1`'s registered plan exactly:

| landscape | initial | batch | rounds | budget |
| --- | --- | --- | --- | --- |
| eqfp611_red | 17 | 16 | 4 | 81 |
| creilov | 25 | 25 | 4 | 125 |
| gb1 | 25 | 25 | 4 | 125 |

HANDOFF rule 2.7: `veritas.preflight` runs at the budget actually about to run, on every
landscape, and a landscape whose gate fails is **not quoted**.

HANDOFF rule 2.12: the ablated arm's manifest claim is verified against the module, not
the CLI flag. `workers/math_worker.py` has **no module-level global** gating
`spatial_filter` (the only such globals are `DPP_DIVERSITY_ENABLED` and
`TWO_STAGE_GP_ENABLED`, both False), so the single gate is `if "spatial_filter" in
ablated` at `math_worker.py:2233`. The manifest's `components.states.spatial_filter` is
therefore honest for this component — but it is still checked empirically on one seed
end-to-end before the campaign burns, by the presence of the
`spatial viability filter ABLATED` log line and the absence of veto calls.

## 2. The published baseline this is measured against

From `docs/HANDOFF.md` §1, re-derived from the prior sweep's `round_metrics.csv` on the
same metric (figures reproduce to four decimals):

| landscape | run dir | panoptes_default | ridge | Δ (published) |
| --- | --- | --- | --- | --- |
| eqfp611_red | `run_20261101_024225_8fa1d86d` | **0.0000** | 15.3542 | −15.3542 |
| creilov | `run_20261101_024544_e82e6214` | 6.9375 | 57.2292 | **−50.2917** |
| gb1 | `run_20261101_025014_5e79e3ae` | 20.3958 | 33.2917 | **−12.8958** |

## 3. The reasoning the predictions come from: is the ceiling BINDING?

`scripts/analysis/spatial_veto_profile.py`, re-run here over all three full universes
before this document was written, gives the veto's **ceiling on reachable hits** — the
count of each landscape's top-1% hits that can pass the filter at all, before any ranking
happens. The engine's achieved hits are set against that ceiling:

| landscape | universe | hits | hits that CAN pass | engine achieved | **fraction of its own ceiling used** |
| --- | --- | --- | --- | --- | --- |
| eqfp611_red | 8,191 | 82 | **0** (0.00%) | 0.00 | ceiling exactly binding |
| creilov | 165,427 | 1,655 | 250 (15.11%) | 6.94 | **2.8%** |
| gb1 | 149,360 | 1,494 | 545 (36.48%) | 20.40 | **3.7%** |

**A ceiling that is not binding cannot be the cause.** On eqFP611-red the ceiling is
*exactly* 0 and the engine scores *exactly* 0, so the veto is sufficient on its own to
explain that result. On CreiLOV and GB1 the engine reaches under 4% of the hits the veto
already lets through, so the veto is not what is stopping it.

A second measurement cuts the same way. The veto does not only remove hits, it changes the
**hit density** of the pool, and on two of three landscapes it *raises* it:

| landscape | hit density, full universe | hit density, veto-passing pool | effect of the veto |
| --- | --- | --- | --- |
| eqfp611_red | 1.00% (82 / 8,191) | **0.00%** (0 / 211) | destroys the pool |
| creilov | 1.00% (1,655 / 165,427) | **2.81%** (250 / 8,889) | **enriches 2.8×** |
| gb1 | 1.00% (1,494 / 149,360) | **1.78%** (545 / 30,594) | **enriches 1.8×** |

So on CreiLOV and GB1, ablating the veto *dilutes* the candidate pool, and the naive
"ablation lifts the ceiling, therefore helps" expectation has the sign in doubt.

Third, the mechanism is **not the same filter on all three**. `filter_spatially_viable_variants`
is three rules, and `--ablate-components spatial_filter` disables all three. Their measured
contribution to rejections differs completely:

| landscape | cohesion (>16 Å mean / >24 Å max) | buried-core | steric | first biting depth |
| --- | --- | --- | --- | --- |
| eqfp611_red | 98.4% of rejections | 1.6% | 0 | d1 (cohesion alone: d2) |
| creilov | **100.0%** | 0 | 0 | d2 |
| gb1 | **0.0% — never fires at any depth** | 63.4% | 36.6% | **d1** |

**On GB1 the depth-monotone cohesion veto described in the task framing contributes
literally zero rejections.** All 118,766 GB1 rejections are single-site chemistry calls
(`V39D` → "Electrostatic desolvation penalty") that fire at depth 1 and are not monotone
in depth. GB1 therefore tests a different mechanism from the other two, and a result there
must not be reported as evidence about the cohesion veto.

## 4. Numeric predictions — per landscape, committed

Metric is `cum_top1_hits_excluding_initial`, 48 paired seeds.

### eqfp611_red (budget 81)

| quantity | prediction |
| --- | --- |
| unablated mean | **exactly 0.0000** (ceiling is 0; forced) |
| ablated mean | **strictly > 0**; point **2.5**, interval **[0.5, 8.0]** |
| Δ(ablated − unablated) | **> 0**, seeds better ≥ **30/48**, seeds worse **0/48** |
| Δ(ablated − ridge) | still a loss: **[−15.0, −7.0]**, t < **−5**, seeds worse ≥ **40/48** |

Rationale: with the veto off the pool is the full 8,191 at 1.00% hit density, where
`random` already scored 0.6875. Any ranker above chance must exceed 0. Reaching ridge's
15.3542 additionally requires ridge-quality ranking, which is a separate question the veto
cannot answer.

### creilov (budget 125)

| quantity | prediction |
| --- | --- |
| unablated mean | 6.94 ± 1.5 (reproduces) |
| Δ(ablated − unablated) | small: **[−5.0, +8.0]**, point **+1.0**, **\|t\| < 4** |
| Δ(ablated − ridge) | **[−55, −38]**, **\|t\| > 15**, seeds worse ≥ **46/48** |
| gap to ridge surviving ablation | **> 38 hits, i.e. > 75% of the 50.29** |

### gb1 (budget 125)

| quantity | prediction |
| --- | --- |
| unablated mean | 20.40 ± 2 (reproduces) |
| Δ(ablated − unablated) | small: **[−6.0, +6.0]**, point **−1.0** (dilution) |
| Δ(ablated − ridge) | **[−18, −7]**, seeds worse ≥ **38/48** |
| gap to ridge surviving ablation | **> 7 hits** |

## 5. What would FALSIFY H-VETO, and what would confirm it

The falsification condition is stated on the **fraction of each landscape's published gap
to ridge that ablation closes**, because "the veto is the cause of the loss" is a claim
about the loss, not about the direction of some movement. Let

    closed = Δ(ablated − ridge) − Δ(unablated − ridge)
    fraction_closed = closed / |Δ(unablated − ridge)|

**H-VETO is CONFIRMED on a landscape** iff `fraction_closed ≥ 0.50` **and** the paired test
on `ablated − unablated` has seeds-better ≥ **36/48** with **|t| > 3**. Half the gap, with
a paired test, or it is not a cause.

**H-VETO is FALSIFIED on a landscape** iff `fraction_closed < 0.50`. Thresholds, computed
now so they cannot be moved later:

| landscape | published Δ vs ridge | H-VETO needs ablated Δ vs ridge ≥ | equivalently, ablated mean ≥ |
| --- | --- | --- | --- |
| eqfp611_red | −15.3542 | **−7.6771** | **7.68** |
| creilov | −50.2917 | **−25.1458** | **32.08** |
| gb1 | −12.8958 | **−6.4479** | **26.84** |

**H-VETO as stated — the veto is the cause of all three losses — is falsified if it fails
on any one of the three.**

A third verdict is available and is pre-defined rather than invented afterwards:
**CONTRIBUTING CAUSE** iff `0.10 ≤ fraction_closed < 0.50` with seeds-better ≥ 30/48, and
**NOT A CAUSE** iff `fraction_closed < 0.10` or the paired test is not separable
(|t| < 2). A landscape where ablation makes things *worse* (`fraction_closed < 0`) is
reported as **NOT A CAUSE, and the veto was net beneficial there**.

### My own pre-registered expectation, which is itself falsifiable

I expect **H-VETO to be FALSIFIED**, in this specific shape:

1. **eqfp611_red — the veto is the complete and sufficient cause of the 0-hit result**
   (ceiling exactly 0, engine exactly 0; this is arithmetic, not a measurement), but
   **NOT the cause of the full loss to ridge**: I predict the ablated mean lands near 2.5,
   below the 7.68 confirmation bar, so `fraction_closed ≈ 0.16` → **CONTRIBUTING CAUSE**.
2. **creilov — NOT a cause.** `fraction_closed < 0.10`, gap to ridge stays > 38.
3. **gb1 — NOT a cause**, and the cohesion rule under test does not even fire there.

**What falsifies my counter-expectation** (recorded so I cannot claim it either way):
if CreiLOV's ablated mean reaches **≥ 32.08**, or GB1's reaches **≥ 26.84**, my
ceiling-is-not-binding analysis is wrong, H-VETO is confirmed on that landscape, and I
will say so plainly.

## 6. Structural defects that limit what an ablation on CreiLOV can mean

`969dd78` found and did **not** fix two structure-mapping defects. Both are re-measured
here. They are recorded in this pre-registration because they change the *interpretation*
of a CreiLOV result, and that interpretation must be fixed in advance:

- **CreiLOV reads `8QI9`, which is phototropin/CrPhotLOV1, a non-aligning homolog.**
  Re-measured: 109 C-α parsed, chain A; **11 of 15 assayed positions mapped, and all 11
  mismatch the library's own parent residue** (29:L≠A, 34:L≠G, 47:Y≠Q, 60:L≠R, 61:Q≠D,
  92:K≠I, …). Positions 3, 4, 5, 7 are absent, so pairs touching them are silently dropped
  from the cohesion mean, which *weakens* the veto. Production passes no
  `reference_sequence`, so `SpatialContactMap`'s own 15% identity guard
  (`spatial_filter.py:221-243`) **never runs** and the mismatch is not caught at runtime.
- **HER2 reads chain A of `1N8Z`**, the trastuzumab **light** chain, while its library
  positions are CDR-**H3** on the heavy chain; 9 of 10 mapped positions disagree with the
  parent. HER2 is not one of my landscapes and its gate fails on informativeness anyway
  (HANDOFF §1), so it is reported as a defect, not as a result.

**Pre-registered interpretation rule.** On CreiLOV, any ablation result is a measurement of
*"the veto computed from the wrong protein's geometry"*, not of *"the veto"*. The two
conclusions must be kept apart in the report:

- **"the veto is wrong here"** — the distances are a different protein's, so whatever the
  veto removes on CreiLOV it removes for a reason with no biophysical content. This is
  supportable from the structure measurement alone, and is already established above.
- **"the veto is harmful in principle"** — needs a landscape where the geometry is
  *correct*. Of my three, only **GB1** has all assayed positions mapped with every residue
  identity matching the parent — and on GB1 the cohesion rule never fires. **TrpB (8EH0,
  62.6% ceiling, correctly mapped) is the only landscape that could test the cohesion veto
  on correct geometry, and it is not in this campaign.** That limitation is recorded now,
  before the numbers, so the report cannot overclaim.

eqFP611-red's `3M24` is in between: 11 of 13 mapped, **1** identity mismatch (174:I≠A).
Its 0-of-82 ceiling is therefore a near-correct-geometry result, which is what makes it the
strongest of the three.

## 7. What is NOT being claimed

- Nothing here is a tuning decision. If ablation closes a gap, the honest default is a
  separate question answered in the report, because the veto encodes a real biophysical
  constraint and turning it off by default is a scientific claim.
- GB1 and CreiLOV are development data (HANDOFF rule 2.4). Nothing here is generalization
  evidence.
- No outcome is permitted to reach selection (rule 2.9). The veto profile and the ceiling
  arithmetic above are retrospective analysis of full universes, computed outside any
  campaign, and no campaign reads them.
