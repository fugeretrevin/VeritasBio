# Pre-registration — the one-hot surrogate arm on CreiLOV, GB1 and TrpB

**Written and committed before any result from the runs it describes was read.** Code
commit under test: `e6c11c7`. Branch `worktree-agent-a14c69a80e0afa44a`.

The arm: `panoptes_onehot_surrogate` — `panoptes_default` with the ablation component
`biophysical_features` removed, so the GP surrogate's fit and predict inputs are one-hot
over the library's declared per-site alphabets (the same encoder `_one_hot_matrix` hands
the baselines) instead of `featurize_variant_on_positions`. Everything else is held fixed:
the calibrated biophysical prior, `AssayAdaptedHead`, the k-DPP similarity kernel, the
candidate generator, the acquisition function and the Pareto stage are untouched, and the
default arm's GP is handed the same array object it always was.

---

## 1. What this is NOT a prediction about

The motivating evidence is **retrospective Spearman rho over a fixed design matrix**,
produced by `scripts/analysis/surrogate_headroom.py`'s protocol: train on a *uniform
random* sample of the library the size of a real budget, predict 3,000 *uniform random*
held-out variants, report rho against truth, 10 repeats on shared splits.

Two things follow, and both constrain what may be quoted:

1. **The absolute rho figures are not campaign results and must never be reported as
   such.** That protocol uses **per-landscape budgets** (19 for CR6261-H1, 81 for
   eqFP611, 125 elsewhere) and a Matern kernel with a **fitted `WhiteKernel`**, where
   `build_gp_model()` uses a fixed `alpha=1e-2`. The noise model differs, so the levels do
   not transfer even though the paired comparison does. The measured pairs — one-hot
   against the engine encoding, same splits, model held fixed — were
   eqfp611_red 0.485/0.419, trpb 0.174/0.149, gb1 0.306/0.199, cr6261_h1 0.743/0.426,
   creilov 0.914/0.837. One-hot better on 8 of 9 landscapes, `|t|>2` on 7, never
   significantly worse.

2. **A better surrogate over a fixed matrix is not a better campaign.** The engine selects
   its own training distribution. The headroom protocol's training set is uniform over the
   library; a campaign's is the engine's own concentrated selection — HANDOFF records the
   engine selecting at 0.68–0.71× the universe's mean depth, and on eqFP611-red the
   candidate set it proposed had a visibly narrower fitness spread than the library, which
   is why an in-run rho of −0.0069 was uninformative. Correlations measured inside a
   narrow self-selected shell are attenuated, and the attenuation need not be symmetric
   between the two encodings: one-hot's advantage is resolving *which residue at which
   site*, which needs the training sample to vary across residues per site. A shallow,
   self-selected sample varies less. A previous investigation in this repository reached
   exactly this conclusion from the other direction — the featurisation adapter it blamed
   had never run, and the cause was the training distribution the engine selects for
   itself.

So **a campaign result that contradicts the rho ranking is a plausible and interesting
outcome, not a bug.** Section 4 says in advance which.

---

## 2. Protocol, frozen

- Landscapes and budgets exactly as `scripts/analysis/generality_sweep.ps1` registers them:
  **creilov** initial 25 / batch 25 / 4 rounds (budget 125), **gb1** 25 / 25 / 4 (125),
  **trpb** 25 / 25 / 4 (125).
- `--seeds 48 --base-seed 0`, held-out seed block, seed is the unit of independence.
- `--methods random,ridge,panoptes_default,panoptes_onehot_surrogate` — all four in **one
  run folder per landscape**, so every arm sees the same shared starting sample per seed
  and every comparison is paired **within** seed. No cross-run pairing.
- `--provenance-mode acknowledge-unflagged`, and the seven `--acknowledge-derived`
  reasons the sweep script carries.
- **Pre-flight gate (rule 2.7) at the budget about to run, before any number from that
  landscape is quoted.** A failed gate makes the landscape un-quotable, not the sweep
  worthless.
- Primary endpoint: `cum_top1_hits_excluding_initial` at each (method, seed)'s final
  round, read out by `scripts/analysis/paired_stats.py`. Verified against the published
  CreiLOV run `run_20261101_024544_e82e6214`, which reproduces HANDOFF's −50.2917 at
  0/48/0 on this metric.
- Verdicts by **paired test** (rule 2.11): paired mean difference, 95% whole-seed
  bootstrap interval, Wilcoxon signed-rank, sign test, and seeds better/worse/tied. A
  difference of means is not a verdict.
- Δ is always written as `arm − comparator`.

---

## 3. The predictions

Point predictions with 90% intervals, on the primary endpoint, 48 paired seeds.

| landscape | Δ vs `panoptes_default` | Δ vs `ridge` | predicted verdict |
| --- | --- | --- | --- |
| **creilov** | **+1.0**, [−3, +8] | **−49**, [−54, −41] | not separable from `panoptes_default`; **the gap to ridge is NOT closed** |
| **gb1** | **−1.0**, [−6, +3] | **−14**, [−19, −9] | not separable, sign weakly negative |
| **trpb** | **0.0**, [−5, +5] | ridge loses on TrpB; Δ vs ridge stays positive | not separable; **a new arm, reported as such** |

### Why these

**CreiLOV — the gap is too large to be a surrogate problem.** Ridge scores 57.23 against
the engine's 6.94 at a 125 budget, on every one of 48 seeds. At most 100 measurements are
selected, so ridge converts well over half of its selections into top-1% hits and the
engine converts about one in fifteen. A ranking improvement worth +0.077 rho on a uniform
sample cannot move ~50 hits. Whatever binds here is *what is on the ballot* or *how the
non-GP terms rank it*, not the GP's column basis. HANDOFF already records CreiLOV as
reachability-limited at 13.3%. The small positive point estimate reflects that CreiLOV is
where one-hot's fixed-matrix advantage is largest in absolute rho (0.914 vs 0.837) and
where the declared alphabets are narrow (35 one-hot columns over 15 sites), so the
estimation problem is mild — if the effect transfers anywhere, it transfers here.

**GB1 — predicted sign is negative.** GB1 is four saturated sites, so one-hot is **80
columns** against the biophysical **26**; at a 125 budget, 80 columns under an isotropic
Matern is a materially harder estimation problem. More importantly, GB1 is the one
landscape where the biophysical basis is recorded as carrying genuine signal: under a
*linear* readout those features beat one-hot there by **+0.103 rho** (t=+6.6, 10 of 10
repeats). Removing them from the GP therefore removes real GB1 information, and the
one-hot fixed-matrix gain (0.306 vs 0.199) has to pay for that.

**TrpB — the negative control, and a different question.** Same geometry as GB1 (4
saturated sites, 80 one-hot columns), but it is the landscape where the engine *wins* and
where acquisition, not generation, is the recorded bottleneck — so it is the one place a
surrogate change could plausibly matter. It is also the landscape with the least headroom:
both encodings sit near the floor (0.174 vs 0.149, a +0.025 difference). Prediction:
indistinguishable.

**On TrpB this arm is a NEW ARM, not a correction.** The TrpB figures that survive in
HANDOFF (45.21 raw top-1%, +7.93 diversity-adjusted at t=+17.1) were measured with the GP
on `engine_biophysical`. `panoptes_onehot_surrogate` on TrpB is a third arm beside
`panoptes_default` and `ridge`, and its number must never be substituted for, averaged
with, or described as a revision of any published TrpB figure. The negative control's job
here is to show the *default* arm is unchanged, which the flag-absent evidence does; it is
not to supply a corrected TrpB number.

---

## 4. Falsification conditions

Each is written so it can be checked mechanically and cannot be renarrated afterwards.

**F1 — primary. The CreiLOV gap closes.** `Δ(onehot_surrogate − ridge)` on CreiLOV
**> −35**, or the arm beats ridge on **≥ 24 of 48** seeds. If F1 fires, the featurisation
was the dominant cause of the single worst engine result in the project, and that is the
most important engine finding here. It must be stated with the paired test and with the
flag-absent evidence beside it, because the comparison only means anything if
`panoptes_default` in the same run folder still reproduces −50.

**F2 — a large campaign effect short of closing the gap.** CreiLOV
`Δ(onehot_surrogate − panoptes_default) ≥ +10` with Wilcoxon p < 0.01.

**F3 — GB1 moves the wrong way for the linear-readout evidence.** GB1
`Δ(onehot_surrogate − panoptes_default) ≥ +5` with p < 0.01.

**F4 — rho does transfer after all.** `Δ(onehot_surrogate − panoptes_default)` is
positive with Wilcoxon p < 0.05 on **all three** landscapes. That contradicts §1.2: it
would mean the fixed-matrix ranking advantage survives the engine's own training
distribution, and the "a better surrogate is not a better campaign" argument above would
be wrong as stated.

**F5 — the inverse, and equally informative.** `Δ(onehot_surrogate − panoptes_default)
≤ −5` with p < 0.01 on **any** landscape. A strictly better fixed-matrix surrogate making
the campaign strictly worse is direct positive evidence for the training-distribution
mechanism, and would be reported as such rather than as a failed arm.

**Not a falsification:** any result on a landscape whose pre-flight gate fails. That
landscape is un-quotable (rule 2.7) and the prediction for it is simply untested.

---

## 5. What a null result licenses — written now, so it is not a consolation prize later

If F1–F5 all fail — the arm is not separable from `panoptes_default` on any of the three
landscapes — the conclusion is specific and usable:

> On CreiLOV and GB1 the engine's deficit against one-hot ridge is **not attributable to
> the surrogate's feature basis.** The arm holds the model class, the prior, the
> acquisition, the Pareto stage and the candidate generator fixed and hands the GP the
> baselines' own encoding, and the deficit does not move.

That leaves exactly two live explanations, and they are separable by other work:

1. **The candidate pool / generation** — what is on the ballot at all. CreiLOV is
   reachability-limited at 13.3%; ridge selects from the enumerated untested eligible pool
   while the engine ranks what it generated. A separate audit of candidate-pool depth
   speaks directly to this, and a null result here *strengthens* its hypothesis by
   eliminating a competing one.
2. **The acquisition and the non-GP ranking terms** — the calibrated biophysical prior
   `m(x)` that the posterior mean is built on top of, the Pareto filter, and greedy
   EI+Thompson. HANDOFF already records greedy acquisition beating EI+Thompson by +11.27
   top-1% on TrpB (t=+12.20), and round 1 selecting **zero** of 193 offered depth-4
   candidates.

A null result here does not weaken either; it removes the surrogate basis from the list,
which is the point of running it as a paired arm instead of flipping a default.

---

## 6. Deliberately out of scope

- `controller/orchestrator.py` is untouched: it takes the engine default and therefore
  keeps the biophysical encoding. Whether its regret-stopping threshold should follow the
  arm is measured and reported separately; no decision is taken here.
- The other six landscapes. Three is what the deadline affords, and they were chosen as
  the two worst losses plus the negative control, not for a favourable result.
- Any aggregate across the three landscapes. Per-landscape only (HANDOFF §5.6).
