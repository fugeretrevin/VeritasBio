# VeritasBio Engine — Per-Stage Improvement Questions, and How the Baselines Beat Us

**Purpose and audience.** Written for a reviewer expert in Bayesian optimisation, active
learning, and ML for protein engineering. It is the algorithmic companion to
`PIPELINE_FOR_BIO_REVIEW.md`, which covers the biology. Read §1 first: we have a
*verified* mechanism for how a plain ridge regression beats us on the hardest part of
the task, and it inverts the usual intuition about regularisation.

Everything numeric is measured from our own 48-seed paired runs. §2 lists five
interventions we tried and the reason each failed, so you do not re-propose them.

**Current baseline (48 paired seeds, 125 measurements, library of 159,999 variants):**

| method | best fitness | SD | top-1% found | SD | P(≥0.70) | P(≥0.90) |
|---|---|---|---|---|---|---|
| random | 0.4443 | 0.146 | 1.1 | 0.9 | 6% | 0% |
| ridge (80 one-hot) | 0.6959 | 0.232 | 31.3 | 15.7 | 67% | **17%** |
| random_forest | 0.6989 | 0.236 | 31.6 | 19.0 | 65% | **19%** |
| ALDE (published) | 0.7134 | 0.182 | 32.8 | 15.3 | 58% | **15%** |
| **ours — reliability mode** | **0.7457** | **0.035** | **45.2** | **4.8** | **98%** | 0% |
| **ours — frontier mode** (§1A) | 0.7570 | 0.092 | 39.9 | 11.3 | **98%** | **17%** |

Reliability mode dominates on breadth and consistency and never reaches the extreme
tail. Frontier mode (§1A) keeps the identical ≥0.70 success curve and reaches the tail
at parity with ridge, at a significant cost in top-1% breadth.

---

## 1. VERIFIED: how ridge reaches the tail without representing the interaction

This is the most useful thing we learned and we want it stress-tested.

### 1.1 The landscape's trap

The five best variants (of 159,999) all require **K at position 227**. But K227's
*marginal* effect, averaged over all 8,000 matched backgrounds at the other three sites,
is **−0.0206**, and it beats the parent residue in only **38.9%** of them. Its benefit
is conditional on **G at 228**: measured across all eight sub-backgrounds, K227 is
positive in every G228 context (+0.299 to +0.819) and negative in every S228 context
(−0.140 to −0.502).

So any model fitted on marginals correctly learns "K227 is bad" and never samples it.
Our engine samples K227 in **1.00%** of picks (1.77% of picks that mutate 227, against
5.26% under uniform choice).

### 1.2 What ridge actually does — observed, not inferred

The decisive fact: **G228 has the strongest positive marginal of any substitution in
the library** — mean **+0.0496**, better than parent in **69.9%** of backgrounds, ranked
first at its position by a wide margin (next best is −0.021).

So G228 is trivially discoverable by any additive model. We then replayed ridge's own
trajectory on seeds where it reached ≥0.90, refitting its coefficients round by round:

**seed 5 (reaches 1.0000):**

| round | K227 picks | K227+G228 picks | round max | fitted coef K@227 | rank | coef G@228 | rank |
|---|---|---|---|---|---|---|---|
| 0 | 1 | 0 | 0.5409 | −0.0071 | 13/20 | +0.1437 | **1/20** |
| 1 | 0 | 0 | 0.7078 | −0.0082 | 11/20 | +0.0890 | **1/20** |
| 2 | 0 | 0 | 0.6830 | **+0.0514** | **2/20** | +0.2937 | **1/20** |
| 3 | 4 | 4 | **1.0000** | **+0.3079** | **1/20** | +0.2532 | 1/20 |
| 4 | 25 | 25 | 0.8951 | +0.2596 | 1/20 | +0.1771 | 2/20 |

**seed 1 (reaches 0.9320):** same shape — G228 ranked 1/20 from round 0; K227 climbs
from rank 3 to rank 1 as KG-bearing variants accumulate.

### 1.3 The mechanism

1. G228's marginal is strongly positive, so ridge ranks it first immediately.
2. Ridge concentrates its sampling on G228-bearing backgrounds.
3. **Within that sub-population, K227's marginal effect is genuinely positive.**
4. Ridge's additive coefficient for K227 therefore flips positive — note seed 5 going
   rank 13 → rank 2 → rank 1 with *zero* K227 picks in rounds 1–2.
5. It then exploits, and reaches the tail.

**Ridge never represents the K227×G228 interaction. It cannot — it is purely additive.
It escapes the epistasis trap because greedy exploitation performs implicit
conditioning: concentration shifts its own training distribution to a slice in which
the conditional effect is marginally visible.**

Our k-DPP + Pareto selection deliberately spreads picks across feature space, so our
training distribution never concentrates, so the conditional never becomes marginally
visible. **The diversity that gives us 98% reliability at ≥0.70 is the same property
that costs us the tail.** Ridge's 17% tail rate and its 33% failure rate at ≥0.70 are
two faces of one high-variance policy.

### 1.4 Questions on this

1. Is this a known phenomenon with a name? "Implicit conditioning through greedy
   exploitation" is our description; we assume this has been studied in the BO or
   bandit literature and would like the reference.
2. Is there a principled way to get both — diversified acquisition for reliability plus
   enough concentration to expose conditional structure? Everything we tried (§2)
   failed because the acquisition function re-imposed its marginal preference on
   whatever slots we reserved.
3. Is the right answer a **two-population** design — one sub-batch under a diversified
   policy, one under a deliberately greedy policy, never mixed — rather than trying to
   tune a single acquisition function?
4. Our GP features are physicochemical property vectors (5 factors per position plus
   pairwise products). A categorical/one-hot kernel could in principle represent
   identity-specific interactions, but with 125 observations and C(4,2)×19² = 2,166
   identity pairs that seems hopeless. Is there a representation that captures
   identity-specific epistasis at this sample size, or is the only viable route the
   distributional trick ridge stumbles into?

---

## 1A. SOLVED: qualify-then-optimize closes the tail gap

A reviewer proposed, instead of our two-population design, a **lexicographic campaign
objective**: run the baseline policy unchanged until a measured variant clears a
qualifying threshold, then spend the remaining budget on frontier search. Under
deterministic reveal this preserves P(reach threshold) and time-to-threshold
*pathwise*, because every trajectory is identical to the baseline until the protected
event has already happened.

Implemented as a selector-level arm (`panoptes_qualify_then_optimize`, threshold 0.70;
frontier phase is greedy Hamming-1 local search around the incumbent, ranked by a ridge
surrogate, with no diversity, Pareto or prior terms). Core engine untouched.

**Result, 48 paired seeds:**

| | baseline | qualify-then-optimize |
|---|---|---|
| P(≥0.50) at 25/50/75/100/125 | 4 / 100 / 100 / 100 / 100 % | **identical** |
| P(≥0.70) at 25/50/75/100/125 | 0 / 88 / 92 / 94 / 98 % | **identical** |
| P(≥0.90) at 25/50/75/100/125 | 0 / 0 / 0 / 0 / 0 % | **0 / 0 / 10 / 15 / 17 %** |
| seeds reaching ≥0.80 | 3/48 | **13/48** |
| seeds reaching ≥0.90 | **0/48** | **8/48** |
| best variant found anywhere | 0.8326 | **1.0000** (global optimum) |
| top-1% hits | 45.23 | **39.90** (−5.33, t = −3.20, significant) |
| best fitness | 0.7457 | 0.7570 (+0.011, t = 0.94, **not** significant, 13–25 on seeds) |

The reliability curves are identical at every budget, as predicted by construction. The
tail gap closes to parity with ridge (17%).

**The trade is explicit.** We gain tail access and lose 5.33 top-1% variants
significantly. Mean best-fitness gain is *not* significant and loses 25 of 38 decided
seeds — the gain is concentrated in a few campaigns that reach much higher. So the
honest description is "threshold attainment at the top improves; average peak does not."

This argues for shipping **two operating modes** rather than one policy:

- **Reliability mode** — 98% reach ≥0.70, ~45 top-1% variants, SD 0.035, 0% tail.
- **Frontier mode** — identical ≥0.70 curve, 17% reach ≥0.90, ~40 top-1% variants, SD 0.092.

**Caveats.** Single landscape. The 0.70 threshold is favourable here specifically
because 88% of campaigns qualify by measurement 50, leaving three batches for the
frontier phase; on a landscape with slower qualification there would be no budget left
to switch into. And 47/48 and 8/48 both carry wide confidence intervals (Wilson 95% for
47/48 is ≈89–99.6%).

**Why this worked when our five attempts did not** (§2): all five tried to make a single
acquisition function serve two objectives, and the acquisition re-imposed its marginal
preference on every reserved slot. Switching policy at an objective boundary avoids the
conflict entirely rather than tuning through it.

---

## 2. Five interventions we tried, and why each failed

Recorded to save you proposing them. All were reverted; the committed baseline is the
best configuration we have measured.

| # | intervention | result | why it failed |
|---|---|---|---|
| 1 | Widen the planner's suggested chemistry set | no change | premise was false — the excluded residues were already being proposed at 10–16% |
| 2 | Exempt annotated active sites from the prior's charge penalties | K227 sampling 1.00% → **1.00%**; no tier hits | the prior was not the binding constraint. The GP's *learned posterior* suppresses K227 independently, and correctly, from marginal evidence |
| 3 | Depth-stratified candidate sampling toward triples/quadruples | top-1% 45.2 → **42.1** (regression) | base-rate error on our part. P(good \| depth) peaks at depth 2–3 (0.092%) and is 3.2× *worse* at depth 4 (0.028%). The tier looks quadruple-heavy only because quadruples outnumber doubles 60:1 |
| 4 | Identity-coverage floor (reserve batch slots for (position, residue) identities never measured) | completely inert; K227 1.00% → 1.00% | a 25-variant random start already yields ~95 identity observations across only 80 possible identities, so nothing is ever "novel". Coverage was never the problem — K227 *is* observed early, measured as marginally bad, and correctly avoided |
| 5 | Conditional scan: reserve 25% of the batch for carriers of the best-marginal identity | +0.0071 best fitness, **t=1.21, not significant**; mechanism unchanged (K227 1.00%, KG pair measured once) | the reserved slots were filled *by acquisition value* among carriers, so the slice inherited the exact preference it was meant to bypass. Fixing G228 and ranking its carriers yields G228-without-K227. A real conditional scan must **enumerate** a second position, not rank it |

**The pattern across all five:** our acquisition function re-imposes its marginal
preference on any slot we reserve, unless the mechanism forces enumeration. We regard
that as the central algorithmic finding, and we do not currently know the right fix.

---

## 3. The gateway finding — generation is fine, one measurement is missing

We traced every tier variant through candidate generation, varying only which leads
were supplied (the real caller passes leads built from evidence):

| leads supplied | pool depths | ≥0.70 in pool | ≥0.90 in pool |
|---|---|---|---|
| none | {1: 72, 2: 888} | 1/64 | 0/5 |
| top 8 of a random 25-variant start | {1, 2, 4} | 1/64 | 0/5 |
| **+ the double `V227K:S228G`** | {1, 2, 3, 4} | **26/64** | **5/5** |
| + the triple `F184I:V227K:S228G` | {1, 2, 4} | 7/64 | 2/5 |

**Once the double `V227K:S228G` is a lead, lead-stacking generates all five top
variants.** Its fitness is 0.6098, against a parent of 0.408. A random 25-variant start
produces leads at 0.13, 0.045, 0.043 — nowhere near it.

So the machinery works end to end. The entire tail hinges on measuring **one specific
depth-2 variant**, which *is* in the candidate pool every round, and which our
acquisition declines to pick because its K227 component has a negative marginal.

**Question 5.** Given that the unlock is a single cheap measurement that is already in
the pool, is there a decision-theoretic framing that would select it — value of
information on an *interaction term* rather than expected improvement on fitness? Our
acquisition optimises E[improvement]; the gateway's value is almost entirely that it
makes other things learnable.

---

## 4. Per-stage questions

### Step 1 — Target initialisation (UniProt/PDB retrieval)
Supplies the sequence, catalytic invariant, annotated active sites, structure.
**Open:** on a fixed library this is nearly the whole decision and it is fully
determined, so there is little to improve. **Question 6:** for a *new* library, is
there a principled way to choose which positions to saturate, given a structure and a
reaction? That is the decision with real leverage and we do not currently attempt it.

### Step 2 — LLM strategic planning
**Measured contribution: zero.** Full pipeline with live planning 0.7519/45.8; the same
engine with a deterministic planner reading the annotated positions 0.7457/45.2. Not
separable. Before we removed leaked content from the prompt templates the LLM arm
scored 0.7912 — the advantage was the leak.

We also showed the planner's output barely reaches the engine. Varying only the
strategy object: raising `max_mutation_depth` 2→4 moved quadruple generation 0%→3%, and
the free-text `epistatic_coupling` field was **provably inert** (two configurations
differing only in that field produced numerically identical output).

**Question 7:** is there any value in an LLM at the *within-library selection* stage, or
is our conclusion right that a fixed library leaves no decision for it? We think its
plausible roles are upstream (library design) or lateral (assay troubleshooting,
designing discriminating experiments), not here.

**Question 8:** we have 160k labelled fitness measurements. Is fine-tuning a model on
them plausibly competitive with ridge on 80 one-hot features at a 125-observation
budget? Our prior is no, and we would like to be wrong with reasons.

### Step 3 — Quantitative selection (the GP, acquisition, Pareto, k-DPP)
Where our advantage and our ceiling both live. **Question 9:** we have never run
selector ablations, so we cannot attribute our advantage to any component — we
previously claimed it came from Pareto + k-DPP and have retracted that as unsupported.
What is the minimal ablation set that would attribute it? Our instinct: GP vs ridge
surrogate × diversified vs greedy acquisition, 2×2, same 48 seeds.

**Question 10:** our prior mean is a hand-built biophysical score with penalties to
−4.5 against a fitness range of −0.16 to 1.0 — a scale error independent of sign. It is
fitted to nothing. Should it be (a) rescaled into fitness units, (b) replaced by a
fitted prior, (c) moved from the mean into the *kernel* as a similarity measure, or
(d) deleted? We suspect its true content is better expressed as uncertainty than as a
mean offset: K227's honest description is "conditionally beneficial, condition unknown,"
which is high variance, not low mean.

### Step 3b — Structural pre-filter
Prunes implausible candidates using the structure. **Question 11:** the burial and
cofactor-proximity terms are computed on the *parent* structure. For a substitution to
a long side chain (V→K) the parent geometry says nothing about where the new terminal
group sits. Is there a cheap way to get this right short of rotamer packing?

### Step 4 — Gatekeeper (spend authorisation)
Deterministic cost model, atomic ledger, fail-closed. A 25-variant physically measured
round costs **$3,314**; the GPU inference beside it costs ~$0.006. Works as intended.
**Question 12:** is per-variant cost the right budget unit, or should the optimiser be
aware that plate formats make certain batch sizes much cheaper per variant?

### Step 5 — ESM-2 scoring
Mean Spearman vs measured fitness **ρ = 0.344** over 20 rounds, range −0.242 to +0.682,
negative in 2 of 20. **Architecturally it is not a prior at all** in our current design —
it is computed *after* batch commitment and only recorded, so it influences nothing.
**Question 13:** at a 125-observation budget, is a frozen-embedding readout plausibly
better than 80 one-hot features, or is the dimensionality fatal? Note the tension:
two of the four optimum substitutions are reversions to wild-type residues, which a
sequence-likelihood model should favour, yet correlation is weak.

### Step 6 / 6b — Ledger, audit chain, reflection
Merkle chain verified across all 48 campaigns. Reflection's contribution is not
separable from zero. **Question 14:** token accounting is now instrumented (2.28M
tokens per 48-seed run, 384 calls, completion exceeding prompt by 40%). Is there a
reason to keep a reflection stage that cannot be shown to help?

---

## 5. Differentiators — which of these are real?

We would like these ranked by how defensible they are to a technical buyer.

| claimed differentiator | evidence | our confidence |
|---|---|---|
| **98% of campaigns reach ≥0.70 vs 58–67%** | 48 paired seeds, censoring-aware | **high** |
| **88% reach ≥0.70 at half budget vs 4–6%** | same | **high** |
| **~40% more top-1% variants** (+12.4 vs ALDE, t=5.07) | paired t-test | **high** |
| **5–7× lower variance** (SD 0.035 vs 0.18–0.24) | same | **high** |
| Real spend authorisation + Merkle audit chain (21 CFR Part 11) | verified across 48 campaigns | high, but is it *valued*? |
| Finds better single variants | **t=1.2–1.7, not separable** | **none — we have retracted this** |
| LLM-driven biophysical reasoning | contribution zero once leaks removed | **none** |
| Multi-objective Pareto (developability, synthesizability, stability) | objectives are *proxies*, never validated against measurement | **unknown — untested** |

**Question 15:** the auxiliary Pareto objectives are hand-built proxies with no
validation. Are they a differentiator, a liability, or noise? They participate in
selection, so if they are wrong they are actively harmful.

**Question 16:** for a wet-lab buyer, which is worth more — one variant at 1.0 found in
17% of campaigns, or 45 variants above the 99th percentile found in 98%? We believe the
latter and have built toward it, but we have never tested that belief against a buyer.

---

## 6. What we think the highest-impact changes are

In our order, for your critique:

1. **A two-population batch** — a diversified sub-batch for reliability plus a
   deliberately greedy sub-batch that concentrates, mimicking ridge's implicit
   conditioning, with the two never sharing an acquisition function. This follows
   directly from §1 and from why all five interventions in §2 failed.
2. **Value-of-information acquisition for interaction terms**, so a measurement whose
   value is "makes other things learnable" can be selected (§3).
3. **Selector ablations** to attribute our existing advantage before trying to extend it.
4. **Rescale or delete the hand-built prior** (Question 10).
5. **Validate or remove the auxiliary Pareto objectives** (Question 15).

Items 3–5 are cleanup with known payoff. Items 1–2 are the only ones we believe could
close the tail gap, and both came from observation rather than theory — which matters,
because the five theory-driven interventions in §2 all failed and two of them actively
degraded the system.

---

## 7. Method notes

- 48 paired seeds; identical 25-variant random start per seed shared by every method;
  4 rounds × 25; 125 unique measurements including the start; outcomes revealed only
  after a batch is committed.
- ALDE is the authors' implementation (jsunn-y/ALDE, commit `d0b3593`) in its own
  environment at their pinned versions (torch 2.1.1, gpytorch 1.11, botorch 0.9.4).
- Paired t-tests, critical |t| = 2.01 at n=48.
- Up to 871 of 159,999 fitness values (0.54%) may be KNN-imputed; our CSV does not flag
  which. All five ≥0.90 variants share the same `K227+G228` motif, so if any are
  imputed, §1 and §3 weaken. We are recovering the flags.
- Single target, single four-site library. Nothing here is claimed to generalise.
