# VeritasBio / PanoptesBio â€” Pipeline, Data Flow, and Open Biological Questions

**Purpose of this document.** It is written for an expert reviewer with strong
molecular-biology and protein-engineering knowledge who has no prior context on this
system. It describes exactly what our pipeline does to data at each step, what each
step is supposed to add, what we have measured, and the specific questions where we
believe domain biology can tell us something our benchmarking cannot.

Everything numeric here is measured from our own runs unless marked `[literature]`.
Where we are uncertain, it says so. We would rather be corrected than agreed with.

---

## 1. What the system is

An automated **Designâ€“Buildâ€“Testâ€“Learn (DBTL)** loop for directed protein evolution.
The intended value proposition is *experiment efficiency*: choose which small batch of
variants to synthesise and assay next, so that a useful protein is found in fewer
wet-lab measurements. A measurement is expensive (our cost model prices one variant at
â‰ˆ **$133**, comprising DNA synthesis plus triplicate robotic plate assay), so the
product is fundamentally about allocating a scarce measurement budget well.

The loop has six steps, described in Â§4.

**What it is not:** it is not a structure predictor, not a folding engine, and not a
de-novo designer in the configuration under discussion here. It selects variants from a
defined combinatorial library.

---

## 2. The biological system under test

### 2.1 Target

**Tryptophan synthase Î²-subunit (TrpB)** from *Thermotoga maritima*.

- TrpB catalyses the final step of tryptophan biosynthesis: condensation of **indole**
  with **L-serine** to give **L-tryptophan**, via a **pyridoxal 5â€²-phosphate (PLP)**
  cofactor. `[literature]`
- The PLP is held as an internal aldimine (Schiff base) to a catalytic **lysine**. In
  our annotation that is **Lys83**, which we treat as an inviolable catalytic invariant
  and never mutate.
- Structure used: **PDB 8EH0**. Sequence length 397 aa, of which 1â€“389 share UniProt
  **P50909** numbering and 390â€“397 are an `LEHHHHHH` purification tag.

### 2.2 The parent is already engineered â€” and this matters

The parent is **Tm9D8\***, a laboratory-evolved thermostable TrpB variant, not
wild-type. Relative to wild-type TmTrpB (P50909) it carries ten substitutions:

```
P19G, E30G, I69V, K96L, P140L, N167D, I184F, L213P, G228S, T292S
```

**Note `I184F` and `G228S`.** The parent has moved *away* from the wild-type residues
at two of the four positions this library varies.

### 2.3 The library

A **combinatorially complete four-site library** at active-site positions
**183, 184, 227, 228**, from Johnston, Iqbal et al. (2024), bioRxiv
`10.1101/2024.06.23.600144`. `[literature]`

- Parent residues at those sites: **V183, F184, V227, S228** (`VFVS`).
- 20â´ = 160,000 combinations; our eligible set is **159,999** (the parent itself is
  excluded as it is the reference).
- Fitness is a normalised activity measure with the landscape maximum scaled to 1.0.
  Higher is better.

### 2.4 The published global optimum, and a striking observation

The best variant is **V183A : F184I : V227K : S228G** (fitness 1.0).

Two of its four substitutions â€” **F184I** and **S228G** â€” are **reversions to the
wild-type TmTrpB residues** that Tm9D8\* had replaced. The optimum partially undoes the
parent's prior engineering while adding two new changes (V183A, V227K).

We do not know how to interpret this and would value your view. Candidate readings:
(a) the earlier engineering traded activity for thermostability and this library
recovers activity; (b) 184/228 are epistatically coupled to 183/227 such that the
wild-type identities are only optimal in the presence of the new substitutions;
(c) coincidence of a rugged landscape.

### 2.5 Shape of the landscape â€” it is overwhelmingly dead

Measured over all 159,999 eligible variants:

| statistic | value |
|---|---|
| minimum fitness | âˆ’0.164 |
| **median fitness** | **0.0124** |
| mean fitness | 0.0211 |
| maximum fitness | 1.000 |
| **fraction â‰¤ 0.05** | **92.15%** |
| **parent (`VFVS`) fitness** | **0.408** |

So: **92% of the library is essentially inactive**, and the *parent itself scores 0.408* â€”
well above the median variant. Most mutation at these four sites is strongly
deleterious. Reaching 0.70 is ~1.7Ã— parent activity; reaching 1.0 is ~2.45Ã—.

Tier sizes:

| threshold | variants | fraction |
|---|---|---|
| â‰¥ 0.50 | 556 | 0.348% |
| â‰¥ 0.70 | 64 | 0.040% |
| â‰¥ 0.90 | **5** | 0.003% |
| â‰¥ 0.95 | 1 | 0.001% |

### 2.6 Hit rate by mutation depth â€” counterintuitive and important

"Depth" = number of the four sites mutated away from parent.

| depth | variants in library | n â‰¥ 0.70 | **P(â‰¥0.70 \| depth)** | n â‰¥ 0.90 | P(â‰¥0.90 \| depth) |
|---|---|---|---|---|---|
| 1 | 76 | 0 | 0.0000% | 0 | 0.0000% |
| 2 | 2,166 | 2 | **0.0923%** | 0 | 0.0000% |
| 3 | 27,436 | 25 | **0.0911%** | 1 | 0.0036% |
| 4 | 130,321 | 37 | **0.0284%** | 4 | 0.0031% |

**Doubles and triples have â‰ˆ3.2Ã— the per-sample hit rate of quadruples.** The â‰¥0.70
tier *looks* quadruple-dominated (58%) only because quadruples outnumber doubles 60:1.
We initially misread this and degraded our own system by chasing depth; see Â§7.3.

All five â‰¥0.90 variants are triples or quadruples, so the extreme tail *is*
depth-dependent even though the broad "good" region is not.

---

## 3. The data that enters the system

1. **Target metadata** â€” UniProt/PDB-derived: sequence, annotated catalytic position
   (83), annotated binding/active-site positions (183/184/227/228), PDB id (8EH0).
   This is the only *a priori* biological knowledge the system is given.
2. **The measured fitness table** â€” 159,999 variantâ†’fitness values, held by an
   evaluator that reveals a value **only after** a method has committed to measuring
   that variant. This simulates the wet lab.
3. **Structure** â€” 8EH0 coordinates, used for contact distances, burial estimates and
   cofactor proximity.

**Provenance caveat, stated up front:** Johnston et al. report that 159,129 variants had
sufficient sequencing coverage to be quantified and that **871 fitness values were
imputed** with a KNN imputer. `[literature]` Our prepared CSV does not flag which, so up
to 871 of our values (0.54%) may be model-filled rather than measured. Every run is
labelled `UNVERIFIED_INCLUDES_IMPUTED`.

---

## 4. The six steps, and what each is supposed to add

### Step 1 â€” Target initialisation (RAG over UniProt/PDB)
Fetches the metadata in Â§3.1. **Adds:** the annotated active-site positions and the
catalytic invariant. **Our note:** on a fixed library this is nearly the whole
"decision" â€” the library already determines which positions can be varied.

### Step 2 â€” Strategic planning (LLM)
An LLM (Meta/Muse `muse-spark-1.3`) is prompted with the protein context plus the
top-scoring evidence so far, and returns strict JSON: residue clusters to explore,
forbidden positions, a maximum mutation depth, and per-position "suggested
chemistries", plus free-text mechanistic hypotheses and an "epistatic coupling"
statement. **Intended to add:** biophysical reasoning about which residues and which
chemistries to try, and why.

**Measured contribution: none.** See Â§6.2. This is the central open problem.

### Step 3 â€” Quantitative selection
The engine builds a candidate pool and scores it:
- **Candidate generation:** singles, lead-stacked variants, exhaustive within-cluster
  and cross-cluster **pairs**, plus GFlowNet-sampled multi-mutants.
- **Featurisation:** per-position 5-factor amino-acid property deltas (hydropathy,
  charge, volume, polarity, aromaticity) plus global and pairwise epistatic
  interaction terms (charge products, volume products, hydrophobic clustering,
  3D-distance-attenuated coupling).
- **Surrogate model:** a heteroscedastic Gaussian Process with a deep-kernel-learning
  feature map, fitted on measured evidence, with a hand-built **biophysical prior
  mean** (see Â§5).
- **Acquisition:** Thompson sampling / expected improvement, then a **4-objective
  Pareto front** (fitness, developability, synthesizability, stability) and a **k-DPP**
  diversity selection to pick the final batch.

**Adds:** the actual batch choice, and this is where our measured advantage arises.
*Correction from an earlier draft:* we previously attributed the advantage specifically
to Pareto + k-DPP diversification. We have run no selector ablation, so that
attribution is unsupported. The measured result is a whole-system advantage with the
mechanism unattributed.

### Step 3b â€” Structural pre-filter
Prunes or penalises variants using 8EH0: steric clashes, buried-core cavity creation,
charge burial, cofactor-shell proximity, and ANM/GNM normal-mode hinge residues.
**Adds:** rejection of physically implausible candidates before spend.

### Step 4 â€” Gatekeeper (spend authorisation)
Deterministic cost model and atomic spend authorisation against a SQLite ledger,
fail-closed. Prices GPU inference in GPU-seconds and physical measurement with a
DNA-synthesis + robotic-assay model. **Adds:** a hard budget guarantee and an audit
record. A 25-variant physically measured round costs **$3,314**.

### Step 5 â€” Execution / measurement
Two things happen:
- **ESM-2 650M masked-marginal scoring** on Modal A10G GPUs â€” an *in-silico* score.
- **The assay** â€” in the benchmark, a lookup of the held-out measured fitness; in
  production, a real plate readout.

The GP learns from the **measured** values. ESM is recorded alongside, not used as the
learning target. **Adds (intended):** a zero-shot prior on sequence plausibility.
**Measured:** ESM-2 vs measured fitness, **Spearman Ï = 0.344** (mean over 20 rounds),
range âˆ’0.242 to +0.682, **negative in 2 of 20 rounds**.

### Step 6 â€” Evidence and audit
Appends `(sequence, measured_score, instrument)` to the ledger and extends an
HMAC-SHA256 Merkle chain for 21 CFR Part 11-style integrity. **Adds:** provenance â€”
which planner proposed each round, which instrument produced each number.

### Step 6b â€” Reflection (LLM)
The LLM is shown prediction residuals (observed âˆ’ GP-predicted) for positive and
negative outliers and asked for a mechanistic diagnosis plus updated/forbidden
chemistries, which are merged into the next round's strategy. **Intended to add:**
closed-loop learning from surprises. **Measured contribution:** not separable from
zero.

---

## 5. The biophysical prior â€” the part most in need of expert review

The GP's prior mean is a hand-built score. Its own docstring states it is **not** a
free energy: "a sum of hand-chosen penalties for volume change, charge burial and
proline introduction; there is no force field, no solvation model, no torsional term,
and no parameterization against any experimental ddG dataset."

Components (per mutation unless stated):

| term | magnitude | condition |
|---|---|---|
| catalytic invariant mutated | âˆ’4.5 | position in annotated catalytic set |
| buried-core cavity creation | âˆ’1.5 | volume delta < âˆ’60 Ã…Â³ at a buried position |
| buried-core overcrowding | âˆ’1.0 | volume delta > +50 Ã…Â³ at a buried position |
| **buried charge desolvation** | **âˆ’2.0** | buried, `AILVMF â†’ KRDENQ` |
| cofactor-shell steric crowding | âˆ’3.5 | < 4.5 Ã… from PLP, volume delta > +40 Ã…Â³ |
| **cofactor-shell charge introduction** | **âˆ’2.5** | < 4.5 Ã… from PLP, `AILVMF â†’ KRDENQ` |
| cofactor-shell proline | âˆ’2.0 | < 4.5 Ã… from PLP |
| entrance-channel decongestion bonus | +0.45 | 4.5â€“12 Ã… from PLP, `F/W/Y â†’ I/V/A/L` |
| entrance-channel flexibility bonus | +0.25 | 4.5â€“12 Ã… from PLP, `S/T/A â†’ G` |
| pairwise contact synergy | +0.35 | pair â‰¤ 8.5 Ã… apart |
| pairwise spatially-uncoupled | âˆ’0.5 | pair > 16 Ã… apart |
| pairwise steric clash | âˆ’2.5 | pair < 7.5 Ã…, both volume delta > +35 Ã…Â³ |
| pairwise like-charge repulsion | âˆ’1.5 | pair < 7.5 Ã…, both K/R or both D/E |
| buried core set | â€” | computed from 8EH0 (â‰¥7 neighbours within 7 Ã…) |

### 5.1 The prior ranks the true optimum as bad

Measured directly:

| variant | prior score |
|---|---|
| **V183A:F184I:V227K:S228G** (the global optimum, fitness 1.0) | **âˆ’4.35** |
| V183A:F184I:S228G (227 left as parent V) | **+0.30** |
| V183A:F184I:V227A:S228G | +0.15 |
| V183A:F184I:V227L:S228G | +0.15 |

Position **227 is buried and within 4.5 Ã… of the PLP cofactor**, so `V227K` takes both
charge penalties (âˆ’2.0 and âˆ’2.5). Single substitutions at 227 split cleanly:
`A,C,H,I,L,M,S,T â†’ 0.0`; `D,E,K,N,Q,R,W,Y â†’ âˆ’4.5`. **K ranks 14th of 19.**

Empirically, **100% of the â‰¥0.90 tier carries K227.** Our engine proposes K at 227 in
**1.0%** of its picks and leaves 227 as parent V in **44%**.

**This is our sharpest question for you.** Is "do not introduce a positive charge into
the PLP cofactor shell" sound general biophysics that happens to be wrong for this
substrate, or is it wrong as stated? Plausible mechanisms we can think of for why K227
might be *beneficial* here:
- Lysine Îµ-amino group participating in substrate positioning or proton shuttling for
  non-natural indole analogues;
- charge compensation for the PLP phosphate;
- a salt bridge or cationâ€“Ï€ interaction with the indole ring;
- opening the channel for bulkier substituted indoles (K is long and flexible).

We tried simply exempting annotated active sites from these charge penalties. It did
**not** change K227 sampling (still 1.0%) and did not improve results, so the penalty
was not the binding constraint â€” see Â§7.2.

### 5.2 RESOLVED: K227 is an obligate pair with G228

Following reviewer guidance we computed the matched contrast
Î”(X) = f(bâ‚â‚ˆâ‚ƒ, bâ‚â‚ˆâ‚„, X, bâ‚‚â‚‚â‚ˆ) âˆ’ f(bâ‚â‚ˆâ‚ƒ, bâ‚â‚ˆâ‚„, V, bâ‚‚â‚‚â‚ˆ) over all 8,000 backgrounds.

**K227's effect is determined entirely by residue 228 â€” in 8 of 8 backgrounds:**

| background (183, 184, 228) | f(V227) | f(K227) | Î” |
|---|---|---|---|
| V F **G** | 0.1904 | 0.6098 | **+0.419** |
| V I **G** | 0.5833 | 0.8825 | **+0.299** |
| A F **G** | 0.0696 | 0.4034 | **+0.334** |
| A I **G** | 0.1807 | **1.0000** | **+0.819** |
| V F **S** | parent (0.408) | 0.0629 | strongly negative |
| V I **S** | 0.5988 | 0.0969 | **âˆ’0.502** |
| A F **S** | 0.2836 | 0.1436 | **âˆ’0.140** |
| A I **S** | 0.6730 | 0.2923 | **âˆ’0.381** |

Positive in every G228 background, negative in every S228 background. Textbook sign
epistasis. All five â‰¥0.90 variants carry **both** K227 and G228:

```
AIKG 1.0000   CLKG 0.9320   VLKG 0.9226   ALKG 0.9166   CIKG 0.9098
```

With `KG` fixed, 183 âˆˆ {A, C, V} and 184 âˆˆ {I, L} are permissive â€” 5 of those 6
combinations exceed 0.90. The tail is a narrow motif with tolerant flanks.

**It is lysine-specific, not charge-general.** Arginine â€” the obvious positive-charge
control â€” has a *worse* mean effect (âˆ’0.0287 vs âˆ’0.0206) and a ceiling of only 0.699
against lysine's 1.000. So "positive charge is favourable here" does not explain it. The
repacking reading is more consistent with the data: G228 removes a side chain and adds
backbone flexibility, and only lysine's specific geometry exploits that.

**Why our engine correctly learns the wrong thing.** Across all 8,000 matched
backgrounds K227 has **mean Î” = âˆ’0.0206** and is beneficial in only **38.9%** of them.
"K227 is bad" is therefore *true on average* â€” our prior and our smoothing GP both
learned a correct marginal fact. But K227 also has the **largest upside of any residue
at that position** (max Î” **+0.819**; next best is M at +0.416). We are averaging away a
high-variance conditional. This also explains why removing the prior penalty changed
nothing: the GP's learned posterior suppresses K227 independently, and correctly, from
the marginal evidence.

This suggests the prior's error is in the **wrong moment**: K227 should carry high
*variance*, not low *mean*. A prior expressing "uncertain" rather than "bad" would leave
the GP free to explore it.

### 5.3 Depth is penalised
Because per-mutation and pairwise penalties are summed, mean prior by depth is
**âˆ’0.695 / âˆ’1.602 / âˆ’2.494 / âˆ’3.371** for depths 1â€“4, and at depth 4 the best
achievable prior was âˆ’0.05 while singles and doubles could reach +0.6. Normalising the
pairwise terms by pair count reduced but did not remove this.

### 5.4 Target-specific content we removed
Two entrance-channel bonuses rewarded `Fâ†’I` and `S/T/Aâ†’G` â€” i.e. **two of the four
optimum substitutions** â€” inside a rule whose comment named TrpB. Several prompt
templates likewise named this library's positions and known-good residues. All removed;
see Â§7.1.

---

## 6. What we have measured

### 6.1 Benchmark protocol
48 paired random seeds. Each seed: an identical 25-variant random starting sample given
to every method, then **4 rounds Ã— 25 newly chosen variants** â€” **125 unique
measurements total, including the starting 25**. Outcomes are revealed only after a
method commits a batch. Competitors:

- **random** â€” uniform sampling of untested variants
- **ridge** â€” ridge regression on 4Ã—20 = 80 one-hot features, refit each round
- **random_forest** â€” same encoding, 300 trees
- **ALDE** â€” the published implementation (jsunn-y/ALDE, Nat Commun 2025), run in its
  own isolated environment at the authors' pinned versions, commit `d0b3593`
- **panoptes** â€” our engine

Results (mean Â± SD across 48 seeds; "top-1%" = count of variants found at or above the
99th percentile of the eligible set):

| method | best fitness | SD | top-1% found | SD | seconds |
|---|---|---|---|---|---|
| random | 0.4443 | 0.146 | 1.1 | 0.9 | 0.0 |
| ridge | 0.6959 | 0.232 | 31.3 | 15.7 | 0.2 |
| random_forest | 0.6989 | 0.236 | 31.6 | 19.0 | 2.6 |
| ALDE | 0.7134 | 0.182 | 32.8 | 15.3 | 9.5 |
| **panoptes** | **0.7457** | **0.035** | **45.2** | **4.8** | 2.2 |

Paired t-tests, n=48, critical |t| = 2.01:

| comparison | metric | mean diff | t | verdict |
|---|---|---|---|---|
| vs ridge | top-1% | +13.9 | 5.93 | **significant** |
| vs random_forest | top-1% | +13.6 | 4.55 | **significant** |
| vs ALDE | top-1% | +12.4 | 5.07 | **significant** |
| vs ridge | best fitness | +0.050 | 1.48 | not separable |
| vs random_forest | best fitness | +0.047 | 1.38 | not separable |
| vs ALDE | best fitness | +0.032 | 1.22 | not separable |

**Reliability â€” P(target reached) at each cumulative budget, over all 48 seeds.**

An earlier version of this document reported median measurements-to-threshold with a
separate censored count. That was wrong: the median was taken over only the seeds that
succeeded, so it described the lucky subset and is not comparable across methods with
different failure rates. Success probability is defined for every seed and needs no
censoring correction. Corrected:

*target â‰¥ 0.50*

| method | 25 | 50 | 75 | 100 | 125 |
|---|---|---|---|---|---|
| random | 4% | 8% | 23% | 31% | 33% |
| ridge | 4% | 42% | 75% | 85% | 88% |
| random_forest | 4% | 33% | 67% | 79% | 85% |
| ALDE | 4% | 50% | 79% | 90% | 90% |
| **panoptes** | 4% | **100%** | **100%** | **100%** | **100%** |

*target â‰¥ 0.70*

| method | 25 | 50 | 75 | 100 | 125 |
|---|---|---|---|---|---|
| random | 0% | 0% | 2% | 6% | 6% |
| ridge | 0% | 4% | 23% | 48% | 67% |
| random_forest | 0% | 4% | 21% | 54% | 65% |
| ALDE | 0% | 6% | 23% | 40% | 58% |
| **panoptes** | 0% | **88%** | **92%** | **94%** | **98%** |

*target â‰¥ 0.90*

| method | 25 | 50 | 75 | 100 | 125 |
|---|---|---|---|---|---|
| ridge | 0% | 0% | 4% | 12% | **17%** |
| random_forest | 0% | 2% | 6% | 12% | **19%** |
| ALDE | 0% | 2% | 6% | 12% | **15%** |
| **panoptes** | 0% | **0%** | **0%** | **0%** | **0%** |

**Summary of our position.** At 50 measurements â€” half the budget â€” we reach â‰¥0.70 in
88% of campaigns against 4â€“6% for every competitor, and â‰¥0.50 in 100% against 33â€“50%.
We find ~40% more top-1% variants with 5â€“7Ã— lower variance. We do **not** find better
single variants than any competitor, and we reach â‰¥0.90 in **0 of 48** campaigns where
they reach it in 15â€“19%.

### 6.2 The LLM contributes nothing measurable
Full pipeline with live LLM planning and reflection: **0.7519 / 45.8**. The same engine
with a deterministic planner that simply reads the annotated active-site positions from
metadata: **0.7457 / 45.2**. Not separable.

Before we removed the leaked content described in Â§7.1, the LLM arm scored 0.7912 â€” the
apparent advantage was the leak.

We also established that the planner's output barely reaches the engine. Holding
everything else constant and varying only the strategy object:

| strategy configuration | depth-1 | depth-2 | depth-3 | depth-4 |
|---|---|---|---|---|
| max_depth = 2 | 16% | 79% | 5% | 0% |
| max_depth = 4 | 8% | 71% | 19% | 3% |
| max_depth = 4, all 20 AA offered | 0% | 92% | 8% | 0% |
| max_depth = 4, explicit "all four MUST be co-mutated" | 11% | 76% | 11% | 3% |

Raising the depth limit from 2 to 4 moved quadruples from 0% to 3%. The free-text
`epistatic_coupling` instruction was **provably inert** â€” two configurations differing
only in that field produced numerically identical output.

So the LLM is currently writing instructions to a component that largely ignores them,
at a cost of ~2.3M tokens per 48-seed run.

---

## 7. Three hypotheses we tested and falsified

Recorded so you do not spend effort on them.

**7.1 "The prompt leaks the answer."** True, and we fixed it â€” a few-shot example
tailored to the target named 3 of 4 optimum residues; a static example named avGFP's
S65T and T203Y; the system prompt cited this library's own loop positions on every
call. Worth â‰ˆ +0.045 of illusory performance. But removing it only *reduced* the LLM
arm to parity; it did not explain the ceiling.

**7.2 "The charge penalty on K227 blocks the optimum."** The penalty is real (Â§5.1), but
exempting annotated active sites from it left K227 sampling at exactly 1.00% and
produced no â‰¥0.90 hits. The penalty is not the binding constraint.

**7.3 "The candidate pool is depth-starved."** We increased deep sampling, which moved
quadruples from 0% to 17% of selections â€” and **made results worse** (top-1% 45.2 â†’
42.1). This was a base-rate error on our part: we reasoned from the tier's composition
(58% quadruples) rather than P(good | depth), which actually *peaks* at depth 2â€“3
(Â§2.6). All changes reverted.

---

## 7A. RESOLVED: the tail deficit is a generation failure, not a scoring failure

Reviewer guidance was to trace tier variants through each pipeline stage before
changing anything. Doing so answered the question immediately.

For a mid-campaign round using the deterministic planner (clusters `[183,184]`,
`[227,228]`, `max_mutation_depth = 4`):

| stage | result |
|---|---|
| candidate pool size | 960 variants |
| **pool depth composition** | **72 singles, 888 doubles, 0 triples, 0 quadruples** |
| of the 64 variants â‰¥0.70, present in pool | **1 (1.6%)** |
| of the 5 variants â‰¥0.90, present in pool | **0 (0.0%)** |

**The top five variants are never generated.** They cannot be scored, ranked, Pareto-
dominated or diversity-rejected, because they never enter the pool. All three of our
earlier hypotheses concerned scoring and selection, and were therefore testing stages
that never see these variants.

Cause: the generator emits exhaustive within- and cross-cluster **pairs** plus singles.
Triples and quadruples arrive only via lead-stacking onto previously measured
high-fitness leads, or via GFlowNet sampling â€” neither of which contributed here. The
â‰¥0.70 tier is 97% triples and quadruples, so its absence from a pairs-only pool is
structural.

**But an accessible path exists.** The double `V227K:S228G` *is* in the pool and scores
**0.6098**, well above the parent's 0.408. From there:

```
V227K:S228G            0.6098   (in pool, above parent)
  + F184I   ->         0.8825
  + V183A   ->         1.0000   (global optimum)
```

Two lead-stacking steps. So the tail is two hops from a variant the generator already
proposes.

The single weak link is step one: we sample K227 in only **1.77%** of picks that mutate
position 227 (uniform over non-parent residues would be 5.26% â€” a 3Ã— depletion), so we
rarely *measure* the gateway double, so it never becomes a lead, so stacking never fires
on it.

This reframes the problem entirely. It is not "our optimiser cannot reach the tail". It
is "our optimiser under-samples one residue whose marginal effect is genuinely negative,
and that residue is the entrance to the only productive region."

**A note on prior scale, which the reviewer raised and we confirm.** The prior's
penalties run to âˆ’4.5 while measured fitness spans âˆ’0.16 to 1.0. As a GP prior *mean*
this is a scale error independent of sign: a âˆ’4.5 prior on a variant whose true fitness
is 1.0 is not a mild preference, it is an assertion three orders of magnitude outside
the observable range. Separately, the `KRDENQ` penalty set lumps charged K/R/D/E with
neutral-polar N/Q and gives opposite charges an identical penalty; and the buried and
cofactor-shell terms appear to double-count one desolvation event (both fire on `V227K`,
totalling âˆ’4.5). We regard all three as defects.

---

## 8. What we are asking you

**Round 1 questions 1 and 3 are now answered** â€” see Â§5.4 and Â§7A. K227 is an obligate
pair with G228, lysine-specific rather than charge-general, and the tail deficit is a
candidate-generation failure. The questions below are what remain, plus new ones those
answers raised.

1. **Mechanism of the K227/G228 pair.** Our data establish *that* K227 requires G228
   (8/8 backgrounds, Î” +0.30 to +0.82 with G, âˆ’0.14 to âˆ’0.50 with S) and that arginine
   does not substitute (ceiling 0.699 vs 1.000). We cannot distinguish *why*. Does the
   repacking reading â€” G228 removing a side chain and adding backbone flexibility that
   lysine's specific length and geometry exploit â€” fit TrpB's active-site architecture?
   What would distinguish it from an electrostatic or conformational-equilibrium
   explanation using measurements we could actually make?

2. **Is "K is good only with G at the adjacent position" generalisable, or idiosyncratic?**
   If adjacent-position flexibility gating a long charged side chain is a recurring
   motif in enzyme active sites, that is a prior worth encoding. If it is specific to
   this geometry, encoding it would be overfitting to one landscape. Which is it?

3. **The parent-reversion observation (Â§2.4).** Two of four optimum substitutions revert
   Tm9D8\* to wild-type TmTrpB residues (F184I, S228G). We now know S228G is also the
   enabling partner for K227, which complicates the reversion reading â€” the ancestral
   residue may be favoured because of what it permits, not because it is ancestral. Does
   that change your view on whether a "favour ancestral identities" prior is justified?

4. **Is the â‰¥0.90 tail worth chasing at a 125-measurement budget?** Five variants in
   160,000, but we now know a two-hop path exists from a double the generator already
   proposes (Â§7A), so the cost of reaching it may be far lower than the rarity suggests.
   Does that change the answer?

5. **How should a prior express "high variance" rather than "low mean"?** Our central
   error is that K227 received a large negative prior mean when the biologically honest
   statement is "this substitution is conditionally beneficial and we do not know the
   condition." Is there an established way to encode conditional or
   context-dependent biophysical expectations that does not collapse to a marginal
   penalty?

5. **Should ESM-2 be in this loop at all?** Ï = 0.344, negative in 2 of 20 rounds. Note
   the tension: two optimum substitutions are reversions *to* wild-type residues, which
   a sequence likelihood model should favour â€” yet overall correlation is weak. Is a
   masked-marginal score simply the wrong readout for an engineered thermostable parent
   assayed on a non-natural substrate? Would a fitness-supervised readout head on ESM
   embeddings be expected to do better?

6. **Where could an LLM genuinely add value in a DBTL loop?** Our result is that on a
   *fixed* library it adds nothing, which we think is structural â€” the library already
   determines the search space, so there is no decision left to make. Candidates we see:
   choosing which positions a *new* library should saturate; proposing off-library
   allosteric or second-shell residues (our planner spontaneously proposed R51, P52,
   L77, which this four-site library cannot measure); or interpreting assay anomalies.
   Which of these is most likely to pay, biologically?

7. **What are we not measuring that we should be?** We score fitness, top-1% recall,
   measurements-to-threshold, and consistency. For a real campaign, what else matters â€”
   expression, solubility, thermostability retention, substrate-scope breadth, mutational
   robustness of the final lead?

---

## 8A. The assay definition we are missing â€” please tell us what to ask for

The reviewer correctly noted that we never specify the substrate or the readout, and
that this gates mechanistic interpretation. We do not have it: our metadata carries the
sequence, structure and position annotations but no assay description, and
"non-natural substrate" in an earlier draft was our inference from the engineering
context, not something we know. We are retrieving the Johnston et al. methods section.

Specifically we will look for: the exact indole analogue and product; whether the
readout is purified enzyme, lysate, whole-cell or growth selection; initial rate versus
endpoint conversion versus enrichment-derived fitness; temperature, pH, substrate
concentrations, PLP supplementation and duration; whether TrpA is present; whether
activity is normalised to enzyme abundance; replicate uncertainty and the meaning of
negative normalised values.

We also accept the point that our claim "â‰¤0.05 is essentially inactive" is not
justified without knowing the detection limit â€” 0.05 is ~12% of parent fitness, which
may or may not be distinguishable from background.

**Highest-priority data-quality action:** all five â‰¥0.90 variants share the `*KG` motif.
If any of them are among the 871 KNN-imputed values, the entire K227/G228 finding
weakens. We are recovering the imputation flags before placing further weight on it.
The five rows to check are listed in Â§5.4.

---

## 9. Known limitations

- Up to 871 of 159,999 fitness values (0.54%) may be KNN-imputed rather than measured.
- A single target, a single four-site library. No claim generalises beyond it. We also
  hold avGFP, GB1 and trastuzumab datasets but have not run this protocol on them.
- The "assay" is a table lookup. Real plate noise, replicate variance, edge-well
  evaporation and 24â€“120 h turnaround are absent, though the cost and authorisation
  are modelled as if performed.
- The hand-built prior (Â§5) is fitted to nothing. We regard this as the weakest
  scientific component of the system.
- `[literature]` claims here are from the cited sources; everything else is measured by
  us and reproducible from the repository.

