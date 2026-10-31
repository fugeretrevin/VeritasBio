# Prior art: depth-adaptive allocation, and the four supporting questions

**Status:** literature and prior-art review. Written 2026-10-31, before the cross-protein
evaluation of `docs/preregistration-cross-protein.md` has been run, and before the §13
intervention has been implemented.

**Scope.** This document answers five questions put to it, in priority order. It does not
re-derive any measurement from this repository; the TrpB/GB1 numbers it refers to are
carried in from `docs/preregistration-cross-protein.md` §6.1 and are treated as established.

**Verification convention.** Every claim about an external paper carries a marker:

- **[full text]** — the relevant section was read.
- **[abstract/snippet]** — only an abstract, a repository README, a press release, or a
  search-engine extract was read. The claim may be right and is not confirmed.
- **[not read]** — cited for identification only. Nothing is asserted about its contents
  beyond its title.

Where no prior art was found, this document says so rather than listing adjacent work as
if it were a hit. Section 7 collects everything that **contradicts or endangers** our
findings or our plan; it is the part worth reading twice.

---

## 0. The verdict, first and unambiguously

> **No prior art found.** There is no published method that estimates per-mutation-depth
> return from a campaign's own measurements and reallocates subsequent batches across
> depth shells. Depth-adaptive allocation, as specified in
> `docs/preregistration-cross-protein.md` §13.1, is novel as far as this review can
> establish.

That verdict is narrow and it needs its boundaries stated, because three of its four
ingredients are individually old:

| ingredient | state of the art | gap |
|---|---|---|
| mutation depth as the thing that matters | **established** (Drummond & Arnold 2005; MULTI-evolve 2026) | treated as a library-design or protocol constant, never as a per-round decision |
| a distance-from-parent knob in the acquisition function | **established** (Zong 2023; ODBO 2022; AntBO 2022) | the knob is a fixed hyperparameter, explicitly left to the user |
| adapting a search radius from the algorithm's own success/failure record | **established** (TuRBO, via ODBO 2022) | the radius lives in an encoded feature space or a hyper-rectangle, **not** in depth-from-wild-type |
| a bandit whose arms are strata, with per-stratum allocation | **established** (MC-UCB, Carpentier & Munos 2011) | applied to Monte-Carlo variance reduction, never to sequence-space depth shells |

So the honest novelty claim is not "nobody has thought about mutation depth". It is:

> **Mutation depth has never been used as an action space.** It appears in this literature
> as a library-construction parameter, as a static regularizer coefficient, as an
> evaluation axis, and as a hand-designed experimental schedule — but no published policy
> treats "which depth shell do I spend the next batch in" as a decision it updates from its
> own returns.

Two consequences follow immediately, and both are uncomfortable:

1. **The strongest competitor to depth-adaptive allocation is not another adaptive method;
   it is a fixed depth schedule.** MULTI-evolve (Science 2026) is a hand-designed one and it
   is the highest-profile protein-engineering result of the past year. Zong's
   `f_λ(s) = f(s) − λ·d(s, s₀)` is a one-line fixed one. This is exactly kill criterion
   §13.5.4, and the literature has already supplied the baselines that would fire it. They
   must be in the arm list from the start, not added after the adaptive version wins.
2. **The methodological ancestor to cite is MC-UCB, not anything in protein engineering.**
   Framing the contribution as "known stratified-bandit machinery, applied to a
   stratification of sequence space that has not been used as an action space before" is
   both more defensible and more likely to survive review than a novelty claim about the
   bandit itself.

---

## 1. Question 1 — has this been done?

Searched under: mutation-depth-aware acquisition; depth-stratified or Hamming-shell batch
allocation; stratified and hierarchical bandits over sequence neighbourhoods;
mutation-distance-constrained acquisition; exploration-radius and trust-region scheduling in
sequence space; epistasis-order-aware acquisition; mutation-order budget allocation.

### 1.1 The closest hits, ranked by how close they actually are

**(a) MULTI-evolve — a fixed, hand-designed depth schedule. Closest in intent.**

Tran V.Q., Nemeth M., Bartie L.J., Chandrasekaran S.S., Fanton A., Moon H.C., Hie B.L.,
Konermann S., Hsu P.D. (2026). *Rapid directed evolution guided by protein language models
and epistatic interactions.* Science. DOI [10.1126/science.aea1820](https://doi.org/10.1126/science.aea1820).
Code: https://github.com/ArcInstitute/MULTI-evolve

**[abstract/snippet — read via the Arc Institute summary, https://www.arcinstitute.org/news/multi-evolve; the Science full text was not accessible to this review]**

What it does: identify ~15–20 single beneficial substitutions (language models or screen),
then measure **all pairwise combinations** of them (~100–200 variants), train on
singles + doubles, and predict combinations up to **depth 7**. As few as 9 predicted
multi-mutants are then tested. One round, against the 5–10 cycles of conventional directed
evolution. Reported: APEX up to 256-fold over wild type, dCasRx 9.8-fold, anti-CD122
antibody 2.7-fold binding.

**How close:** it allocates measurement effort by mutation order — depth 1, then depth 2,
then a jump to high depth. That is a depth schedule and it is the closest thing in the
literature to our intervention's *purpose*.

**How it differs, specifically:** the schedule is **fixed in advance and identical for every
target**. There is no estimate of per-depth return, no reallocation, and no round in which
the method could decide that depth 2 is not paying and move. The depth-7 jump is
extrapolation from a surrogate, not an allocation decision informed by depth-7 returns.

**Why this matters to us more than a citation:** MULTI-evolve is, in effect, the fixed-depth
heuristic of our own kill criterion §13.5.4, published in Science, with wet-lab results. If
depth-adaptive allocation cannot beat a singles-then-pairs-then-extrapolate schedule, the
intervention is dead and the field already has the answer.

**(b) Zong (2023) — a static distance-from-wild-type penalty. Closest in mechanism.**

Zong C. (2023). *Protein Sequence Design with Batch Bayesian Optimisation.*
arXiv:[2303.10429](https://arxiv.org/abs/2303.10429). DOI
[10.48550/arXiv.2303.10429](https://doi.org/10.48550/arXiv.2303.10429)

**[abstract/snippet]** Defines a regularized objective `f_λ(s) = f(s) − λ·d(s, s₀)` where
`d` is Hamming distance to the starting sequence, and restricts search to the wild-type
neighbourhood on the grounds that fitness is often gained through few substitutions. Larger
λ keeps the search shallow; smaller λ explores further.

**How close:** this *is* a depth knob in the acquisition function, aimed at the wild-type
reference, with the right sign.

**How it differs, specifically:** λ ≥ 0 is a **fixed hyperparameter**, and the paper says in
as many words that the right value "depends on the problem and the desired balance". It is
set once, by hand, before the campaign. Nothing updates it. There is no per-depth return
estimate. One λ meeting two different neighbourhood structures is precisely the failure we
measured on TrpB versus GB1 — so this paper documents the knob and hands the problem of
setting it to the user, which is the gap we propose to close.

**(c) ODBO — an adaptive radius, in the wrong space. Closest in machinery.**

Cheng L., Yang Z., Liao B., Hsieh C., Zhang S. (2022). *ODBO: Bayesian Optimization with
Search Space Prescreening for Directed Protein Evolution.*
arXiv:[2205.09548](https://arxiv.org/abs/2205.09548)

**[full text, via ar5iv]** Adapts TuRBO-1 to directed evolution. Appendix C.3: the trust
region's base length `L` is driven by the algorithm's own success and failure counters —
`L` doubles (capped at `L_max`) after `τ_s` consecutive successes and halves after `τ_f`
consecutive failures. Datasets: GB1 4-site (149,361 variants), GB1 55-site, Ube4b, avGFP.
Initial set 40 on GB1(4); q = 1 sequential for the headline results, with q = 5 and q = 10
batch ablations; 50 iterations. Metrics: average maximum fitness, and **top 1%, 2% and 5%
hit ratios**. XGBOD prescreening is scored with accuracy/precision/recall/FNR/FPR.

**How close:** this is an exploration radius that the algorithm adapts from its own
observations, with no labels and no metadata — the same control loop we propose.

**How it differs, specifically:** Section 3.5 puts the trust region in the **encoded feature
space**, as a hyper-rectangle centred on the current best solution. The paper never mentions
mutational distance. So the radius is not depth-from-wild-type, the reference point moves
with the incumbent rather than being the library's wild-type shell structure, and the
quantity adapted is a region size rather than an allocation across shells. Initial sampling
balances residue coverage per site but **does not allocate by mutation count**.

**Secondary value to us:** ODBO's top-1%/2%/5% hit-ratio reporting is a direct precedent for
our raw top-tier endpoint (Question 4), and its GB1 4-site universe of 149,361 is an
independent corroboration of our 149,360 (the one-row gap is pre-registered as discrepancy
D-3).

**(d) Trust regions defined by sequence distance — the mechanism exists, nobody adapts it
per depth.**

- **AntBO**: Khan A. et al. (2022), *Toward real-world automated antibody design with
  combinatorial Bayesian optimization*, Cell Reports Methods.
  https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(22)00276-4
  **[abstract/snippet]** Builds a trust region in combinatorial sequence space for CDRH3
  design, sized as the distance of the best-seen sequence from the random start, shaped by
  biophysical constraints. Distance-based and adaptive in size, but the arm is still a
  region, not a depth, and the reference is the start sequence.
- **PepCompass**: arXiv:[2510.01988](https://arxiv.org/abs/2510.01988)
  **[abstract/snippet]** Trust region = sequences within a fixed edit distance of the current
  best, with BO inside it. The distance is **fixed**.

**(e) Tree Search-Based Evolutionary Bandits — a bandit, and a depth, that never meet.**

Qiu J., Yuan H., Zhang J., Chen W., Wang H., Wang M. (2024). *Tree Search-Based Evolutionary
Bandits for Protein Sequence Optimization.* AAAI-24.
arXiv:[2401.06173](https://arxiv.org/abs/2401.06173);
proceedings: https://ojs.aaai.org/index.php/AAAI/article/view/29386

**[abstract + third-party summary; the AAAI full text was behind bot protection and the
arXiv HTML render does not exist, so the handling of N is *not* confirmed from the paper
itself]** Nodes are sequences, edges are mutations or recombinations. The depth `N` is
described as the maximum allowed mutations from wild type and appears as a **parameter in
the Bayesian regret bound**. The bandit (TreeNeuralUCB / TreeNeuralTS) scores candidate
sequences by predicted fitness plus uncertainty. Batch parameters are `n₁` mutants, `n₂`
recombinants, `n₃` retained. Datasets: AAV (~20²⁸), TEM (~20²⁸⁶), AAYL49 (~20¹¹⁸). Reported:
top AAV sequences typically found within 3–5 mutations of wild type.

**How close:** it is the only protein-sequence work found that puts a bandit and a
mutation-depth bound in the same algorithm.

**How it differs, specifically:** **the arms are sequences, not depths.** `N` is a bound on
the search tree, not an arm, and no per-depth allocation appears anywhere in the described
pseudo-code (`n₁`/`n₂`/`n₃` split mutants from recombinants, not shallow from deep). The
"3–5 mutations" observation is a reported outcome, not a controlled variable.

**Caveat on this entry:** because the full text was not readable, there remains a residual
risk that the paper schedules `N` across rounds in a way the abstract does not advertise.
**Recommended action: obtain the AAAI PDF and check how `N` is set before any novelty claim
is published.** This is the single largest unresolved risk to the §0 verdict.

**(f) FolDE — diversity-aware batching, not depth-aware.**

*Low-N Protein Activity Optimization with FolDE*, PMC12636759.
**[abstract/snippet]** Adds a diversity-aware batch selector (constant-liar). The authors
report it "improves batch diversity; this is important in multi-mutation campaigns but had
limited effect in our benchmarks". Diversity in the batch, not allocation across depths.

### 1.2 The methodological ancestor, outside protein engineering

Carpentier A., Munos R. (2011). *Finite-Time Analysis of Stratified Sampling for Monte
Carlo.* NeurIPS 2011.
https://papers.neurips.cc/paper_files/paper/2011/file/872488f88d1b2db54d55bc8bba2fad1b-Paper.pdf
Journal version: Carpentier, Munos, Antos, JMLR 16 (2015), carpentier15a.
Hierarchical extension: MC-ULCB, allocating samples over a hierarchical partition.

**[abstract/snippet]** Stratified Monte-Carlo integration framed as a bandit in which
**each stratum is an arm**; MC-UCB allocates samples in proportion to an upper confidence
bound on each stratum's standard deviation, with Bernstein-based confidence intervals and a
finite-time pseudo-regret analysis (Õ(K n^{−4/3}), later tightened in K, with a matching
lower bound).

This is the correct citation for the *estimator*. Note the objective differs: MC-UCB
allocates to reduce integration variance, whereas we allocate to maximize per-depth hit
return, which makes our problem a plain stratified-UCB / best-arm-allocation problem rather
than MC-UCB itself. Cite it as the lineage, not as the method.

### 1.3 What was searched and not found

Stated as clean negatives:

- **No prior art found** for a policy that estimates per-depth (per-Hamming-shell) return
  from its own measurements and reallocates batch slots across depths.
- **No prior art found** for a stratified or hierarchical bandit whose arms are
  mutation-depth shells of a protein library.
- **No prior art found** for scheduling an exploration radius in sequence space across
  rounds from observed returns (as opposed to TuRBO-style success/failure counters in
  feature space, which does exist — see 1.1(c)).
- **No prior art found** for epistasis-order-aware acquisition in the sense of allocating
  budget across mutation orders by measured return. Order-aware **modelling** exists
  (e.g. a transformer that caps epistasis order by attention depth, eLife reviewed preprint
  https://elifesciences.org/reviewed-preprints/108005 **[not read]**), and order-stratified
  **evaluation** exists (CombinGym's k-vs-rest, §3), but not order-aware allocation.

---

## 2. Question 2 — the three papers already in our plan

### 2.1 ALDE — Yang et al., Nature Communications 2025

Yang J., Lal R.G., Bowden J.C., Astudillo R., Hameedi M.A., Kaur S., Hill M., Yue Y.,
Arnold F.H. (2025). *Active learning-assisted directed evolution.* Nat Commun 16:714.
DOI [10.1038/s41467-025-55987-8](https://doi.org/10.1038/s41467-025-55987-8).
Full text: PMC11739421.

**[full text]**

**The published protocol, exactly.**

| setting | ALDE as published | ours |
|---|---|---|
| simulated campaigns (GB1, TrpB, 4 sites) | random batch of **96**, then **4 rounds of 96** = **480** measurements; averaged over **70** random initializations | **25**, then **4 rounds of 25** = **125**; 48 seeds |
| wet lab (ParPgb, 5 sites) | 384 variants screened for the initial library, **216** usable as training data; rounds 1 and 2 each tested the top **90**; three rounds total; "six 96-well plates" | n/a |
| landscapes | GB1 and TrpB only, both 20⁴ = 160,000 | same two, plus the CombinGym corpus |
| models | boosting ensemble, GP, DNN ensemble, DKL | — |
| encodings | AAIndex, Georgiev, one-hot, ESM2 | — |
| acquisitions | greedy, UCB (β = 4), Thompson sampling. Batch expected improvement was tried and dropped as slow without clear gain | — |
| primary metric | **maximum fitness at end of campaign, normalized to the best variant in the design space** | diversity-adjusted yield at d ≥ 2 |
| top-1% | **not reported as a metric** | our secondary |
| diversity | discussed only qualitatively ("expanded mutational diversity" in round 2); **no quantitative measure** | quantified |
| calibration | MAE, miscalibration area, Spearman(uncertainty, error) | — |

**How much the adaptation could matter: a lot, and in the direction that flatters us.**

1. **The round structure is identical; the data volume is not.** Our 125 measurements are
   **26%** of ALDE's 480. Round counts match (4), so the comparison is not "fewer rounds" —
   it is **the same number of model refits, each on 3.84× less data.**
2. **The regime difference lands squarely on ALDE's model class.** ALDE's published
   configuration space is DNN ensembles and deep kernel learning over ESM2 embeddings, and
   its wet-lab campaign used a **DNN ensemble with 216 training points**. At 25 initial
   observations over a 160,000-variant space, a DNN ensemble and a DKL surrogate are being
   run at roughly a tenth of the smallest training-set size the paper ever demonstrates.
   One-hot ridge degrades far more gracefully there. **A deep-surrogate baseline evaluated
   below its demonstrated data regime will look worse than it is, and our engine's
   advantage over it will look larger than it is.**
3. **The paper reports no budget or batch-size sensitivity analysis.** **[full text]** It
   shows full per-iteration optimization trajectories (Fig. S73) and states that how samples
   and rounds should scale with the number of simultaneously explored sites is an **open
   question for future work**. Fig. S73 is a partial substitute: a trajectory lets one read
   off performance at a smaller **cumulative** budget (e.g. after 1 round = 192
   measurements). It does **not** answer the question we need, which is performance at a
   smaller **batch size** — 25 rather than 96 — because the trajectory holds batch size
   fixed at 96 throughout. Nothing in the paper licenses extrapolation to a 25-wide batch.
4. **Direction of the one available external datum is against the simple reading.** A 2026
   preprint on far-red fluorescent protein evolution reports ProteinGym simulations in which
   "more experimental batches with fewer samples per batch consistently accelerated
   convergence" (bioRxiv [10.64898/2026.08.12.744534](https://doi.org/10.64898/2026.08.12.744534),
   **[abstract/snippet]**, not peer reviewed). That is a claim about **reallocating a fixed
   total** across more, smaller rounds. Our adaptation does the opposite thing — it holds
   rounds fixed and **cuts the total** — so this preprint does not rescue the comparison.

**What this obliges us to do.** `docs/preregistration-cross-protein.md` §1.1 already requires
the deviation sentence on every external claim, and §7.2 already forbids silent substitution
of a different ALDE configuration. Both are correct and sufficient on process. The
substantive addition this review recommends:

> The development headline "best arm 50.44 top-1% hits against ALDE's 30.42" is measured at
> **26% of the budget at which ALDE was demonstrated**, with ALDE's deep surrogates far
> below their demonstrated training-set size. It should be reported with ridge as the
> reference baseline (as §7 already designates) and **should not be the sentence that
> carries the claim**. If an ALDE comparison is to carry weight, the honest version is to
> also run the arm grid at ALDE's own 96 + 4×96 schedule on TrpB and GB1 — which are
> development data, so this costs nothing in held-out evidence — and report both.

### 2.2 Greenman, Amini & Yang — uncertainty quantification benchmark

Greenman K.P., Amini A.P., Yang K.K. (2025). *Benchmarking uncertainty quantification for
protein engineering.* PLOS Comput Biol 21(1):e1012639. DOI
[10.1371/journal.pcbi.1012639](https://doi.org/10.1371/journal.pcbi.1012639).
Full text: PMC11741572. Earlier versions: bioRxiv 2023.04.17.536962; ICLR 2022 workshop.

**[full text]**

**Setup.** Seven UQ methods — Bayesian ridge, GP, and five CNN variants (ensemble, dropout,
SVI, MVE, evidential). Eight of the 15 FLIP tasks: AAV/Random, Meltome/Random, GB1/Random,
AAV/Random-vs-Designed, GB1/1-vs-Rest, AAV/7-vs-Rest, GB1/2-vs-Rest, GB1/3-vs-Rest. One-hot
and protein-language-model embeddings. Metrics: RMSE, MAE, R², Spearman ρ; miscalibration
area (AUCE), coverage (% within ±2σ, ideal 95%), width (4σ normalized by training range),
and ρ_unc = Spearman(uncertainty, |error|). Retrospective active learning: start from a
random 1%, 5% or 10% of the training data; **five acquisition rounds with training-set sizes
equally spaced on a log scale**; three seeds; strategies are "explorative greedy" (highest
uncertainty), "explorative sample" (uncertainty-weighted), and random; headline metric is
Spearman ρ.

**Does it support our proposed mechanism — that poor calibration produces the shallow-
exploitation signature we observe? No. Stated plainly: it does not, and we must not cite it
as if it did.**

- **[full text]** The paper **does not link miscalibration to over-exploitation or
  under-exploration.** It reports that calibration quality varies with domain shift in a
  model-dependent and non-monotone way: some models are well calibrated on the hardest
  splits while others are poorly calibrated even on random splits.
- **[full text]** It **does not stratify any result by mutation count or distance from wild
  type.** Domain shift is summarized only at the split level. The GB1 k-vs-Rest splits are
  described as "designed primarily to test generalization from few- to many-mutation
  sequences" — which is the right axis — but no performance-versus-depth curve is reported.
  So the paper establishes that the few-to-many-mutation axis is a recognized hard
  generalization axis, and nothing more than that.
- **[full text]** The one statement in the direction we want is weak and is the authors'
  own suggestion, not a demonstration: when uncertainties are poorly calibrated early on,
  **some randomness may help**. CNN SVI is noted to show low coverage and low width, which
  they read as possible over/under-confidence; MVE degrades under high domain shift. None of
  this is tied to selection depth.

**The finding that actually matters to us, and it cuts against our engine.**
**[full text]** Across their benchmark, uncertainty-based acquisition **almost always beat
random but never beat greedy**. UCB ≈ greedy; Thompson sampling sat between greedy/UCB and
random; evidential UQ on AAV was **worse than random**. The authors offer two readings:
better UQ methods are needed, or these landscapes are simple enough that pure exploitation
is optimal.

That is an independent, peer-reviewed result saying that on protein fitness landscapes an
uncertainty-driven acquisition should not be expected to beat a greedy baseline. Our GB1
result — ridge-greedy beating the reliability policy by 22.49 raw top-1% hits — is
**consistent with their benchmark**, and our TrpB result is the exception that needs
explaining, not the rule. See §7.2.

**Usable conclusion.** Cite Greenman et al. for (a) the greedy-is-hard-to-beat baseline
result, which strengthens the case for ridge-greedy as our reference arm, and (b) the
few-to-many-mutation split as a recognized generalization axis. **Do not cite it for the
calibration→shallow-exploitation mechanism.** If we want that mechanism, we have to
demonstrate it ourselves, from the per-round per-depth rows that
`docs/preregistration-cross-protein.md` §10 already requires — specifically by showing that
predicted-uncertainty rank and realized-error rank decouple as depth rises. No external
paper will carry that argument for us.

### 2.3 CombinGym — Chen et al., bioRxiv 2026

Chen Y., Fu L., Lu X., Li W., Gao Y., Wang Y., Ruan Z., Si T. (2026). *CombinGym: a
benchmark platform for machine learning-assisted design of combinatorial protein variants.*
bioRxiv, posted 2026-03-25. DOI path
[10.64898/2026.03.24.714074](https://doi.org/10.64898/2026.03.24.714074).
Repository: https://github.com/sitonglab/CombinGym Website: https://www.combingym.org/
**Not peer reviewed.**

**[abstract via ScienceCast https://bio.sciencecast.org/casts/d4fvyphgo79k; repository README
and GitHub commit API; the bioRxiv full text was not reachable by this review]**

**Landscapes and proteins — the answer, and the reconciliation.**

- Abstract: **14 curated combinatorial mutagenesis datasets spanning 9 proteins**, covering
  binding affinity, fluorescence and enzymatic activity.
- README: 14 datasets, **9 proteins**, **4 to 16 mutated positions**, **>400,000 data
  points**.
- This matches `docs/preregistration-cross-protein.md` §3.2 exactly, including its
  resolution that landscape ≠ protein (three proteins carry multiple phenotypes: 3 + 3 + 2 +
  6×1 = 14 over 9). **That reconciliation is confirmed independently here and needs no
  revision.**

**The version to pin.** CombinGym has **no releases and no tags** (Releases and Packages
sections empty), so a commit is the only pin available. Confirmed from the GitHub commits
API:

| field | value |
|---|---|
| tip of `main` | `c5dccd514a9acfe1089a5f1555faea98afd78019` |
| date (UTC) | 2026-09-16T07:33:53Z |
| message | "Add files via upload" |
| previous commits | `ebbd2c43…` (2026-09-16, deletes `Data/Structure/SaCas9(SaCas9-KKH).pdb`), then a gap back to `e635a959…` / `0851af0c…` / `a0814441…` (all 2025-03-19) |
| total | 18 commits, 0 tags, 0 releases |
| license | MIT |
| citation | README says "currently under submission"; paper link "to be added upon publication" |

**This is the same commit `docs/preregistration-cross-protein.md` §3.1 already pins.
Confirmed correct.** Two additional notes for the freeze:

- The 2026-09-16 pair of commits is the **first change in eighteen months** (previous
  activity 2025-03-19). The corpus is effectively static but has demonstrated that it
  changes without a version bump — and one of those two commits **deleted a data file**.
  A commit pin is therefore not merely the best available option, it is load-bearing.
- No Zenodo DOI and no versioned data release exist. The hash table in §3.4 of the
  pre-registration is the only defence against a silent upstream edit. Keep it.

**Baselines it ships — and they are not what we need.** README lists nine, which map to the
abstract's five methodological categories (alignment-based, protein language, structure-based,
sequence-label, substitution-based):

**EVmutation, DeepSequence, ESM-1b, ESM-1v, GVP-Mut, MAVE-NN, CNN, Ridge, BLOSUM62.**

Metrics: **Spearman ρ and NDCG**, in zero-shot and supervised settings. Benchmark artifacts
at the pin: `Benchmarks/Leaderboard_{NDCG,Spearman}.xlsx` and
`Benchmarks/{NDCG,Spearman}/{0,1,2,3}-vs-rest_*.xlsx`.

> **The critical limitation for us: every one of those nine is a *fitness predictor*, scored
> on correlation and ranking. CombinGym ships no active-learning loop, no sequential design
> protocol and no batch-selection utility.** This is confirmed from both the README and the
> abstract. There is nothing to inherit: our arm list (random / ridge-greedy / RF-greedy /
> ALDE / three engine configurations) is the campaign baseline set, and CombinGym supplies
> only the landscapes. `docs/preregistration-cross-protein.md` §7 is therefore correctly
> specified, but the paper cannot be cited as supplying campaign baselines, and a reviewer
> who assumes otherwise should be corrected. One genuine overlap: CombinGym's `Ridge` is a
> predictor on the same footing as our A2 surrogate, which is a small point in favour of
> ridge as the reference baseline.

**The k-vs-rest splits, and the one place CombinGym bears directly on H1.** 0-vs-rest through
3-vs-rest; the stated conclusion is that **lower-order mutation data is valuable for
predicting higher-order mutant properties**. Validated in silico (brightness of an
oxygen-independent fluorescent protein) and experimentally (enzyme substrate specificity,
"substantial increase in specific activity"). **[abstract]**

The abstract does **not** say whether this varies across proteins or assays — but the
per-split leaderboards at the pin (`{0,1,2,3}-vs-rest_Spearman.xlsx`, 14 keys each) **contain
exactly that breakdown**. This is an opportunity, and §4.3 of this document picks it up.

---

## 3. Question 3 — does near-wild-type shell enrichment vary predictably by protein family,
assay type or library design?

### 3.1 The answer

> **No prior art found.** There is no established literature establishing that near-wild-type
> shell enrichment varies systematically and *predictably* across protein families or assay
> types. No cross-protein meta-analysis of hit-mass-versus-mutation-depth, stratified by
> protein class or assay type, was located. Three independent searches were run on this
> question and all three came back empty on the specific claim.

This is a clean negative, and it is the commercially interesting one. **It means the
prospective version of our covariate — "predict shallow-shell richness from protein class,
assay type or library design, then choose targets" — is unclaimed.** It also means we cannot
borrow a prior for it; it would have to be built.

### 3.2 What does exist, and exactly how far short it falls

**(a) The functional fraction decays with mutation count. Established, aggregate, not
stratified.**

Multiple studies report that the probability a protein remains functional after *m*
mutations falls off **exponentially** in *m*, with a review-level model of W ≈ e^(−0.36n),
sometimes steeper because of epistasis, and with documented deviation over the first few
mutations when the starting protein's stability or activity is unusually high or low.
**[abstract/snippet]** Sources: "Stability effects of mutations and protein evolvability",
Curr Opin Struct Biol (https://www.sciencedirect.com/science/article/abs/pii/S0959440X09001249);
Bloom et al., *Evolution favors protein mutational robustness in sufficiently large
populations*, BMC Biol 5:29 (arXiv:[0704.1885](https://arxiv.org/abs/0704.1885)).

**Why this is not our covariate.** It describes the **functional fraction** (does the protein
still work), aggregated over proteins. Ours is **top-tier hit density within a depth shell**
(is this shell enriched in the best variants) for a **specific 4–16 site library**. These can
point in opposite directions: a shell can be rich in functional variants and poor in
improved ones, which is the ordinary case for a well-adapted wild type. Nothing in this
literature stratifies by family or assay.

**(b) Drummond & Arnold 2005 — the one source that says the optimum is protein-dependent.
Closest thing to a prospective rule, and it is about library construction.**

Drummond D.A., Iverson B.L., Georgiou G., Arnold F.H. (2005). *Why high-error-rate random
mutagenesis libraries are enriched in functional and improved proteins.* J Mol Biol
350:806–816. PMID [15939434](https://pubmed.ncbi.nlm.nih.gov/15939434/).
Preprint: arXiv:[q-bio/0411041](https://arxiv.org/abs/q-bio/0411041).
PDF: https://cheme.caltech.edu/groups/fha/publications/drummondJMB2005.pdf

**[abstract/snippet]** An optimal mutation rate exists, balancing uniqueness (lost if the
rate is too low) against retained function (lost if too high). Error-prone PCR's non-Poisson
mutation distribution explains the excess of functional clones at high error rates. For
their highest-rate scFv parameters the predicted optimum is ~3.1 nucleotide mutations per
gene, yielding >35,000 unique functional sequences. **Crucially: optimal rates depend on both
the protein and the mutagenesis protocol, and the more tolerant the protein, the higher the
optimal rate.**

**How close:** this is the field's canonical statement that **the right mutation depth is a
property of the target**, with a model and an experiment behind it. It is the right ancestor
for the *idea* that depth is target-dependent, and it should be the first citation in any
write-up of our covariate.

**How it falls short of what we need:** (i) it is about **library construction** — the
distribution of mutations you build in — not about allocating measurements across depths
during a campaign; (ii) the optimum is worked out for **one protein and one protocol**, with
the explicit warning not to transfer the number; (iii) "tolerance" is not operationalized
into anything measurable before a campaign, so there is **no prospective rule** here. It
establishes the dependence and leaves the prediction problem open.

**(c) Family-wide comparison is explicitly described as scarce.**

A B1 β-lactamase study states that "the experimental study of functional constraints at a
family-wide level by comparing residue-level mutational epistasis across homologous
sequences remains scarce" **[abstract/snippet]**
(bioRxiv 2023.11.27.568921). The nearest multi-protein comparison located is a 2025 preprint
over **three** proteins (TEM-1, YAP1 WW, Pab1 RRM2) reporting a clear but nonlinear
single/double-mutant fitness relationship across all three (bioRxiv 2025.09.14.676130
**[abstract/snippet]**). Three proteins is not a meta-analysis.

**(d) Assay-class intuitions exist as qualitative claims, not measured rates.**

- Active-site-focused libraries are described in review as having "the highest success rate"
  for altering substrate or product specificity, especially iteratively; one source argues
  this is because random methods mostly hit residues far from the active site.
  **[abstract/snippet]** — both are **qualitative judgments**, not measured
  improved-variant fractions.
- GFP chromophore-pocket work reports that most active-site mutations destabilize or
  disable, and functional multipoint mutants are "exceptionally rare", while a designed
  library nonetheless recovered >16,000 functional designs with up to eight active-site
  mutations (PMC10199939 **[abstract/snippet]**). This shows **library design** moves shell
  structure a great deal — which supports the covariate being real — but gives no
  cross-family rate.
- No source was found reporting improved-variant rates for protein–protein binding
  interfaces on a footing comparable to enzyme active sites. **The enzyme-versus-interface
  comparison our H1 would most like to lean on does not exist in the literature.**

### 3.3 The structure already in our own corpus, and the one external paper to check

Our pre-registered corpus statistics (`docs/preregistration-cross-protein.md` §6.5) contain a
pattern that, if it held up, would be the beginning of the missing literature:

- **Antibody affinity-maturation landscapes are deep-shell enriched, monotonically.** CR6261
  H9: 0% hit rate up to depth 5, 22.9% at depth 9, 40.0% at depth 10. CR9114 H1 likewise
  rises with depth. These are structurally the inverse of TrpB.
- **Enzyme and photosensor landscapes are shallow-favouring.** CreiLOV falls from 9.66% at
  depth 2 to 0.12% at depth 12; RhlA's hits sit at mean depth 1.89–2.89.
- **The same protein can go either way depending on the readout.** eqFP611 blue falls
  monotonically with depth (30.8% at depth 1); eqFP611 **red** has *zero* hits at depths 0–3
  and peaks at depth 8 — same protein, same 8,192 genotypes, opposite shell structure.

That last observation is the important one and it **cuts against the commercial framing**:
if one library's two readouts have opposite depth structure, then **shallow-shell richness
cannot be a property of the protein or its family alone.** Assay and readout matter at least
as much, and possibly more. A target-selection tool keyed on protein class would mis-call
eqFP611 red.

**One external paper should be read before anything is claimed here.** The CR9114/CR6261
landscapes originate in:

Phillips A.M. et al. *Binding affinity landscapes constrain the evolution of broadly
neutralizing anti-influenza antibodies.* eLife 10:e71393.
DOI [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393).
https://elifesciences.org/articles/71393

**[not read — the full-text fetch was refused by the summarization model used in this
review, twice, so nothing is asserted about its contents.]** Its title indicates it analyses
how affinity is **constrained** across its combinatorial library, which is exactly the
depth-structure question for those two landscapes. If it already reports that high affinity
is reachable only at high substitution counts, then the deep-shell structure of our two
antibody targets is **published, known in advance, and not our finding** — and our H1 ranking
of those targets is partly anticipated. **This must be read before the H1 result is written
up.** It is listed as an open item in §8.

### 3.4 The prospective version, and whether it is reachable

Our covariate is currently retrospective: it is computed from a **fully measured** landscape
(`docs/preregistration-cross-protein.md` §6.2 is explicit and correct that this is what
makes it pre-registerable). Converting it to a prospective target-selection tool requires
predicting shell enrichment **before** the landscape is measured. Three candidate routes,
with honest assessments:

1. **From protein class / assay type.** No literature exists (§3.1) and our own eqFP611
   blue-versus-red split argues the signal is not primarily in protein class (§3.3).
   **Weak.** Do not build a product thesis on this without first establishing the
   cross-family relationship that §3.1 says is missing.
2. **From library design.** Better grounded: Drummond & Arnold establishes
   protocol-dependence, and the GFP designed-library result shows design moves shell
   structure a lot. Also directly actionable, since library design is a parameter we would
   control. **Most promising, and least explored.**
3. **From CombinGym's own k-vs-rest leaderboards.** See §4.3 — this is the cheapest
   available test and it uses already-published numbers.

And one route that is **already ruled out, correctly, by our own pre-registration**: §13.4
shows an upfront probe pre-screen is hopeless (expected hit counts of 0.6 versus 0.26 at 8
probes for TrpB versus GB1). That arithmetic is sound and this review found nothing that
contradicts it.

---

## 4. Question 4 — what endpoints does this field use?

### 4.1 What the field actually reports

| source | primary | diversity / separation endpoint | verification |
|---|---|---|---|
| **ALDE** (Nat Commun 2025) | max fitness at campaign end, normalized to the design space's best | **none quantitative**; diversity discussed qualitatively only. Top-1% not reported | [full text] |
| **ODBO** (2022) | average max fitness over 10 runs | none; but reports **top 1% / 2% / 5% hit ratios** | [full text] |
| **FLEXS / AdaLead** (2020) | max fitness `y = max φ(x)` | **yes** — `\|S\|`, the count of sequences above a threshold `y_τ`; and, where the landscape is enumerable, **`\|M'_{y_τ}\|`, the number of distinct local maxima found above `y_τ``. Explicitly "diversity in a narrow sense". No top-K rule, no Hamming-separation rule. Fitness and diversity reported **separately**; no combined endpoint | [full text, ar5iv] |
| **Jain et al.** GFlowNets (ICML 2022) | mean oracle score over top-K (K = 100 AMP, K = 128 TFBind8/GFP) | **yes** — diversity = mean pairwise distance within the top-K set; novelty = mean distance from each top-K candidate to its nearest neighbour in the initial dataset. Distance measure left unnamed | [full text, ar5iv] |
| **ProSpero** (2025) | max and mean fitness of top 100 | **yes** — diversity = mean pairwise **Hamming** distance among top 100; novelty = mean Hamming distance from top 100 to the start sequence; plus validity (biophysical percentile band) and structure (pLDDT, pTM, scPerplexity) | [full text, arXiv HTML] |
| **Kirjner et al. / Jain lineage** | max fitness; top-k (k = 128) | diversity = **median** pairwise distance; novelty = median of each candidate's min distance to the training set; **Levenshtein**, because sequences can differ in length | [abstract/snippet] |
| **CombinGym** (2026) | Spearman ρ | **none** — plus NDCG for the mutant-design benchmark | [abstract + README] |
| **Greenman et al.** (2025) | Spearman ρ (AL); RMSE/MAE/R² | none | [full text] |
| antibody AL preprint | **Recall@2%** | none | [abstract/snippet] |

### 4.2 The verdict on our choice

**There is precedent for diversity-adjusted and portfolio-style endpoints, but no standard,
and nothing exactly matching ours.** Specifically:

- **Our endpoint is not unprecedented.** FLEXS's `|M'_{y_τ}|` — distinct local maxima above a
  fitness threshold — is a **threshold-plus-separation count**, which is the same shape as
  diversity-adjusted yield at d ≥ 2: it counts *separated* hits rather than raw hits, for
  the same reason. Cite it. ODBO's top-1% hit ratio is the precedent for our mandatory
  secondary. Together these cover both of our endpoints.
- **Our separation rule is stricter and better specified than any of them.** The
  Jain/ProSpero/Kirjner family reports **mean or median pairwise distance over a top-K set**,
  which is a *descriptive summary of a fixed-size set*, not a count. It has a known failure
  mode: a method can raise mean pairwise distance by adding one far-flung candidate while
  the rest stay clustered. A greedy minimum-separation walk (ours) cannot be gamed that way.
  The sequential-filter definition in `docs/preregistration-cross-protein.md` §2 — sort by
  fitness, keep a variant only if Hamming ≥ 2 from everything already kept, applied strictly
  **after** measurement — is a cleaner object than the field's default and should be
  presented as such rather than apologized for.
- **Reporting fitness and diversity *separately* is the field's dominant convention**
  (FLEXS, Jain, ProSpero all do this). Our design — one primary diversity-adjusted endpoint
  plus a **mandatory** raw secondary reported even when it contradicts the primary
  (§2, §9.2) — is more disciplined than the convention and is defensible. But see the
  reviewer expectation below.

**What reviewers will expect, and where we will take fire:**

1. **Max fitness achieved.** ALDE, ODBO, FLEXS, Jain and ProSpero all report it; it is the
   one endpoint universal across this literature. `docs/preregistration-cross-protein.md`
   §9.2.2 already requires "best-so-far fitness against budget (the full curve)". **Good —
   this is the single most important reviewer expectation and it is already covered.** Do not
   let it be dropped to an appendix.
2. **A novelty metric, distinct from diversity.** Jain, ProSpero and the Kirjner lineage all
   report distance-to-training-set alongside within-set diversity. We report redundancy and
   cluster count (§9.2.9) but **no explicit novelty endpoint**. A reviewer from the
   generative-design community will ask. Mean depth of selected variants — which §10 already
   captures per round — is a natural answer, and it happens to be the mechanism variable for
   H1. Consider naming it as a reported secondary rather than leaving it as diagnostic
   instrumentation.
3. **Justifying a non-standard primary.** Changing the primary endpoint away from max fitness
   invites the charge of endpoint-shopping. Our defences are strong and should be stated
   together and up front: the choice was made by the project owner **before the corpus was
   characterised** (§2); the raw count is mandatory and reported even when it contradicts
   (§2); there is **no post-hoc endpoint promotion** (§11.3); and the development data
   already contains a case where the two disagree in sign (frontier beats ridge on raw count
   and is significantly worse on distinct hits). That last fact is the real argument — an
   empirical demonstration that the raw count can be won by collecting near-duplicates — and
   it should lead, not follow.
4. **Comparability of d ≥ 2 across library geometries.** §9.3 of the pre-registration
   already identifies this (the filter keeps 9.8% of GB1's hits but 57% of CR9114's) and
   already fixes it by aggregating ceiling shares. **This is handled, and it is the kind of
   thing a reviewer finds and a benchmark paper usually misses.** No external precedent for
   ceiling-share aggregation was located, so expect to have to explain it from first
   principles.

**No prior art found** for a single scalar that *combines* fitness and diversity into one
portfolio-style endpoint in protein campaign benchmarking. Every source located reports them
separately. One paper on smoothed fitness landscapes explicitly warns that higher diversity
and novelty are **not** equivalent to better performance and should be read descriptively
(arXiv:[2307.00494](https://arxiv.org/abs/2307.00494) **[abstract/snippet]**) — which is an
argument for our structure (one primary count, one mandatory raw secondary) over a weighted
composite.

### 4.3 A cheap, concrete opportunity

CombinGym's pinned `Benchmarks/{Spearman,NDCG}/{0,1,2,3}-vs-rest_*.xlsx` give, for each of
the 14 landscapes, how well each of nine predictors transfers from low-order to high-order
mutants. That is a **published, per-landscape, prediction-transfer analogue of our shell
covariate**, computed by someone else, on the same corpus, at the same pin.

Correlating CombinGym's per-landscape k-vs-rest transfer scores against our shallow-shell
enrichment is a few hours of work on already-frozen numbers, requires no campaign, and would
either (a) give a second, independently published covariate that correlates with ours —
materially strengthening H1's external validity at n = 14 rather than n = 4 — or (b) show
they are unrelated, which is itself informative about what the shell covariate measures.

**Caution, and it is not minor:** this is an analysis over the independent targets' label
structure. `docs/preregistration-cross-protein.md` §4.4 reclassifies a final-test protein as
development data if inspection is used to modify anything. Reading published Spearman values
is not modifying the algorithm, and §6.2 already establishes that label-only landscape
statistics are published for every target precisely so that computing them is not a hidden
look — so this is probably admissible on the same footing as §9.3's ceilings. **But it should
be declared as an amendment in §16 before it is run, not after**, and it must not be allowed
to influence any arm, prior or configuration. If there is any doubt, run it on the reserve
level first.

---

## 5. Question 5 — detecting candidate-generator collapse

### 5.1 The answer

> **No prior art found** for a standardized coverage-based diagnostic for candidate-generator
> collapse in generative sequence models used for protein library design. No paper was
> located that defines distinct-variant coverage, per-round unique fraction, or
> fraction-of-landscape-explored as a **benchmark metric** for this failure mode.

The pieces exist in four separate places, none of them assembled:

### 5.2 What exists

**(a) The failure has been observed and counted in a real MLDE campaign. Closest hit.**

Wu Z., Kaltenbach H.-M., Arnold F.H. et al. (2019). *Machine learning-assisted directed
protein evolution with combinatorial libraries.* PNAS 116(18).
DOI [10.1073/pnas.1901979116](https://doi.org/10.1073/pnas.1901979116).
Preprint: arXiv:[1902.07231](https://arxiv.org/abs/1902.07231)

**[abstract/snippet]** The authors sequenced a model-proposed set **because the predicted
libraries could contain repeats**, and report that only **47 and 39 of 90** proposals were
unique in two position sets. The supplement also records the usual ~95% library-coverage
target, commonly approximated by ~3× oversampling.

**How close:** this is generator redundancy, measured, in exactly our setting
(combinatorial-library MLDE), from the same lab as ALDE. 47/90 and 39/90 are unique
fractions of 52% and 43%.

**How it differs:** it is an **observation inside one campaign**, not a metric, not a
diagnostic, and not proposed for reuse. Nobody defined a threshold or tracked it across
rounds. **This is the single best citation available for "generator collapse is a real,
documented failure mode in this exact setting" — and the absence of any follow-up metric is
the gap our coverage module fills.**

**(b) Mode counting in the GFlowNet literature — the nearest methodological relative, and
the attribution is commonly wrong.**

Jain M. et al. (2022). *Biological Sequence Design with GFlowNets.* ICML, PMLR 162.
arXiv:[2203.04115](https://arxiv.org/abs/2203.04115)

**[full text, ar5iv]** Important correction to a widespread belief: **Jain et al. 2022 does
not report a "number of modes found" metric**, and gives no Hamming threshold for declaring
two sequences distinct modes. "Modes" appear only qualitatively, as modes of the oracle. Its
actual metrics are performance (mean oracle score over top-K), diversity (mean pairwise
distance within top-K) and novelty (distance to nearest neighbour in `D₀`). It
**does not analyse or detect mode collapse in its own generator**; it cites Bengio et al.
2021 that RL generators "can cause candidate generation to get stuck in a single mode" as
*motivation*, and uses a uniform exploration mixture δ to prevent "getting stuck in one or a
few modes". Protocol: AMP N = 10 rounds at b = 1000; TFBind8 and GFP N = 1 round at b = 128.
It does report that spurious proxy modes led to poor de novo GFP generation.

Mode **counting** appears later, in LS-GFN (arXiv:[2310.02710](https://arxiv.org/abs/2310.02710)
**[abstract/snippet]**), which reports >5,000 modes found versus <1,000 for standard GFlowNet
methods on sEH. That is where to look for how mode counts are defined and reported.

**Direct relevance:** our engine uses a GFlowNet for candidate enrichment, and the
configuration whose generator collapsed to 294 distinct variants (0.18% of the landscape)
across 96 campaigns is a mode-collapse event in precisely the sense this literature worries
about. The literature supplies the **vocabulary** (modes, coverage) and the **concern**, but
not a measurement protocol for a fixed-pool combinatorial library.

**(c) Coverage-style metrics in adjacent generative-model evaluation.**

- Graph generative models: many classical metrics respond poorly to mode collapse, while
  **recall, coverage and some MMD scores track diversity better**
  (arXiv:[2201.09871](https://arxiv.org/abs/2201.09871) **[abstract/snippet]**).
- Generative drug discovery: **uniqueness (fraction of distinct designs) is only a basic
  sanity check**; counting structurally distinct clusters or unique substructures is
  preferred (arXiv:[2501.05457](https://arxiv.org/abs/2501.05457) **[abstract/snippet]**).
  Transferring this to sequences would need a sequence-identity or embedding clustering step.
- Vendi Score: a general diversity metric that captures more than a mode count
  (arXiv:[2210.02410](https://arxiv.org/abs/2210.02410) **[abstract/snippet]**).
- A conditional-generative protein-design preprint proposes two heuristic diversity
  estimates plus a **duality gap** that can flag mode collapse
  (bioRxiv 2021.11.10.467885 **[abstract/snippet]**).
- ProteinZero adds an embedding-space diversity regularizer against mode collapse in online
  RL, with Hamming distance as a separate sequence-level check
  (arXiv:[2506.07459](https://arxiv.org/abs/2506.07459) **[abstract/snippet]**).

**(d) Unique-sequence counts reported in passing.** A 2026 generator benchmark reports
"71,032 sequences remained after filtering, of which 67,664 were unique" and names filtering
as an under-explored bottleneck (bioRxiv 2026.05.08.723784 **[abstract/snippet]**).

### 5.3 Usable conclusion

The coverage metrics we built — distinct variants proposed, fraction of the landscape
reached, attainable-fitness ceiling under the realized proposal set — have **no direct
precedent as a diagnostic**, and the concrete observation that makes them necessary (a
generator reaching 0.18% of the landscape while plain ridge reaches 2.52%, capping attainable
fitness at 0.7036 against a global optimum of 1.0) has one close ancestor: Wu et al.'s 47/90
and 39/90. The framing with the best chance of landing:

> Generator redundancy in combinatorial-library MLDE was **observed and counted** by Wu et
> al. (2019) but never turned into a metric. The GFlowNet literature **names** mode collapse
> as the failure mode and regularizes against it, but measures only within-top-K diversity
> and novelty, not landscape coverage. We report coverage and the **attainable-fitness
> ceiling it implies**, which converts a diversity observation into a bound on what the
> campaign could possibly have achieved.

That last clause is the part with no precedent at all: no located source converts a coverage
measurement into a **ceiling on attainable fitness**. If anything in this document is
publishable on its own, it is that.

---

## 6. Pinned citation list

Everything asserted above, with a stable identifier. Verification markers as in §0.

| # | work | identifier | read |
|---|---|---|---|
| 1 | Yang et al. 2025, ALDE, Nat Commun 16:714 | [10.1038/s41467-025-55987-8](https://doi.org/10.1038/s41467-025-55987-8); PMC11739421 | full text |
| 2 | Greenman, Amini & Yang 2025, PLOS Comput Biol 21(1):e1012639 | [10.1371/journal.pcbi.1012639](https://doi.org/10.1371/journal.pcbi.1012639); PMC11741572 | full text |
| 3 | Chen et al. 2026, CombinGym, bioRxiv | [10.64898/2026.03.24.714074](https://doi.org/10.64898/2026.03.24.714074); repo @ `c5dccd51…`; combingym.org | abstract + repo |
| 4 | Tran et al. 2026, MULTI-evolve, Science | [10.1126/science.aea1820](https://doi.org/10.1126/science.aea1820) | press summary only |
| 5 | "Additive baselines furnish no evidence for epistasis learning by MULTI-evolve", bioRxiv | [10.64898/2026.04.23.719915](https://doi.org/10.64898/2026.04.23.719915) | abstract |
| 6 | Zong 2023, Protein Sequence Design with Batch BO | arXiv:[2303.10429](https://arxiv.org/abs/2303.10429) | snippet |
| 7 | Cheng et al. 2022, ODBO | arXiv:[2205.09548](https://arxiv.org/abs/2205.09548) | full text |
| 8 | Qiu et al. 2024, Tree Search-Based Evolutionary Bandits, AAAI-24 | arXiv:[2401.06173](https://arxiv.org/abs/2401.06173); https://ojs.aaai.org/index.php/AAAI/article/view/29386 | abstract only |
| 9 | Khan et al. 2022, AntBO, Cell Reports Methods | https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(22)00276-4 | snippet |
| 10 | PepCompass 2025 | arXiv:[2510.01988](https://arxiv.org/abs/2510.01988) | snippet |
| 11 | Kmicikiewicz, Fortuin & Szczurek 2025, ProSpero | arXiv:[2505.22494](https://arxiv.org/abs/2505.22494); code github.com/szczurek-lab/ProSpero | full text |
| 12 | Sinai et al. 2020, AdaLead / FLEXS | arXiv:[2010.02141](https://arxiv.org/abs/2010.02141) | full text |
| 13 | Jain et al. 2022, Biological Sequence Design with GFlowNets, ICML | arXiv:[2203.04115](https://arxiv.org/abs/2203.04115); PMLR v162 jain22a | full text |
| 14 | LS-GFN | arXiv:[2310.02710](https://arxiv.org/abs/2310.02710) | snippet |
| 15 | Drummond, Iverson, Georgiou & Arnold 2005, J Mol Biol 350:806 | PMID [15939434](https://pubmed.ncbi.nlm.nih.gov/15939434/); arXiv:[q-bio/0411041](https://arxiv.org/abs/q-bio/0411041) | snippet |
| 16 | Wu, Kaltenbach & Arnold et al. 2019, PNAS 116(18) | [10.1073/pnas.1901979116](https://doi.org/10.1073/pnas.1901979116); arXiv:[1902.07231](https://arxiv.org/abs/1902.07231) | snippet |
| 17 | Carpentier & Munos 2011, MC-UCB, NeurIPS | https://papers.neurips.cc/paper_files/paper/2011/file/872488f88d1b2db54d55bc8bba2fad1b-Paper.pdf; JMLR 16 (2015) carpentier15a | snippet |
| 18 | Phillips et al., eLife 10:e71393 (source of CR9114/CR6261) | [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393) | **not read** |
| 19 | FolDE | PMC12636759 | snippet |
| 20 | Far-red FP active-learning preprint, bioRxiv | [10.64898/2026.08.12.744534](https://doi.org/10.64898/2026.08.12.744534) | snippet |
| 21 | Bloom et al., BMC Biol 5:29, mutational robustness | arXiv:[0704.1885](https://arxiv.org/abs/0704.1885) | snippet |
| 22 | "Improving protein optimization with smoothed fitness landscapes" | arXiv:[2307.00494](https://arxiv.org/abs/2307.00494) | snippet |
| 23 | "How Evaluation Choices Distort the Outcome of Generative Drug Discovery" | arXiv:[2501.05457](https://arxiv.org/abs/2501.05457) | snippet |
| 24 | "On Evaluation Metrics for Graph Generative Models" | arXiv:[2201.09871](https://arxiv.org/abs/2201.09871) | snippet |
| 25 | Vendi Score | arXiv:[2210.02410](https://arxiv.org/abs/2210.02410) | snippet |
| 26 | ProteinZero | arXiv:[2506.07459](https://arxiv.org/abs/2506.07459) | snippet |
| 27 | ProteinGym, NeurIPS 2023 D&B | [10.1101/2023.12.07.570727](https://doi.org/10.1101/2023.12.07.570727); proceedings.neurips.cc/paper_files/paper/2023/hash/cac723e5ff29f65e3fcbb0739ae91bee-Abstract.html | snippet |

---

## 7. Things that contradict our findings or endanger our plan

Ordered by how much damage each could do.

**7.1 — Our ALDE comparison is run at a budget where ALDE's model class was never
demonstrated, and ALDE publishes no sensitivity analysis to cover it.** [§2.1] ALDE's
simulated protocol on these exact two landscapes is 96 + 4×96 = 480; ours is 125 (26%). The
round count matches, so the gap is entirely in data volume per refit. ALDE's published
configurations are DNN ensembles and DKL; its wet-lab run used a DNN ensemble on 216 points.
At 25 initial observations these are an order of magnitude below any demonstrated
training-set size, where one-hot ridge degrades far more gracefully. **The 50.44-vs-30.42
headline is therefore partly an artefact of the adaptation, and the artefact runs in our
favour.** Fig. S73's trajectories do not cover this: they vary cumulative budget at fixed
batch size 96. Mitigation: run the arm grid at ALDE's own 96 + 4×96 schedule on TrpB and GB1
(development data, so free) and report both.

**7.2 — A peer-reviewed benchmark says uncertainty-based acquisition does not beat greedy on
protein landscapes. Our GB1 loss is its expected outcome; our TrpB win is the anomaly.**
[§2.2] Greenman et al. 2025, across eight FLIP tasks and seven UQ methods: uncertainty-based
acquisition "never outperform[s] greedily sampling the sequences with the highest predicted
values". This reframes our two development results. The narrative we would naturally write —
"we win on TrpB, we have a fixable problem on GB1" — is backwards relative to the
literature's prior, which is "greedy wins; a win over greedy needs explaining". Writing it
the literature's way is both more honest and more interesting: the TrpB win becomes a
specific claim about rich shallow shells rather than a general claim about the engine.

**7.3 — Greenman et al. does not support the calibration mechanism we wanted to cite it
for.** [§2.2] It makes no link between miscalibration and over-exploitation, and stratifies
nothing by mutation count. Citing it for that mechanism would be a misrepresentation that a
reviewer who knows the paper would catch immediately. We must demonstrate the mechanism from
our own §10 per-depth rows or drop the claim.

**7.4 — The literature has already published the fixed-depth baseline that fires our own kill
criterion, and it is in Science.** [§1.1(a), §0] MULTI-evolve's singles → all-pairs →
extrapolate-to-depth-7 is a hand-designed depth schedule with wet-lab results. Kill criterion
§13.5.4 says to abandon depth-adaptive allocation if a simpler fixed-depth heuristic matches
it. **That heuristic must be in the arm list from the start.** Adding it only after the
adaptive version wins would be indefensible, and a reviewer will ask for it by name.

**7.5 — eqFP611 blue versus red argues against protein class as a prospective predictor of
shell richness.** [§3.3] Same protein, same library, same 8,192 genotypes; blue falls
monotonically with depth while red has zero hits at depths 0–3 and peaks at depth 8. If
readout alone can invert shell structure, the commercial framing in Question 3 — predict
shallow richness from protein class and use it for target selection — is in trouble as
stated. Library design (§3.4 route 2) is the better-grounded route; protein class is the
weakest.

**7.6 — The depth structure of our two antibody targets may already be published.** [§3.3]
Phillips et al., eLife 10:e71393, titled *Binding affinity landscapes **constrain** the
evolution of broadly neutralizing anti-influenza antibodies*, is the source of CR9114 and
CR6261. If it reports that high affinity requires many substitutions, then the deep-shell
structure of two of our seven independent landscapes is prior knowledge, our H1 ranking is
partly anticipated, and the ordering of ranks 5–7 in §6.5 is less of a blind prediction than
it reads. **This review could not read the paper** (the fetch was refused twice). It must be
read before the H1 result is written up.

**7.7 — Unresolved risk to the novelty verdict: Qiu et al. 2024.** [§1.1(e)] The AAAI full
text was unreadable (bot protection; no arXiv HTML render). Its depth parameter `N` is
described as the maximum allowed mutations from wild type and appears in the regret bound;
the abstract and a third-party summary give no indication that it is scheduled or learned,
but this is **not confirmed from the paper**. Obtain the PDF and check before publishing a
novelty claim.

**7.8 — A preprint disputes MULTI-evolve's central epistasis claim, which cuts both ways.**
[§1.1(a)] *Additive baselines furnish no evidence for epistasis learning by MULTI-evolve*
(bioRxiv [10.64898/2026.04.23.719915](https://doi.org/10.64898/2026.04.23.719915),
**[abstract/snippet]**, not peer reviewed) reports that MULTI-evolve's multimutant
predictions are near-perfectly correlated with an additive model's across all three
applications, and that its networks do not beat an additive model on held-out prediction.
A related eLetter on the Science paper argues the result depends on low-order interactions
dominating. **For us:** (i) it weakens MULTI-evolve as evidence that *epistasis-order-aware
modelling* works, while leaving the depth *schedule* intact as a baseline — so 7.4 stands;
(ii) it is a warning about our own GB1 result, where ridge-greedy beats the engine. If
additive models explain most of the usable signal on these landscapes, then an engine that
loses to ridge is not failing to model epistasis — it is failing to be as simple as the
problem. §13.6 of the pre-registration already names this possibility and is right to.

**7.9 — CombinGym supplies no campaign baselines, so a reviewer expecting the benchmark's
own baselines will be disappointed.** [§2.3] Its nine baselines are fitness predictors scored
on Spearman and NDCG. There is no AL loop in the repository. Our arm list is entirely our
own, which is correct but must be stated, not assumed.

**7.10 — CombinGym changed without a version bump, and one of those commits deleted a data
file.** [§2.3] Eighteen months of no activity, then two commits on 2026-09-16, one of which
removed `Data/Structure/SaCas9(SaCas9-KKH).pdb`. The commit pin and the SHA table in §3.4 of
the pre-registration are doing real work. Keep them, and do not let anyone "update to
latest".

**7.11 — No prior art is not the same as no competition.** [§0] The gap is a gap in
*published method*. All four ingredients of depth-adaptive allocation are individually
standard, and the composition is a natural next step from Zong's static λ plus ODBO's
adaptive radius. The novelty window is narrow and the right framing is "known stratified-
bandit machinery applied to an action space nobody has used", not a claim of methodological
invention.

---

## 8. Open items this review could not close

1. **Qiu et al. 2024 (AAAI-24) full text** — how is `N` set? Largest unresolved risk to the
   §0 novelty verdict. [§1.1(e), §7.7]
2. **Phillips et al., eLife 10:e71393** — does it already publish the depth structure of
   CR9114 / CR6261? Must be read before the H1 write-up. [§3.3, §7.6]
3. **MULTI-evolve, Science full text** — confirm from the paper (not the press summary) that
   the depth schedule is fixed, and extract the exact variant counts per order so it can be
   implemented as a baseline arm. [§1.1(a), §7.4]
4. **CombinGym bioRxiv full text** — confirm from the paper that the k-vs-rest conclusion is
   not reported per-protein or per-assay, and recover the nine algorithms' mapping to the
   five categories. [§2.3]
5. **Zong 2023 full text** — confirm the λ formulation and that λ is never adapted. Currently
   snippet-level only, and it is load-bearing for §0's "the knob exists but is static". [§1.1(b)]
6. **Drummond & Arnold 2005 full text and DOI** — the DOI was not verified by this review;
   PMID 15939434 and the Caltech PDF are the identifiers used here. Confirm before citing a
   DOI. [§3.2(b)]
