# Is the H1 ordering already in the literature?

**Status:** literature check. Written 2026-10-31, before any cross-protein arm has been run.
It answers one question about `docs/preregistration-cross-protein.md` §6.5: for each landscape
in the pre-registered per-protein enrichment ordering, how much of that landscape's position is
recoverable from published figures by a reader who has never seen our data.

**Scope and non-scope.** This document does not re-derive any number from this repository. The
enrichments, shell cuts and depth profiles it quotes are carried in from
`docs/preregistration-cross-protein.md` §6.3 and §6.5 and are treated as established. It does
not modify the pre-registration, the registry, or any code. It changes nothing about what will
be run; it changes what may honestly be claimed afterwards.

**The ordering under audit** (`docs/preregistration-cross-protein.md` §6.5, per-protein rows):

    CreiLOV 1.875x  >  eqFP611/mTagBFP2 1.439x  >  CR9114 0.584x  >  CR6261 0.105x

n = 4, so the smallest attainable one-sided p is 1/24 = 0.042 and only a perfect ordering
reaches significance. A position that was already published is not a prediction, so the count
of genuinely blind positions is the thing that decides whether the test means anything.

**Verification convention**, carried over from `docs/prior-art-depth-allocation.md`:

- **[full text]** — the relevant section was read.
- **[figure read]** — the published figure image itself was obtained and inspected, not just
  its caption.
- **[caption]** — only the figure legend was read.
- **[abstract/snippet]** — only an abstract or search extract was read.

Every primary source below was obtained as a PDF, a publisher XML full text, or a figure image
and read locally. Nothing in §2–§6 rests on an abstract.

---

## 0. The verdict, first

> **Three of the four H1 positions are recoverable from published figures. One — eqFP611's —
> is recoverable in direction but not in rank, because its rank is an artefact of averaging
> three readouts in a way no paper reports. So the per-protein ordering is, at best,
> **one position blind out of four**, and arguably **zero**.**
>
> **The H1 ordering is therefore not a discovery. H1 is not thereby void — what is unpublished
> is the other half of the correlation, our optimizer's per-landscape effect sizes — but any
> sentence presenting the ordering itself as a finding of this project must be struck.**

| landscape | verdict | the thing that settles it |
|---|---|---|
| eqFP611 blue vs red | **ANTICIPATED** | Poelwijk Fig 6a: bright genotypes are blue for mutational steps 0–6 and red for steps 7–13; the colour switch sits at the neck, i.e. exactly at our floor(13/2) = 6 cut |
| eqFP611 protein-mean rank (2nd of 4) | **PARTLY ANTICIPATED** | direction published; the 1.439x value is a mean over blue/red/combined, a construct with no published counterpart |
| CreiLOV | **ANTICIPATED** | Chen Fig 3c publishes the fraction of each mutation-count shell that beats WT, falling monotonically from ~0.52 at depth 2 to ~0 by depth 13 |
| CR6261 | **PARTLY ANTICIPATED, leaning anticipated** | Phillips Fig 1E/1F: the CR6261 germline binds *neither* antigen, and H9's distribution reaches its ceiling only at depth 9–10 — while CR9114's germline already binds H1 at −logK_D 8.4 |
| CR9114 | **ANTICIPATED** (carried in) — **but the citation behind it is the wrong figure; see §5** | Phillips Fig 1D, not Fig 5A |
| TrpB (development) | **NOT ANTICIPATED** | Johnston 2024 publishes no fitness-versus-mutation-count object of any kind |
| GB1 (development) | **NOT ANTICIPATED** | Wu 2016 publishes a global fitness distribution, never stratified by mutation count |

The §6 result is the one worth reading twice: **H1's own derivation was not available from the
literature.** The covariate was found on two landscapes whose shell structure nobody has
published, and is being tested on four whose shell structure everybody could have read. That is
an awkward shape for a pre-registered test and §7 says what to do about it.

---

## 1. What counts as anticipation, and the distinction everything turns on

Our covariate is a **hit density per depth shell**: the fraction of the shell `1 ≤ d ≤ cut`
that falls in the landscape's top tier, divided by the base rate. Three different published
objects get confused with it, and they are not the same claim:

1. **Maximum fitness per mutation count.** Says where the single best variant at each depth
   sits. Says nothing about how many variants at that depth are good. A shell can contain the
   global optimum and still be hit-poor, and — the case that actually bites here — a shell can
   contain a variant near the ceiling while holding a vanishing fraction of the global top 1%.
2. **Functional fraction per mutation count.** The fraction of a shell that works at all.
   Closer, but the threshold is at the floor, not at the top tier. A well-adapted parent makes
   its shallow shell rich in *functional* variants and poor in *improved* ones.
3. **Hit density per shell at a high threshold.** Ours.

The rule applied throughout: an object of type 3 at any high threshold counts as anticipating
our covariate. An object of type 2 counts as anticipating its *direction* where the published
threshold is near the parent's own value, because the parent is the reference the optimizer
starts from. An object of type 1 counts for nothing, and §5 shows a case where reading type 1
as type 3 gets the sign wrong.

A second rule, applied against ourselves: **a published figure anticipates a position only if a
reader could reach that position from it without being told our answer first.** Where that
requires knowing our shell cut, it is recorded as a condition, not waved past.

---

## 2. eqFP611 / mTagBFP2 — Poelwijk, Socolich & Ranganathan 2019

Poelwijk FJ, Socolich M, Ranganathan R (2019). *Learning the pattern of epistasis linking
genotype and phenotype in a protein.* Nature Communications 10:4213.
DOI [10.1038/s41467-019-12130-8](https://doi.org/10.1038/s41467-019-12130-8).
Full text: [PMC6746860](https://pmc.ncbi.nlm.nih.gov/articles/PMC6746860/). **[full text,
publisher PDF; Fig 6 read as an image]**

**Verdict: the blue-versus-red contrast is ANTICIPATED. The protein's rank in the H1 ordering
is PARTLY ANTICIPATED.**

### 2.1 The library design alone carries most of it

**[full text]** The library is not "a fluorescent protein with 13 mutated sites". It is the
complete 2^13 = 8,192-genotype space **between two named parents**: mTagBFP2, a bright blue
variant, and mKate2, a bright deep-red variant of eqFP611, "separated by thirteen mutations"
(Fig 1a and the opening of Results). The frozen corpus takes the blue parent as the depth-0
reference — the pinned wild-type file is literally named
`Data/WT Sequence/eqFP611(mTagBFP2).fasta` — so **depth 0 is the blue protein and depth 13 is
the red protein, by construction of the published library.**

That single published fact predicts the sign of both readouts' enrichment before any
measurement. Red fluorescence is the phenotype of the depth-13 genotype; it cannot be
shallow-enriched. Blue fluorescence is the phenotype of the depth-0 genotype; a monotone
decline away from it is the default expectation, not a finding.

This disposes of the question the task asked about the parent being blue-optimised. The paper
does not need to *discuss* mTagBFP2 as a poor starting point for red — it defines the library
as the path from mTagBFP2 to mKate2, which says the same thing more strongly.

### 2.2 Fig 6a is a published depth-resolved colour map, and the switch lands on our cut

**[figure read]** Figure 6a lays out every above-threshold genotype (brightness y > 0.73, the
value of the red parent) in columns **by mutational step from mTagBFP2 (step 0) to mKate2
(step 13)**, with each node coloured by its fluorescence colour. The published picture is
unambiguous: **every node from step 0 through about step 6 is blue, every node from about step
7 through step 13 is red**, with one or two mixed nodes at the boundary. The text states it in
words: the solution space is "shaped like a dumbbell, with two densities of functionally bright
sequences near to the parental genotypes connected by a narrow neck", and "**the colour
switches at the neck**".

Our pre-specified shallow shell for a 13-site library is `1 ≤ d ≤ floor(13/2) = 6`. The
published colour boundary is at step 6/7. **The cut falls exactly on the published switch.**

So a reader with Fig 6a and our cut rule predicts: blue's hits are entirely inside the shallow
shell (enrichment well above 1), red's are entirely outside it (enrichment well below 1). Our
measurements are blue 1.854x and red 0.634x, with blue falling monotonically from 30.8% at
depth 1 and red holding **zero** hits at depths 0–3 and peaking at depth 8. The depth-8 peak is
one step past the published neck.

### 2.3 Two further published statements in the same direction

- **[full text]** "the data show that L63M — the only chromophore mutation — is necessary but
  insufficient on its own to produce the red chromophore. Instead, **red fluorescence requires
  several specifically-ordered mutational steps** after L63M". This is an explicit published
  statement that the red phenotype needs multiple mutations.
- **[full text, Methods]** The red parent's apparent brightness in this assay is **0.172** of
  the blue parent's (0.772 spectral brightness ratio × 0.222 for the assay's excitation and
  emission bands), and the channels are normalised to that ratio. A reader who works through
  the Methods therefore knows that the top of the *combined* brightness scale sits on the blue
  side of the path — i.e. shallow. Our combined enrichment is 1.830x, essentially blue's.

### 2.4 What the paper does **not** publish

- **No per-channel fitness.** The published phenotype is a single scalar ȳ, the quadratic
  combination of the normalised blue and red enrichments (Methods, Eq. for `x_a`). The
  per-channel values `E_a,blue` and `E_a,red` exist in the deposited data and are what
  CombinGym exposes as the `blue` and `red` columns, but **no figure in the paper plots blue
  fitness or red fitness against mutation count.**
- **No hit-density curve at a high threshold.** Fig 6b plots `log(fraction func)` against
  mutational step and is symmetric and U-shaped — 0 at step 0, minimum ≈ −2.1 at step 7, back
  to 0 at step 13. Its threshold is "functional and in the connected network", which is type 2
  of §1, not type 3. Read alone, Fig 6b predicts a *combined* enrichment near 1.0 by symmetry;
  the asymmetry that gives us 1.830x comes from the 0.172 normalisation, which is in the
  Methods rather than the figure.
- **No numbers resembling 1.854x / 0.634x.**

### 2.5 What a knowledgeable reader could have said without our data

> "The eqFP611 library runs from a blue parent at depth 0 to a red parent at depth 13, and
> Poelwijk Fig 6a shows the bright genotypes are blue up to about step 6 and red after it. So
> the blue readout's best variants sit in the shallow half and the red readout's sit in the
> deep half. Split the library at the midpoint and blue will be hit-enriched, red hit-depleted,
> by a wide margin."

That is H2's primary instance, stated in advance, from one published figure. The contrast
`eqfp611_blue` vs `eqfp611_red` — described in `docs/preregistration-cross-protein.md` §6.7.1
as "the cleanest available within-protein manipulation of the covariate" — **is clean, and it
is also already in print.**

What a reader could *not* have said is where the mean of blue, red and combined lands relative
to CreiLOV. 1.439x is the average of three readouts of one library, two of which (blue 1.854,
combined 1.830) are near-duplicates of each other by construction. No publication reports that
average, and its value depends on our decision to admit all three channels as separate
landscapes and then average them. **That is the one part of eqFP611's position that is ours.**

---

## 3. CreiLOV — Chen et al. 2023

Chen Y, Hu R, Li K, Zhang Y, Fu L, Zhang J, Si T (2023). *Deep Mutational Scanning of an
Oxygen-Independent Fluorescent Protein CreiLOV for Comprehensive Profiling of Mutational and
Epistatic Effects.* ACS Synthetic Biology.
DOI [10.1021/acssynbio.2c00662](https://doi.org/10.1021/acssynbio.2c00662).
PMID [37066862](https://pubmed.ncbi.nlm.nih.gov/37066862/); full text
[PMC10204710](https://pmc.ncbi.nlm.nih.gov/articles/PMC10204710/). **[full text, publisher XML;
Fig 3 read as an image]**

**Verdict: ANTICIPATED. This is the strongest case of the four, and it is not close.**

### 3.1 The paper publishes a hit-density-versus-depth curve

**[figure read]** Figure 3 contains two panels that are, between them, our covariate at a lower
threshold:

- **Fig 3b** — violin plots of fluorescence against **number of mutated residues**, from `wt`
  through 15, with median and quartile lines. The distributions shift downward monotonically;
  at depths 1–3 the bulk sits at or above the WT level, and by depth 10 the bulk is roughly an
  order of magnitude below it.
- **Fig 3c** — "**Fraction of sampled genotypes** observed to be **>WT**" plotted against number
  of mutated residues, 2 through 15. The observed curve (red solid) falls monotonically from
  **≈0.52 at depth 2** to **≈0 by depth 13–14**, with the complementary <WT curve (blue solid)
  rising to 1.

**Fig 3c is a type-3 object** in the sense of §1 — a fraction of each depth shell exceeding a
threshold, plotted against depth — with the threshold at the parent's own fitness rather than
at the top 1%. It is the closest published analogue to our covariate found anywhere in this
audit.

### 3.2 The text states the conclusion outright

**[full text]** "When separating the CreiLOV variants into different groups based on mutation
numbers, we noticed **the overall fluorescence decreased monotonically with the increasing
mutation numbers** (Figure 3b,c)." And, reinforcing it: "**negative epistasis became more
prevalent with increasing numbers of mutations**" (Fig 3e, which plots the fraction of each
shell showing strong negative epistasis rising from ~0 at depth 2 to ~0.5 at depth 15).

### 3.3 The ProteinGym curation adds nothing and conceals nothing

**[full text of the reference file]** `PHOT_CHLRE_Chen_2023` in ProteinGym's
`reference_files/DMS_substitutions.csv` is this assay: 167,529 total mutants, 2,122 singles,
165,407 multiples, binarization cutoff 0.6317 by median, pointing at the same DOI. The
reference file carries **no depth-stratified fitness statistic** — its depth fields are counts
of singles and multiples, nothing more. ProteinGym distributes the full label table, so anyone
can compute the covariate; but computing it is doing our analysis, not reading a published
result, and it is Chen Fig 3b,c that does the anticipating here.

### 3.4 The honest caveat

Fig 3c's threshold is "brighter than WT", not "in the top 1%". These can come apart: a
landscape can show a monotone decline in the above-WT fraction while its very best variants
require several mutations. They do not come apart here — Fig 3b shows the *upper* tail of each
violin collapsing with depth too, so the top tier moves with the bulk. But the distinction is
real and the write-up should state it rather than claim Fig 3c *is* our covariate.

### 3.5 What a knowledgeable reader could have said without our data

> "Chen Fig 3c shows the fraction of CreiLOV genotypes beating wild type falling monotonically
> from about half at two mutations to nothing by thirteen. CreiLOV's shallow shell is the
> richest part of its landscape, and by a large factor."

Our measured profile — hit rate 9.66% at depth 2 falling to 0.12% at depth 12, enrichment 9.18x
at the depth ≤ 2 cut and 1.875x at the primary floor(15/2) = 7 cut — is that statement with our
threshold substituted. **CreiLOV's first place is published.**

---

## 4. CR6261 — Phillips et al. 2021

Phillips AM, Lawrence KR, Moulana A, Dupic T, Chang J, Johnson MS, Cvijovic I, Mora T,
Walczak AM, Desai MM (2021). *Binding affinity landscapes constrain the evolution of broadly
neutralizing anti-influenza antibodies.* eLife 10:e71393.
DOI [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393).
**[full text, publisher PDF; Figs 1 and 5 read as images]** — obtained by downloading
`https://cdn.elifesciences.org/articles/71393/elife-71393-v3.pdf` from the asset list at
`https://api.elifesciences.org/articles/71393` and extracting text locally. The HTML route is
refused by a content filter, as `docs/prior-art-depth-allocation.md` §3.3 recorded; the PDF
route worked first try.

**Verdict: PARTLY ANTICIPATED, leaning anticipated. The paper's prose points away from our
ordering; the paper's figures point straight at it, and the figures win.**

### 4.1 The reconciliation holds, and it is not the whole story

The established nuance — that the paper says "improvements in binding can occur early on for
both antigens" with 92% and 81% of variants binding, while we measure 0.210x and 0.000x — is
correct and is not re-derived here. It reconciles because Fig 5B is a **maximum** per shell
(type 1 of §1) and we measure **top-tier density** per shell (type 3).

But resting CR6261's "unanticipated" status on that reconciliation would be wrong, because the
paper publishes a type-3-adjacent object too, and it is in Figure 1.

### 4.2 Figure 1F is a per-shell distribution, not a maximum

**[figure read]** Figure 1F plots, for CR6261, **boxplots of −logK_D against number of somatic
mutations, one panel per antigen**. A boxplot per shell shows median, interquartile range and
whiskers — i.e. where the *mass* of each shell sits, not merely its extreme. Reading the two
panels:

- **H1**: the germline (depth 0) sits pinned at the ≤7 floor. Medians climb through depths 2–5
  and the box reaches the ceiling band (≈9.4–9.5) from roughly depth 6–7 onward.
- **H9**: depths 0–1 are at the floor; the distribution reaches its own ceiling (≈10.1) only at
  roughly **depth 9–10**, visibly later than H1 reaches its ceiling.

**H9's top band arrives later than H1's, in a published figure.** That is our H1 > H9 ordering
(0.210x > 0.000x), and the text states it in words as well: among sequential selection
scenarios, those beginning with H1 are more likely "as the first two mutational steps can
improve affinity to H1 more than H9, and **mutations late in maturation can improve affinity to
H9 more than H1** (Figure 1F, Figure 5B)". The paper itself cites Fig 1F for this, not only
Fig 5B.

### 4.3 The germline comparison is the part that puts CR6261 below CR9114

This is the piece that matters for the H1 rank, and it is published in the same figure:

- **[figure read]** Fig 1E and Fig 1F: the **CR6261 germline binds neither H1 nor H9** — both
  coordinates are pinned at the ≤7 detection floor, and Fig 5B's depth-0 maxima are ≈6.9–7.0.
- **[figure read]** Fig 1C and Fig 1D: the **CR9114 germline already binds H1** at −logK_D ≈
  8.4, and 97% of the CR9114 library has measurable H1 affinity.

A library whose depth-0 genotype does not bind at all cannot have a top-tier-rich shallow
shell. A library whose depth-0 genotype already binds respectably can. **Side by side in one
published figure, CR6261's shallow shells are visibly poorer than CR9114's**, which is our
0.105x < 0.584x.

Supporting this from the Methods: **[full text]** "For CR6261, **all mutations at the first
mutational step are neutral** (with the exception of one mutation that improves affinity for H1
only)", and from Results: T29P "is the only single mutation that can improve H1 affinity,
albeit rather modestly". Eleven sites, one modestly-improving single. Our depth ≤ 2 shell for
CR6261 H1 holds 58 variants and **zero** hits; the published sentence predicts that.

### 4.4 What points the other way, and why it loses

The paper's framing — abstract, introduction and discussion — contrasts CR9114 as the
constrained antibody against CR6261 as the unconstrained one: "the acquisition of breadth to
more similar antigens for CR6261 is less constrained", "most variants can bind both antigens
(92% for H1, 81% for H9)", "many likely paths that continually improve in binding to both
antigens". A reader who worked from the prose would rank CR6261 **above** CR9114 and get the
adjacency wrong.

That is a genuine source of blindness, and it is why this verdict is PARTLY rather than fully
anticipated. But it is blindness produced by reading the summary instead of the figure, and the
figure is Figure 1 of the paper. We should not claim credit for a position that Figure 1
displays merely because the abstract is phrased in a different currency.

### 4.5 What a knowledgeable reader could have said without our data

> "Phillips Fig 1E/1F shows the CR6261 germline binding neither antigen while the CR9114
> germline already binds H1 at 8.4, and Fig 1F shows H9's affinity distribution reaching its
> ceiling only around nine or ten mutations. CR6261's shallow shells will be the poorest in the
> set, with H9 poorer than H1 — below CR9114, and last overall."

Conditional on being told our floor(S/2) cut and our 1%/2% tier definitions. Without those, the
reader gets the direction and the H1 > H9 split but not the numeric rank.

---

## 5. CR9114 — carried in, with a correction to the citation that supports it

**Verdict: ANTICIPATED. Established previously and not re-derived here.** But this audit read
the figures, and the figure currently cited for the verdict does not support it — a different
figure in the same paper does. Recording this because the write-up will cite one of them.

**The problem with Fig 5A.** **[figure read]** Figure 5A plots **maximum** −logK_D against
number of somatic mutations for each of CR9114's three antigens. For **H1** — which is the only
CR9114 landscape in the frozen corpus (`docs/preregistration-cross-protein.md` §3.2: "`Clean/`
at this pin holds CR9114 H1 only") — the maximum is ≈8.4 at depth 0, reaches ≈9.4 by **depth 1**
and ≈9.5 by depth 2, then flat to depth 16. Read as a statement about our covariate, **Fig 5A
says CR9114/H1's best binders are available immediately** — the opposite sign to our 0.584x.
The deep structure visible in Fig 5A is in the **H3 and influenza B** traces, and the "eight
mutations required for influenza B binding" result is likewise about influenza B. Neither
landscape is in the frozen corpus. The main text reinforces the hazard: "**For H1, only a few
mutations are needed to improve from the germline affinity.**"

**The figure that does support it is Fig 1D.** **[figure read]** Figure 1D plots **boxplots of
H1 −logK_D against number of somatic mutations**, 0 through 16. The medians rise monotonically
from ≈8.5 at depth 0 to the ≈9.5 ceiling band, and the **interquartile range only sits inside
the ceiling band from roughly depth 9–10 onward**. The upper whiskers reach near-ceiling values
early — which is exactly what Fig 5A reports — while the *mass* of each shallow shell stays a
full log unit below. **The fraction of a shell that is top-tier rises monotonically with depth**,
which is our 0.584x with 227 hits in the 38,835-variant depth 1–8 shell and zero in the
134-variant depth ≤ 2 shell.

This is the cleanest illustration in the whole audit of the §1 rule, and it cuts in our favour
on the facts and against us on the conclusion: Fig 5A would have told a reader the wrong thing,
and Fig 1D tells them the right thing. The CR9114 verdict stands; **the supporting citation
should be Phillips Fig 1D (per-shell affinity distributions), not Fig 5A (per-shell maxima),
and the claim should be scoped to H1 rather than to the influenza B result.**

---

## 6. TrpB and GB1 — the development landscapes, and the reason to check them

The question asked of these two: if their shallow-shell richness is also published, then H1
itself was derivable from the literature rather than discovered from our data.

**It is not. Both are clean negatives.**

### 6.1 TrpB — Johnston et al. 2024: NOT ANTICIPATED

Johnston KE, Almhjell PJ, Watkins-Dulaney EJ, Liu G, Porter NJ, Yang J, Arnold FH (2024).
*A combinatorially complete epistatic fitness landscape in an enzyme active site.* PNAS
121:e2400439121. DOI [10.1073/pnas.2400439121](https://doi.org/10.1073/pnas.2400439121).
Full text [PMC11317637](https://pmc.ncbi.nlm.nih.gov/articles/PMC11317637/). **[full text,
publisher XML; all figure captions read]**

The paper publishes **no fitness-versus-mutation-count object of any kind** — no distribution
per mutation count, no functional fraction per mutation count, no maximum per mutation count.
Its figures are: assay schematic and site selection (Fig 1); pairwise epistasis types split by
fitness quartile and by site pair (Fig 2); path maps and path-accessibility ECDFs from the
parent `VFVS` to the top variant `AIKG`, plus directed-evolution simulations (Fig 3); kinetic
and structural characterisation (Fig 4); off-the-shelf predictor classification (Fig 5).
Mutation depth appears once, as the maximum step count in the direct-path analysis.

What a reader *could* take from it points the wrong way: Fig 3A's path map runs from the parent
`VFVS` to the top variant `AIKG`, which differs at all four sites. **The published global
optimum is a quadruple mutant** — a type-1 fact that, read as a density claim, suggests TrpB is
deep-favouring. Our measurement is the opposite: 7.49x shallow enrichment, 168 of the top 1% in
the 2,242-variant depth 1–2 shell.

### 6.2 GB1 — Wu et al. 2016: NOT ANTICIPATED

Wu NC, Dai L, Olson CA, Lloyd-Smith JO, Sun R (2016). *Adaptation in protein fitness landscapes
is facilitated by indirect paths.* eLife 5:e16965.
DOI [10.7554/eLife.16965](https://doi.org/10.7554/eLife.16965).
Full text [PMC4985287](https://pmc.ncbi.nlm.nih.gov/articles/PMC4985287/). **[full text,
publisher XML; all figure captions read]**

The paper reports that **2.4% of mutants are beneficial** (fitness > 1) and that 19.7% are
lethal, and Figure 1—figure supplement 4 shows the fitness distribution of all non-lethal
genotypes — **globally, never stratified by mutation count**. Figure 1C analyses the WT
neighbourhood, but for epistasis *types*, not hit density. Nothing in the paper plots fitness
or hit fraction against number of substitutions.

Adjacent curation does not fill the gap either: FLIP's and ProteinGym's GB1 `1-vs-rest`,
`2-vs-rest` and `3-vs-rest` splits, and CombinGym's `{0,1,2,3}-vs-rest` leaderboards, are
stratifications of **predictor transfer**, not of hit density — the point
`docs/prior-art-depth-allocation.md` §2.3 and §4.3 already make.

### 6.3 The consequence, which is not the one the question expected

The hypothesis was **not** derivable from the literature. H1 was built from a TrpB-versus-GB1
contrast (7.49x against 3.32x) that exists in no publication, and this audit found no
cross-protein meta-analysis of shell enrichment either — consistent with
`docs/prior-art-depth-allocation.md` §3.1's clean negative.

So the shape of the problem is the reverse of the worry: **the covariate's derivation is
original and its test set is pre-published.** The ordering H1 predicts was knowable; the
relationship H1 asserts was not.

---

## 7. Overall assessment

### 7.1 How many of the four H1 positions are genuinely blind

**At most one, and that one only in a bookkeeping sense.**

| position | blind? | why |
|---|---|---|
| CreiLOV, 1st | **no** | Chen Fig 3b,c publish a monotone hit-fraction decline with depth. The only published object of our statistic's type in the whole set. |
| eqFP611, 2nd | **direction no, rank partly** | Poelwijk Fig 6a fixes blue-shallow/red-deep with the colour switch on our cut. The 1.439x figure is a mean over three readouts that no paper treats as three landscapes. |
| CR9114, 3rd | **no** | Phillips Fig 1D. Established before this audit; the supporting citation needs the correction in §5. |
| CR6261, 4th | **no, with one honest caveat** | Phillips Fig 1E/1F: a germline that binds neither antigen, H9 reaching its ceiling at depth 9–10, "all mutations at the first mutational step are neutral". The paper's *prose* would mislead; its Figure 1 would not. |

Three published figures — Chen Fig 3c, Poelwijk Fig 6a, Phillips Fig 1D/1F — between them fix
all four positions in direction, and fix three of the four in rank once our shell cut is
supplied. **The pre-registered ordering is, to a good approximation, a transcription of
published material into our units.**

### 7.2 Does H1 remain a meaningful test?

**Yes, but it is a narrower test than §6.1 of the pre-registration describes, and the
difference must be written down.**

H1 is a Spearman correlation between two columns. Only one of them is now in question:

- **The covariate column is published.** Nothing in it is a prediction.
- **The outcome column is not.** Nobody has run the VeritasBio reliability policy, or ridge, on
  these seven landscapes at a 125-measurement budget. The per-landscape effect sizes are
  genuinely unknown, and the §6.5 prediction built on them — that reliability will **lose to
  ridge on every independent landscape**, most severely on `cr6261_h9` and `cr6261_h1`, least
  on `creilov` and `eqfp611_blue` — is a real, falsifiable, pre-registered forecast. It
  survives this audit intact.

So the correct description of what H1 tests is:

> **Given a landscape ordering that was already available in the literature, does our
> optimizer's measured advantage track it?**

That is a legitimate pre-registered question. What it is *not* is evidence that we discovered
the ordering, and §6.1's framing of the covariate as "the one structural difference found"
should be read as applying to TrpB-versus-GB1 (where it is true and unpublished, §6.3 above)
and not to the independent set.

### 7.3 The sharper problem this creates for H2

H2 is the better-powered hypothesis (48 seeds against n = 4 proteins) and the pre-registration
calls `eqfp611_blue` vs `eqfp611_red` "the cleanest available within-protein manipulation of
the covariate". It is. It is also the *most* anticipated single item in this audit: Poelwijk
Fig 6a draws the manipulation directly, with the colour switch at our cut.

This does **not** weaken H2 as a mechanism test — the mechanism claim is about our policy's
behaviour, and that is unpublished. It does mean H2's design section must stop presenting the
blue/red shell-structure contrast as something we found. The sentence in §6.7.1 that the
comparison "varies the covariate and holds the rest fixed" is correct and should stay; the
implication that the variation was discovered here should go. The correct phrasing is that
**Poelwijk et al. built a library in which the covariate varies between two readouts, and we
are using it as a natural experiment** — which is a better argument, because it means the
manipulation was not chosen by us after seeing our own results.

### 7.4 Recommended dispositions

Offered, not applied — this document modifies nothing.

1. **Amend `docs/preregistration-cross-protein.md` §16** with a dated entry recording that the
   §6.5 ordering is substantially recoverable from Poelwijk Fig 6a, Chen Fig 3b,c and Phillips
   Fig 1D/1F, naming this document. Before the runs, so it is an amendment and not a defence.
2. **Add the three citations to §6.5 itself**, alongside each protein's enrichment row, so the
   ordering is never again reported as if it were measured here first.
3. **Fix the CR9114 citation** per §5: Phillips Fig 1D, scoped to H1, not Fig 5A and not the
   influenza B eight-mutation result.
4. **Re-word §6.7.1** per §7.3: the blue/red covariate contrast is a published property of the
   Poelwijk library, used as a natural experiment.
5. **Do not weaken the §6.5 across-the-board-failure prediction.** It is the part of H1 that is
   genuinely blind and it is the part worth defending.
6. **Note for the reserve promotion** (§6.5 of the pre-registration): if the set is expanded to
   n = 7 proteins, the same audit should be run on HIV-1 protease, SaCas9, SpCas9 and RhlA
   *before* promotion, for the same reason. Promoting four more pre-published positions would
   take the exact p-value from 1/24 to 1/5,040 without adding a single blind prediction.

---

## 8. Sources

All obtained and read locally. Journal HTML was not used for eLife 71393 (content-filter
refusal, as previously recorded); the PDF was taken from the asset list at
`https://api.elifesciences.org/articles/71393`.

| source | DOI / stable URL | what was read |
|---|---|---|
| Poelwijk FJ, Socolich M, Ranganathan R (2019) Nat Commun 10:4213 | [10.1038/s41467-019-12130-8](https://doi.org/10.1038/s41467-019-12130-8) · [PMC6746860](https://pmc.ncbi.nlm.nih.gov/articles/PMC6746860/) | full text (publisher PDF); Fig 6 image |
| Chen Y, Hu R, Li K, Zhang Y, Fu L, Zhang J, Si T (2023) ACS Synth Biol | [10.1021/acssynbio.2c00662](https://doi.org/10.1021/acssynbio.2c00662) · [PMC10204710](https://pmc.ncbi.nlm.nih.gov/articles/PMC10204710/) | full text (publisher XML); Fig 3 image |
| Phillips AM, Lawrence KR, Moulana A, Dupic T, Chang J, Johnson MS, Cvijovic I, Mora T, Walczak AM, Desai MM (2021) eLife 10:e71393 | [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393) | full text (publisher PDF); Figs 1 and 5 images |
| Johnston KE, Almhjell PJ, Watkins-Dulaney EJ, Liu G, Porter NJ, Yang J, Arnold FH (2024) PNAS 121:e2400439121 | [10.1073/pnas.2400439121](https://doi.org/10.1073/pnas.2400439121) · [PMC11317637](https://pmc.ncbi.nlm.nih.gov/articles/PMC11317637/) | full text (publisher XML); all figure captions |
| Wu NC, Dai L, Olson CA, Lloyd-Smith JO, Sun R (2016) eLife 5:e16965 | [10.7554/eLife.16965](https://doi.org/10.7554/eLife.16965) · [PMC4985287](https://pmc.ncbi.nlm.nih.gov/articles/PMC4985287/) | full text (publisher XML); all figure captions |
| ProteinGym `DMS_substitutions.csv`, row `PHOT_CHLRE_Chen_2023` | https://github.com/OATML-Markslab/ProteinGym | reference file row and header |
| CombinGym @ `c5dccd51`, `Data/WT Sequence/eqFP611(mTagBFP2).fasta` | https://github.com/sitonglab/CombinGym | filename only, to confirm the depth-0 reference is the blue parent |

**Could not obtain:** nothing. Every source this audit needed was retrieved in full.
