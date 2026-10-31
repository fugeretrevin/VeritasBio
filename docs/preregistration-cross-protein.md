# Pre-registration: independent cross-protein evaluation

**Status:** frozen protocol. Written and committed **before any cross-protein arm has been
run and before any cross-protein performance figure has been examined.** That ordering is
the only thing that makes this document worth anything; if it is edited after results are
seen, the edit must appear in §16 as a dated amendment with its reason, and the original
text must remain readable in git history.

**Frozen on:** 2026-10-31.
**Corpus pin:** CombinGym @ `c5dccd514a9acfe1089a5f1555faea98afd78019` (§3).
**Machine-checkable companion:** `scripts/freeze_cross_protein_corpus.py` plus the lock file
it writes at `scripts/data/cross_protein_corpus.lock.json`.

This is not an architecture decision record, so it does not live in `docs/adr/`. It sits
alongside `docs/545.md`, `docs/NEXT_STEPS.md` and `docs/ENGINE_IMPROVEMENT_QUESTIONS.md`,
which are the repository's measurement-and-state documents. Every figure attributed to this
repository is reproducible from it; regenerate rather than cite.

---

## 1. What is being tested, and what is not

The engine's results to date are confined to two four-site saturated libraries, and
algorithm selection, component ablation and configuration choices were all made against
those two. The question here is narrow:

> On protein families that were never used to select an algorithm, a component, a prior or
> a hyperparameter, does the VeritasBio reliability policy produce more **distinct**
> high-fitness variants per measurement than simple baselines at a fixed small budget?

Three things are explicitly **not** being tested, and no result from this study may be
reported as evidence for them: that the engine is a general-purpose optimizer; that any
individual component contributes to the result; that the measured advantage transfers to
wet-lab campaigns. The second is the subject of the 2^5 joint ablation
(`docs/NEXT_STEPS.md` §3.1) and the third of nothing in this repository.

### 1.1 Protocol deviation that must accompany every external claim

The protocol is **25 shared initial observations, then 4 rounds of 25, for 125 unique
measurements**. This is an **adapted** protocol. It is **not** ALDE's published 96-variant
batch schedule. Any comparison against ALDE numbers — ours or the published ones — must
carry that sentence. A result quoted without it is a misquotation.

### 1.2 Established facts carried in (not to be re-derived)

Measured on held-out seeds on the development landscapes. Listed so the predictions below
are anchored, not so they are repeated as findings about anything else.

| fact | value |
|---|---|
| TrpB eligible universe | 159,999 variants; sites 183/184/227/228; parent `VFVS` |
| GB1 eligible universe | 149,360 variants (prepared CSV); sites 39/40/41/54; parent `VDGV` |
| HER2 trastuzumab CDR-H3 | 11,277 labelled variants, 10 sites (99–108), against a 20^10 product |
| headline, TrpB held-out seeds | best arm 50.44 top-1% hits against ALDE's 30.42 |
| GB1 | ridge, random forest, ALDE and frontier statistically indistinguishable |
| front-loading, 50-measurement budget | 6.3x the best baseline's top-1% hits, t = +25.9 |
| front-loading, 125-measurement budget | 1.7x, t = +8.3 — and the best *configuration* also changes with budget |
| raw count vs diversity | reliability spreads (2.77x redundancy, 4.52 clusters); frontier concentrates (7.84x) and is significantly **worse** than ridge on distinct hits while beating it on raw count |
| reliability vs ridge, TrpB | +13.43 raw top-1%, +9.30 diversity-adjusted |
| reliability vs ridge, GB1 | −22.49 raw top-1% |

The last three rows are why the primary endpoint is not a raw count.

---

## 2. Primary endpoint and the reason for it

**Primary endpoint: diversity-adjusted yield at Hamming d ≥ 2.**

Definition, fixed here: take every measured variant in the target's top tier (§9.1), sort
by fitness descending, and walk the list keeping a variant only when it is at Hamming
distance ≥ 2 from every variant already kept. The count of kept variants is the yield. The
filter is applied **strictly after measurement** — it never influences selection, and no
arm is given the distance rule as an objective.

**Mandatory secondary: raw top-tier count.** It is reported for every arm, target and
budget without exception, including where it contradicts the primary. The TrpB/GB1
divergence above is the whole reason both exist: an arm that wins on raw count by
collecting near-duplicates has not done the thing a campaign needs.

This choice was made by the project owner before the corpus was characterised and is not
revisited by this document.

---

## 3. The corpus freeze

### 3.1 The pin

CombinGym has **no tags and no releases**, so the pin is a commit.

| field | value |
|---|---|
| repository | `https://github.com/sitonglab/CombinGym` |
| commit | `c5dccd514a9acfe1089a5f1555faea98afd78019` |
| commit date (UTC) | `2026-09-16T07:33:53Z` |
| tree SHA | `911f19c7ba82e112056f154d3e8b2073f70dc7d1` |
| state at pin time | tip of `main`; 18 commits; 0 tags; 0 releases |

The tree SHA is recorded so the pinned file list can be re-derived without trusting this
document: `GET /repos/sitonglab/CombinGym/git/trees/911f19c7.../?recursive=1` returns
exactly the blob SHAs below.

### 3.2 Reconciling "14 landscapes over 9 proteins"

The README states "14 combinatorial variant datasets, featuring 4 to 16 mutations across 9
proteins, and encompassing over 400,000 meticulously curated data points." Two different
counts are in circulation because **landscape ≠ protein**, and three proteins carry more
than one phenotype. `Data/Data_summary.xlsx` has **9 rows, one per protein**;
`Benchmarks/Spearman/*.xlsx` has **14 distinct (protein, phenotype) keys**. They agree:

| protein | landscapes at this pin | sites | variant space | measured (summary) |
|---|---|---|---|---|
| SaCas9 (KKH) | activity | 8 | 1,296 | 1,296 |
| SpCas9 | activity | 8 | 952 | 948 |
| HIV-1 protease | activity | 11 | 4,608 | 1,176 |
| RhlA | activity_RLs, activity_Rha-(C8-C10), activity_%Rha-(C8-C10) | 11 | 663,552 | 2,346 |
| CreiLOV | fluorescence | 15 | 184,320 | 165,428 |
| eqFP611 / mTagBFP2 | fluorescence_blue, fluorescence_red, fluorescence_combined | 13 | 8,192 | 8,192 |
| CR9114 | binding_H1 | 16 | 65,536 | 48,841 |
| CR6261 | binding_H1, binding_H9 | 11 | 2,048 | 1,917 |
| GB1 | binding | 4 | 160,000 | 149,361 |

3 + 3 + 2 + 6×1 = **14 landscapes over 9 proteins**. Any source citing a different
collection size is describing a different version of CombinGym and **must not be mixed with
this pin**; reconcile it against the table above first.

Note that `Clean/` at this pin holds **CR9114 H1 only**. The CR9114 H3 and B landscapes
from the upstream publication are not in the frozen corpus, so this protocol cannot use
them.

### 3.3 Curation, and what this protocol reads

Only the files the protocol actually reads are pinned: the curated DMS tables, the
wild-type sequences, and the dataset summary. MSAs, predicted structures and the upstream
baseline code are not part of this evaluation and are deliberately not pinned.

Eligible-universe rule, applied identically to every target: a row enters the eligible
universe only if its variant identifier, its mutation depth and its phenotype value are all
present and parseable, and the value is finite. Rows failing that are dropped and the drop
count is recorded. No value is imputed, back-filled or inferred. This is the same rule the
HER2 curation already applies in `benchmark_panoptes_isolated.py`, where 25,114 held-out
rows carrying `DMS_score` exactly 0.0 are an **absent label**, not a measurement of no
binding, and are excluded.

Direction is **higher is better** for every landscape listed.

### 3.4 Hashes

`git_blob_sha1` is authoritative: it is what GitHub's tree API reports for the pinned commit
and what `git hash-object` recomputes. `sha256` is the ordinary digest of the bytes,
independent of git's framing. A file matching both is the file this protocol was frozen
against.

| path | bytes | git blob SHA-1 | SHA-256 |
|---|---|---|---|
| `Data/Data_summary.xlsx` | 13,055 | `45d6067d821f338df36045f56c4959fc3cffc86d` | `a3b873e4eb16bae12b41de5f2d46a77d3fecd43e81cd81ef931d51a1c2b38c8a` |
| `Data/DMS/Clean/eqFP611_clean.xlsx` | 616,607 | `0d7423fe554453675fc2af10ed0320ffdfc21904` | `27e26e9d198e04f81740349399a8e5fb4687f4f2af7b27b0f9d374e04b0a8d63` |
| `Data/DMS/Clean/CreiLOV_clean.xlsx` | 26,263,997 | `9ad0d8ff3cf7c742280c84befcd33780a575b31d` | `5b1bfaa6788808dfd2692a3757e5d5aa30eb6097653543c1e1c60caa079ad27f` |
| `Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx` | 4,403,699 | `d9895e1f7e6a25fba4c6e2aced9843f7ad4a7000` | `1a1053862608769f5c7dd34cb93e003d26bcf8404a57dd8edcaa5b6463640084` |
| `Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx` | 164,703 | `90864f0970ab50b6b9ee0623f1acecd8aac81d3b` | `43778937140871704228516976eccf5d5fac7ac79dcc1de8e94007803da9a0f5` |
| `Data/DMS/Clean/bnAbs_CR6261_H9_clean.xlsx` | 149,858 | `97f2a084427e67b4d017654bb7ad158f55a2c9cc` | `9ddf1989ee1a4a9449b67bba9ff14ab8b0190ca827611e9196542bb4c2de7063` |
| `Data/DMS/Clean/GB1_clean.xlsx` | 4,350,703 | `c4b9a7dd4eedceb10c952ebc897c69362f215a4d` | `1cbbb4554a4cd0a1c4d564cf37965e9a8a30cab360fab351b707774b017b842f` |
| `Data/DMS/Clean/HIV-1_clean.xlsx` | 47,401 | `ad96e7c6b000bb9f1033d91fe297fb41ced62c5c` | `a0ed221928a81e6e3d5b57571867dccfc4771bc7f4db6629c2da2cf8a99a6023` |
| `Data/DMS/Clean/RhlA_clean.xlsx` | 287,350 | `0cb393531ab6fff5fa8f037221d50cf543429dcd` | `d63fdd66be2a9394a89e8321a4e4553d30b2f8b553773d6faddfb5c114c0f876` |
| `Data/DMS/Clean/SaCas9_clean.xlsx` | 126,881 | `76b854e6ec16c09c99efd2e153f6f2a10154f280` | `8f74ec4bfce62a15e53bde8c381803e6f29bad6603614073233c47bc6f5d8648` |
| `Data/DMS/Clean/SpCas9_clean.xlsx` | 113,908 | `00ae9e52824b835ddbd09c569630d813df7414b0` | `48283f3ea546aecf67ed00e64831c69a62f6ba858ca5df018a7ce5d49f74c103` |

Nine wild-type FASTA files are pinned on the same terms; their hashes are in the lock file
rather than inline here because nothing in the analysis depends on reading them by hand.

All **20 pinned files were fetched and verified on 2026-10-31**: every one matches both the
git blob SHA-1 published by GitHub's tree API for the pinned commit and the SHA-256 recorded
in the lock file. `verify` exits 0 on the frozen corpus and exits 1 on a single altered byte.

### 3.5 Making the freeze checkable

```
python scripts/freeze_cross_protein_corpus.py fetch    # download the pinned files
python scripts/freeze_cross_protein_corpus.py record   # write the lock file
python scripts/freeze_cross_protein_corpus.py verify   # exit 0 iff the corpus matches the pin
python scripts/freeze_cross_protein_corpus.py stats    # the pre-run landscape statistics of §9.3 and §6.2
```

`stats` writes `scripts/data/cross_protein_corpus_stats.json`, which holds every number in
§6.5, §9.1, §9.3 and §9.4 together with the full depth-stratified hit-rate profile and Wilson
intervals for each of the 14 landscapes. That file is committed with this document and is the
machine-readable form of the freeze.

`verify` recomputes both hashes for every pinned file and exits non-zero on any mismatch,
missing file or size disagreement. It also refuses to run if the lock file's commit
disagrees with the script's pin, so a silent re-pin cannot pass unnoticed. The script
depends on the standard library only — including its own minimal `.xlsx` reader — because a
freeze that needs `pandas` and `openpyxl` to be verifiable is a freeze that changes when
those packages change. The corpus itself is **not** committed (it lands in the gitignored
`downloads/combingym/`); the pin and the hashes are what is committed.

### 3.6 Known discrepancies in the upstream corpus

These are recorded now so they cannot be discovered later and treated as results. Only D-6
is resolved; the rest stand.

| # | discrepancy | status |
|---|---|---|
| D-1 | `Data_summary.xlsx` reports 948 measured SpCas9 variants; `SpCas9_clean.xlsx` has 5,748 rows of which only **597** carry a parseable `mean`. The eligible universe used here is 597. | unresolved; summary and clean file disagree |
| D-2 | `Data_summary.xlsx` reports 1,917 measured CR6261 variants; the two clean files hold 1,887 (H1) and 1,842 (H9). The summary's own footnote says that where multiple landscapes exist the measured count "would be possibly different", so 1,917 is probably a union. | probable explanation, unverified |
| D-3 | `Data_summary.xlsx` and `GB1_clean.xlsx` both give 149,361 GB1 variants; this repository's prepared CSV has **149,360**. A one-row difference, most likely wild-type inclusion or a de-duplication. | unresolved; see §4.3 |
| D-4 | `Data_summary.xlsx` reports 1,176 measured HIV-1 variants; `HIV-1_clean.xlsx` has 1,175 rows. | unresolved, one row |
| D-5 | RhlA's `Data_summary.xlsx` row cites no publication. | unresolved; provenance unestablished |
| D-6 | RhlA's `mutant` column is not fixed-width (observed widths 2–17), so it cannot be used as the variant identifier the way eqFP611's `mutant` column can. Its `genotype` column *is* fixed-width (11 residues) and is used instead. | resolved by column choice; recorded so the wrong column is not picked later |
| D-7 | Column conventions are inconsistent across the corpus: the variant identifier is a 13-residue combo (eqFP611), an 11-residue combo (RhlA), or a full 56-, 99-, 119-, 121-, 1,053- or 1,368-residue sequence (GB1, HIV-1, CreiLOV, the antibodies, the Cas9s); the depth column is `n_mut` or `num_mut`; the value column is `mean`, `h1_mean`, `h9_mean`, `fitness`, `Mean`, `blue`/`red`/`combined` or three named activity columns. | handled by the per-landscape registry in `freeze_cross_protein_corpus.py`; a silent mis-mapping would be invisible, so each landscape's resolved column names are printed in the stats artifact |

Where the summary and the clean file disagree, **the clean file wins**, because the clean
file is what is read.

---

## 4. Two-level validation structure

### 4.1 Levels

There are **two validation levels**, development and final test. The reserve below is not a
third level; it is unallocated holdback that has not been assigned to either, and anything
promoted out of it joins the final test level under §6.5.

**Development level.** Free to inspect, tune against, and iterate on. Results from it are
never evidence of generalization.

- **TrpB** (four-site, 159,999 variants) — the landscape every configuration choice was made on.
- **GB1** (four-site, 149,360 variants) — selected as a second landscape and then used for
  mechanism analysis; it is development data, not a held-out check.
- **HER2 trastuzumab CDR-H3** — **designated development data.** See §4.2.

**Final test level.** Protein families never used to select an algorithm, a component, a
prior, a hyperparameter or a configuration. Declared in §5 before any arm is run.

**Reserve level.** Frozen targets held back entirely, to be promoted into the final test
level only under the reclassification rule of §4.4 or to raise statistical power (§6.5):
HIV-1 protease, SaCas9, SpCas9, RhlA. Nothing about them has been inspected beyond the
label-only statistics in §9.3, which are published here for every target precisely so that
computing them is not a hidden look.

### 4.2 HER2 is development data, not held-out — a correction

HER2 appeared on the candidate list for the independent set. It cannot be used that way.
Commit `38e98be` added the pool-restricted selection path to `select_next_micro_batch`
**because HER2 needed it**: a 10-site library with 11,277 labelled variants against a 20^10
product makes free generation propose variants nobody measured, so the batch came up short.
That same commit also decided, by looking at the HER2 file, that the 25,114 `test` rows
carrying `DMS_score` exactly 0.0 are absent labels and must be excluded.

Both are algorithm and curation changes made in response to inspecting HER2. Under the rule
in §4.4 — which this document would otherwise be applying to someone else — HER2 is
reclassified as development data. It stays in the study as a development target and its
results carry no weight in any generalization claim.

This also means the pool-restricted code path was **built on a development target** and will
be exercised on the independent ones. That is a confound and is registered as such in §7.3.

### 4.3 Development-level cross-check

GB1 appears in both this repository's prepared CSV and the frozen corpus. The two will be
compared once, for file-level agreement only (eligible universe size, fitness distribution,
depth profile), to establish that the corpus loader and the existing harness agree on a
landscape both can read. This is a plumbing check on development data. It is not an
endpoint and it is not reported as a result. D-3 above is the one known disagreement.

### 4.4 Reclassification rule

> If a final-test protein is inspected and that inspection is used to modify the algorithm,
> a component, a prior, a hyperparameter, a configuration, the candidate-generation path or
> the curation of any target, that protein is **reclassified as development data** from that
> moment. It may still be reported, labelled as development. A new independent protein
> family must then be reserved from §4.1's reserve level — or, if the reserve is exhausted,
> acquired — before any generalization claim is made again.

"Inspected" includes looking at an arm's performance on it. The rule has already fired once,
on HER2 (§4.2), before any cross-protein arm ran; it is not hypothetical.

Reclassification events are appended to §16 with a date and the change that caused them.

---

## 5. The independent target set, and an honest count

Declared before any performance on any of these targets has been examined.

### 5.1 Final test level

| landscape key | protein | protein family | phenotype | sites | eligible universe |
|---|---|---|---|---|---|
| `eqfp611_blue` | eqFP611 / mTagBFP2 | GFP-like fluorescent protein (Anthozoa) | blue fluorescence | 13 | 8,192 |
| `eqfp611_red` | eqFP611 / mTagBFP2 | GFP-like fluorescent protein (Anthozoa) | red fluorescence | 13 | 8,192 |
| `eqfp611_combined` | eqFP611 / mTagBFP2 | GFP-like fluorescent protein (Anthozoa) | combined | 13 | 8,192 |
| `creilov` | CreiLOV | LOV / flavin-binding photosensor | fluorescence | 15 | 165,428 |
| `cr9114_h1` | CR9114 | immunoglobulin heavy-chain variable domain | binding, H1 | 16 | 48,841 |
| `cr6261_h1` | CR6261 | immunoglobulin heavy-chain variable domain | binding, H1 | 11 | 1,887 |
| `cr6261_h9` | CR6261 | immunoglobulin heavy-chain variable domain | binding, H9 | 11 | 1,842 |

### 5.2 The count

Seven landscapes. That is the largest number in this section and the least meaningful one.

- **7 landscapes.**
- **4 distinct proteins**: eqFP611/mTagBFP2, CreiLOV, CR9114, CR6261.
- **3 independent protein families**: GFP-like fluorescent protein; LOV photosensor;
  immunoglobulin heavy-chain variable domain.

The reductions are not bookkeeping; each one removes a kind of evidence that is not there:

1. **CR6261 against H1 and H9 is one protein family assayed against two antigens.** It is
   one protein, one library, one set of 11 sites, two binding targets. Paired assays on one
   protein are not independent evidence of cross-protein generalization. The same applies
   to eqFP611's blue, red and combined channels: **one protein, one library of 8,192
   genotypes, three readouts.**
2. **CR9114 and CR6261 are both VH1-69-derived anti-influenza broadly neutralizing antibody
   heavy chains**, characterised in the same publication for that reason. Counting them as
   two independent families would overstate the breadth of the test. They are counted as one
   family, with two proteins inside it.
3. Trastuzumab's CDR-H3 library is also an immunoglobulin VH domain, so even if HER2 had
   remained eligible (§4.2) it would have added a protein to that family, not a family.

So the honest headline is: **three independent protein families, four proteins, seven
landscapes.** Any claim of generalization from this study is a claim about three families.
At n = 4 proteins the across-protein statistics in §6.4 and §12 have very little power, and
§6.5 sets out what expanding the set would buy.

For the primary per-protein analyses the four rows are: eqFP611/mTagBFP2 (3 landscapes),
CreiLOV (1), CR9114 (1), CR6261 (2). All four are present and characterised; nothing in the
final test level is outstanding.

---

## 6. Pre-registered primary hypothesis

### 6.1 The hypothesis

> **H1.** VeritasBio's reliability policy outperforms strong baselines on landscapes whose
> near-wild-type shell is enriched in top-tier hits, and underperforms on landscapes where
> it is not.

**Why it is a hypothesis and not a summary.** The mechanism analysis on TrpB and GB1 found
that the *policy* does not change between landscapes and the *location of the hits* barely
changes either; what changes is the richness of the shallow shell.

| | TrpB | GB1 |
|---|---|---|
| top-tier hits at depth 3–4 | ~90–95% | ~90–95% |
| mean depth of a hit | 3.38 | 3.59 |
| deep-shell (depth ≥ 3) hit rate | 0.91% | 0.97% |
| **shallow-shell (depth ≤ 2) hit rate** | **168/2,242 = 7.49%** (7.5x base rate) | **72/2,167 = 3.32%** (3.3x) |
| reliability's mean selection depth | 2.63 | 2.56 |
| baselines' mean selection depth | 3.48–3.60 | 3.48–3.60 |
| round-1 hits, reliability vs ridge | 10.44 vs 2.69 | 0.94 vs 3.35 |
| reliability vs ridge, raw top-1% | +13.43 | −22.49 |

The policy is a shallow exploiter on both targets. The payoff is not. TrpB's shallow shell
is 2.26x richer than GB1's, and that is the one structural difference found.

**Current support: n = 2 landscapes, one in each direction.** That is suggestive and
nothing more. Two points cannot distinguish this explanation from several others, and the
two landscapes that produced it are its training data (§6.6).

### 6.2 Why it can be pre-registered

Shallow-shell enrichment is computable from a fully-measured landscape using **labels and
the library definition only** — no model, no surrogate, no optimizer, no campaign. So it
can be declared as a covariate in advance rather than fitted to results afterwards. The
protocol therefore **requires** it, for every target, before any arm runs, published
alongside the ceilings of §9.3. `scripts/freeze_cross_protein_corpus.py stats` computes it.

### 6.3 Defining depth, and the shell, across libraries of different width

**Depth** = the number of library positions at which a variant differs from the wild-type
reference. Not the number of mutations in the whole protein, and not a sequence distance:
only library positions count. This definition is width-independent and applies unchanged to
a 4-site 20-way library, a 13-site binary library and HER2's 10-site CDR-H3 library.

**The shell cut cannot be a fixed absolute depth, and the reason is resolution, not taste.**
The candidate libraries are mostly *binary* — one alternative residue per site, not 19 — so
the depth ≤ 2 shell holds between 59 and 197 variants on them, against 2,242 on TrpB. At a
base rate near 1% such a shell yields an expected hit count of one or two, and the enrichment
estimate is unusable. Measured on the frozen corpus:

| landscape | sites | depth ≤ 2 shell size | hits in it | depth ≤ 2 enrichment (95% CI) | resolvable? |
|---|---|---|---|---|---|
| TrpB (dev) | 4 (20-way) | 2,242 | 168 | 7.49x | yes |
| GB1 (dev) | 4 (20-way) | 2,168 | 72 | 3.32x (2.64–4.16) | yes |
| eqFP611 blue | 13 (binary) | 92 | 14 | 15.20x (9.28–23.91) | no |
| eqFP611 combined | 13 (binary) | 92 | 14 | 15.20x (9.28–23.91) | no |
| eqFP611 red | 13 (binary) | 92 | 0 | 0.00x (0.00–4.00) | no |
| CreiLOV | 15 (restricted) | 197 | 18 | 9.13x (5.86–13.97) | no |
| CR9114 H1 | 16 (binary) | 134 | 0 | 0.00x (0.00–2.78) | no |
| CR6261 H1 | 11 (binary) | 59 | 0 | 0.00x (0.00–3.04) | no |
| CR6261 H9 | 11 (binary) | 66 | 0 | 0.00x (0.00–2.74) | no |

Four of the seven independent landscapes tie at exactly zero with intervals spanning the
whole plausible range. A rank test on that column is a test on ties.

**Pre-specified primary shell definition: depth ≤ floor(S/2), where S is the number of
library sites.** On a four-site library floor(4/2) = 2, so it **reproduces the development
definition exactly** on TrpB and GB1 without redefinition, and it is resolvable on every
independent target. This is the covariate the primary test of §6.4 uses.

**Pre-specified secondary definitions,** reported for every target, used as sensitivity
analyses and never as a substitute when the primary is unwelcome: absolute depth ≤ 2
(mechanistically faithful — the candidate generator proposes singles and pairs — but
unresolvable where the shell holds fewer than 500 variants), and depth ≤ round(S/4).

**Degeneracy guard, pre-specified.** The covariate is declared *undefined* for a target when
the shell cut does not discriminate: specifically when the cut is at or above the maximum
observed depth, or the shell contains more than 80% of the eligible universe. RhlA fails
this — its measured variants reach depth 3 only, so a floor(11/2) = 5 cut admits 100% of the
universe and returns an enrichment of exactly 1.00 by construction. **All three RhlA
landscapes are therefore excluded from the H1 test**, though they remain usable for every
other endpoint. HIV-1 protease is a borderline case — its floor(11/2) = 5 shell holds 72.5%
of the universe and all 59 of its hits — and is flagged rather than excluded; its enrichment
of 1.38x at the primary cut against 8.85x at depth ≤ 2 is the most severe cut-sensitivity in
the corpus and must be reported with both values.

**For HER2's 10-site library** the primary cut is floor(10/2) = 5. HER2 is development data
(§4.2) so this does not enter the test, but the number is fixed here so that it cannot be
chosen later.

### 6.4 The pre-specified test

**Statistic.** Spearman rank correlation `rho_s` between a landscape's primary shallow-shell
enrichment (§6.3) and the reliability arm's effect size versus ridge on the primary endpoint
at the 125-measurement budget, computed one landscape per row. Per-landscape effect sizes
are the seed-clustered paired mean contrasts of §11.

**Unit of the test.** Two versions, both pre-specified and both reported:

- **Per-protein (primary).** One row per *protein*, so that eqFP611's three channels and
  CR6261's two antigens do not inflate n. Where a protein contributes several landscapes,
  the row is the mean of its landscapes' enrichments and the mean of its effect sizes, with
  the within-protein spread reported. **n = 4 proteins.**
- **Per-landscape (secondary, explicitly non-independent).** n = 7. Reported because it is
  more sensitive and labelled as pseudo-replicated because it is.

**p-value.** Exact permutation over all orderings, not an asymptotic approximation. At n = 4
there are 24 orderings, so the smallest attainable one-sided p is 1/24 = 0.042: a *perfect*
ordering is only just significant and anything short of perfect is not. This is stated in
advance: *at the current corpus size the primary test is underpowered by construction, and a
non-significant result is uninformative about H1.* The per-landscape version is not a repair
for that; it is the same four proteins counted seven times.

**Alpha.** 0.05, one-sided, for the single primary test. Predicted direction: positive.

**Falsification.** H1 is treated as disconfirmed if any of the following holds:

1. `rho_s` ≤ 0 on the per-protein primary test;
2. a landscape in the upper half of the enrichment ordering shows a material regression
   (§12.3) against ridge;
3. a landscape in the lower half shows an improvement exceeding the minimum scientifically
   meaningful threshold (§12.3);
4. the primary and the depth ≤ 2 secondary cut disagree in the *sign* of `rho_s` — in which
   case H1 is recorded as **unresolved under the shell definition**, not as supported. The
   favourable definition is not selected afterwards.

### 6.5 The pre-registered ordering, and the prediction it forces

The covariate is computed and frozen here. The independent landscapes, ranked by primary
shallow-shell enrichment, with the development anchors shown for scale:

| rank | landscape | level | shell cut | shell size | hits | shell hit rate | enrichment (95% CI) |
|---|---|---|---|---|---|---|---|
| — | TrpB | development | ≤2 | 2,242 | 168 | 7.49% | 7.49x |
| — | GB1 | development | ≤2 | 2,168 | 72 | 3.32% | 3.32x (2.64–4.16) |
| 1 | `eqfp611_blue` | final test | ≤6 | 4,096 | 77 | 1.880% | **1.880x** (1.51–2.34) |
| 2 | `creilov` | final test | ≤7 | 62,791 | 1,178 | 1.876% | **1.876x** (1.77–1.98) |
| 3 | `eqfp611_combined` | final test | ≤6 | 4,096 | 76 | 1.856% | **1.856x** (1.48–2.31) |
| 4 | `eqfp611_red` | final test | ≤6 | 4,096 | 26 | 0.635% | **0.635x** (0.43–0.93) |
| 5 | `cr9114_h1` | final test | ≤8 | 32,445 | 202 | 0.623% | **0.623x** (0.54–0.71) |
| 6 | `cr6261_h1` | final test | ≤5 | 946 | 4 | 0.423% | **0.210x** (0.08–0.54) |
| 7 | `cr6261_h9` | final test | ≤5 | 918 | 0 | 0.000% | **0.000x** (0.00–0.21) |

**Ranks 1–3 are a near-tie and ranks 4–5 are a near-tie,** at 1.880/1.876/1.856 and
0.635/0.623, with overlapping intervals. The ordering within each cluster is not meaningfully
resolved and the rank test should be read as separating three groups — {blue, CreiLOV,
combined} > {red, CR9114} > {CR6261 H1} > {CR6261 H9} — rather than seven distinct levels.
This is stated now so that a `rho_s` driven by an unresolved within-cluster ordering is not
reported as a clean monotone result.

The per-protein rows the primary test uses, from the means of the above:

| protein | landscapes | mean primary enrichment |
|---|---|---|
| CreiLOV | 1 | 1.876x |
| eqFP611 / mTagBFP2 | 3 | 1.457x (spread 0.635–1.880) |
| CR9114 | 1 | 0.623x |
| CR6261 | 2 | 0.105x (spread 0.000–0.210) |

Four proteins, a strict ordering, no ties.

**The prediction this forces, stated plainly.** Every independent enrichment measured
(≤ 1.88x) sits **below GB1's 3.32x**, where reliability already lost to ridge by 22.49 raw
top-1% hits. H1 therefore predicts that **the current reliability configuration will lose to
ridge on every independent landscape in this corpus**, most severely on `cr6261_h9` and
`cr6261_h1`, and least severely on `eqfp611_blue`, `creilov` and `eqfp611_combined`. If
reliability instead wins broadly, H1 is wrong and so is the mechanism story behind it — and
that is a more interesting result than the one predicted.

Registering a prediction of across-the-board failure is the point of doing this before the
runs. It is also the result that would be hardest to publish and easiest to re-narrate
afterwards, which is why it is written down first.

Two features of this corpus are worth naming before the numbers arrive, because they make
the test sharper than n = 3 families usually allows:

- **`eqfp611_blue` against `eqfp611_red` holds the protein, the library and the 8,192
  genotypes fixed and varies only the shell structure** (blue: hit rate falls monotonically
  with depth, 30.8% at depth 1, 0.08% at depth 8; red: *zero* hits at depth 0–3, peaking at
  1.63% at depth 8). It is the cleanest available within-protein manipulation of the
  covariate, and it is not independent cross-protein evidence. Both statements are true at
  once.
- **CR9114 and CR6261 are affinity-maturation landscapes where hit rate rises monotonically
  with depth** (CR6261 H9: 0% up to depth 5, 22.9% at depth 9, 40.0% at depth 10). They are
  structurally the opposite of TrpB and are the hardest case for a shallow exploiter.

One caveat on the covariate itself, recorded because it cuts against the clean reading above:
**CreiLOV's depth profile is shallow-favouring but its primary shell cut hides that.** Its
hit rate falls from 9.66% at depth 2 to 0.12% at depth 12, yet floor(15/2) = 7 admits 62,791
variants — 38% of the universe, spanning the whole productive range — so the enrichment comes
out at 1.88x while the depth ≤ 2 cut gives 9.13x and the depth ≤ 4 cut 4.64x. The primary cut
is the resolvable one, not the most discriminating one, and on wide libraries it compresses
real structure. This is why the full depth profile is a required published object (§10) and
why §6.4 falsification condition 4 exists.

**Expanding the set.** Promoting the four reserve proteins (HIV-1 protease, SaCas9, SpCas9,
RhlA) would take the endpoint tables from 4 proteins to 8 and add three families — Cas9
endonuclease, retroviral aspartyl protease, rhamnosyltransferase — for six families in all.
Recommended, with three caveats that reduce what it buys the H1 test specifically:

- **SaCas9 and SpCas9 are two orthologs of one family**, so the family count goes to six,
  not seven.
- **RhlA is excluded from the H1 test** by the degeneracy guard (§6.3). It adds a protein to
  the endpoint tables only.
- **SpCas9's 597-variant universe** makes the absolute-budget protocol unusable (§8.3), and
  under the fractional schedule it yields a single budget point (§8.2), so it cannot carry
  budget sensitivity.

So the H1 test would go from n = 4 to **n = 7 proteins**, where a perfect ordering gives an
exact one-sided p of 1/5,040 instead of 1/24. That is the real gain, and it is large.

Their covariates are already computed and frozen alongside the others, so promoting them
later is not a fresh look: HIV-1 protease 1.38x, SaCas9 1.59x, SpCas9 2.59x, RhlA undefined.
Note that SpCas9's 2.59x is the highest independent value in the whole corpus and sits
between GB1's 3.32x and the rest — which is precisely why its 20.9% universe consumption
must not be allowed to carry an H1 result.

### 6.6 H1's training data

H1 was derived from TrpB and GB1. Those two landscapes are its training data and **cannot be
counted as evidence for it**. They appear in §6.5 as anchors so the independent values can be
read against a known scale, not as two of the test's observations. The test is run on the
independent families only.

---

## 7. Arms

All arms run on **identical measured candidate pools and identical starting samples**, per
seed, per target. Pairing is by seed: arm A and arm B on seed s see the same 25 initial
observations drawn from the same eligible universe.

| # | arm | description |
|---|---|---|
| A1 | random | uniform sampling from the eligible universe |
| A2 | ridge, greedy | one-hot features over library positions, refit each round, greedy selection, no diversity term |
| A3 | random forest, greedy | same features, `RandomForestRegressor`, refit each round, greedy |
| A4 | ALDE | the published implementation (`jsunn-y/ALDE`), where compatible — see §7.2 |
| A5 | VeritasBio **reliability** | the shipped default configuration (`panoptes_default`) |
| A6 | VeritasBio **frontier** | qualify-then-optimize (`panoptes_qualify_then_optimize`) |
| A7 | VeritasBio **lightweight** | one predefined low-cost configuration, fixed in §7.1 before any run |

Ridge (A2) is the designated reference baseline for every effect size, because it is the
simplest arm that is competitive and because the development results show it beating the
engine on distinct hits in one of the two regimes already.

### 7.1 The lightweight configuration must be fixed before running

A7 exists to answer "does the expensive machinery earn its cost", so it cannot be tuned
after seeing results. It is defined as: the reliability policy's candidate generation and
acquisition with the DKL surrogate, the GFlowNet enrichment and the LLM planner all
disabled, retaining the structural prior and the diversity (DPP/Pareto) stage. Its exact
flag set is recorded in the run manifest (§10) on first use and may not change afterwards;
if it must change, that is an amendment in §16 and the earlier runs are reported separately.

### 7.2 Incompatibility is documented, never substituted

Where an arm cannot be run on a target, the protocol **records the incompatibility and
reports the cell as absent**. Silently substituting a different method for ALDE, or a
different ALDE configuration, is prohibited; so is dropping the arm from the table without
a note. Each absent cell carries: the target, the arm, the exact failure, the probe output,
and whether the cause is the arm, the harness or the dataset.

ALDE is expected to be the common case. `probe_alde` already records the checkout, the
interpreter, the import result and any deviation from the pinned `alde.yml`
(`python=3.11`, `pytorch=2.1.1`, `gpytorch==1.11`, `botorch==0.9.4`). Additional expected
causes on this corpus: the binary per-site alphabets, and the 121- and 1,053-residue
genotype representations, which are not the saturated-site encoding ALDE's published code
assumes. If ALDE cannot run on a target, the comparison on that target is against A1–A3
and A5–A7 only, and the headline ALDE figure from §1.2 is not restated in that context.

### 7.3 Registered confound: the pool-restricted selection path

Every landscape in the frozen corpus has a restricted per-site alphabet — mostly a single
alternative residue per site — so free candidate generation would propose variants that were
never assayed. All independent targets must therefore run under the **pool-restricted**
selection path, in which the selector is handed the untested measured variants and ranks
them.

That path **skips GFlowNet enrichment** (sampling fresh trajectories would propose
out-of-pool variants), while scoring, Pareto and DPP still run. So the component set under
test on the independent targets is *not* the component set that produced the development
results on TrpB and GB1, which ran saturated and unrestricted.

Consequences, registered in advance:

1. Any cross-protein result is a result about the **pool-restricted** configuration. It must
   be described that way, including in the title of any table.
2. For comparability, the reliability and frontier arms are **also** run pool-restricted on
   TrpB and GB1, so that a development-to-independent comparison differs in the target and
   not in the code path. Those runs are development results.
3. The pool-restricted path was built in response to HER2, a development target (§4.2).

---

## 8. Budgets

Budget sensitivity is a **primary analysis**, not a robustness appendix, because the
development result is front-loaded: 6.3x the best baseline at 50 measurements (t = +25.9)
falling to 1.7x at 125 (t = +8.3), with the best *configuration* also changing with budget.
A single-budget headline would be a misreport.

### 8.1 Absolute budgets and starting observations

| total budget | starting observations | selected rounds | per-round batch |
|---|---|---|---|
| 25 | 25 | 0 | — |
| 50 | 25 | 1 | 25 |
| 75 | 25 | 2 | 25 |
| 100 | 25 | 3 | 25 |
| 125 | 25 | 4 | 25 |

The five budgets are **prefixes of one campaign per seed**, not five separate campaigns:
the 50-measurement figure is the state of the 125-measurement campaign after round 1. This
is what makes the paired per-seed contrasts of §11 valid across budgets.

**Budget 25 is an identity check, not a comparison.** With 25 shared initial observations
and zero selected measurements, every arm has measured exactly the same 25 variants. The
contrast between any two arms at budget 25 must be **exactly zero on every endpoint for
every seed**. A non-zero contrast there is an implementation leak — shared-start divergence,
seed misalignment, or an arm peeking — and the run is void, not interesting. This check runs
first on every target.

### 8.2 Fixed-fraction sensitivity analysis

Required for every target, and load-bearing on the small ones: a 125-measurement budget
consumes a very different share of the eligible universe from target to target.

| landscape | eligible universe | 125 as a share |
|---|---|---|
| `creilov` | 165,428 | 0.08% |
| `cr9114_h1` | 48,841 | 0.26% |
| `eqfp611_*` | 8,192 | 1.53% |
| `cr6261_h1` | 1,887 | 6.62% |
| `cr6261_h9` | 1,842 | 6.79% |
| `sacas9` (reserve) | 1,296 | 9.65% |
| `spcas9` (reserve) | 597 | **20.94%** |

A result that holds only where the campaign has consumed a fifth of the eligible universe is
a result about exhausting a library, not about search. The fixed-fraction analysis therefore
runs budgets at **0.25%, 0.5%, 1%, 2% and 4% of the eligible universe**, with the 1:4
structure preserved: starting observations = round(B/5), then four equal rounds of
round(B/5), with any remainder added to the final round.

Two pre-specified bounds, with their consequences:

- **Floor.** Fractional budgets below 15 total measurements are not run; the cell is reported
  absent. This removes 0.25% and 0.5% on both CR6261 landscapes and SaCas9, and 0.25%, 0.5%,
  1% and 2% on SpCas9 — **leaving SpCas9 a single fractional budget point (4% = 24
  measurements)**, which is enough for a level comparison and not enough for budget
  sensitivity. Said here rather than discovered later.
- **Cap.** Fractional budgets are capped at **625 total measurements** (five times the
  absolute protocol) for tractability. The cap binds on `creilov` above 0.25% and on
  `cr9114_h1` above 1%; where it binds, the cell is reported as capped with the fraction
  actually reached, never as the nominal fraction. The small libraries — which are the reason
  this analysis exists — are unaffected by the cap.

### 8.3 Targets where the absolute protocol is not usable

Pre-specified, from the universe sizes above rather than from any result:

- **`spcas9` (reserve, 597 variants).** 125 measurements is 20.9% of the universe. The
  absolute-budget protocol is **not** run; only the fixed-fraction analysis is, and the
  target is excluded from the equal-protein-weighted aggregate of §12.2.
- **`cr6261_h1` and `cr6261_h9` (1,887 and 1,842).** Run under both schedules. The absolute
  results are reported with the share-of-universe figure attached in the same table.
- **`sacas9` (reserve, 1,296).** Same treatment as CR6261, at 9.65%.

---

## 9. Endpoints

### 9.1 Tier selection, fixed before running

The "top 1%" tier does not survive contact with small libraries: 1% of CR6261 H1's 1,887
variants is 19 members, and a diversity-adjusted count over 19 hits has no dynamic range.
The tier is therefore selected per target by a rule fixed here, never per result:

> Walk the ladder 1%, 2%, 5%, 10% and take the **first tier holding at least 30 members.**

Resulting assignment on the frozen corpus:

| landscape | eligible | tier | tier members | boundary ties |
|---|---|---|---|---|
| `eqfp611_blue` / `_red` / `_combined` | 8,192 | 1% | 82 | 1 |
| `cr9114_h1` | 48,841 | 1% | 489 | 1 |
| `cr6261_h1` | 1,887 | **2%** | 38 | 1 |
| `cr6261_h9` | 1,842 | **2%** | 37 | 1 |
| `creilov` | 165,428 | 1% | 1,655 | 1 |
| `hiv1_protease` (reserve) | 1,175 | **5%** | 59 | 1 |
| `rhla_*` (reserve) | 2,346 | **2%** | 47 | 1 |
| `sacas9` (reserve) | 1,296 | **5%** | 65 | 1 |
| `spcas9` (reserve) | 597 | **5%** | 30 | 1 |
| TrpB, GB1 (development) | 159,999 / 149,360 | 1% | 1,600 / 1,494 | — |

Because the tier differs across targets, **raw counts are not comparable across targets**
and the cross-protein aggregate of §12.2 uses yield as a share of the achievable ceiling,
never a raw count. The tier used is printed in every table.

Boundary ties: where the tier boundary falls inside a tie, both the strictly-above count and
the tie-inclusive count are reported. The tie mass is 1 on every target characterised so
far, so this is currently a formality; it is fixed in advance because a ceiling-censored
target would make it decisive.

### 9.2 The endpoint list

**Primary.** Diversity-adjusted yield at d ≥ 2 (§2).

**Mandatory secondary — all of these, for every arm, target and budget:**

1. raw top-tier count;
2. best-so-far fitness against budget (the full curve, not its endpoint);
3. probability of reaching the top 0.1% and top 0.01%, **where resolvable** — see §9.5;
4. measurements to first qualifying hit (right-censored at the budget; censored campaigns
   reported as censored, never dropped and never imputed);
5. normalized improvement over random, `(arm − random) / (ceiling − random)`, per seed;
6. yield as a share of the **achievable ceiling** (§9.3), for both the primary and the raw
   endpoint — this is the only cross-protein-comparable quantity;
7. worst-case and lower-quantile campaign performance: minimum, 10th and 25th percentile
   over seeds, because a policy that is good on average and catastrophic in a tenth of
   campaigns is not usable;
8. candidate-generation and acquisition runtime, per round, wall-clock and per-measurement;
9. redundancy and cluster count, matching the development instrumentation (2.77x/4.52 for
   reliability, 7.84x for frontier), so that a raw-count win driven by near-duplicates is
   visible in the same table.

### 9.3 Ceilings are computed before anything runs

**Raw counts compared without ceilings are meaningless.** For every target, the maximum
achievable top-tier count and the diversity-adjusted ceiling at d ≥ 2 are computed from
labels alone and published here, before any arm runs.

| landscape | tier | raw hit ceiling | d≥2 ceiling | d≥3 ceiling | mean depth of a hit | max depth |
|---|---|---|---|---|---|---|
| `eqfp611_blue` | 1% | 82 | **38** | 14 | 3.96 | 13 |
| `eqfp611_red` | 1% | 82 | **28** | 10 | 7.11 | 13 |
| `eqfp611_combined` | 1% | 82 | **38** | 15 | 3.98 | 13 |
| `cr9114_h1` | 1% | 489 | **280** | 118 | 8.91 | 15 |
| `cr6261_h1` | 2% | 38 | **18** | 9 | 7.18 | 11 |
| `cr6261_h9` | 2% | 37 | **21** | 8 | 8.30 | 11 |
| `creilov` | 1% | 1,655 | **647** | 188 | 6.50 | 15 |
| `hiv1_protease` (res.) | 5% | 59 | 22 | 7 | 2.24 | 9 |
| `rhla_rls` (res.) | 2% | 47 | 23 | 7 | 1.89 | 3 |
| `rhla_rha` (res.) | 2% | 47 | 27 | 10 | 2.81 | 3 |
| `rhla_pct_rha` (res.) | 2% | 47 | 22 | 7 | 2.89 | 3 |
| `sacas9` (res.) | 5% | 65 | 35 | 17 | 3.66 | 8 |
| `spcas9` (res.) | 5% | 30 | 15 | 9 | 3.20 | 8 |
| GB1 (dev) | 1% | 1,494 | 146 | 24 | 3.59 | 4 |

Achievable ceiling at a budget B is `min(B, ceiling)`. Two consequences are already visible
and must be carried into interpretation:

- **On CR6261 the primary endpoint saturates below the budget.** The d ≥ 2 ceiling is 18
  (H1) and 21 (H9), so a 125-measurement campaign cannot score above 18 or 21 no matter how
  good it is. The dynamic range of the primary endpoint there is roughly 0–18, and ceiling
  share is the only sensible reading.
- **The d ≥ 2 filter bites far harder on four-site 20-way libraries than on wide binary
  ones.** It keeps 146 of GB1's 1,494 hits (9.8%) but 38 of eqFP611's 82 (46%), 647 of
  CreiLOV's 1,655 (39%) and 280 of CR9114's 489 (57%). Absolute diversity-adjusted counts are
  therefore **not comparable across library geometries**, which is the second reason §12.2
  aggregates ceiling shares.

### 9.4 Censoring, measured per target

| landscape | floor / min | mass at min | ceiling / max | mass at max | reading |
|---|---|---|---|---|---|
| `eqfp611_blue` | 0.0909 | 0.012% | 1.608 | 0.012% | effectively uncensored |
| `eqfp611_red` | 0.0252 | 0.012% | 1.692 | 0.012% | effectively uncensored |
| `creilov` | 620.3 | 0.0006% | 15,690 | 0.0006% | uncensored; 165,428 distinct values in 165,428 rows |
| `cr9114_h1` | 7.00 | **1.91%** | 9.835 | 0.002% | mild floor at the −log10(Kd) detection limit |
| `cr6261_h1` | 7.00 | **11.34%** | 9.729 | 0.053% | censored from below |
| `cr6261_h9` | 7.00 | **20.20%** | 10.15 | 0.054% | censored from below |
| GB1 (dev) | 0.00 | **19.74%** | 8.762 | 0.0007% | censored from below |
| HER2 (dev) | 0.10 | **51.7%** | — | — | heavily censored from below |

**All observed censoring in this corpus is from below.** That matters for which endpoints
survive:

- The **top-tier endpoints are unaffected** by a floor. Censoring at the bottom does not
  disturb which variants are in the top 1% or 2%.
- **Normalized improvement over random is affected**, because the random arm's mean is
  dragged toward the floor and the normalizer shrinks. On targets with ≥10% mass at the
  floor (`cr6261_h1`, `cr6261_h9`, GB1, HER2) this endpoint is reported with the floor mass
  printed beside it, and the alternative normalizer `(arm − random) / (ceiling − random)`
  computed over **ranks rather than values** is reported as well.
- **Surrogate fit is affected.** A ridge or random-forest regression against 20% of targets
  pinned at a detection limit is fitting a censored response. This is a property of the
  target, not a defect of the arm, and it is not corrected for; it is noted so that a poor
  baseline on CR6261 H9 is not misread as the engine doing well.
- **Pre-specified alternative endpoints for heavily censored targets** (floor mass ≥ 25%,
  which covers HER2 at 51.7% and no independent target at present): the primary endpoint
  becomes diversity-adjusted yield restricted to variants **strictly above the floor**, the
  raw secondary becomes the count strictly above the floor, and a third endpoint is added —
  the fraction of measurements spent on floor-valued variants, which is the cost a censored
  assay actually imposes. These are declared now so HER2 results are not re-specified later.

### 9.5 Where an endpoint cannot be supported

Said plainly rather than reported as a zero:

- **Top 0.1% and top 0.01% are not resolvable on the small targets.** 0.1% of CR6261 H1 is
  2 variants and 0.01% is 0. These endpoints are computed only where the tier holds at least
  10 members — on this corpus that is `cr9114_h1` (49 at 0.1%, 5 at 0.01%, the latter at the
  limit), `creilov` (166 and 17), and the development landscapes. Elsewhere the cell reads
  **"not resolvable at this universe size"**, never 0.0. The `eqfp611_*` landscapes have 9
  variants at 0.1% and 1 at 0.01%, so both are unresolvable there; at 0.01% the endpoint
  would be "did the campaign find the single global maximum", which is reported under
  best-so-far instead.
- **The shallow-shell covariate is undefined on all three RhlA landscapes** (§6.3 degeneracy
  guard): their measured variants reach depth 3 only, so no cut discriminates. Every other
  endpoint is computable there.
- **Statistical power for the H1 test is inadequate at n = 4 proteins** (§6.4). Not a missing
  endpoint, but the same category of honesty.

---

## 10. Data capture requirements

The analyses above and the intervention in §13 cannot be run from summary rows. The
benchmark must persist, per run, a manifest and a per-round table. If these rows are not
captured, testing H1 or the §13 intervention later requires repeating the entire study.

**Run manifest (one per arm × target × seed).** Corpus commit; dataset path and both
hashes; eligible universe size and the drop count; library key, site list and wild-type
reference; tier used and tier size; the ceilings of §9.3; arm name and full flag set
(including A7's, §7.1); pool-restricted yes/no; seed; initial-sample variant ids; code
commit; wall-clock and per-stage runtime.

**Per-round, per-depth rows** — the extension being built into the coverage module in
parallel, which this protocol depends on:

| field | why |
|---|---|
| round index, depth | the keys |
| candidates **proposed** at that depth | separates "the generator never offered it" from "acquisition declined it" — the distinction §13 turns on |
| candidates **selected** at that depth | gives the mean selection depth (2.63 / 2.56 on the development landscapes) |
| hits **returned** at that depth | gives per-depth realised return, which is what a depth-adaptive policy would act on |
| mean and max **revealed fitness** at that depth | the signal a sequential allocator would see |
| cumulative distinct-hit count at d ≥ 2 | the primary endpoint's trajectory, not just its endpoint |

Rows are written for every arm, including the baselines, because the baselines' selection
depth (3.48–3.60) is the comparison that makes the reliability arm's 2.6 meaningful.

---

## 11. Statistical analysis

### 11.1 The unit of independence is the seed

**The seed is the unit, not the observation, not the variant, not the round, not the
background pair.** A real error of this kind was found and fixed on 2026-10-31: 768
background-pair differences were pooled as independent when each of 48 seeds had contributed
16 of them. That inflated t-statistics by **1.3x to 2.6x**. Nothing in this protocol may
pool within-seed sub-units as independent.

### 11.2 Procedure

1. **Per-seed contrasts.** For each (target, budget, endpoint, arm pair), compute one
   difference per seed: `d_s = Y_s(arm) − Y_s(ridge)`, with both arms on the same seed's
   shared initial observations. The seed-level contrast vector is the only input to
   inference.
2. **Paired test over seeds.** Paired t-test on `{d_s}`, df = n_seeds − 1, reported with the
   mean contrast, its standard error, a 95% CI and n_seeds. Wilcoxon signed-rank is reported
   alongside wherever the contrast distribution is visibly skewed or the endpoint is a small
   count — which is most of this corpus (§9.3).
3. **Cluster bootstrap.** Resample **whole seeds** with replacement, B = 10,000, recompute
   the mean contrast, report the percentile interval. Never resample within a seed.
4. **Across targets.** The cluster is (target, seed) nested in target. The aggregate is
   computed with **equal weight per protein** (not per landscape, not per seed), and the
   bootstrap is stratified: resample seeds within each protein, and resample proteins as a
   second level where the across-protein interval is reported. Seeds from different targets
   are never treated as exchangeable.
5. **Seed count.** 48 paired seeds per arm × target × budget, matching the development runs.
   Any target run with fewer carries the count in every table.
6. **Reporting.** Every interval states what was resampled. Every t-statistic states its df
   and n_seeds. A figure without n_seeds attached is not reportable.

### 11.3 Multiplicity

The full grid is 7 arms × 5 budgets × 9-plus endpoints × 7 landscapes. Treating that as one
family would guarantee both false positives and, with Bonferroni, no power at all. The
pre-specified split:

- **Primary family.** `{reliability vs ridge}` × `{5 budgets}` × `{primary endpoint}`, per
  protein. That is 5 tests per protein and 20 across four proteins. **Holm–Bonferroni** at
  family-wise alpha 0.05, applied **within each protein** (5 tests), with the per-protein
  adjusted p-values reported; and separately across the 20, reported as the stricter figure.
  Both are shown; neither is chosen afterwards.
- **The H1 test** (§6.4) is a single pre-specified test with its own alpha of 0.05,
  one-sided, outside the primary family.
- **Everything else is exploratory.** All other arm pairs, all secondary endpoints, all
  subgroup and sensitivity analyses are reported with unadjusted CIs and **no significance
  claims**. "Exploratory" appears in the table caption, not in a footnote.
- **No post-hoc endpoint promotion.** An endpoint that was secondary when this document was
  frozen stays secondary whatever it shows.
- **No post-hoc budget selection.** The headline is the whole budget curve. A single budget
  may be quoted only with the other four beside it.

---

## 12. Success criterion

### 12.1 What the criterion is not

It is **not** "beats ridge on 3 of 5 datasets". A majority-of-datasets rule counts a
0.3-hit win on one target and a 22-hit collapse on another as 1–1, which is the opposite of
how a practitioner would read it. It rewards trivial wins and conceals severe losses.

### 12.2 What is reported, in full, every time

1. **The full effect versus ridge for every independent protein** — mean per-seed contrast,
   95% seed-cluster CI, n_seeds, at every budget, on the primary and the raw endpoint. Every
   protein appears whether it helps the case or not.
2. **Effects versus the stronger baselines** — random forest, ALDE where compatible (§7.2),
   and the lightweight configuration A7. An improvement over random alone is not a result.
3. **The fraction of independent proteins improving**, with the count and the denominator
   written out (for example "2 of 4 proteins; 3 families"), never as a bare percentage.
4. **The distribution and magnitude of regressions** — every regression, its size, its CI,
   and the worst one named in the summary sentence. Not a count of regressions: their sizes.
5. **An equal-protein-weighted aggregate effect**, computed on **yield as a share of the
   achievable ceiling** (§9.3) because tiers and library geometries differ across targets
   (§9.1, §9.3). One protein, one weight.
6. **Uncertainty intervals on all of the above**, from the seed cluster bootstrap, and the
   pre-specified sensitivity analyses: both shell definitions (§6.3), both budget schedules
   (§8), rank-based and value-based normalizers on censored targets (§9.4), and
   pool-restricted development runs for comparability (§7.3).

### 12.3 Thresholds, fixed in advance

Set from the development magnitudes, which is legitimate because TrpB and GB1 are
development data. Both are defined on the primary endpoint, per protein, on the per-seed
paired contrast against ridge.

**Minimum scientifically meaningful improvement (MSMI).** All three must hold:

- mean contrast ≥ **+20%** of ridge's diversity-adjusted yield on that protein; **and**
- mean contrast ≥ **+2 distinct qualifying hits** in absolute terms; **and**
- the 95% seed-cluster CI excludes zero.

Rationale for the absolute floor: at a 125-measurement budget a one-hit difference is inside
the noise of a single plate and does not change a campaign decision.

**Material regression threshold (MRT).** Either condition suffices:

- mean contrast ≤ **−10%** of ridge's diversity-adjusted yield on that protein with a 95% CI
  excluding zero; **or**
- mean contrast ≤ **−2 distinct qualifying hits** with a 95% CI excluding zero.

**The asymmetry is deliberate.** A regression is material at half the relative magnitude
required for an improvement to count. A tool that sometimes loses badly is harder to deploy
than one that sometimes wins mildly, and the development data already contains a −22.49
regression. The asymmetry is stated so that it cannot later be read as a bug.

On targets where the primary endpoint saturates below the budget (CR6261, §9.3) the relative
thresholds are computed against ridge's ceiling share rather than its count, and the
absolute floor of 2 hits is retained unchanged.

### 12.4 How the gate is framed

The criterion is **guidance for continued development**. It is explicitly **not** proof of
universal superiority, and a pass must never be reported as one.

- **Continue** this development direction if: no independent protein shows a material
  regression, **and** at least one independent *family* shows an improvement reaching MSMI,
  **and** the equal-protein-weighted aggregate ceiling-share effect is positive with a CI
  excluding zero.
- **Narrow** it — reframe as a specialist tool for an identified landscape class (§13.5) —
  if improvements and regressions both occur and sort in the direction H1 predicts.
- **Redirect** it to §13 if *every* independent protein regresses **and** the regressions sort
  in the order H1 predicts. That is the outcome §6.5 actually predicts. It means the current
  configuration is not shippable as a general tool and that the mechanism is understood well
  enough to try one specific fix. It is a stop for the current configuration and a start for
  depth-adaptive allocation, under §13.3's ordering and §13.5's kill criterion. It is not a
  licence for a second fix if the first fails.
- **Stop** it, under §14, otherwise — in particular if every independent protein regresses
  and the regressions do **not** sort as H1 predicts, which would mean the losses are real
  and unexplained.

A pass licenses one sentence: *on three independent protein families, under a pool-restricted
adapted 125-measurement protocol, this configuration produced more distinct top-tier variants
than one-hot ridge, with these per-protein effects and these regressions.* Nothing broader
follows from n = 3 families.

---

## 13. The one planned intervention: depth-adaptive allocation

Declared in advance, before any cross-protein result exists, so that it is a prediction
rather than a post-hoc fix.

### 13.1 The intervention

**Depth-adaptive allocation.** A policy that estimates per-depth return from **its own
measurements** each round and allocates later batches toward the depth that is paying. No
labels, no landscape metadata, no prior knowledge of the target, no outcome information the
campaign has not itself revealed.

This is the single planned evidence-driven change to the policy. It follows from the §6.1
mechanism: the policy selects at mean depth ~2.6 on both development landscapes while the
hits sit at depth 3.4–3.6, so it is misaligned with the hit mass on both — and happens to be
rescued on TrpB by a shallow shell rich enough to pay anyway.

### 13.2 Predicted direction, and what would disconfirm it

**Prediction:** depth-adaptive allocation improves results **more** on landscapes with poor
shallow-shell enrichment, where the fixed-depth policy is most misaligned, and is **roughly
neutral** on landscapes like TrpB where shallow exploitation already pays.

**Both of these disconfirm the mechanism even if mean performance rises:**

- a **uniform** improvement across all landscapes — that would indicate something other than
  depth misalignment was fixed;
- an improvement **concentrated on rich-shell** landscapes — the opposite of the predicted
  gradient.

Stating this now is the point. "It got better on average" will not be accepted as support
for the mechanism.

### 13.3 When it may be implemented

**After the cross-protein evaluation of this document is complete and its results frozen.**
Not before, and not in parallel. It is then validated on the **independent families only**;
TrpB and GB1 generated the hypothesis and cannot test it. Implementing it earlier would
convert the independent targets into development data under §4.4 and destroy the only
held-out evidence available.

### 13.4 A one-shot upfront pre-screen is not viable — recorded so it is not re-proposed

Estimating a target's shallow-shell richness from a handful of upfront probe measurements is
badly underpowered. At 8 probe measurements in the shallow shell, the expected hit counts are
roughly **0.6** at TrpB's 7.5% rate and **0.26** at GB1's 3.3% rate. Distinguishing those two
distributions from 8 Bernoulli draws is not possible at any useful confidence, and a
125-measurement campaign cannot spare meaningfully more than that without destroying the
budget the intervention is supposed to spend better.

Therefore the intervention **must be sequential allocation by observed return across
rounds**, not an upfront estimate. The one-shot pre-screen is ruled out here, with its
arithmetic, so that it is not re-proposed as a fresh idea later.

### 13.5 Kill criterion

Abandon depth-adaptive allocation if any of the following holds. These are instances of §14
and are not softer than it.

1. It improves **only TrpB** after repeated tuning.
2. Its benefit **vanishes on the independent families**.
3. It fails to help **specifically on the poor-shallow-shell landscapes** its mechanism
   predicts it should help — this is the mechanism-specific kill condition, and it fires even
   if average performance improved (§13.2).
4. A **simpler fixed-depth heuristic matches it** — for example, a fixed allocation to depth
   3–4, or simply selecting at the depth where hits have so far been found, with no adaptive
   estimator at all.
5. It needs hidden outcome information to choose its schedule.

### 13.6 The uncomfortable alternative, named in advance

If depth allocation is what actually matters, then the engine's advantage may reduce to
**"allocate budget by observed per-depth returns"** — needing no DKL surrogate, no GFlowNet,
no DPP diversity stage and no structural prior.

This is not a speculative worry. The 2^5 joint ablation already points that way: **removing
`dkl` and `gflownet` improved top-1% yield on both development landscapes.** That result
arrived after `docs/NEXT_STEPS.md` §3.1 was written describing the 2^5 design as the
instrument still to be built; it is carried here as reported and is not re-derived by this
document. Section 7.3 adds that GFlowNet enrichment is skipped anyway under the
pool-restricted path every independent target requires, so the cross-protein study will not
even exercise it.

**Pre-registered position: that outcome is a legitimate and publishable result implying a
substantially simpler product. It is not a failure, and it is not to be explained away.** If
the data supports it, the write-up says that the measured advantage is attributable to budget
allocation across mutation depth, that the heavier machinery was not necessary to obtain it,
and that the product should shrink accordingly. Naming it here removes the option of
re-narrating the data later.

---

## 14. Stopping rules

A direction is abandoned, or deferred with its reason recorded, when any of the following
holds. These are commitments, not warnings; each one has a measurement behind it that
already exists in this repository or is specified above.

1. **It improves only TrpB after repeated tuning.** The candidate-generation work is already
   0 for 6 (`docs/NEXT_STEPS.md` §3.2); a seventh attempt on the same target is not evidence.
2. **Its benefit vanishes on independent proteins.**
3. **It loses to simpler approaches with no compensating benefit** — one-hot ridge, a greedy
   random forest, or the lightweight configuration A7.
4. **It needs hidden outcome information to select configurations.** If a configuration can
   only be chosen by looking at the answer, it cannot be chosen in a real campaign.
5. **It depends on unverifiable labels.** The TrpB imputation argument
   (`docs/NEXT_STEPS.md` §2.3) is sound inference but is not the upstream flags; a result
   resting on which rows were imputed is not a result until the flags exist.
6. **It causes repeated data-insensitive selection** — the same variants chosen regardless of
   what the measurements said. Detectable from the §10 per-round rows.
7. **Its advantage disappears under cluster-aware uncertainty analysis** (§11). The 1.3x–2.6x
   inflation found on 2026-10-31 is the reason this rule exists.
8. **It improves raw counts only by collecting redundant variants.** Detectable from the
   primary endpoint against the raw secondary, and from the redundancy and cluster-count
   columns (§9.2 item 9).
9. **Its cost outweighs its measured benefit** — runtime, GPU spend or LLM calls against
   measured distinct-hit gain.

### 14.1 Negative results stay in the repository

Every negative result stays, with a configuration reproducible from the repository: the arm's
full flag set, the corpus commit and hashes, the seed list and the harness commit. Deleting a
negative result, or leaving it only in a commit message, is how `docs/545.md` came to record a
48-campaign live run as "never measured" (`docs/NEXT_STEPS.md` §1, §2.2). Negative results are
listed in the results document under their own heading, not in an appendix.

This includes the possibility, registered in §6.5, that **every** independent landscape shows
a regression. If that is what the data says, that is what gets written down.

---

## 15. Order of operations

The sequence matters, because some of these steps destroy the value of later ones if run out
of order.

1. Freeze this document and `scripts/freeze_cross_protein_corpus.py`. **Commit before any
   cross-protein arm runs.** ← this step
2. `fetch`, `record`, `verify` the corpus; commit the lock file. **Done 2026-10-31, 20/20
   files verified.**
3. `stats` on every target; the covariates of §6.5 and the ceilings of §9.3 are published
   above. Label-only, no arms. **Done 2026-10-31.** Recompute TrpB's anchor with the same
   code path (§16 O-4).
4. Run the budget-25 identity check (§8.1) on every target. A failure voids the run.
5. Run the development-level plumbing cross-check on GB1 (§4.3).
6. Run all arms on all targets, both budget schedules, 48 seeds. Capture §10 rows.
7. Analyse per §11, report per §12, with no changes to §§2, 6, 8, 9, 11 or 12.
8. Only then consider §13.

---

## 16. Open items this document does not resolve

| # | item | consequence |
|---|---|---|
| O-1 | ~~`CreiLOV_clean.xlsx` not fully retrieved.~~ **Resolved on 2026-10-31** before the freeze: all 20 pinned files fetched, verified against the commit tree and characterised. Nothing in §§3.4, 5.1, 6.5, 9.1, 9.3 or 9.4 is outstanding. | none |
| O-2 | **RhlA's library is measured only to depth 3** of its 11 sites (2,346 of 663,552 variants, 0.35% saturation), so the shallow-shell covariate is undefined there under the §6.3 degeneracy guard and no cut can be made to discriminate. An earlier draft of this document recorded RhlA as blocked for the primary endpoint; that was wrong — its `genotype` column is fixed-width (D-6) and the endpoint computes. | RhlA is excluded from the H1 test only. Its other endpoints are reported. |
| O-3 | **ALDE compatibility on this corpus is unprobed.** The binary per-site alphabets and the 121-/1,053-residue genotype representations are not the saturated-site encoding its published code assumes. | Likely absent cells, to be documented per §7.2 rather than substituted. |
| O-4 | **TrpB's shallow-shell figures (168/2,242 = 7.49%) were not recomputed here.** They are carried from the mechanism analysis. GB1's were independently reproduced by `freeze_cross_protein_corpus.py stats` against the frozen corpus — 72/2,168 = 3.3210% and a mean hit depth of 3.59, against the reported 72/2,167 = 3.32% and 3.59 — which validates the method but not TrpB's number. | Recompute TrpB's with the same code path before the H1 test is run, so both anchors come from one implementation. |
| O-5 | **The upstream corpus has seven recorded discrepancies or inconsistencies** (§3.6, D-1 to D-7), two of them unexplained one-row differences and one — SpCas9, 948 claimed against 597 usable — a factor of 1.6. | The clean file wins where they disagree; D-5 means RhlA's provenance is unestablished. |
| O-6 | **The H1 test is underpowered at the present corpus size** (§6.4): with n = 4 proteins the smallest attainable exact one-sided p is 0.042, so only a perfect ordering can reach alpha. Three of the seven landscape-level enrichments are also a near-tie (§6.5). | A non-significant result is uninformative. Expanding to the reserve proteins (§6.5) is the only fix, and it has its own caveats. |
| O-7 | **Whether the pool-restricted path changes the engine's behaviour enough to invalidate comparison with the development results is unknown** (§7.3). The mitigation — running TrpB and GB1 pool-restricted too — has not been run. | Until it is, development-to-independent comparisons confound the target with the code path. |
| O-8 | **No wet-lab validation of anything here.** Every figure is a retrospective simulation over a fully-measured library. | Nothing in this study is evidence about a real campaign. |

---

## 17. Amendments

Appended, dated, with a reason. The original text above is not edited in place. Any amendment
written after cross-protein results have been examined must say so in its first sentence.

| date | section | change | reason |
|---|---|---|---|
| 2026-10-31 | — | initial freeze: corpus pinned and verified (20/20 files), covariates and ceilings computed for all 7 final-test and 6 reserve landscapes, H1 and the §13 intervention registered | — |

No cross-protein arm had been run and no cross-protein performance figure had been examined
when this table was written.
