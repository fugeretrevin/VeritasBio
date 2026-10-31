# Cross-protein landscape inventory

Acquisition, verification and label-only characterisation of the landscapes available for
the cross-protein generalization study. Every number here is model-free: it is computed
from labels and genotypes alone, before any arm is run, so it can be pre-registered.

- **Status:** 4 of 5 requested target proteins obtained. One requested landscape
  (CR9114 against H9) **does not exist** — see [What could not be obtained](#what-could-not-be-obtained).
- **Corpus pin:** CombinGym @ `c5dccd514a9acfe1089a5f1555faea98afd78019`, verified file by file.
- **Reproduce:** `scripts/data/combingym_corpus.lock.json` holds the pin and every hash.

---

## 1. Provenance and verification

### Source

| field | value |
|---|---|
| collection | CombinGym (14 landscapes over 9 proteins) |
| repository | `https://github.com/sitonglab/CombinGym` |
| pin | commit `c5dccd514a9acfe1089a5f1555faea98afd78019` |
| commit date | 2026-09-16T07:33:53Z |
| tree | `911f19c7ba82e112056f154d3e8b2073f70dc7d1` |
| paper | Chen *et al.*, "CombinGym: a benchmark platform for machine learning-assisted design of combinatorial protein variants", bioRxiv, DOI `10.64898/2026.03.24.714074` |
| retrieved | 2026-10-31 |

The repository publishes **no tags and no releases**, so a commit SHA is the strongest
available pin. It is not a floating download: the commit is immutable and the file
contents are pinned by hash below.

### How the pin was verified

The expected hashes were **not** taken from any list written by hand. They were read from
GitHub's own tree API for the pinned commit:

```
https://api.github.com/repos/sitonglab/CombinGym/git/trees/911f19c7ba82e112056f154d3e8b2073f70dc7d1?recursive=1
```

That listing returned `truncated: false` and 1,087 blobs. For each of the 20 files used
here, the git blob SHA-1 of the bytes that arrived on disk was recomputed
(`SHA1("blob <len>\0" + content)`, which is what `git hash-object` produces) and compared
to the tree's value, along with the byte size. **20/20 matched.** A SHA-256 of each file
is recorded as well, independent of git's framing.

Downloads were written to their own fresh directory (`downloads/combingym-c5dccd5/`,
gitignored), all analysis scripts were kept in a separate directory, every path was passed
as a command-line argument, and every interpreter invocation used `python -I`.

### Independent cross-checks

Three checks were run that could each have caught a silent parsing error, and all three
passed on numbers that were fixed before this work began:

1. **Reference landscape reproduction.** The implementation reproduces the pre-confirmed
   shallow-shell figures exactly: TrpB 168/2,242 = **7.4933 %** shallow and **0.9077 %**
   deep; GB1 72/2,167 = **3.3226 %** shallow and **0.9661 %** deep. Agreement to four
   decimal places on both landscapes.
2. **ProteinGym agreement.** ProteinGym's reference file reports GB1 `SPG1_STRSG_Wu_2016`
   as **149,360** mutants — exactly the eligible universe computed here — and its CreiLOV
   entry `PHOT_CHLRE_Chen_2023` as **165,407** multi-mutants, which is exactly the 165,427
   variants computed here minus the 20 depth-1 singles.
3. **Round-trip through two unrelated canonicalisation paths.** Statistics were computed
   once from the full protein-sequence column in the workbooks and again from the
   mutation-string column of the exported CSVs. Both paths agree exactly on every field
   including the diversity-adjusted ceilings, which is the quantity a truncation bug would
   corrupt.

### Registry cross-validation

The library positions were **not** trusted from any single source. For every landscape the
positions listed in CombinGym's `Data_summary.xlsx` were checked against the positions that
actually vary across the data table, and against the wild-type FASTA:

- the set of varying sequence indices equals the registry's position list, **exactly**, for
  all six new landscapes (so no library position is missing and none is spurious);
- the wild-type FASTA residue at every library position equals the registry's residue;
- the depth computed as Hamming-distance-to-parent equals the table's own `n_mut` column on
  **all 232,376 rows**, with zero mismatches.

For eqFP611 this is a closed loop: the registry's sites `D20, V45, L63, …, K231` spell
`DVLTFNSAAYNNK`, which is byte-identical to the `mutant` string on the table's own
`n_mut = 0` row.

### Note for the coordinator

An untracked file `scripts/freeze_cross_protein_corpus.py` exists in the main checkout from
a parallel session and pins the same commit. That pin was verified here independently and
is **correct**. Its `PINNED_FILES` list also matches the tree. The two efforts should be
reconciled rather than both landing; this document's lock file is
`scripts/data/combingym_corpus.lock.json`, that script's is
`scripts/data/cross_protein_corpus.lock.json`.

---

## 2. Honest protein-family count

This is the number that governs what the study can claim, so it is stated before any
performance figure.

**9 landscapes → 7 distinct proteins → 5 structural families.**

| family | proteins | landscapes |
|---|---|---|
| Tryptophan synthase β-subunit (PLP fold type II) | TrpB | TrpB |
| Protein G domain B1 | GB1 | GB1 |
| **Immunoglobulin heavy-chain variable domain** | trastuzumab, CR9114, CR6261 | HER2, CR9114-H1, CR6261-H1, CR6261-H9 |
| GFP-like fluorescent protein (11-strand β-barrel) | eqFP611 / mTagBFP2 | eqFP611-blue, eqFP611-red |
| LOV / flavin-binding photosensor (PAS fold) | CreiLOV | CreiLOV |

What this acquisition actually adds:

- **+2 new families**, not +4: GFP-like fluorescent protein and LOV photosensor. CR9114 and
  CR6261 are new *proteins* but they fall in the immunoglobulin V-domain family that the
  existing HER2/trastuzumab landscape already occupies.
- **CR6261-H1 and CR6261-H9 are one protein assayed against two antigens.** They share a
  genotype set, a wild-type reference and 11 library positions; only the antigen differs.
  They are not independent evidence of cross-protein generalization and must not be counted
  as two families or as two proteins.
- **eqFP611-blue and eqFP611-red are one protein with two phenotypes**, on a single shared
  8,191-variant genotype set. Same caveat.
- CR9114 and CR6261 are genuinely different antibodies (different germlines, different
  epitopes), so they are two proteins — but within one family.

Counting generously (distinct proteins) gives 7; counting by family gives 5; counting
*independent* families newly added by this work gives 2.

---

## 3. Comparison table

### 3a. Structure

| landscape | variants | sites | WT reference | alphabet/site | comb. product | saturation |
|---|---:|---:|---|---|---:|---:|
| TrpB | 159,999 | 4 | `VFVS` | 20 (uniform) | 160,000 | 99.999 % |
| GB1 | 149,360 | 4 | `VDGV` | 20 (uniform) | 160,000 | 93.35 % |
| HER2/trastuzumab | 11,276 | 10 | `WGGDGFYAMD` | 20,20,19,19,18,13,12,10,9,5 | 182,463,840,000 | 6.2e-06 % |
| eqFP611 blue | 8,191 | 13 | `DVLTFNSAAYNNK` | 2 (uniform, binary) | 8,192 | 99.99 % |
| eqFP611 red | 8,191 | 13 | `DVLTFNSAAYNNK` | 2 (uniform, binary) | 8,192 | 99.99 % |
| CreiLOV | 165,427 | 15 | `GLRTAGQRDIDRVVT` | 2 (12 sites), 3 (2 sites), 5 (1 site) | 184,320 | 89.75 % |
| CR9114-H1 | 48,840 | 16 | `FSSITANTKSTSSRYY` | 2 (uniform, binary) | 65,536 | 74.52 % |
| CR6261-H1 | 1,886 | 11 | `TSANQKSTSAL` | 2 (uniform, binary) | 2,048 | 92.09 % |
| CR6261-H9 | 1,841 | 11 | `TSANQKSTSAL` | 2 (uniform, binary) | 2,048 | 89.89 % |

The observed combinatorial product equals CombinGym's independently stated "Variant Space"
for every landscape, including CreiLOV's non-uniform 184,320 = 5 × 3 × 3 × 2¹².

**The new libraries are binary, not saturated.** Four of the six carry exactly 2 residues
per position: each site is wild-type or one specific alternative. This is a structural
difference from TrpB and GB1, where every site carries all 20 amino acids. It changes what
"combinatorial library" means for a selector — the candidate space is 2ⁿ, not 20ⁿ — and it
is why these libraries reach 75–100 % saturation with far fewer measurements.

**The eligible universe excludes the wild type.** This is the benchmark's own rule
(`benchmark_panoptes_isolated.py:542`, *"the parent itself is not a library member"*): a
record canonicalising to the all-parent genotype yields an empty variant id and is dropped.
Applying it is what makes GB1 149,360 rather than 149,361 and the GB1 shallow shell 2,167
rather than 2,168, and it is corroborated by ProteinGym reporting the same 149,360. All six
new landscapes ship an explicit depth-0 row, and all six have it excluded.

### 3b. Labels, completeness and censoring

| landscape | rows in file | label completeness | unlabelled split | floor | mass at floor | ceiling | mass at ceiling | distinct values |
|---|---:|---|---|---:|---:|---:|---:|---:|
| TrpB | 159,999 | complete | none | -0.164071 | 1 (0.0006 %) | 1.0 | 1 | 159,630 |
| GB1 | 149,361 | complete | none | **0.0** | **29,477 (19.7 %)** | 8.76197 | 1 | 28,138 |
| HER2/trastuzumab | 36,391 | **31.0 % labelled** | **yes — `stage` != `train`; 25,114 rows carry `DMS_score` exactly 0.0, an absent label, not a measured zero** | **0.1** | **5,832 (51.7 %)** | 15.0842 | 1 | 1,918 |
| eqFP611 blue | 8,192 | complete | none | 0.0908746 | 1 (0.012 %) | 1.60767 | 1 | 8,073 |
| eqFP611 red | 8,192 | complete | none | 0.025218 | 1 (0.012 %) | 1.69244 | 1 | 8,085 |
| CreiLOV | 165,428 | complete | none | 2.79259 | 1 (0.0006 %) | 4.19552 | 1 | 165,427 |
| CR9114-H1 | 48,841 | complete | none | **7.0** | **934 (1.91 %)** | 9.83451 | 1 | 47,907 |
| CR6261-H1 | 1,887 | complete | none | **7.0** | **213 (11.3 %)** | 9.72861 | 1 | 1,674 |
| CR6261-H9 | 1,842 | complete | none | **7.0** | **371 (20.2 %)** | 10.15 | 1 | 1,471 |

**Fitness direction: higher is better for all nine landscapes.** No sign flip is needed
anywhere.

Units and transformation:

| landscape | units |
|---|---|
| TrpB | normalised enzymatic activity, linear |
| GB1 | relative fitness (enrichment ratio), linear |
| HER2 | binding classification score, linear, floored at 0.1 |
| eqFP611 blue / red | log-scale fluorescence brightness, each channel separately |
| CreiLOV | log10 mean fluorescence intensity (3 replicates averaged in log space) |
| CR9114, CR6261 | −log10(Kd) binding affinity, i.e. already a log scale |

**Censoring, by landscape.**

- **Not censored:** CreiLOV (every one of 165,427 values is distinct — no ties at all),
  eqFP611 blue and red (single variant at each extreme), TrpB.
- **Censored from below at exactly 0.0:** GB1, 29,477 variants (19.7 %).
- **Censored from below at exactly 0.1:** HER2, 5,832 variants (51.7 %) — *and* separately
  burdened by the unlabelled split above.
- **Censored from below at exactly 7.0:** all three bnAb landscapes. This is the assay's
  detection limit on the −log10(Kd) scale (Kd = 100 nM). CR6261-H9 is the most heavily
  censored at 20.2 %, CR6261-H1 at 11.3 %, CR9114-H1 at 1.9 %.
- **No landscape is censored from above:** the maximum value is held by exactly one variant
  everywhere.

> **Caveat on the bnAb floor.** CombinGym's *clean* tables put non-binders at **7.0**, but
> the *raw* tables (`Data/DMS/Raw/bnAbs_*_raw.csv`) encode them at **6.0**. CombinGym's
> cleaning step therefore changes the censoring value, and the clean tables also drop rows
> (CR9114: 48,841 clean rows against a 65,536 variant space). Anything that compares against
> the published Phillips *et al.* numbers must use one convention consistently. The clean
> tables are what is used here.

---

## 4. The two pre-computed quantities

### 4a. Maximum achievable top-1 % hit count and diversity-adjusted ceiling

Rule: the hit set is the `ceil(0.01 × N)` highest-valued variants of the eligible universe.
The ceiling is the same greedy rule the protocol scores with, applied to the *whole* hit
set: walk hits in fitness-descending order and keep one only when it is at least `d` Hamming
from everything already kept, over the **full-length** library genotype. No method can report
more diversity-adjusted hits than this, at any budget.

| landscape | universe | **max top-1 % hits** | **d≥2** | **d≥3** | **d≥4** | 125 as % of universe | hits reachable in 125 |
|---|---:|---:|---:|---:|---:|---:|---:|
| TrpB | 159,999 | 1,600 | 166 | 20 | 4 | 0.078 % | 125 |
| GB1 | 149,360 | 1,494 | 146 | 24 | 5 | 0.084 % | 125 |
| HER2/trastuzumab | 11,276 | 113 | 113 | 111 | 94 | 1.11 % | 113 |
| eqFP611 blue | 8,191 | 82 | 39 | 14 | 6 | 1.53 % | 82 |
| eqFP611 red | 8,191 | 82 | 28 | 10 | 6 | 1.53 % | 82 |
| CreiLOV | 165,427 | 1,655 | 647 | 188 | 71 | 0.076 % | 125 |
| CR9114-H1 | 48,840 | 489 | 280 | 118 | 54 | 0.26 % | 125 |
| **CR6261-H1** | 1,886 | **19** | 10 | 6 | 4 | **6.63 %** | **19** |
| **CR6261-H9** | 1,841 | **19** | 13 | 4 | 2 | **6.79 %** | **19** |

Reading this table:

- **The d≥2 ceiling varies 65-fold across landscapes (10 to 647).** Raw hit counts compared
  across these landscapes without normalising by the ceiling are not interpretable, which is
  the reason these were computed in advance.
- **HER2's ceiling is nearly unconstrained** (113 of 113 at d≥2, 94 at d≥4) because a 10-site
  library sampled at 6×10⁻⁸ saturation has no two hits close together. The diversity
  adjustment does essentially nothing on HER2 and will not discriminate between arms there.
- **The 4-site saturated libraries are the opposite:** TrpB collapses 1,600 raw hits to 4 at
  d≥4, because with only 4 sites `d≥4` means "differs everywhere". On a 4-site library d=4 is
  a near-degenerate constraint and should be read with that in mind.
- **eqFP611 and CR6261 cannot use the full budget productively.** Their entire top-1 % tier
  holds 82 and 19 variants, fewer than the 125-measurement budget, so an oracle exhausts the
  hit set before the budget is spent and the endpoint saturates.

### 4b. Shallow-shell enrichment

Depth = number of library positions differing from the wild-type reference. Base rate is the
1 % tier rate by construction.

Columns 3–7 are the fixed `depth≤2` cut; columns 8–12 are the fractional
`depth≤⌊sites/2⌋` cut.

| landscape | sites | d2 shell | d2 hits | d2 rate | **d2 enrich** | d2 deep rate | frac cut | frac shell | frac hits | frac rate | **frac enrich** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TrpB | 4 | 2,242 | 168 | 7.4933 % | **7.49x** | 0.9077 % | ≤2 | 2,242 | 168 | 7.4933 % | **7.49x** |
| GB1 | 4 | 2,167 | 72 | 3.3226 % | **3.32x** | 0.9661 % | ≤2 | 2,167 | 72 | 3.3226 % | **3.32x** |
| HER2 | 10 | *9* | 0 | — | *n/e* | 1.003 % | ≤5 | 2,602 | 23 | 0.8839 % | **0.88x** |
| eqFP611 blue | 13 | *91* | 13 | 14.29 % | *14.3x* | 0.8519 % | ≤6 | 4,095 | 77 | 1.880 % | **1.88x** |
| eqFP611 red | 13 | *91* | 0 | 0 % | *0x* | 1.012 % | ≤6 | 4,095 | 26 | 0.6349 % | **0.63x** |
| CreiLOV | 15 | *196* | 18 | 9.184 % | *9.18x* | 0.9907 % | ≤7 | 62,790 | 1,178 | 1.876 % | **1.88x** |
| CR9114-H1 | 16 | *133* | 0 | 0 % | *0x* | 1.004 % | ≤8 | 32,444 | 202 | 0.6226 % | **0.62x** |
| CR6261-H1 | 11 | *58* | 0 | 0 % | *0x* | 1.039 % | ≤5 | 945 | 3 | 0.3175 % | **0.32x** |
| CR6261-H9 | 11 | *65* | 0 | 0 % | *0x* | 1.070 % | ≤5 | 917 | 0 | 0 % | **0x** |

*Italic* fixed-cut figures rest on shells of 9–196 variants and are reported for
completeness only; see the resolvability argument below. `n/e` = not estimable.

#### The headline finding: the covariate's sign is not shared

The brief's two reference landscapes both show strong positive shallow enrichment (7.5x,
3.3x) with near-equal deep shells (0.91 %, 0.97 %). **That pattern does not generalize.**
Among the six new landscapes:

- **CreiLOV behaves like TrpB** — 9.18x shallow enrichment, deep rate 0.99 %. Same sign,
  comparable magnitude.
- **eqFP611-blue behaves like TrpB, more strongly** — 14.3x.
- **eqFP611-red, CR9114-H1, CR6261-H1 and CR6261-H9 have *zero* hits in the shallow shell.**
  Their hit rate *rises monotonically with depth*. CR6261-H9 is the extreme case: 0 hits
  through depth 6, then 4/297, 5/149, 7/48, 2/10, and the single depth-11 variant is a hit.
- HER2 is essentially flat in depth (0.88x fractional).

The cause is mechanical and must be stated in any pre-registration that uses this covariate:
**for CR9114 and CR6261 the depth-0 reference is the germline antibody precursor, not a
functional parent.** The 11 and 16 library positions are the somatic mutations that separate
germline from the matured broadly-neutralizing antibody, so the optimum sits at *maximum*
depth, not minimum. CR6261's germline sits at the assay floor (its wild-type −log10(Kd) is
exactly 7.0 for both H1 and H9 — it does not detectably bind either antigen). eqFP611-red is
the same situation in a different guise: the parent `DVLTFNSAAYNNK` is the blue-optimised
mTagBFP2, so red fluorescence requires moving far from it.

Consequently shallow-shell enrichment is **not a comparable covariate across these nine
landscapes without conditioning on whether the depth-0 reference is a functional parent or a
non-functional ancestor.** Four of nine are the latter. Pooling them into one regression
would mix two different quantities.

#### Proposed width-appropriate definition

A fixed `depth≤2` cut does not transfer past 4 sites, for a reason that is quantitative
rather than aesthetic: **the shell becomes too small to estimate a ~1 % rate.**

| landscape | sites | fixed `d≤2` shell | shell as % of universe | fractional shell | as % of universe |
|---|---:|---:|---:|---:|---:|
| TrpB | 4 | 2,242 | 1.40 % | 2,242 | 1.40 % |
| GB1 | 4 | 2,167 | 1.45 % | 2,167 | 1.45 % |
| HER2 | 10 | **9** | 0.08 % | 2,602 | 23.1 % |
| eqFP611 | 13 | **91** | 1.11 % | 4,095 | 50.0 % |
| CreiLOV | 15 | **196** | 0.12 % | 62,790 | 38.0 % |
| CR9114-H1 | 16 | **133** | 0.27 % | 32,444 | 66.4 % |
| CR6261-H1 | 11 | **58** | 3.08 % | 945 | 50.1 % |
| CR6261-H9 | 11 | **65** | 3.53 % | 917 | 49.8 % |

**Recommendation: use `depth ≤ ⌊sites/2⌋` as the primary pre-registered definition, and
report the fixed `depth ≤ 2` figures alongside it.** Justification:

1. **It is exactly backward-compatible.** On a 4-site library `⌊4/2⌋ = 2`, so TrpB and GB1
   keep their pre-registered values *unchanged to the digit* — 7.4933 % / 7.49x and
   3.3226 % / 3.32x. No re-registration of the reference landscapes is needed.
2. **It is the only one of the two that is estimable.** Estimating a rate of order 1 %
   usefully needs on the order of 500+ shell members. Under the fixed cut, **not one** of the
   six new landscapes reaches 500 (9 to 196). Under the fractional cut every landscape clears
   it (917 to 62,790). HER2's fixed shell of 9 variants cannot distinguish a 0 % rate from a
   10 % rate at all.
3. **"Shallow" is a claim about the fraction of the designed library perturbed**, not an
   absolute mutation count. Two mutations out of 4 designed sites and 2 out of 16 are not the
   same experiment; 2-of-4 and 8-of-16 are.

Two honest caveats on this choice:

- **Fraction-of-sites does not preserve fraction-of-library, and the fixed cut is actually
  the better-matched one on that axis.** The fixed `d≤2` shell holds 0.08–3.5 % of each
  universe, broadly comparable across libraries; the fractional shell holds 23–66 %. The
  reason is alphabet size: on a binary library, `depth ≤ ⌊n/2⌋` captures about half of 2ⁿ by
  binomial symmetry, whereas on a 20-letter 4-site library `depth ≤ 2` captures only 1.4 %.
  The fractional cut therefore buys statistical power at the cost of no longer selecting a
  *rare* shell. Which property matters depends on what the covariate is for; both columns are
  reported so either can be used, and the full depth profiles below let any third cut be
  derived without re-reading the data.
- A 50 %-of-library shell is not "shallow" in any lab sense. For narrative purposes
  `depth ≤ ⌈sites/4⌉` is closer to the spirit, but it drops CR6261 to 149–162 members and
  becomes unresolvable again.

#### Full depth profiles (hits / variants at each depth)

Published so that any alternative cut is derivable without re-reading the corpus.

| landscape | depth profile |
|---|---|
| TrpB (4) | 1: 12/76 · 2: 156/2,166 · 3: 639/27,436 · 4: 793/130,321 |
| GB1 (4) | 1: 6/76 · 2: 66/2,091 · 3: 463/26,019 · 4: 959/121,174 |
| HER2 (10) | 1: 0/3 · 2: 0/6 · 3: 2/66 · 4: 7/591 · 5: 14/1,936 · 6: 23/3,118 · 7: 33/2,945 · 8: 27/1,743 · 9: 5/688 · 10: 2/180 |
| eqFP611 blue (13) | 1: 4/13 · 2: 9/78 · 3: 20/286 · 4: 18/715 · 5: 17/1,287 · 6: 9/1,716 · 7: 4/1,716 · 8: 1/1,287 · 9–13: 0/1,093 |
| eqFP611 red (13) | 1: 0/13 · 2: 0/78 · 3: 0/286 · 4: 2/715 · 5: 5/1,287 · 6: 19/1,716 · 7: 25/1,716 · 8: 21/1,287 · 9: 7/715 · 10: 3/286 · 11–13: 0/92 |
| CreiLOV (15) | 1: 1/20 · 2: 17/176 · 3: 57/978 · 4: 145/3,565 · 5: 274/9,603 · 6: 325/19,193 · 7: 359/29,255 · 8: 263/34,230 · 9: 127/30,925 · 10: 67/21,167 · 11: 14/11,022 · 12: 5/4,093 · 13: 1/1,031 · 14: 0/159 · 15: 0/10 |
| CR9114-H1 (16) | 1: 0/16 · 2: 0/117 · 3: 1/538 · 4: 3/1,706 · 5: 10/3,946 · 6: 38/6,936 · 7: 55/9,367 · 8: 95/9,818 · 9: 110/7,984 · 10: 80/5,001 · 11: 52/2,365 · 12: 32/819 · 13: 9/196 · 14: 4/29 · 15: 0/2 |
| CR6261-H1 (11) | 1: 0/9 · 2: 0/49 · 3: 0/149 · 4: 2/308 · 5: 1/430 · 6: 3/430 · 7: 4/302 · 8: 5/150 · 9: 4/48 · 10: 0/10 · 11: 0/1 |
| CR6261-H9 (11) | 1: 0/11 · 2: 0/54 · 3: 0/151 · 4: 0/290 · 5: 0/411 · 6: 0/419 · 7: 4/297 · 8: 5/149 · 9: 7/48 · 10: 2/10 · 11: 1/1 |

CR9114's maximum observed depth is 15, not 16: the fully-matured all-mutation genotype is
absent from the curated table.

---

## 5. What could not be obtained

| requested | status | reason |
|---|---|---|
| eqFP611 **or** mTagBFP2 | **obtained** | These are the **same dataset**, not two. CombinGym's registry names the protein `mTagBFP2` and the file `eqFP611_clean.xlsx`; it is the Poelwijk *et al.* 2019 13-site binary library spanning the blue and red parents. One landscape, two phenotype columns. |
| CreiLOV | **obtained** | 165,427 variants. Also independently present in ProteinGym as `PHOT_CHLRE_Chen_2023`, which is the full single+multi DMS; CombinGym's table is the combinatorial subset and is the right source here. |
| CR9114 vs H1 | **obtained** | 48,840 variants. |
| **CR9114 vs H9** | **DOES NOT EXIST** | Phillips *et al.* 2021 never assayed CR9114 against H9. The raw table's own column header is `h1_*, h3_*, fluB_*` — CR9114 was measured against **H1, H3 and influenza B**. H9 was measured for **CR6261 only** (`h1_*, h9_*`). The brief's target list conflates the two antibodies' antigen panels. This cannot be acquired from any source because the experiment was not performed. |
| CR6261 vs H1 | **obtained** | 1,886 variants. |
| CR6261 vs H9 | **obtained** | 1,841 variants. |
| CR9114 vs H3 / influenza B | *available, not prepared* | Present in `Data/DMS/Raw/bnAbs_CR9114_raw.csv` (columns `h3_mean`, `fluB_mean`) but **not curated into a clean table** by CombinGym. Obtainable with modest extra work. It would add antigens, **not** protein families, and the raw table uses the 6.0 floor rather than the clean tables' 7.0. Flagged as an option, deliberately not built. |

Also checked and **not** a source for these: **ProteinGym** (217 substitution assays, 69 with
multiple mutants) carries CreiLOV and GB1 but contains **no** CR9114, CR6261, eqFP611 /
mTagBFP2, TrpB or HER2/trastuzumab entry. CombinGym is the only source for the bnAbs and
eqFP611.

Nothing was unobtainable for access, licensing or network reasons. The single gap is an
experiment that was never run.

---

## 6. Datasets judged unsuitable for the 125-measurement protocol

The protocol is 25 initial + 4 rounds of 25 = 125 unique measurements.

### Unsuitable as-is: CR6261-H1 and CR6261-H9

| problem | CR6261-H1 | CR6261-H9 |
|---|---|---|
| universe | 1,886 | 1,841 |
| 125 measurements as a share of the library | **6.63 %** | **6.79 %** |
| max achievable top-1 % hits | **19** | **19** |
| d≥4 ceiling | 4 | **2** |
| shallow-shell hits (fixed and fractional) | 0 and 3 | **0 and 0** |

Three independent reasons:

1. **The budget exhausts a large share of the library.** 125 of 1,886 is 6.6 %, against
   0.08 % on TrpB — an 85-fold difference in how much of the landscape an arm gets to see.
   A fixed 125-measurement budget is not the same experiment on these two libraries, so
   results are not comparable to the reference landscapes. **These need a fixed-fraction-of-library
   sensitivity analysis** (e.g. hold the budget at a constant ~0.08 % of the universe, which
   is 1–2 measurements here and is clearly infeasible; or accept the fraction and report it
   as a separate stratum rather than pooled).
2. **The endpoint cannot resolve.** The entire top-1 % tier is 19 variants. A count out of 19
   carries Poisson noise of roughly ±4, and the budget (125) exceeds the whole hit set, so an
   oracle and a good method score identically. The d≥4 ceiling on H9 is **2** — two arms
   differing by one distinct hit span the full dynamic range.
3. **The covariate is degenerate.** CR6261-H9 has **zero** top-1 % hits under both the fixed
   and the fractional shallow-shell definition, so its enrichment is 0x under both and
   carries no information for a covariate regression.

Recommendation: retain CR6261 as a **qualitative** demonstration that the shallow-shell
covariate can invert, which is genuinely valuable, but **exclude it from any pooled
performance comparison** at a 125-measurement budget, or analyse it in a separate
fixed-fraction stratum.

### Borderline: eqFP611 blue and red

8,191 variants; 125 measurements = 1.53 % of the library; top-1 % tier = 82 variants, which
is **smaller than the budget**. An arm can in principle measure the entire hit set, so the
raw endpoint saturates and only the diversity-adjusted figure (39 and 28 at d≥2) retains
headroom. Usable, but the raw top-1 % count should not be the primary endpoint here.

### Suitable without reservation

**CreiLOV** (165,427 variants, 0.076 % budget share, 1,655-hit tier, 647 at d≥2, uncensored,
positive 9.18x shallow enrichment) is the strongest new landscape and the closest structural
analogue to TrpB and GB1 in every respect except alphabet. **CR9114-H1** (48,840 variants,
0.26 % budget share, 489-hit tier, 280 at d≥2) is also well-sized, with the caveat that its
depth structure is inverted.

### Pre-existing caveat worth re-flagging

**HER2** carries both a 51.7 % mass at the 0.1 floor and a 25,114-row unlabelled split whose
`DMS_score` is forced to exactly 0.0. Admitting the latter would put 69 % of the universe at
an invented value; it is correctly excluded by `labelled_stage="train"`. Its fixed `depth≤2`
shell holds **9 variants** and supports no shell estimate at all — the concrete instance of
the brief's concern about libraries wider than 4 sites.

---

## 7. Files written

### Prepared datasets — `scripts/data/`

Column format matches the three CSVs already present
(`mutant,mutated_sequence,DMS_score,DMS_score_bin,stage`).

| file | rows | source sheet / column |
|---|---:|---|
| `EQFP611_BLUE_Poelwijk_2019.csv` | 8,191 | `eqFP611_clean.xlsx` / `genodata` / `blue` |
| `EQFP611_RED_Poelwijk_2019.csv` | 8,191 | `eqFP611_clean.xlsx` / `genodata` / `red` |
| `CREILOV_CHLRE_Chen_2023.csv` | 165,427 | `CreiLOV_clean.xlsx` / `Sheet1` / `mean_log` |
| `CR9114_H1_Phillips_2021.csv` | 48,840 | `bnAbs_CR9114_H1_clean.xlsx` / `CR91114-h1` / `h1_mean` |
| `CR6261_H1_Phillips_2021.csv` | 1,886 | `bnAbs_CR6261_H1_clean.xlsx` / `CR6261-h1` / `h1_mean` |
| `CR6261_H9_Phillips_2021.csv` | 1,841 | `bnAbs_CR6261_H9_clean.xlsx` / `CR6261` / `h9_mean` |

`scripts/data/combingym_corpus.lock.json` holds the upstream pin, the git blob SHA-1 and
SHA-256 of all 20 source files, and the SHA-256 and row count of each prepared CSV.

Encoding decisions, stated rather than implied:

- `mutant` uses ProteinGym notation over the library positions **with wild-type positions
  omitted**, matching the existing CSVs and `_combo_from_notation`. So `D20G:V45A` on a
  13-site library is a 13-site genotype whose other 11 sites are wild-type.
- The wild-type row is **dropped** (parent is not a library member, per §3a).
- `stage` is `train` on every row: these landscapes are fully labelled and ship no unlabelled
  split, unlike HER2.
- `DMS_score_bin` is 1 when `DMS_score >=` the wild-type value, else 0. **This is a new rule,
  not an inherited one** — the three existing CSVs each use a different, undocumented cutoff
  (TrpB 0.408085, GB1 1.0, HER2 0.1), so there was no repo-wide convention to follow. The
  benchmark does not read this column. Note that on CR6261 the wild type sits at the assay
  floor (7.0), so `DMS_score_bin` is 1 for every row and carries no information there.

### Untrusted-download directories (gitignored)

`downloads/combingym-c5dccd5/` — the 20 pinned source files.
`downloads/proteingym-refs/` — ProteinGym's `DMS_substitutions.csv`, for the coverage check.
`downloads/combingym-raw-probe/` — first 500 bytes of the two raw bnAb tables, used only to
read their column headers and settle the CR9114/H9 question.

### Out of scope, deliberately untouched

No dataset adapter was written and the benchmark was not modified: `veritas/datasets/` and
`benchmark_panoptes_isolated.py` belong to other agents. This work is acquisition,
verification and characterisation only.
