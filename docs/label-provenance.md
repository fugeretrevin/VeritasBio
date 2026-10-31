# Label provenance: which fitness values in our corpus were measured, and which were not

**Status:** literature and provenance research, written 2026-10-31. Companion to
`docs/prior-art-depth-allocation.md`; this document closes open items **1** and **2** of that
document's §8, and resolves two label-provenance questions that block external scientific
claims.

**Scope.** Four items. Items 1 and 2 are label provenance for the TrpB and GB1 landscapes.
Items 3 and 4 are risks to the novelty claim of `docs/preregistration-cross-protein.md` §13.1.
No source code was modified in producing this document.

**Verification convention** (same markers as `docs/prior-art-depth-allocation.md`):

- **[full text]** — the relevant section was read.
- **[abstract/snippet]** — only an abstract, README or search extract was read.
- **[verified in data]** — the claim was checked by computation against the primary data file,
  not just read in prose.
- **[could not obtain]** — the source was not retrieved. Nothing is asserted about its contents.

---

## 0. The verdicts, first

### 0.1 Item 4 — the novelty verdict, stated first and unambiguously

> **The §0 novelty verdict of `docs/prior-art-depth-allocation.md` SURVIVES. Qiu et al. 2024
> is not prior art for depth-as-an-action-space.** The AAAI-24 full text was obtained and
> read. The bandit's arms are **individual protein sequences**, not depths or shells. The
> mutation-depth quantity `N` is a **fixed constant of the regret analysis** that defines the
> search region; it is not an algorithm input, not scheduled across rounds, and never updated
> from observed returns. No per-depth allocation appears anywhere in the paper.

Open item §8.1 and risk §7.7 of `docs/prior-art-depth-allocation.md` are **closed**. The
entry at §1.1(e) should be upgraded from *[abstract only]* to *[full text]* and its caveat
paragraph ("there remains a residual risk that the paper schedules `N` across rounds") should
be struck.

### 0.2 The other three, in one line each

| item | verdict |
|---|---|
| **1. TrpB imputed identifiers** | **FOUND, and recovered in machine-readable form.** All 871 are published. Our file carries all 871, values matching to float precision. `docs/trpb_imputed_variants_871.csv`. |
| **2. GB1 label provenance** | **CLEAN.** Our 149,361 rows are exactly Wu et al.'s measured set. Zero overlap with their 10,639 imputed values. The absent ~10,639 combinations are genuinely absent, not imputed. |
| **3. Phillips et al. eLife 71393** | **PARTLY ANTICIPATED, and the paper also contradicts us.** It does plot max affinity against somatic mutation count (Figure 5A,B). For CR9114 the deep-shell structure is prior knowledge. For **CR6261 it concludes the opposite of our measurement**, which is the single most valuable thing in this document. |

---

## 1. Item 1 — TrpB's imputed-variant identifiers

**The landscape.** `scripts/data/TRPB_THEMA_Johnston_2024.csv`, sites 183/184/227/228, wild
type VFVS, stamped `UNVERIFIED_INCLUDES_IMPUTED`.

**The source.**

Johnston K.E., Almhjell P.J., Watkins-Dulaney E.J., Liu G., Porter N.J., Yang J., Arnold F.H.
(2024). *A combinatorially complete epistatic fitness landscape in an enzyme active site.*
Proc Natl Acad Sci U S A 121(32):e2400439121.
DOI [10.1073/pnas.2400439121](https://doi.org/10.1073/pnas.2400439121). PMID 39074291.
PMCID PMC11317637. Preprint: bioRxiv
[10.1101/2024.06.23.600144](https://doi.org/10.1101/2024.06.23.600144).

**[full text — bioRxiv PDF read in full; PNAS version read via PMC11317637]**

### 1.1 The published measured/imputed split

The counts are published in the Results, verbatim from the bioRxiv full text:

> "From the selection performed on the library at these four sites, 159,129 variants (99.45%
> of the total library) had sufficient sequencing coverage for quantification in both
> replicate experiments."

and

> "The missing 871 fitness values were imputed for downstream analyses (Methods, Fitness score
> calculations, SI Appendix, Fig S16)."

The method is named in the Methods, verbatim:

> "To enable the analyses in this paper, fitness values for the missing variants in the 4-site
> landscape were imputed with the sklearn KNN imputer (Fig S16), and the imputed scores
> reported previously by the authors were used for GB1 (26)."

So the split is **159,129 measured / 871 imputed**, out of 160,000. Imputation is
`sklearn.impute.KNNImputer`. (The second clause of that sentence is the thread that resolves
Item 2 — see §2.)

The paper also tracks the imputed set internally rather than forgetting it, which is why the
identifiers exist: *"No imputed variants were used as starting points, but they were used as
intermediate variants when necessary"* (path analyses), and *"Local optima were determined by
looking at all non-imputed variants classified as active."*

### 1.2 The identifiers exist, and here they are

**The identifiers are published.** They are not in the SI Appendix PDF and not in ProteinGym
(see §1.4). They are a file inside the CaltechDATA deposit named in the paper's Data,
Materials, and Software Availability statement:

| field | value |
|---|---|
| deposit | CaltechDATA, DOI [10.22002/h5rah-5z170](https://doi.org/10.22002/h5rah-5z170) |
| landing page | https://data.caltech.edu/records/h5rah-5z170 |
| archive | `data.zip`, 3,321,330,601 bytes, 559 members |
| **member** | **`data/figure_data/4-site_imputed/20230828_KNN_imputed_TrpB.csv`** |
| member size | 38,055 bytes uncompressed |
| rows | **871** data rows + header |
| columns | `,AAs,AA1,AA2,AA3,AA4,fitness,active` |

Raw sequencing is separately at SRA accession
[PRJNA1127511](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1127511).

The deposit is 3.3 GB, but the server honours HTTP range requests (`accept-ranges: bytes`), so
the member was retrieved by reading the ZIP central directory from the tail and then
range-fetching 15,174 compressed bytes. No bulk download is needed. The reproduction recipe is
in §1.6.

First rows as published:

```
,AAs,AA1,AA2,AA3,AA4,fitness,active
0,MKKG,M,K,K,G,0.4415112029968438,False
1,FQMG,F,Q,M,G,0.2949203141193624,False
2,MQKG,M,Q,K,G,0.26097919334107833,False
```

`AAs` is the four-residue genotype in site order 183, 184, 227, 228. All 871 rows have
`active = False`.

**Machine-readable form for us:** `docs/trpb_imputed_variants_871.csv`, 871 rows, columns
`mutant, aa_code_183_184_227_228, imputed_fitness, label_source`, where `mutant` is in the
same `V183M:F184K:V227K:S228G` form our CSVs use (sites equal to wild type are omitted, as in
our files) and `label_source` is the constant `KNN_IMPUTED`.

### 1.3 Verification against our own file

**[verified in data]** Computed by rebuilding the full four-site genotype for every row of
`scripts/data/TRPB_THEMA_Johnston_2024.csv` from wild type VFVS plus its `mutant` tokens, then
intersecting with the published `AAs` set:

| check | result |
|---|---|
| published imputed variants | 871 |
| distinct genotypes in our file | 159,999 |
| published imputed variants present in our file | **871 / 871** |
| value agreement (`DMS_score` vs published `fitness`) | **871 / 871 within 1e-6** (max residual 5.0e-10, i.e. float round-trip only) |
| wild type VFVS present in our file | **no** |

So our file is **159,999 variants = 159,128 measured mutants + 871 imputed**, with the
wild-type genotype itself excluded. The arithmetic closes exactly: Johnston's 159,129 measured
variants minus the wild type is 159,128, plus 871 is 159,999.

**The imputed fraction of our TrpB file is 871 / 159,999 = 0.544%.**

Two consequences worth acting on:

1. **The `UNVERIFIED_INCLUDES_IMPUTED` stamp can be replaced by an exact row-level flag.** We
   no longer need to say "some values are imputed and we do not know which". We know which.
2. **The 871 imputed rows are currently spread across all three of our `stage` splits:
   685 `train`, 100 `valid`, 86 `test`** **[verified in data]**. Model-derived labels are
   therefore present in our evaluation sets. At 0.544% this is unlikely to move an aggregate,
   but it is a real and now-removable defect, and it is the kind of thing a reviewer asks
   about once they know the dataset has imputed values at all. Note that all 871 are
   `active = False` **[verified in data: 871/871]**,
   so they are concentrated in the low-fitness mass and cannot be top-1% hits — the top-tier
   endpoints are not contaminated. The exposure is to anything that uses the full
   distribution (correlation metrics, calibration, the null model).

### 1.4 ProteinGym does not curate this landscape — clean negative

**[verified in data]** `reference_files/DMS_substitutions.csv` at the tip of
[OATML-Markslab/ProteinGym](https://github.com/OATML-Markslab/ProteinGym) lists **217**
substitution DMS assays. Searching it for `TRPB`, `TrpB`, `Tm9D8`, `Phillips`, `CR9114` and
`CR6261` returns **nothing**.

So:

- **ProteinGym has no curation of the TrpB landscape.** There is no ProteinGym-side
  measured/imputed annotation to inherit, and no ProteinGym assay ID for it. Our filename
  `TRPB_THEMA_Johnston_2024.csv` follows the ProteinGym naming convention but is **not** a
  ProteinGym assay ID.
- **ProteinGym does not curate the Phillips antibody landscapes either.** Same conclusion for
  `CR6261_H1`, `CR6261_H9` and `CR9114_H1`.

### 1.5 What ALDE does with these values, which matters for our comparison

Yang J., Lal R.G., Bowden J.C., Astudillo R., Hameedi M.A., Kaur S., Hill M., Yue Y.,
Arnold F.H. (2025). *Active learning-assisted directed evolution.* Nat Commun 16:714.
DOI [10.1038/s41467-025-55987-8](https://doi.org/10.1038/s41467-025-55987-8). PMC11739421.

**[full text, via PMC]** ALDE uses GB1 and TrpB as its two in-silico landscapes and describes
them as four-residue combinatorial datasets for which *"fitnesses have been measured
experimentally for nearly all of the 20⁴ = 160,000 variants"*. **ALDE does not mention
imputation anywhere**, and does not say how it handles the 871. Data availability:
https://github.com/jsunn-y/ALDE and https://zenodo.org/records/12196802.

**For us:** the imputed values are propagated silently by the most prominent prior user of
this landscape. That is mildly reassuring for comparability — we and ALDE are evaluating on
the same 160,000-row object — and it is also a small, free, honest differentiator: we can
state the split and they do not. It does **not** license us to leave imputed labels in our
test split; it only means the baseline comparison is not invalidated by them.

### 1.6 Reproduction recipe

The member can be recovered without downloading 3.3 GB:

1. `GET https://data.caltech.edu/records/h5rah-5z170/files/data.zip?download=1` with
   `Range: bytes=<size-1MiB>-<size-1>` where size = 3,321,330,601. Parse the end-of-central-
   directory record (the archive is **not** ZIP64: 559 entries, central directory at offset
   3,321,249,346, length 81,233).
2. Range-fetch the central directory and locate the entry
   `data/figure_data/4-site_imputed/20230828_KNN_imputed_TrpB.csv` — local header offset
   **1,805,298,899**, compressed size **15,174**, uncompressed **38,055**, method 8 (deflate).
3. Range-fetch 30 bytes at the local header offset to read the name and extra-field lengths,
   then range-fetch the compressed payload and `zlib.decompress(..., -15)`.

Note the deposit also contains `data/figure_data/4-site_merged_replicates/20230827/four-site_simplified_AA_data.csv`
(8,494,938 bytes), which is the measured four-site amino-acid table. It was not needed, because
the imputed list is the complement and is published directly — but it is the file to pull if we
ever want to re-derive the measured values independently of our own CSV.

---

## 2. Item 2 — GB1 label provenance

**The landscape.** `scripts/data/GB1_STRSG_Wu_2016.csv`, 149,361 data rows, sites 39/40/41/54,
wild type VDGV.

**The source.**

Wu N.C., Dai L., Olson C.A., Lloyd-Smith J.O., Sun R. (2016). *Adaptation in protein fitness
landscapes is facilitated by indirect paths.* eLife 5:e16965.
DOI [10.7554/eLife.16965](https://doi.org/10.7554/eLife.16965).
https://elifesciences.org/articles/16965

**[full text]**

### 2.1 Wu et al. did impute — but published the imputed values in a separate file

This is the thing worth knowing, and Johnston et al. flagged it for us: *"the imputed scores
reported previously by the authors were used for GB1 (26)"*, where reference 26 is Wu et al.
2016. So GB1 **does** have published imputed values, and the question is only whether we are
carrying them.

From the Wu et al. full text:

- Measured: *"The fraction of all possible variants that passed this filter was 93.4% (149,361
  out of 160,000 all possible variants)."*
- Missing: *"The fitness values for 10,639 variants (6.6% of the entire sequence space) were
  not directly measured."* Causes are non-sampling in the library and an input-count filter
  (variants with input count below 10 were removed).
- Imputation method: **lasso-regularized regression on log fitness** with main effects, pairwise
  terms and three-way interactions among sites 39, 41 and 54 — 9,102 coefficients, λ = 1e-4 by
  10-fold cross-validated MSE, Matlab R2014b, trained on the 119,884 non-lethal measured
  variants, Pearson r = 0.93 on measured variants. *"We then used the fitted model to impute the
  fitness of the 10,639 missing variants and complete the entire fitness landscape."*
- **The two sets are published as two separate files:**
  - **Supplementary file 1** — *"The fitness of each profiled variant."*
    DOI [10.7554/eLife.16965.024](https://doi.org/10.7554/eLife.16965.024),
    `elife-16965-supp1-v4.xlsx`
  - **Supplementary file 2** — *"Imputed fitness values for missing variants."*
    DOI [10.7554/eLife.16965.025](https://doi.org/10.7554/eLife.16965.025),
    `elife-16965-supp2-v4.xlsx`

Raw sequencing: BioProject
[PRJNA278685](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278685). Scripts:
github.com/wchnicholas/ProteinGFourMutants.

### 2.2 Verification: our file is the measured set, exactly

Both supplementary files were downloaded and parsed, and compared against our CSV.

**[verified in data]**

| check | result |
|---|---|
| Supplementary file 1 sheet dimension | `A1:E149362` → **149,361** data rows |
| Supplementary file 1 columns | `Variants, HD, Count input, Count selected, Fitness` |
| Supplementary file 2 sheet dimension | `A1:B10640` → **10,639** data rows |
| Supplementary file 2 columns | `Variants, Imputed fitness` |
| measured ∩ imputed | **0** |
| measured ∪ imputed | **160,000** |
| distinct 4-site genotypes in our file | **149,361** |
| our file ∩ imputed (supp. 2) | **0** |
| our file ∩ measured (supp. 1) | **149,361** |
| in our file but not in supp. 1 | **0** |
| in supp. 1 but not in our file | **0** |

> **Our GB1 file is Wu et al.'s measured set, exactly and completely: every row is measured,
> no row is imputed, and not one of the 10,639 imputed variants is present.** The ~10,639
> absent combinations are genuinely absent. Nothing model-derived is in this landscape.

### 2.3 This also closes pre-registered discrepancy D-3

`docs/prior-art-depth-allocation.md` §1.1(c) records ODBO's GB1 count of 149,361 against our
149,360 and notes "the one-row gap is pre-registered as discrepancy D-3".

**There is no discrepancy. All three numbers are the same object counted two ways:**

| source | number | what it counts |
|---|---|---|
| Wu et al. Supplementary file 1 | 149,361 | measured **variants**, including the wild-type genotype VDGV |
| ODBO (arXiv:[2205.09548](https://arxiv.org/abs/2205.09548)) | 149,361 | same |
| ProteinGym `SPG1_STRSG_Wu_2016` | 149,360 | **mutants** — 76 singles + 149,284 multiples, excluding wild type |
| our `GB1_STRSG_Wu_2016.csv` | 149,361 rows | 149,360 mutants **plus** an explicit `WT` row |

**[verified in data]** The ProteinGym reference row for `SPG1_STRSG_Wu_2016` gives
`DMS_total_number_mutants = 149360`, `DMS_number_single_mutants = 76`,
`DMS_number_multiple_mutants = 149284`, `raw_DMS_phenotype_name = Fitness`,
`raw_DMS_mutant_column = Variants` — the column names of Supplementary file **1**, confirming
ProteinGym's chain also runs through the measured file only.

D-3 can be closed as a definitional artefact, not a data defect. Note also that our GB1 file is
the one landscape in the corpus whose ProteinGym assay ID differs from our filename: ProteinGym
calls it **`SPG1_STRSG_Wu_2016`**, not `GB1_STRSG_Wu_2016`.

### 2.4 One incidental corroboration

**[verified in data]** The CaltechDATA deposit for Johnston et al. contains
`data/figure_data/GB1_data/GB1_missing_data.xlsx` at **340,223 bytes** — byte-identical in size
to the `elife-16965-supp2-v4.xlsx` we downloaded (340,223 bytes). Johnston et al.'s GB1 imputed
values are Wu et al.'s Supplementary file 2, unmodified. The chain is clean end to end.

---

## 3. Item 3 — Phillips et al., eLife 10:e71393

Phillips A.M., Lawrence K.R., Moulana A., Dupic T., Chang J., Johnson M.S., Cvijovic I.,
Mora T., Walczak A.M., Desai M.M. (2021). *Binding affinity landscapes constrain the evolution
of broadly neutralizing anti-influenza antibodies.* eLife 10:e71393.
DOI [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393).

**[full text]** — obtained as `https://cdn.elifesciences.org/articles/71393/elife-71393-v2.pdf`
and text-extracted locally. **Note for future reviews: the two earlier failures were not the
paper's fault.** Fetching the eLife HTML through a summarizing fetch tool is refused by the
summarizer's content filter on this article. Downloading the PDF directly and extracting the
text locally works on the first attempt. Open item §8.2 of
`docs/prior-art-depth-allocation.md` is **closed**.

### 3.1 Antigens — our belief is confirmed exactly

> **Confirmed.** CR9114 was assayed against **H1, H3 and influenza B**. CR6261 was assayed
> against **H1 and H9 only**. There is **no CR9114/H9 experiment** in the paper.

Verbatim: *"For CR6261, we chose two divergent group 1 HA subtypes (H1 and H9 …), while for
CR9114, we chose the three highly divergent subtypes present in the vaccine (H1 from group 1,
H3 from group 2, and influenza B …)."* The figure source data confirm it: *"Source data 1.
CR9114 library -logK_D to H1, H3, and influenza B. Source data 2. CR6261 library -logK_D to H1
and H9."*

Library construction, for the record: CR9114 = **16 mutations, 65,536 variants**; CR6261 =
**11 mutations, 2,048 variants**. Both are heavy-chain-only scFv yeast display, assayed by
Tite-Seq. Reported `N_var` after filtering: **65,091** for CR9114 and **1,812** for CR6261.

### 3.2 The assay floor — confirmed, and it is exactly our 7.0

> **Confirmed.** Verbatim: *"Inferred affinities outside our titration boundaries (10⁻¹¹–10⁻⁶ M
> for H3 and influenza B, 10⁻¹²–10⁻⁷ M for H1 and H9) are pinned to the boundary."* Methods
> repeat it: *"Inferred K_D outside of the titration boundaries were then pinned to the
> boundaries (10⁻¹² and 10⁻⁷ for H1 and H9; 10⁻¹¹ and 10⁻⁶ for H3 and FluB)."*

So for H1 and H9 the floor is 10⁻⁷ M, i.e. **−log₁₀(K_D) = 7.0 exactly**, which is precisely
where our CR6261 germline sits for both antigens. **[verified in data]** Our files bear this
out: minimum `DMS_score` is exactly 7.0 in both `CR6261_H1_Phillips_2021.csv` (213 rows at the
floor, of 1,886) and `CR6261_H9_Phillips_2021.csv` (371 rows at the floor, of 1,841). This is
consistent with the paper's binder fractions for CR6261 — **1,673/1,812 (92%) for H1** and
**1,470/1,812 (81%) for H9** — i.e. ~8% and ~19% non-binders.

Also note a sequence-fitting detail that explains the floor pile-up: *"Sequences leading to a
failed optimization were deemed 'non-binding'."*

**On the specific question "does it state that the germline precursor does not detectably
bind":** **No — not as an explicit sentence, for either antibody.** What the paper does say is
weaker and indirect:

- For CR9114 against H1: *"For H1, only a few mutations are needed to improve from the germline
  affinity"* — which presupposes a germline affinity worth improving on, not a non-binder.
- For CR6261 against H1: *"T29P is highly likely to occur first in scenarios that begin with
  H1, as this is the only single mutation that can improve H1 affinity, albeit rather
  modestly."* That is a strong depth-1 statement (1 of 11 single mutations does anything at
  all) and it is consistent with a germline at or near the floor, but it is not a statement
  that the germline does not bind.
- The germline is marked in Figure 1C–F and the per-antibody binder percentages are given, but
  no germline −log₁₀(K_D) value is quoted in the text.

**Recommendation:** we may state that *our* CR6261 germline measurement sits at the assay floor
of 7.0 for both H1 and H9, and cite the paper for the floor being 7.0 and for 8%/19% of
variants being non-binders. We should **not** attribute to the paper a claim that the germline
precursor does not detectably bind. If we need that claim, it has to be made as our own
reading of the deposited source data.

### 3.3 Does it report hit rate or affinity rising with mutation count? Yes, and it quantifies it

> **Yes. This is the part that bears on the novelty claim, and the answer is that for CR9114
> the depth structure is published, plotted, and quantified.**

Three distinct places:

**(a) Figure 5A,B is literally an affinity-versus-depth plot.** Caption, verbatim:

> *"(A, B) Maximum binding affinity achievable for sequences with a given number of mutations.
> For each antigen for (A) CR9114 and (B) CR6261, the maximum observed (circles) and
> model-predicted (triangles) affinity for each number of somatic mutations is shown."*

And in the body: *"Considering the maximum −logK_D achieved by sequences with a given number of
mutations (a proxy for time), we see that improvements in H1 binding can be realized early on,
whereas improvements in H3 binding are not possible until later, and even later for influenza
B (Figure 5A)."*

**(b) Figure 1D,F stratifies the affinity distribution by depth.** Caption: *"Distributions of
library binding affinities for (D) CR9114 and (F) CR6261 for each antigen (gray histogram,
right) separated by number of somatic mutations (boxplots, left)."*

**(c) It gives a count of required mutations per antigen, with a frequency criterion.**

> *"We observe three specific mutations that are required for binding to H3 (present at over
> 90% frequency in the set of binding sequences) … We also observe eight specific mutations
> that are required for binding to influenza B … Notably, these sets of required mutations in
> CR9114 exhibit a nested structure: mutations beneficial for H1 are required for H3, and
> mutations required for H3 are required for influenza B."*

Summary claim: *"For CR9114, binding affinity to H1 can be achieved through different sets of
few mutations with complex interactions. In contrast, a specific set of many mutations with
strong synergistic interactions is required to bind H3, and to an even greater extent,
influenza B."* Binder fractions: **97% H1, 11% H3, 0.3% influenza B** (63,419 / 7,174 / 198
of 65,091).

**Demonstrated vs claimed.** The required-mutation counts (3 for H3, 8 for influenza B), the
nested structure, the binder fractions and Figures 1D/F and 5A/B are **demonstrated** from
measurement. The evolutionary narrative built on top — that affinity to diverse antigens was
*likely acquired sequentially* in vivo — is explicitly **inference under a model**, and the
authors say so: *"we cannot conclusively determine how CR9114 and CR6261 evolved in vivo."*

### 3.4 Consequence for our hypothesis ordering — and the contradiction that matters more

**For CR9114: the deep-shell structure is prior knowledge. Risk 7.6 is REAL.** An eight-mutation
requirement for influenza B binding, a nested requirement structure, and a published
max-affinity-versus-mutation-count plot together mean that "CR9114's hits live deep" is
published. Our H1 ranking of CR9114 is **anticipated, not predicted**, and must be written that
way.

**For CR6261: the paper reaches the OPPOSITE conclusion from our measurement. This is the most
valuable finding in this document.**

The paper's verbal conclusion for CR6261 is that depth structure is *absent*:

> *"This hierarchical structure is in striking contrast to the CR6261 library, in which most
> variants can bind both antigens (92% for H1, 81% for H9), variants have a similar K_D
> distribution, and many variants display intermediate affinity to both antigens."*

> *"In contrast, for CR6261, improvements in binding can occur early on for both antigens
> (Figure 5B) and the number of uphill paths is more similar across single-antigen and
> sequential selection pressures."*

> *"In contrast, CR6261 could have acquired affinity to H1 and H9 in a continuous and
> simultaneous manner."*

Our own CR6261 H9 measurement is the reverse: **zero top-1% hits through depth 6**, then 4/297
at depth 7, 5/149 at 8, 7/48 at 9, 2/10 at 10, and 1/1 at depth 11.

These are not necessarily in conflict, and the distinction is the whole point:

| | Phillips Figure 5B | our measurement |
|---|---|---|
| statistic | **maximum** −log₁₀(K_D) in each depth shell | **density of top-1% variants** in each depth shell |
| what it shows | affinity *can* start improving at low depth | the *best* variants are concentrated at high depth |

A shell can contain a variant that improves measurably on the germline while containing no
member of the global top 1%. Both statements can hold at once. **But the paper's prose says
CR6261 has no hierarchical depth structure, and our top-1% curve says it has a strong one.**

Three things follow, and they are all in our favour if handled correctly:

1. **Our CR6261 H9 result is not anticipated — it is arguably contrary to the published
   conclusion.** The paper's own summary of CR6261 is "unconstrained, early, simultaneous". If
   we show a monotone top-1% depth gradient there, that is a genuinely new reading of a
   published landscape, not a restatement of it.
2. **It therefore has to be defended on exactly that distinction**, explicitly, in the write-up:
   maximum-per-shell versus top-tier density. A reviewer who knows this paper will otherwise
   read our claim as contradicting it, and we will lose the exchange by not having anticipated it.
3. **The paper supplies one piece of corroboration we should use.** Among sequential scenarios
   for CR6261, *"those beginning with H1 are more likely than those beginning with H9, as the
   first two mutational steps can improve affinity to H1 more than H9, and mutations late in
   maturation can improve affinity to H9 more than H1."* **"Mutations late in maturation can
   improve affinity to H9 more than H1" is a published, measured, depth-asymmetric statement
   pointing in our direction, and it is specific to H9** — the exact antibody-antigen pair
   where we see the deep-shell gradient. This should be quoted.

**Net effect on `docs/prior-art-depth-allocation.md` §7.6:** it should be rewritten rather than
deleted. It is **confirmed for CR9114** and **inverted for CR6261**. The §6.5 corpus statement
that "antibody affinity-maturation landscapes are deep-shell enriched, monotonically" is
partly published prior knowledge for one of the two and a novel, mildly contrarian result for
the other. Reporting it that way is both more accurate and more interesting than either
extreme.

### 3.5 One row-count discrepancy to note, unresolved

**[verified in data]** The paper reports `N_var` = **1,812** for CR6261 (both antigens) and
**65,091** for CR9114. Our files hold **1,886** rows (CR6261 H1), **1,841** rows (CR6261 H9)
and **48,840** rows (CR9114 H1).

- The CR6261 counts *exceed* 1,812 and differ between antigens, which is consistent with the
  paper's 1,812 being the variants measured for **both** antigens while our files were filtered
  per antigen. Plausible, not verified.
- The CR9114 H1 count of 48,840 is **25% below** the paper's 65,091 for H1, which the paper
  states has 97% measurable binding. This gap is not explained by anything found in the paper
  and is **not** the floor-pinning rule.

This is a provenance gap in our own ingestion, not in the source. It is a smaller issue than
Items 1 and 2 but it is a real one: **we should establish which filter produced 48,840 from
65,091 before any per-depth CR9114 count is published**, since a 25% unexplained attrition
could itself be depth-correlated and would then bias exactly the shell statistics we are
claiming. Recommended as a new open item.

---

## 4. Item 4 — Qiu et al., AAAI-24

Qiu J., Yuan H., Zhang J., Chen W., Wang H., Wang M. (2024). *Tree Search-Based Evolutionary
Bandits for Protein Sequence Optimization.* AAAI-24, Proc. AAAI Conf. Artif. Intell.
38(13):14686–14694. DOI [10.1609/aaai.v38i13.29386](https://doi.org/10.1609/aaai.v38i13.29386).
arXiv:[2401.06173](https://arxiv.org/abs/2401.06173).

**[full text]** — obtained as `https://arxiv.org/pdf/2401.06173` (19 pages incl. appendices)
and text-extracted locally. An arXiv HTML render **does** exist at
https://arxiv.org/html/2401.06173v1 and was also read; the earlier report that it does not is
incorrect. The AAAI PDF is additionally open at
https://ojs.aaai.org/index.php/AAAI/article/download/29386/30618.

### 4.1 Are the arms sequences, or depths/shells? — **Sequences.**

> **The arms are individual protein sequences.** The paper never uses the word "arm" in its
> method section at all; the bandit layer is a scoring function over candidate sequences.

The mechanism: a scoring function `F_t(x) = f_t(φ(x))` where `φ` is a pretrained protein
embedding and `f_t` a fitness model, with **UCB** (TreeNeuralUCB) or **Thompson Sampling**
(TreeNeuralTS) applied to that score. The tree is rooted at the wild type; each **node is a
protein sequence**; mutation children are one edit from a parent, so depth happens to equal the
accumulated mutation count; recombination children combine two active-set sequences and can
jump several layers. Candidates are generated by `TREE SEARCH`, scored, and the top-scoring ones
enter the active set.

`Algorithm 2: TREE SEARCH(A, F(·))` takes parameters `n₁, n₂, n₃, ρ`:

- `n₁` — number of random-mutation candidates per round
- `n₂` — number of recombination candidates per round
- `n₃` — active-set size
- `ρ` — fraction of the previous active set retained

**None of these is a per-depth budget.** `n₁`/`n₂` split mutation from recombination, not
shallow from deep. The paper does say of `ρ` that it *"can also be viewed as balancing depth and
width in tree search"* — but that is an interpretive remark about a fixed constant, not an
allocation mechanism, and nothing in the algorithm updates `ρ`.

### 4.2 Is the mutation-depth bound fixed, a regret term, or learned? — **A fixed constant of the analysis.**

> **`N` is a fixed quantity of the theoretical analysis. It is not an algorithm input, is not
> scheduled, and is never updated from observed returns.**

Verbatim, where `N` is introduced:

> *"we focus on a local searching region X̄ := {x : d(x, x^wt) ≤ N} ⊂ X and aim to identify the
> optimal sequences within X̄. Here d(·,·) denotes the hamming distance between any two
> arbitrary sequences in X and N controls the mutation counts (the depth of tree search)."*

`N` then appears only in the theory:

- It defines the region over which Bayesian regret is measured
  (`BayesRGT(T) = E[Σ_t (max_{x ∈ X̄} F(x) − F(x_t))]`).
- It sets the embedding-space radius `R := L_φ N`, via Lipschitz `φ`.
- It enters Theorem 3.6's bound as a `(1 + 4L_φ²N²/r²)` factor, and the confidence width `β_T`
  through a `d·ln(d·L_φ·N·‖f*‖_k·M·T²)` term. The expectation bound later simplifies and `N`
  drops out.

**`N` does not appear in Algorithm 1 or Algorithm 2 as an input or a variable.** It is not
enforced as a cap in the experiments. Mutation counts in the results are **reported outcomes,
not controlled allocations**: *"able to find top designs using reasonably small mutation
counts"*, *"single-digit mutation counts"*, top AAV sequences within ~5 mutations. The paper's
own depth claim is a complexity statement, not an allocation one: *"If we just use pure tree
search without any bandit learning, tree search has a complexity that grows exponentially with
depth. With the use of bandit learning, the sample complexity reduces to quadratic with tree
depth."*

Experimental protocol, for the record: **10 rounds, 100 queries per round**; datasets AAV, TEM,
AAYL49. The learned object `f_t` is a GP posterior resampled under a rarely-switching rule
(resample when `det(V_t)/det(V_pre(t)) > 2`) — **the fitness model is updated from returns, the
depth is not.**

### 4.3 Verdict

> **Qiu et al. 2024 is not prior art for mutation depth as an action space.** It is a bandit
> over sequences, inside a tree whose depth happens to be mutation count, with a
> Hamming-radius constant in the regret analysis. No published policy in this paper allocates
> across depth shells from its own measurements.

The §1.1(e) characterization in `docs/prior-art-depth-allocation.md` — *"a bandit, and a depth,
that never meet"* and *"the arms are sequences, not depths"* — was **correct**, and is now
confirmed from the paper rather than inferred from the abstract. Risk **7.7** and open item
**8.1** are closed.

### 4.4 Residual due diligence on the wider novelty claim

One further extended search was run for a policy that allocates budget across mutation-depth
shells from observed returns. **No prior art found.** The nearest items, and why each falls
short:

- **BoBa**, arXiv:[2606.26657](https://arxiv.org/abs/2606.26657), *Target-Aware Bandit
  Allocation for Scalable Surrogate Optimization in Chemical Space* **[abstract/snippet]** — a
  bandit that *does* shift sampling effort across **partitions of a library**, which is the
  right shape. But the partitions are structure-aware clusters of a small-molecule docking
  library, not Hamming-distance shells of a protein library, and the domain is not protein
  engineering. **This is the closest structural analogue found so far and should be added to
  the citation list as such.** It arguably displaces MC-UCB as the nearest "bandit over strata
  with per-stratum allocation" precedent in a design setting, and it should be read in full
  before the novelty claim is published.
- An out-of-distribution lab-in-the-loop active-learning study reports that selecting antigens
  by **average Hamming distance to already-labeled points** beat random
  **[abstract/snippet]** — Hamming distance as a diversity criterion, not as an arm set.
- A bandit over **batch sizes** as arms, in bird-sound classification **[abstract/snippet]** —
  adaptive budget, wrong action space, unrelated domain.

**Recommended new open item:** obtain BoBa (arXiv:2606.26657) in full. It is now the largest
residual risk to the §0 verdict, replacing Qiu et al.

---

## 5. Actions this document recommends

Ordered by how much they protect an external claim.

1. **Flag the 871 TrpB rows, and remove them from `valid`/`test`.** The identifiers are in
   `docs/trpb_imputed_variants_871.csv`. Replace the `UNVERIFIED_INCLUDES_IMPUTED` stamp with
   a row-level `label_source` column (`MEASURED` / `KNN_IMPUTED`). All 871 are `active = False`,
   so top-tier endpoints are unaffected; distributional metrics are. Report the proportion
   (0.544%) wherever the landscape is described. *(Requires a code change — not made here.)*
2. **Mark GB1 as fully measured.** No caveat is needed. Record the Supplementary file 1 vs 2
   split so the question is not reopened.
3. **Close discrepancy D-3** as a definitional artefact (§2.3), and record that ProteinGym's
   assay ID for our GB1 file is `SPG1_STRSG_Wu_2016`.
4. **Rewrite `docs/prior-art-depth-allocation.md` §7.6** as confirmed-for-CR9114 and
   inverted-for-CR6261 (§3.4). Pre-register the maximum-per-shell versus top-tier-density
   distinction *before* the H1 result is written up.
5. **Open a new item for the CR9114 H1 attrition**: 65,091 published → 48,840 in our file, 25%
   unexplained, potentially depth-correlated (§3.5).
6. **Close open items §8.1 and §8.2** of `docs/prior-art-depth-allocation.md`; upgrade entries
   8 and 18 of its §6 citation table to **[full text]**.
7. **Open a new item for BoBa** (arXiv:2606.26657) as the largest remaining novelty risk (§4.4).

---

## 6. Sources, with stable identifiers

| # | work | identifier | read |
|---|---|---|---|
| 1 | Johnston et al. 2024, TrpB 4-site landscape, PNAS 121(32):e2400439121 | [10.1073/pnas.2400439121](https://doi.org/10.1073/pnas.2400439121); PMC11317637; preprint [10.1101/2024.06.23.600144](https://doi.org/10.1101/2024.06.23.600144) | full text |
| 2 | Johnston et al. data deposit, CaltechDATA | [10.22002/h5rah-5z170](https://doi.org/10.22002/h5rah-5z170) | member `20230828_KNN_imputed_TrpB.csv` retrieved and verified |
| 3 | Johnston et al. raw reads, SRA | [PRJNA1127511](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1127511) | identifier only |
| 4 | Wu et al. 2016, GB1 4-site landscape, eLife 5:e16965 | [10.7554/eLife.16965](https://doi.org/10.7554/eLife.16965) | full text |
| 5 | Wu et al. Supplementary file 1 — measured fitness, 149,361 variants | [10.7554/eLife.16965.024](https://doi.org/10.7554/eLife.16965.024) | downloaded and verified |
| 6 | Wu et al. Supplementary file 2 — imputed fitness, 10,639 variants | [10.7554/eLife.16965.025](https://doi.org/10.7554/eLife.16965.025) | downloaded and verified |
| 7 | Wu et al. raw reads, SRA BioProject | [PRJNA278685](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278685) | identifier only |
| 8 | Phillips et al. 2021, CR9114/CR6261 landscapes, eLife 10:e71393 | [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393) | **full text (previously [not read])** |
| 9 | Qiu et al. 2024, AAAI-24 38(13):14686–14694 | [10.1609/aaai.v38i13.29386](https://doi.org/10.1609/aaai.v38i13.29386); arXiv:[2401.06173](https://arxiv.org/abs/2401.06173) | **full text (previously [abstract only])** |
| 10 | Yang et al. 2025, ALDE, Nat Commun 16:714 | [10.1038/s41467-025-55987-8](https://doi.org/10.1038/s41467-025-55987-8); PMC11739421 | full text (via PMC) |
| 11 | Cheng et al. 2022, ODBO | arXiv:[2205.09548](https://arxiv.org/abs/2205.09548) | carried from `prior-art-depth-allocation.md` |
| 12 | ProteinGym reference file, `DMS_substitutions.csv` (217 assays) | https://github.com/OATML-Markslab/ProteinGym | downloaded and searched |
| 13 | BoBa, *Target-Aware Bandit Allocation … Chemical Space* | arXiv:[2606.26657](https://arxiv.org/abs/2606.26657) | abstract/snippet — **new open item** |

### Could not obtain

- **PNAS publisher HTML** for Johnston et al. (`https://www.pnas.org/doi/10.1073/pnas.2400439121`)
  returned HTTP 403. Worked around via PMC11317637 and the bioRxiv PDF; no content is asserted
  from the publisher page.
- **Nature publisher HTML** for ALDE returned a 303 to an identity provider. Worked around via
  PMC11739421.
- **Johnston et al. SI Appendix Fig. S16** (the imputation figure) was not read; it is inside the
  ~7.3 MB SI PDF. Nothing is asserted about it beyond the Results and Methods sentences quoted
  in §1.1. It was not needed, because the imputed identifiers were recovered directly.
- **`four-site_simplified_AA_data.csv`** (the measured four-site table in the CaltechDATA
  deposit) was located but not downloaded; the imputed complement was sufficient.
