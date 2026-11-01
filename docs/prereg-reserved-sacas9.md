# Pre-registration: the reserved KKH-SaCas9 landscape `sacas9_sg5`

**Status: SEALED and UNSCORED.** Written and committed 2026-11-01. As of this commit no
method, arm or campaign has been run on this landscape, at any budget, on any seed. No
`benchmark_comparison_results/` run directory scores it, `load_landscape('sacas9_sg5')`
reports no manifest threshold, and `scripts/analysis/` holds no output derived from it.

**The first person to score it must do so against this document unchanged.** Any later
change to the label, the tier, the budget, the endpoint or the comparator set is an
**amendment**, and an amendment must be dated, justified and committed *before* the arm
it affects runs — the §5.4 condition that `docs/HANDOFF.md` records cannot be met once an
outcome is known, and a false blindness claim is worse than no amendment.

**Why this document exists.** Rule 2.4 has now fired twice: on HER2, and on 2026-11-01 on
**CreiLOV**, because the `wt_anchor` intervention was built against CreiLOV's measured
failure (+31.771 at 48 paired seeds, t=+35.88, 48/0/0). A seven-arm baseline panel had
just established CreiLOV as the only non-development landscape it tested. Of the nine
registered landscapes, TrpB, GB1 and HER2 are development data by §2.4; CR9114 is, because
the generator fix was designed against its measured failure; CreiLOV now is; both eqFP611
channels and both CR6261 landscapes have been measured and re-measured repeatedly, and
CR6261's endpoint is degenerate at 19 hits. **The only generalization evidence that exists
anywhere in this project is eqFP611-blue, +1.625 at t=+7.83 on 38 of 48 seeds.**

This document reserves one more. Choosing the endpoint now is legitimate precisely because
**no outcome on this landscape exists yet**; it is the one moment in the project's life
when that choice can be made without contaminating it, and it is being made in writing,
in advance, with every number it rests on published below.

**Reproduce everything here with:**

```
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_sacas9.py validate
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_sacas9.py readouts
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_sacas9.py prepare
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_sacas9.py endpoint
```

`--corpus-dir` defaults to the **real repository root's** `downloads\combingym`, because
`downloads/` is gitignored and a git worktree cannot see it (rule 2.8). If you pass it
explicitly it must name the **corpus directory itself** — the folder holding `Data/` — not
the repository root, which still fails.

---

## 0. The verdict, first

| question | answer |
|---|---|
| Is the library definition validated? | **Yes, four ways, all agreeing.** §2 |
| Is the label unambiguous? | **Yes: one named readout, sg5 / PAM TGCCGT.** §3 |
| Does the registered 125-measurement endpoint transfer? | **No. It FAILS the gate.** §4.1 |
| Does an amended endpoint survive the gate? | **Yes: 5% tier, budget 50 — PASS.** §4.2 |
| Is the endpoint degenerate? | **No** at the 5% tier (18.46% removed at d>=2). §4.3 |
| Is `sacas9_sg5` a usable reserved landscape? | **Yes**, with the two limits in §8. |

The one thing it is **not**: a landscape that can carry the registered protocol unchanged.
It is 1,295 variants. §9 surveys what a landscape that could would require.

---

## 1. What is reserved

| field | value |
|---|---|
| registry key | `sacas9_sg5` (`veritas/calibration/landscapes.py::REGISTRY`) |
| protein | KKH-SaCas9 — PAM-relaxed *Staphylococcus aureus* Cas9, UniProt `J7RUA5` |
| family | **Cas9 endonuclease — a sixth structural family, new to this corpus** |
| library | 8 sites, mixed alphabets 2·2·3·3·3·2·2·3 |
| positions | 887, 888, 889, 985, 986, 988, 989, 991 (1-based on a 1,053-residue parent) |
| parent residues | `LNANNLLR` |
| design space | **1,296, and all 1,296 are present — combinatorially complete** |
| eligible universe | **1,295** (the parent is not a library member) |
| depths present | 1–8, with 12 / 62 / 180 / 321 / 360 / 248 / 96 / 16 variants |
| label | sg5 log2 enrichment, PAM `TGCCGT` (§3) |
| prepared table | `scripts/data/SACAS9_KKH_SG5_Thean_2022.csv`, 1,295 rows, 73,055 bytes |
| sha256 (LF) | `3d9722c3cb3ce628e80e691e11565fa74b0bdf397ba78e011782ce82e9b85a00` |
| upstream source | `Data/DMS/Raw/SaCas9_raw.xlsx`, CombinGym @ `c5dccd514a9acfe1089a5f1555faea98afd78019` |
| citation | Thean D.G.L., Chu H.Y., Fong J.H.C. *et al.*, "Machine learning-coupled combinatorial mutagenesis enables resource-efficient engineering of CRISPR-Cas9 genome editor activities", *Nature Communications* **13**(1):2219 (2022) — transcribed verbatim from CombinGym's `Data/Data_summary.xlsx`. **The DOI was not verified in this task and is therefore not quoted**: no network lookup was made, and "a caption is not a citation" (HANDOFF §3). |

**Why SaCas9 and not one of the other reserve landscapes.** `docs/corpus-audit.md` §3 lists
five `reserve` entries. HIV-1 protease is 1,175 rows at **25.5%** saturation and falls from
100% of depth 1 to 0.6% of depth 9; RhlA ×3 is 2,346 rows at **0.35%** saturation and holds
nothing at all beyond depth 3, which is why §6.3's degeneracy guard already declares its
covariate undefined; SpCas9 is 597 of 952 rows at **3.44%** saturation with **depth-correlated
measurement attrition** (95.8% retention at depth 2 falling monotonically to 47.2% at depth
8, §5.1). SaCas9 is the only one that is combinatorially complete, and completeness is what
makes it worth reserving: *"every genotype in the design space is present, at every depth,
in both tables. There is no room for a filter to hide."*

**The pin did not move and the lock did not change.** `SaCas9_raw.xlsx` is already a pinned
entry of `scripts/data/cross_protein_corpus.lock.json`, so reading it changes *what is
read*, not the pin — the same move `docs/cr9114-reprep.md` made. Verified at the time of
writing:

```
freeze_cross_protein_corpus.py verify  ->  CORPUS VERIFIED: 29/29 files
                                           match CombinGym @ c5dccd514a9a
```

A `FAIL corpus.lock` inside a worktree is an environment artefact, not drift:
`freeze_cross_protein_corpus.py:55` resolves `downloads/combingym` against its own
`__file__` and `downloads/` is gitignored. Pass `--corpus-dir` pointed at the corpus
directory itself.

---

## 2. The library definition, validated four ways

A silent mis-mapping here is the failure mode that made HER2 read the trastuzumab light
chain, so the definition is not trusted from any single source.
`docs/dataset-inventory.md` §1 cross-validated the six CombinGym landscapes against three
sources; this is the same check against four. **All four agree, and
`scripts/data/prepare_sacas9.py validate` refuses to write anything if they do not.**

| # | source | what it says |
|---|---|---|
| 1 | CombinGym's registry, `Data/Data_summary.xlsx` | `L887, N888, A889, N985, N986, L988, L989, R991`; 8 sites; variant space 1,296; DMS rows 1,296; sequence length 1,053 |
| 2 | the positions that actually **vary** across all 1,296 clean `genotype` strings | `[887, 888, 889, 985, 986, 988, 989, 991]` — **exactly** source 1, none missing and none spurious |
| 3 | the wild-type FASTA `SaCas9(SaCas9-KKH).fasta`, read at those 1-based positions | spells `LNANNLLR` — **exactly** source 1 |
| 4 | the **raw** table's own column header and parent row | header `887-888-889-985-986-988-989-991`; parent row `L-N-A-N-N-L-L-R` — **exactly** source 1 |

A **closed loop**, which is what eqFP611 had and HER2 did not: exactly **one** of the 1,296
`genotype` strings is byte-identical to the FASTA, and it is the row whose `n_mut` is 0. The
library's parent and the deposited wild type are the same sequence, not two sequences that
happen to agree at eight positions.

Four further invariants, all measured, none assumed:

| invariant | result |
|---|---|
| `n_mut` equals Hamming-distance-to-parent | **0 of 1,296 rows disagree** |
| rows varying at any of the other 1,045 sequence positions | **0 of 1,296** |
| per-site alphabets multiply to the design space | `2·2·3·3·3·2·2·3 = 1,296` = distinct genotypes = rows = the registry's own 1,296 |
| `mutant`-string round trip through `veritas.calibration.genotypes` | **0 of 1,296 fail** |

Per-site alphabets, published because they are the library's definition and because a
20-residue block per site would encode 140 columns that cannot occur (§6):

| site | residues | | site | residues |
|---|---|---|---|---|
| `L887` | `LR` | | `N986` | `LNT` |
| `N888` | `NQ` | | `L988` | `DL` |
| `A889` | `ARS` | | `L989` | `LR` |
| `N985` | `DLN` | | `R991` | `IQR` |

Genotypes are canonicalised to full library length with the parent residue filled at
omitted positions, and **equal length is asserted** — nothing here zips or truncates
(rule 2.1), and no sixth canonicalisation was written: the module is
`veritas.calibration.genotypes`, imported. Before any of it,
`prepare_sacas9.py validate` reproduces the two pre-registered parser anchors from the
committed tables (rule 2.2): **TrpB 168/2,242 = 7.4933%** and **GB1 72/2,167 = 3.3226%**,
depth 0 excluded from both shells and retained in the universe. Both PASS.

`load_landscape` re-checks the positions and the parent residues against the data on every
load, so a later mis-mapping raises rather than passing.

---

## 3. The label — addressed, not papered over

`docs/corpus-audit.md` §5.2 records the defect: *"SaCas9's clean phenotype silently averages
three of the six sgRNAs the raw table reports."* Re-derived here: the clean `Mean` equals
`mean(sg1, sg2, sg3)` to **1.110e-16** on all 1,296 rows and differs from the mean over all
six on **1,296 of 1,296**. Clean `sg1`–`sg3` equal raw `sg1`–`sg3` at **exactly 0.0**, so
nothing turns on which table is read — only on the provenance being the one table that
carries all six columns.

Two measurements decide the label. Both are assay structure, model-free, and could not
differ between two methods.

### 3.1 The six readouts are close to rank-independent

Spearman rho over the 1,295 mutants:

|  | sg1 | sg2 | sg3 | sg5 | sg6 | sg7 |
|---|---|---|---|---|---|---|
| **sg1** | 1.000 | 0.151 | 0.143 | 0.011 | 0.077 | 0.075 |
| **sg2** | 0.151 | 1.000 | 0.126 | −0.011 | −0.006 | 0.016 |
| **sg3** | 0.143 | 0.126 | 1.000 | −0.005 | 0.079 | 0.034 |
| **sg5** | 0.011 | −0.011 | −0.005 | 1.000 | 0.018 | −0.016 |
| **sg6** | 0.077 | −0.006 | 0.079 | 0.018 | 1.000 | 0.018 |
| **sg7** | 0.075 | 0.016 | 0.034 | −0.016 | 0.018 | 1.000 |

"SaCas9 activity" is **not one latent phenotype**. It is six nearly independent fitness
landscapes over one genotype set, and |rho| <= 0.018 within sg5/sg6/sg7. An average over
them is a composite that no assay measured, correlated only **0.32–0.52** with each of its
own constituents (`mean(all 6)` row) — so neither CombinGym's three-of-six `Mean` nor a
mean over all six is a label a reserved landscape can be defined by. **A reserved landscape
whose endpoint is ambiguous is not reserved, it is merely unused.**

### 3.2 On sg1, sg2 and sg3 there is nothing to optimise

| sgRNA | PAM | parent | library max | headroom | mutants above the parent | parent's percentile | distinct values |
|---|---|---|---|---|---|---|---|
| sg1 | `ACAAGT` | 1.2211 | 1.5423 | 0.3212 | **3** | 99.77% | 1,295 |
| sg2 | `GGTGGT` | 1.3108 | 1.2712 | **−0.0395** | **0** | **100.00%** | 1,295 |
| sg3 | `TGGAGT` | 1.1737 | 1.4246 | 0.2509 | **5** | 99.61% | 1,295 |
| **sg5** | **`TGCCGT`** | **0.0530** | **0.9267** | **0.8737** | **495** | **61.78%** | 1,294 |
| sg6 | `TCTTGT` | −0.1339 | 0.6071 | 0.7409 | 1,012 | 21.85% | 1,295 |
| sg7 | `GGGTGT` | −0.0823 | 0.4530 | 0.5353 | 940 | 27.41% | 1,295 |

On the three readouts CombinGym's `Mean` covers, the KKH parent already sits at the
**99.6th–100.00th percentile** of the 1,295 mutants, and on sg2 it exceeds *every* mutant
in the library. Every member of a top-tier defined on those readouts is **worse than the
starting enzyme**: 3, 0 and 5 of 65 tier members beat the parent at the 5% tier, against
**65 of 65 on each of sg5, sg6 and sg7**. A benchmark that rewards finding variants worse
than wild type is not a benchmark. **So the clean `Mean` phenotype is unusable as an
optimisation endpoint independently of being undocumented** — which is a stronger reason
than the audit's, and it is the reason recorded here.

sg5, sg6 and sg7 are the readouts the library is actually for: KKH-SaCas9 was engineered to
relax the PAM requirement, and these are target sites it still handles poorly. Making it
active on them is the engineering problem the data was generated to answer.

### 3.3 The registered label

> **`DMS_score` is the log2 enrichment of the KKH-SaCas9 variant on sgRNA `sg5`
> (PAM `TGCCGT`), read verbatim from the `sg5 (PAM: TGCCGT)` column of
> `Data/DMS/Raw/SaCas9_raw.xlsx` at CombinGym commit `c5dccd514a9a`. It is one named
> measured readout. It is not an average, not a composite, and not CombinGym's `Mean`.**

**Why sg5 of the three usable readouts:** the largest parent-relative headroom, **0.8737**
log2 units against sg6's 0.7409 and sg7's 0.5353. That is the clearest optimisation problem
of the three; it is a property of the assay; and it is declared here with the landscape
unscored. It is *not* the readout with the most favourable endpoint geometry — sg7 removes
27.69% at d>=2 against sg5's 18.46% — and it is recorded that way so the choice cannot
later be read as having been made for the endpoint's benefit.

**sg6 and sg7 are the pre-declared second and third reserve.** They are deliberately NOT
registered in `REGISTRY`, no prepared table exists for either, and they must not be scored
when sg5 is scored. `prepare_sacas9.py --readout sg6|sg7 prepare` produces them without new
code, and their endpoint geometry is published in §4.3 so a later promotion needs no new
inspection. **Caveat, carried verbatim from `docs/dataset-inventory.md` §2:** they are *one
protein with three phenotypes over one shared genotype set*. Like eqFP611-blue/red and
CR6261-H1/H9, they are **not three proteins and not three families**, and must never be
counted as three independent generalization data points.

Sign: **higher is better.** No transformation is applied; the value is the published log2
enrichment. **No censoring:** exactly one variant sits at the minimum (−0.814488) and one at
the maximum (0.926728), so there is no floor or ceiling mass — unlike CR9114's 934 rows at
7.0 or HER2's 5,832 at 0.1. **1,294 distinct values among 1,295 variants** (one tied pair,
not at the tier boundary).

**Label provenance:** every value is a direct measurement from the published enrichment
table. **No value is imputed** — the library is combinatorially complete and the raw table
reports all six sgRNAs for all 1,296 genotypes with no blanks, so the TrpB problem (871
KNN-imputed rows spread across all three splits, `docs/label-provenance.md` §1) does not
arise here and no `label_source` column is needed. `stage` is `train` on every row; there
is no held-out split and therefore no absent-label placeholder, so rule 2.3's HER2 hazard
does not arise either.

---

## 4. The endpoint — the crux, and it required an amendment

### 4.1 The registered endpoint does not transfer, and this is measured

The protocol's 125-measurement budget is **9.65%** of a 1,295-variant universe, and a 1%
tier is **13 variants** — the same degeneracy that makes CR6261 unusable at 19.
`veritas/preflight.py::check_endpoint_viability` is the instrument, and its verdict on the
registered protocol unchanged is:

```
FAIL  endpoint.sacas9_sg5   a budget of 125 measurements is not smaller than the
                            13-variant hit tier: an arm can reach the whole tier, so
                            the endpoint cannot discriminate
      measured {'raw': 13, 'ceiling': 12, 'informativeness': 0.0769, 'budget': 125,
                'tier_size': 13}
```

The full gate at protocol defaults exits **1**. This is not a preference for an amendment;
it is what forces one.

### 4.2 The amendment, and its verdict

Two parameters change, and nothing else. Both are pre-registered here, before any arm runs.

**The tier: 5%, chosen by the protocol's own ladder, not by me.**
`scripts/freeze_cross_protein_corpus.py` already defines `TIER_LADDER = (0.01, 0.02, 0.05,
0.10)` and `MIN_TIER_MEMBERS = 30`, and deepens the tier to the shallowest rung holding at
least 30 members — the mechanism that already deepened CR6261 to 2%. On 1,295 variants: 1%
is 13 members and 2% is 26, both under 30, so the ladder picks **5% = 65 members**. The
`MIN_TIER_MEMBERS` rule is load-bearing here and not merely satisfied: the 2% tier would
*also* have been **NEAR-DEGENERATE**, removing only 3.85% of the hit set at d>=2, under the
5% informativeness floor.

The hit set is fixed by an **order statistic**, not an interpolated percentile, so it is
unambiguous:

> **threshold = the 65th largest `DMS_score` = `0.37560747100293446`, at-or-above
> semantics, giving exactly 65 hits.** The 66th largest is `0.37438035248487095`, so there
> is **no boundary tie mass** and the cut is not arbitrary. (For the record, `numpy`'s
> interpolated 95th percentile is `0.37474848804028993` and yields the same 65 hits.)

**The budget: 50 = 10 initial + 4 rounds × 10.** A uniform **0.4×** scaling of the
registered `25 + 4 × 25 = 125`, which preserves both the 20% initial fraction and the
five-round structure exactly. The scale was fixed by two requirements: the budget must be
smaller than the hit tier (the gate's second condition), and it must not approach
enumeration of a 1,295-variant universe.

| ratio | `sacas9_sg5` @ 50 | TrpB @ 125 | CR6261-H1 @ 125 (degenerate) |
|---|---|---|---|
| budget / hit tier | **76.9%** | 7.8% | 658% |
| budget / universe | **3.86%** | 0.078% | 6.6% |

```
PASS  endpoint.sacas9_sg5   sacas9_sg5: d>=2 ceiling 53 of 65 (18.5% removed) and a
                            budget of 50 under a 65-variant tier - usable as the
                            primary endpoint
      measured {'raw': 65, 'ceiling': 53, 'informativeness': 0.1846, 'budget': 50,
                'tier_size': 65}
```

> **The registered primary endpoint is: diversity-adjusted top-5% yield at Hamming
> d >= 2, over a 50-measurement budget of 10 initial + 4 × 10, excluding the initial
> sample, at 48 paired seeds.** Secondary: raw top-5% yield, and the same contrast at
> d >= 3 (ceiling 29 of 65, 55.38% removed). Tertiary, descriptive only: best-so-far
> trajectory, and whether the global optimum was found.

Uniform-random expectation at this budget is **2.510** raw hits against a per-arm ceiling of
`min(50, 65) = 50` raw and 50 at d>=2 — a ~20× dynamic range, so the endpoint has room to
separate arms.

### 4.2.1 What this endpoint does about redundancy, which is why d>=2 is the primary

The raw count can be satisfied by re-finding one basin, and that is not a hypothetical. On
TrpB the one-hot GP surrogate adds **+5.00 raw top-1% hits at t=+5.33** — large and highly
significant — while the same contrast at **d>=2 is +0.10 at t=+0.29, not separable**. The
entire raw gain was **redundant hits within Hamming 1 of a hit already counted**. A raw-count
primary would therefore have reported a decisive surrogate improvement that the diversity
adjustment shows to be nothing.

On this landscape the pressure is far worse, not better. TrpB's 125 measurements are 0.078%
of 159,999 variants; here 50 measurements are **3.86%** of 1,295, and the whole library sits
within Hamming 8 of itself with 575 of 1,295 variants inside depth 4. The denser the
universe relative to the budget, the easier it is to fill a batch with neighbours of one
good variant. The instrument that measures this is
`veritas/calibration/endpoint.py::diversity_ceiling`, and on the registered tier it says the
adjustment has real work to do:

| | raw | d>=2 | d>=3 |
|---|---|---|---|
| 5% tier ceiling | 65 | **53** | 29 |
| hits the adjustment removes | — | **12** | **36** |
| share of the hit set removed | — | **18.46%** | **55.38%** |

The adjustment is **not vacuous here, and that is the measurement that matters.** The tier's
2,080 unordered pairs have distance histogram `{1: 18, 2: 100, 3: 310, 4: 518, 5: 593,
6: 359, 7: 155, 8: 27}`, so the **minimum pairwise Hamming distance within the hit set is 1**
and 18 pairs sit at it. Contrast HER2, where the minimum pairwise distance across the whole
top-1% set is **already 2**, so a d>=2 filter can act on nothing and the ceiling equals the
raw count (113 == 113) for every arm. On `sacas9_sg5` there is redundancy in the hit set for
the rule to find, and it finds 12 hits' worth.

So **the primary endpoint cannot be satisfied by re-finding the same basin.** Concretely: a
selection that includes a hit and then its Hamming-1 neighbours scores every one of them raw
and only the first at d>=2, because each neighbour is at distance 1 from the kept hit. (The
rule is not a blunt instrument either — two Hamming-1 neighbours of a common centre that
differ at *different* sites are distance 2 apart and both survive if the centre is not kept.
It penalises re-finding a *point*, not exploring a region.) The registered reporting
requirement in §5.1 keeps raw, d>=2 and d>=3 together for this reason — the TrpB case is a
worked example of a raw gain that reverses under the adjustment — and §5.3 requirement 5
(same sign at all three) exists so that this landscape cannot produce the same artefact.

One honest limit on the instrument: at d>=3 the ceiling is 29 of 65, so an arm spending 50
measurements can at most report 29 — the endpoint stops being budget-limited and becomes
ceiling-limited. That is why **d>=2 is the primary and d>=3 is a secondary**, and it is
stated here rather than discovered when two arms both post 29.

### 4.3 The degeneracy verdict, at every rung and for every readout

`INFORMATIVENESS_FLOOR` is 5% of the hit set removed; `DEGENERATE` means the diversity
adjustment is a numerical no-op for every arm, as on HER2 (113 == 113).

| readout | 1% (13) | 2% (26) | **5% (65)** | 10% (130) |
|---|---|---|---|---|
| **sg5** (registered) | 7.69% | **3.85% NEAR-DEGENERATE** | **18.46% USABLE** | 35.38% |
| sg6 (2nd reserve) | **0.00% DEGENERATE** | 11.54% | 24.62% USABLE | 40.00% |
| sg7 (3rd reserve) | 15.38% | 26.92% | 27.69% USABLE | 39.23% |
| sg1 | 46.15% | 46.15% | 36.92% | 42.31% |
| sg2 | 61.54% | 65.38% | 38.46% | 43.85% |
| sg3 | 53.85% | 61.54% | 41.54% | 43.85% |

At the registered 5% tier the endpoint is **USABLE, not degenerate and not
near-degenerate**, on sg5 and on both reserve readouts. Note that it is degenerate or
near-degenerate at the shallower rungs on exactly the three readouts that pose an
optimisation problem — which is why the tier ladder, and not the 1% default, had to be used.

### 4.4 Two endpoints that were weighed and rejected

**A deepened 10% tier** (130 hits, 35.38% removed at d>=2) survives the gate and would
allow a larger budget. **Rejected**, because the protocol's own ladder stops at the first
rung with >=30 members and taking a later rung would be a preference rather than a rule; and
because a top-10% "hit" on a 1,295-variant library is a weaker claim than a top-5% one.
Recorded here so that moving to it later is visibly an amendment.

**P(find the global optimum)** — rejected as the primary, kept as a tertiary descriptive
statistic. The optimum is unique and well defined: `L887R:A889S:L988D:R991I` at
`0.9267276666300931`, **depth 4**, one variant, no tie. But a Bernoulli endpoint at 48
paired seeds resolves only very large differences in success probability (roughly >= 20
percentage points by McNemar), against a count endpoint whose measured MDE elsewhere in this
project has been **0.42–2.67 hits**. Using it as the primary would substitute a low-power
test for a high-power one. It was also checked for a specific trap that the sg3 readout does
have and sg5 does not: under sg3 the global optimum sits at **depth 1**, so exhausting the
12 single mutants finds it for free; under **sg5 the 5% tier holds zero hits at depth 1**
(depth histogram of the tier: 2:3, 3:10, 4:15, 5:22, 6:8, 7:6, 8:1) and the optimum is at
depth 4.

**A reduced budget with the 1% tier** (e.g. 12 measurements, so that 12 < 13) — rejected as
absurd: it is below the initial sample of any protocol variant and cannot support a
multi-round campaign at all.

### 4.5 The threshold, and rule 2.3

Rule 2.3: thresholds come from manifests; never recompute one. **No run has scored this
landscape, so there is no manifest, and that is the correct and expected state of a reserved
landscape — it is not a defect.** `veritas/calibration/landscapes.py::_manifest_threshold`
was fixed on `main` on 2026-11-01 in `017de36`: it previously knew only the legacy flat
schema and, on no match, **silently recomputed the 99th percentile**, which rule 2.3 forbids
outright. It now reads both schemas, **raises** when runs exist for a landscape but their
threshold cannot be read, and returns `None` only when no run has scored it. For
`sacas9_sg5` it returns `None`, and the loader records a bootstrap percentile whose source
string says so:

```
threshold        = 0.54567545643975834     # the 99th percentile - a BOOTSTRAP value
threshold_source = "computed (N manifest(s) read, none scored sacas9_sg5)
                    - NOT a measured threshold"
# or, from inside a git worktree, where the gitignored results root is invisible:
threshold_source = "computed (no results directory at ...benchmark_comparison_results)"
```

Both strings mean the same thing and both are correct here. **That value is a bootstrap
percentile, not a measured threshold, and it is NOT the registered threshold of this
landscape.** The loader's percentile is hard-coded at 99.0 and therefore describes the 1%
tier that §4.1 just failed. **The registered threshold is the 5%-tier order statistic
`0.37560747100293446` fixed in §4.2**, and it too is a bootstrap value until a run records
one.

**The first scored run's manifest becomes the authority.** It must record the registered
value as `results.top1.threshold_value`, after which rule 2.3 governs, the loader reads it
from the manifest, and this document's number becomes a cross-check rather than the source.
If the two ever disagree, the eligible universe moved between them, which is a finding and
not something to average.

This is a wart and it is written down rather than hidden: `tier_fraction` on
`LandscapeSpec` drives `tier_size` but not the loader's percentile, so a reader who takes
`load_landscape('sacas9_sg5').threshold` at face value gets the wrong tier. Do not do that.
`prepare_sacas9.py endpoint` computes both and prints both verdicts for exactly this reason.

---

## 5. The predictions, written down in advance

### 5.1 The comparison

48 paired seeds, seed as the unit of independence, budget gate-verified matched at 50 on
every arm and every seed, 10,000-sample whole-seed bootstrap. **A verdict needs a paired
test, not a mean difference** (rule 2.11): `eqFP611-blue 9.50 vs 8.79` was published as a
WIN and the paired test on the same numbers gives t=0.77, better on 25 of 48 seeds.

Carrying the post-hoc reporting requirement that `docs/HANDOFF.md` §1 records, any
comparative claim on this landscape must report **all four** of:

1. Δ against `ridge` (the registered reference, unchanged);
2. Δ against the **highest-mean arm of the seven-arm baseline panel** (`random_forest`,
   `gbm_onehot`, `knn_hamming`, `gp_onehot`, `greedy_measured`, `ridge`, `random`) — because
   ridge is 1st of seven on CreiLOV, 2nd on TrpB and **6th of seven on CR9114**, so a
   single-reference design reports advantages that are not comparable across landscapes;
3. `ridge − random` with **both** the t statistic and the sign test;
4. the contrast at d>=2 (the primary) **and** at raw and d>=3.

Or be stated as *"advantage over ridge, comparator strength unmeasured"* — never as an
advantage over "the baselines".

`alde` must **not** be entered. 140 of its 160 one-hot columns are constant on this library
(§6), which is the condition under which it silently returns index 0 forever. Availability is
not compatibility (rule 2.6).

### 5.2 What the covariates predict, and it is a loss

| covariate | value | comparison |
|---|---|---|
| shallow-shell enrichment, `1 <= depth <= floor(8/2) = 4`, at the 5% tier | **0.9702×** (28/575 = 4.8696%) | CreiLOV 1.875×, eqFP611-blue 1.854×, **GB1 ≈ 3.32×**, eqFP611-red 0.634×, CR9114 0.584×, CR6261-H1 0.210× |
| mean depth of the hit set | **4.677** of 8 | the deep-shell regime |
| hits at depth 1 | **0 of 12** | nothing is reachable by exhausting the singles |
| hits reachable within the generator's documented radius-4 ceiling | **28 of 65 = 43.1%** | published as a covariate; the depth-ceiling hypothesis is itself falsified (HANDOFF) |
| parent functionality for the measured phenotype | **not functional** (0.0530, 61.78th percentile) | the 9-of-9 regularity ties this to near-zero shallow enrichment |
| primary shell size | 575 variants | above the >=500 resolvability floor |

`docs/HANDOFF.md` §8 pre-registers that the engine **loses**, on the grounds that every
independent landscape sits below GB1's shallow-shell enrichment. `sacas9_sg5` sits at
**0.9702× against GB1's ≈3.32×**, its hits are deep (mean 4.677 of 8), and its parent is
non-functional for the measured phenotype. Every covariate points the same way.

> **Pre-registered prediction: the engine loses or is not separable on `sacas9_sg5`.**
> Point prediction for Δ = `panoptes_default − ridge` on the primary endpoint:
> **Δ ∈ [−8, +2]** diversity-adjusted hits at d>=2. The most likely overall outcome is
> **mechanism supported, generalization undetermined**, which §8 of `HANDOFF.md` already
> records in advance as a legitimate reportable result specifically so that it cannot be
> re-narrated as failure later.

### 5.3 What would count as the engine generalizing

**All five** of the following, or it does not count:

1. Δ > 0 against `ridge` on the primary endpoint;
2. Δ > 0 against the **highest-mean arm of the seven-arm panel** — not only against ridge;
3. the bootstrap 95% **lower** bound above 0;
4. the sign test over 48 seeds at p < 0.05, i.e. a seed majority and not a coin flip;
5. the **same sign** at raw, d>=2 and d>=3.

Requirement 2 exists because "no engine win against the strongest available baseline on any
landscape except TrpB" is the current state. Requirement 3 exists because the CreiLOV
confirmation set a dominance threshold in advance and cleared it on the lower bound.
Requirement 4 exists because eqFP611-blue failed exactly it. Requirement 5 exists because
`knn_hamming` beat the engine on CR9114 at all three endpoints while greedy's and RF's leads
dissolved at d>=3.

### 5.4 What would count as the engine not generalizing

- **A decisive loss:** Δ <= 0 against the highest-mean panel arm with the bootstrap 95%
  **upper** bound below 0. Combined with §5.2 this is the predicted outcome.
- **A tight null:** |Δ| below the measured MDE *and* the MDE itself smaller than the
  smallest effect this project has claimed elsewhere — the CR9114 shape (+0.042 against an
  MDE of 0.924), which is a null and not a power failure.
- **Undetermined, and reportable as such:** |Δ| below the measured MDE with the MDE too
  large to exclude a claimed-size effect. *"Not separable" is not "identical"*: the point
  estimate, the interval and the MDE must all be quoted.

### 5.5 One sharp, falsifiable side prediction

The fabricated wild-type training row that `select_next_micro_batch` injects carries the
literal label `y = 0.0` whatever the assay's scale. On `sacas9_sg5`:

| quantity | value |
|---|---|
| label mean / SD | −0.005590 / 0.230225 |
| label range | 1.741215 |
| **`y = 0.0` in SDs from the centre of the distribution** | **0.0243** |
| the parent's **measured** sg5 value in SDs from the centre | 0.2544 |
| `|fabricated − measured|` in label SDs | **0.2301** |
| does `y = 0.0` lie inside `[min, max]`? | **yes** |

**0.0243 SDs is the smallest violation of any landscape in this corpus** — smaller than
GB1's 0.20, TrpB's 0.35, eqFP611-red's 0.80, blue's 1.19, CR9114's 14.14 and CreiLOV's
15.68. The fabricated row is almost exactly where the real parent sits.

> **Prediction: `--ablate-components wt_anchor` must be a NULL on `sacas9_sg5`.** If it is
> not, the WT-row account needs revising again. And a loss on this landscape **cannot** be
> attributed to that defect, which closes off the one explanation that was available for
> CreiLOV. Recall that the label-scale theory is already falsified — CR9114 violates the
> scale by 14.14 SDs and 2.47 label-ranges and still measures Δ=+0.042 — so this landscape
> is a second, independent test of it, from the opposite end.

### 5.6 What this landscape can and cannot support if the result is positive

**Can:** *"on one independent protein from a sixth structural family, over a
combinatorially complete 1,296-variant library at a 50-measurement budget, the advantage
over the strongest of seven baselines is X."*

**Cannot:** a generality claim. One landscape is one landscape; 48 campaigns on one
biological system are not 48 independent replicates (§5.6 of `HANDOFF.md`); the budget and
tier are **amended**, so nothing measured here is directly comparable to a 125-measurement
figure on another landscape; and sg6/sg7 are the same protein, so they cannot be counted as
further families.

---

## 6. Arm compatibility, recorded before anyone is surprised by it

`veritas/preflight.py::REFERENCE_CONSTANT_COLUMNS["sacas9_sg5"] = (140, 160)`. The ALDE
family allocates a column per amino acid per site regardless of what the library contains,
so 8 sites × 20 = 160 columns of which the mixed 2·2·3·3·3·2·2·3 design can occupy only 20.
**140 constant columns**, each of which normalises to `0/0 = NaN`, after which `np.argmax`
on an all-NaN array returns index 0 forever. This is `INVALID_NUMERICAL` under §5.8's arm
validity statuses, not `UNAVAILABLE`: ALDE will report `available: true` here and will not
function. Do not enter it.

The declared per-site alphabets are legitimate to take from the observed data **only
because the library is combinatorially complete** — observed and designed alphabets are the
same set, which is not true of a sparse library like HER2 and is why
`allowed_aa_per_position` is documented as never inferred from observed sequences.

---

## 7. What the first person to score it must do, and nothing more

The library, the label, the threshold, the tier, the budget, the endpoint, the comparators
and the predictions are all fixed above. One bounded engineering step remains, and **its
field values are pre-registered here so that it is a transcription and not a decision**:

```python
# benchmark_panoptes_isolated.py :: LIBRARIES
"sacas9_sg5": LibrarySpec(
    key="sacas9_sg5",
    label="KKH-SaCas9 eight-site PAM-interacting library, sgRNA sg5 / PAM TGCCGT "
          "(Thean et al. 2022)",
    dataset=REPO_ROOT / "scripts" / "data" / "SACAS9_KKH_SG5_Thean_2022.csv",
    accession="SACAS9_KKH",                      # UniProt J7RUA5
    positions=(887, 888, 889, 985, 986, 988, 989, 991),
    parent_residues="LNANNLLR",
    structure_pdb=None,                          # deliberately None - see below
    published_optimum="L887R:A889S:L988D:R991I",
    published_measured_count=1_296,
    published_imputed_count=0,
    provenance_note=(
        "1,296 of 1,296 genotypes, combinatorially complete over a mixed "
        "2.2.3.3.3.2.2.3 design; every value a direct measurement, none imputed. The "
        "label is one named readout (sg5, PAM TGCCGT) read from the pinned RAW table, "
        "NOT CombinGym's clean `Mean`, which averages three of six sgRNAs. See "
        "docs/prereg-reserved-sacas9.md."
    ),
    labelled_stage=None,                         # no held-out split; stage is 'train'
    pool_restricted=False,                       # 100% saturated: generation always
                                                 # proposes a variant that exists
    allowed_aa_per_position=("LR", "NQ", "ARS", "DLN", "LNT", "DL", "LR", "IQR"),
    assay_units="log2 enrichment on sgRNA sg5 (PAM TGCCGT)",
    detection_floor=None,                        # one variant at the min, no floor mass
    detection_ceiling=None,                      # one variant at the max
    source_publication="Thean et al. 2022, Nat Commun 13(1):2219",
),
```

The wild-type sequence is already recorded in `scripts/data/library_wt_sequences.json`
under `SACAS9_KKH_SG5_Thean_2022.csv` (1,053 residues, offset 1, 1,295 rows), byte-identical
to the pinned FASTA and to the library's own `n_mut = 0` row.

**`structure_pdb` is registered as `None` on purpose, and raising it is an amendment.** No
SaCas9 structure is pinned in the corpus lock, and the spatial veto has already been
measured reading **a different protein's geometry** on CreiLOV — `8QI9` maps 11 of 15
positions and **all 11 disagree with the library's parent residue** — where it nonetheless
*helped*, by raising candidate-pool hit density 2.8×, so the benefit was pool density and
not biophysics. An unpinned PDB chosen at scoring time is an unregistered degree of freedom
on a reserved landscape. If a structure is wanted, pin it, verify that all eight library
positions map and that the parent residue agrees at every one, and date the amendment.

Also required before the run, and not optional:

- `.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library sacas9_sg5 --arms <arms>
  --run-dir <dir> --budget 50` — **with `--run-dir`**. `generality_sweep.ps1` invokes the
  gate *without* it, which makes both `threshold` and `budget` report NOT-APPLICABLE, and
  that is how an un-passed `threshold.<lib>` gate check travelled with three published runs.
- The seven `--acknowledge-derived KEY=WHY` flags for the non-reserve landscapes that read a
  clean table with a raw counterpart, citing `docs/corpus-audit.md` §3.1–§3.2 (0 recoverable
  rows everywhere but CR9114). With them, all three corpus checks PASS.
- The row-order invariance check with a real runner. The prepared table is sorted
  **fitness-descending**, so row 0 is the global optimum and an arm whose surrogate is
  entirely NaN would post a perfect best-fitness.

---

## 8. The two honest limits

1. **It is 1,295 variants, so it can never carry the registered 125-measurement endpoint.**
   Nothing measured here is directly comparable to a 125-measurement figure on another
   landscape. It is a fast, honest, independent test at an amended budget — not a
   replacement for a landscape at 10⁴–10⁵. §9.
2. **The budget is 50 with a 10-variant initial sample**, which is thin against a 20-column
   one-hot design space: the round-1 surrogate is fitted on fewer points than it has
   features. That is a deliberate consequence of keeping the budget below the hit tier, it
   applies **identically to every arm**, and it is the kind of thing that must be stated
   with the number rather than discovered afterwards.

---

## 9. What a landscape supporting the registered protocol unchanged would require

The registered protocol needs `budget 125 < 1% tier`, i.e. **N > 12,500**, a **fixed-site**
combinatorial library (so that `LandscapeSpec`'s positions-and-parent-residues contract
holds), and d>=2 informativeness >= 5% at the 1% tier. Nothing was downloaded for this
survey.

### 9.1 The already-pinned commit cannot supply one — a definitive negative

Adding lock entries at commit `c5dccd514a9a` is explicitly permitted, so this was checked
first. CombinGym at that commit is **14 landscapes over 9 proteins**, all 10 clean workbooks
and all 9 raw tables already pinned (29/29). Their universes, from
`docs/corpus-audit.md` §3:

| still unmeasured at the pin | universe | verdict |
|---|---|---|
| `rhla_rls` / `rhla_pct_rha` / `rhla_rha` | 2,346 | far below 12,500; 0.35% saturated; nothing beyond depth 3 |
| `sacas9` | 1,296 | reserved here |
| `hiv1_protease` | 1,175 | below; 25.5% saturated, falls to 0.6% of depth 9 |
| `spcas9` | 952 (of 597 real rows) | below; 3.44% saturated, depth-correlated attrition |

Everything at or above 12,500 at that commit is already spent: CreiLOV 165,428 and GB1
149,361 and CR9114 65,536 are development data, both eqFP611 channels (8,192) are below
threshold *and* measured.

The 14 landscapes are the whole of CombinGym's DMS content at that commit: all **10** clean
workbooks and all **9** raw tables are already pinned, and the remainder of the commit's
1,087 blobs are MSAs, structures and baseline code —
`scripts/freeze_cross_protein_corpus.py` states that scope explicitly ("MSAs, structures and
baseline code are not part of this evaluation and are not pinned"), and
`docs/corpus-audit.md` §1 enumerates the raw counterpart of every clean workbook with its
git blob SHA-1. There is no unpinned DMS table left to add.

> **Verdict: no additional file at the already-pinned commit would serve.** Every landscape
> there is either below the 12,500-variant floor that the registered 1% tier requires, or
> already development data, or already measured. Growing the lock at `c5dccd514a9a` is
> permitted and would cost almost nothing — but there is nothing there to grow it with.
> **The one exception is a different phenotype on an already-pinned file, not a new file:
> `h3_mean` in `bnAbs_CR9114_raw.csv` (candidate #2 below), which needs no download and no
> lock change at all.** The pin itself must not move.

### 9.2 Ranked candidates

| # | candidate | N | independent? | cost | verdict |
|---|---|---|---|---|---|
| 1 | **New fixed-site combinatorial library from the primary literature**, new family | target 1.3×10⁴–10⁵ | **yes** | 10–18 h | **Recommended.** The only route to genuinely independent evidence. |
| 2 | **CR9114 vs H3**, from the already-pinned `bnAbs_CR9114_raw.csv` | 65,535 | **no** | 3–5 h | Cheapest by far, needs no download and no lock change — but see below. |
| 3 | **ProteinGym screen** | unknown | maybe | 2–3 h to screen | Do the screen; expect a low yield. |
| 4 | `GFP_AEQVI_Sarkisyan_2016.csv`, already on disk | 51,714 | n/a | — | **Ruled out.** |
| 5 | CR9114 vs influenza B, same raw table | 65,534 | no | — | **Ruled out:** 65,336 of 65,534 values sit at the 6.0 floor (99.7%). |

**#1 — the recommendation.** The search is precisely specified, which is what makes it
tractable: a combinatorially complete or near-complete library over a **fixed** site set,
N between 1.3×10⁴ and 10⁵, per-variant values published machine-readably, **one** named
measured phenotype (not an average), from a family **not already in this corpus** — so not
PLP fold type II (TrpB), protein G domain B1 (GB1), immunoglobulin V-domain (HER2, CR9114,
CR6261), GFP-like β-barrel (eqFP611), LOV/PAS (CreiLOV), or Cas9 (SaCas9, SpCas9). Cost:
4–8 h to find and verify candidates against §9.3, 4–8 h to prepare (prepared table,
four-way validation, digest pin, lock entry for the upstream file), 2 h for the gate and the
pre-registration. **No candidate is named here**, deliberately: two of this project's
literature claims were wrong and both were caught by reading a figure's axes rather than its
caption, and asserting a dataset I have not opened would be the same error.

**#2 — CR9114 against H3: cheap, real, and not independent.** `bnAbs_CR9114_raw.csv` is
already pinned and already on disk, and it carries `h3_mean` for 65,535 genotypes
(`docs/corpus-audit.md` §5.3) with **no** clean table curated for it, so nothing downstream
has ever read it. The 1% tier would be 655 variants and the 125-measurement budget passes
the non-saturation condition comfortably. Three things disqualify it as *generalization*
evidence, and only the third is fatal to that use:

- 58,361 of 65,535 values sit at the **6.0** assay floor — **89.1%** of the library carries
  no gradient. The floor is 6.0 and not 7.0 because the H3 titration boundary is 10⁻¹¹–10⁻⁶ M
  against H1's 10⁻¹²–10⁻⁷ M; it is a real censoring difference, not a cleaning artefact.
- It is the **same protein and the same 16-site genotype set** as CR9114-H1. Per
  `dataset-inventory.md` §2 and `HANDOFF.md` §5.6, that is one protein with two readouts —
  the eqFP611 and CR6261 caveat — not a second protein and not a second family.
- The **generator fix was designed against CR9114's measured failure**, and that fix targets
  the *library geometry* (16 binary sites), which H3 shares exactly. Rule 2.4 therefore
  reaches it: CR9114-H3 is development data for that intervention.

So it is worth preparing, but as **§5.6's paired within-library comparison across biological
readouts** (the H2 shape) — not as independent generalization evidence, and it must never be
reported as the latter. Cost 3–5 h, because the library validation is already done for
CR9114 and only the value column, the prepared table, the censoring characterisation and the
endpoint check are new.

**#3 — ProteinGym.** Screen rather than assume: download `reference_files/DMS_substitutions.csv`
(217 substitution assays) and keep only rows with `DMS_number_multiple_mutants` >= 12,500,
then for each survivor count the **distinct mutated positions** in its data file and keep
only those with a small fixed site set. ProteinGym is predominantly single-substitution DMS,
and it is already verified **not** to curate TrpB, CR9114 or CR6261
(`docs/label-provenance.md` §1.4), so the prior is low — but the screen is 2–3 hours and
answers the question definitively.

**#4 — ruled out, and it is the cautionary precedent.** `GFP_AEQVI_Sarkisyan_2016.csv` is
already present at **51,714 variants** and is **not** a fixed-site combinatorial library:
its mutations scatter across the whole protein, so it does not fit `LandscapeSpec` at all —
there is no site set to declare and no parent combo to canonicalise against. Separately,
commit `d85dca2` records that an avGFP benchmark had already **scored the model on its own
training set**. Size is not suitability, and presence on disk is not provenance.

### 9.3 What a new landscape must pass before it is trusted

Each of these is written against a defect that already cost this project a wrong result.

| # | check | the defect it is written against |
|---|---|---|
| 1 | **Raw-vs-clean diff: zero recoverable rows**, or the filter identified and the table re-prepared from raw | CR9114's clean workbook silently dropped **16,695** rows — a whole `(29, 106)` quadrant with depth-correlated attrition — and **moved a published covariate by −6.01% relative**. Recoverable rows (missing from clean, measured in raw) is the only quantity that can move a covariate. |
| 2 | **Four-way library validation** (§2): site list, varying indices, WT FASTA, and the source table's own position header — with a byte-identical depth-0 row where one exists | HER2 read the **trastuzumab light chain**. |
| 3 | **Parser anchors first** (rule 2.2): TrpB 168/2,242 = 7.4933% and GB1 72/2,167 = 3.3226% | a parser that cannot reproduce them makes nothing downstream trustworthy |
| 4 | **Row-level label provenance**: measured vs imputed, and the imputed fraction stated | TrpB carries **871 KNN-imputed** values spread across all three splits (685/100/86) |
| 5 | **No overlap with any zero-shot or pretrained model's training data**, or the exposure recorded | `d85dca2`: an avGFP benchmark **scored the model on its own training set** |
| 6 | **Censoring characterised**, floor and ceiling mass measured on the **raw** table | the inventory asserted a 6.0 bnAb floor that was false for every scored column (§5.3) |
| 7 | **`check_endpoint_viability` at the intended tier AND budget, before any arm runs** | HER2's d>=2 ceiling is 113 against 113 raw — the adjustment is a no-op and the primary endpoint silently becomes the secondary; CR6261's 19-variant tier is saturated by a 125 budget |
| 8 | **Encoding compatibility per arm, kept separate from availability** (rule 2.6) | ALDE reports `available: true` on **seven** landscapes where it returns index 0 forever |
| 9 | **Row-order invariance**, one seed, shuffled table | prepared tables are fitness-descending, so **CR9114's row 0 is the global optimum** |
| 10 | **Both files hash-pinned**: the upstream file in the corpus lock (entries added, never altered, pin never moved) and the prepared table in `PREPARED_TABLE_DIGESTS` | the lock hashed only upstream files, so the tables the arms actually read had no coverage |
| 11 | **A pre-registration committed before any arm runs**, fixing the label, tier, budget, endpoint, comparators, predictions and the falsification condition — and the landscape left unscored until then | rule 2.4 has fired twice; this document is the eleventh check applied to itself |

---

## 10. Changelog

| date | change |
|---|---|
| 2026-11-01 | Created and sealed. Landscape `sacas9_sg5` reserved and **unscored**. Library validated four ways; label fixed as sg5/`TGCCGT`; tier 5% by the protocol's own ladder; budget 50; primary endpoint diversity-adjusted top-5% yield at d>=2; `check_endpoint_viability` PASS at the amendment and FAIL at the registered protocol; predictions in §5. sg6 and sg7 declared second and third reserve. |

Any further row in this table is an amendment and must say what changed, why, and that it
was committed before the arm it affects ran.
