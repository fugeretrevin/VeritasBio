# Corpus-wide audit: is the CR9114 defect anywhere else?

**Status:** closed. Written 2026-10-31. Companion to [`docs/cr9114-reprep.md`](cr9114-reprep.md),
which named the defect in one landscape. This document rules it in or out in **all fourteen**.

**Reproduce:**

```
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\genotype_canon.py
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\corpus_raw_fetch.py fetch
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\corpus_raw_fetch.py verify
.venv\Scripts\python.exe -X utf8 -I scripts\freeze_cross_protein_corpus.py verify
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\corpus_attrition.py
```

---

## 0. The verdict, first

**Clean bill of health. The D-8 defect is unique to CR9114 and exists in no other landscape.
Not one pre-registered covariate moves. The freeze holds.**

Of the nine clean workbooks never previously diffed:

- **seven have a row-for-row identical set** to their raw counterpart — zero rows missing in
  either direction, and in six of them the values agree to the last bit;
- **two drop rows** (CR6261, SpCas9), and in both the attrition is **entirely accounted for by
  rows the raw table never measured**. The rule *"kept iff the raw measurement is present"*
  reproduces each kept set **exactly — 0 false positives, 0 false negatives**;
- **zero landscapes carry a structural exclusion**. A sweep of every single-site and
  site-pair rule of the CR9114 form — `NOT(mutated at sites ...) AND the measurement is
  present` — finds **no rule reaching 0 FP / 0 FN on any landscape but CR9114**, where it
  recovers `(29, 106)` at exactly 0 and 0;
- **no covariate can move**, for the strongest possible reason: the count of **recoverable
  rows — missing from clean yet carrying a real measurement in raw — is 0 everywhere except
  CR9114** (16,253 there). A row the assay never measured cannot enter the eligible universe
  under any convention, so it cannot shift a statistic however its absence is distributed.

Two genuine but *different* findings fell out of the audit and are recorded in §5: SpCas9's
measurement attrition is strongly depth-correlated (§5.1), and SaCas9's clean phenotype
silently averages three of the six sgRNAs the raw table reports (§5.2). Neither is D-8 class
and neither touches a pre-registered number. A third, §5.3, corrects
[`docs/dataset-inventory.md`](dataset-inventory.md) §3.

---

## 1. Provenance: every landscape has a raw table at the already-pinned commit

`scripts/data/cross_protein_corpus.lock.json` pins the ten clean workbooks and exactly **one**
raw table (`bnAbs_CR9114_raw.csv`, added by the re-preparation). The other nine clean
workbooks had never been compared against anything.

They all have a raw counterpart at the **same pinned commit**
`c5dccd514a9acfe1089a5f1555faea98afd78019`. The blob SHA-1s were read from GitHub's tree API
for the pinned tree `911f19c7ba82e112056f154d3e8b2073f70dc7d1` — the same authority the freeze
script cites for its own pins — and every download was verified against them before it was
read. A download whose hash disagrees is deleted, not kept.

| clean workbook | raw counterpart | git blob SHA-1 | size |
|---|---|---|---:|
| `bnAbs_CR9114_H1_clean.xlsx` | `Data/DMS/Raw/bnAbs_CR9114_raw.csv` | `8ac7352f…a094d` | 12,258,680 |
| `bnAbs_CR6261_H1_clean.xlsx` | `Data/DMS/Raw/bnAbs_CR6261_raw.csv` | `1158a5dc…a06c36` | 358,244 |
| `bnAbs_CR6261_H9_clean.xlsx` | *(same file — one raw table serves both antigens)* | | |
| `CreiLOV_clean.xlsx` | `Data/DMS/Raw/CreiLOV_raw.xlsx` | `300672c7…a9df0` | 22,619,510 |
| `eqFP611_clean.xlsx` | `Data/DMS/Raw/eqFP611_raw.xlsx` | `0028edb5…dd4875` | 737,089 |
| `GB1_clean.xlsx` | `Data/DMS/Raw/GB1_raw.xlsx` | `08b80e43…0c1581` | 6,633,543 |
| `HIV-1_clean.xlsx` | `Data/DMS/Raw/HIV-1_raw.tsv` | `4611f984…147ec6` | 68,024 |
| `RhlA_clean.xlsx` | `Data/DMS/Raw/RhlA_raw.xlsx` | `a8704f81…2c87a8` | 287,545 |
| `SaCas9_clean.xlsx` | `Data/DMS/Raw/SaCas9_raw.xlsx` | `f7fafde2…1cd35` | 242,262 |
| `SpCas9_clean.xlsx` | `Data/DMS/Raw/SpCas9_raw.xlsx` | `072d1c00…859776` | 112,221 |

`freeze_cross_protein_corpus.py verify` reports **CORPUS VERIFIED: 21/21**, and
`corpus_raw_fetch.py verify` reports **RAW CORPUS VERIFIED: 9/9**. **The pin is preserved:**
nothing here introduces a new upstream, and no lock file or existing pin was altered. §7
records the one recommendation that is the coordinator's call.

---

## 2. Parser validation, before any number downstream

Mutation strings omit wild-type positions, so a 16-site genotype may present as three tokens.
`scripts/analysis/genotype_canon.py` expands every genotype to exactly one character per
library position, filling the parent residue at omitted positions, and **asserts equal length
— it never zips or truncates**. `corpus_attrition.py` runs that module's gate first and
refuses to audit anything if it fails.

| landscape | universe | tier | shallow shell | expected | result |
|---|---:|---:|---|---|---|
| TrpB | 159,999 | 1,600 | 168/2,242 = **7.4933 %** | 168/2,242 = 7.4933 % | **PASS** |
| GB1 | 149,361 | 1,494 | 72/2,167 = **3.3226 %** | 72/2,167 = 3.3226 % | **PASS** |

Deep-shell rates also land on the inventory's independently recorded values (TrpB 0.9077 %,
GB1 0.9661 %), under the convention that depth 0 is excluded from both shells and retained in
the universe.

**Three further validations** make the audit's recomputations trustworthy rather than merely
self-consistent. The same code, run over the clean tables, reproduces:

1. **every frozen covariate in [`docs/preregistration-cross-protein.md`](preregistration-cross-protein.md) §6.5**
   — CreiLOV 1.8753x, blue 1.8541x, combined 1.8297x, red 0.6343x, CR9114 0.5845x, CR6261-H1
   0.2102x, CR6261-H9 0.0000x — with the shell sizes and hit counts matching digit for digit
   (62,790/1,178; 4,095/76; 4,095/75; 4,095/26; 38,835/227; 945/4; 917/0);
2. **CR9114's filter, independently re-derived**: the composite sweep recovers sites
   `(29, 106)` at **0 FP and 0 FN**, which is `docs/cr9114-reprep.md` §2.2's result reached by
   a differently written sweep;
3. **the reweighting bound at 0.5622x** on CR9114, the figure `cr9114-reprep.md` §6 records,
   against a truth of 0.5845x.

A pipeline that reproduces the known defect and all seven frozen covariates is a pipeline
whose *negative* findings mean something.

---

## 3. The result, per landscape

Design space is computed exactly: `C(n,k)` for a binary library, and for the mixed alphabets
here the coefficient of `x^k` in `prod_i (1 + (a_i - 1) x)`, which is the proper
generalisation and not an approximation to it.

| landscape | level | raw | clean | diff | cause | depth-corr. | covariate |
|---|---|---:|---:|---:|---|---|---|
| `creilov` | final test | 165,428 | 165,428 | **0** | none — identical sets | — | **does not move** |
| `eqfp611_blue` | final test | 8,192 | 8,192 | **0** | none — identical sets | — | **does not move** |
| `eqfp611_red` | final test | 8,192 | 8,192 | **0** | none — identical sets | — | **does not move** |
| `eqfp611_combined` | final test | 8,192 | 8,192 | **0** | none — identical sets | — | **does not move** |
| `cr9114_h1` | final test | 65,536 | 48,841 | 16,695 | **D-8 structural exclusion** | **YES** | **moved −0.0374** (already fixed) |
| `cr6261_h1` | final test | 1,917 | 1,887 | 30 | absent measurement (0 FP, 0 FN) | yes, **vacuous** | **does not move** |
| `cr6261_h9` | final test | 1,917 | 1,842 | 75 | absent measurement (0 FP, 0 FN) | no | **does not move** |
| `hiv1_protease` | reserve | 1,175 | 1,175 | **0** | none — identical sets | — | not pre-registered |
| `rhla_rls` | reserve | 2,346 | 2,346 | **0** | none — identical sets | — | not pre-registered |
| `rhla_pct_rha` | reserve | 2,346 | 2,346 | **0** | none — identical sets | — | not pre-registered |
| `rhla_rha` | reserve | 2,346 | 2,346 | **0** | none — identical sets | — | not pre-registered |
| `sacas9` | reserve | 1,296 | 1,296 | **0** | none — identical sets | — | not pre-registered |
| `spcas9` | reserve | 952 | 597 | 355 | absent measurement (0 FP, 0 FN) | yes, **vacuous** | not pre-registered |
| `gb1_combingym` | development | 149,361 | 149,361 | **0** | none — identical sets | — | not pre-registered |

**`clean − raw = 0` on every single landscape.** No clean table contains a genotype its raw
counterpart does not, which would have been just as alarming as the reverse.

### 3.1 The seven identical landscapes

| landscape | rows | design space | saturation | max \|raw − clean\| |
|---|---:|---:|---:|---:|
| CreiLOV | 165,428 | 184,320 (mixed 2·2·2·3·3·2·5·2⁸) | 89.75 % | **0** |
| eqFP611 blue / red / combined | 8,192 | 8,192 (2¹³) | **100.00 %** | **0** |
| GB1 (CombinGym) | 149,361 | 160,000 (20⁴) | 93.35 % | **0** |
| HIV-1 protease | 1,175 | 4,608 (mixed) | 25.50 % | **0** |
| RhlA ×3 | 2,346 | 663,552 (mixed) | 0.35 % | **0** |
| SaCas9 | 1,296 | 1,296 (mixed 2·2·3·3·3·2·2·3) | **100.00 %** | 1.1e−16 |

eqFP611 and SaCas9 are **combinatorially complete** — every genotype in the design space is
present, at every depth, in both tables. There is no room for a filter to hide.

Where saturation is below 100 % the shortfall is in the **raw table itself** — library
synthesis and sequencing dropout — and is identical in clean and raw, so it is not a cleaning
artefact and cannot have been introduced after the freeze. It is recorded because it is a real
limit on these landscapes, and in two cases it is severe and depth-structured: **HIV-1 falls
from 100 % of depth-1 to 0.6 % of depth-9**, and **RhlA holds nothing at all beyond depth 3**
(which is why §6.3's degeneracy guard already declares the covariate undefined on all three
RhlA landscapes). CreiLOV declines gently, 100 % at depth 1 to 62.5 % at depth 15. **None of
this is attrition the cleaning step caused, and none of it moves a frozen number** — the
covariate was computed from the clean table, which is the raw table.

### 3.2 CR6261 H1 and H9 — benign, and provably unable to move

The two highest-stakes unchecked landscapes, both pre-registered final-test targets.

The raw `bnAbs_CR6261_raw.csv` holds 1,917 of 2,048 genotypes and serves **both** antigens. Its
`genotype` column is an 11-bit site indicator and `som_mut` its popcount; both were asserted
consistent on every row, and the clean workbooks' 121-residue sequences map onto the same 11
library positions (28, 30, 58, 59, 62, 74, 75, 76, 77, 79, 104; parent `TSANQKSTSAL`).

| | H1 | H9 |
|---|---:|---:|
| raw genotypes | 1,917 | 1,917 |
| raw rows with a measurement | **1,887** | **1,842** |
| clean rows | **1,887** | **1,842** |
| missing from clean | 30 | 75 |
| of which **unmeasured in raw** | **30 (all)** | **75 (all)** |
| **recoverable rows** | **0** | **0** |
| rule *kept iff measured* | **0 FP, 0 FN** | **0 FP, 0 FN** |
| best composite D-8-form rule | (74, 75): 423 FP, 0 FN | (75, 77): 407 FP, 0 FN |
| value fidelity on shared rows | 9.77e−15 | 9.59e−14 |
| enrichment from clean | **0.2102x** | **0.0000x** |
| enrichment over the raw universe | **0.2102x** | **0.0000x** |

**The measured sets are identical, so the covariates are identical.** Not "close" — the same
rows, the same tier, the same shell, the same number. `0.210x` and `0.000x` stand exactly as
frozen.

H1's membership difference does trip the crude depth-correlation threshold, at an 18.2 pp
spread, but the whole of it is **2 rows out of 11 at depth 1**, and both are rows the assay
never quantified. With 0 recoverable rows the depth profile of the difference is **vacuous**:
there is nothing to reweight. The best site-pair rule misclassifies 423 kept genotypes; the
quadrant signature that gives CR9114 away — one pair at 0 error, the next-best off by
thousands — is absent.

### 3.3 CR9114 — the negative control, reproduced

The audit re-derives the known defect from scratch: 65,536 raw against 48,841 clean, 16,695
missing, **16,253 of them carrying a real measurement**, retention falling 100 % → 12.5 % → 0 %,
and the composite rule at sites **(29, 106)** reproducing the kept set at **0 FP and 0 FN**
while the next-best rule is off by 16,139. The covariate moves 0.6219x → 0.5845x, −6.01 %
relative. All of this matches `docs/cr9114-reprep.md`. **No action — already corrected.**

---

## 4. Would a covariate move? The decisive quantity

The audit reports a field the CR9114 investigation did not need, because there it was large:
**recoverable rows**, meaning rows missing from clean that carry a real measurement in raw.

This is the only quantity that can move a covariate. A missing row the assay never measured
could not have entered the eligible universe under any convention — the freeze script drops
blank and unparseable values before it computes anything — so its absence cannot shift a
statistic no matter how depth-correlated that absence looks.

| landscape | recoverable rows | covariate can move? |
|---|---:|---|
| `cr9114_h1` | **16,253** | yes — and it did, by −6.01 % relative |
| every other landscape | **0** | **no** |

Recomputing each covariate over the raw universe with identical code confirms it: a delta of
**exactly 0.0000** everywhere but CR9114.

### 4.1 The reweighting bound, as a fallback and as a check

Where no raw table can be obtained, the fallback is the depth-blind reweighting bound: hold
each depth's observed hit rate fixed and reweight to the exact design-space population. **Both
the shell and the base rate must be reweighted** — a table that under-represents deep shells
also understates the overall hit rate whenever hit rate varies with depth, and dividing a
reweighted shell rate by the table's own unreweighted base rate is what makes the bound look
falsely optimistic.

Implemented that way it reproduces `cr9114-reprep.md` §6's **0.5622x** against a truth of
0.5845x — the right direction, 3.8 % relative. It was not needed as a fallback anywhere: a raw
table was obtained for every landscape. Its values on the complete landscapes (eqFP611,
SaCas9) equal the measured enrichment exactly, as they must when saturation is 100 %, which is
a further check on the implementation.

---

## 5. Three findings that are real but not D-8 class

### 5.1 SpCas9's measurement attrition IS depth-correlated

The clean `SpCas9` sheet holds 5,748 rows of which **5,151 are wholly empty padding**; the real
content is **597 rows** (596 variants plus the parent). The raw table holds **952**.

The 355 missing rows are explained exactly: **every one has a blank or `NA` value in at least
one of the two on-target columns the clean table keeps** (`RFPsg5 ON`, `RFPsg8 ON`), and **not
one kept row has either missing**. The rule is `RFPsg5 ON present AND RFPsg8 ON present`, at
**0 FP and 0 FN**. The clean `mean` was confirmed to equal `mean(RFPsg5 ON, RFPsg8 ON)` on all
597 rows, and `target_RFPsg5` to equal `RFPsg5 ON` exactly.

But the retention is strongly depth-correlated — **95.8 % at depth 2 falling monotonically to
47.2 % at depth 8**, a 48.6 pp spread:

| depth | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| raw | 7 | 24 | 64 | 148 | 252 | 267 | 153 | 36 |
| clean | 7 | 23 | 54 | 115 | 162 | 138 | 80 | 17 |
| retention | 100 % | 95.8 % | 84.4 % | 77.7 % | 64.3 % | 51.7 % | 52.3 % | 47.2 % |

**Why this is not a D-8-class defect, and what it still costs.** The cause is absent
measurement, not an undocumented structural exclusion: there is no quadrant, no value-blind
rule, and the 355 rows carry **no recoverable signal** — 0 of them have both on-target values,
so none could enter the universe. SpCas9 carries no pre-registered covariate (it is `reserve`),
and its enrichment over the raw universe is identical to the clean one, 2.5000x, because the
measured sets coincide.

What it does mean is that **SpCas9's per-depth statistics describe a depth-biased sample of its
library, and the bias is in the assay rather than the curation.** Deeper SpCas9 variants were
likelier to fail one of the two on-target readouts. If SpCas9 is ever promoted from reserve to
a scored target, this belongs in its write-up, and the depth bias should be stated alongside
the number rather than discovered afterwards. Note also that SpCas9's library is only **3.44 %
saturated** (952 of 27,648) and is a designed subset, not a combinatorially complete library,
so per-depth rates there are weakly grounded regardless.

### 5.2 SaCas9's clean phenotype averages three of six sgRNAs, undocumented

SaCas9's rows are perfect: 1,296 raw against 1,296 clean, zero difference, and a
**combinatorially complete** library. But the raw table reports **six** sgRNAs — sg1, sg2, sg3,
sg5, sg6, sg7 — and the clean workbook keeps only sg1–sg3. Its `Mean` column was verified to
equal `mean(sg1, sg2, sg3)` on all 1,296 rows, and to differ from the mean over all six on
**all 1,296**.

This is a **phenotype-definition** choice, not row attrition: no genotype is lost, no depth is
touched, and no covariate can move. It is recorded because the clean table gives no indication
that half the measured sgRNA panel was dropped, and anyone comparing SaCas9 against the
upstream paper will not reproduce its numbers without knowing. The same shape of column
reduction appears in SpCas9, where the clean table drops the two off-target readouts
(`RFPsg5 OFF5-2`, `RFPsg8 OFF5`) — there at least the kept columns are named for what they are.

### 5.3 The `dataset-inventory.md` §3 floor caveat is wrong, and now corrected

The inventory stated that the raw bnAb tables encode non-binders at **6.0** and that
*"CombinGym's cleaning step therefore changes the censoring value"*. **That is false for every
scored bnAb column, on both raw tables.** Measured directly, column by column:

| raw table | column | n | min | rows at 6.0 | rows at 7.0 | scored? |
|---|---|---:|---:|---:|---:|---|
| `bnAbs_CR9114_raw.csv` | `h1_mean` | 65,094 | **7.0** | **0** | **1,675** | yes |
| | `h3_mean` | 65,535 | 6.0 | 58,361 | 0 | no — no clean table |
| | `fluB_mean` | 65,534 | 6.0 | 65,336 | 0 | no — no clean table |
| `bnAbs_CR6261_raw.csv` | `h1_mean` | 1,887 | **7.0** | **0** | **214** | yes |
| | `h9_mean` | 1,842 | **7.0** | **0** | **372** | yes |

Three corrections follow:

1. **No cleaning step changes any censoring value, on any of the three bnAb landscapes.** Raw
   and clean agree on both the floor (7.0) and the mass at it: CR6261-H1 214 and 214,
   CR6261-H9 372 and 372. CR9114's mass differs (1,675 raw against 934 clean) **only because
   the excluded quadrant was modestly enriched for floored rows** — the floor *value* is 7.0 in
   both, and that gap is now closed by the re-preparation.
2. **The 6.0 floor applies to `h3_mean` and `fluB_mean` only** — exactly as the paper's
   titration boundaries require (10⁻¹¹–10⁻⁶ M for H3 and influenza B, 10⁻¹²–10⁻⁷ M for H1 and
   H9). Neither column is curated into a clean table and neither is scored anywhere.
3. **The plural `bnAbs_*_raw.csv` is wrong twice over**: `bnAbs_CR6261_raw.csv` has **no**
   `h3_mean` or `fluB_mean` column at all, so there is no 6.0 anywhere in it. The claim had no
   support in that file whatsoever.

The caveat has been narrowed in `dataset-inventory.md` §3 accordingly. The inventory's
CR6261 masses of 213 and 371 are also correct as written: they count the prepared CSVs, which
drop the germline row, and CR6261's germline sits at the floor.

---

## 6. What this means for the freeze

**Nothing moves. The pre-registration stands as frozen**, with CR9114's already-amended 0.5845x
and no further change:

| landscape | §6.5 frozen | audited | verdict |
|---|---:|---:|---|
| `creilov` | 1.875x | 1.8753x | **unchanged** |
| `eqfp611_blue` | 1.854x | 1.8541x | **unchanged** |
| `eqfp611_combined` | 1.830x | 1.8297x | **unchanged** |
| `eqfp611_red` | 0.634x | 0.6343x | **unchanged** |
| `cr9114_h1` | 0.584x | 0.5845x | **unchanged** (amended before any arm ran) |
| `cr6261_h1` | 0.210x | 0.2102x | **unchanged** |
| `cr6261_h9` | 0.000x | 0.0000x | **unchanged** |

The per-protein ordering the primary test uses — CreiLOV 1.875x > eqFP611 1.439x > CR9114
0.584x > CR6261 0.105x — is untouched, and so is the §6.5 prediction it forces.

The honest summary of the exercise: **CR9114 was the anomaly, not the first instance of a
pattern.** It is the only landscape CombinGym curates from a derived workbook that silently
disagrees with its own raw table about which genotypes exist. That was worth proving rather
than assuming, and the cost of proving it was one afternoon and nine hash-verified downloads.

---

## 7. Recommendations (the coordinator's call — nothing here was acted on)

1. **Extend the lock to the nine raw tables.** They are now the evidence that the clean tables
   are faithful, and that evidence should be hash-pinned like everything else. The blob SHA-1s
   are in §1 and in `scripts/analysis/corpus_raw_fetch.py`. This only ever adds entries; it
   alters none. `freeze_cross_protein_corpus.py verify` already guards against a landscape
   being repointed at a file the lock does not cover, which is the drift the CR9114
   re-preparation exposed — that guard is why this is worth doing rather than optional.
2. **Record SpCas9's depth-correlated measurement attrition (§5.1) before SpCas9 is ever
   promoted** from reserve to a scored target. No landscape needs re-preparing for it.
3. **Note SaCas9's three-of-six sgRNA phenotype (§5.2)** in the landscape registry comment,
   next to the CR9114 note that already lives there.
4. **No landscape needs re-preparing.** This was the outcome worth hoping for, and the numbers
   in §3 are what justify saying so.

---

## 8. What was produced

| artefact | role |
|---|---|
| `scripts/analysis/genotype_canon.py` | canonicalisation with an asserted length invariant, exact mixed-alphabet design-space expectations, shell statistics, the reweighting bound, and the TrpB/GB1 validation gate |
| `scripts/analysis/corpus_raw_fetch.py` | fetches and hash-verifies the nine raw tables at the pinned commit; deletes anything that fails |
| `scripts/analysis/corpus_attrition.py` | the audit: set differences, cause attribution, depth-correlation, the single-site and site-pair sweeps, and the covariate recomputation |
| `scripts/analysis/cr9114_attrition.py` | pre-existing; left unchanged. Its output is now the negative control described in `cr9114-reprep.md` §8 |

`corpus_attrition.py` supersedes `cr9114_attrition.py` in scope: it covers all fourteen
landscapes rather than three, computes the exact design space for mixed alphabets rather than
assuming binary sites, and diffs against the raw table rather than against a binomial alone.
