# CR9114 H1: the 25% attrition, named exactly, and the landscape re-prepared

**Status:** closed. Written 2026-10-31. Closes open item **5** of
[`docs/label-provenance.md`](label-provenance.md) §5 and the unresolved discrepancy at its §3.5.

**Reproduce:**

```
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_cr9114_h1.py validate
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_cr9114_h1.py prepare --source <elife-71393-fig1-data1-v3.csv>
.venv\Scripts\python.exe -X utf8 -I scripts\data\prepare_cr9114_h1.py stats
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\cr9114_attrition.py
```

---

## 0. The verdicts, first

1. **The filter is named, exactly, with zero residual.** It was **not** a binder threshold, not a
   both-antigen requirement, not a confidence or replicate filter, and not the 6.0-vs-7.0 floor
   divergence. CombinGym's clean workbook excludes **an entire quadrant of the design space**:
   every genotype carrying **both** mutation 1 (`F29S`) **and** mutation 16 (`Y106S`) — all
   16,384 of them, regardless of affinity — plus the 311 unquantified rows outside that quadrant.
   `65,536 − 16,384 − 311 = 48,841` clean rows, i.e. our 48,840 mutants. The rule reproduces the
   kept set **exactly**: 0 false positives, 0 false negatives.
2. **The full table was obtained**, and it is not a new upstream. Phillips *et al.*'s deposited
   Figure 1 — source data 1 holds all 65,536 genotypes and is **byte-identical** to
   `Data/DMS/Raw/bnAbs_CR9114_raw.csv` at the already-pinned CombinGym commit. **The corpus pin is
   preserved.** We were reading the derived clean workbook where the raw table was already on disk.
3. **Retention is now flat across depth: 98.3 % to 100.0 %**, against 100.0 % → 12.5 % before.
   `cr9114_attrition.py`'s verdict flips from *DEPTH-CORRELATED* to *approximately depth-blind*.
4. **The pre-registered 0.622x shallow enrichment was an overestimate, as suspected. The true
   value is 0.5845x** — a 6.0 % relative reduction. The hit count moves 489 → 651 and the d≥2
   ceiling 280 → 365.
5. **CR9114 can now carry a pre-registered per-depth ranking**, with two caveats that this work
   does not touch (§7).

---

## 1. Parser validation, before anything downstream

Mutation strings omit wild-type positions, so a 16-site genotype may present as three tokens.
Genotypes are canonicalised to full library length with wild type filled in at omitted positions,
and **equal length is asserted, never zipped or truncated** — the bug that inverted a conclusion
earlier in this project.

The implementation was validated by reproducing the two pre-registered reference figures from the
committed TrpB and GB1 files *before* any CR9114 number was computed, under the convention that
depth 0 is excluded from both shells but retained in the eligible universe:

| landscape | universe | tier | shallow shell | expected | result |
|---|---:|---:|---|---|---|
| TrpB | 159,999 | 1,600 | 168/2,242 = **7.4933 %** | 168/2,242 = 7.4933 % | **PASS** |
| GB1 | 149,361 | 1,494 | 72/2,167 = **3.3226 %** | 72/2,167 = 3.3226 % | **PASS** |

Deep-shell rates also match the inventory's independently recorded values (TrpB 0.9077 %,
GB1 0.9661 %). As a third check, the same code reproduces **every** committed CR9114 figure from
the old file — 489 hits, threshold 9.6077550213894, 0.6226 % shallow rate, 0.6218x enrichment,
d≥2 ceiling 280, d≥3 ceiling 118, mean hit depth 8.9100 — and its whole depth profile, digit for
digit against [`docs/dataset-inventory.md`](dataset-inventory.md) §4b and §4c.

---

## 2. What the 48,840-row file actually was

### 2.1 Provenance chain

`scripts/data/CR9114_H1_Phillips_2021.csv` was derived from CombinGym
`Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx`, sheet `CR91114-h1`, column `h1_mean`, at commit
`c5dccd514a9acfe1089a5f1555faea98afd78019` — 48,841 rows, less the wild-type row.

The clean workbook keys rows by the **121-residue protein sequence**, not by Phillips' 16-bit
genotype string, which is why the two tables could not previously be lined up. They were mapped by
recovering the 16 varying sequence indices from the clean table itself and cross-checking every
parent residue against `Data/WT Sequence/bnAbs_CR9114.fasta`:

| | |
|---|---|
| library positions (1-based) | 29, 30, 31, 52, 57, 58, 59, 71, 74, 75, 76, 77, 84, 87, 95, 106 |
| parent (germline) residues | `FSSITANTKSTSSRYY` |
| mutant residues | `SNNSSTASIFSNNTFS` |
| germline row | present, `n_mut = 0`, `h1_mean` 8.424622015254, sequence identical to the FASTA |

Under that mapping all 48,841 clean rows map to distinct bitstrings, and **all 48,841 `h1_mean`
values match the raw table exactly**. These positions and parent residues are identical to those
already hard-coded in `veritas/calibration/landscapes.py`, so nothing downstream needs to change.

### 2.2 The filter, named

Per-site retention is the tell. If the filter were site-blind every site would be kept at the
overall 74.5 %:

| site | seq. position | mutation | kept when bit = 1 |
|---:|---:|---|---|
| **1** | 29 | F→S | **16,139 / 32,768 = 49.25 %** |
| 2–15 | 30 … 95 | — | 24,390–24,576 / 32,768 = **74.43 %–75.00 %** |
| **16** | 106 | Y→S | **16,374 / 32,768 = 49.97 %** |

Two sites, and only those two, are depleted by a third. The joint contingency closes it:

| site 1 | site 16 | genotypes | kept | kept % | of which `h1_mean` blank |
|---|---|---:|---:|---:|---:|
| 0 | 0 | 16,384 | 16,328 | 99.66 % | 56 |
| 0 | 1 | 16,384 | 16,374 | 99.94 % | 10 |
| 1 | 0 | 16,384 | 16,139 | 98.50 % | 245 |
| **1** | **1** | **16,384** | **0** | **0.00 %** | 131 |

> **The rule is `NOT(site 1 AND site 16) AND h1_mean is present`.**
> It predicts **48,841** rows and matches CombinGym's kept set **exactly — 0 FP, 0 FN.**

A brute-force sweep over all 120 site pairs confirms the pair is unique: `(1, 16)` gives symmetric
difference **0**; the next-best pair, `(1, 6)`, is off by 16,177 rows.

### 2.3 Why that filter is mechanically depth-correlated

The fraction of depth-*k* genotypes carrying two *specified* mutations is

```
C(14, k-2) / C(16, k)  =  k(k-1) / (16 x 15)  =  k(k-1) / 240
```

so retention is `1 − k(k−1)/240`, which falls monotonically and quadratically in depth. That is
the observed curve, not an approximation to it:

| depth | predicted `1 − k(k−1)/240` | observed retention |
|---:|---:|---:|
| 1 | 100.0 % | 100.0 % |
| 4 | 95.0 % | 93.7 % |
| 8 | 76.7 % | 76.3 % |
| 12 | **45.0 %** | **45.0 %** |
| 15 | **12.5 %** | **12.5 %** |
| 16 | 0 % | 0 % |

The small shortfalls at low depth are the 442 unquantified rows. So the entire "25 % unexplained
attrition" was one two-site exclusion, and its depth signature is arithmetic.

### 2.4 What the filter was *not* — four hypotheses, all falsified

- **Not a binder / non-binder threshold.** 934 *kept* rows sit exactly at the 7.0 assay floor, and
  741 *dropped* rows do too. Both sides contain floored non-binders. The dropped set reaches
  `h1_mean` 9.7936, deep inside the top tier.
- **Not a both-antigen or cross-antigen requirement.** The clean workbook carries no H3 or
  influenza-B column at all, and 58,361 raw rows sit at the H3 floor while most are kept.
- **Not a confidence or replicate filter.** Dropped rows are, if anything, *more* precise:
  `h1_sem` median 0.0262 dropped against 0.0361 kept (p90 0.0794 vs 0.0772). And 672 kept rows rest
  on a **single** replicate, while 15,041 dropped rows have all three.
- **Not the 6.0-vs-7.0 floor divergence.** This correction matters:
  > [`docs/dataset-inventory.md`](dataset-inventory.md) §3 states that *"the raw tables
  > (`Data/DMS/Raw/bnAbs_*_raw.csv`) encode them at 6.0"* and that *"CombinGym's cleaning step
  > therefore changes the censoring value"*. **That is false for H1.** The raw `h1_mean` column's
  > minimum is **7.0**, with 1,675 rows at the floor — the same convention as the clean table. The
  > 6.0 floor applies to `h3_mean` and `fluB_mean` only, exactly as the paper's titration
  > boundaries require (10⁻¹¹–10⁻⁶ M for H3 and influenza B, 10⁻¹²–10⁻⁷ M for H1 and H9). No floor
  > reconciliation was needed, and the caveat should be narrowed to H3/FluB.

**16,253 of the 16,695 dropped rows carry a genuine, measured H1 affinity.** The cleaning step
discarded real measurements, not low-quality ones. The mechanism inside CombinGym is undocumented
and not recoverable from the published artefacts; what is established here is the rule it applied,
and that the rule is value-blind.

---

## 3. The authoritative table

| field | value |
|---|---|
| citation | Phillips A.M. *et al.* (2021) *eLife* **10**:e71393, Figure 1 — source data 1 |
| title | "CR9114 library -logK_D to H1, H3, and influenza B" |
| URL | `https://cdn.elifesciences.org/articles/71393/elife-71393-fig1-data1-v3.csv` |
| size | 12,258,680 bytes |
| sha256 | `ee31f8fc26fce2730d2ae772e7a095aa334b4c6c5d8e7fbfe4b4a5a1c9a794fe` |
| rows | **65,536** — the complete 2¹⁶ design space, one row per genotype |
| columns | `genotype`, 3 replicates + mean + sem for each of H1/H3/fluB, `pos1…pos16`, `som_mut` |

The asset list was taken from eLife's own API (`https://api.elifesciences.org/articles/71393`)
rather than the article HTML, which a summarizing fetch tool refuses on this article.

**The pin is preserved.** This file is **byte-identical** (same sha256, same 12,258,680 bytes) to
CombinGym's `Data/DMS/Raw/bnAbs_CR9114_raw.csv` at the pinned commit. The re-preparation reads a
file the pinned corpus already ships; it does not introduce a new upstream.

### Reconciling with the paper's `N_var`

- 65,536 genotypes; **442** have a blank `h1_mean` (not quantified) → **65,094 measured**.
- Of those, **1,675** sit at the 7.0 floor → **63,419 binders**.
- The paper states **63,419** H1 binders (97 %). **The binder count reconciles exactly.**
- The paper's `N_var = 65,091` is 3 fewer than 65,094; the gap sits entirely in the non-binder
  count (paper 1,672 vs source data 1,675). The deposited source data is authoritative and is what
  is used here. This 3-row discrepancy is noted, not resolved.

---

## 4. The re-prepared file

| | old | new |
|---|---|---|
| rows | 48,840 | **65,093** |
| bytes | 3,142,855 | 4,296,350 |
| sha256 | — | `1d9a57745636f27d5fdea5f8798a852f4852b9dc7fa2e60884d49bd930972d83` |

65,093 = 65,094 measured − 1 wild-type row (dropped, per the repo convention that the parent is not
a library member).

Conventions match the committed CR6261 files exactly, all verified by assertion:

- columns `mutant,DMS_score,DMS_score_bin,stage`; `mutated_sequence` **dropped** (the wild type
  lives in `library_wt_sequences.json` and still matches the CombinGym FASTA byte for byte);
- `mutant` in ProteinGym notation over the 16 library positions, wild-type positions **omitted**,
  tokens in ascending position order (0 violations in 65,093 rows);
- `DMS_score_bin` = 1 iff `DMS_score >= 8.424622015254`, the germline value (0 violations) — the
  same rule the old file used;
- `stage` = `train` on every row; rows sorted by `DMS_score` descending, ties broken by `mutant`.

**Fidelity against the old file:**

- **0 of the 48,840 old rows were lost.** 16,253 rows recovered.
- All 65,093 `DMS_score` strings are **byte-identical to the published eLife source strings**. The
  old file carried `.xlsx`-rounded 15-digit values; max deviation 9.77 × 10⁻¹⁵, i.e. float
  round-trip only. The new file is strictly more faithful.
- **The fully-matured antibody is back.** The depth-16 genotype — all 16 somatic mutations, the
  actual CR9114 — was **entirely absent** from the old landscape, since it necessarily carries both
  excluded mutations. It is present at `DMS_score` 9.592294332007965. The inventory's remark that
  *"CR9114's maximum observed depth is 15 … the fully-matured all-mutation genotype is absent from
  the curated table"* is now obsolete.

`33/33` tests in `tests/test_calibration_landscapes.py` and `tests/test_calibration_genotypes.py`
pass against the new file, unchanged.

---

## 5. The new retention curve — flat

```
  CR9114_H1   16 binary sites, design space 65,536
    our rows 65,093   paper reports 65,091   at floor 7.0: 1,675
    depth    expected  observed  retention            (was)
    1              16        16    100.0%            100.0%
    2             120       118     98.3%             97.5%
    3             560       552     98.6%             96.1%
    4           1,820     1,796     98.7%             93.7%
    5           4,368     4,302     98.5%             90.3%
    6           8,008     7,915     98.8%             86.6%
    7          11,440    11,339     99.1%             81.9%
    8          12,870    12,797     99.4%             76.3%
    9          11,440    11,394     99.6%             69.8%
    10          8,008     7,984     99.7%             62.5%
    11          4,368     4,363     99.9%             54.1%
    12          1,820     1,820    100.0%             45.0%
    13            560       560    100.0%             35.0%
    14            120       120    100.0%             24.2%
    15             16        16    100.0%             12.5%
    16              1         1    100.0%              0.0%
    retention range over depths with >=10 expected: 98.3% to 100.0%
    => approximately depth-blind - shell statistics unbiased, counts only shrink
```

**Retention is flat: a 1.7-percentage-point spread, against 87.5 points before.** The residual
443 absent rows are the 442 genuinely unquantified genotypes plus the dropped wild type. What
remains is not merely flatter but **non-monotone in the opposite direction** — slightly *more*
missingness at low depth (98.3 % at depth 2) than high (100.0 % at depths 12–16) — so it cannot
manufacture a shallow-to-deep gradient. Depths 9 through 16 are complete or within 0.4 %.

---

## 6. Recomputed statistics, and what moves

Identical code and conventions on both files; the old file's column is a verified reproduction of
the committed values.

| quantity | committed (48,840) | **re-prepared (65,093)** | change |
|---|---:|---:|---|
| variant count | 48,840 | **65,093** | **+16,253 (+33.3 %)** |
| top-1 % tier size = hit count | **489** | **651** | **+162 (+33.1 %)** |
| top-1 % threshold | 9.6077550213894 | **9.6088080625474** | +0.00105 |
| base rate | 1.0012 % | 1.0001 % | — |
| shallow shell, `depth ≤ 8` | 202/32,444 = 0.6226 % | **227/38,835 = 0.5845 %** | −0.0381 pp |
| **enrichment, `depth ≤ 8` (primary)** | **0.6218x** | **0.5845x** | **−0.0373 (−6.0 % rel.)** |
| deep shell, `depth > 8` | 287/16,396 = 1.7504 % | 424/26,258 = 1.6147 % | −0.136 pp |
| shallow shell, fixed `depth ≤ 2` | 0/133 = 0 % | **0/134 = 0 %** | **unchanged** |
| **enrichment, fixed `depth ≤ 2`** | **0.0000x** | **0.0000x** | **unchanged** |
| mean depth of a top-1 % hit | 8.9100 | **9.2657** | **+0.3557** |
| max depth observed | 15 | **16** | **+1** |
| diversity-adjusted ceiling `d ≥ 2` | **280** | **365** | **+85 (+30.4 %)** |
| diversity-adjusted ceiling `d ≥ 3` | 118 | **143** | +25 (+21.2 %) |
| budget share at 125 measurements | 0.26 % | **0.192 %** | −0.068 pp |

**Per-depth hit rate** (the quantity the filter was biasing):

| depth | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| new | 0/16 | 0/118 | 1/552 | 3/1,796 | 10/4,302 | 39/7,915 | 60/11,339 | 114/12,797 | 139/11,394 | 110/7,984 | 89/4,363 | 57/1,820 | 20/560 | 8/120 | 1/16 | 0/1 |
| rate | 0 % | 0 % | 0.18 % | 0.17 % | 0.23 % | 0.49 % | 0.53 % | 0.89 % | 1.22 % | 1.38 % | 2.04 % | 3.13 % | 3.57 % | 6.67 % | 6.25 % | 0 % |

### Which pre-registered values move

**Move, and must be re-recorded:** the universe (48,841 → 65,093), the hit count
(**489 → 651**), the top-1 % threshold, the primary shallow enrichment
(**0.622x → 0.5845x**), the deep-shell rate (1.750 % → 1.615 %), the mean hit depth
(8.910 → 9.266), the max observed depth (15 → 16), both diversity-adjusted ceilings
(**d≥2 280 → 365**, d≥3 118 → 143) and the budget share.

**Does not move:** the fixed `depth ≤ 2` cut stays at **0 hits and 0.0000x enrichment**, and stays
below the 500-variant resolvability floor (shell of 134). The *sign and qualitative shape* of the
finding also do not move: shallow enrichment remains far below 1 and the hit rate still rises
monotonically with depth. **CR9114's inverted depth structure was real, not an artefact of the
filter — but its magnitude was overstated by 6.0 %.**

### Cross-check: the reweighting bound predicted this

The bias bound available *without* the full table — hold each depth's observed hit rate fixed and
reweight the shells to the exact C(16,k) population — gives:

| | enrichment, `depth ≤ 8` |
|---|---:|
| as committed (biased) | 0.6218x |
| depth-blind reweighting **bound** | **0.5622x** |
| **truth, from the full table** | **0.5845x** |

The bound called the direction correctly (the committed value was an overestimate) and landed
within 0.022, 3.8 % relative. That the two independent routes agree is the strongest available
check that the re-preparation did not introduce a new error.

---

## 7. Can CR9114 carry a pre-registered ranking?

> **Yes.** Every depth shell is now ≥ 98.3 % complete, the residual incompleteness is not
> depth-monotone and runs the opposite way, and no per-depth statistic is shaped by a filter any
> more. The covariate ordering may use CR9114's primary-cut enrichment of **0.5845x**, replacing
> 0.622x.

Two caveats survive this work unchanged, and neither is a data-quality problem:

1. **Report the primary `depth ≤ ⌊16/2⌋ = 8` cut, not the fixed `depth ≤ 2` cut.** The fixed shell
   holds 134 variants and 0 hits; a 1 %-order rate is not estimable there, and `0.0000x` carries no
   information. This was already the pre-registered recommendation and it still binds.
2. **CR9114's deep-shell structure is anticipated prior knowledge, not a prediction.** Phillips
   Figure 5A publishes max affinity against somatic-mutation count, and the paper quantifies the
   nested requirement structure. Per [`docs/label-provenance.md`](label-provenance.md) §3.4, risk
   7.6 is **real for CR9114** and the H1 ranking must be written as *anticipated*. Fixing the
   attrition makes the number trustworthy; it does not make it novel. (The genuinely
   non-anticipated result remains CR6261-H9, which runs *contrary* to the paper's prose.)

Also worth retaining: 1,675 rows (2.57 %) are censored at the 7.0 assay floor, slightly up from
1.91 % in the old file, because the excluded quadrant was modestly enriched for floored rows.

---

## 8. Recorded values that are now stale

**No lock file was altered.** `scripts/data/cross_protein_corpus.lock.json` pins only the 20
upstream CombinGym blobs (10 `.xlsx`, 9 `.fasta`, 1 summary) and records **no hash of any derived
CSV**, so re-preparing the CSV changes nothing in it. Verified: no entry's path contains `csv`.

| artefact | stale value | should be | action |
|---|---|---|---|
| `scripts/data/library_wt_sequences.json` | `rows: 48840` | `65093` | **updated here** (it describes this file directly; the WT sequence itself is unchanged) |
| `scripts/data/cross_protein_corpus_stats.json`, `cr9114_h1` | `eligible_universe` 48841, `hit_count_ceiling` 489, `diversity_adjusted_ceiling_d2` 280, `_d3` 118, `enrichment_vs_base_rate` 0.62186, `mean_depth_of_hit` 8.91002, `max_depth_observed` 15, `tier_boundary_value` 9.60775502138942 | §6 | **not touched** — it is the freeze script's output over the `.xlsx`; needs deliberate regeneration, and that script reads the clean workbook, so it must be repointed at the raw table first |
| `docs/preregistration-cross-protein.md` | 48,841 at lines ~123, 342, 1114, 1171; 489 hits at ~1171 | §6 | **not touched** — pre-registered; the coordinator's call |
| `docs/dataset-inventory.md` | §3 table 48,840/48,841 and 74.52 % saturation; §3 floor caveat (wrong for H1, §2.4); §4a 489/280/118/54; §4c 0.6226 %/0.62x; §4c depth profile; §4d "max depth is 15"; §5 and §6 row counts | §2.4, §5, §6 | **not touched** |
| `docs/alde-compatibility.md` | 48,840 at lines 39, 64 | 65,093 | **not touched** |
| `docs/label-provenance.md` §3.5, §5 item 5 | "25 % unexplained", open item | resolved by this document | **not touched** |
| `scripts/analysis/cr9114_attrition.py` docstring | "our prepared file holds 48,840" | historical | leave as the record of the finding; its output is now the *negative* control |

The diversity-adjusted ceiling change (280 → 365 at d≥2) is the one with the most downstream
reach: it is the denominator the inventory uses to argue that raw hit counts are not comparable
across landscapes, and CR9114's share of a 125-measurement budget drops from 0.26 % to 0.192 %.

---

## 9. Sources

| # | item | identifier | read |
|---|---|---|---|
| 1 | Phillips A.M., Lawrence K.R., Moulana A., Dupic T., Chang J., Johnson M.S., Cvijovic I., Mora T., Walczak A.M., Desai M.M. (2021) *Binding affinity landscapes constrain the evolution of broadly neutralizing anti-influenza antibodies.* eLife 10:e71393 | [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393) | full text (per `label-provenance.md` §3) |
| 2 | Figure 1 — source data 1, CR9114 library −logK_D to H1, H3, influenza B | `elife-71393-fig1-data1-v3.csv`, sha256 `ee31f8fc…94fe` | **downloaded, parsed, 65,536 rows verified** |
| 3 | eLife asset manifest | `https://api.elifesciences.org/articles/71393` | downloaded |
| 4 | CombinGym, raw bnAb table | `Data/DMS/Raw/bnAbs_CR9114_raw.csv` @ `c5dccd51` | **downloaded; byte-identical to #2** |
| 5 | CombinGym, clean bnAb workbook | `Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx` @ `c5dccd51` | parsed; filter reverse-engineered exactly |
| 6 | Chen *et al.*, CombinGym | DOI `10.64898/2026.03.24.714074` | carried from `dataset-inventory.md` |
