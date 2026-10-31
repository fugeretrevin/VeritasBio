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

A fourth thing is **not a finding of this project at all**: the shallow-shell enrichment
**ordering** of §6.5. It is substantially a transcription of published figures into this
protocol's units (§6.5.1) and is a legitimate covariate, never a result.

### 1.1 Protocol deviation that must accompany every external claim

The protocol is **25 shared initial observations, then 4 rounds of 25, for 125 unique
measurements**. This is an **adapted** protocol. It is **not** ALDE's published schedule,
which is 96 initial plus 4 rounds of 96 = **480 measurements** averaged over 70
initializations — so this protocol runs at **26% of ALDE's demonstrated budget** at roughly a
quarter of its batch size. Any comparison against ALDE numbers — ours or the published ones —
must carry both sentences. A result quoted without them is a misquotation. §7.6 registers this
as a named threat to the baseline comparison and specifies the mitigation.

### 1.2 Established facts carried in (not to be re-derived)

Measured on held-out seeds on the development landscapes. Listed so the predictions below
are anchored, not so they are repeated as findings about anything else.

| fact | value |
|---|---|
| TrpB eligible universe | 159,999 variants; sites 183/184/227/228; parent `VFVS` |
| GB1 eligible universe | 149,360 variants (prepared CSV); sites 39/40/41/54; parent `VDGV` |
| HER2 trastuzumab CDR-H3 | 11,277 labelled rows, of which one is a literal `WT`, so the eligible universe of mutants is **11,276** and the assay-floor count **5,832** (not 11,277 and 5,833); 10 sites (99–108), against a 20^10 product — see §6.3 |
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
| CR9114 | binding_H1 | 16 | 65,536 | 48,841 **(clean workbook; the raw table holds 65,094 — see D-8)** |
| CR6261 | binding_H1, binding_H9 | 11 | 2,048 | 1,917 |
| GB1 | binding | 4 | 160,000 | 149,361 |

3 + 3 + 2 + 6×1 = **14 landscapes over 9 proteins**. Any source citing a different
collection size is describing a different version of CombinGym and **must not be mixed with
this pin**; reconcile it against the table above first.

Note that `Clean/` at this pin holds **CR9114 H1 only**. The CR9114 H3 and influenza B
landscapes from the upstream publication are not in the frozen corpus, so this protocol
cannot use them.

**There is no CR9114/H9 landscape, and there never was.** Phillips et al. assayed CR9114
against H1, H3 and influenza B; **H9 was measured for CR6261 only**. This is recorded
explicitly because the CR6261 H1/H9 pair invites the symmetric assumption, and a missing
landscape that was never measured is a different thing from one that exists upstream and is
absent from the pin. Nothing in this protocol may list a CR9114/H9 cell as unobtainable,
pending or absent: the experiment was not performed.

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
All **nine raw tables** under `Data/DMS/Raw/` are pinned as well, on the same terms:

| raw table | bytes | git blob SHA-1 | SHA-256 |
|---|---|---|---|
| `bnAbs_CR9114_raw.csv` | 12,258,680 | `8ac7352f539a9b4d33330dc3764b0ff5e86a094d` | `ee31f8fc26fce2730d2ae772e7a095aa334b4c6c5d8e7fbfe4b4a5a1c9a794fe` |
| `bnAbs_CR6261_raw.csv` | 358,244 | `1158a5dcfcc54a66467056ce096ecfbe13a06c36` | `6acc7dc8352d694db5f3aedb36e8c4216f61d400b4d85828c056f08580e9c61e` |
| `CreiLOV_raw.xlsx` | 22,619,510 | `300672c721f1fc1a5398f689c9c961dd6eaa9df0` | `02b314ba88bec25a6a814d04294ad415f4b150b6285b6d5a4ec0a848df5f8b70` |
| `eqFP611_raw.xlsx` | 737,089 | `0028edb5b41202512bfd26ca37bb14be25dd4875` | `d893d36645d10334b326696e98edbb7222b15bbeb946176482ed83637120904b` |
| `GB1_raw.xlsx` | 6,633,543 | `08b80e4367dcbe8ea8bfa70dff5d0908380c1581` | `ac9a16295cd24110ad3ef4bcc57741f023ff9f3f5c3d9f5fdc0a94b996dca9b7` |
| `HIV-1_raw.tsv` | 68,024 | `4611f98454408f6edfdacb570be57e3b34147ec6` | `0a174f172af42840463c94c80f29fe52c29926cab03524e9db9f1688f781f386` |
| `RhlA_raw.xlsx` | 287,545 | `a8704f81527fbd089bac63847fb734c0fa2c87a8` | `8302186312185d94062872198a144e493879875442d3b53bc6a33b347128ad0d` |
| `SaCas9_raw.xlsx` | 242,262 | `f7fafde24a570a11143452de5e808d933ea1cd35` | `ebdcc3edaf51821e6160db6760d1262145e66b43257dcebb46c18e1d123c7686` |
| `SpCas9_raw.xlsx` | 112,221 | `072d1c0064647fda885652c788dd82d899859776` | `2832e7f18f82b5b0a5769c0152505b2fdf15f8c0b720824d31c26603b89d80a3` |

Nine wild-type FASTA files are pinned on the same terms; their hashes are in the lock file
rather than inline here because nothing in the analysis depends on reading them by hand.

All **29 pinned files were fetched and verified**: every one matches both the git blob SHA-1
published by GitHub's tree API for the pinned commit and the SHA-256 recorded in the lock
file. `verify` exits 0 on the frozen corpus and exits 1 on a single altered byte.

**Why the raw tables are pinned, and why only one of them is read.** CR9114 is **scored**
against its raw table because the clean workbook is defective (D-8). The other eight are
pinned for a different reason: they are **the evidence that the clean tables are faithful**.
`docs/corpus-audit.md` establishes `clean − raw = 0` recoverable rows on every landscape but
CR9114, and that negative result is what licenses scoring the other landscapes from their
clean workbooks at all (§3.7). An unpinned raw table leaves that evidence unpinned — a later
audit could run against different bytes while `verify` still exited 0, which is precisely the
hole the CR9114 episode exposed. Pinning the evidence costs nothing and closes it.

**Both extensions were additive and the pin never changed.** The raw CR9114 table was pinned
first (21 files), then the remaining eight (29 files). On each occasion:

- the **commit, date and tree SHA of §3.1 are untouched** — every raw table's blob SHA-1 was
  cross-checked against the tree SHA `911f19c7…` recorded there;
- the lock's `corpus` block is **byte-identical** to the version committed at the freeze;
- **no existing entry changed and none was removed** — all were recomputed from the same local
  files and came out identical;
- the lock grew **20 → 21 → 29**, by addition only.

Stated this precisely so nobody later suspects the freeze was quietly broken: CR9114 being
scored against the raw table changes **which pinned file is read**, not what the corpus is.

Leaving the raw tables out was the alternative and it would have been worse. The lock as
originally frozen recorded only the 20 clean-and-metadata blobs and **no hash of any raw or
derived table**, so a landscape repointed at an unlisted file would have been scored against
bytes outside the freeze while `verify` still passed, because `verify` only walked the lock.
That hole is closed in two ways — the raw tables are pinned, and `verify` now fails if
`PINNED_FILES` or any landscape's source file is absent from the lock (§3.5).

One honest note on process: "no lock change" was the instruction on both occasions and was
wrong on both. The rule that survives is narrower and worth keeping — **the pin must not
change; the lock may only grow.** A frozen *pin* is the scientific commitment; a lock that
refuses to cover a file someone is about to read is not protecting anything.

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
disagrees with the script's pin, so a silent re-pin cannot pass unnoticed, and — added
2026-10-31 after the CR9114 episode showed the gap — it **fails if any pinned path, or any
file a landscape in the registry reads, is absent from the lock**. Before that guard, a
landscape repointed at an unlisted file would have been scored against unhashed bytes with
`verify` still exiting 0. The script
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

| D-8 | **`bnAbs_CR9114_H1_clean.xlsx` silently drops a whole quadrant of the library.** The exclusion filter is exactly `NOT(site 1 AND site 16) AND h1_mean is present`: every genotype carrying **both F29S and Y106S** is removed — **16,384 of them, one quadrant of the 2^16 space** — regardless of affinity, with zero false positives and zero false negatives, and a sweep over all 120 site pairs confirms (1,16) is the unique such pair. The loss is **depth-correlated in closed form**: a depth-k genotype carries two specified mutations with probability C(14,k−2)/C(16,k) = k(k−1)/240, so retention is 1 − k(k−1)/240, predicting **45.0% / 12.5% / 0%** at depths 12 / 15 / 16 against observed **45.0% / 12.5% / 0%**. The consequence worth stating in one sentence: **the fully matured antibody — the actual CR9114, with all 16 somatic mutations — was entirely absent**, because it necessarily carries both excluded mutations. The clean table omitted the antibody the landscape is about. | **Resolved, 2026-10-31.** CR9114 is scored against `Data/DMS/Raw/bnAbs_CR9114_raw.csv` at the same pinned commit (65,094 eligible rows against the workbook's 48,841). See §6.9 and `docs/cr9114-reprep.md`. The clean workbook is no longer read for CR9114. |

Where the summary and the clean file disagree, **the clean file wins**, because the clean
file is what is read — **except for CR9114**, where D-8 showed the clean file to be the
defective artifact and the raw table is read instead. That exception is the only one, and it
is a reminder that "the file we read wins" is a tie-breaking rule, not a guarantee that the
file we read is correct.

### 3.7 The clean-versus-raw audit: D-8 is unique, and no covariate moves

O-14 asked whether a second silent exclusion of the D-8 kind lurked in the other nine clean
workbooks, four of them final-test targets whose covariates this protocol has already
pre-registered. It does not. Full detail in `docs/corpus-audit.md`.

> **The D-8 defect is unique to CR9114. Not one pre-registered covariate moves. The freeze
> holds. No landscape needs re-preparing.**

**The reasoning, which is what makes the negative meaningful.** A verdict of "we looked and
found nothing" is worth little; what makes this one load-bearing is the quantity it turns on.

The decisive quantity is **recoverable rows** — rows missing from a clean table yet carrying a
**real measurement** in the raw table. That count is **16,253 for CR9114 and exactly 0 for
every other landscape.** The argument is then closed rather than statistical: a row with no
measurement **can never enter the eligible universe** under the curation rule of §3.3, so its
absence cannot shift a covariate **however depth-correlated the absence looks**. Recomputing
every covariate over the raw universe gives a delta of **exactly 0.0000** except for CR9114.

So the clean workbooks are faithful **on row membership, value fidelity and depth structure**
for all nine, and the eight non-CR9114 landscapes may be scored from them.

**Why the audit's negatives are credible: the same code demonstrably finds the positive.** A
null from an instrument never shown to detect anything is not evidence. This one:

- reproduced all **seven frozen §6.5 covariates digit for digit**;
- **re-derived CR9114's filter independently** at 0 false positives and 0 false negatives;
- reproduced the documented **reweighting bound at 0.5622x** (§6.9.1).

A composite sweep over every single-site and site-pair rule of the CR9114 form found nothing
at 0 FP / 0 FN anywhere else: CR6261's best candidate misses by **423 kept rows** and
SpCas9's by **361**. And **eqFP611 and SaCas9 are combinatorially complete** — every genotype
of their design space is present — so there is no room for a filter to hide in either.

#### 3.7.1 Scope caveat, which bounds the clean bill of health

The audit examined **row membership, value fidelity and depth structure**. It did **not**
examine each clean table's **column choices** against its upstream publication.

That distinction is not hypothetical: the SaCas9 phenotype divergence recorded in §6.10
surfaced **incidentally** rather than by design. So more of that kind may exist, particularly
in the reserve landscapes nobody has had reason to read closely.

> **A clean bill on row membership is not a clean bill on phenotype definition.**

Phenotype-column provenance is therefore a **precondition on promotion**, alongside the
anticipation audit (§6.5, O-15). For the four final-test targets the columns in use are the
ones the registry names and the §9.4 censoring figures describe, which is a weaker check than
a publication comparison and is recorded as such in O-16.

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
**because HER2 needed it**: a 10-site library with 11,277 labelled rows — 11,276 mutants plus
a literal `WT`, §1.2 — against a 20^10 product makes free generation propose variants nobody
measured, so the batch came up short.
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

### 4.5 Suitability finding: CR6261 is excluded from the pooled aggregate

A target can be legitimately held out and still be the wrong instrument. CR6261 is
**unsuitable for pooled performance comparison** on three independent grounds, all computable
before any arm runs and all recorded here rather than discovered in the results:

1. **Budget-to-library ratio.** 125 measurements is **6.63% of its 1,886 mutants**, against
   **0.08%** on TrpB — a factor of 83. A difference measured while a campaign has consumed
   one variant in fifteen is not the same quantity as one measured at one in twelve hundred.
2. **The budget exceeds the hit set.** Its entire top-1% tier is **19 variants**, so a
   125-measurement campaign could in principle measure every hit several times over. At that
   ratio **an oracle scores identically to a merely good method**, and the endpoint stops
   discriminating between them. This is why the tier ladder of §9.1 deepens CR6261 to 2%, and
   deepening mitigates but does not remove the problem: at 2% the tier is 38 and the d ≥ 2
   ceiling is 18 (§9.3), still below the budget.
3. **Its H9 covariate is degenerate at zero** under both cut definitions — 0.000x at
   floor(S/2) and 0.000x at depth ≤ 2 (§6.3) — so it contributes a tied value to the H1
   ordering rather than a measurement.

**Disposition, fixed in advance.** CR6261 H1 and H9 are **excluded from the
equal-protein-weighted pooled aggregate** of §12.2 item 5, and are analysed instead as a
**separate fixed-fraction stratum** (§8.2). They are **retained in full** for everything else,
and specifically retained as:

- the **qualitative demonstration that the covariate can invert** — H9 has zero top-1% hits
  anywhere in its shallow shell while its hit rate rises to 22.9% at depth 9 and 40.0% at
  depth 10, which is the most extreme depletion in the corpus and is genuinely valuable as a
  structural contrast;
- the **replication instance for H2** (§6.7.5), where its weaknesses are already stated.

Excluding it from the pooled aggregate reduces that aggregate to **3 proteins** — eqFP611,
CreiLOV, CR9114 — and that figure, not 4, is what any pooled number reports. The H1 ordering
still uses all four proteins, because an ordering tolerates a tied endpoint where a
pooled mean does not; §6.4's power statement already accounts for the tie.

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
| `cr9114_h1` | CR9114 | immunoglobulin heavy-chain variable domain | binding, H1 | 16 | **65,094** (65,093 mutants; raw table, D-8) |
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

### 5.3 The whole study: seven proteins, five structural families

The count above is about **held-out evidence**. The count across the whole study, development
and final test together, is different and is recorded separately so neither can be quoted for
the other:

| structural family | proteins | level |
|---|---|---|
| tryptophan synthase β subunit | TrpB | development |
| protein G domain B1 | GB1 | development |
| **immunoglobulin V-domain** | **trastuzumab, CR9114, CR6261** — three proteins, **one family** | trastuzumab development (§4.2); CR9114 and CR6261 final test |
| GFP-like fluorescent protein | eqFP611 / mTagBFP2 | final test |
| LOV / flavin-binding photosensor | CreiLOV | final test |

**Seven proteins, five families.** The CombinGym acquisition contributes **+2 families** the
project previously had no data in at all — GFP-like fluorescent protein and LOV photosensor —
taking the study from three families to five.

**And one caveat that the three-family figure in §5.2 hides.** Of those three independent
families, the immunoglobulin V-domain is **already represented in development data** by
trastuzumab (§4.2). CR9114 and CR6261 are different proteins from trastuzumab, so they are
held out at the **protein** level; they are **not** held out at the **family** level. The
practical reading:

- **2 families** are genuinely new to the project — eqFP611 and CreiLOV.
- **1 family** — immunoglobulin V-domain — carries held-out proteins but has a development
  sibling whose inspection drove an algorithm change (the pool-restricted path, §4.2).

Any claim of transfer to a *new* protein family therefore rests on **two** families, not
three. The antibody results remain held-out evidence at the protein level and are reported as
such, with this caveat attached.

---

## 6. Pre-registered hypotheses

Two, testing different things on different designs:

- **H1 (§6.1–6.6)** — *generalization.* Does the advantage track shallow-shell richness
  **across protein families**? Unit: protein, n = 4. Power marginal by construction (§6.4).
- **H2 (§6.7)** — *mechanism.* Is shallow-shell richness what drives the advantage **within
  one protein**? Unit: seed, n = 48. Adequately powered.

§6.8 fixes the order in which they are interpreted and states that neither substitutes for
the other.

### 6.1 H1 — the cross-protein hypothesis

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
| deep-shell (depth ≥ 3) hit rate | 1,432/157,757 = **0.9077%** | 1,422/147,193 = **0.9661%** |
| **shallow-shell (1 ≤ depth ≤ 2) hit rate** | **168/2,242 = 7.4933%** (7.5x base rate) | **72/2,167 = 3.3226%** (3.3x) |
| hit rate by depth, TrpB | d1 15.79%, d2 7.20%, d3 2.33%, d4 0.608% | — |
| best variant at depth 4 | 1.0000 | 8.7620 |

Both rows of shell figures have now been computed by **three independent implementations** —
the original mechanism analysis, a parallel agent's own genotype canonicalisation and shell
profiler, and `freeze_cross_protein_corpus.py stats` against the frozen corpus — and they
agree to the digit (§16 O-4).
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

**Depth 0 — the wild type — is excluded from both shells and retained in the eligible
universe.** The shells partition the *mutants*: shallow is `1 ≤ depth ≤ cut`, deep is
`depth > cut`, and the wild-type row belongs to neither while still counting toward the tier
and the base rate. A wild type is not a search target.

This is stated because it had already produced a one-row divergence between implementations.
Several datasets carry a literal `WT` row and several do not, so the same definition gives
different denominators depending on whether depth 0 is swept in:

| | with depth 0 | excluding depth 0 (**adopted**) |
|---|---|---|
| GB1 shallow | 72/2,168 = 3.3210% | **72/2,167 = 3.3226%** |
| GB1 deep | — | **1,422/147,193 = 0.9661%** |
| TrpB shallow | 168/2,242 (no depth-0 row exists) | **168/2,242 = 7.4933%** — convention is a no-op |

GB1's universe of 149,361 rows splits as 1 + 2,167 + 147,193; its CSV carries a literal `WT`
row and TrpB's 159,999 does not. The same cause was confirmed independently on HER2, whose
`train` split carries an explicit `WT` row: its eligible universe is **11,276**, not 11,277,
and its assay-floor count is **5,832**, not 5,833. The adopted convention reproduces the
development anchors exactly and is what `freeze_cross_protein_corpus.py` implements. The
difference is immaterial in size — and two anchors computed under different conventions is
precisely the drift a freeze exists to prevent, which is why it is written down rather than
left to whichever implementation runs last.

**The shell cut cannot be a fixed absolute depth, and the reason is resolution, not taste.**
The candidate libraries are mostly *binary* — one alternative residue per site, not 19 — so
the depth ≤ 2 shell holds between 58 and 196 variants on them, against 2,242 on TrpB. At a
base rate near 1% such a shell yields an expected hit count of one or two, and the enrichment
estimate is unusable. Measured on the frozen corpus:

| landscape | sites | depth ≤ 2 shell size | hits in it | depth ≤ 2 enrichment (95% CI) | resolvable? |
|---|---|---|---|---|---|
| TrpB (dev) | 4 (20-way) | 2,242 | 168 | 7.49x | yes |
| GB1 (dev) | 4 (20-way) | 2,167 | 72 | 3.32x (2.64–4.16) | yes |
| eqFP611 blue | 13 (binary) | 91 | 13 | 14.27x (8.53–22.90) | no |
| eqFP611 combined | 13 (binary) | 91 | 13 | 14.27x (8.53–22.90) | no |
| eqFP611 red | 13 (binary) | 91 | 0 | 0.00x (0.00–4.05) | no |
| CreiLOV | 15 (restricted) | 196 | 18 | 9.18x (5.89–14.04) | no |
| CR9114 H1 | 16 (binary) | 134 | 0 | 0.00x (0.00–2.78) | no |
| CR6261 H1 | 11 (binary) | 58 | 0 | 0.00x (0.00–3.08) | no |
| CR6261 H9 | 11 (binary) | 65 | 0 | 0.00x (0.00–2.78) | no |

Four of the seven independent landscapes tie at exactly zero with intervals spanning the
whole plausible range. A rank test on that column is a test on ties.

**Pre-specified primary shell definition: 1 ≤ depth ≤ floor(S/2), where S is the number of
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
universe and returns an enrichment of 0.98–1.00 — it leaves the deep shell empty, which is
the degeneracy the guard is for. **All three RhlA landscapes are therefore excluded from the
H1 test**, though they remain usable for every other endpoint. HIV-1 protease is a borderline
case — its floor(11/2) = 5 shell holds 72.4% of the universe and 58 of its 59 hits, leaving a
deep shell with zero — and is flagged rather than excluded; its enrichment of 1.357x at the
primary cut against 8.71x at depth ≤ 2 is the most severe cut-sensitivity in the corpus and
must be reported with both values.

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

**H1's power is marginal by construction, and here is the demonstration.** The test needs a
perfect ordering of four items, and the covariate's own top three are separated by margins
that move on bookkeeping: adopting the depth-0 convention of §6.3 — excluding **one**
wild-type row from a denominator — swapped landscape ranks 1 and 2, CreiLOV 1.875x against
eqFP611 blue 1.854x, where before the fix it had been 1.880x against 1.876x in the other
order. A 0.02x margin moved on a single row. The per-protein ordering survived that change,
which is why it and not the landscape ordering is the primary unit, but the episode sets a
realistic expectation: **H1 will most likely be recorded as unresolved.** It is kept because
an unresolved pre-registered test is still worth more than a post-hoc explanation, not
because it is expected to decide anything. The mechanism itself is tested by H2 (§6.7), which
is where the power actually is.

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

| rank | landscape | level | shell cut | shell size | hits | shell hit rate | enrichment (95% CI) | published anticipation |
|---|---|---|---|---|---|---|---|---|
| — | TrpB | development | 1–2 | 2,242 | 168 | 7.4933% | 7.49x | **none** — Johnston 2024 publishes no fitness-versus-mutation-count object |
| — | GB1 | development | 1–2 | 2,167 | 72 | 3.3226% | 3.32x (2.64–4.16) | **none** — Wu 2016 never stratifies by mutation count |
| 1 | `creilov` | final test | 1–7 | 62,790 | 1,178 | 1.8761% | **1.875x** (1.77–1.98) | **anticipated** — Chen Fig 3c |
| 2 | `eqfp611_blue` | final test | 1–6 | 4,095 | 76 | 1.8559% | **1.854x** (1.48–2.31) | **anticipated** — Poelwijk Fig 6a |
| 3 | `eqfp611_combined` | final test | 1–6 | 4,095 | 75 | 1.8315% | **1.830x** (1.46–2.29) | **anticipated** — Poelwijk Fig 6a |
| 4 | `eqfp611_red` | final test | 1–6 | 4,095 | 26 | 0.6349% | **0.634x** (0.43–0.93) | **anticipated** — Poelwijk Fig 6a |
| 5 | `cr9114_h1` | final test | 1–8 | 38,835 | 227 | 0.5845% | **0.584x** (0.51–0.67) | **anticipated** — Phillips Fig 1D, scoped to H1 (§6.9.2) |
| 6 | `cr6261_h1` | final test | 1–5 | 945 | 4 | 0.4233% | **0.210x** (0.08–0.54) | **partly anticipated** — Phillips Fig 1E/1F |
| 7 | `cr6261_h9` | final test | 1–5 | 917 | 0 | 0.0000% | **0.000x** (0.00–0.21) | **partly anticipated** — Phillips Fig 1E/1F |

**Ranks 1–3 are a near-tie and ranks 4–5 are a near-tie,** at 1.875/1.854/1.830 and
0.634/0.584, with heavily overlapping intervals. The ordering within each cluster is not
meaningfully resolved and the rank test should be read as separating four groups — {CreiLOV,
blue, combined} > {red, CR9114} > {CR6261 H1} > {CR6261 H9} — rather than seven distinct
levels. Concretely: adopting the depth-0 convention of §6.3 swapped ranks 1 and 2 (CreiLOV
1.8753 against blue 1.8541, previously 1.8761 against 1.8799), which is how little separates
them. This is stated now so that a `rho_s` driven by an unresolved within-cluster ordering is
not reported as a clean monotone result.

The per-protein rows the primary test uses, from the means of the above:

| protein | landscapes | mean primary enrichment | position anticipated by | status |
|---|---|---|---|---|
| CreiLOV | 1 | 1.875x | **Chen Fig 3c** — "fraction of sampled genotypes observed to be >WT" against number of mutated residues, falling monotonically from ~0.52 at depth 2 to ~0 by depth 13 | anticipated |
| eqFP611 / mTagBFP2 | 3 | 1.439x (spread 0.634–1.854) | **Poelwijk Fig 6a** — above-threshold genotypes coloured blue for steps 0–~6 and red for ~7–13 | anticipated |
| CR9114 | 1 | 0.584x | **Phillips Fig 1D, scoped to H1** — per-shell affinity boxplots whose interquartile range enters the ceiling band only from ~depth 9–10 (§6.9.2) | anticipated |
| CR6261 | 2 | 0.105x (spread 0.000–0.210) | **Phillips Fig 1E/1F** — the CR6261 germline binds *neither* antigen, both pinned at the ≤7.0 floor, while Fig 1C/1D shows the CR9114 germline already binding H1 at −logK_D ≈ 8.4 | partly anticipated, leaning anticipated |

Citations in full:

- **Poelwijk FJ, Socolich M, Ranganathan R (2019).** *Learning the pattern of epistasis linking
  genotype and phenotype in a protein.* Nat Commun 10:4213. DOI
  [10.1038/s41467-019-12130-8](https://doi.org/10.1038/s41467-019-12130-8).
- **Chen Y, Hu R, Li K, Zhang Y, Fu L, Zhang J, Si T (2023).** *Deep mutational scanning of an
  oxygen-independent fluorescent protein CreiLOV…* ACS Synth Biol 12(5):1461–1473. DOI
  [10.1021/acssynbio.2c00662](https://doi.org/10.1021/acssynbio.2c00662).
- **Phillips AM, Lawrence KR, Moulana A, et al. (2021).** *Binding affinity landscapes constrain
  the evolution of broadly neutralizing anti-influenza antibodies.* eLife 10:e71393. DOI
  [10.7554/eLife.71393](https://doi.org/10.7554/eLife.71393).
- Clean negatives: **Johnston et al. (2024)** PNAS 121:e2400439121, DOI
  [10.1073/pnas.2400439121](https://doi.org/10.1073/pnas.2400439121) — no fitness-versus-mutation-count
  object of any kind; **Wu NC et al. (2016)** eLife 5:e16965, DOI
  [10.7554/eLife.16965](https://doi.org/10.7554/eLife.16965) — a global fitness distribution,
  never stratified by mutation count.

Four proteins, a strict ordering, no ties — and, unlike the landscape-level ordering, one
that is **unchanged** by the depth-0 convention.

#### 6.5.1 The ordering is pre-published. The relationship is not.

**At most one of H1's four positions is genuinely blind, and arguably zero.** Three published
figures fix all four positions in direction and three of four in rank, once this protocol's
floor(S/2) cut is supplied. The full audit is `docs/ordering-anticipation.md`; the verdicts
are in the tables above. Note what this does *not* rescue: it is no use pleading that our cut
rule was needed to produce the exact numbers, because the cut rule is in this document and a
reader has it.

The distinction that saves H1, in exactly these terms:

> **The covariate column is published. The outcome column is not.**

H1 asserts a **relationship** between shallow-shell enrichment and this project's
**per-landscape effect sizes**. The ordering was knowable from the literature; the effect
sizes are not, and §6.5's prediction of across-the-board loss to ridge is **entirely blind**.
So H1 is not void. The compact formulation, to be used verbatim rather than paraphrased:

> **H1's derivation was original and unpublished; its test set is pre-published.**

**What must be struck.** Any sentence presenting **the ordering itself** as a finding of this
project. The enrichment ordering is a transcription of published figures into this protocol's
units. It is a legitimate covariate and an illegitimate result. Every table reporting the
ordering carries the anticipation column above.

**What survives untouched.** The per-landscape effect sizes, the §12.2 reporting set, the
material-regression and minimum-improvement findings, H2 (§6.7), and the prediction below.

**Two of this corpus's five characterised proteins are clean negatives — and they are the
development pair.** TrpB and GB1, where the mechanism was derived, have no published
depth-stratified object at all. So the derivation genuinely could not have been read off the
literature; only the test set could. That is the right way round for honesty about the
derivation and the wrong way round for the strength of the test.

**Required alongside this ordering, every time it is reported: each target's library site
count** — 16 for CR9114, 15 for CreiLOV, 13 for eqFP611, 11 for CR6261. This is disclosure,
not a correction: the Spearman correlation between site count and enrichment across the four
proteins is only **rho_s = +0.400 (exact one-sided p = 9/24 = 0.3750)**, because CR9114 is the
widest library in the corpus yet third of four in enrichment, and that single inversion
decorrelates the two. **Width is not a serious confound for H1.** The reason width still
matters is different and more serious — under a proportional depth bias the policy would not
be a shallow exploiter on any of these libraries at all, which would make the mechanism
inapplicable to every independent target. See §16 O-12, resolved from the §10.1 read-out.

**The prediction this forces, stated plainly.** Every independent enrichment measured
(≤ 1.88x) sits **below GB1's 3.32x**, where reliability already lost to ridge by 22.49 raw
top-1% hits. H1 therefore predicts that **the current reliability configuration will lose to
ridge on every independent landscape in this corpus**, most severely on `cr6261_h9` and
`cr6261_h1`, and least severely on `creilov`, `eqfp611_blue` and `eqfp611_combined`. If
reliability instead wins broadly, H1 is wrong and so is the mechanism story behind it — and
that is a more interesting result than the one predicted.

Registering a prediction of across-the-board failure is the point of doing this before the
runs. It is also the result that would be hardest to publish and easiest to re-narrate
afterwards, which is why it is written down first.

**This prediction is not weakened by §6.5.1 and is not to be softened in response to it.** The
anticipation audit struck the *ordering* as a finding; it leaves this prediction entirely
blind, because no published figure says anything about how the VeritasBio reliability policy
performs against one-hot ridge on any of these landscapes. With the ordering demoted to a
covariate, **this prediction is now the main thing H1 tests**, and it is the part of H1 that
was always the real content.

Two features of this corpus are worth naming before the numbers arrive, because they make
the test sharper than n = 3 families usually allows:

- **`eqfp611_blue` against `eqfp611_red` holds the protein, the library and the 8,192
  genotypes fixed and varies only the shell structure** (blue: hit rate falls monotonically
  with depth, 30.8% at depth 1, 0.08% at depth 8; red: *zero* hits at depth 0–3, peaking at
  1.63% at depth 8). It is the cleanest available within-protein manipulation of the
  covariate, and it is not independent cross-protein evidence. Both statements are true at
  once — which is why it is now a pre-registered hypothesis of its own, **H2 (§6.7)**, tested
  separately from H1 and licensing a mechanism claim rather than a generalization claim.
- **CR9114 and CR6261 are affinity-maturation landscapes where hit rate rises monotonically
  with depth** (CR6261 H9: 0% up to depth 5, 22.9% at depth 9, 40.0% at depth 10). They are
  structurally the opposite of TrpB and are the hardest case for a shallow exploiter.

One caveat on the covariate itself, recorded because it cuts against the clean reading above:
**CreiLOV's depth profile is shallow-favouring but its primary shell cut hides that.** Its
hit rate falls from 9.66% at depth 2 to 0.12% at depth 12, yet floor(15/2) = 7 admits 62,790
variants — 38% of the universe, spanning the whole productive range — so the enrichment comes
out at 1.875x while the depth ≤ 2 cut gives 9.18x and the depth ≤ 4 cut 4.6x. The primary cut
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
exact one-sided p of 1/5,040 instead of 1/24.

#### Standing prohibition: no n-expansion without an anticipation audit first

That apparent gain is a trap, and it is recorded here as a precondition rather than left as a
temptation:

> **No reserve protein may be promoted into the H1 test until the §6.5.1 anticipation audit
> has been run on it and its verdict recorded in this document.**

The reasoning. Promoting four more **pre-published** positions would move the exact attainable
one-sided p from **1/24 to 1/5,040** while adding **not a single blind prediction**. It would
read as a decisive power gain and be nothing of the kind: a longer transcription of published
figures is still a transcription. The audit is cheap — it is a literature read, not a compute
job — and it must precede promotion, because after promotion the number looks authoritative
and the question stops being asked.

Each promoted protein's row carries its anticipation verdict in the §6.5 tables, and the
promoted test reports **how many of its n positions are blind** alongside the p-value. A
p-value from an ordering with zero blind positions is reported as what it is: a measure of how
well published figures agree with each other.

**A promoted protein needs two clearances, not one:** the anticipation audit above, **and** a
provenance check of the kind §3.7.1 shows is still outstanding. The two reserve findings in
§6.10 are why.

Their covariates are already computed and frozen alongside the others, so promoting them
later is not a fresh look: HIV-1 protease 1.357x, SaCas9 1.560x, SpCas9 2.500x, RhlA
undefined. Note that SpCas9's 2.50x is the highest independent value in the whole corpus and
sits between GB1's 3.32x and the rest — which is precisely why its 20.9% universe consumption
must not be allowed to carry an H1 result.

### 6.6 H1's training data

H1 was derived from TrpB and GB1. Those two landscapes are its training data and **cannot be
counted as evidence for it**. They appear in §6.5 as anchors so the independent values can be
read against a known scale, not as two of the test's observations. The test is run on the
independent families only.

### 6.7 H2 — the within-protein mechanism test

> **H2.** On a single protein and a single library, the reliability policy's advantage over
> ridge is larger on the readout whose shallow shell is hit-enriched than on the readout whose
> shallow shell is depleted.

H2 tests the **mechanism**. H1 tests **generalization**. They are different hypotheses with
different designs and different power, and neither substitutes for the other (§6.8).

#### 6.7.1 Design — why this is the cleaner causal comparison

The primary instance is **eqFP611 blue against eqFP611 red**.

| held constant | varies |
|---|---|
| protein (eqFP611 / mTagBFP2) | shallow-shell enrichment: **1.854x (blue) against 0.634x (red)**, a 2.9x ratio |
| the genotype set — the identical 8,192 sequences, row for row | the depth at which hits sit: blue's mean hit depth 3.96, red's 7.11 |
| library design, 13 sites, binary alphabet, 100% saturation | |
| assay platform and publication (Poelwijk et al. 2019) | |
| tier (1%), tier size (82 hits), base rate (1.0010%) | |
| diversity ceiling geometry (d≥2: 38 blue, 28 red) | |
| censoring (both effectively uncensored, 0.012% at each extreme) | |

**This is a natural experiment Poelwijk built, not a manipulation this protocol selected —
and that is stronger provenance, not weaker.** The library is not "13 mutated sites" in the
abstract. It is the **complete 2^13 space between two named parents**: **mTagBFP2** (blue, at
depth 0) and **mKate2** (red, at depth 13). The pinned wild-type file is literally
`eqFP611(mTagBFP2).fasta`, so depth *is* the interpolation coordinate between the two
phenotypes, by the original authors' design. Blue's hits sit near its own parent and red's
near the other one because that is what the library is.

Poelwijk Fig 6a makes it explicit: every above-threshold genotype is coloured by mutational
step, **blue from 0 to ~6 and red from ~7 to 13**. This protocol's cut is floor(13/2) = **6**,
so the **published colour switch lands exactly on it.** The cut was derived from the
resolvability argument of §6.3, not from that figure, and it coincides.

Why this is better provenance than a contrast of our own construction: the blue/red split
**cannot have been chosen after seeing our results**, because it was fixed by a 2019
publication designing an interpolation between two fluorescent proteins. A within-protein
covariate manipulation that predates the hypothesis by seven years is about as clean a causal
handle as an observational corpus can offer. It is still not cross-protein generalization
evidence (§6.7.5), and it is now also a position the literature anticipates (§6.5.1) — but
H2's content was never the ordering; it is the **effect-size difference**, which no figure
publishes.

The cross-protein comparison of H1 varies protein, fold, library width, alphabet, saturation,
tier, censoring and ceiling geometry **all at once**, and then asks which of them explains a
difference. This comparison varies the covariate and holds the rest fixed. It is a better
causal design by a wide margin, and it was already in this document as a remark (§6.5) rather
than a test. That was an error of nerve, corrected here.

What it is **not**: evidence of cross-protein generalization. One protein cannot speak to
transfer across families. §6.7.5 states that limit without hedging.

#### 6.7.2 Power

Paired by seed over **48 seeds**, against n = 4 proteins for H1. That is the whole point.

#### 6.7.3 The statistic

A difference of differences, with the seed as the unit of independence throughout:

```
for each seed s:
    D_s = [ Y_s(reliability, blue) − Y_s(ridge, blue) ]
        − [ Y_s(reliability, red)  − Y_s(ridge, red)  ]
```

on the primary endpoint (diversity-adjusted yield at d ≥ 2) at the 125-measurement budget.
Inference on `{D_s}` is a paired t-test over seeds, df = 47, plus a **cluster bootstrap over
whole seeds** (B = 10,000, percentile interval) and a Wilcoxon signed-rank, exactly as §11
requires of everything else. No within-seed sub-unit is pooled. The same contrast is reported
at all five budgets and on the raw secondary endpoint; the 125-measurement primary-endpoint
cell is the pre-specified test.

Seeds are shared across the two readouts: seed s draws **the same 25 initial genotypes** on
blue and on red, which is possible because the genotype set is identical. That makes `D_s` a
within-seed, within-genotype-set contrast and removes initial-sample variation from it
entirely. Alpha 0.05, one-sided, predicted direction positive.

#### 6.7.4 Two predictions of unequal strength

They are separated here because they do not rest on the same premises, and conflating them
would let one fail and be recorded as though the other had.

**P1 — the directional prediction. This is H2's primary claim and it is robust.**

> `D_s` is **positive**: reliability's edge over ridge is larger on blue than on red.

P1 holds for **any** selection depth that sits below red's hit mass. Red has **zero top-1%
hits at depths 0, 1, 2 and 3** and all 82 of its hits at depth 4 and deeper; blue's hit rate
falls monotonically from depth 1. So any policy biased shallower than red's hit mass does
relatively better on blue, whatever that bias turns out to be numerically. P1 does not depend
on resolving the question in P2 and is scored independently of it.

Blue's shell is also the richest at depth ≤ 2 in the entire independent corpus — 14.27x, 13
hits in 91 variants — so if the mechanism holds, blue is where reliability has its best chance
on any independent landscape.

**P2 — the point prediction. Contingent, and labelled so.**

> **Conditional on the policy's depth bias being absolute rather than proportional to library
> width:** on eqFP611 red the reliability arm finds **essentially nothing** — a
> diversity-adjusted yield at or near zero at every budget, because the region it samples
> contains no qualifying variants at all.

The premise that cannot be assumed: *"reliability selects at mean depth ~2.6"* was measured on
**four-site** libraries (TrpB 2.63, GB1 2.56). eqFP611 has **thirteen** binary sites. Nothing
measured so far establishes which of two readings transfers:

| reading | selection depth on a 13-site library | where red's hits are | P2 |
|---|---|---|---|
| **absolute** — the bias is ~2.6 mutations | ~2.6 (20% of sites) | zero hits at depth ≤ 3 | holds: near-zero yield |
| **proportional** — the bias is ~65% of sites (2.6/4) | ~8.5 (65% of sites) | **peak hit rate, 1.63% at depth 8** | fails: red may do well |

The two readings make **opposite** predictions about red, and the proportional reading puts the
policy exactly on red's hit mass. The question is genuinely open: our only evidence is two
libraries of **identical width**, which cannot distinguish an absolute constant from a
proportion. No prior is asserted here beyond that.

**Interpretation rule, fixed in advance.** A healthy diversity-adjusted yield on red
**falsifies P2**. It does **not** by itself falsify H2 if the measured mean selected depth on
eqFP611 turns out to be proportional rather than absolute — that outcome is a statement about
how the depth bias scales with library width, not about whether shallow-shell richness drives
the advantage. In that case P2 is recorded as falsified, the width-scaling finding is reported
as such (§16 O-12), and H2 is scored on P1 alone. The reverse is also fixed: if the measured
selection depth is absolute (~2.6 on a 13-site library) and red nonetheless yields well, that
**is** a refutation of the mechanism, because the policy sampled a region containing no hits
and found hits anyway, which cannot happen.

**P2 is resolved from measurement, not from assumption.** §10 requires mean selected depth to
be recorded per round per target, in both absolute and fractional terms, and §10.1 requires it
to be read out **before** P2 is scored.

Note the honest tension with H1: blue's floor(S/2) enrichment of 1.854x is still well below
GB1's 3.32x, so H1 predicts reliability loses to ridge on blue in absolute terms. H2's P1
predicts only that it loses *less* on blue than on red. **Both can hold at once**, and if they
do, the reportable statement is "the mechanism is confirmed and the policy still loses
everywhere", which is a coherent and uncomfortable result, not a contradiction.

#### 6.7.5 Replication instance, and falsification

**Replication: CR6261 H1 against CR6261 H9** (0.210x against 0.000x). Same structure — one
protein, one 11-site library, two antigens, two readouts of the same genotypes. Pre-registered
as a replication of H2 with the same statistic.

It is the **weaker instance** and that is stated in advance, not discovered later: both
enrichments sit near zero, so the contrast is small (0.21x against 0.00x, against blue-red's
1.854x against 0.634x); the d ≥ 2 ceilings are only 18 and 21, so the endpoint saturates
(§9.3); and the two antigens' genotype sets are **not** identical (1,887 rows against 1,842),
so the shared-initial-sample construction of §6.7.3 applies only to the genotypes present in
both. The size of that intersection has not been computed and is not asserted here; it must be
computed and recorded before the contrast is run, and the contrast is computed on the
intersection with its size stated in the table. A null on this instance is weak evidence
against H2; a reversal is not.

**H2 is disconfirmed if any of the following holds.** All four are conditions on P1, which is
H2's claim; P2's failure is scored separately per §6.7.4 and is not on this list.

1. no blue-minus-red difference — the 95% seed-cluster CI on `mean(D_s)` contains zero at the
   125-measurement budget;
2. a difference in the **wrong direction** — reliability's edge is larger on red;
3. reliability returns a materially non-zero diversity-adjusted yield on **red** *and* the
   measured mean selected depth on eqFP611 is **absolute** (near 2.6 of 13 sites, i.e. below
   red's hit mass), while a blue-red difference nonetheless appears. That combination is
   incoherent under the mechanism — the policy would have sampled a region containing no
   qualifying variants and found them anyway — so it means the difference is real but not
   caused by shallow-shell depletion. **If instead the measured selection depth is
   proportional (near 8.5 of 13), a healthy red yield is not a disconfirmation of H2 at all**;
   it falsifies P2 only, and the width-scaling result is reported in its own right (§16 O-12).
   An earlier draft of this clause said a healthy red yield "kills the mechanism regardless of
   what blue does" — that was too strong, and it would have let a width-transfer effect be
   misrecorded as a refutation of the mechanism;
4. the CR6261 replication reverses (not merely nulls).

**The limit, stated plainly.** A positive H2 establishes that **shallow-shell richness drives
the advantage within one protein**. It says **nothing on its own** about transfer across
protein families. It cannot be quoted as cross-protein evidence, it cannot be aggregated with
H1, and a sentence that cites H2 for generalization is a misreport. H2's own training data
problem is also recorded: the covariate was derived from TrpB and GB1 (§6.6), so eqFP611 is
legitimately held out for H2, but the blue-red contrast was *noticed* while characterising the
corpus for this document — the covariate values in §6.5 were computed from labels before any
arm ran, which is what makes H2 prospective, and no arm performance on eqFP611 had been
examined when H2 was written.

### 6.8 Interpretation order: H1 and H2 answer different questions

Fixed in advance so the result cannot be re-framed afterwards.

| | H1 (§6.1–6.6) | H2 (§6.7) |
|---|---|---|
| question | does the advantage transfer across protein families? | is shallow-shell richness what drives the advantage? |
| design | 4 proteins, everything varies | 1 protein, 1 genotype set, only the covariate varies |
| unit, n | protein, n = 4 | seed, n = 48 |
| power | marginal by construction; min attainable p = 0.042 | adequate |
| what it licenses | a generalization claim | a mechanism claim |

**Neither substitutes for the other.** H2 cannot rescue a failed H1 and H1 cannot rescue a
failed H2. Reported in this order: H2 first, because it is the better-powered test of the
narrower question; then H1, with its power limitation restated; then the per-protein effects
of §12.2, which stand on their own regardless of either hypothesis.

**The reporting order matters more than it did when it was fixed.** §6.5.1 established that
H1's covariate ordering is substantially a transcription of published figures into this
protocol's units — at most one of four positions blind, arguably zero. H1 retains real content
(the relationship, and the blind across-the-board-loss prediction of §6.5), but **H2 now
carries more of the study's weight** than when this order was written. H2's own content was
never an ordering: it is the per-seed effect-size difference between two readouts, which no
publication reports. H2 first was the right call for power; it is now also the right call for
originality.

**The most likely outcome is pre-named: H2 confirmed, H1 unresolved.** That combination is
**mechanism supported, generalization undetermined**, and it is a legitimate and reportable
result. It is not a failure and must not be narrated as one; nor may it be narrated as a
generalization success. The sentence it licenses is: *within one protein and one library, the
policy's advantage tracks shallow-shell richness as predicted; whether that transfers across
protein families is undetermined at n = 4 proteins, and the per-protein effects are as
tabulated.*

The symmetric warning: **H1 confirmed with H2 falsified** would mean the across-protein
ordering tracks the covariate while the within-protein manipulation does not — which would
point to a confound correlated with the covariate rather than to the covariate itself, and
must be reported that way rather than as support for the mechanism.

### 6.9 CR9114's covariate was revised before any arm ran, and what that does to its standing

The CR9114 value in §6.5 changed from **0.622x to 0.584x** on 2026-10-31, before any arm ran
on it, because the curated workbook this protocol had been reading excludes a whole quadrant
of the library (D-8). The full correction is in `docs/cr9114-reprep.md`; the moved values are
in §5.1, §6.3, §6.5, §8.2, §9.1, §9.3, §9.4 and §9.5, and the regenerated
`cross_protein_corpus_stats.json`. Two things follow that matter more than the numbers.

#### 6.9.1 The old figure was an overestimate, and that was predicted quantitatively first

Before the full table was obtained, a depth-blind reweighting bound predicted the corrected
enrichment at **0.5622x**. The truth is **0.5845x** — the right direction, **3.8% relative**.
The error's direction was forced: the workbook's filter removes deep genotypes
preferentially, deep genotypes are where CR9114's hits live (hit rate rises monotonically with
depth to 13.79% at depth 14), so a clean-only shell ratio had to *understate* the deep shell's
richness and therefore *overstate* the shallow shell's enrichment. It did, by 6.0% relative.

**Methodological note worth keeping.** Two independent routes to the same quantity agreeing
to 3.8% is the best available evidence that the re-preparation introduced no new error, since
a mistake in either would have to be matched by a compensating mistake in the other. It also
establishes **reweighting as a usable bound when a full table cannot be obtained** — not a
substitute for the full table, but a check that is informative rather than decorative, and one
this protocol should reach for the next time a corpus file is suspect.

#### 6.9.2 CR9114 can now carry a pre-registered ranking — but its position is anticipated, not predicted

**It can carry a ranking.** Every shell in the re-prepared table is **≥98.3% complete**, and
the residual incompleteness is **non-monotone in the opposite direction**, so it cannot
manufacture a depth gradient of the kind the covariate measures. The §6.5 ordering uses
**0.584x** at the primary cut and CR9114 stays at rank 5 of 7 landscapes and rank 3 of 4
proteins — the position is unchanged, now on trustworthy numbers.

**Its position is anticipated, not predicted** — but **the citation given for that in the
first version of this section was the wrong figure, and read literally it points the other
way.** The correction, recorded rather than quietly swapped:

- **What was cited: Phillips Fig 5A. Withdrawn.** Fig 5A plots the **maximum** −logK_D per
  mutation count. For **H1 — the only CR9114 landscape in this pin** — that maximum is ≈8.4 at
  depth 0 and already ≈9.4 by **depth 1**, then flat. Read as a statement about this
  protocol's covariate, Fig 5A says CR9114's *shallow* shell is where the good variants are,
  which is the **opposite** of the measured 0.584x depletion. The deep structure visible in
  Fig 5A belongs to the **H3 and influenza B** traces, and the widely quoted "eight mutations
  required" figure is about **influenza B**. Neither landscape is in this corpus. Phillips'
  main text is explicit: *"For H1, only a few mutations are needed to improve from the
  germline affinity."*
- **What is cited now: Phillips Fig 1D, scoped to H1.** Fig 1D plots **boxplots of the
  affinity distribution per mutation count**. Its upper whiskers reach near-ceiling early —
  which is exactly what Fig 5A reports, and why the two are easy to confuse — while the
  **interquartile range only enters the ceiling band from ~depth 9–10**. A *per-shell
  distribution* is the object whose shape this protocol's hit-rate statistic depends on; a
  *per-shell maximum* is not, because one variant cannot set a shell's hit rate.

**The verdict stands and the evidence was misattributed.** CR9114's position is still
anticipated, on Fig 1D.

**Methods note, because this is the second citation error in this chain and both failed the
same way.** The earlier one was the claim that CombinGym's cleaning changed the H1 censoring
value from 6.0 (§9.4); this one was Fig 5A. **Both were caught by reading the figure rather
than the caption or the summary**, and both would have survived any number of re-readings of
the surrounding prose. Where this protocol cites a published figure as anticipating a
covariate, the citation is only valid if someone has **read that figure's axes and scoping**
— which trace, which antigen, maximum or distribution — and recorded what they read. A caption
is not a citation. Applied to §6.5.1's audit: its three load-bearing citations are each marked
in `docs/ordering-anticipation.md` as figure-reads rather than caption-reads, which is why
they are usable.

Both statements hold at once and neither cancels the other:

- fixing the attrition makes the number **trustworthy**; it does not make it **novel**;
- CR9114 remains legitimate held-out evidence about **arm performance**, which no Phillips
  figure speaks to — Fig 1D describes the landscape, not any method's behaviour on it;
- but **H1's ordering may not be presented as though all four positions were blind**. Where
  the ordering is reported, CR9114's is flagged as literature-anticipated. With n = 4 proteins
  and one position anticipated, the ordering's evidential value is weaker still than §6.4's
  power statement alone implies.

Whether the other three positions are likewise recoverable from published depth-versus-fitness
figures was recorded as O-13 and is now **answered — §6.5.1**. All four are anticipated, so
the caveat above understated the problem rather than overstating it.

### 6.10 Two reserve-landscape findings, recorded as preconditions on promotion

Both come from the clean-versus-raw audit (§3.7). Neither affects any pre-registered covariate
today, because neither landscape is in the final test level. Both would matter the moment one
is promoted, and both would be easy to miss then — which is why they are written down beside
the promotion prohibition of §6.5 rather than left in the audit document.

#### 6.10.1 SpCas9's attrition *is* strongly depth-correlated — and harmless

| | |
|---|---|
| retention at depth 2 | 95.8% |
| retention at depth 8 | 47.2% |
| spread | **48.6 points, monotone** |
| cause, at 0 FP / 0 FN | `RFPsg5 ON present AND RFPsg8 ON present` |
| enrichment effect | **none** |

This is the same *shape* as D-8 and emphatically **not** the same *defect*. All 355 missing
rows **lack an on-target value**, so **none are recoverable** (§3.7), and an unmeasured row can
never enter the eligible universe. The covariate is unchanged.

The reason to record it anyway: a depth-correlated retention curve with a 48.6-point spread is
exactly what would be noticed after a surprising SpCas9 result and read as an explanation. It
is not one, and that is settled here rather than argued later.

Two further SpCas9 facts belonging with it, both already load-bearing in §8.3 and §9.1: its
clean sheet is **597 real rows plus 5,151 wholly empty padding rows**, and its library is only
**3.44% saturated**. The padding is why the drop count matters in §3.3, and the saturation is
why §8.3 refuses to run the absolute-budget protocol on it.

#### 6.10.2 SaCas9's clean phenotype is a 3-of-6 average, undocumented upstream

The clean table's `Mean` column was verified as `mean(sg1, sg2, sg3)` on **all 1,296 rows** —
it silently averages **three of the six sgRNAs** the raw table reports.

Row-complete and depth-blind, so there is **no covariate effect**. But it is a
**phenotype-definition divergence from the upstream publication**, not a row-membership
question, and it is a different species of problem from D-8: D-8 changed *which variants* were
present; this changes *what the number means*.

It is also the only instance of its kind anyone has found, and it was found **incidentally**
(§3.7.1). That is the uncomfortable part: the audit was not looking for it, which gives no
information about how many more there are.

**Precondition on promotion, fixed here:** a reserve protein requires **both** the anticipation
audit of §6.5 **and** a check that each clean phenotype column means what its upstream
publication says it means. SaCas9 specifically may not be promoted until the 3-of-6 average is
either justified against the publication or replaced with the full 6-sgRNA definition, with the
choice stated in advance of any run.

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
| A8 | **MULTI-evolve fixed depth schedule** | the published singles → all-pairs → extrapolate schedule, re-implemented to this protocol's budget — see §7.5 |

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

### 7.4 The cross-protein aggregate is invalid until the anchor runs exist

The mitigation in §7.3 item 2 is a **precondition, not a nice-to-have**. Until the
pool-restricted development anchors have been measured on TrpB and GB1, any cross-protein
number confounds two things that moved at once:

- the **target effect** — a different protein, library geometry and shell structure;
- the **mode effect** — pool-restricted selection with GFlowNet enrichment skipped.

A single difference cannot be attributed to the first while the second is unmeasured. So:

> **No equal-protein-weighted aggregate (§12.2 item 5), no H1 test (§6.4) and no gate decision
> (§12.4) may be computed or reported before the pool-restricted TrpB and GB1 anchor runs
> exist.** Per-target, per-arm tables may be produced and circulated in the interim, each
> carrying the words "pool-restricted; mode effect not yet separated from target effect".

**H2 (§6.7) is the one exception, and for a reason, not as a convenience.** Its statistic is a
difference of differences within a single library: blue and red run under the identical
pool-restricted mode, on the identical 8,192 genotypes, with the identical per-seed initial
sample. Whatever the mode effect is, it enters both halves of `D_s` and cancels. H2 is
therefore valid before the anchor runs land — which is a further reason it, and not H1, is
reported first (§6.8). The same cancellation does **not** apply to any cross-target contrast,
because there the mode and the target move together.

**Anchors.** Those runs are in progress on held-out seeds 48–95, same 25 + 4×25 protocol. The
GFlowNet-enabled held-out figures they will be differenced against, raw top-1% hits, are
fixed here so the comparison cannot be re-baselined afterwards:

| landscape | reliability | frontier |
|---|---|---|
| TrpB | 45.21 | 40.50 |
| GB1 | 15.02 | 37.21 |

The mode effect is `pool-restricted − GFlowNet-enabled` on each of those four cells, estimated
with the seed-clustered contrasts of §11. It is reported **separately and before** any target
effect, and every cross-protein table states its magnitude in the caption. If the mode effect
on either development landscape is itself larger than the material-regression threshold of
§12.3, then the pool-restricted configuration is a different method from the one the
development results describe, and the write-up must say so rather than treating the two as
one engine.

Note what the GB1 column already shows: reliability 15.02 against frontier 37.21. The two
VeritasBio modes are not close to each other on GB1 even before the mode effect is measured,
which is a further reason the two arms are reported separately throughout and never merged
into "the engine".

### 7.5 A8 — the published fixed depth schedule, added now rather than later

**MULTI-evolve** (Tran et al., *Science* 2026, DOI
[10.1126/science.aea1820](https://doi.org/10.1126/science.aea1820); code at
`github.com/ArcInstitute/MULTI-evolve`) is a **published, hand-designed fixed depth schedule**:
identify ~15–20 single beneficial substitutions, measure **all pairwise combinations** of
them, train on singles plus doubles, and predict combinations out to **depth 7**, testing as
few as 9 predicted multi-mutants.

It is, in effect, the simpler fixed-depth heuristic that kill criterion §13.5 item 4 already
names — published, in Science, with wet-lab results, in the past year. The honest thing is to
put it in the arm grid **before** the adaptive version is measured, not after it has won. A
fixed schedule that matches depth-adaptive allocation would kill the adaptive version under
§13.5.4, and that comparison should be available the first time the question is asked.

Registered constraints on A8, so the comparison is not quietly rigged in either direction:

1. **It is a re-implementation to this protocol's budget, not the published method.** The
   published schedule assumes a singles screen plus ~100–200 pairwise measurements; this
   protocol's whole budget is 125. The adaptation is: spend round 1 on singles ranked by the
   revealed initial sample, round 2 on pairwise combinations of the best singles found, and
   rounds 3–4 on surrogate-predicted deeper combinations. The exact rule is fixed in the run
   manifest on first use and may not change afterwards, under the same lock as A7 (§7.1).
2. **This is not a reproduction of MULTI-evolve's published results** and must never be
   described as one. It is a budget-matched fixed depth schedule in the spirit of the
   published one. Any table carrying A8 says so.
3. **Where the library forbids the schedule, the cell is absent, not substituted** (§7.2).
   Pool-restricted binary libraries may not contain the required pairwise combinations of the
   chosen singles; where the enumerated pool lacks them, the shortfall is recorded.
4. A8 is a **baseline**, not a VeritasBio arm. It enters the "effects versus the stronger
   baselines" row of §12.2 item 2 alongside random forest, ALDE and A7.

### 7.6 Named threat to the baseline comparison: ALDE is run far below its own budget

ALDE as published uses a random batch of **96** then **4 rounds of 96** = **480**
measurements, averaged over 70 random initializations, on GB1 and TrpB. This protocol's 125
measurements are **26% of that demonstrated budget**, and the batch size is 25 rather than 96.

Three consequences, registered as a threat rather than defended:

1. **ALDE's deep surrogates run far below any training size it demonstrates.** A DKL or DNN
   ensemble fitted on 25, 50 or 75 revealed measurements is outside the regime the paper
   evidences. A poor ALDE result at these budgets is at least partly a statement about the
   budget, not only about the method.
2. **ALDE publishes no budget sensitivity analysis.** Its reported trajectories hold batch
   size fixed at 96 throughout, so they answer "performance at a smaller *cumulative* budget",
   not "performance at a smaller *batch size*". Nothing in the paper licenses extrapolation
   to a 25-wide batch.
3. Therefore **the ALDE comparison must not be the sentence that carries the claim.** ALDE is
   a reference baseline (§7); ridge (A2) is the designated reference for every effect size,
   and that is deliberate.

**Mitigation, pre-specified.** The full arm grid is **also** run at ALDE's own
**96 + 4×96 = 480** schedule on **TrpB and GB1**. Both are development data (§4.1), so this
costs nothing in held-out evidence, and it is the only way to say whether an ALDE deficit at
125 is the method or the budget. Both schedules are reported side by side wherever ALDE
appears. Any external ALDE comparison carries §1.1's protocol-deviation sentence and, where
the 125-measurement schedule is quoted, the 26%-of-demonstrated-budget figure as well.

This threat applies in reduced form to every arm with a deep surrogate, including the
VeritasBio arms, and that is noted rather than treated as an ALDE-specific excuse.

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
| `cr9114_h1` | 65,094 | 0.19% |
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
  `cr9114_h1` above 0.5% (1% of its 65,094 variants is 651); where it binds, the cell is
  reported as capped with the fraction
  actually reached, never as the nominal fraction. The small libraries — which are the reason
  this analysis exists — are unaffected by the cap.

### 8.3 Targets where the absolute protocol is not usable

Pre-specified, from the universe sizes above rather than from any result:

- **`spcas9` (reserve, 597 variants).** 125 measurements is 20.9% of the universe. The
  absolute-budget protocol is **not** run; only the fixed-fraction analysis is, and the
  target is excluded from the equal-protein-weighted aggregate of §12.2.
- **`cr6261_h1` and `cr6261_h9` (1,887 and 1,842 rows; 1,886 and 1,841 mutants).** Run under
  both schedules, with the share-of-universe figure attached in the same table. Under §4.5
  they are **excluded from the pooled aggregate** and analysed as a separate fixed-fraction
  stratum: at 6.63% of the library the absolute schedule's budget exceeds CR6261's entire
  top-1% hit set of 19 variants, so the endpoint cannot separate an oracle from a good method.
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
| `cr9114_h1` | 65,094 | 1% | 651 | 1 |
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
| `cr9114_h1` | 1% | 651 | **365** | 143 | 9.27 | 16 |
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
  CreiLOV's 1,655 (39%) and 365 of CR9114's 651 (56%). Absolute diversity-adjusted counts are
  therefore **not comparable across library geometries**, which is the second reason §12.2
  aggregates ceiling shares.

### 9.4 Censoring, measured per target

| landscape | floor / min | mass at min | ceiling / max | mass at max | reading |
|---|---|---|---|---|---|
| `eqfp611_blue` | 0.0909 | 0.012% | 1.608 | 0.012% | effectively uncensored |
| `eqfp611_red` | 0.0252 | 0.012% | 1.692 | 0.012% | effectively uncensored |
| `creilov` | 620.3 | 0.0006% | 15,690 | 0.0006% | uncensored; 165,428 distinct values in 165,428 rows |
| `cr9114_h1` | 7.00 | **2.57%** (1,675 rows) | 9.835 | 0.002% | mild floor at the −log10(Kd) detection limit |
| `cr6261_h1` | 7.00 | **11.34%** | 9.729 | 0.053% | censored from below |
| `cr6261_h9` | 7.00 | **20.20%** | 10.15 | 0.054% | censored from below |
| GB1 (dev) | 0.00 | **19.74%** | 8.762 | 0.0007% | censored from below |
| HER2 (dev) | 0.10 | **51.7%** | — | — | heavily censored from below |

**Correction to a claim about these floors — now settled per column.** It has been stated
elsewhere that "the raw bnAb tables" encode non-binders at 6.0 and that CombinGym's cleaning
step therefore changes the censoring value. **The plural was wrong twice over, and so was the
claim:**

| column | 6.0 floor? | minimum | raw vs clean | curated / scored here? |
|---|---|---|---|---|
| `h1_mean`, CR9114 | **no — zero rows at 6.0** | 7.0 | **agree** on floor and mass | scored |
| `h1_mean`, CR6261 | **no — zero rows at 6.0** | 7.0 | **agree** on floor and mass | scored |
| `h9_mean`, CR6261 | no | 7.0 | agree | scored |
| `h3_mean`, CR9114 | yes | 6.0 | — | **not curated, not scored** |
| `fluB_mean`, CR9114 | yes | 6.0 | — | **not curated, not scored** |
| any such column, CR6261 | **column does not exist** — `bnAbs_CR6261_raw.csv` has none | — | — | — |

So: the 6.0 floor exists **only** in CR9114's `h3_mean` and `fluB_mean`, neither of which is
curated or scored, and `bnAbs_CR6261_raw.csv` has no such column at all. **No cleaning step
changes any censoring value on any bnAb landscape**, and the H1 and H9 figures in the table
above are the raw tables' own. The "cleaning changed the censoring value" caveat applies to
nothing this protocol scores, and it is one of the two citation failures the methods note in
§6.9.2 was written about.

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
  10 members — on this corpus that is `cr9114_h1` (66 at 0.1%, 7 at 0.01%, the latter below
  the limit and therefore **not** resolvable), `creilov` (166 and 17), and the development
  landscapes. Elsewhere the cell reads
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

### 10.1 Mean selected depth is a reported quantity, in two units

Library width varies from **4 to 16 sites** across this corpus. A depth bias measured on
four-site libraries cannot be assumed to transfer as an absolute number rather than as a
proportion of sites, and the two readings diverge sharply as libraries widen (§6.7.4 P2). So
the question is settled by measurement:

> **For every arm, every target and every round, mean selected depth is reported twice: as an
> absolute depth, and as a fraction of that library's site count.** The same two figures are
> reported for mean *proposed* depth and for the depth distribution, not only its mean.

Pre-specified read-out, done **before** P2 is scored and before any cross-protein comparison
is interpreted:

1. Regress the reliability arm's mean selected absolute depth on library site count across all
   targets. A slope near zero supports the **absolute** reading; a slope near 0.65 — the
   development value, 2.6 of 4 sites — supports the **proportional** reading. Report the slope
   with its seed-clustered interval and the per-target values, not just the fit.
2. State the verdict, or state that the data does not separate the two readings, before
   scoring P2.
3. The same read-out is produced for every other arm, because the baselines' depth behaviour
   (3.48–3.60 on four-site libraries, i.e. 0.87–0.90 of sites) is subject to exactly the same
   question and is the comparison that gives the reliability figure meaning.

This costs no extra compute: the rows are already required above, and the analysis is a
regression over values the runs will have produced anyway. It is written down because an
unmeasured premise that quietly becomes an assumption is how §6.7.4 nearly acquired a point
prediction resting on nothing.

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

The full grid is 8 arms × 5 budgets × 9-plus endpoints × 7 landscapes. Treating that as one
family would guarantee both false positives and, with Bonferroni, no power at all. The
pre-specified split:

- **Primary family.** `{reliability vs ridge}` × `{5 budgets}` × `{primary endpoint}`, per
  protein. That is 5 tests per protein and 20 across four proteins. **Holm–Bonferroni** at
  family-wise alpha 0.05, applied **within each protein** (5 tests), with the per-protein
  adjusted p-values reported; and separately across the 20, reported as the stricter figure.
  Both are shown; neither is chosen afterwards.
- **The H1 test** (§6.4) is a single pre-specified test with its own alpha of 0.05,
  one-sided, outside the primary family.
- **The H2 test** (§6.7) is a second single pre-specified test, also alpha 0.05 one-sided,
  also outside the primary family: one contrast, `mean(D_s)` on the primary endpoint at the
  125-measurement budget, on the eqFP611 blue-red instance. The CR6261 replication and the
  other four budgets are reported alongside it and are **not** additional confirmatory tests;
  they are labelled as replication and sensitivity. H1 and H2 are not corrected against each
  other: they test different hypotheses on different designs, each gets one pre-specified
  test, and both results are reported whatever they show (§6.8).
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
2. **Effects versus the stronger baselines** — random forest, ALDE where compatible (§7.2,
   with the budget threat of §7.6 attached), the lightweight configuration A7, and the
   published fixed depth schedule A8 (§7.5). An improvement over random alone is not a result,
   and an improvement that A8 also achieves is not a result for the engine.
3. **The fraction of independent proteins improving**, with the count and the denominator
   written out (for example "2 of 4 proteins; 3 families"), never as a bare percentage.
4. **The distribution and magnitude of regressions** — every regression, its size, its CI,
   and the worst one named in the summary sentence. Not a count of regressions: their sizes.
5. **An equal-protein-weighted aggregate effect**, computed on **yield as a share of the
   achievable ceiling** (§9.3) because tiers and library geometries differ across targets
   (§9.1, §9.3). One protein, one weight. **CR6261 is excluded from this aggregate** on the
   suitability grounds of §4.5, so the aggregate is over **3 proteins** — eqFP611, CreiLOV,
   CR9114 — and every pooled number states that denominator. CR6261 is reported separately as
   a fixed-fraction stratum.
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
4. A **simpler fixed-depth heuristic matches it** — a fixed allocation to depth 3–4, simply
   selecting at the depth where hits have so far been found with no adaptive estimator at all,
   or **arm A8, the published MULTI-evolve singles → all-pairs → extrapolate schedule
   (§7.5)**. A8 is in the grid from the start precisely so this criterion can fire on the
   first measurement rather than after the adaptive version has been declared a win.
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
2. `fetch`, `record`, `verify` the corpus; commit the lock file. **Done; 29/29 files verified
   — 20 at the freeze, plus the raw CR9114 table, plus the remaining eight raw tables, all at
   the same unchanged pin (§3.4).**
3. `stats` on every target; the covariates of §6.5 and the ceilings of §9.3 are published
   above. Label-only, no arms. **Done 2026-10-31.** Recompute TrpB's anchor with the same
   code path (§16 O-4).
4. Run the budget-25 identity check (§8.1) on every target. A failure voids the run.
5. Run the development-level plumbing cross-check on GB1 (§4.3).
6. Run all arms on all targets, both budget schedules, 48 seeds. Capture §10 rows.
7. **Measure the pool-restricted TrpB and GB1 anchors** (§7.4) and report the mode effect. In
   progress on seeds 48–95. Steps 9 and 10 are blocked until this exists.
8. **Test H2** (§6.7). Not blocked by step 7 — the mode effect cancels in its
   difference-of-differences (§7.4) — so it can be reported as soon as step 6 covers eqFP611.
9. Run the arm grid at ALDE's own **96 + 4×96** schedule on TrpB and GB1 (§7.6). Development
   data, so this costs nothing held-out.
10. Analyse per §11, report per §12 in the order fixed by §6.8, with no changes to §§2, 6, 8,
    9, 11 or 12.
11. Only then consider §13.

Steps 6, 7 and 9 are independent of each other and may run concurrently. Step 8 needs only
step 6's eqFP611 cells.

---

## 16. Open items this document does not resolve

| # | item | consequence |
|---|---|---|
| O-1 | ~~`CreiLOV_clean.xlsx` not fully retrieved.~~ **Resolved on 2026-10-31** before the freeze: all 20 pinned files fetched, verified against the commit tree and characterised. Nothing in §§3.4, 5.1, 6.5, 9.1, 9.3 or 9.4 is outstanding. | none |
| O-2 | **RhlA's library is measured only to depth 3** of its 11 sites (2,346 of 663,552 variants, 0.35% saturation), so the shallow-shell covariate is undefined there under the §6.3 degeneracy guard and no cut can be made to discriminate. An earlier draft of this document recorded RhlA as blocked for the primary endpoint; that was wrong — its `genotype` column is fixed-width (D-6) and the endpoint computes. | RhlA is excluded from the H1 test only. Its other endpoints are reported. |
| O-3 | **ALDE compatibility on this corpus is unprobed.** The binary per-site alphabets and the 121-/1,053-residue genotype representations are not the saturated-site encoding its published code assumes. | Likely absent cells, to be documented per §7.2 rather than substituted. |
| O-4 | ~~TrpB's shallow-shell figures were not recomputed here.~~ **Resolved on 2026-10-31.** TrpB has now been reproduced by a **second, independent implementation** — a parallel agent's own genotype canonicalisation and shell profiler, run against the real landscape files — giving shallow 168/2,242 = **7.4933%**, deep 1,432/157,757 = **0.9077%**, per-depth 15.79% / 7.20% / 2.33% / 0.608% at depths 1–4 and a depth-4 maximum of 1.0000. The same implementation gives GB1 shallow 72/2,167 = **3.3226%**, deep 1,422/147,193 = **0.9661%**, depth-4 maximum 8.7620. `freeze_cross_protein_corpus.py stats` reproduces GB1's figures **exactly** under the §6.3 depth-0 convention. Both anchors now come from agreeing implementations; the earlier 2,168 / 3.3210% was the depth-0 convention difference and nothing else. | none |
| O-5 | **The upstream corpus has seven recorded discrepancies or inconsistencies** (§3.6, D-1 to D-7), two of them unexplained one-row differences and one — SpCas9, 948 claimed against 597 usable — a factor of 1.6. | The clean file wins where they disagree; D-5 means RhlA's provenance is unestablished. |
| O-15 | **Whether the four reserve proteins' positions are also pre-published is unchecked**, and §6.5's standing prohibition blocks promoting any of them into the H1 test until it is. Promoting four pre-published positions would move the exact attainable one-sided p from 1/24 to 1/5,040 while adding **zero** blind predictions. | The only sanctioned route to n = 7 runs the §6.5.1 audit first and records each verdict. A promoted test must report how many of its n positions are blind alongside its p-value. **A promoted protein now needs two clearances**: that audit *and* the phenotype-provenance check of O-16 (§6.10). SaCas9 specifically is blocked until its 3-of-6 sgRNA average is justified or replaced. |
| O-6 | **H1's power is marginal by construction and it will probably be recorded as unresolved** (§6.4). With n = 4 proteins the smallest attainable exact one-sided p is 0.042, so only a perfect ordering reaches alpha; the covariate's top three landscape values are a near-tie; and the demonstration is that excluding one wild-type row from a denominator swapped ranks 1 and 2 on a 0.02x margin. H1 is additionally ranking the *magnitudes of four losses*, since every independent value sits below GB1's 3.32x — and, per O-13, at most one of those four positions is blind. | H1 is kept but is no longer the only formal test of the mechanism: **H2 (§6.7)** tests the mechanism at n = 48 seeds on a design where only the covariate moves, and now carries more of the study's weight (§6.8). Expanding to the reserve proteins would take H1 to n = 7 **only subject to the §6.5 anticipation prohibition** (O-15). |
| O-9 | **A8 is a budget-matched re-implementation, not MULTI-evolve** (§7.5). The published schedule assumes a singles screen plus ~100–200 pairwise measurements; this protocol's entire budget is 125. The adaptation rule is fixed on first use but has not been written or tested. | A8 results are evidence about a fixed depth schedule at a 125-measurement budget, never a reproduction of the published result. |
| O-10 | **ALDE runs at 26% of its own demonstrated budget here** (§7.6): 125 against 96 + 4×96 = 480, and at batch 25 rather than 96. ALDE publishes no batch-size sensitivity analysis, so a deficit at these budgets is partly a statement about the budget. | ALDE must not carry the claim; ridge is the reference baseline. Mitigation is the 480-measurement development runs on TrpB and GB1, which have not been run. |
| O-11 | **The CR6261 H1/H9 genotype intersection has not been computed** (§6.7.5). H2's replication instance needs it, because the two antigens have 1,887 and 1,842 rows and the shared-initial-sample construction applies only to genotypes present in both. | Compute and record before running the replication contrast. |
| O-12 | **Whether the policy's selection depth scales with library width is unmeasured — and this matters because it governs whether the mechanism means anything beyond four sites at all, not because it confounds H1.** The depth bias (reliability 2.63 / 2.56; baselines 3.48–3.60) was measured on **two libraries of identical four-site width**, which cannot distinguish an absolute constant from a proportion of sites. **The serious branch: under the proportional reading** (~65% of sites) the policy selects near depth 8.5 of 13 on eqFP611 and 10.4 of 16 on CR9114, i.e. it is **not a shallow exploiter on any library wider than four sites**. The shallow-exploitation mechanism would then describe **only TrpB and GB1 — both development data — and would not apply to a single independent target.** That is a more serious threat to the project than the generalization result itself, because it would void the mechanism's relevance rather than merely weaken the evidence for it. Under the **absolute** reading the policy instead gets effectively shallower as libraries widen (2.6 of 4 sites is 65% of positions, 2.6 of 16 is 16%) and performance should degrade monotonically with site count. **Correction to an earlier draft of this item:** it claimed site count covaries with the enrichment ordering strongly enough to confound H1. Checked against the frozen per-protein values — CR9114 16 sites / **0.584x**, CreiLOV 15 / 1.875x, eqFP611 13 / 1.439x, CR6261 11 / 0.105x — the Spearman correlation between site count and enrichment is **rho_s = +0.400, exact one-sided p = 9/24 = 0.3750**. (The CR9114 re-preparation of §6.9 moved that value from 0.622x to 0.584x and left the rank pattern, and therefore rho_s and p, **unchanged**.) The two orderings differ at one position and it is the decisive one: **CR9114 is the widest library in the corpus (16 sites) yet only third of four in enrichment.** That single inversion decorrelates width from the covariate enough that **H1 is not seriously confounded by it**. The confound framing was an overstatement and is withdrawn; the premise-undercutting branch above is the reason this item is a priority. | **Testable from the §10 rows at zero extra compute**, by the read-out pre-specified in §10.1, and it must be done before P2 (§6.7.4) is scored. Site count is still reported alongside enrichment in every H1 table (§6.5) — as disclosure, not as a correction for a confound that the data does not support. |
| O-7 | **The size of the mode effect is still unknown** (§7.3). The mitigation — running reliability and frontier pool-restricted on TrpB and GB1 — is **now in progress on held-out seeds 48–95**, against the fixed GFlowNet-enabled anchors of §7.4 (TrpB 45.21 / 40.50, GB1 15.02 / 37.21). | Until those runs land, the aggregate, the H1 test and the gate decision are all blocked by §7.4. |
| O-13 | ~~How much of the H1 ordering is recoverable from published literature is unchecked.~~ **Answered, before any run** (`docs/ordering-anticipation.md`, §6.5.1). **At most one of H1's four positions is genuinely blind, and arguably zero.** eqFP611 — Poelwijk Fig 6a, anticipated. CreiLOV — Chen Fig 3c, which publishes the fraction of each depth shell beating WT and is the only published object of this protocol's own statistic type found in the audit, anticipated. CR6261 — Phillips Fig 1E/1F, partly anticipated. CR9114 — anticipated, on **Fig 1D scoped to H1**, after the Fig 5A citation was withdrawn as wrong and pointing the other way (§6.9.2). TrpB and GB1 are **clean negatives**, so the derivation could not have been read off the literature; only the test set could. | H1 is **not void**: the covariate column is published, the outcome column is not, so the relationship and the §6.5 across-the-board-loss prediction remain blind. What is struck is any presentation of **the ordering itself** as a finding of this project. Consequences recorded: anticipation columns in §6.5, the standing prohibition on n-expansion without an audit (§6.5), and H2 carrying more weight (§6.8). |
| O-14 | ~~The clean/raw divergence was only checked for CR9114.~~ **Closed, before any run** (`docs/corpus-audit.md`, §3.7). **The D-8 defect is unique to CR9114. Not one pre-registered covariate moves. The freeze holds. No landscape needs re-preparing.** The decisive quantity is **recoverable rows** — missing from clean yet carrying a real measurement in raw — which is **16,253 for CR9114 and exactly 0 everywhere else**; an unmeasured row can never enter the eligible universe (§3.3), so its absence cannot shift a covariate however depth-correlated it looks, and recomputing every covariate over the raw universe gives a delta of **exactly 0.0000** bar CR9114. A composite sweep over every single-site and site-pair rule of the CR9114 form found nothing at 0 FP / 0 FN elsewhere (CR6261's best candidate off by 423 kept rows, SpCas9's by 361), and **eqFP611 and SaCas9 are combinatorially complete**, so no filter can hide in either. The negatives are credible because the same code finds the positive: it reproduced all seven frozen §6.5 covariates digit for digit, re-derived CR9114's filter at 0 FP / 0 FN, and reproduced the 0.5622x reweighting bound. | None for the final-test set. Two reserve findings became promotion preconditions (§6.10). **Scope caveat: the audit covered row membership, value fidelity and depth structure, not phenotype-column choices** — see O-16. |
| O-16 | **No clean table's phenotype *column* choice has been checked against its upstream publication** (§3.7.1). The audit covered row membership, value fidelity and depth structure only. The one divergence found — SaCas9's `Mean` silently averaging 3 of the 6 sgRNAs the raw table reports (§6.10.2) — surfaced **incidentally**, which gives no information about how many more exist. For the four final-test targets the columns in use are the ones the registry names and the §9.4 figures describe, which is a weaker check than a publication comparison. | **A clean bill on row membership is not a clean bill on phenotype definition.** A phenotype-provenance check is a precondition on promoting any reserve protein (§6.5, §6.10.2), and should also be run on the four final-test columns before the results are written up. It would not move a covariate, but it could change what a covariate *means*. |
| O-8 | **No wet-lab validation of anything here.** Every figure is a retrospective simulation over a fully-measured library. | Nothing in this study is evidence about a real campaign. |

---

## 17. Amendments

Appended, dated, with a reason. The original text above is not edited in place. Any amendment
written after cross-protein results have been examined must say so in its first sentence.

| date | section | change | reason |
|---|---|---|---|
| 2026-10-31 | — | initial freeze: corpus pinned and verified (20/20 files), covariates and ceilings computed for all 7 final-test and 6 reserve landscapes, H1 and the §13 intervention registered | — |
| 2026-10-31 | §6.3 | fixed the depth-0 convention: the wild type is retained in the eligible universe and the tier, excluded from both shells. All shell figures recomputed under it; landscape-level ranks 1 and 2 swapped (CreiLOV now above eqFP611 blue), the per-protein ordering unchanged | two implementations had disagreed by one row on GB1 (2,167 vs 2,168), the cause being a literal `WT` row; convention drift between implementations is what the freeze exists to prevent |
| 2026-10-31 | §16 O-4 | closed: TrpB's shell figures reproduced by a second independent implementation; both development anchors now agree | the gap flagged at the initial freeze was closed by parallel work, not by this document |
| 2026-10-31 | §7.4 | added: the cross-protein aggregate, the H1 test and the gate decision are blocked until the pool-restricted TrpB/GB1 anchor runs exist, and the mode effect must be attributed separately from the target effect | the §7.3 confound needed an enforceable precondition, not just a note; the anchor runs are now in progress on seeds 48–95 |
| 2026-10-31 | §6.7, §6.8, §6.4, §11.3 | added **H2**, the within-protein mechanism hypothesis (eqFP611 blue against red, difference of differences paired over 48 seeds, with CR6261 H1/H9 as a weaker replication), a fixed interpretation order for H1 against H2, and the statement that H1's power is marginal by construction with the rank-1/2 flip as its demonstration | H1 needs a perfect ordering of four items whose top three move on one-row conventions, so it will probably be unresolved either way; the within-protein contrast was already identified as the cleanest covariate manipulation in the corpus and was a remark rather than a test. H2 is exempt from the §7.4 block because the mode effect cancels in its difference of differences |
| 2026-10-31 | §7, §7.5, §13.5, §12.2 | added arm **A8**, the published MULTI-evolve fixed depth schedule (Tran et al., Science 2026, 10.1126/science.aea1820), budget-matched to this protocol, with the registered caveat that it is not a reproduction of the published result | §13.5 item 4 kills depth-adaptive allocation if a simpler fixed-depth heuristic matches it, and the strongest such heuristic is already published; it belongs in the grid before the adaptive version is measured, not after it wins |
| 2026-10-31 | §6.7.4, §6.7.5, §6.5, §10.1 | split H2's prediction into **P1** (directional, robust, H2's actual claim) and **P2** (the near-point prediction on red, now explicitly **conditional** on the depth bias being absolute rather than proportional to library width); softened the falsification clause that had made a healthy red yield fatal to the mechanism regardless of blue; required mean selected depth to be reported in both absolute and fractional units for every arm and target, with a pre-specified read-out that settles absolute-versus-proportional before P2 is scored; required site count alongside the H1 enrichment ordering | the premise "the policy selects at mean depth ~2.6" was measured only on two **four-site** libraries. eqFP611 has thirteen. Absolute (~2.6, 20% of sites) and proportional (~8.5, 65% of sites) readings make **opposite** predictions about red, whose hits peak at depth 8, and nothing measured distinguishes them. As written, a width-transfer effect would have been misrecorded as a refutation of the mechanism |
| 2026-10-31 | §16 O-12 | added: whether selection depth scales with library width is unmeasured, testable from the §10 rows at no extra compute | if the bias is absolute the policy gets effectively shallower as libraries widen; if proportional, it is not a shallow exploiter on any library wider than four sites. Either answer changes the interpretation |
| 2026-10-31 | §16 O-12, §6.5 | **corrected**: the claim that site count covaries with the enrichment ordering strongly enough to confound H1 was an overstatement and is withdrawn. Measured on the frozen per-protein values, rho_s = **+0.400**, exact one-sided p = **9/24 = 0.3750**, because CR9114 is the widest library (16 sites) yet third of four in enrichment. O-12 now leads with the premise-undercutting branch instead: under a proportional depth bias the policy is not a shallow exploiter on any library wider than four sites, so the mechanism would describe only TrpB and GB1 — both development data — and no independent target | the confound risk and the premise-undercutting risk are not the same risk, and the item led with the weaker one. The first is not supported by the data; the second is unmeasured and would void the mechanism's relevance rather than merely weaken its evidence |
| 2026-10-31 | §3.2 | recorded that **no CR9114/H9 landscape exists** — Phillips et al. assayed CR9114 against H1, H3 and influenza B, and H9 for CR6261 only — and that it may not be listed as unobtainable, pending or absent | the CR6261 H1/H9 pair invites the symmetric assumption; an experiment never performed is not a missing file |
| 2026-10-31 | §5.3 | added the whole-study count: **seven proteins, five structural families**, with trastuzumab, CR9114 and CR6261 as three proteins in **one** immunoglobulin V-domain family, and the acquisition contributing **+2** families the project had no data in. Also recorded that of §5.2's three independent families, the Ig V-domain is **already represented in development data** by trastuzumab, so transfer to a *new* family rests on **two** families, not three | the earlier count did not place trastuzumab in the same family as the two bnAbs, which both understated the study-wide family count and overstated the independence of the antibody results at the family level |
| **2026-10-31** | §3.6 D-8, §6.9 | **CR9114 re-prepared and its covariate revised, before any arm ran on it.** `bnAbs_CR9114_H1_clean.xlsx` excludes every genotype carrying both F29S and Y106S — 16,384 genotypes, one quadrant of the 2^16 space — regardless of affinity, with retention falling in closed form as 1 − k(k−1)/240 and predicting 45.0% / 12.5% / 0% at depths 12 / 15 / 16 against those observed exactly. The fully matured 16-mutation antibody was entirely absent. CR9114 is now scored against the raw table at the same pinned commit | a depth-correlated filter on the covariate's own axis is the one defect that could manufacture the gradient H1 measures. Found and fixed before any arm ran, so the covariate is revised rather than the result reinterpreted |
| **2026-10-31** | §5.1, §6.3, §6.5, §8.2, §9.1, §9.3, §9.4, §9.5, §16 O-12, `cross_protein_corpus_stats.json` | **values moved**: eligible universe 48,841 → **65,094** rows (65,093 mutants); top-1% hits 489 → **651**; tier boundary 9.6077550214 → **9.6088080625**; primary enrichment 0.622x → **0.584x** (−6.0% relative); d≤2 enrichment **0.000x, unchanged** (shell of 134); mean hit depth 8.91 → **9.27**; max observed depth 15 → **16**; d≥2 ceiling 280 → **365** (+30.4%); d≥3 ceiling 118 → **143**; 125-measurement budget share 0.26% → **0.19%**; floor mass 1.91% → **2.57%** (1,675 rows); top 0.1% / 0.01% tiers 49 / 5 → **66 / 7**. The fractional-budget cap now binds on `cr9114_h1` above 0.5% rather than 1%. Regenerated the stats artifact; **only the `cr9114_h1` block changed**, the other 13 landscapes are byte-identical. The per-protein H1 ordering and `rho_s = +0.400, p = 9/24` are unchanged | the stats artifact is generated from the corpus, so it had to be regenerated from the raw table rather than hand-edited; its generator read the clean workbook and would otherwise have reproduced the hole |
| **2026-10-31** | §3.4, §3.5, `cross_protein_corpus.lock.json`, `freeze_cross_protein_corpus.py` | **the pin did not change; the lock gained one entry.** `Data/DMS/Raw/bnAbs_CR9114_raw.csv` is now pinned (git blob `8ac7352f…`, sha256 `ee31f8fc…`, 12,258,680 bytes) at the same commit `c5dccd51`. The `corpus` block is byte-identical and **no existing entry changed** — all 20 were recomputed and came out identical — taking the lock from 20 to 21 files. `verify` now also **fails if any pinned path, or any file a landscape reads, is absent from the lock**; it reports 21/21 | the lock recorded only upstream blobs and no derived table, so a landscape repointed at an unlisted file would have been scored against unhashed bytes while `verify` still exited 0. Pinning the raw table closes that for CR9114 and the new guard closes the class |
| **2026-10-31** | §6.9.1 | recorded that the pre-fix **depth-blind reweighting bound predicted 0.5622x against the measured 0.5845x** — right direction, 3.8% relative — and that two independent routes agreeing is the best available evidence the re-preparation introduced no new error | it validates reweighting as a usable bound when a full table cannot be obtained, which is a reusable method rather than a one-off observation |
| **2026-10-31** | §6.9.2, §16 O-13 | recorded that CR9114's position in the H1 ordering is **anticipated, not predicted** — ~~Phillips Figure 5A already publishes max affinity against depth for it~~ **(citation superseded two rows below: the correct figure is Fig 1D; Fig 5A was wrong and pointed the other way)** — so fixing the attrition makes the number trustworthy but not novel, and the ordering may not be presented as though all four positions were blind | an ordering recoverable from published figures is not a test; whether the other three positions are likewise recoverable is unchecked (O-13) |
| **2026-10-31** | §6.5, §6.5.1, §16 O-13 | **O-13 answered, before any run: H1's covariate ordering is substantially pre-published.** At most one of its four positions is genuinely blind, arguably zero. eqFP611 — **Poelwijk Fig 6a** (Nat Commun 10:4213, 10.1038/s41467-019-12130-8); CreiLOV — **Chen Fig 3c** (ACS Synth Biol 12(5):1461–1473, 10.1021/acssynbio.2c00662), the only published object of this protocol's own statistic type in the audit; CR6261 — **Phillips Fig 1E/1F** (eLife 10:e71393, 10.7554/eLife.71393); CR9114 — Phillips **Fig 1D**. TrpB (Johnston 2024) and GB1 (Wu 2016) are clean negatives. Anticipation columns added to both §6.5 tables, with the citations in full | an ordering recoverable from published figures is not a finding. H1 survives because **the covariate column is published and the outcome column is not** — it asserts a relationship, and the per-landscape effect sizes plus §6.5's across-the-board-loss prediction remain entirely blind. **H1's derivation was original and unpublished; its test set is pre-published.** What is struck is any sentence presenting the ordering itself as this project's finding |
| **2026-10-31** | §6.9.2 | **withdrew the Phillips Fig 5A citation for CR9114 and replaced it with Fig 1D scoped to H1.** Fig 5A plots per-shell **maxima**, and for H1 — the only CR9114 landscape in the pin — that maximum is ≈8.4 at depth 0 and ≈9.4 by **depth 1**, flat thereafter, so read literally it says the shallow shell is rich, the **opposite** of the measured 0.584x depletion. Its deep structure belongs to the H3 and influenza B traces and the "eight mutations" figure is about influenza B; neither is in this corpus. Fig 1D's per-shell **distributions** have an interquartile range entering the ceiling band only from ~depth 9–10, which is the object a hit-rate statistic depends on. Added a methods note: this is the second citation error in this chain, both were caught by reading the figure rather than the caption, and a caption is not a citation | the verdict was right and the evidence was misattributed. An anticipation claim is valid only if someone has read the figure's axes and scoping — which trace, which antigen, maximum or distribution — and recorded it |
| **2026-10-31** | §6.7.1 | **re-framed eqFP611 blue-versus-red as a natural experiment Poelwijk built**, not a manipulation this protocol selected: the library is the complete 2^13 space **between two named parents**, mTagBFP2 (blue, depth 0) and mKate2 (red, depth 13) — the pinned WT file being `eqFP611(mTagBFP2).fasta` — and Fig 6a's published colour switch at step ~6 lands exactly on the independently derived floor(13/2) = 6 cut | stronger provenance, not weaker: a covariate manipulation fixed by a 2019 publication **cannot have been chosen after seeing our results**. H2's content was never the ordering but the effect-size difference, which no figure publishes |
| **2026-10-31** | §6.5, §16 O-15 | added a **standing prohibition**: no reserve protein may be promoted into the H1 test until the §6.5.1 anticipation audit has been run on it and its verdict recorded; a promoted test reports how many of its n positions are blind alongside its p-value | promoting four more pre-published positions would move the exact attainable one-sided p from **1/24 to 1/5,040 while adding not a single blind prediction**. It would read as a decisive power gain and be nothing of the kind. The audit is a literature read, not a compute job, and must precede promotion because afterwards the number looks authoritative and the question stops being asked |
| **2026-10-31** | §6.5, §6.8 | recorded that §6.5's across-the-board-failure prediction is **not weakened by the anticipation audit and is not to be softened in response to it** — it is now the main thing H1 tests — and that **H2 carries more of the study's weight** than when the reporting order was fixed | no published figure says anything about how the reliability policy performs against one-hot ridge, so that prediction stays blind. H2 first was already right for power; it is now also right for originality |
| **2026-10-31** | §3.7, §16 O-14 | **O-14 closed, before any run: a clean bill of health.** The D-8 defect is **unique to CR9114**; not one pre-registered covariate moves; the freeze holds; no landscape needs re-preparing. Recorded the reasoning and not just the verdict: the decisive quantity is **recoverable rows** (missing from clean, real measurement in raw) = **16,253 for CR9114, exactly 0 everywhere else**, and an unmeasured row can never enter the eligible universe, so its absence cannot shift a covariate however depth-correlated it looks — recomputing every covariate over the raw universe gives a delta of **exactly 0.0000** bar CR9114. Composite sweep over all single-site and site-pair rules found nothing at 0 FP / 0 FN elsewhere (CR6261 off by 423 kept rows, SpCas9 by 361); eqFP611 and SaCas9 are combinatorially complete, so no filter can hide. Negatives credible because the same code reproduced all seven frozen covariates digit for digit, re-derived CR9114's filter at 0 FP / 0 FN, and reproduced the 0.5622x bound | a null from an instrument never shown to detect anything is not evidence; this instrument demonstrably finds the positive. And the recoverable-rows argument closes the question deductively rather than statistically, which is what makes a negative result load-bearing |
| **2026-10-31** | §3.4, §3.5, `cross_protein_corpus.lock.json`, `freeze_cross_protein_corpus.py` | **extended the lock to all nine raw tables, 21 → 29 files.** Additive only: `corpus` block byte-identical, **no existing entry changed or removed**, and every raw table's blob SHA-1 independently cross-checked against the tree SHA `911f19c7…` of §3.1. The commit, date and tree SHA are untouched. `verify` reports 29/29 | the nine raw tables are now **the evidence that the clean tables are faithful** (§3.7). Unpinned, that evidence is unpinned, and a later audit could run against different bytes while `verify` still exited 0 — the same hole the CR9114 episode exposed. Also recorded the narrower rule that survives: **the pin must not change; the lock may only grow** |
| **2026-10-31** | §3.7.1, §16 O-16 | recorded the audit's **scope caveat**: it covered row membership, value fidelity and depth structure, **not** each clean table's phenotype-**column** choices against its upstream publication. **A clean bill on row membership is not a clean bill on phenotype definition** | the one divergence found, SaCas9's 3-of-6 sgRNA average, surfaced **incidentally** rather than by design, which gives no information about how many more exist |
| **2026-10-31** | §6.10, §6.5, §16 O-15 | recorded the two reserve findings as **preconditions on promotion**. **SpCas9**: attrition *is* strongly depth-correlated — 95.8% retention at depth 2 falling monotonically to 47.2% at depth 8, a 48.6-point spread, cause named exactly as `RFPsg5 ON present AND RFPsg8 ON present` at 0 FP / 0 FN — yet **harmless**, because all 355 missing rows lack an on-target value so none are recoverable; also 597 real rows plus **5,151 wholly empty padding rows** against a library 3.44% saturated. **SaCas9**: the clean `Mean` is `mean(sg1,sg2,sg3)` on all 1,296 rows, silently averaging **3 of the 6** sgRNAs the raw table reports — row-complete and depth-blind, so no covariate effect, but a phenotype-definition divergence from the upstream publication. A promoted protein now needs **two** clearances, anticipation and provenance; SaCas9 is blocked until the 3-of-6 average is justified or replaced | both would matter only on promotion and would be easy to miss then. SpCas9's 48.6-point curve is exactly what would be noticed after a surprising result and read as an explanation — it is not one, and that is settled now rather than argued later. SaCas9's is a different species from D-8: D-8 changed which variants were present, this changes what the number means |
| **2026-10-31** | §9.4 | **settled the bnAb floor claim per column**, and the plural was wrong twice over: `h1_mean` has **zero rows at 6.0** on both CR9114 and CR6261 with a minimum of 7.0, and raw agrees with clean on floor and mass; the 6.0 floor exists **only** in CR9114's `h3_mean` and `fluB_mean`, neither curated nor scored; and `bnAbs_CR6261_raw.csv` has **no such column at all**. **No cleaning step changes any censoring value on any bnAb landscape** | an earlier version of this row said only "false for H1", which left the impression that some bnAb landscape somewhere was affected. None is |
| **2026-10-31** | §9.4, §16 O-14 | **corrected an inventory claim**: the raw bnAb tables do **not** floor H1 at 6.0. Raw `h1_mean` has a minimum of 7.0 with 1,675 rows at the floor, identical to clean; the 6.0 floor applies to `h3_mean` and `fluB_mean` only. The "cleaning changed the censoring value" caveat is false for H1 and applies to no landscape this protocol scores. Also opened O-14: the other nine clean workbooks have not been diffed against their raw counterparts | a false caveat about censoring would have been carried into the censored-endpoint reasoning of §9.4; and D-8 was found in the one landscape anyone looked at, which says nothing about the other nine |
| 2026-10-31 | §4.5, §8.3, §12.2 | added: **CR6261 is excluded from the pooled aggregate** on three pre-computable grounds — 125 measurements is 6.63% of its 1,886 mutants against 0.08% on TrpB; its top-1% tier is 19 variants so the budget exceeds the hit set and an oracle scores like a good method; and its H9 covariate is degenerate at zero under both cut definitions. Retained in full as the qualitative demonstration that the covariate can invert, and as H2's replication instance. The pooled aggregate is now over **3 proteins** | a target can be legitimately held out and still be the wrong instrument; the ratio and the tier size are properties of the dataset, not of any result |
| 2026-10-31 | §7.6 | added the named threat that ALDE runs at 26% of its own demonstrated budget (125 against 96 + 4×96 = 480, batch 25 rather than 96, no published batch-size sensitivity), with the mitigation of also running the grid at 96 + 4×96 on TrpB and GB1 | both are development data, so the mitigation costs nothing held-out, and without it an ALDE deficit cannot be attributed to the method rather than the budget |

Every row above is dated **2026-10-31**: the initial freeze and all amendments to date fall on
the same day, in the table's order, which is why the order and not the date distinguishes
them. The rows from "CR9114 re-prepared" onward are the CR9114 round and came last.

**No cross-protein arm had been run and no cross-protein performance figure had been examined
when any row of this table was written.** Every amendment so far revises a *covariate, a
convention, a hash, an arm or a threshold* — never an interpretation of a result, because there
are no results yet.
