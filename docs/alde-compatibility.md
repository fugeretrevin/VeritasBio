# ALDE compatibility on the cross-protein corpus

Resolves open question **O-3** of [`preregistration-cross-protein.md`](preregistration-cross-protein.md)
("ALDE compatibility on this corpus is unprobed") and supplies the HER2 ALDE arm that was
wrongly reported unavailable in run `run_20261031_051736_45ae3dd8`.

- **Headline:** ALDE as published runs **only on fully site-saturated libraries**. It does not
  fail loudly anywhere else — it completes, returns a full-length trajectory, and the
  trajectory is **a single variant selected once per query for the entire budget**.
- **Scope of the damage:** this affects **all six new landscapes *and* HER2**. Only TrpB and
  GB1 — the two 20-residues-per-site libraries — are unaffected.
- **The HER2 ALDE arm now exists and must not be quoted as a performance figure.** It scores
  0.292 mean top-1% hits against reliability's 3.125, loses 48–0 over 48 paired seeds, and
  **spends 26 of its 125 measurements**: it is the shared starting sample plus one variant
  measured 100 times. §5.
- The literature review's stated mechanism (long genotype representations) is **not** the
  cause. The cause is a divide-by-zero in ALDE's design-space normalisation. See §3.
- **Best estimate of what a repaired ALDE would have scored on HER2: 1.50** (48 seeds, a
  documented deviation) — between random and random forest, not a contender. §4a.

Everything below was produced with the authors' code at
`C:\Users\stoof\OneDrive\Documents\GitHub\ALDE`, commit `d0b3593dd179`, running in its own
interpreter (`C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe`,
python 3.11.9, torch 2.1.1+cpu, gpytorch 1.11, botorch 0.9.4, numpy 1.26.4). **No file in the
ALDE checkout was modified**; `git status` in that checkout is clean.

Reproduce with:

```powershell
# the HER2 re-run with ALDE included (48 seeds)
.venv\Scripts\python.exe -X utf8 benchmark_panoptes_isolated.py --library her2 `
  --methods random,ridge,random_forest,alde,panoptes_default,panoptes_qualify_then_optimize `
  --seeds 48 --base-seed 0 --initial-size 25 --rounds 4 --batch-size 25 `
  --provenance-mode acknowledge-unflagged `
  --alde-path   "C:\Users\stoof\OneDrive\Documents\GitHub\ALDE" `
  --alde-python "C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe"

# paired per-seed statistics over any run folder
.venv\Scripts\python.exe -X utf8 scripts\analysis\paired_stats.py <run_dir> --reference panoptes_default

# compatibility probes - MUST use ALDE's own interpreter, not .venv
C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe -X utf8 `
  scripts\analysis\alde_compat_probe.py encode --out encode.json
C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe -X utf8 `
  scripts\analysis\alde_compat_probe.py run --landscape cr9114_h1 --seeds 3 --out run.json
```

Path notes, because they have cost this project time before: `benchmark_comparison_results/` is
gitignored and does not exist in a git worktree, ALDE is **not** in this repository, and the main
`.venv` (python 3.14) has no torch — so a torch import failure there says nothing about ALDE's
availability. Every path above is absolute or repo-relative from the main checkout for that
reason.

### Cross-check on the probe itself

`alde_compat_probe.py` reconstructs each genotype independently of `veritas/datasets/` — from
the `mutant` column plus the parent in `scripts/data/library_wt_sequences.json` — so its numbers
are a second path to quantities `docs/dataset-inventory.md` computed by a different route. All of
them agree, which is what licenses the rest of this document:

| quantity | probe | `dataset-inventory.md` |
|---|---|---|
| eligible universe, all nine landscapes | 165,427 / 48,840 / 8,191 / 8,191 / 1,886 / 1,841 / 11,276 / 159,999 / 149,360 | identical |
| top-1% tier size | 1,655 / 489 / 82 / 82 / 19 / 19 / 113 / 1,600 / 1,494 | identical |
| sites and per-site alphabets | 15 / 16 / 13 / 13 / 11 / 11 / 10 / 4 / 4 | identical |
| parent residues | `GLRTAGQRDIDRVVT`, `FSSITANTKSTSSRYY`, `DVLTFNSAAYNNK`, `TSANQKSTSAL` | identical |

The probe also drops zero rows on all six new landscapes and recovers the published ceiling
exactly where a run reached it (CR9114-H1 9.83451, CR6261-H9 10.15000, eqFP611-red 1.69244).
The parent residue read from the WT sequence was additionally asserted against the wild-type
letter in every `mutant` token; no mismatch occurred on any landscape.

One departure from the inventory worth noting: the inventory lists HER2's per-site alphabets as
`20,20,19,19,18,13,12,10,9,5` (sorted descending), while the probe reports them **in position
order** as `9,20,19,18,13,12,5,19,10,20`. Same multiset, different presentation.

---

## 1. Per-landscape verdict

The protocol asks for RUNS / FAILS / RUNS-WITH-DEVIATION per landscape. The honest answer needs
a fourth category, because the interesting outcome is none of those three: the code **completes
without raising** and returns a trajectory that is not a search.

| landscape | variants | sites | alphabet/site | `generate_onehot` | encoding dim | constant one-hot columns | end-to-end | verdict |
|---|---:|---:|---|---|---:|---:|---|---|
| CREILOV_CHLRE_Chen_2023 | 165,427 | 15 | 2×12, 3×2, 5×1 | OK | 300 | **265 / 300 (88%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| CR9114_H1_Phillips_2021 | 48,840 | 16 | 2 (binary) | OK | 320 | **288 / 320 (90%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| EQFP611_BLUE_Poelwijk_2019 | 8,191 | 13 | 2 (binary) | OK | 260 | **234 / 260 (90%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| EQFP611_RED_Poelwijk_2019 | 8,191 | 13 | 2 (binary) | OK | 260 | **234 / 260 (90%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| CR6261_H1_Phillips_2021 | 1,886 | 11 | 2 (binary) | OK | 220 | **198 / 220 (90%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| CR6261_H9_Phillips_2021 | 1,841 | 11 | 2 (binary) | OK | 220 | **198 / 220 (90%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| *HER2/trastuzumab (already in the study)* | 11,276 | 10 | 9,20,19,18,13,12,5,19,10,20 | OK | 200 | **55 / 200 (28%)** | completes; 1 distinct pick per 100 queries | **FAILS (silently)** |
| *TrpB (control)* | 159,999 | 4 | 20 (uniform) | OK | 80 | **0** | completes; 100 distinct picks | **RUNS** |
| *GB1 (control)* | 149,360 | 4 | 20 (uniform) | OK | 80 | **0** | completes; 100 distinct picks | **RUNS** |

**There is no exception and no partial case.** Every library with at least one constant one-hot
column degenerates; every library with none works. The threshold is not "binary", not "many
sites", and not "long parent sequence" — it is **"is every one of the 20 amino acids observed at
every library position"**, which is true only of TrpB and GB1.

### What `generate_onehot` actually does with a binary alphabet

The literature review's concern was that ALDE's one-hot encoder assumes a saturated-site
alphabet and would reject a 2-residue site. **It does not reject it.**
`src/encoding_utils.py:46` builds its index from the module-level
`ALL_AAS` tuple of all 20 amino acids and allocates `len(seqs[0]) × 20` columns regardless of
what the library actually contains:

```python
one_hot_dict = {aa: i for i, aa in enumerate(ALL_AAS)}
onehot_array = torch.zeros([len(seqs), len(seqs[0]), len(ALL_AAS)])
```

So the call succeeds on all six landscapes, every row is a valid one-hot (verified:
`per_site.sum(dim=2) == 1` everywhere), and no character in any of the six genotype sets falls
outside `ALL_AAS`. The encoding is merely **wasteful**: on a binary library 18 of every 20
columns are structurally zero. That waste is what later becomes fatal.

### The dimensionalities are not the problem either

`[d, 30, 30, 1]` with `d` of 220–320 is a larger input layer than the `d = 80` of ALDE's
published TrpB/GB1 runs, but it is well within what the architecture handles — the guarded
diagnostic in §4 trains and searches normally at exactly these dimensionalities. Nothing failed
for reasons of size, memory or runtime. Per-seed wall time was 5–27 s on CPU for all nine
landscapes, CreiLOV's 165,427 × 300 design space included.

### Correction to the pre-registration's statement of the risk

§7.2 and open question **O-3** attribute the expected incompatibility to "the binary per-site
alphabets and the 121-/1,053-residue genotype representations". Two corrections:

1. **No genotype representation of these six landscapes ever reaches ALDE.** The benchmark hands
   its `alde` arm `view.combos`, which is one residue per library position — 11 to 16 characters.
   The full-length parent is used only to fill silent positions during ingestion. Full-length
   sequence length is therefore irrelevant to ALDE compatibility.
2. **The 1,053-residue figure does not belong to any of the six.** Per
   `scripts/data/library_wt_sequences.json` the parents are 119 (CreiLOV), 121 (all three
   antibody landscapes) and 233 (both eqFP611 landscapes) residues. 1,053 is the SpCas9
   representation from risk D-7, and SpCas9 is not in this corpus.
3. The binary-alphabet half of the statement is **right about the cause but wrong about the
   mechanism**: binary alphabets break ALDE, but through the normalisation step, not the encoder,
   and the same break hits HER2, which is not binary.

---

## 2. What the degenerate trajectory looks like

Configuration is exactly the benchmark's `alde` arm: `DNN_ENSEMBLE`, architecture
`[d, 30, 30, 1]`, Thompson sampling, `lrelu`, lr 1e-3, `train_iter` 300, dropout 0,
25 shared initial observations, `budget` 100, `batch_size` 25.

| landscape | seeds | status | distinct picks / 100 queries | duplicate picks | top-1% hits in the 25 starting variants | top-1% hits after 100 more queries |
|---|---:|---|---:|---:|---:|---:|
| creilov | 2 | completed | **1** | 99 | 1, 2 | 2, 3 |
| cr9114_h1 | 2 | completed | **1** | 99 | 1, 0 | 2, 1 |
| eqfp611_blue | 2 | completed | **1** | 99 | 0, 0 | 1, 1 |
| eqfp611_red | 2 | completed | **1** | 99 | 0, 0 | 0, 0 |
| cr6261_h1 | 3 | completed | **1** | 99 | 0, 0, 0 | 0, 0, 0 |
| cr6261_h9 | 2 | completed | **1** | 99 | 0, 0 | 0, 0 |
| her2 | 3 | completed | **1** | 99 | 0, 1, 1 | 0, 1, 1 |
| **trpb (control)** | 2 | completed | **100** | 0 | 0, 0 | **29, 28** |
| **gb1 (control)** | 2 | completed | **100** | 0 | 0, 0 | **52, 43** |

The picked index is **always literal index 0** of the design space, on every seed of every
affected landscape. The first ten picks of every affected run are `0,0,0,0,0,0,0,0,0,0`.

### The benchmark's own output already recorded this

This is visible in artefacts already on disk without any new instrumentation, which is worth
saying because it means the signal was available and was read as a performance result:

- `selection_trace.csv`, HER2, seed 0, rounds 1–4: the variant `A106T:M107H:D108V` appears
  **100 times**. One row is `accepted` (revealed fitness **0.2032**, barely off the 0.1 assay
  floor); the other 99 are `rejected / duplicate_within_batch`.
- `round_metrics.csv` reports `n_unique_measurements = 26` for every ALDE seed on HER2, against
  125 for every other arm. **The arm spends 99 of its 125 measurements re-measuring one
  variant.**

A reviewer reading only `summary.csv` sees `mean_measurements_used` and an ALDE top-1% figure
and has no reason to suspect either. The `n_unique_measurements` column is the tell, and it is
the column that should gate any future ALDE cell.

---

## 3. Root cause

Not the encoder. The design-space normalisation two steps later.

`src/optimize.py:128`, in `BayesianOptimization.__init__`:

```python
self.disc_X = botorch.utils.transforms.normalize(disc_X, self.domain)
```

where `self.domain` comes from `src/utils.py:39`:

```python
def domain_discrete(x):
    lower = x.min(dim=0)[0]
    upper = x.max(dim=0)[0]
    return lower, upper
```

and botorch's `normalize` is `(X - bounds[0]) / (bounds[1] - bounds[0])`.

A one-hot column for an amino acid that **never appears at that position anywhere in the design
space** is constant zero, so `upper == lower == 0` and the column normalises to `0 / 0 = NaN`.
From there, mechanically:

1. Every row of `disc_X` contains at least one NaN (verified: 100% of rows on all seven
   affected landscapes).
2. The first `Linear` layer of the DNN propagates NaN across all units, so **every** prediction
   over the design space is NaN — for every member of the 5-model ensemble, in every round.
3. `src/acquisition.py:51` masks already-queried points with
   `self.preds[samp_indices] = min(self.preds)`. `min` of an all-NaN array is NaN, so the mask
   writes NaN over NaN and excludes nothing.
4. `src/acquisition.py:53`, `ind = np.argmax(self.preds)`. NumPy's `argmax` on an all-NaN array
   returns **0**.

Every query of every round therefore returns index 0. Verified directly:

```
normalize(X, domain) ->  [[0., nan, nan], ...]        # constant columns become NaN
Linear(3,1)(N)       ->  [nan nan nan nan nan]        # NaN propagates
min(preds)           ->  nan                          # the "exclude queried" mask is a no-op
np.argmax(preds)     ->  0                            # the same index, forever
```

**Why the authors never hit this.** ALDE's README describes the package as being for
"combinatorial **site-saturation** libraries", and `generate_domain.py:27` builds its design
space as `generate_all_combos(nsites)` = `itertools.product(ALL_AAS, repeat=nsites)` — the
complete 20^n product. Under that construction every one-hot column varies by definition, so
`upper == lower` is unreachable and no guard is needed. Their published simulations use two
combinatorially complete 4-site libraries (TrpB and GB1). The assumption is sound for their
scope and silently violated by every library in this corpus except those same two.

### This is a correctness bug in ALDE, not a configuration mismatch

Worth stating plainly for the pre-registration: this is not "ALDE performs badly outside its
domain". ALDE **produces a NaN-valued surrogate and does not notice**. There is no warning, no
non-zero exit, no assertion. On a library it was not designed for it returns a well-formed
trajectory file of the correct length whose contents are one index repeated.

### The row-order accident that makes it dangerous

Because the degenerate pick is always index 0, the arm's reported `best_measured_fitness` is
determined by **the row order of the prepared CSV**, which is an arbitrary artefact of how the
corpus was exported:

| landscape | fitness of row 0 | landscape max | row 0's percentile |
|---|---:|---:|---:|
| **CR9114_H1** | **9.83451** | **9.83451** | **100.00 — row 0 IS the global optimum** |
| CREILOV | 4.13941 | 4.19552 | 99.60 |
| EQFP611_BLUE | 1.54427 | 1.60767 | 99.44 |
| EQFP611_RED | 0.07564 | 1.69244 | 14.44 |
| CR6261_H1 | 7.00000 | 9.72861 | 11.29 (at the assay floor) |
| CR6261_H9 | 7.00000 | 10.15000 | 20.15 (at the assay floor) |
| HER2 | 0.10000 | 15.08420 | 51.72 (at the assay floor) |

**On CR9114-H1 the broken arm finds the global maximum on every seed with probability 1**,
because the clean table is sorted by affinity descending. A `best_measured_fitness` endpoint
would have shown ALDE achieving a perfect result on CR9114-H1 while its surrogate was entirely
NaN. Had CR9114-H1 been run before this probe, that cell would have looked like the strongest
ALDE result in the study.

This is the concrete reason the protocol's "document the absent cell" rule matters more than it
looks: the failure does not announce itself, and on one of the six landscapes it masquerades as
a win.

---

## 4. Is the normalisation the *sole* cause? (diagnostic deviation)

To establish that nothing else about these landscapes troubles ALDE, the loop was re-run with
one change and nothing else: for columns where `upper == lower`, the reported upper bound is
widened from 0 to 1, so `normalize` divides by 1 and leaves the constant column at its constant
value. On a library with no constant columns the change is a no-op.

**This is a deviation from the published implementation and the results below are not ALDE as
published.** It was not applied to `benchmark_panoptes_isolated.py`, to the ALDE checkout, or to
any run folder. It exists to answer one question. Flag:
`alde_compat_probe.py run --guard-constant-columns`.

Three seeds per landscape, same configuration as §2, 100 queries after 25 shared starting
observations. "hits" is top-1% hits in the whole 125-measurement trajectory.

| landscape | top-1% tier size | distinct picks: published → guarded | hits, published code | hits, guarded |
|---|---:|---|---|---|
| her2 | 113 | 1 → **100** | 0, 1, 1 | **2, 5, 2** |
| creilov | 1,655 | 1 → **100** | 2, 3 | **45, 67, 41** |
| cr9114_h1 | 489 | 1 → **100** | 2, 1 | **10, 0, 0** |
| eqfp611_blue | 82 | 1 → **100** | 1, 1 | **5, 36, 37** |
| eqfp611_red | 82 | 1 → **100** | 0, 0 | **42, 0, 43** |
| cr6261_h1 | 19 | 1 → **100** | 0, 0, 0 | **13, 13, 13** |
| cr6261_h9 | 19 | 1 → **100** | 0, 0 | **11, 19, 11** |

Three things in that table are worth reading carefully:

- **CR6261-H9 seed 1 returned 19 of 19** — the entire top-1% tier. This independently confirms
  `dataset-inventory.md`'s warning that the CR6261 endpoint saturates at a 125-measurement budget
  and cannot discriminate between arms.
- **eqFP611-red (0 against 42 and 43) and CR9114-H1 (10, 0, 0) have enormous between-seed
  variance**, consistent with both landscapes having zero hits in their shallow shells and their
  optima sitting far from the depth-0 reference. These are landscapes where 48 seeds will be the
  minimum, not a luxury.
- **On CR9114-H1 the repaired arm scores *worse* on best-fitness than the broken one** (9.699
  against 9.83451). That is not a regression — it is §3's row-order accident, and it is the
  cleanest possible demonstration of why the broken arm must never be reported as a result.

With the single-line guard ALDE becomes a functioning searcher at every one of these
dimensionalities and alphabet sizes. **The only incompatibility is the `0/0`.** Encoding,
architecture, acquisition, ensemble training and the 220–320-dimensional inputs are all fine.

### 4a. How good is a repaired ALDE on HER2? Not as good as three seeds suggested

The three guarded HER2 seeds above (2, 5, 2) average 3.0, which is suspiciously close to
reliability's 3.125. Because that number would have changed the recommendation in §7, it was
re-run at **48 seeds** rather than quoted:

| arm | mean top-1% hits (incl. start) | sd | distinct measurements |
|---|---:|---:|---:|
| random | 1.188 | 0.867 | 125 |
| **ALDE, guarded deviation, 48 seeds** | **1.500** | **1.399** | **125** |
| random_forest | 1.688 | 1.728 | 125 |
| ridge | 1.958 | 1.663 | 125 |
| reliability (`panoptes_default`) | 3.125 | 0.981 | 125 |
| frontier (`panoptes_qualify_then_optimize`) | 3.062 | 1.040 | 121.5 |
| *ALDE, published code* | *0.292* | *0.544* | *26* |

**So the honest conclusion is the opposite of what three seeds implied.** A repaired ALDE on
HER2 scores **1.50**, which sits between random (1.19) and random forest (1.69), and is well
short of reliability's 3.125. Range over the 48 seeds was 0 to 5 hits; all 48 seeds made 100
distinct picks, so none of them degenerated.

Two caveats on that 1.50:

- It is **not paired** with the benchmark's arms. The probe draws its own 25 starting variants
  with a local RNG rather than the benchmark's shared starting sample, so this is a comparison of
  means, not a within-seed comparison. The sanity check is that the probe's mean hits in the
  starting sample is **0.25**, against the analytic expectation 25 × 113/11,276 = **0.2506** and
  the benchmark ALDE arm's 0.2917 — the starting distributions line up.
- The guard is the minimal fix, not a principled one (§8).

This matters for how the absent cell is read: the HER2 ALDE cell is absent **because of a bug**,
and the best available estimate of what it would have contained is a baseline-tier result, not a
contender. Both halves of that sentence need saying — the first stops the 0.292 being quoted as
ALDE's performance, the second stops the absence being treated as having hidden a winner.

---

## 5. HER2 with ALDE included — the re-run

Authoritative run folder: **`benchmark_comparison_results/run_20261031_132813_ff1a3423`**
(a first run, `run_20261031_130814_4ffecf18`, produced bit-identical numbers against an earlier
version of the benchmark file; see §5.4).
48 paired seeds (0–47), 25 initial + 4 × 25 = 125 measurements, `--library her2`,
`--provenance-mode acknowledge-unflagged`. Universe 11,276 variants; top-1% threshold
**2.4882**; 113 hits in the tier; achievable ceiling 113. ALDE resolved by absolute path to the
checkout at `C:\Users\stoof\OneDrive\Documents\GitHub\ALDE` (`d0b3593dd179`) and the interpreter
at `C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe`;
`probe_alde` reported `available: true`, no missing files, no missing dependencies.

### 5.1 Endpoint, with the reference run for comparison

| arm | mean top-1% hits (incl. start) | sd | reference run (48 seeds) | **distinct measurements** | mean best fitness |
|---|---:|---:|---:|---:|---:|
| random | 1.188 | 0.867 | 1.188 ✓ | 125 | 4.220 |
| ridge | 1.958 | 1.663 | 1.958 ✓ | 125 | 4.054 |
| random_forest | 1.688 | 1.728 | 1.688 ✓ | 125 | 3.613 |
| reliability (`panoptes_default`) | 3.125 | 0.981 | 3.125 ✓ | 125 | 4.473 |
| frontier (`panoptes_qualify_then_optimize`) | 3.062 | 1.040 | 3.062 ✓ | 121.5 | 4.444 |
| **alde (published code)** | **0.292** | **0.544** | *absent — wrongly skipped* | **26** | **2.059** |

The five pre-existing arms reproduce the reference figures **to four decimal places on all five**,
so adding ALDE perturbed nothing and the pairing is intact.

### 5.2 Paired statistics against reliability

Seed is the unit of independence; each difference is taken within a seed, over all 48 shared
seeds. Produced by `scripts/analysis/paired_stats.py`.

Metric: cumulative top-1% hits **including** the shared starting sample.

| arm vs `panoptes_default` | paired n | mean diff | 95% bootstrap CI | Wilcoxon p | sign-test p | W/L/T |
|---|---:|---:|---|---:|---:|---|
| **alde** | **48** | **−2.833** | **[−3.062, −2.604]** | **9.67e−10** | **7.11e−15** | **0 / 48 / 0** |
| frontier | 48 | −0.063 | [−0.188, +0.021] | 0.257 | 0.625 | 1 / 3 / 44 |
| random | 48 | −1.938 | [−2.229, −1.625] | 4.92e−09 | 2.61e−12 | 1 / 44 / 3 |
| random_forest | 48 | −1.438 | [−1.875, −0.958] | 4.18e−06 | 1.64e−06 | 6 / 37 / 5 |
| ridge | 48 | −1.167 | [−1.625, −0.688] | 6.49e−05 | 4.19e−05 | 8 / 35 / 5 |

ALDE loses to reliability on **48 of 48 seeds with zero ties**. It also loses to random on the
same basis.

### 5.3 Why those numbers are not a result, and must not be published as one

Metric: cumulative top-1% hits **excluding** the starting sample — what the arm itself found:

| arm | mean top-1% hits found by the arm |
|---|---:|
| reliability | 2.833 |
| ridge | 1.667 |
| random_forest | 1.396 |
| random | 0.896 |
| **alde** | **0.000 on every one of the 48 seeds** |

ALDE's endpoint on HER2 is **exactly its shared starting sample**. The per-round trace for every
seed is:

| round | n_unique_measurements | requested | accepted | rejected | status |
|---:|---:|---:|---:|---:|---|
| 0 | 25 | 25 | 25 | 0 | completed |
| 1 | 26 | 25 | 1 | 24 | short_batch |
| 2 | 26 | 25 | 0 | 25 | short_batch |
| 3 | 26 | 25 | 0 | 25 | short_batch |
| 4 | 26 | 25 | 0 | 25 | short_batch |

**0.292 is the hit count of 25 random variants** (expectation 25 × 113/11,276 = 0.251). The arm
measured one additional variant — `A106T:M107H:D108V`, fitness 0.2032, just off the 0.1 assay
floor — and then re-proposed that same variant 99 more times. Cause in §3.

**So the arm that was "wrongly skipped" would not have been informative even if it had been
included.** The reported skip in `run_20261031_051736_45ae3dd8` was wrong about the *reason*
(ALDE is on disk and fully available; it was a path-resolution failure, per the task brief) but
the cell it left empty is the cell this document says should be empty anyway — for a completely
different and much worse reason.

### 5.4 Process caveat: the repository moved under both runs, and it did not matter

`benchmark_panoptes_isolated.py` is being actively rebased by another agent, so neither run has
`repository_unchanged: true`. This was handled by running twice and checking what actually moved.

| | first run | authoritative re-run |
|---|---|---|
| folder | `run_20261031_130814_4ffecf18` | `run_20261031_132813_ff1a3423` |
| HEAD before → after | `a47664a` → `abbe430` | `abbe430` → `6aea8f0` |
| manifest schema | 1.0.0 | 2.0.0 (`run.state: PARTIAL`) |
| did `benchmark_panoptes_isolated.py` change during the run? | **yes** (4 commits, incl. "legible ALDE resolution") | **no** |
| did `workers/math_worker.py` or `veritas/datasets/` change? | yes (`math_worker.py`) | **no** |
| what changed during the run | the benchmark and selector | only a new, unrelated `veritas/calibration/` package and its tests |
| alde mean top-1% hits | 0.2917 | **0.2917** |
| alde distinct measurements | 26.00 | **26.00** |
| all five other arms | reference figures to 4 dp | **identical** |

The re-run's `repository_unchanged: false` is therefore a source-tree-hash flag, **not** evidence
that anything the run executed changed: `git diff abbe430 6aea8f0` touches only
`veritas/calibration/` and `tests/test_calibration_*`, neither of which the benchmark imports.
And the two runs agreeing to the last digit across two different versions of the benchmark file
is itself the strongest available check that the ALDE result is a property of ALDE and not of the
harness.

`run.state` is `PARTIAL` on the re-run, with `state_reason` "199 round(s) committed a short
batch". **192 of those 199 short rounds are ALDE's** (4 rounds × 48 seeds); the remainder belong
to the frontier arm, which also runs short on HER2 and did so in the reference run too. The
`PARTIAL` state is, in this run, very largely a report of the bug this document describes.

---

## 6. Budget: ALDE at 25 initial versus its own 96 + 4×96

This is the measurement §7.6 of the pre-registration pre-specified as its mitigation ("run the
arm grid at ALDE's own 96 + 4×96 schedule on TrpB and GB1") and recorded as not yet done. Here
is the ALDE arm of it, on the two landscapes where ALDE actually runs, 3 seeds each, published
code unmodified.

| landscape | schedule | total measurements | distinct picks | top-1% hits per seed | mean hits | **hits per measurement** |
|---|---|---:|---:|---|---:|---:|
| TrpB (tier 1,600) | 25 + 4×25 (ours) | 125 | 100 | 29, 28, 17 | 24.7 | **0.198** |
| TrpB | 96 + 4×96 (ALDE's) | 480 | 384 | 277, 234, 257 | 256.0 | **0.533** |
| GB1 (tier 1,494) | 25 + 4×25 (ours) | 125 | 100 | 52, 43, 29 | 41.3 | **0.331** |
| GB1 | 96 + 4×96 (ALDE's) | 480 | 384 | 259, 192, 199 | 216.7 | **0.451** |

**Direct answer to the fairness question: no, ALDE is not behaving pathologically at 25 initial
observations.** On TrpB and GB1 it searches normally — 100 distinct picks, no duplicates, no NaN,
and it finds 25–41 top-1% hits. What it does show is a **real efficiency penalty**: its
hits-per-measurement is **2.7× lower** on TrpB (0.198 vs 0.533) and **1.4× lower** on GB1
(0.331 vs 0.451) at our budget than at its own. Some of that is inherent — any learner improves
with training data, so a per-measurement rate rising with budget is expected and is not by
itself a handicap — but the magnitude on TrpB is large enough that it should not be waved away.

Three qualifications, so this is not over-read:

1. **Three seeds.** These are indicative figures, not the pre-registered mitigation run. The
   mitigation as specified is the *full arm grid* at 480, which still has not been run; this
   covers only the ALDE arm.
2. **Not comparable to the repo's published 30.42 TrpB figure.** These probe runs draw their own
   initial sample with a local RNG rather than the benchmark's shared starting variants, and
   count hits over the whole trajectory. They establish a within-probe 125-vs-480 contrast, not a
   restatement of the benchmark's number.
3. **The 480 schedule is only meaningful on TrpB and GB1 anyway.** On HER2 and on five of the six
   new landscapes, 480 measurements is 4.3% to 25% of the entire universe, and on CR6261 it
   exceeds it. The published schedule is not transferable to this corpus even where the code
   works.

So the existing TrpB and GB1 comparisons are **fair in the sense that matters** — ALDE is doing
real work there, not failing — while carrying a quantified budget handicap that §7.6 is right to
insist be reported alongside any ALDE figure.

---

## 7. Is ALDE usable as the reference baseline for the cross-protein study?

**No — not as published, on any landscape in this study except TrpB and GB1.** My honest read:

### 7.1 What must change in the pre-registration, before the study runs

1. **O-3 is resolved, and the answer is worse than "likely absent cells".** ALDE is an **absent
   cell on seven of the nine landscapes**: all six new ones and HER2. It is present only on
   TrpB and GB1. The rule in the **pre-registration's** §7.2 ("the protocol reports the cell as
   absent") applies, and should be applied now rather than discovered per-run.

   *All "§7.x" references below are to `preregistration-cross-protein.md`, not to this document.*
2. **The absence is not detectable by `probe_alde`.** That function checks checkout, interpreter,
   dependencies and imports, and on every one of these landscapes it correctly reports
   `available: true`. Availability is not compatibility. An arm can be fully available and
   return a null trajectory. **The pre-registration's "document the absent cell" discipline is
   only enforceable if something actually detects the failure**, and today nothing does.
3. **A detection rule exists and is cheap.** Either of these is sufficient and neither requires
   touching ALDE:
   - *pre-flight:* refuse the ALDE arm when the design space has any constant one-hot column
     (`(X.max(0) == X.min(0)).any()`), which is a one-line check over data already in memory;
   - *post-hoc:* refuse the result when `n_unique_measurements` is below the nominal budget — the
     benchmark **already computes this** and already wrote `26` into `round_metrics.csv`.
   I have not added either, because `benchmark_panoptes_isolated.py` is another agent's file.
   This is the single most valuable follow-up and should be a ticket.
4. **Correct the stated mechanism** in §7.2 and O-3 (see §1 of this document): drop the
   "121-/1,053-residue genotype representations" reasoning, which is not the cause and in the
   1,053 case not even this corpus.
5. **The HER2 ALDE figure from this re-run must be published as a null arm, with its
   `n_unique_measurements = 26`, not as ALDE's performance on HER2.** Reporting it as a
   performance number would be the most misleading row in the study.

### 7.2 On keeping ridge as the reference baseline

§7.6 already concludes "ALDE must not carry the claim; ridge is the reference baseline." That
conclusion is **correct and this work strengthens it**, but for a different reason than §7.6
gives. §7.6's reason is budget. The stronger reason is **coverage**: ALDE can be computed on
2 of 9 landscapes, so it cannot be the common comparator of a cross-protein study at all. A
reference baseline that exists on 22% of the corpus is not a reference baseline.

Ridge, random forest and random run on all nine. They should carry the comparison.

### 7.3 On whether to repair ALDE

§4 of this document shows a one-line guard makes ALDE work at every dimensionality and alphabet here. On HER2 at
48 seeds the repaired arm scores **1.50**, i.e. baseline-tier, not a contender (§4a). That cuts
the urgency but not the case:

- **Arguing for repair:** leaving it broken means the study's strongest baseline is absent from
  seven of nine cells. An absent strongest-baseline column is a real weakness in the evidence,
  and "their code divides by zero" is a thin reason to omit a method rather than fix an obvious
  bug in it. The guarded runs also produce genuinely informative numbers on landscapes where no
  ALDE figure otherwise exists at all (CreiLOV 41–67 hits of a 1,655 tier; CR6261-H9 reaching
  the full 19-hit ceiling).
- **Arguing against repair:** the moment the guard is applied, the arm is no longer "the
  published ALDE implementation run unmodified", which is the entire basis on which
  `AldeSelector` claims faithfulness, and the pre-registration explicitly prohibits "a different
  ALDE configuration".

**My recommendation:** report **both**, labelled unambiguously, and let neither carry the
headline.

- `alde` — published code, unmodified. **Absent** on seven of nine landscapes, with this
  document cited as the reason. Present on TrpB and GB1.
- `alde_constant_column_guard` — a **clearly-named separate arm**, not a configuration of the
  first, carrying the exact diff in its `config` block, run on all nine. It is the honest answer
  to "what would ALDE have scored here", and it is a strictly better artefact than an empty
  column.

This keeps the faithfulness claim intact where it is made, documents the incompatibility as the
protocol requires, and does not throw away the information. What must **not** happen is the
guard being folded into the `alde` arm, or the published-code null results being reported as
ALDE's performance.

### 7.4 On the fairness of the existing TrpB and GB1 comparisons

Those two landscapes are the only ones where ALDE was ever doing real work, and there it is
**not pathological at 25 initial observations** (§6 of this document): full distinct-pick trajectories, 25–41
top-1% hits. The 125-measurement budget costs it a factor of **2.7 in hits per measurement on
TrpB** and **1.4 on GB1** relative to its own 480-measurement schedule. So the published
TrpB (30.42) and GB1 (35.65) ALDE figures are real measurements of ALDE, obtained under a budget
handicap that §7.6 correctly insists be stated. They are not artefacts. The newly-measured
landscapes are the artefacts, and they are artefacts of a bug rather than of a budget.

---

## 8. What this document does not establish

- **Seed counts on the compatibility probes are 2–3**, except the 48-seed guarded HER2 run.
  They are sufficient for a binary compatibility verdict (the degeneracy is deterministic —
  index 0 on every seed) and **not** sufficient for any performance claim. The guarded hit counts
  in §4 and §6 are indicative only. §4a is a worked example of why: the 3-seed guarded HER2 mean
  was 3.0 and the 48-seed mean was 1.50, which would have inverted the recommendation in §7.3.
- **The guarded probe runs are not paired with the benchmark's arms** (own RNG for the starting
  sample), so §4a is a comparison of means. A paired figure needs the guarded arm run inside
  `benchmark_panoptes_isolated.py`.
- **The guard was not validated against the authors.** It is the minimal change that removes the
  `0/0`; it is not necessarily what the authors would do. A principled fix would more likely
  restrict the one-hot alphabet to the residues observed per position, which would also reduce
  `d` from 320 to 32 on CR9114 and is a larger change.
- **No 480-measurement full arm grid.** §6 covers the ALDE arm only.
- **ALDE's other encodings and models were not probed** (georgiev, DKL, GP, boosting ensemble).
  The failure is in `BayesianOptimization.__init__`, upstream of model choice, so it should hit
  every model; the georgiev encoding has no constant columns and would plausibly be unaffected.
  Neither was tested, and neither is what the benchmark's `alde` arm uses.
