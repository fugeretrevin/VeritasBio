# Handoff — state as of 2026-10-31

For a Claude session picking this up cold. Read this before touching anything. The
rules in **§2 are load-bearing**: each one is written down because violating it
already produced a wrong result in this project, in one case silently inverting a
published conclusion.

Suite: **674 passing, 2 skipped, 0 failures**. Everything described here is committed and
pushed to `origin/main`; run `git log --oneline -40` for the day's work rather than
trusting a commit hash quoted in prose, which goes stale the moment this file is edited.

> **Newer than this file:** [`HANDOFF-2026-10-31-reachability.md`](HANDOFF-2026-10-31-reachability.md)
> and its analysis [`day1-reachability.md`](day1-reachability.md). The candidate generator
> was hardcoded to one library geometry — four saturated sites, 20 residues each, depth
> capped at 4 by a schema validator — and was therefore blind on seven of the nine
> landscapes. Fixed and measured: the engine now beats one-hot ridge on CR9114 (5.00 vs
> 1.833 top-1%, 48 paired seeds, global optimum found) and on eqFP611-blue (9.50 vs 8.79),
> and scores **zero** on eqFP611-red against ridge's 15.35. TrpB is bit-identical, as the
> negative control requires. Three further rules are added there, §7 — the first being
> that a task-stop is not a process kill, which voided ten run directories.

---

## 1. Where things stand in one page

### The cross-protein sweep has now RUN, and the engine loses

**Read this before the TrpB table below.** 2026-11-01, commit `406ef58`, clean tree, six
landscapes, 48 paired seeds each, budgets matched per arm, seed as the unit of
independence, 10,000-sample whole-seed bootstrap. Δ is `panoptes_default − ridge` on
top-1% hits excluding the initial sample. Reproduce with
`scripts/analysis/generality_sweep.ps1` and read out with the recipe in §4.

| landscape | Δ vs ridge | t | seeds better/worse/tied | verdict |
| --- | --- | --- | --- | --- |
| CreiLOV | **−50.29** | **−32.01** | **0 / 48 / 0** | loss on every seed |
| GB1 | **−12.90** | **−5.63** | 6 / 41 / 1 | decisive loss |
| eqFP611-red | 0 hits vs ridge's 15.35 | — | — | loss; 7 short batches |
| eqFP611-blue | +0.73 | +0.77 | 25 / 18 / 5 | **not separable** |
| CR9114-H1 | **+3.17** | **+6.10** | 41 / 5 / 2 | win — but development data |
| HER2 | — | — | — | **gate FAILED, un-quotable** |

**One win on development data, one tie, three losses.** CR9114 is development data under
rule 2.4 — the intervention was designed against its measured failure — so it is
mechanism evidence, not generalization evidence. The day1 §7a pre-registration predicted
improvement on the reachability-limited landscapes; CreiLOV is reachability-limited at
13.3% and lost on **48 of 48 seeds**. That prediction is **falsified**.

This is close to what §8 of this document pre-registered ("that we lose… across-the-board
regression against ridge"), which is the one good thing about it: it was written down in
advance and so cannot be re-narrated now.

> **The veto hypothesis for these losses is FALSIFIED.** Measured 2026-11-01 at 48
> paired seeds with `--ablate-components spatial_filter`, branch
> `worktree-agent-a038e3b42e920227a` (`3804237`), pre-registered in
> `docs/prereg-spatial-veto-ablation.md` **before** any ablated arm ran. Ablating the
> veto **costs** hits on CreiLOV (−2.04, t=−2.33) and GB1 (−4.27, t=−3.22): on both it
> was *helping*, by raising the candidate pool's hit density 2.8× and 1.8×. On
> eqFP611-red it fully explains the **zero** — the ceiling is 0 of 82 and the arm scored
> exactly 0.0000 on all 48 seeds, which is arithmetic, not a measurement — but ablated it
> reaches only 4.21 against ridge's 15.35, so **73% of that loss is not the veto**. Every
> figure in the table above reproduced bit-for-bit, paired statistics included. **The
> CreiLOV and GB1 losses are therefore unexplained.** See
> `docs/spatial-veto-ablation-results.md`.
>
> Two cautions carried from that campaign: the veto helped CreiLOV *while reading a
> different protein's geometry* (`8QI9` maps 11 of 15 positions and all 11 disagree with
> the library's parent residue), so the benefit is **pool density, not biophysics** — a
> random filter enriching density 2.8× would have scored the same. And on GB1 the
> cohesion rule contributes **zero** rejections at any depth (the library fits in an
> ~8.2 Å sphere); all 118,766 GB1 rejections are single-site buried-core/steric calls
> firing at depth 1. **No landscape measured so far tests the cohesion veto on correct
> geometry.** TrpB is the only one that could, and it was not in scope.

> **THE CREILOV LOSS IS ROOT-CAUSED, AND IT IS A FABRICATED TRAINING ROW.** Measured
> 2026-11-01, branch `worktree-agent-a83565bcde59a81f3` (`bb925c9`), **12 paired seeds**
> with `configuration.seeding.seed_block_sha256` verified identical across arms and ridge
> held constant at 59.167. `select_next_micro_batch` injects one row per round for the
> wild type carrying the literal label `y = 0.0` whatever the assay's scale, at variance
> `1e-4` against a real measurement's `1e-2`. CreiLOV's real labels occupy
> **2.793–4.196**, so that row sits **11.6 units below the minimum — 8.3× the entire
> label range, 15.68 label SDs from centre** — and every round's surrogate is fitted
> through it.
>
> | arm | engine | Δ vs default | t | B/W/T |
> | --- | --- | --- | --- | --- |
> | default | 9.083 | — | — | — |
> | `prior` | 13.333 | +4.250 | 2.30 | 11/1/0 |
> | `physchem_features` | 9.250 | +0.167 | 0.19 (**p=0.74**) | 4/5/3 |
> | `prior,physchem_features` | 17.000 | +7.917 | 3.75 | 10/0/2 |
> | **`wt_anchor`** | **38.833** | **+29.750** | **12.17** | **12/0/0** |
> | `wt_anchor,prior` | 43.083 | +34.000 | 12.44 | 12/0/0 |
> | **`wt_anchor,prior,physchem_features`** | **50.167** | **+41.083** | **19.48** | **12/0/0** |
>
> −50.08 closes to **−9.00**, and best fitness reaches **4.1880 ± 0.0001 against ridge's
> 4.1851 ± 0.0043** — above it. **TrpB (0.35 label SDs) and GB1 (0.20) are the only
> landscapes genuinely centred near zero, i.e. the only two where this defect is
> invisible — and they are exactly the two every quotable engine figure here came from.**
> CR9114 violates the scale by 14.14 SDs and the engine *wins* there, so this is not a
> complete theory of the result set.
>
> **Not a default change, deliberately.** `wt_anchor` is an ablatable component; dropping
> or correcting the row moves TrpB and GB1 and would re-base every submitted number.
> **Next:** `--ablate-components wt_anchor` on CreiLOV at **48** paired seeds, then
> cross-landscape, where the scale-violation ordering predicts TrpB/GB1 unaffected as a
> built-in negative control.
>
> **The candidate-pool depth hypothesis is also FALSIFIED, on all nine landscapes.**
> `max_mutation_depth` does not cap the pool — `recombine_lead_variants` carries it past
> the declared depth from round 1. The pre-veto ballot spans depths **1–12** on CreiLOV
> and reaches **99.9%** of its hit set; the depth ceiling costs **1 hit of 1,655** there
> and **0** on the other eight. CreiLOV binds at **acquisition**: `ceil_post` is 70
> against ridge's 69, so the ballot sufficed, and stage-wise Spearman rho degrades
> monotonically (GP residual **+0.148** → composite µ(x) **−0.032** → post-acquisition
> **−0.104** → post-Pareto **−0.121**) with the prior at **−0.249, negative in 48 of 48
> campaigns**. The committed batch scores **0.81×** its own ballot's base rate — worse
> than drawing uniformly from the ballot it built. See
> `scripts/analysis/pool_hit_reachability.py` and `ranking_decomposition.py`.
>
> Two instrument defects found alongside, **not fixed**: `GenerationCollapseWarning`'s
> `depth_reach` divides by `max_library_depth`, and CreiLOV's depth 15 holds 10 variants
> and **zero hits**, so it demands proposals into an empty region, fires on every
> landscape including the winning `ridge` arm, and **would have corroborated the depth
> hypothesis that is now falsified**. And `generation_coverage.json`'s depth histograms
> are **post-veto** (`record_round` runs after the veto has replaced `novel_candidates`),
> so that artifact cannot distinguish a narrow ballot from a filtered one.

> **`threshold.<lib>` FAILS the pre-flight gate on the three runs this table is quoted
> from** — `run_20261101_024225_8fa1d86d`, `…_024544_e82e6214`, `…_025014_5e79e3ae` —
> because `qualification.quantile` and `qualification.source` are not recorded in the
> manifest. Nobody had seen it because `scripts/analysis/generality_sweep.ps1` invokes
> the gate **without** `--run-dir`, which makes both `threshold` and `budget` report
> NOT-APPLICABLE. The bar it concerns is consumed only by the
> `qualify_then_optimize`/`conditional_probe` arms and the failure is symmetric across
> arms, so it **cannot bias any paired contrast above**. It is nonetheless an un-passed
> gate check that the published figures carry, and it must be stated wherever they are.
> Run the gate with `--run-dir` from now on. Not fixed.

**Every previously reported figure reproduced exactly** — CR9114 +3.167 against the
reported 5.000−1.833, eqFP611-blue +0.729 against 9.500−8.792, red's ridge 15.354 to three
decimals. The measurements were sound. Two **verdicts** were wrong:

- **eqFP611-blue was called a WIN. It is not.** +0.73 at t=0.77, better on 25 of 48
  seeds — a coin flip. The win was asserted from a mean difference with no paired test,
  and it was one of only two generalization data points.
- **GB1 was called "parity at every budget". It is a decisive loss**, ridge 33.29 against
  the engine's 20.40, better on 6 of 48 seeds.

**Both controls are unsatisfactory, so the geometry diagnosis is NOT cleanly confirmed:**

- **GB1 (negative control) held, weakly.** 20.90 documented → 20.40 measured, so the
  generator fixes are a no-op on a saturated 4-site library as predicted. But only one
  GB1 48-seed run exists on disk, so the "before" is a number quoted in a document, not a
  re-analysable artifact. No per-seed data, no t, no interval.
- **HER2 (null control) MOVED**: 2.833 → 1.188, t=−8.29, paired on 48 seeds. HER2 is
  `pool_restricted`, so `generate_candidate_variants` is never called and nothing in the
  fix set can reach it. That is either a side effect or the dpp/dkl default change, and
  it **cannot be attributed from the artifacts** because both manifests wrongly claim dpp
  and dkl were enabled (fixed in `e083333`; see §2.11).
- **HER2's gate fails on INFORMATIVENESS, which no budget fixes.** There is also a real
  off-by-one — `generality_sweep.sh` asks for `25 + 4×22` = 113 against §3C's registered
  `budget = min(125, tier − 1)` = **112** — and it should be corrected. But correcting it
  does **not** open the gate. Verified by running the gate at 112: it still FAILS, because
  the d≥2 ceiling is 113 against a raw count of 113, i.e. **0.00% of the hit set removed
  against a 5% floor**. §3C separates the two checks and says so explicitly:
  informativeness is "a property of the library; no budget change affects it."

  An earlier version of this section claimed HER2 "was one measurement over the line" and
  that 112 would pass. That was wrong — it conflated the saturation check with the
  informativeness check. **§3A's "permanently untestable at this endpoint" stands.**

**A lead on HER2 that does NOT survive scrutiny.** Between the two runs k-DPP went
on → off and HER2's yield fell by 1.65 (t=−8.29), which is the direction §4's sign flip
predicts on a diffuse landscape. It is **not** usable, for three independent reasons, any
one of which is sufficient:

1. Confounded across a commit range, not a single-component change.
2. HER2's endpoint cannot discriminate two methods at all (above), so no arm run there is
   a registered test — and no budget choice changes that.
3. §3F: HER2's diffuse band is itself almost entirely **marginal skew**, not tail
   structure. Even a clean sign flip there would be confounded with library sparsity.

So do **not** spend a run on a controlled HER2 DPP arm expecting to close the
registration. It would be mechanism evidence on an un-gated landscape at best.

> **Careful:** `--ablate-components dpp` and `dkl` are currently **no-ops**, since both
> default off. An arm so labelled is identical to the default arm. To measure either, set
> the global True for the baseline arm. Pinned by
> `test_ablating_an_already_disabled_component_is_a_no_op`.

### What survives and can be quoted

On TrpB, held-out seeds 48–95, 25 initial + 4×25 = 125 measurements:

| result | figure |
| --- | --- |
| `panoptes_default` ("reliability"), raw top-1% hits | 45.21 (50.44 with `dkl,gflownet` ablated) |
| vs ridge, diversity-adjusted at d≥2 | **+7.93, t=+17.1** — larger than the raw t of +8.4 |
| at a 50-measurement budget | **6.3×** the best baseline, decaying to 1.7× at 125 |
| GB1 | ~~parity at every budget~~ — **WITHDRAWN**, measured as a loss; see above |

That arm was budget-matched at exactly 125, and TrpB is the one landscape where the
corrected qualification bar accepts the identical variant set, so none of the
2026-10-31 fixes touch it.

> **`panoptes_default` changed meaning on 2026-10-31.** `DPP_DIVERSITY_ENABLED` and
> `TWO_STAGE_GP_ENABLED` now default **False**, because both were measured to *hurt*
> top-1% yield. So the 45.21 above and any figure produced after that commit are
> **different arms under one name**. Always carry the flag state with the number; a
> bare "panoptes_default" is no longer a sufficient description of a configuration.
>
> **GB1's parity claim is withdrawn, and the re-measurement is done.** It was suspect
> because removing k-DPP is worth +6.25 top-1% hits on GB1 (t=3.92, 36 of 48 seeds,
> [`preregistration-motif-narrowness.md`](preregistration-motif-narrowness.md) §3B) —
> more than twice its effect on TrpB. The sweep settles it: on today's defaults the
> engine scores 20.40 against ridge's 33.29, better on 6 of 48 seeds. GB1 was never
> parity; it is a loss.

**Mechanistic findings, all intact:**

> **Read the first two scoped to TrpB.** They were measured on TrpB only, and on TrpB
> they hold. They do **not** generalise, and stating them unscoped inverts the diagnosis
> everywhere else: on seven of the nine landscapes the generator was *structurally
> blind* and generation is the confirmed bottleneck
> ([`HANDOFF-2026-10-31-reachability.md`](HANDOFF-2026-10-31-reachability.md)). Both
> pictures are true because TrpB is 4 saturated sites — depth 4 is full saturation, so
> the generator's radius-4 ceiling is invisible there and only acquisition can bind.

- **On TrpB**, the candidate generator is healthy — 2.80% landscape coverage at mean
  depth 3.28 in four campaigns, where ridge needs 48 campaigns to reach 2.52%.
- **On TrpB, acquisition is the bottleneck.** Round 1 was offered 193 depth-4 candidates
  across four campaigns and selected **zero** (9→0, 18→0, 139→0, 27→0). Consistent with
  greedy acquisition later measuring +11.27 top-1% (t=+12.20) against EI+Thompson.
- `_AnnotationDrivenPlanner` hard-codes `max_mutation_depth = 2` in round 1, and the
  GFlowNet sampler's `max_gfn_depth = max(strategy.max_mutation_depth, 3)` is the only
  thing that lifts it. That is why ablating `gflownet` shifts the engine shallower and
  *improved* TrpB yield.
- Depth bias is **proportional, not absolute**: reliability selects at 0.68–0.71× the
  universe's mean depth across library widths 4 and 10. `random` validates the metric
  at 1.002 and ridge sits at 0.92–1.03.
- **Parent functionality predicts enrichment on 9 of 9 landscapes.** Parent functional
  for the measured phenotype → shallow enrichment above 1×; parent not functional →
  near zero. This needs no labels, unlike shell enrichment itself.

### What is void — do not quote these

- **Every** `frontier` / `panoptes_qualify_then_optimize` number on every landscape
  except TrpB. The qualification bar was an absolute 0.70 that gated nothing.
- Every engine arm mean on GB1 and TrpB computed at 117–124 measurements against
  baselines' 125.
- All HER2 frontier results. The reported 3.062 against reliability's 3.125 is void.
- Every `alde` row on the seven landscapes where ALDE cannot function.

### What has never been done

**Zero independent protein validation.** HER2 is development data — the
pool-restricted code path exists *because* that library required it. The
cross-protein study has not run.

---

## 2. Rules you must not break

Each of these cost us a wrong result.

### 2.1 Genotype canonicalisation

Mutation strings **omit wild-type positions**. `V183A:F184A:V227A` is a *four-site*
genotype whose fourth site is WT. Comparing such strings directly lets `zip()`
truncate a 4-tuple to 3 and understate Hamming distances.

**This bug silently inverted a conclusion about hit diversity** — it made a real
advantage appear to vanish. Always expand to full library length filling WT at
omitted positions, and **assert equal length**. Never zip, never truncate.

Reuse `veritas/calibration/genotypes.py`. Do not write a sixth implementation.

### 2.2 Validate any genotype code before trusting its output

Before any new number, reproduce these two:

```
TrpB shallow shell (1 <= depth <= 2): 168 / 2,242 = 7.4933%
GB1  shallow shell (1 <= depth <= 2):  72 / 2,167 = 3.3226%
```

Convention: **depth 0 is excluded from both shells but retained in the universe.** A
wild type is not a search target. If you cannot reproduce both figures, your parser
is wrong and nothing downstream is trustworthy.

### 2.3 Thresholds come from manifests

Read `top1.threshold` from each run's own `manifest.json`. Do not recompute it. The
eligible universe depends on library-specific eligibility rules — HER2 excludes
25,114 unlabelled rows whose `DMS_score` is exactly 0.0, which are absent labels and
not measurements of zero.

### 2.4 TrpB, GB1 and HER2 are all development data

Do not tune anything against them and do not present a result on them as
generalization evidence. If you inspect a held-out target and change the algorithm
because of what you saw, that target is reclassified as development data and a new
one must be reserved. This rule has already fired once, on HER2.

### 2.5 The pin must not change; the lock may only grow

`scripts/data/cross_protein_corpus.lock.json` pins CombinGym at commit
`c5dccd514a9acfe1089a5f1555faea98afd78019`. The commit, date and tree SHA are a
scientific commitment and must not move. Adding entries is fine and sometimes
required — a lock that refuses to cover a file someone is about to read protects
nothing. `scripts/freeze_cross_protein_corpus.py verify` must report **29/29**.

### 2.6 Availability is not compatibility

An external arm can be installed, importable, and still unable to function on a
given library. ALDE reports `available: true` on seven landscapes where it silently
returns index 0 forever. The manifest keeps `availability_probe` and `compatibility`
as separate keys. Never collapse them.

### 2.7 Run the pre-flight gate before any scored run

```
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library <key> [--run-dir <dir>]
```

Six checks, each written against a defect that actually happened. Non-zero exit means
do not trust the run. Two checks need an actual campaign and report NOT-APPLICABLE
with the missing input named rather than passing quietly.

### 2.8 Path resolution

A git worktree **cannot see** gitignored paths or anything outside the repository.
This cost three full agent runs in one night. `benchmark_comparison_results/` is
gitignored (`.gitignore:246`) and does not exist in a worktree. ALDE is **not in this
repo** — it is at `GitHub\ALDE` with its own interpreter at `GitHub\.venv-alde`.

Use `veritas/paths.py`. A failure to resolve must **name every path it checked** — a
bare "not found" is what let our strongest baseline get silently dropped from a
48-seed run. See `docs/agents/paths.md`.

### 2.9 Never let outcomes reach selection

Anything reading true fitness is retrospective only and gated behind an explicit
access token in `veritas/diagnostics/retrospective.py`. Configuration fingerprints
must identify a run **without consulting outcomes** — there is a test asserting a
32-arm sweep yields 32 distinct fingerprints with no outcome leakage.

### 2.10 A check validated only against its own fixture proves nothing

The pre-flight threshold check's first version tested `p >= 1 - 1e-12`. On real HER2
the figure is 0.908, so **it passed the exact defect it was written for.** Always
validate a new check against the real artifact, not the fixture you built beside it.

---

### 2.11 A verdict needs a paired test, not a mean difference

`eqFP611-blue: 9.50 vs 8.79` was published as a **WIN**. The paired test on the same
numbers gives t=0.77, better on 25 of 48 seeds — a coin flip. It was one of only two
generalization data points this project had, and it was not a result at all.

Before any arm is called better or worse: compute the paired statistic over seeds, report
seeds-better / worse / tied, and state the interval. A difference of means with no test
is not a verdict. Two numbers side by side invite one anyway, so write the test next to
them.

The same rule killed the redundancy argument in `docs/545.md` §4 and inverted two of five
ablation signs. It is the most repeated error in this project.

### 2.12 The manifest's component states are only as honest as what they read

`component_states()` derived state from `--ablate-components` alone, so `dpp` and `dkl` —
gated by module globals — were recorded as **enabled** in every run where they were in
fact off. The `configuration_fingerprint`, whose documented job is telling two
configurations apart, was blind to the change: a dpp-on and a dpp-off run of one library
hashed identically.

That cost a real attribution. The HER2 null control moved by t=−8.29 and the cause could
not be read from the artifacts, because both manifests claimed the same component state.
Fixed in `e083333`, which **changes the fingerprint of every future run** — so a
fingerprint comparison across that commit is meaningless. Quote the commit, not the
fingerprint.

The general form: when a behaviour is gated in two places, a record that reads one of
them is worse than no record, because it will be trusted.

## 3. Known traps

- **Our prepared CSVs are sorted by fitness descending.** Combined with ALDE's
  index-0 bug, **CR9114's row 0 is the global optimum** — a wholly broken arm would
  have posted a perfect best-fitness. Any method must be invariant to row order;
  `veritas/preflight.py` has the shuffle check.
- **`build_synthetic_universe` is 4-site only**, so `--smoke-test --library her2`
  fails on the 10-site geometry. Pre-existing.
- **The frontier policy's Hamming-1 shell is exhaustible.** On TrpB it holds ≤76
  untested neighbours, so the arm cannot fill 125 with `rejected = 0`. On HER2 the
  median library member has **zero** in-library Hamming-1 neighbours. This is an open
  design problem, not a bug to patch — see §5.
- **A caption is not a citation.** Two of our literature claims were wrong and both
  were caught by reading the figure's axes and scoping rather than its caption. One
  cited figure, read literally, said the opposite of what we claimed.
- **`--qualify-threshold` 0.70 is below the global minimum** on four landscapes
  (bnAb assays are −log₁₀(K_D) floored at 7.0). Superseded by `--qualify-quantile`,
  default 99.96.
- **Reading eLife HTML through a summarizing fetch tool is refused by a content
  filter** on some articles. Download the PDF and extract locally; asset lists are at
  `api.elifesciences.org/articles/<id>`.

---

## 4. Where things live

| what | where |
| --- | --- |
| Benchmark harness | `benchmark_panoptes_isolated.py` (~4,800 lines) |
| Selection engine | `workers/math_worker.py::select_next_micro_batch` |
| Dataset schema + ingestion | `veritas/datasets/` |
| Coverage, trace, depth returns | `veritas/diagnostics/` |
| Threshold / endpoint calibration | `veritas/calibration/` |
| Pre-flight gate | `veritas/preflight.py` |
| Path resolution | `veritas/paths.py` |
| Analysis scripts + their outputs | `scripts/analysis/` (read its README first) |
| Frozen protocol | `docs/preregistration-cross-protein.md` |
| Corpus characterization | `docs/dataset-inventory.md` |
| Audit: clean vs raw tables | `docs/corpus-audit.md` |
| Label provenance | `docs/label-provenance.md` |
| ALDE incompatibility | `docs/alde-compatibility.md` |
| Literature anticipation of H1 | `docs/ordering-anticipation.md` |
| CR9114 re-preparation | `docs/cr9114-reprep.md` |
| Prior art on depth allocation | `docs/prior-art-depth-allocation.md` |

Run artifacts live in the **gitignored** `benchmark_comparison_results/`. Resolve by
absolute path from the real repo root.

---

## 5. The open problems — now answered

All six were resolved by external review on 2026-10-31. These are decisions, not
suggestions; implement to them.

### 5.1 Candidate source → `frontier_pool_v2`

None of the three options we proposed. **Locality means proximity to promising
*measured* sequences within the permitted candidate pool — not exactly one
substitution.** Build an adaptive neighbourhood: multiple measured leads as centres,
include unmeasured candidates in increasing distance order, expand the radius until a
shortlist capacity is reached (~250 for a batch of 25 — an engineering starting point,
not an optimum), deduplicate across leads, deterministic tie-breaking.

Reserve a **fixed, documented** fraction of capacity for candidates outside the local
neighbourhoods. This is part of the policy, *not* an undocumented fallback that fires
after failure.

**Three concepts we had conflated must be separate in code and docs:**

| concept | meaning |
| --- | --- |
| `design_space` | sequences the library specification permits |
| `assay_available_pool` | sequences selectable in the experiment |
| `oracle_evaluable_pool` | sequences with recorded outcomes, for retrospective replay |

They coincide in a fully characterised library and **must not be treated as
conceptually identical**. Restricting selection to historically measured variants is
acceptable for retrospective benchmarking only if disclosed as an oracle-coverage
constraint and applied identically to every method. Missing historical labels do not
mean a sequence is biologically impossible.

Log the effective search radius. **If most selections come from large distances or the
broader-search component, the result must not be presented as demonstrating the
benefit of local optimisation.** Keep `frontier_legacy_v1` alongside `v2`; never
silently modify the old one and reuse its results.

### 5.2 Qualification → remove landscape quantiles from online decisions entirely

**Our quantile fix is still an oracle leak.** This is a sharper diagnosis than we had.
A prospective lab cannot know the 99.96th percentile before measuring most of its
library, and no finite-library adjustment or minimum-count floor fixes that
information problem. The quantile is fine as a *retrospective evaluation* threshold
and unacceptable as an *online decision* threshold.

**First implementation — a fixed, oracle-free phase schedule.** Under the
125-measurement protocol: 25 shared random, then 25 reliability-oriented, then 75
frontier-oriented. Switch at 50 total measurements. It is independent of the unknown
fitness distribution, works across assay scales, cannot be satisfied by the initial
sample, and needs no global-optimum estimate.

It is **deliberately less sophisticated** than what it replaces, and that is the
point: it removes one source of dataset-specific behaviour so the frontier
acquisition strategy itself can be evaluated.

Rename it: this is a **`phase_transition`**, not a fitness qualification — reserve
"qualification" for an actual user-defined experimental target.

Only after benchmarking the fixed schedule, investigate observed-progress switching
(switch when a measurement improves on the best initial measurement by a predefined
margin, with a forced checkpoint so it cannot stay in phase one forever). Treat it as
a separate arm, not an assumed improvement.

Note this makes `v2` a genuinely new arm **on every landscape including TrpB**, since
the old switch was data-dependent and the new one is positional.

### 5.3 Stopping → out of the primary study

Keep fixed-budget benchmarking as the primary protocol. Changing candidate generation,
qualification *and* stopping at once makes attribution impossible. There is no
universally valid oracle-free rule that proves an optimum has been found from a small
sample.

Two specific warnings. **Do not stop because best fitness has not improved for two
rounds** — a flat curve can mean local exploitation failure rather than convergence.
And **budget saturation is a metric problem, not a budget problem**: that 125
measurements exceeds CR6261's 19-variant hit tier makes raw top-1% count a poor
endpoint there, to be fixed with measurement-efficiency curves and normalisation
against the eligible hit-set size. Do not use early stopping to paper over a badly
defined metric.

### 5.4 Depth covariate → keep the primary, add a secondary

Preserve `depth <= floor(S/2)` as the primary H1 covariate. Add a
**distribution-relative secondary** — enrichment below the *median eligible-library*
mutation depth — plus a continuous depth analysis rather than relying on a binary
split, which can distinguish a smooth relationship from a threshold artifact.

Amending is acceptable because no outcomes have been inspected, provided the original
is archived, the amendment is dated before analysis, and results are reported
whichever covariate is more favourable. Do **not** describe the secondary as blind —
published figures already reveal depth-related patterns.

Limitation to document: a depth distribution computed from historically measured
variants is a property of the **sampled library**, not of feasible sequence space.

### 5.5 Depth bias → frozen-state counterfactuals

**A correction to our own hypothesis, and it matters.** We suggested that posterior
variance growing with mutational distance would make an acquisition function prefer
shallow candidates. **That has the sign backwards.** If the acquisition function
rewards uncertainty with a positive coefficient, higher variance at depth would
*favour* deep candidates, not suppress them. Read the actual acquisition equation
before hypothesising.

Method: freeze one decision state (measured data and outcomes, fitted surrogate,
candidate pool, prediction vector, structural-prior components, acquisition inputs,
random state), then evaluate alternate configurations **on that same state**. Trace
the candidate depth distribution at every stage — generation, prediction, acquisition
ranking, Pareto, DPP, final batch — and distinguish filtering from score reordering
from diversity selection.

Four counterfactuals on the frozen state: prior neutralised; uncertainty term varied
with predictions held fixed; Pareto off; DPP off. Then **one small paired
full-campaign run** to confirm, because frozen-state interventions cannot show
long-term effects — changed selections alter future training data.

Primary evidence is controlled ranking differences under one frozen state.
Correlations are secondary.

**And a reframe:** lower selected depth is not intrinsically bad. If TrpB's shallow
region is genuinely hit-enriched, shallow preference is *rational* there. The goal is
not to increase depth; it is to determine when the engine's depth preference aligns
with the scientific objective.

### 5.6 H2 → narrow claims only

What H2 can support if positive: *"the optimization advantage changes predictably
between two fitness landscapes defined over the same combinatorial genotype space."*

What it cannot: general superiority across families, a universal depth relationship,
causality without an intervention isolating depth preference, prospective wet-lab
effectiveness, or reliable adaptation to new targets. **48 campaigns on one biological
system are not 48 independent replicates.**

**Our natural-experiment framing was over-read.** Better described as a *paired,
within-library comparison across biological readouts* than a causal natural
experiment — the two assays may differ in measurement noise, dynamic range and
censoring, not only in depth-related geometry. Those alternatives must be evaluated.

Three proteins across two new families supports a **preliminary cross-protein
characterisation**, not a generalization claim. Report per-protein prominently; any
aggregate is exploratory.

### 5.7 Genotype architecture → consolidate after the freeze

Separate the production implementation from independent test oracles. **Do not have
every test call the canonical implementation and compare it against itself.** Retain
the TrpB and GB1 anchors, the CR9114 reweighting cases, and the WT-omission example as
golden fixtures, plus property-based tests (round-trip conversion, equivalent mutation
strings, Hamming symmetry, strict equal-length distance, per-site alphabets, invalid
residues).

Timing: **not during the frozen study** unless a correctness defect demands it. Do it
in a branch, run old and new against the full reference corpus, verify identical
benchmark inputs and eligibility decisions, then remove duplication. Resolve the two
`LibrarySpec` classes — do not keep two incompatible classes with the same name
indefinitely.

### 5.8 Arm validity statuses — implement this immediately

Every benchmark arm carries an explicit status, and **only `VALID_COMPLETE` arms at
equal budget enter the primary matched-budget comparison**:

`VALID_COMPLETE` · `VALID_BUDGET_LIMITED` · `INVALID_NUMERICAL` ·
`INVALID_PROTOCOL` · `UNAVAILABLE` · `PROVENANCE_UNRESOLVED`

This is how the void results from 2026-10-31 get represented rather than deleted. A
legacy policy that cannot fill its budget is `VALID_BUDGET_LIMITED` and must not have
its mean compared directly against fully budgeted methods. ALDE on the seven
incompatible landscapes is `INVALID_NUMERICAL`. Any numerical fix to ALDE becomes a
separately named `ALDE_stabilized` that must be validated against the authors'
behaviour on compatible datasets — a stabilised implementation is not automatically
the published method.

Full reasoning: the expert review of 2026-10-31, and the "Open Problems for Expert
Review" doc it responds to.

---

## 6. Decisions — all three that were open are now closed

Recorded because an earlier version of this file listed them as blocking, and a
session reading a stale copy would stop for no reason:

- **The CR6261 qualification floor — dissolved, not decided.** The question was
  whether to floor the qualifying set at a minimum count or declare the stage
  inapplicable below some library size. §5.2 removes the global-quantile bar from
  online decisions entirely, so there is no floor to set. Do not re-raise it.
- **The frontier policy's exhausted shell — answered.** §5.1, `frontier_pool_v2`.
- **Whether to redefine the covariate — answered.** §5.4: keep the pre-registered
  primary, add a distribution-relative secondary and a continuous depth analysis.

### What genuinely still needs a human

- **Hackathon framing and timing.** The demo narrative is already fully supported by
  what exists and needs no new runs; see §1. Given that the pre-registration predicts
  across-the-board regression (§8), a result is unlikely to improve the story.
- **Whether to publish or raise on a within-protein result.** The external review
  declined to promise publication from a paired result alone and named what a strong
  outcome needs instead: the H2 contrast, frozen-state evidence identifying the stage,
  a full-campaign intervention demonstrating it, **and** confirmation on at least one
  further independent family. That is a strategy call, not a technical one.
- **Whether to report the ALDE constant-column bug upstream.** It is a genuine
  correctness defect in published code that makes the method silently return index 0
  on any library with an unobserved residue at a position. We have not contacted the
  authors and would not without being asked.

---

## 7. The order of work, as revised by the expert

Candidate-space correctness comes **before** any performance comparison. Running the
development grid first would spend compute characterising a frontier arm that
structurally cannot fill its budget, then compare it to baselines that fill theirs.

**Phase 0 — freeze the specifications before any large run.** Define
`frontier_pool_v2`, the fixed phase schedule, the `design_space` /
`assay_available_pool` / `oracle_evaluable_pool` semantics, primary and secondary
depth covariates, and the arm validity statuses. Preserve the legacy policy IDs.
Commit all of it before viewing any outcome.

**Phase 1 — candidate-space correctness.** Adaptive feasible-pool neighbourhoods,
multi-lead generation, the documented broader-search component, deduplication,
budget-fill validation, candidate-source provenance, policy versioning. Required
tests: saturated 4-site; sparse 10-site HER2-like; a **disconnected Hamming-1 graph**;
mixed per-position alphabets; a small library with few unmeasured candidates left;
duplicate lead neighbourhoods; exhaustion; deterministic tie-breaking; exact budget
fulfilment when feasible; no hidden-outcome access. *Acceptance is correctness, not
performance — better fitness is a separate hypothesis.*

**Phase 2 — replace qualification** with the fixed phase schedule and
observed-data-only transition logic, with transition provenance recorded. Keep the
oracle-quantile behaviour only as a clearly named retrospective legacy analysis.

**Phase 3 — mechanistic diagnostics.** Frozen-state interventions for prior,
uncertainty term, Pareto and DPP; stage-by-stage depth tracing; then one small paired
full-campaign confirmation.

**Phase 4 — the independent benchmark.** Validated ridge, random forest, random,
other compatible baselines, VeritasBio reliability, legacy frontier where applicable,
and the new adaptive frontier. Full budget accounting, best-so-far, top-1% yield,
diversity-adjusted yield, tail discovery where resolvable, depth distribution,
candidate-source behaviour, paired intervals. **Do not tune the new policy on the
final independent landscapes.**

**Phase 5 — stopping, as a separate pre-registered study.** Not in the primary.

**Phase 6 — consolidate the genotype infrastructure**, after the scientific freeze.

**Explicitly out of scope for now:** a new GP architecture, a foundation model, a new
interaction-probing algorithm, or an automated policy-learning controller. The next
contribution comes from understanding and correcting the existing decision process,
not from adding complexity.

Two smaller items that are not blocked:

- **Publish digests for the eight prepared CSVs.** The lock hashes only upstream
  CombinGym files, never the tables the arms actually read. Only CR9114 has a digest.
- **Acknowledge the eight non-CR9114 derived tables** in the gate via
  `--acknowledge-derived KEY=WHY`, citing `docs/corpus-audit.md`. The audit proved
  clean ≡ raw for all of them, so the check is now over-strict by evidence — but it
  should be acknowledged explicitly, not silenced.

---

## 8. What the pre-registration predicts

**That we lose.** Every independent landscape sits below GB1's shallow-shell
enrichment, where the engine already lost by 22.49 hits. The protocol predicts
across-the-board regression against ridge, and that prediction is deliberately
**not** softened — with the ordering demoted to a published covariate, it is the main
thing H1 still tests.

The most likely overall outcome is pre-registered as **mechanism supported,
generalization undetermined**. That is recorded in advance as a legitimate reportable
result, specifically so it cannot be re-narrated as failure later.
