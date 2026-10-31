# Handoff — state as of 2026-10-31

For a Claude session picking this up cold. Read this before touching anything. The
rules in **§2 are load-bearing**: each one is written down because violating it
already produced a wrong result in this project, in one case silently inverting a
published conclusion.

Main is at `6ba8d1d`, pushed. Suite: **664 passing, 0 failures**.

---

## 1. Where things stand in one page

### What survives and can be quoted

On TrpB, held-out seeds 48–95, 25 initial + 4×25 = 125 measurements:

| result | figure |
| --- | --- |
| `panoptes_default` ("reliability"), raw top-1% hits | 45.21 (50.44 with `dkl,gflownet` ablated) |
| vs ridge, diversity-adjusted at d≥2 | **+7.93, t=+17.1** — larger than the raw t of +8.4 |
| at a 50-measurement budget | **6.3×** the best baseline, decaying to 1.7× at 125 |
| GB1 | parity at every budget |

That arm was budget-matched at exactly 125, and TrpB is the one landscape where the
corrected qualification bar accepts the identical variant set, so none of the
2026-10-31 fixes touch it.

**Mechanistic findings, all intact:**

- The candidate generator is healthy — 2.80% landscape coverage at mean depth 3.28 in
  four campaigns, where ridge needs 48 campaigns to reach 2.52%.
- **Acquisition is the bottleneck.** Round 1 was offered 193 depth-4 candidates across
  four campaigns and selected **zero** (9→0, 18→0, 139→0, 27→0).
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

## 5. Open problems — do not start work without an answer

These went to an external expert on 2026-10-31. **Do not implement a fix to 1 or 2
without the answer**; both were deliberately left open because the obvious fixes risk
either changing what the algorithm is or leaking oracle information.

1. **Candidate-source exhaustion.** Is a locality-constrained acquisition policy
   *supposed* to degrade to global ranking when its locality is exhausted, or is
   exhaustion the signal that the policy does not apply to that library? Determines
   whether `frontier` appears in the cross-protein study at all. The agreed shape of
   any fix is a **new named policy version**, so existing results stay interpretable
   as results of the old policy.
2. **Oracle-free qualification and stopping.** A quantile transfers the *rarity* of a
   bar but not its *attainability* — on CR6261 the top 0.04% is one variant. Every
   floor we can construct needs tier information we should not have. And there is **no
   stopping rule at all**, on targets where the 125 budget already exceeds the entire
   hit set (HER2 113, eqFP611 82, CR6261 19).
3. **Covariate and mechanism are on different scales.** The mechanism is relative to
   the universe's mean depth; the covariate cut is relative to library width.
4. **Which depth-discount mechanism to test first** — uncertainty calibration, the
   structural prior, k-DPP, or Pareto. All four correlate with depth, so four
   correlations will all come back significant and prove nothing.
5. **What H2 alone can claim.** H1's ordering is substantially pre-published; H2 is
   one protein.

Full context: the "Open Problems for Expert Review" doc.

---

## 6. Decisions only the user can make

- **The CR6261 qualification floor.** Blocks any CR6261 frontier arm. Either floor
  the qualifying set at a declared minimum count, or declare the stage inapplicable
  below a stated library size. Must be fixed *before* the arm runs.
- **The frontier policy's exhausted shell** (open problem 1) — a policy change, not a
  correctness fix.
- **Whether to redefine the covariate** (open problem 3), which means amending a
  frozen pre-registration before any results exist.

---

## 7. The order of work, as revised by the expert

Candidate-space correctness comes **before** any performance comparison. Running the
development grid first would spend compute characterising a frontier arm that
structurally cannot fill its budget, then compare it to baselines that fill theirs.

1. Fix candidate-source exhaustion via a versioned policy *(blocked on answer)*
2. Define qualification and stopping rules that work at both library scales without
   hidden-oracle information *(blocked on answer)*
3. Test both mechanisms on constrained CombinGym libraries
4. Only then re-run the development grid and compare performance
5. Then H2 on eqFP611 — it is exempt from the anchor-run block and reportable as soon
   as eqFP611 runs
6. Then the full cross-protein evaluation

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
