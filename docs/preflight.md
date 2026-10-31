# The pre-flight gate

Five defects in the last day each invalidated or nearly invalidated a result. Every one of
them was **already recorded in an artefact we wrote ourselves**, and nothing asserted on
it:

| the defect | the artefact that already held the number |
|---|---|
| ALDE ran as a null arm on 7 of 9 landscapes, picking the same variant 100 times | `round_metrics.csv`: `n_unique_measurements = 26` against a budget of 125 |
| our own arms measured 120–124 variants where baselines measured exactly 125 | `round_metrics.csv` |
| CR9114 was missing 25% of its variants, with depth-correlated attrition | the depth profile in `docs/dataset-inventory.md` |
| the primary endpoint was degenerate on HER2: the d>=2 ceiling equalled the raw ceiling at 113 | both numbers published |
| the qualification threshold gated nothing on 8 of 9 landscapes | every HER2 seed logged "qualified before round 1" |

So this is not a detection problem, it is an **assertion** problem.
[`veritas/preflight.py`](../veritas/preflight.py) is the assertions.

Nothing in it is a bare boolean. Every check returns a `CheckResult` with a status of
`PASS` / `FAIL` / `NOT-APPLICABLE` / `SKIPPED`, the measured value, the expectation and —
for anything that is not a pass — a reason. `bool(result)` raises `TypeError`, because
`if check(...):` quietly reading "not applicable" as "fine" is the failure mode the gate
exists to prevent. A skipped check is a recorded result, never an absence.

## How to run it

```powershell
# everything that can be answered from the repository alone
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library her2 --arms alde,ridge

# with a run folder, which adds the budget accounting and the qualification bar
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library gb1 --arms alde,panoptes_default `
  --run-dir benchmark_comparison_results\run_20261031_132813_ff1a3423

# one check at a time; a skip is recorded in the report
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library trpb --arms alde --only encoding
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library trpb --arms alde --skip corpus,row_order

# machine-readable
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library her2 --arms alde --json
```

Exit status is **0** when nothing applicable failed and **1** when any check failed.
`NOT-APPLICABLE` and `SKIPPED` do not fail the gate, and both always state why.

In Python:

```python
from veritas import preflight

report = preflight.run_preflight(library="her2", arms=("alde",), run_dir=run)
report.exit_code          # 1 if anything failed
report.failures           # the CheckResults that failed
report.to_dict()          # JSON-safe, for a run manifest
```

Each check is also individually runnable: `check_corpus_integrity`,
`check_encoding_compatibility`, `check_budget`, `check_endpoint_viability`,
`check_threshold_transfer`, `check_row_order_invariance`.

## The six checks

### 1. `corpus` — corpus integrity

Every file a landscape reads is present in `scripts/data/cross_protein_corpus.lock.json`
and both hashes match, and a derived table is flagged when a raw one exists at the same
pin.

* **Motivating defect:** CombinGym's clean CR9114 workbook silently excluded every
  genotype carrying both F29S and Y106S — 16,384 of them, a whole quadrant of the 2^16
  space — with attrition that is depth-correlated in closed form and reaches 0% at depth
  16, removing the matured antibody the landscape is about. See
  [`cr9114-reprep.md`](cr9114-reprep.md).
* **What it emits:** three results.
  * `corpus.lock` — *calls* `scripts/freeze_cross_protein_corpus.py verify` rather than
    reimplementing it, and reports its exit status and output.
  * `corpus.derived_table` — what `verify` cannot see: the clean file it hashes is
    exactly the file that was read, so a hash match says nothing about whether a raw
    table at the same commit holds more of the design space. Acknowledge a reviewed case
    with `--acknowledge-derived KEY=WHY`; the acknowledgement is recorded in the result.
  * `corpus.prepared` — the sha256 of a prepared table under `scripts/data/`. The lock
    pins the **upstream** corpus; the arms read these derived CSVs, which it never
    hashes. `PREPARED_TABLE_DIGESTS` pins `cr9114_h1` to the digest published in
    `cr9114-reprep.md` §4; unpinned landscapes are reported `NOT-APPLICABLE` with the gap
    named, not passed over.

### 2. `encoding` — encoding compatibility

For any arm handed an encoded design space, `(X.max(0) == X.min(0)).any()` must be false.

* **Motivating defect:** a constant one-hot column is a residue never observed at that
  position. botorch's `normalize` divides by `upper - lower`, so the column becomes
  `0/0 = NaN`, NaN propagates through the first `Linear` layer, the "exclude already
  queried" mask writes NaN over NaN, and `np.argmax` over an all-NaN array returns `0` —
  the same index, forever. See [`alde-compatibility.md`](alde-compatibility.md) §3.
* **Availability is not compatibility.** `probe_alde` checks the checkout, the
  interpreter, the dependencies and the imports, and correctly reports
  `available: true` on all nine landscapes — including the seven where the arm returns
  one index for its entire budget. The result type therefore keeps `detail["available"]`
  (what the probe said) and `detail["compatible"]` (what this check measured) as separate
  keys.
* **Reference figures**, from ALDE's own code: TrpB 0, GB1 0, HER2 55/200, CreiLOV
  265/300, CR9114 288/320, eqFP611 234/260, CR6261 198/220. `REFERENCE_CONSTANT_COLUMNS`
  holds them and the test suite asserts that `encode_landscape` + `constant_columns`
  reproduces the measured counts on CR6261 H1/H9 and eqFP611-blue from the real tables.
* **Relationship to `benchmark_panoptes_isolated.constant_encoding_columns`.** That is the
  same one-line condition enforced *inside* a run. This one is independent of it on
  purpose: the pre-flight gate has to answer before a run starts and without importing a
  4,800-line harness that is under active edit. Consolidating the two is a follow-up; the
  formula, not an implementation, is what they share.
* The check applies to arms whose pipeline normalises the design space by column range —
  `ENCODED_DESIGN_SPACE_ARMS`, the ALDE family. Any other arm is `NOT-APPLICABLE` with
  that stated; `--encoded-arms` extends the set.

### 3. `budget` — distinct measurements against the declared budget

Distinct measurements must equal the declared budget, and a shortfall is **classified**.

* **Motivating defects:** two of them, which must not collapse into one word.
  * **DEGENERATE** — the arm could not propose a distinct candidate; unrecoverable. HER2
    ALDE spent **26 of 125**: the 25 shared starting variants plus one variant measured
    100 times.
  * **RECOVERABLE** — the arm's proposals were ineligible and the round ran short. GB1
    reliability averaged **123.23** (min 120) and TrpB frontier **120.23** (min 116) this
    way. The arm searched; it spent slightly less than the stated budget.
  * **UNCLASSIFIED** — a shortfall with no rejection-reason evidence. It still fails,
    because the two above require opposite responses, and the result names the evidence
    it would need (`selection_trace.csv`'s `reason` column).
* **Evidence, in order of preference.** `round_metrics.csv` now carries
  `backfill_added`, `distinct_proposed` and `degenerate_proposals`, which are the run's
  own determination, so the classifier reuses them instead of inferring. Where they are
  absent — older runs — it falls back to the `reason` column of `selection_trace.csv`,
  then to `requested`/`accepted`/`rejected`, and says UNCLASSIFIED rather than guessing.
  That fallback matters: on the authoritative HER2 run the frontier arm has **no rejected
  rows at all** in the trace, because it offered fewer candidates than the batch size
  rather than having proposals refused.
* **Not static.** `budget_observations_from_run(run_dir)` reads those two files. Without a
  `--run-dir` the check is `NOT-APPLICABLE` and says that the numbers only exist once a
  run has written them.
* **Validated against the real artefacts.** On
  `benchmark_comparison_results/run_20261031_132813_ff1a3423` (HER2, 48 seeds) the check
  reports `alde` **26 of 125, DEGENERATE**; `panoptes_qualify_then_optimize` min 101,
  mean 121.50, **RECOVERABLE**; and `random`, `ridge`, `random_forest`,
  `panoptes_default` all at 125, PASS.

### 4. `endpoint` — endpoint viability

The diversity-adjusted ceiling at the intended `d` must differ materially from the raw
top-1% count, **and** the measurement budget must be smaller than the hit tier.

* **Motivating defect:** on HER2 the d>=2 ceiling equals the raw ceiling at 113, so the
  adjustment is a no-op and the pre-registered primary endpoint silently becomes the
  secondary one. HER2 is the only `DEGENERATE` target in this corpus.
* The degeneracy rule is **imported** from
  [`veritas/calibration/endpoint.py`](../veritas/calibration/endpoint.py), not duplicated:
  `DEGENERATE` when the ceiling equals the raw count, `NEAR-DEGENERATE` under 5% removed,
  `USABLE` otherwise.
* The second condition catches budget saturation, which is a separate way for the endpoint
  to stop discriminating: at 125 measurements both eqFP611 channels (tier 82) and both
  CR6261 landscapes (tier 19) are saturated. A `USABLE` verdict with a saturated budget
  still fails, and `detail["budget_saturates_tier"]` says which condition bit.

### 5. `threshold` — threshold transfer

The qualification bar must be expressed as a **quantile**, must resolve from the same
values array that defines `top1.threshold`, and the manifest must record the quantile, the
resolved absolute and the source.

* **Motivating defect:** `--qualify-threshold` is an absolute, defaulting to 0.70. An
  absolute does not transfer: 0.70 is TrpB's 99.96th percentile — a bar in the extreme
  upper tail — and lies below the global minimum on four landscapes, where it qualifies
  every variant instantly. Every HER2 seed logged "qualified before round 1".
* The conversion is **imported** from
  [`veritas/calibration/qualification.py`](../veritas/calibration/qualification.py).
* Four things are checked against the data rather than believed: the quantile is recorded;
  the recorded absolute really resolves from *this* landscape's values at that quantile;
  `top1.threshold` recomputes from the same array, so the two bars are on one scale; and
  the bar actually gates something — it is above the global minimum, and the exact
  hypergeometric probability that a 25-variant initial draw already qualifies is below
  `GATING_P_MAX`.
* **`GATING_P_MAX` is 0.5, not 1.0, and that is load-bearing.** On HER2 the absolute 0.70
  qualifies 1,027 of 11,276 variants, so the probability that the starting sample already
  clears the bar is **0.908** — not 1. A `p == 1.0` test would have passed the exact defect
  it was written for. It is a judgement bar rather than a derived one, so it is a named,
  overridable constant: a stage satisfied before round 1 on more than half of seeds does
  not gate.
* **Manifest shapes.** `QUALIFICATION_MANIFEST_PATHS` and `TOP1_MANIFEST_PATHS` list where
  the bar is looked for, best shape first. Real manifests today record
  `configuration.components.qualify_threshold: 0.7` and
  `results.top1.threshold_value` — a bare absolute with no quantile and no source, which
  is the defect. A failure names **every path it checked**, never just "not found".

### 6. `row_order` — row-order invariance

Shuffle the input table's row order, re-run one seed, and assert identical results.

* **Motivating defect:** any correct method is invariant to input row order. ALDE was not:
  it returned index 0 forever. Our prepared tables are sorted by fitness descending, so
  **CR9114's row 0 IS the global optimum** — a wholly broken arm with an all-NaN surrogate
  would have posted a perfect `best_measured_fitness` there and beaten a working one. This
  check catches that entire class of index bug, including instances we have not hit.
* **Cheap on purpose:** one seed, one extra run, and a fixed shuffle seed
  (`DEFAULT_SHUFFLE_SEED = 20261031`) recorded in the result so a failure reproduces
  exactly. `shuffled_order(n, shuffle_seed)` regenerates the permutation.
* **Not static.** It needs a runner it can call twice:
  `runner(rows, seed) -> result`, passed as `runner=` or
  `--invariance-runner module:function`. Without one the result is `NOT-APPLICABLE` with
  that stated — and it still reports `detail["row_zero_hazard"]`, which is a property of
  the table: whether row 0 is the global optimum, and what percentile it sits at.

## Genotype parsing

Mutation strings **omit wild-type positions**, so `V183A:F184A:V227A` is a four-site
genotype whose fourth site is wild type. This bug inverted a conclusion earlier in this
project. Everything here works on canonical full-length combos from
[`veritas/calibration/genotypes.py`](../veritas/calibration/genotypes.py) — imported, not
rewritten — and `one_hot` **raises** on genotypes of unequal length rather than zipping or
truncating them.

## What is static and what is not

| check | static? | what it needs |
|---|---|---|
| `corpus` | yes | the corpus mirror and the lock file |
| `encoding` | yes | the prepared table for the landscape |
| `endpoint` | yes | the prepared table for the landscape |
| `threshold` | yes, given a manifest | a `manifest.json`, or `manifest=` |
| `budget` | **no** | a run's `round_metrics.csv` (+ `selection_trace.csv` for the reason split) |
| `row_order` | **no** | a runner the check can call twice on one seed |

`CHECKS_REQUIRING_A_RUN` names the last two in code, so a report can say which part of the
gate a static invocation did not cover.

## Tests

```powershell
.venv\Scripts\python.exe -m pytest tests\test_preflight.py -q
```

Every check has a test that **fails on a constructed bad input** and passes on a good one,
and the bad inputs reuse the real figures — ALDE's 26 of 125, HER2's 113 == 113, CR6261's
198 of 220 constant one-hot columns, 0.70 below a bnAb assay floor. A check that cannot
catch the defect it was written for is worthless, so the figures are the fixtures.

## Follow-up: wiring the gate into the benchmark

Not done here, deliberately — `benchmark_panoptes_isolated.py` and
`workers/math_worker.py` are another agent's files. The integration is three calls:

1. **Before the arm grid runs**, call `run_preflight(library=..., arms=...)` with
   `only={"corpus", "encoding", "endpoint", "threshold"}` and refuse to start on a
   failure, or record the report in the manifest under a `preflight` key. The `encoding`
   result alone would have refused the ALDE arm on seven of nine landscapes before a
   single measurement was spent.
2. **After the run**, call `check_budget(budget_observations_from_run(run_dir), ...)` and
   let the classification drive `run.state`: a DEGENERATE arm is not a null result, while
   a RECOVERABLE shortfall is the protocol working as designed on an unsaturated library.
   `round_metrics.csv` already holds both numbers.
3. **Occasionally**, pass the benchmark's own single-seed entry point as
   `--invariance-runner` so the row-order check runs against the real campaign rather than
   a stand-in.
