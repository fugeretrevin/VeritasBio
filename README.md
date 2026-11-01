# VeritasBio

**An autonomous design agent that refuses to spend your budget on a campaign that cannot
work — because we built one that did, and our own metrics could not tell.**

Meta Global AI Developer Hackathon (with Cerebral Valley), track **Push the Frontier**.
First commit `cd313c4`, 2026-10-26 12:17 ET, the day the event opened. Everything in this
repository was built during the event; `git log --reverse` is the record.

VeritasBio is two things in one repository. The first is an autonomous protein-design
agent: each round a Meta model (`muse-spark-1.3` on `api.meta.ai`) reads the target's
structure and annotations, the assay results so far, and the residuals from its own
previous predictions, then emits a search strategy as schema-validated JSON — which
positions to mutate, which residue clusters, how deep — and the campaign executes that
plan against a real measured fitness landscape. The second is the verification layer that
decides whether the plan is worth paying for. The second exists because of what we found
while measuring the first.

## The finding

When a generative system proposes candidates that are filtered *before* they are scored, a
low score is **ambiguous**. It can mean the model ranked badly, or it can mean nothing the
model produced was ever eligible to be ranked. Aggregate metrics cannot separate those two.
They need opposite fixes. The second is invisible unless you instrument for it.

We hit it. On the CR9114 antibody landscape our own engine's arm:

| what the arm did | figure |
|---|---|
| candidate proposals built | **757,829**, none above mutation depth 4 |
| offered to the eligibility judge | **25,712 rejected, 32 accepted** |
| rejection reasons | **100% `not_in_eligible_universe`**, no other reason once |
| measurement budget spent | **26 of 125** |
| rounds that committed a short batch | **192 of 192** |
| top-1% variants found | **0 of 651** |
| what the aggregate metrics said | mean top-1% hits and best fitness placed it exactly where a weak-but-working method sits |
| what this harness said | **`Run state: PARTIAL`**, arm **`DEGENERATE`** |

Source: [`docs/day1-reachability.md`](docs/day1-reachability.md) §6, reproduced in the
summary block of [`scripts/demo/gate_catches_planner.py`](scripts/demo/gate_catches_planner.py).

The last two rows are the point, and they need stating carefully. **This harness did catch
it** — a run demotes to `PARTIAL` on any degenerate arm, and the per-arm `DEGENERATE`
classification is `veritas/preflight.py` check 3. Both exist *because of this episode*.
The code comment on the demotion says so: "A degenerate run must never leave the state
SUCCESS. The HER2 ALDE run completed cleanly, wrote every file, and was a null arm."

What did not distinguish it is the thing everyone reports: mean top-1% yield and best
fitness. Those two numbers are identical in shape for an arm that searched and lost and an
arm that never searched. The rejection reason that would have separated them was recorded
in the run's own artifacts and used by no reported metric.

The cause was not the model or the surrogate. The candidate generator was hardcoded to one
library geometry — four saturated sites, 20 residues each, depth capped at 4 by a pydantic
validator in `controller/prompts.py`. Six of the nine standard landscapes in this corpus
are binary or mixed-alphabet. The two the generator was validated on, TrpB and GB1, are the
only two of its native shape — and they are also the two the MLDE/ALDE literature
benchmarks on. The `le=4` bound alone put 647 of CR9114's 651 top-1% hits out of reach at
any budget. No metric in the results table distinguished that from a method that searched
and simply lost.

## Verify the Meta model is real

One command, one live call, no mock path:

```powershell
.venv\Scripts\python.exe -X utf8 scripts/demo/gate_catches_planner.py
```

It prints, from the live `Orchestrator` object:

* `strategy_engine` and `strategy_engine_is_mock` — set per round in
  `controller/orchestrator.py::_generate_strategy`, and reset at the top of every round so
  one successful call cannot label a whole campaign live.
* the **provider's own token counters** (`token_usage`), accumulated from the API's `usage`
  block in `_call_llama_api`. Not an estimate.
* the strategy the model actually returned: hypothesis, `max_mutation_depth`, and each
  target cluster with its positions and rationale.

Where the call is made: `controller/orchestrator.py::_call_llama_api` POSTs to
`https://api.meta.ai/v1/chat/completions` with model `muse-spark-1.3` (overridable via
`META_API_URL` / `META_MODEL_ID`), `response_format: {"type": "json_object"}`,
temperature 0.2. The response is parsed into `SearchStrategy`
([`controller/prompts.py`](controller/prompts.py)), a pydantic model — a malformed plan is a
validation error, not a silently degraded round.

Three further ways to check the claim rather than take it:

1. **The arm cannot fall back to a mock.** `require_live_llm=True` makes a failed live call
   fail the arm. Artifact: in the CR9114 run below, `panoptes_muse` reports
   `FAILED, 0 of 12 seeds completed — the arm did not run, which is not a null result`,
   and the log carries the provider's error body. An arm that quietly substituted a mock
   would have posted a plausible number instead.
2. **Token counts are measured.** The TrpB run below records
   `panoptes_muse: 568,189 tokens over 96 model call(s), 5,919/call` in `report.txt` and
   `summary.csv`. Before commit `9583708` the client discarded the provider's usage block
   and every report said "token usage: unknown"; that commit also fixed an arm that
   reported "live Muse" while making **zero** calls.
3. **The model's plan is enforced, not decorative.** `restrict_strategy_to_assayable`
   clips the returned plan to the positions the ordered library actually varies, and the
   demo withholds that scope on the first call so you can see what the planner proposes
   unconstrained.

Requires `LLAMA_API_KEY` (or `META_API_KEY`) in `.env`, and `VERITAS_USE_OS_TRUSTSTORE=1`
on a network that terminates TLS. Costs one call, roughly 5k tokens.

## The four checks

Four questions have to be answered before a campaign's number means anything. Each is code
in this repository, each was written against a defect that actually happened here, and each
names the artifact it read.

**1. Scope — can the assay read what the planner proposed?**
`scripts/demo/gate_catches_planner.py` gate check 1, enforced in production by
`controller/prompts.py::restrict_strategy_to_assayable`. The model sees a protein; the lab
has a plate. A distal allosteric residue may genuinely matter and still be unmeasurable
with the library in hand, and a proposal there cannot be ordered, measured or scored. This
check passes on TrpB today because an earlier fix (`8099060`) derived the planner's
positions from the library spec. It stays, because nothing guarantees the next target's
annotation matches the library someone actually ordered.

**2. Reach — can that plan reach where the answers live?**
Gate check 2. It takes the `max_mutation_depth` the model asked for and asks, per
landscape, what fraction of the top-tier variants sit at or below that depth — using the
landscape's own depth and hit arrays from `veritas/calibration/landscapes.py`, not a second
implementation. One model call, one plan, seven libraries: sound on some, **unable to
succeed** on others. It is the check the model cannot make for itself: where a library's
good variants sit in mutation depth is a property of the measured landscape, not of the
sequence, the structure or the annotation.

> **Be precise about what this check binds on, because we were not.** A reach refusal says
> *a proposer that honours this depth cannot succeed.* That is a statement about the plan,
> and about any system that enforces a depth bound — the `le=4` pydantic bound did enforce
> one, which is how it put 647 of CR9114's 651 hits out of reach. It is **not** a
> prediction about this engine's yield today: `recombine_lead_variants` stacks measured
> leads and carries the candidate pool *past* the declared depth from round 1. Measured
> 2026-11-01 over the real selection path on all nine landscapes, the pre-veto ballot
> spans depths 1–12 on CreiLOV and reaches **99.9%** of its hit set, and the depth ceiling
> costs **1 hit of 1,655** there and **0** on the other eight — instrumentation in
> `scripts/analysis/pool_hit_reachability.py` on branch
> `worktree-agent-a83565bcde59a81f3` (`bb925c9`), not yet merged. An earlier version of
> `veritas/check.py` asserted the stronger reading. It was wrong, the tool now prints the
> distinction itself, and a test pins the specific wrong sentence so it cannot return.

**3. Budget — did the arm actually spend its budget, and if not, why?**
`veritas/preflight.py::check_budget`, reading a run's `round_metrics.csv` and
`selection_trace.csv`. A shortfall is **classified**, not summarised: `DEGENERATE` (the arm
could not propose a distinct candidate — ALDE on HER2 spent 26 of 125, being the 25 shared
starting variants plus one variant measured 100 times), `RECOVERABLE` (proposals were
ineligible and the round ran short), or `UNCLASSIFIED`, which still fails, because the first
two require opposite responses. A degenerate arm is not a null result.

**4. Endpoint — can the metric tell two methods apart at all?**
`veritas/preflight.py::check_endpoint_viability`, with the degeneracy rule imported from
`veritas/calibration/endpoint.py`. On HER2 the diversity-adjusted d≥2 ceiling equals the
raw top-1% ceiling at 113 against a raw count of 113 — 0.00% of the hit set removed against
a 5% floor — so the pre-registered primary endpoint silently becomes the secondary one. No
budget change fixes that; it is a property of the library. **HER2 results are therefore
not quotable from this project**, and we say so rather than reporting them.

Checks 3 and 4 are two of the six in the pre-flight gate. The other four are `corpus`
(every table a landscape reads is hash-pinned), `encoding` (a constant one-hot column makes
botorch's min-max normaliser divide by zero, which is what made ALDE return index 0
forever), `threshold` (the qualification bar must be a quantile that resolves from the same
array as the top-1% threshold), and `row_order` (shuffle the input table, re-run one seed,
assert identical results — our prepared tables are sorted by fitness descending, so
CR9114's row 0 *is* the global optimum and a wholly broken arm would have posted a perfect
best-fitness). Nothing in the gate is a bare boolean: `bool(result)` raises `TypeError`, so
`if check(...):` cannot quietly read NOT-APPLICABLE as "fine". Full documentation:
[`docs/preflight.md`](docs/preflight.md).

```powershell
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library her2 --arms alde,ridge
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library trpb --arms panoptes_default --run-dir <run>
```

## `veritas check` — the gate as a command

`veritas.preflight` grades a campaign that has already run. `veritas check` grades a
**plan**, before anything is spent. Three ceilings, each reported separately and then
composed, because collapsing them would hide which stage binds:

```powershell
# the live Meta planner proposes; the gate decides
.venv\Scripts\python.exe -X utf8 -m veritas.check --library trpb --from-model

# a plan you have on disk, or one typed by hand
.venv\Scripts\python.exe -X utf8 -m veritas.check --library creilov --plan plan.json
.venv\Scripts\python.exe -X utf8 -m veritas.check --library creilov --depth 3 --all-landscapes
```

Exit **0** the plan can succeed · **1** it cannot, *or a ceiling could not be measured* ·
**2** the check itself could not run. An unmeasured ceiling is not an open one, so
`UNKNOWN` never composes to `PASS` and `--skip-viability` cannot clear a plan. A gate that
exits 0 on "I don't know" is not a gate.

`--from-model` makes a real Meta Model API call and has **no mock path**, on purpose: a
mock would prove nothing about what a model proposes. With the assayable scope withheld by
default, because supplying it lets `restrict_strategy_to_assayable` clip the plan before
anything downstream sees it — correct for a scored run, and it makes the catch invisible.

**It audits its own input.** A 100% spatial-veto pass rate has two unrelated meanings:
the veto rejected nothing, or no assayed position mapped and
`filter_spatially_viable_variants` returned its input unchanged. Worse, a structure can
load and be the *wrong* structure, and production validates neither, because
`reference_sequence` is never passed at the call site so `SpatialContactMap`'s identity
guard never fires. So `grade_structure` has three bands rather than two, and they are not
interchangeable:

| library | structure | mapped positions disagreeing with the library's own parent | grade |
|---|---|---|---|
| TrpB | `8EH0` | 1 of 4 (`184: P` vs `F`) — right protein, right numbering, one substitution | **CAVEAT** |
| HER2 | `1N8Z` | 9 of 10 — chain A is the trastuzumab *light* chain; the library's positions are CDR-H3 on the *heavy* chain | **UNTRUSTED** |
| CreiLOV | `8QI9` | 11 of 11 — phototropin / CrPhotLOV1, a non-aligning homolog | **UNTRUSTED** |

Grading TrpB as badly compromised as the other two would be its own kind of dishonesty.
The audit can only ever **downgrade** a verdict, so a bad input can never manufacture a
pass.

## Measured results

Protocol, shared by every arm: one evaluator owns the complete measured fitness table and
selectors never see it; per seed the same 25 starting variants for every method, then 4
rounds of 25 newly selected variants, 125 unique measurements per method per seed including
the start; a method commits its batch before outcomes are revealed. The seed is the unit of
independence and every comparison below is a **paired** test over seeds with
seeds-better/worse/tied reported. A difference of means with no paired test is not a verdict
here — that rule cost this project two wrong verdicts
([`docs/HANDOFF.md`](docs/HANDOFF.md) rule 2.11).

### The live agent on TrpB — 12 paired seeds

`run_20261101_043406_8662bdf0`, commit `9583708`, clean tree, TrpB Tm9D8\* four-site active
site (Johnston et al. 2024), 159,999 eligible variants, top-1% threshold 0.349399 (1,600
qualifying variants).

| arm | peak fitness | sd | top-1% found | sd |
|---|---|---|---|---|
| random | 0.4455 | 0.1591 | 0.92 | 0.79 |
| one-hot ridge | 0.6968 | **0.3012** | 29.17 | 17.02 |
| engine, non-LLM planner (`panoptes_default`) | 0.7288 | 0.0242 | 50.33 | 6.71 |
| **engine, live Meta planner (`panoptes_muse`)** | 0.7421 | **0.0316** | 49.42 | 6.04 |

Paired tests, n=12 (the delta is identical whether or not the shared starting sample is
counted, since every arm starts from it):

| comparison | metric | delta | t | better/worse/tied | verdict |
|---|---|---|---|---|---|
| live agent vs ridge | top-1% | **+20.25** | **+3.77** | 11/1/0 | agent better |
| live agent vs random | top-1% | +48.50 | +27.55 | 12/0/0 | agent better |
| live agent vs ridge | peak fitness | +0.045 | +0.50 | 6/6/0 | **not separable** |
| live agent vs non-LLM planner | top-1% | −0.92 | −0.44 | 6/5/1 | **not separable** |
| live agent vs non-LLM planner | peak fitness | +0.013 | +1.20 | 7/3/2 | **not separable** |

State the fourth row honestly: **the Meta model plans as competently as our tuned
heuristic planner, and on this landscape it does not beat it.** What it does is plan from
the target's annotations and its own residuals rather than from a per-target rule we wrote.

The variance is the result that matters for a real campaign. Ridge's per-seed peak fitness
ranged **0.0646 to 1.0000** (sd 0.3012) — it found the global optimum of the eligible set
on one seed and essentially nothing on two others. The live agent ranged **0.7036 to
0.7903** (sd 0.0316), roughly a tenth of the spread. The engine's advantage here is
**reliability, not peak performance**, which is what matters when you run a campaign once
with real synthesis budget rather than twelve times with different seeds.

Cost, from the provider's own usage block: **568,189 tokens over 96 model calls,
5,919/call**, for 12 full campaigns. No dollar figure is hardcoded; multiply by your rate.

### CR9114-H1 antibody — 12 paired seeds, engine without the live planner

`run_20261101_051322_4bb22522`, same protocol and commit. CR9114 VH 16-site library against
H1 antigen (Phillips et al. 2021), 65,093 eligible variants, 651 qualifying.

| arm | peak fitness | sd | top-1% found | sd |
|---|---|---|---|---|
| random | 9.6373 | 0.0211 | 1.33 | 0.65 |
| one-hot ridge | 9.6264 | 0.0419 | 1.83 | 1.80 |
| **engine (`panoptes_default`)** | **9.6945** | 0.0525 | **6.50** | 2.35 |

| comparison | metric | delta | t | better/worse/tied | verdict |
|---|---|---|---|---|---|
| engine vs ridge | top-1% | **+4.67** | **+7.39** | 12/0/0 | engine better |
| engine vs ridge | peak fitness | **+0.068** | **+4.49** | 11/0/1 | engine better |
| ridge vs random | top-1% | +0.50 | +1.00 | 5/5/2 | **not separable** |
| ridge vs random | peak fitness | −0.011 | −0.90 | 4/5/3 | **not separable** |

One seed reached **9.8345** — the global optimum of the landscape
(`S31N:T57S:T71S:K74I:S75F:T76S:S77N`, mutation depth 7, which the old depth-4 ceiling made
unreachable by construction). Note the third and fourth rows: on this landscape ridge, the
standard baseline, is **no better than guessing**.

Two disclosures on this run. Its state is `PARTIAL`, because the live-planner arm was
requested and failed on all 12 seeds (`target_accession` was `None` for this library — a
configuration defect, fixed afterwards in `c97b4fb`); no figure above comes from that arm.
And the live planner would be refused here anyway — see *Content policy* below.

### Across six landscapes, 48 paired seeds each — the engine loses

The generality sweep, commit `406ef58`, clean tree, budgets matched per arm, 10,000-sample
whole-seed bootstrap. Δ is `panoptes_default − ridge` on top-1% hits excluding the initial
sample. Reproduce with `scripts/analysis/generality_sweep.ps1`.

| landscape | Δ vs ridge | t | better/worse/tied | verdict |
|---|---|---|---|---|
| CreiLOV | **−50.29** | **−32.01** | 0/48/0 | loss on every seed |
| GB1 | **−12.90** | **−5.63** | 6/41/1 | decisive loss |
| eqFP611-red | 0 hits vs ridge's 15.35 | — | — | loss, total |
| eqFP611-blue | +0.73 | +0.77 | 25/18/5 | **not separable** |
| CR9114-H1 | **+3.17** | **+6.10** | 41/5/2 | win — but development data |
| HER2 | — | — | — | **endpoint gate FAILED, not quotable** |

**One win on development data, one tie, three losses.** That is the honest summary of
generality for the selection engine, and it is in this README because a
[pre-registration](docs/preregistration-cross-protein.md) written **before** the sweep ran
predicted across-the-board regression against ridge. It was right. That is why these
numbers are reportable rather than re-narrated: the prediction is in the repository, dated,
and it was not softened afterwards. Two earlier verdicts of ours did not survive the paired
test and are withdrawn in place: eqFP611-blue was published as a win and is a coin flip,
and GB1 was described as "parity at every budget" and is a decisive loss.

### What we then falsified about ourselves

Three explanations for those losses were pre-registered and **measured false**, in one day.
Each had been plausible enough to act on, and two of them were ours.

1. **"The spatial viability veto caused them."** Ablating the veto at 48 paired seeds
   *costs* hits: CreiLOV −2.04 (t=−2.33), GB1 −4.27 (t=−3.22). On both landscapes it was
   **helping**, by raising the candidate pool's hit density 2.8× and 1.8×. On eqFP611-red it
   does fully explain the **zero** — the achievable ceiling there is 0 of 82 hits and the arm
   scored exactly 0.0000 on all 48 seeds, which is arithmetic rather than a measurement —
   but ablated it reaches only 4.21 against ridge's 15.35, so **73% of that loss is not the
   veto**. Pre-registration and results on branch `worktree-agent-a038e3b42e920227a`
   (`3804237`), not yet merged.
2. **"The candidate pool was too shallow."** False on all nine landscapes; see the
   blockquote under check 2 above. The depth ceiling costs 1 hit of 1,655 on CreiLOV and 0
   elsewhere.
3. **"The biophysical featurisation was the problem."** `physchem_features` against one-hot
   is a **null** on CreiLOV standalone: +0.167, p=0.74, better on 4 of 12 seeds.

What CreiLOV binds on is **ranking**, and the dominant cause is a defect in our own training
data: `select_next_micro_batch` injects one row per round for the wild type carrying the
literal label `y = 0.0` *whatever the assay's scale*, at variance `1e-4` against a real
measurement's `1e-2`. CreiLOV's real labels occupy 2.793–4.196, so that row sits 11.6 units
below the minimum — **15.68 label SDs from centre** — and every round's surrogate is fitted
through it. Ablating it is worth **+29.75 hits (t=12.17, 12 of 12 seeds)**; with the prior
and one-hot features as well, **+41.08**, closing −50.08 to −9.00, and best fitness reaches
4.1880 ± 0.0001 against ridge's 4.1851 ± 0.0043.

**TrpB (0.35 label SDs) and GB1 (0.20) are the only landscapes whose assays are genuinely
centred near zero — the only two where this defect is invisible — and they are exactly the
two every quotable engine figure in this repository came from.** They are also the only two
of the candidate generator's native shape, which is a second, independent reason the same
pair flatters this engine. CR9114 violates the scale by 14.14 SDs and the engine *wins*
there, so this is not a complete theory of the result set.

> **Provisional, and marked as such: 12 paired seeds, not 48.** A 48-seed confirmation is
> running, with TrpB and GB1 as built-in negative controls — the theory predicts the
> ablation is a no-op on both, and if it moves TrpB materially then the theory is wrong
> *and* every TrpB figure above changes meaning. The fix is deliberately an **ablatable
> component, not a default**: correcting or dropping that row re-bases every number on this
> page. Do not quote the −9.00 as a result until the confirmation lands.

## Limitations

Read these before quoting anything above.

* **TrpB, GB1 and HER2 are development data** ([`docs/HANDOFF.md`](docs/HANDOFF.md) rule
  2.4). The intervention was designed against failures observed on them, so an improvement
  there is **mechanism evidence, not generalization evidence**. CR9114 is development data
  for the same reason: the generator fix was designed against its measured failure.
* **Zero independent protein validation, and zero wet-lab validation.** Every number here
  is retrospective replay against a published measured landscape. Nothing was synthesised.
* **The live-agent evidence is one landscape and 12 seeds**, not the 48 the cross-protein
  protocol uses. Live planning was not run across the sweep.
* **TrpB's labels are `HISTORICAL_UNKNOWN_PROVENANCE`.** Johnston et al. imputed 871 of the
  159,999 values with a KNN imputer and the prepared CSV carries no per-row flag, so no
  individual label in a TrpB run can be claimed as measured. Do not assume a shared
  contaminated universe makes the comparison between methods unbiased — imputed values sit
  where sequencing coverage was poorest, which is not uniform over the landscape, so the
  bias is differential, not common-mode. Nothing in these runs bounds it.
* **The `threshold` gate check FAILS on the three runs the sweep table is quoted from.**
  `qualification.quantile` and `qualification.source` are not recorded in those manifests.
  Nobody had seen it because `scripts/analysis/generality_sweep.ps1` invokes the gate
  *without* `--run-dir`, which makes both `threshold` and `budget` report NOT-APPLICABLE.
  The bar it concerns is consumed only by the `qualify_then_optimize`/`conditional_probe`
  arms and the failure is symmetric across arms, so it **cannot bias any paired contrast
  above** — but it is an un-passed check that those figures carry, and it is stated
  wherever they are. Not fixed.
* **No landscape measured so far tests the spatial cohesion veto on correct geometry.** On
  GB1 the cohesion rule contributes **zero** rejections at any depth — the library fits in
  an ~8.2 Å sphere — and all 118,766 GB1 rejections are single-site buried-core/steric calls
  firing at depth 1. On CreiLOV the structure is a non-aligning homolog. TrpB is the only
  landscape that could test it, and it has not been ablated. So the veto's measured benefit
  on CreiLOV is **pool density, not biophysics**: a random filter enriching density 2.8×
  would have scored the same.
* **`panoptes_default` changed meaning on 2026-10-31.** k-DPP diversity and the two-stage
  GP now default **off**, because both were measured to hurt top-1% yield. Always carry the
  flag state with the number; every run's `manifest.json` records components enabled *and*
  disabled.
* **Steps 4 and 5 of the production loop are not exercised by the benchmark.** Gatekeeper
  authorisation and Modal ESM-2 execution are skipped, because `Orchestrator.run_stream()`
  offers no interface for substituting an experimental assay for ESM scoring. The benchmark
  measures planning and selection.
* **The compliance and wet-lab paths are architecture, not demonstrated capability.** The
  Opentrons OT-2 protocol generator, Twist/IDT order CSVs, Benchling LIMS adapter and the
  21 CFR Part 11 Merkle hash chain exist in `workers/` and `ledger/` and are tested, but
  none of them was exercised in any result on this page. Treat them as structure.
* **`scripts/run_real_data_benchmark.py` (avGFP) is not a predictive-skill measurement.**
  It sets `all_eval_feats = np.array(train_x)` and correlates the model against its own
  training targets, so its rho of +0.9996 describes fit, not skill, and the same run finds
  0 top-1% variants. Earlier revisions of this README reported it as "high ranking
  fidelity". It is not. Use `benchmark_panoptes_isolated.py`.
* **The ALDE baseline is incompatible, not unavailable, on seven of nine landscapes.** It
  reports `available: true` and then silently returns index 0 forever. Availability is not
  compatibility; `docs/alde-compatibility.md` has the mechanism.
* **Only the cited run artifacts are in this repository.** `benchmark_comparison_results/`
  is gitignored (`.gitignore:246`), because the benchmark writes hundreds of run
  directories. The runs behind every figure on this page are copied verbatim into
  [`evidence/`](evidence/README.md) — manifest, summary, per-round metrics, report and log
  for each — so a reader can recompute any paired statistic here without an API key or a
  re-run. That directory also names three traps in the data, including the fact that one
  run's failed live arm leaves round-0 rows which a naive aggregation will misread as a
  result.

## Content policy, and why we did not work around it

The live planner is **refused on antibody targets**. A round-2 prompt carrying measured
affinity evidence against influenza hemagglutinin returns, from the Meta API:

```
Meta Model API 400 for model 'muse-spark-1.3': {"error":{"code":"content_policy_violation",
"message":"The response was filtered due to the prompt triggering our content management
policy.","param":"prompt","type":"invalid_request_error"}}
```

Recorded verbatim in [`evidence/run_20261101_052018_8755fea6/benchmark.log`](evidence/run_20261101_052018_8755fea6/benchmark.log).
The round-1 prompt, which carries no affinity evidence, succeeds; the refusal arrives when
the evidence block does. We consider this correct dual-use behaviour and the project does
not route around it: `scripts/demo/gate_catches_planner.py` offers enzyme and
fluorescent-protein targets only (`OFFERED = ("trpb", "gb1", "creilov", "eqfp611_blue",
"eqfp611_red")`), and because `require_live_llm=True` the benchmark's live arm **fails
loudly** on an antibody rather than substituting a mock planner and reporting a number.

This is a real constraint on LLM-planned biological design and practitioners should know
about it: the targets where a learned planner would be most valuable commercially are
exactly the targets where the provider will — correctly — decline to plan.

## How to run

```powershell
# environment
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
copy .env.example .env     # then add LLAMA_API_KEY or META_API_KEY

# 1. the live model call and the two agent-plan checks (one call, ~5k tokens)
.venv\Scripts\python.exe -X utf8 scripts/demo/gate_catches_planner.py
.venv\Scripts\python.exe -X utf8 scripts/demo/gate_catches_planner.py --library gb1

# 2. the pre-flight gate, no run required for four of its six checks
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library her2 --arms alde,ridge

# 3. the live agent on TrpB - reproduces run_20261101_043406_8662bdf0
.venv\Scripts\python.exe benchmark_panoptes_isolated.py --library trpb `
  --methods random,ridge,panoptes_default,panoptes_muse `
  --seeds 12 --base-seed 0 --rounds 4 --initial-size 25 --batch-size 25 `
  --enable-live-muse --max-model-calls 12 --provenance-mode acknowledge-unflagged

# 4. the engine on CR9114, no model calls - reproduces run_20261101_051322_4bb22522
.venv\Scripts\python.exe benchmark_panoptes_isolated.py --library cr9114_h1 `
  --methods random,ridge,panoptes_default `
  --seeds 12 --base-seed 0 --rounds 4 --initial-size 25 --batch-size 25 `
  --provenance-mode acknowledge-unflagged

# 5. the six-landscape sweep (long)
pwsh scripts/analysis/generality_sweep.ps1

# tests
.venv\Scripts\python.exe -m pytest -q
```

Each run writes a fresh folder under `benchmark_comparison_results/` containing
`manifest.json` (authoritative, machine-readable, carries a configuration fingerprint and
**no** outcomes), `report.txt`, `summary.csv`, `round_metrics.csv`, `selection_trace.csv`,
`candidate_trace.csv` and the log. `--help` lists every flag.
[`docs/HANDOFF.md`](docs/HANDOFF.md) records the suite at 674 passing, 2 skipped as of
2026-10-31.

## Repository map

| path | what |
|---|---|
| `benchmark_panoptes_isolated.py` | the benchmark harness and the arm definitions, including `panoptes_muse` (live Meta planning + residual reflection). Imports production modules read-only; writes only into a fresh run folder. Start at the module docstring. |
| `controller/orchestrator.py` | the 6-step agent loop. `_call_llama_api` is the Meta Model API client; `_generate_strategy` is round planning; `_generate_reflection` is residual analysis between rounds. |
| `controller/prompts.py` | prompt construction and the `SearchStrategy` / `ReflectionDiagnosis` pydantic schemas, plus `restrict_strategy_to_assayable` (check 1). |
| `controller/gatekeeper.py` | biophysical validation and the fail-closed spend authorisation (architecture; not exercised by the benchmark). |
| `workers/math_worker.py` | the selection engine: GP surrogate, candidate generation, GFlowNet sampler, Pareto and k-DPP stages. `select_next_micro_batch` is the entry point. |
| `veritas/preflight.py` | the six-check pre-flight gate (checks 3 and 4 above). |
| `veritas/calibration/` | landscapes, genotype canonicalisation, thresholds, endpoint degeneracy. `genotypes.py` is the single canonical genotype implementation — mutation strings omit wild-type positions, which once silently inverted a conclusion here. |
| `veritas/datasets/` | the library/dataset schema and validating ingestion. Imports nothing from `workers` or `controller`, so a new library is describable without touching the optimiser. |
| `veritas/diagnostics/` | coverage, candidate traces, depth returns, and the token-gated retrospective path that is the only code allowed to read true fitness. |
| `scripts/demo/gate_catches_planner.py` | the live-model demo and gate checks 1 and 2. |
| `scripts/data/` | the nine prepared landscape CSVs, pinned by `cross_protein_corpus.lock.json`. |
| `scripts/analysis/` | the sweep and analysis scripts. Read its README first. |
| `docs/` | see below. |
| `ledger/`, `ui/`, `web/`, `server.py`, `main.py` | the product surface: SQLite/PostgreSQL ledger with the Part 11 hash chain, Gradio dashboard, React frontend, FastAPI bridge. Architecture; not exercised by any result here. |

Documents a judge should read, in order: [`docs/HANDOFF.md`](docs/HANDOFF.md) (current
state, what is quotable, what is void, and the §2 rules that are load-bearing because each
one already produced a wrong result here) ·
[`docs/HANDOFF-2026-10-31-reachability.md`](docs/HANDOFF-2026-10-31-reachability.md) (the
generator defect, the fix, and the measurement) ·
[`docs/day1-reachability.md`](docs/day1-reachability.md) (depth tables and the void-run
record) · [`docs/preflight.md`](docs/preflight.md) (the gate) ·
[`docs/preregistration-cross-protein.md`](docs/preregistration-cross-protein.md) (the
frozen protocol, which predicted the losses above) ·
[`docs/alde-compatibility.md`](docs/alde-compatibility.md) ·
[`docs/dataset-inventory.md`](docs/dataset-inventory.md) ·
[`docs/corpus-audit.md`](docs/corpus-audit.md). Several root-level documents
(`ARCHITECTURAL_STATE_AND_VALUE_OVERVIEW.md`, `INNOVATION_ROADMAP.md`, `PRODUCT.md`,
`FRONTIER_GAP_RESOLUTION_AND_EXECUTION_PLAN.md`) predate the measurements on this page;
where they disagree with this README or with `docs/HANDOFF.md`, those two are correct.

## Provenance

* **First commit** `cd313c4`, 2026-10-26 12:17 ET — the day the hackathon opened. All work
  on this page was built during the event. `git log --reverse --date=iso` is the record,
  and every run's `manifest.json` records the commit it ran at and whether the tree was
  dirty.
* **Landscape data is published third-party data**, not ours: TrpB (Johnston et al. 2024),
  GB1 (Wu et al. 2016), HER2/trastuzumab (Mason et al. 2021), CR9114 and CR6261 (Phillips
  et al. 2021), eqFP611 blue and red (Poelwijk et al. 2019), CreiLOV (Chen et al. 2023),
  avGFP (Sarkisyan et al. 2016).
  `scripts/data/cross_protein_corpus.lock.json` pins the CombinGym mirror at commit
  `c5dccd514a9acfe1089a5f1555faea98afd78019`; that pin is a scientific commitment and does
  not move. CR9114 was re-prepared from the raw eLife table rather than CombinGym's clean
  workbook, which drops 16,384 genotypes with depth-correlated attrition
  ([`docs/cr9114-reprep.md`](docs/cr9114-reprep.md)).
* **Baselines are the published implementations or standard ones**: scikit-learn
  `Ridge(alpha=1.0)` and `RandomForestRegressor` on one-hot features, and the authors' ALDE
  implementation (jsunn-y/ALDE), availability-probed and never silently replaced by a
  homemade method.
* **Outcomes never reach selection.** Anything that reads true fitness is retrospective
  only and gated behind an explicit access token in
  `veritas/diagnostics/retrospective.py`; a configuration fingerprint must identify a run
  without consulting outcomes, and there is a test asserting a 32-arm sweep yields 32
  distinct fingerprints with no outcome leakage.

Built for the Meta Global AI Developer Hackathon. Foundation models: Meta `muse-spark-1.3`
via the Meta Model API, and Meta FAIR ESM-2 (architecture path). Infrastructure: Modal Labs.
