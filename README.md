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

## A second mode: qualify, then optimize

The default policy optimises for **breadth** — how many top-1% variants a campaign finds.
That is the right objective when you want a panel of candidates, and the wrong one when you
want *one very good* variant. `panoptes_qualify_then_optimize` runs the default policy
unchanged until a measured variant clears a qualification bar, then switches the remaining
budget to greedy Hamming-1 search around the incumbent:

```powershell
.venv\Scripts\python.exe -X utf8 benchmark_panoptes_isolated.py --library trpb `
  --methods ridge,panoptes_default,panoptes_qualify_then_optimize --qualify-quantile 99.96
```

Because every trajectory is **identical to the default until the protected event has
occurred**, the reliability curves are preserved by construction — and they measured
byte-identical at every budget. The bar is a quantile of the eligible universe, not an
absolute: an absolute does not transfer across landscapes, which is what the gate's
`threshold` check exists to catch (0.70 is TrpB's 99.96th percentile and sits *below the
global minimum* on four of the nine landscapes).

| | default | qualify-then-optimize |
|---|---|---|
| P(≥0.50) and P(≥0.70) at every budget | — | **identical** |
| seeds reaching ≥0.80 | 3/48 | **13/48** |
| seeds reaching ≥0.90 | **0/48** | **8/48** |
| best variant found | 0.8326 | **1.0000** |
| top-1% hits | 45.23 | 39.90 (**−5.33**, t=−3.20) |
| best fitness | 0.7457 | 0.7570 (t=0.94, **not** significant) |

> **These figures were measured on a superseded tree and need one re-run before they are
> quoted.** `panoptes_default` on this block measures 51.42 today, not the 45.23 above —
> `DPP_DIVERSITY_ENABLED` and `TWO_STAGE_GP_ENABLED` changed default since. The *shape* of
> the result is what is being claimed here, not the decimals.
>
> One more caveat that is a property of the landscape, not the method: the 99.96 bar is
> favourable on TrpB specifically because 88% of campaigns qualify by measurement 50,
> leaving three batches to spend on the frontier phase. On a slower landscape there would
> be no budget left to switch into, and the mode would degenerate to the default.

Honest description: **threshold attainment at the top improves; average peak does not.**
The mean best-fitness gain is not significant and loses 25 of 38 decided seeds. It buys the
tail for a real and measured cost in breadth, which is a trade a campaign owner should make
deliberately rather than a strict improvement. It is shipped as a *mode*, not a default, for
exactly that reason.

This is also the mode the baseline panel independently argues for: `greedy_measured`, pure
local exploitation with **zero learned parameters**, is the strongest arm on CR9114. Where
the assay floors and most of the library carries no gradient, exploiting a measured lead
beats modelling.

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
| engine vs **ridge** | top-1% | **+4.67** | **+7.39** | 12/0/0 | better **than ridge** |
| engine vs **ridge** | peak fitness | **+0.068** | **+4.49** | 11/0/1 | better **than ridge** |
| ridge vs random | top-1% | +0.50 | +1.00 | 5/5/2 | **not separable** |
| ridge vs random | peak fitness | −0.011 | −0.90 | 4/5/3 | **not separable** |

One seed reached **9.8345** — the global optimum of the landscape
(`S31N:T57S:T71S:K74I:S75F:T76S:S77N`, mutation depth 7, which the old depth-4 ceiling made
unreachable by construction).

> **This is a win over `ridge`, and on this landscape that is a low bar. Measured
> 2026-11-01 at 48 paired seeds with a seven-arm panel, budgets gate-verified matched, the
> engine LOSES on CR9114 to four standard baselines:**
>
> | arm | top-1% | engine Δ vs it | t |
> |---|---|---|---|
> | `greedy_measured` — **no model at all** | **11.60** | **−6.60** | −6.25, 10/36/2 |
> | `random_forest` | 11.40 | −6.40 | −4.84 |
> | `knn_hamming` | 10.98 | −5.98 | −7.18 |
> | `gbm_onehot` | 9.00 | −4.00 | −3.90 |
> | `panoptes_default` | 5.00 | — | — |
> | **`ridge`** — the registered reference, **6th of seven** | 1.83 | +3.17 | +6.10 |
> | `random` | 0.88 | | |
>
> `greedy_measured` has **zero learned parameters** — lexicographic Hamming distance to the
> fitness-ordered measured variants — and it beats the engine by 6.60 hits on 36 of 48
> seeds. The obvious objection, that local-exploit arms return one dense basin, was
> **tested rather than argued**: `greedy_measured` and `random_forest` leads do dissolve at
> d≥3, but **`knn_hamming` beats the engine at raw, d≥2 and d≥3** (+1.79, t=+5.6), so the
> conclusion survives all three endpoints.
>
> And the "ridge is no better than guessing" reading below needs **both** tests, not one.
> At 48 seeds `ridge − random` is Δ=+0.958, t=+2.46, Wilcoxon p=0.044 — but **21/19/8
> seeds, sign test p=0.875.** The t is carried by a heavy tail (sd 2.60 on a mean of 1.83);
> on a seed-majority basis it is a coin flip, while every other arm separates from random
> at sign p ≤ 1.8e−6. Ridge's weakness here is **assay censoring, not library size**: the
> −log₁₀K_D readout floors, so most of the library carries no gradient. Branch
> `worktree-agent-a8325f2a02fbcdc90` (`493826b`), not merged.

Note the third and fourth rows: on this landscape ridge, the standard baseline, is **no
better than guessing** on a seed-majority basis — read with the blockquote above, which
gives both tests and the four arms that beat us here.

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

> **Narrowed 2026-11-01 by a seven-arm baseline panel at 48 paired seeds** (branch
> `worktree-agent-a8325f2a02fbcdc90`, `493826b`; budgets gate-verified matched, all arms
> spending all 125 measurements on all 48 seeds). The sweep above uses only `random`,
> `ridge` and the engine, and **`ridge` is not a uniformly fair reference**: it is 1st of
> seven on CreiLOV, 2nd on TrpB, and **6th of seven on CR9114**. Against the strongest
> arm per landscape:
>
> | landscape | strongest baseline | it scores | engine | Δ | t | b/w/t |
> |---|---|---|---|---|---|---|
> | TrpB | `gp_onehot` | 38.54 | 51.42 | **+12.88** | +3.57 | 31/16/1 |
> | CR9114-H1 | `greedy_measured` | 11.60 | 5.00 | **−6.60** | −6.25 | 10/36/2 |
> | CreiLOV | `ridge` | 57.23 | 6.94 | **−50.29** | −32.01 | 0/48/0 |
>
> So the honest summary is **not** "one win, one tie, three losses" but: **no win against
> the strongest available baseline on any landscape except TrpB, and none at all off
> development data.** CreiLOV is the only non-development landscape in the panel and the
> engine loses there to every arm except `random`. The TrpB margin is real, large and
> diversity-robust (+8.69 at d≥2 on 48 of 48 seeds) — and it is on development data *and*
> on the one landscape where the fabricated `y = 0.0` WT row cannot show itself. That is
> evidence the method works on a 4-site saturated library with a zero-centred assay: the
> single configuration in which nothing known to be broken is exercised.
>
> **Settled 2026-11-01 on ONE seed block, all nine arms in one run folder** (branch
> `worktree-agent-aa53d65b7e6a10b8f`, `9158f6a`; TrpB on the designated held-out seeds
> 48–95, GB1 on 0–47; every arm spent all 125 measurements on all 48 seeds, gate-verified).
> The two results above lived on *different* seed blocks and could not be combined. Here
> they are measured together, and every premise figure reproduces **exactly** — TrpB's
> +12.875/t=+3.57/31-16-1, the d≥2 +7.750/t=+13.40, GB1's −12.896/t=−5.63/6-41-1, the
> redundancy ratios 7.06×/4.54×.
>
> **On TrpB the engine beats all seven baselines on all three endpoints, every contrast
> separable:**
>
> | endpoint | strongest baseline | it scores | engine + one-hot | Δ | t | b/w/t |
> |---|---|---|---|---|---|---|
> | raw | `gp_onehot` | 38.54 | **56.42** | **+17.88** | +5.08 | 37/11/0 |
> | **d≥2** | `ridge` | 6.40 | **14.25** | **+7.85** | **+15.41** | **48/0/0** |
> | d≥3 | `ridge` | 2.06 | 4.79 | **+2.73** | +10.75 | 43/2/3 |
>
> The strongest baseline is a *different arm under each endpoint* — `gp_onehot` on raw,
> `ridge` at d≥2 and d≥3 — and `gp_onehot` buys its raw lead with redundancy, 7.06 raw hits
> per distinct one against ridge's 4.54. Pre-registered as F5 and confirmed.
>
> **The splice was worth 2.7 hits.** Naively combining the two earlier runs gives
> 51.42 + 7.71 − 38.54 = **+20.59**; measured on one block it is **+17.88** — inside the
> ±6 band the pre-registration declared in advance. The splice got the sign and the order
> of magnitude right and the number wrong, which is precisely what cannot be known without
> measuring, and is why the withdrawn 6.3× claim is not rehabilitated by this.
>
> **On GB1, `ridge` was never the strongest raw baseline. `gp_onehot` is, at 46.02 against
> ridge's 33.29.** The panel never ran GB1, so **every GB1 contrast in this README and in
> HANDOFF §1 is against the project's third-strongest raw arm, and the loss is worse than
> recorded**: engine+one-hot is **−21.31** against `gp_onehot` (t=−9.10, 7/41/0). At d≥2 it
> is +0.79 (not separable) and at d≥3 **+0.75** (t=+4.23) — so GB1 is a decisive *raw* loss
> and a shallow *diversity-adjusted* win, and the endpoint decides the verdict.
>
> One-hot closes **33.4%** of GB1's ridge gap, against a pre-registered 33% with an 80%
> interval of [10%, 60%]. So **GB1's loss is not a feature-space defect** — and together
> with the two earlier campaigns it is not the spatial veto (ablating it *costs* GB1 4.27
> hits, and the cohesion rule contributes **zero** GB1 rejections at any depth) and not the
> fabricated WT row (GB1 sits 0.20 label SDs from centre and its 48-seed ablation is +0.542
> at t=+0.41). **Three mechanisms eliminated; GB1 remains unexplained.**
>
> Validity, demonstrated rather than asserted: after merging two branches, **21,120 cells
> across all 22 behavioural columns** agree exactly against `main` while
> `source_tree_sha256` differs — the code changed and the behaviour did not. The
> `_one_hot_matrix` refactor was checked against `main`'s inline loop as an independent
> oracle over the **full** universes: **12,799,920 TrpB and 11,948,800 GB1 cells bitwise
> identical.**
>
> Also re-measured on the same held-out block: `panoptes_default` on TrpB is **51.42**, not
> the 45.21 recorded earlier, and the d≥2 figure is **+7.75 (t=+13.4)**, not +7.93
> (t=+17.1) — different arms under one name, since `DPP_DIVERSITY_ENABLED` and
> `TWO_STAGE_GP_ENABLED` changed default between the two measurements. And the
> front-loading claim of "6.3× the best baseline at 50 measurements" is **withdrawn as
> incorrect**: it divided by ALDE (2.25) when ridge is the best baseline at that budget
> (2.83 → 5.0×), `panoptes_default` itself gives 3.7×, and the two ends of the quoted range
> use different engine arms. See [`docs/HANDOFF.md`](docs/HANDOFF.md) §1.

**One win on development data, one tie, three losses** — ~~the honest summary of
generality~~, now narrowed by the panel above; read that blockquote in its place. It is in
this README because a
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
3. **"The biophysical featurisation was the problem."** Measured twice, independently, and
   a **null on CreiLOV both times**: +0.167 (p=0.74, 12 seeds) and **+0.79 (p=0.33, 48
   paired seeds)**. The 48-seed arm comes with the default path proven unchanged rather
   than asserted — **15,840 cells across all 22 behavioural columns of `round_metrics.csv`
   agree exactly** against the parent commit over 3 arms × 48 seeds × 5 rounds, the only
   differences being wall clock. Branch `worktree-agent-a14c69a80e0afa44a` (`a779422`).

   **It is not a null everywhere — but the gain is REDUNDANT, and that correction is the
   sharpest result of the day.** On **TrpB** the one-hot GP adds **+5.00 raw top-1% hits
   at t=+5.33** (40/6/2) — and **+0.10 at d≥2, t=+0.29, 22/19/7: NOT SEPARABLE.** The
   engine's redundancy rises 3.63× → 3.96× and the diversity-adjusted count is flat at
   14.15 → 14.25. So it finds five more hits per campaign and **nothing at all** once hits
   within Hamming 1 of an already-counted hit are removed. **Any claim that the one-hot
   surrogate improves the engine must be scoped to the raw count.** On GB1 it is +4.31 raw
   but confounded, because the GFlowNet reward scores candidates *through* this GP and
   deep-enrichment moved 0.690 → 0.801.

   This is a finding that could not be obtained from the earlier, separate measurements —
   it needed all arms in one run folder on one seed block, and it is why the splice below
   was not trusted.

   **Better surrogate rho does not predict a better campaign.** Rho gain ran GB1 (+0.107) >
   CreiLOV (+0.077) > TrpB (+0.025); campaign gain runs the other way round — TrpB > GB1 >
   CreiLOV, Spearman −0.5 over the three. Offered mechanism: the redundant 10-D block is
   33% of the biophysical columns on TrpB and GB1 but only 12% on CreiLOV, where one-hot is
   the *narrower* encoding — a geometry a fixed-matrix rho protocol cannot see.

What CreiLOV binds on is **ranking**, and the dominant cause is a defect in our own training
data: `select_next_micro_batch` injects one row per round for the wild type carrying the
literal label `y = 0.0` *whatever the assay's scale*, at variance `1e-4` against a real
measurement's `1e-2`. CreiLOV's real labels occupy 2.793–4.196, so that row sits **2.793
raw units below the minimum — 11.62 label SDs, or 1.99× the entire label range — and
15.68 label SDs from the centre** of the distribution, and every round's surrogate is
fitted through it.

**Confirmed at 48 paired seeds: +31.771, t=+35.88, better on 48 of 48 seeds, bootstrap
[+30.000, +33.438].** It accounts for **63.2%** of CreiLOV's −50.292 deficit. With the
prior and one-hot features as well, **+42.583** (t=+46.24, 48/0/0). Branch
`confirm-wt-anchor-48` (`90472c2`), pre-registration committed before any arm ran; the
measured value lands inside the pre-registered [+24, +36] and the bootstrap's **lower**
bound clears the dominance threshold of +25.04 that was fixed in advance.

**The engine still loses CreiLOV**: −18.521 with the row ablated, −7.708 with all three
changes. None of this is a CreiLOV win.

**TrpB and GB1 are unaffected, as the negative control required** — +1.458 (t=+1.60,
bootstrap [−0.292, +3.250]) and +0.542 (t=+0.41, [−2.021, +3.105]), both straddling zero.
**Every quotable TrpB and GB1 figure survives.** One precision the agent insisted on, and
it is right: *"not separable" is not "identical"* — TrpB's point estimate sits just under
its minimum detectable effect of 1.826, so a real one-to-three-hit effect is not excluded.

> **The label-scale theory is FALSIFIED, and this was pre-registered.** The theory said
> distance from zero predicts the defect's effect. **CR9114 violates the scale by 14.14
> SDs and by 2.47 label-ranges — *worse* than CreiLOV's 1.99 — and the effect is
> Δ=+0.042, t=+0.09, bootstrap [−0.854, +0.938].** A tight null against an MDE of 0.924,
> not a power failure. The measured ordering by relative effect is CreiLOV (+457.9%) ≫
> eqFP611-blue (+17.1%) > TrpB (+2.9%) ≈ GB1 (+2.7%) > **CR9114 (+0.8%)** — CR9114 at the
> wrong end of *both* candidate covariates. The agent's own primary CR9114 prediction was
> wrong and its pre-registration said in advance what that would mean.
>
> So the honest statement is narrower than the one this section first carried: **the row
> is the dominant cause of the CreiLOV loss specifically. It is not a general explanation
> of the engine's losses, and distance from zero does not predict its effect.** **GB1's
> −12.896 remains unexplained** and has now survived two pre-registered root-cause
> campaigns. eqFP611-red was **REFUSED** — its arm failed the budget gate at 71 of 81.

> **RULE 2.4 HAS FIRED: CreiLOV IS NOW DEVELOPMENT DATA.** `wt_anchor` was built against
> CreiLOV's measured failure, so the +31.771 is **mechanism evidence, not generalization
> evidence**. Combined with the baseline panel's finding that CreiLOV was the only
> non-development landscape in it, **this project now has no non-development landscape
> among its measured results**, and a new independent one must be reserved before any
> generality claim is made again.
>
> The only generalization evidence the confirmation produced is **eqFP611-blue: +1.625 at
> t=+7.83 on 38 of 48 seeds** — small, real, and on a landscape the intervention was not
> designed against.

> **Still an ablatable component, not a default.** Correcting or dropping that row
> re-bases every number on this page, so it ships as a **named policy version** or not at
> all. The recommendation, drafted and deliberately unexecuted: use the measured WT label,
> omit the row where no measured label exists, give it a real measurement's `1e-2`
> variance rather than `1e-4`, and re-measure all nine landscapes at 48 seeds first —
> precisely because "not separable" is not "identical".
>
> **One term of the three-way stack is actively disputed, so the −7.708 figure is weaker
> than it reads.** At 48 paired seeds with the prior ablated, the one-hot increment measures
> **−0.19, CI [−1.04, +0.65], p=0.547** — an interval that **excludes** the +3.67 measured
> at 12 seeds, so it is not an underpowered estimate of the same quantity. The two
> implementations differ in **scope**: one moves all five featurisation call sites including
> `AssayAdaptedHead`, the other only the GP's fit/predict inputs. With the prior ablated
> `zero_shot_scores` is zero and `AssayAdaptedHead` becomes the dominant non-GP term, so
> **"one-hot is worth +7.08 once the prior is gone" may be an `AssayAdaptedHead` effect
> attributed to the GP.** The prior ablation itself does replicate: 6.94 → 10.94, **+4.00**
> at 48 seeds against +4.25 at 12. The decisive experiment is both components on one
> branch and one seed block, and it has not been run.

## Limitations

Read these before quoting anything above.

* **Every landscape here is now development data, and that is the limitation to read
  first.** TrpB, GB1 and HER2 by [`docs/HANDOFF.md`](docs/HANDOFF.md) rule 2.4; CR9114
  because the generator fix was designed against its measured failure; **CreiLOV as of
  2026-11-01**, because `wt_anchor` was built against *its* measured failure. Rule 2.4 has
  now fired twice. An improvement on development data is **mechanism evidence, not
  generalization evidence**.

  **The only generalization evidence in this project is eqFP611-blue: +1.625 at t=+7.83 on
  38 of 48 seeds**, from the `wt_anchor` confirmation, on a landscape the intervention was
  not designed against. Everything else is mechanism.

  And the pinned corpus cannot currently supply a replacement: the four landscapes never
  promoted to the registry are all too small or too sparse for a 125-measurement budget
  and a 1% tier — `sacas9` is combinatorially complete but only 1,296 variants (13 hits at
  the 1% tier), `spcas9` has 597 real rows, `hiv1_protease` drops to 0.6% saturation by
  depth 9, and `rhla` holds nothing beyond depth 3. Reserving a real one needs either a
  documented endpoint amendment for a small complete universe or a new dataset added to the
  lock (which rule 2.5 permits; the pin itself may not move). This is the top item on the
  roadmap, and until it is done **no generality claim should be made from this repository.**
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
