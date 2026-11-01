# Phase 0 — specification freeze

**Status:** frozen specification. Written and committed **before any Phase 1–4 run, and
before any outcome of a `frontier_pool_v2` or `phase_schedule_v2` arm exists.** It
discharges the Phase 0 item of [`HANDOFF.md`](HANDOFF.md) §7:

> *"Phase 0 — freeze the specifications before any large run. Define `frontier_pool_v2`,
> the fixed phase schedule, the `design_space` / `assay_available_pool` /
> `oracle_evaluable_pool` semantics, primary and secondary depth covariates, and the arm
> validity statuses. Preserve the legacy policy IDs. Commit all of it before viewing any
> outcome."*

**Frozen on:** 2026-11-01.
**Authority:** the decisions in [`HANDOFF.md`](HANDOFF.md) §5 are the external reviewer's
and are **not relitigated here**. Where §5 settles something this document cites it and
specifies only how it is computed, recorded and checked. Where §5 leaves a value or a
boundary rule undetermined, this document **fixes** it and marks the provenance column
`this doc`, so that a later reader can tell an expert decision from an engineering choice.
Where a decision genuinely has no answer available, it is in §10, not invented.

**Amendment rule, inherited from
[`preregistration-cross-protein.md`](preregistration-cross-protein.md) §17.** Amendments
are appended to §11 with a date and a reason; the text above is not edited in place. **Any
amendment written after an outcome of the arm it governs has been examined must say so in
its first sentence.** That clause is already load-bearing here — see §5.5.

**What this document is not.** It is not a pre-registration of a prediction. It fixes
*definitions, provenance and admissibility*, not an expected result. The prediction
remains [`HANDOFF.md`](HANDOFF.md) §8: across-the-board regression against ridge,
recorded in advance and not to be softened.

---

## 0. Conformance vocabulary, and the two findings that reweight this document

### 0.1 How to read a requirement

| word | meaning |
|---|---|
| **MUST** | a conforming run satisfies it; a run that does not is non-conforming and its numbers are not admissible evidence for the claim the requirement protects |
| **MUST NOT** | as above, inverted |
| **RECORDS** | the named key appears in `manifest.json` (or the named CSV column) with the stated semantics; absence is non-conformance, not a default |
| **VERIFIED BY** | the artefact-only procedure a later reader runs to decide conformance without re-running anything |

Every section below carries all four. A requirement that cannot be checked from the run's
own artefacts is not a requirement in this document — that is the lesson of
[`preflight.md`](preflight.md): every one of five invalidating defects was already written
down in an artefact we produced, and nothing asserted on it.

### 0.2 Two findings of 2026-11-01 that change the weighting

**(F1) The binding constraint is probably the training distribution, not the surrogate.**
An offline study established that the engine's featurisation is *not* worse than one-hot as
an information carrier (better on GB1 by **+0.103** with a linear read-out), and that
eqFP611-red's in-run Spearman of **−0.0069** reproduces offline from the **training sample
the engine selects for itself** — identically for one-hot features. So representation is
exonerated as the primary suspect and the roadmap points at §5.1's candidate source and
the acquisition stage. **Consequence for this document: §2 (`frontier_pool_v2`) carries
more weight than anything about representation, and no Phase 1 acceptance criterion may be
stated in terms of features.**

**(F2) Every published comparison is against ridge, which is not the hardest baseline.**
`random_forest` beats both ridge and the GP on **3 of 9** landscapes
([`day1-reachability.md`](day1-reachability.md) §4 / `scripts/analysis/surrogate_headroom.py`:
gb1 0.314, eqfp611_blue 0.780, eqfp611_red 0.587, cr9114_h1 0.906). **Consequence: §7.3
requires `random_forest` in any comparison that claims a win.** Whether the
*pre-registered primary contrast* changes from "versus ridge" is a pre-registration
amendment and is in §10.

---

## 1. The three pool semantics

§5.1 ruled that three concepts this project had conflated **must be separate in code and
in prose**, because *"missing historical labels do not mean a sequence is biologically
impossible."* Everything else in this document is defined in terms of these three, so they
come first.

### 1.1 Definitions

| pool | definition | depends on |
|---|---|---|
| **`design_space`** | every genotype the **library specification permits**: the Cartesian product over declared positions of each position's declared allowed residues, including depth 0. A property of the design, not of any dataset. | `LibrarySpec.positions`, `LibrarySpec.allowed_aa_per_position` (or all 20 per site when `None`) |
| **`assay_available_pool`** | the genotypes a selector **may propose** in a given run. Retrospectively this is the oracle-coverage constraint (§1.3); prospectively it is `design_space` ∩ what can actually be ordered and assayed. | the run's selection mode and `oracle_evaluable_pool` |
| **`oracle_evaluable_pool`** | the genotypes with a **usable measured label** under the curation rule of `preregistration-cross-protein.md` §3.3 — identifier, depth and phenotype all present, parseable and finite. This is the set `view.untested(revealed)` draws from. | the prepared table and the eligibility rules |

**The three coincide only in a fully characterised library and MUST NOT be treated as
conceptually identical** (§5.1). HER2 is the standing counter-example: its 25,114
absent-label rows are outside `oracle_evaluable_pool` and **inside** `design_space`, and
its `design_space` exceeds its `oracle_evaluable_pool` by seven orders of magnitude.

### 1.2 Sizes, carried in and asserted

Carried in from existing documents, **not re-derived here**. A conforming run recomputes
both columns and **MUST** reproduce them; a mismatch is a configuration error, not noise.

| library | sites | `design_space` | `design_space` formula | `oracle_evaluable_pool` | source for column 5 |
|---|---|---|---|---|---|
| trpb | 4 | 160,000 | 20⁴ | 159,999 | day1 §4 |
| gb1 | 4 | 160,000 | 20⁴ | 149,361 (prepared CSV 149,360 — D-3) | day1 §4, prereg §1.2 |
| her2 | 10 | **182,463,840,000** | 20·20·19·19·18·13·12·10·9·5 | 11,277 rows, **11,276** mutants | day1 §6b, prereg §1.2 |
| cr9114_h1 | 16 | 65,536 | 2¹⁶ | 65,094 rows, **65,093** mutants | prereg §5.1 |
| cr6261_h1 | 11 | 2,048 | 2¹¹ | 1,886 | day1 §4 |
| cr6261_h9 | 11 | 2,048 | 2¹¹ | 1,841 | day1 §4 |
| eqfp611_blue | 13 | 8,192 | 2¹³ | 8,191 | day1 §4 |
| eqfp611_red | 13 | 8,192 | 2¹³ | 8,191 | day1 §4 |
| creilov | 15 | 184,320 | 2¹²·3²·5 | 165,427 | day1 §4, prereg §3.2 |

HER2's 1.824×10¹¹ reproduces day1 §8c's "1.8e11" from the declared alphabets alone, and
CreiLOV's 184,320 reproduces `preregistration-cross-protein.md` §3.2. Those two agreements
are what make this table checkable rather than decorative.

> **The wild-type row is a convention, not a fact, and the documents differ.** day1 §4
> counts rows; prereg §1.2 and §6.3 count mutants and exclude depth 0 from both shells
> while retaining it in the universe. A conforming artefact **RECORDS**
> `pools.wt_row_convention` ∈ {`rows`, `mutants_excluding_depth0`} and every count beside
> it. Two anchors computed under different conventions is precisely the drift a freeze
> exists to prevent (prereg §17, 2026-10-31).

### 1.3 The oracle-coverage constraint, stated as a constraint

Restricting selection to historically measured variants is acceptable **for retrospective
benchmarking only if it is disclosed as an oracle-coverage constraint and applied
identically to every method** (§5.1).

- In a retrospective run, `assay_available_pool := oracle_evaluable_pool`.
- The disclosure **MUST** appear in the manifest as a limitation string, not only as a
  derived key, and **MUST** state the ratio `|oracle_evaluable_pool| / |design_space|`.
  On HER2 that ratio is 6.18×10⁻⁸; on TrpB it is 0.99999.
- The constraint **MUST** be applied identically to every arm in the run. An arm proposing
  outside `assay_available_pool` and having the surplus rejected is **not** the same
  constraint as an arm handed the pool — see §1.4, which is the whole of the CR9114
  starvation.

### 1.4 Where each is computed, and which one each existing comparison actually used

| code site | pool it realises |
|---|---|
| `LibrarySpec.positions` + `allowed_aa_per_position` + `alphabet_for()` (`benchmark_panoptes_isolated.py:200–291`) | `design_space` — **computable today but never materialised or recorded.** Phase 0 adds the record (§8). |
| `Universe` / `Evaluator` ingestion, `Evaluator.eligible_ids` (`:1620`) | `oracle_evaluable_pool` |
| `CandidateView.variant_ids`, `CandidateView.untested(revealed)` (`:1648–1664`) | `oracle_evaluable_pool` minus the measured set |
| `_judge(...)`, reason `not_in_eligible_universe` (`:3627`) | the judge that enforces `assay_available_pool` **after** proposal |
| `select_next_micro_batch(..., candidate_pool=...)` when `pool_restricted` or `candidate_source == "enumerate"` (`:3398`) | `assay_available_pool` handed in directly |

**What each existing comparison actually used.** This table is a deliverable of Phase 0 in
its own right: the asymmetry in its last column is the mechanism behind
[`day1-reachability.md`](day1-reachability.md) §6a.

| arm | proposed from | judged against | effect |
|---|---|---|---|
| `random` | `oracle_evaluable_pool` − measured (`view.untested`) | same | zero rejections by construction |
| `ridge`, `random_forest` | `oracle_evaluable_pool` − measured | same | zero rejections by construction |
| `panoptes_default`, generative mode | the **generator's reach inside `design_space`** — radius-bounded, alphabet-bounded | `oracle_evaluable_pool` | on CR9114, 25,712 rejected / 32 accepted, **100 % `not_in_eligible_universe`**; the arm reached 26 of 125 |
| `panoptes_default`, `pool_restricted` (HER2) | `assay_available_pool` handed in | same | no rejections; `generate_candidate_variants` never called |
| `panoptes_default --production-candidate-source enumerate` | `assay_available_pool` handed in | same | the baselines' task exactly; "strictly easier than generation" per the arm's own note (`:2790–2797`) |
| `panoptes_qualify_then_optimize`, generative | distance-1 shell of the incumbent **∩ `oracle_evaluable_pool`** (`:3276–3286`) | same | shell exhaustible; ≤76 untested neighbours on TrpB |
| `panoptes_qualify_then_optimize`, `pool_restricted` | `assay_available_pool` (`:3269–3274`) | same | in 41 of 48 HER2 seeds it selected the identical 125 variants as the reliability arm — not a distinct policy |
| `alde` | the encoded `oracle_evaluable_pool`, **plus** the exported truth table (`Evaluator.export_truth_table`) | same | the one documented exception to hidden results |

**So no published engine-versus-baseline comparison on a generative library was a
comparison on one pool.** That is a statement about the task, not about the policy, and it
**MUST** accompany any re-use of a pre-Phase-1 generative number.

### 1.5 What a conforming run records, and how it is verified

**RECORDS** under `configuration.pools`:

```
pools.design_space_size                 int        (exact; may be a big integer)
pools.design_space_formula              str        the per-site alphabet product, written out
pools.oracle_evaluable_pool_size        int
pools.assay_available_pool_size         int
pools.assay_available_pool_rule         enum       {oracle_coverage_retrospective,
                                                    design_space_prospective}
pools.oracle_coverage_ratio             float      oracle_evaluable / design_space
pools.wt_row_convention                 enum       {rows, mutants_excluding_depth0}
pools.identical_across_arms             bool       true iff every arm in the run saw the
                                                   same assay_available_pool
```

**VERIFIED BY:** (a) recomputing `design_space_size` from `LibrarySpec` and comparing with
§1.2; (b) asserting `assay_available_pool_size == oracle_evaluable_pool_size` whenever the
rule is `oracle_coverage_retrospective`; (c) asserting
`pools.identical_across_arms == true` for any run whose arms are compared to each other;
(d) for each arm, cross-reading `selection_trace.csv`'s rejection-reason histogram —
a non-zero `not_in_eligible_universe` count **proves** that arm proposed outside
`assay_available_pool` and therefore did not face the same constraint as the baselines.

---

## 2. `frontier_pool_v2` — the candidate source

§5.1: *"Locality means proximity to promising **measured** sequences within the permitted
candidate pool — not exactly one substitution."* Build an adaptive neighbourhood: multiple
measured leads as centres, unmeasured candidates in increasing distance order, radius
expanding until a shortlist capacity is reached, deduplicated across leads, deterministic
tie-breaking, with a **fixed documented fraction** of capacity reserved for candidates
outside the local neighbourhoods — *"part of the policy, not an undocumented fallback that
fires after failure."*

### 2.1 Inputs, and what it is forbidden to see

`build_frontier_pool_v2` is a **pure function** of:

| input | from |
|---|---|
| `library` | `LibrarySpec` — positions, parent residues, per-site alphabets |
| `assay_available_pool` | §1, as a list in the universe's fixed lexicographic order |
| `measured` | `revealed.measured` — **observed labels only** |
| `batch_size` | the protocol |
| `seed`, `round_index` | the run |
| `params` | §2.3 |

It **MUST NOT** receive the evaluator, the truth table, `top1.threshold`, the qualification
quantile, or any object through which true fitness of an unmeasured variant is reachable
(rule 2.9). **VERIFIED BY** a signature test: the function's parameters cannot name any of
those objects, and a test constructs it with a sentinel evaluator that raises on any
attribute access.

### 2.2 The algorithm

1. **Canonicalise.** Every genotype is converted to a full-length combo with
   `veritas.calibration.genotypes.canonical_combo(variant_id, positions, parent_residues,
   strict=True)`, and every distance is `genotypes.hamming`, which **raises** on unequal
   lengths. Mutation strings omit wild-type positions (rule 2.1); nothing here zips or
   truncates, and **no sixth implementation is written** (§5.7 — there are already five).
2. **Leads.** Rank `measured` by observed value descending; tie-break by canonical combo
   ascending. Take the top `L`. The incumbent best is always lead 0.
3. **Shells.** For `d = 1, 2, …, S` (S = `library.n_sites`; **there is no depth cap** —
   the radius ceiling is the library's own width, never a schema literal: that literal
   `le=4` put 647 of 651 CR9114 hits out of reach at any budget), and for each lead,
   collect the members of `assay_available_pool` at Hamming distance exactly `d` from that
   lead that are not in `measured`.
4. **Fill order.** Ascending `d`; within a `d`, round-robin over leads in lead rank order,
   so one lead cannot monopolise capacity; within a `(lead, d)` cell, ascending canonical
   combo. Deterministic at every level.
5. **Dedup.** By canonical variant id across all leads. A candidate reachable from several
   leads is attributed to the **nearest** lead; ties go to the lowest lead rank. The
   attribution is recorded, because it is what makes the effective radius meaningful.
6. **Local capacity.** Stop the local fill at `ceil((1 − f) · C)` candidates.
7. **Exploration reserve.** Fill `floor(f · C)` from `assay_available_pool` − measured −
   local shortlist, by **uniform sampling without replacement** under
   `random.Random((seed, round_index, "explore"))`. Uniform is chosen because it is the
   only rule that assumes nothing about where hits sit, and because day1 §4b established
   that uniform sampling over the pool is **unhandicapped** where the local operator is
   not. A depth-stratified reserve is a **named alternative arm**
   (`frontier_pool_v2_stratified`), to be measured, never assumed (§10).
8. **Capacity transfer, recorded.** If the local shells exhaust before their capacity, the
   unused local capacity transfers to the exploration component and the transfer is
   recorded as `capacity_transferred_to_exploration`. If the whole pool exhausts, the
   shortlist is returned **short**. It is never padded, and a short shortlist flows into
   the budget accounting of §6, where it becomes a status — not into a silent fallback.
9. **Rebuild cadence.** The pool is rebuilt **once per round**, after that round's
   measurements are revealed, so the leads are current. It is **not** rebuilt inside the
   backfill loop; a backfill request (`:3687–3720`) is served from the same round's
   shortlist, continued, so backfilled picks remain the policy's own ordering.

### 2.3 Frozen parameters

| parameter | symbol | frozen value | provenance |
|---|---|---|---|
| shortlist capacity | `C` | `10 · batch_size` (= **250** at batch 25) | §5.1 ("~250 for a batch of 25 — an engineering starting point, not an optimum") |
| lead count | `L` | **8** | `this doc`. Carried forward unchanged from v1's inline literal (`math_worker.py:1808`) **deliberately**, so the lead count is not confounded with the neighbourhood change. Varying it is a named sensitivity arm, not a silent edit. |
| reserved exploration fraction | `f` | **0.20** | `this doc`. §5.1 requires *a fixed documented fraction* and does not set it. 0.20 gives 50 of 250 at batch 25 — large enough that the reserve is measurable per round, small enough that the arm remains a locality policy. **No evidence selects this value**; it must not be tuned on the final-test landscapes (prereg §4.4). |
| max radius | `d_max` | `library.n_sites` | §5.1 ("expand the radius until a shortlist capacity is reached"); the ceiling is the library, per day1 §2 |
| locality-claim threshold | `λ` | **0.50** | `this doc`, to make §5.1's presentation rule checkable — see §2.6 |

### 2.4 How it differs from the current generator

| axis | current generator (post-2026-10-31 fixes) | `frontier_pool_v2` |
|---|---|---|
| source of candidates | constructed from heuristics inside `design_space`, then filtered by the judge | drawn **from `assay_available_pool`**, so `not_in_eligible_universe` is zero by construction on an enumerable library |
| depth reach | effective depth ∈ {3, 4} (`prompts.py:100` `le=64` now, but `max(strategy.max_mutation_depth, 3)` at `math_worker.py:1178` and the lead-skip at `:1197` still bound it) | `1 … S`, expansion driven by capacity |
| leads | 8, and **a lead at or beyond `max_depth` is skipped entirely** (`:1197–1199`), so a measured depth-9 variant contributes nothing | 8, with no depth predicate on a lead |
| WT reversion | absent; `generate_hamming1_lead_neighbors` implements it correctly and is never called | implicit — a reversion is simply a pool member at distance `d` from the lead |
| position / alphabet truncation | `available_pos[:6]`, `target_aas[:10]`, numerically lowest sites | none; shells are complete at each `d` up to capacity |
| budget | stages demand ~17,954 slots against `max_candidates=2500`; cluster stages get zero whenever a lead exists | one explicit capacity `C`, allocated by the rule in §2.2 |
| exhaustion | silent short batch, or a fallback not visible per round | recorded capacity transfer, then a recorded short shortlist |

### 2.5 How it differs from `--production-candidate-source enumerate`

`enumerate` hands the arm **the entire untested `assay_available_pool`** — the identical
task `RandomSelector` and `RidgeSelector` always perform. Its own note says it: *"strictly
easier than generation; success here does not validate generative proposal"*
(`:2790–2797`).

| | `enumerate` | `frontier_pool_v2` |
|---|---|---|
| shortlist size | the whole untested pool (159,974 on TrpB at round 1) | bounded at `C` ≈ 250 |
| locality | **none**; there is no neighbourhood and no radius to report | an ordered neighbourhood with a measured effective radius |
| what a win would validate | ranking, i.e. surrogate + acquisition — **not** generation | a bounded locality policy, with the locality claim checkable by §2.6 |
| feasibility where the pool is not enumerable | undefined in general (on HER2 it coincides with `pool_restricted` because that library *is* its pool) | defined wherever `assay_available_pool` can be sampled |
| role in the study | the task-symmetry control; keeps surrogate/acquisition separable from generation | the Phase 1/Phase 4 production candidate source |

Both are retained. `enumerate` remains the control that isolates ranking from generation,
and per (F1) it is now the **more** informative of the two for diagnosing the training
distribution.

### 2.6 The presentation rule, made checkable

§5.1: *"If most selections come from large distances or the broader-search component, the
result must not be presented as demonstrating the benefit of local optimisation."*

Define, per (arm, landscape), over **adaptive rounds only**:

```
locality_share = (accepted measurements attributed to the LOCAL component
                  at distance <= 2 from their attributed lead)
                 / (all accepted adaptive measurements)
```

**A result may be presented as evidence for local optimisation only if
`locality_share >= λ` (0.50) on that landscape, and the share MUST be quoted next to the
claim.** Below λ the same numbers are still reportable — as a result about a bounded
mixed-source policy, which is what it would then be.

### 2.7 Policy IDs, and the legacy freeze

§5.1: *"Keep `frontier_legacy_v1` alongside `v2`; never silently modify the old one and
reuse its results."* Phase 0 therefore freezes a **policy ID registry**.

| ID | meaning | status |
|---|---|---|
| `panoptes_default` | production planning + selection, default strategy | **legacy, frozen** |
| `panoptes_muse` | as above with live Muse planning | **legacy, frozen** |
| `panoptes_qualify_then_optimize` | baseline until a measured variant clears an absolute/quantile bar, then the Hamming-1 frontier policy | **legacy, frozen** |
| `panoptes_conditional_probe` | as above with the conditional-probe second phase | **legacy, frozen** |
| `frontier_legacy_v1` | the second-phase candidate source of the two IDs above: `frontier_hamming1` generatively, `frontier_pool_restricted` on a sparse library (`:3274`, `:3286`) | **named here**; behaviour unchanged |
| `frontier_pool_v2` | §2 | **new** |
| `frontier_pool_v2_stratified` | §2.2 step 7 alternative reserve | **new, unmeasured** |
| `phase_schedule_v2` | §3 | **new** |

**A legacy ID's behaviour MUST NOT change. If it must change, a new ID is minted.** The
reason is recorded in [`HANDOFF.md`](HANDOFF.md) §1: `panoptes_default` already changed
meaning under one name when `DPP_DIVERSITY_ENABLED` and `TWO_STAGE_GP_ENABLED` went to
`False`, so *"the 45.21 above and any figure produced after that commit are different arms
under one name."*

**Arm identity is therefore the tuple**, and a bare policy ID is not a configuration
description:

```
arm_identity = (policy_id,
                component_states,              # the EFFECTIVE states, rule 2.12
                candidate_source_id,           # frontier_legacy_v1 | frontier_pool_v2 | enumerate | pool_restricted
                phase_schedule_id,             # legacy_qualification | phase_schedule_v2
                budget_schedule)
```

### 2.8 What a conforming run records

**RECORDS** under `configuration.candidate_source`:

```
candidate_source.id                     enum      see §2.7
candidate_source.version                str       semver; bumped on any behavioural change
candidate_source.params                 {C, L, f, d_max, lambda}
candidate_source.rebuild_cadence        enum      {per_round}
candidate_source.reserve_rule           enum      {uniform, stratified}
candidate_source.outcome_isolation      str       the assertion of §2.1, and the test name
```

**RECORDS** per round, in `round_metrics.csv` (new columns):

```
candidate_source_id          lead_count              shortlist_size
effective_radius_max         effective_radius_p50    local_shortlist_size
reserve_shortlist_size       capacity_transferred_to_exploration
n_accepted_from_local        n_accepted_from_reserve
n_accepted_by_distance       (JSON histogram, distance -> count)
```

**VERIFIED BY:** (a) `candidate_trace.csv`'s `source` column must agree with
`n_accepted_from_local + n_accepted_from_reserve` per round — day1 §9a is the reason this
is checked and not assumed: the ballot was 6,506 gflownet against 909 combinatorial while
the audit was scoped to the combinatorial generator, and *"scope a code audit to the
measured dominant path"* is rule 7.5 of that handoff; (b) recomputing
`locality_share` from the distance histogram; (c) re-running
`build_frontier_pool_v2` from the recorded inputs and asserting a byte-identical
shortlist — determinism is checkable because the function is pure.

---

## 3. The fixed phase schedule

§5.2 removed landscape quantiles from online decisions **entirely**: *"A prospective lab
cannot know the 99.96th percentile before measuring most of its library… The quantile is
fine as a retrospective evaluation threshold and unacceptable as an online decision
threshold."*

### 3.1 The schedule

Under the 125-measurement protocol: **25 shared random, then 25 reliability-oriented, then
75 frontier-oriented. Switch at 50 total measurements.**

| phase | id | measurements | candidate source | selection policy |
|---|---|---|---|---|
| 0 | `shared_initial` | round 0, 25 | n/a | `random.Random(seed).sample` over `assay_available_pool`, identical for every arm in the seed |
| 1 | `reliability` | round 1, 25 | the arm's configured source | the reliability policy: `select_next_micro_batch` with the arm's component states |
| 2 | `frontier` | rounds 2–4, 75 | the arm's configured source, **in frontier configuration** | §3.2 |

**Generalisation to other budgets, fixed now.** §5.2 states the schedule only for 125, and
the sweep actually ran budgets of 81, 113 and 19. The rule is **positional**: the
transition occurs at the end of adaptive round 1, i.e. after `initial_size + batch_size`
measurements. On the 125 protocol that is exactly **50**, as §5.2 fixes. The resolved value
is recorded as `phase_transition_at_measurements`, so no reader has to re-derive it.

### 3.2 What "frontier-oriented" means, concretely

§5.2 names the phases by policy flavour, and the only frontier policies that exist are the
ones §5.1 replaces. The schedule is therefore specified as a **configuration delta**, and
the set of knobs a phase may touch is closed:

| knob | phase 1 (`reliability`) | phase 2 (`frontier`) |
|---|---|---|
| candidate source | arm's configured source, full capacity split `f` | same source, **`f` unchanged** |
| acquisition | the arm's configured acquisition mode | **greedy** (exploitative) |
| k-DPP diversity | the arm's `component_states` | **unchanged** — not a phase knob |
| two-stage GP / DKL | the arm's `component_states` | **unchanged** — not a phase knob |
| prior, Pareto | unchanged | unchanged |

**A phase MUST NOT change any knob outside this table.** The restriction exists because
*"changing candidate generation, qualification and stopping at once makes attribution
impossible"* (§5.3), and because the 2026-10-31 HER2 movement of t = −8.29 could not be
attributed precisely when two things moved across a commit range
([`HANDOFF.md`](HANDOFF.md) §1).

Greedy is named for phase 2 on measured grounds, not taste: on TrpB greedy acquisition
measured **+11.27 top-1 % (t = +12.20)** against EI+Thompson. That is development-landscape
evidence (rule 2.4) and is recorded as the reason for a configuration choice, not as a
result.

### 3.3 The transition rule, and the proof that it is oracle-free

The phase decision is a function

```python
def phase_for(round_index: int,
              measurements_at_round_start: int,
              schedule: PhaseSchedule) -> Phase: ...
```

**whose signature cannot name the evaluator, a fitness value, a threshold or a quantile.**
That arity is the conformance proof, and it is checkable by `inspect.signature`. A test
**MUST** assert that the module defining it imports nothing from
`veritas.calibration.qualification`, nothing from the evaluator, and no landscape table.

Properties this buys, each one a §5.2 requirement:

- independent of the unknown fitness distribution;
- works across assay scales, so the bnAb floor at −log₁₀(K_D) = 7.0 is irrelevant to it;
- **cannot be satisfied by the initial sample** — the defect behind *"every HER2 seed
  logged qualified before round 1"* at an exact hypergeometric probability of 0.908;
- needs no global-optimum estimate.

It is **deliberately less sophisticated** than what it replaces, and that is the point
(§5.2): it removes one source of dataset-specific behaviour so that the frontier
acquisition strategy itself can be evaluated.

**Rename.** This is a **`phase_transition`**, not a fitness qualification. The word
"qualification" is reserved for an actual user-defined experimental target (§5.2).
`panoptes_qualify_then_optimize` keeps its ID (§2.7) and keeps its legacy meaning; the new
arm is a new ID.

### 3.4 The CR6261 qualification floor was dissolved, not decided

[`HANDOFF.md`](HANDOFF.md) §6: *"**The CR6261 qualification floor — dissolved, not
decided.** The question was whether to floor the qualifying set at a minimum count or
declare the stage inapplicable below some library size. §5.2 removes the global-quantile
bar from online decisions entirely, so there is no floor to set. **Do not re-raise it.**"*

Consequences for conformance:

- `configuration.components.qualify_minimum_set_floor.enabled` stays **`False`** under
  `phase_schedule_v2`, and its `note` **MUST** be changed from "a protocol decision that
  must be taken in advance" to a statement of dissolution citing §5.2 and §6. The present
  note invites exactly the re-raise that §6 forbids. (Text change only;
  `benchmark_panoptes_isolated.py` is another agent's file — Phase 2 lands it.)
- A `phase_schedule_v2` arm **MUST NOT** record a `qualify_quantile` or a
  `qualify_threshold_absolute` as an **online** parameter. Both keys remain, under
  `results.retrospective.*`, for the evaluation threshold of §3.5.
- The observation that at the default quantile the qualifying set is **one variant** on the
  CR6261 landscapes — i.e. qualification ≡ finding the global best — is retained as a
  recorded property of those landscapes, not as an argument for a floor.

### 3.5 The quantile survives, renamed and relocated

The oracle quantile is kept **only as a clearly named retrospective legacy analysis**
(§7, Phase 2). It moves to `results.retrospective.qualification_quantile` with the
resolved absolute, the source, and the qualifying count beside it, and it is reachable in
code only through `veritas/diagnostics/retrospective.py` behind that module's explicit
access token (rule 2.9).

### 3.6 Observed-progress switching is a later, separate arm

§5.2: only **after** benchmarking the fixed schedule, investigate observed-progress
switching — switch when a measurement improves on the best initial measurement by a
predefined margin, with a forced checkpoint so a campaign cannot stay in phase 1 forever —
and *"treat it as a separate arm, not an assumed improvement."* It is **out of scope for
Phase 0 and Phase 2**, and its margin is deliberately not fixed here (§10).

### 3.7 `v2` is a new arm on every landscape, TrpB included

§5.2: *"this makes `v2` a genuinely new arm **on every landscape including TrpB**, since
the old switch was data-dependent and the new one is positional."*

So **no TrpB figure carries over to a `phase_schedule_v2` arm**, and the bit-identity
argument that validated the 2026-10-31 generator fixes on TrpB (day1 §9e) does **not**
apply here. A conforming Phase 2 report **MUST NOT** present TrpB as a negative control
for the schedule change.

### 3.8 What a conforming run records

**RECORDS** under `configuration.phase_schedule`:

```
phase_schedule.id                           enum   {legacy_qualification, phase_schedule_v2}
phase_schedule.rule                         enum   {positional}
phase_schedule.phases                       [{phase_id, first_round, last_round,
                                              n_measurements, candidate_source,
                                              acquisition}]
phase_schedule.transition_at_measurements   int    50 on the 125 protocol
phase_schedule.transition_inputs            list   MUST equal ["measurements_at_round_start"]
phase_schedule.oracle_free_assertion        str    the §3.3 signature test's name
```

**RECORDS** per round, in `round_metrics.csv`: `phase`, `phase_rule_id`,
`measurements_at_round_start`, `phase_transition_fired` (bool),
`phase_transition_reason` (a literal such as
`"positional: 50 >= transition_at_measurements"`).

**VERIFIED BY:** (a) every seed of an arm **MUST** transition at the same round — a
positional rule that varies by seed is not positional, and that single assertion
discriminates `phase_schedule_v2` from `legacy_qualification` from the artefacts alone;
(b) `phase_transition_inputs == ["measurements_at_round_start"]`;
(c) no `qualified before round N` note appears in any `phase_schedule_v2` arm's notes.

---

## 4. Depth covariates — primary and secondary

§5.4: *"Preserve `depth <= floor(S/2)` as the primary H1 covariate. Add a
**distribution-relative secondary** — enrichment below the median eligible-library
mutation depth — plus a continuous depth analysis rather than relying on a binary split."*

### 4.1 Primary — unchanged, and already pre-registered

Shell `1 <= depth <= floor(S/2)`; enrichment = shell hit rate ÷ overall hit rate on the
landscape's tier. Depth 0 is **excluded from both shells and retained in the universe and
the tier** (prereg §6.3). Nothing about it changes here. Lives in
`scripts/freeze_cross_protein_corpus.py stats` and `veritas/calibration/`.

Validation gate, rule 2.2 — before any new number, these two **MUST** reproduce:

```
TrpB shallow shell (1 <= depth <= 2): 168 / 2,242 = 7.4933%
GB1  shallow shell (1 <= depth <= 2):  72 / 2,167 = 3.3226%
```

If either fails, the parser is wrong and nothing downstream is trustworthy.

### 4.2 Secondary — distribution-relative, specified

§5.4 names the quantity and not its boundary rule. Fixed here:

| element | specification | provenance |
|---|---|---|
| the pool the median is taken over | **`oracle_evaluable_pool`**, mutants only, depth 0 excluded | `this doc`; §5.4 says "eligible-library", which is §1's `oracle_evaluable_pool` and not `design_space` — the distinction is the point of §1 |
| `m` | the median mutation depth of that pool; reported to 4 decimals, and **not rounded** before use | `this doc` |
| the shell | `{v : 1 <= depth(v) < m}` — **strict** inequality, so a tie at the median falls outside the shell | `this doc`; a strict rule needs no separate tie convention |
| enrichment | shell hit rate ÷ overall hit rate on the same tier as the primary | §5.4 by analogy with §4.1 |
| degeneracy guard | the covariate is **undefined** for a landscape when the shell holds > 80 % of the pool, or when `floor(m) >= max observed depth` | prereg §6.3's guard, applied unchanged |

`m`, the shell size, the shell hit count and the guard verdict **MUST** be published for
every landscape **before** any analysis using the secondary covariate.

**Limitation to document, verbatim from §5.4:** *a depth distribution computed from
historically measured variants is a property of the **sampled library**, not of feasible
sequence space.* In §1's terms: the secondary covariate is a property of
`oracle_evaluable_pool`. On HER2, where that pool is 6.18×10⁻⁸ of `design_space`, the
covariate says almost nothing about the design.

### 4.3 Continuous depth analysis

Reported alongside both binary covariates, for every landscape, so that a smooth
relationship can be distinguished from a threshold artifact (§5.4):

1. the complete **per-depth** table — variants, hits, hit rate, enrichment, max score —
   in the form of day1 §3, which is the object that retracted "hits start at depth 9";
2. a logistic regression of top-tier membership on depth, and separately on `depth / S`,
   each reporting the coefficient, its interval and `n`;
3. a monotonicity statement (CR9114's enrichment rises monotonically with depth; CR6261-H9
   is zero throughout its shallow shell and reaches 22.9 % at depth 9).

Inference remains on the **primary** covariate. The continuous analysis is descriptive and
is labelled so.

### 4.4 Reporting rule

**Results are reported under whichever covariate is more favourable *and* whichever is
less** (§5.4). Both appear in every table. An endpoint or covariate that was secondary at
this freeze stays secondary whatever it shows (prereg §11.3).

**The secondary covariate MUST NOT be described as blind.** Published figures already
reveal depth-related patterns (§5.4), and prereg §16 O-13 settles that at most one of H1's
four positions is genuinely blind, arguably zero.

### 4.5 The amendment licence §5.4 relies on is no longer available — flagged, not resolved

§5.4 permits the amendment *"because no outcomes have been inspected, provided the original
is archived, the amendment is dated before analysis, and results are reported whichever
covariate is more favourable."*

**That precondition is false as of 2026-11-01.** The generality sweep has run and its
per-landscape effect sizes against ridge have been examined and published
([`HANDOFF.md`](HANDOFF.md) §1: CreiLOV −50.29, GB1 −12.90, eqFP611-red 0 vs 15.35,
eqFP611-blue +0.73, CR9114 +3.17, HER2 gate failed). Four of those are final-test
landscapes.

The conditions that **can** still be met are met by this document: the original is archived
(prereg §6.3 is unedited and in git history), and the amendment is dated before any
analysis **using the secondary covariate**. The condition that cannot be met is the first
one. So, per prereg §17, the amendment row for the secondary covariate **MUST** state in
its first sentence that it was written after cross-protein results were examined, and the
H1 rows for the six landscapes already read out **MUST** be labelled exploratory under the
secondary covariate. Whether that is sufficient, or whether a new landscape must be
reserved under prereg §4.4, is a **human decision** and is in §10 — not resolved here.

### 4.6 What a conforming analysis records

**RECORDS** in the stats artefact (`scripts/data/cross_protein_corpus_stats.json`) per
landscape: `depth_covariate.primary.{cut, shell_size, hits, hit_rate, enrichment, ci}`;
`depth_covariate.secondary.{median_depth, pool, cut_rule, shell_size, hits, hit_rate,
enrichment, ci, degeneracy_verdict}`; `depth_covariate.continuous.{per_depth_table,
logit_depth, logit_depth_fraction, monotone}`.

**VERIFIED BY:** reproducing the two rule-2.2 anchors; recomputing `m` from the prepared
table; and confirming that every landscape carries all three blocks, including the ones
where the covariate is `undefined` — an undefined covariate is a recorded result, never an
absence (`preflight.py`'s `CheckResult` invariant).

---

## 5. Arm validity statuses

§5.8, *implement this immediately*: every benchmark arm carries an explicit status, and
**only `VALID_COMPLETE` arms at equal budget enter the primary matched-budget comparison.**

### 5.1 The vocabulary

Frozen as §5.8 states it, six members, no additions:

```
VALID_COMPLETE · VALID_BUDGET_LIMITED · INVALID_NUMERICAL ·
INVALID_PROTOCOL · UNAVAILABLE · PROVENANCE_UNRESOLVED
```

Every non-`VALID_COMPLETE` status **MUST** carry an `invalid_reason` from a closed
enumeration, because the four real instances of §5.3 are distinguishable only at that
level:

```
degenerate_proposals · ineligible_proposal_starvation · eligibility_shortfall ·
endpoint_not_informative · no_adaptive_round_executed · policy_degraded_to_baseline ·
numerical_collapse · arm_absent · source_tree_changed_mid_run ·
component_state_unrecorded · concurrent_run_contention
```

This is how the void results of 2026-10-31 **get represented rather than deleted** (§5.8).

### 5.2 The decision procedure

Computed per **(arm, landscape, budget schedule)** from that run's own artefacts, in this
order. **First match wins**; the precedence is part of the specification, because more than
one condition can hold.

| # | condition | status | reason |
|---|---|---|---|
| 1 | the arm was requested and not run — checkout, interpreter or dependency absent | `UNAVAILABLE` | `arm_absent` |
| 2 | `repository_state_before.source_tree_sha256 != repository_state_after...`; or the manifest's component states are known-blind (any run before `e083333`); or the run directory falls in a recorded contention window | `PROVENANCE_UNRESOLVED` | the matching reason |
| 3 | the run's endpoint gate FAILED on informativeness for this landscape and endpoint | `INVALID_PROTOCOL` | `endpoint_not_informative` |
| 4 | **no adaptive round executed** on any seed (every seed failed at or before round 1), judged over `round_index >= 1` only | `INVALID_PROTOCOL` | `no_adaptive_round_executed` |
| 5 | any seed recorded `degenerate_proposals = true`, or the budget classifier returns `DEGENERATE` | `INVALID_NUMERICAL` | `degenerate_proposals` / `numerical_collapse` |
| 6 | distinct measurements fall below `DEGENERATE_RUN_BUDGET_FRACTION · budget` (0.9) on any seed **and** the rejection histogram is dominated by `not_in_eligible_universe` | `INVALID_PROTOCOL` | `ineligible_proposal_starvation` |
| 7 | the arm's own second-phase policy produced nothing and the round ran the baseline instead, on a material share of seeds | `INVALID_PROTOCOL` | `policy_degraded_to_baseline` |
| 8 | distinct measurements < budget on ≥ 1 seed, shortfall classified `RECOVERABLE` | `VALID_BUDGET_LIMITED` | `eligibility_shortfall` |
| 9 | the pre-flight gate's report is absent from the run | `PROVENANCE_UNRESOLVED` | `component_state_unrecorded` |
| — | otherwise | `VALID_COMPLETE` | — |

**Two rules that make this computable at all:**

- **The status is arm-level and the evidence is seed-level**, so the aggregation is
  specified: conditions 4–8 are evaluated over all seeds of the arm, and **one bad seed is
  enough**. An arm is not rescued by its good seeds, because a mean over a mixture of
  statuses is the error §5.8 exists to prevent.
- **The status MUST NOT be derived from `round_metrics.csv` rows with
  `round_index == 0`.** Round 0 is appended at `benchmark_panoptes_isolated.py:3558` with
  `status="completed"` **before the selector is ever consulted**, so round 0 says
  `completed` for an arm that never ran. This is the §5.3 instance 4 defect, in one line.

### 5.3 The four real instances, classified

| instance | status | reason | discriminator in the artefacts |
|---|---|---|---|
| **An arm at 26 of 125 reporting `Run state: SUCCESS`** — and there are **two different arms at that number, with two different statuses** | see below | see below | `distinct_proposed` and the rejection histogram |
| ↳ HER2 `alde`: 25 shared starting variants plus **one variant measured 100 times** | `INVALID_NUMERICAL` | `degenerate_proposals` | `distinct_proposed = 1`; rejections dominated by `already_observed`; constant one-hot columns 55/200 → NaN → `argmax` returns 0 forever |
| ↳ CR9114 `panoptes_default`: **757,829 distinct candidates**, none above depth 4, 25,712 rejected / 32 accepted, **100 % `not_in_eligible_universe`** | `INVALID_PROTOCOL` | `ineligible_proposal_starvation` | `distinct_proposed` large; rejections 100 % `not_in_eligible_universe`; 192 of 192 rounds short |
| **eqFP611-red at 80.27 of 81 with 7 short-batch seeds** | `VALID_BUDGET_LIMITED` | `eligibility_shortfall` | shortfall classified `RECOVERABLE`; the arm searched |
| **HER2** — the endpoint gate fails on **informativeness** and **no budget fixes it**: the d ≥ 2 ceiling is 113 against a raw count of 113, i.e. **0.00 % of the hit set removed against a 5 % floor** | `INVALID_PROTOCOL` for **every** arm in that run | `endpoint_not_informative` | the gate's `endpoint` result, `detail["budget_saturates_tier"]`, and the equality 113 == 113 |
| **A live-LLM arm that failed all 12 seeds yet left round-0 rows in `round_metrics.csv` with `status=completed`** | `INVALID_PROTOCOL` | `no_adaptive_round_executed` | zero rows with `round_index >= 1`; `outcome.status = "failed"`; and the arm's `model_calls = 0` per round where the ID claims live calls |

Three notes that the classification is specifically designed to force:

1. **`VALID_BUDGET_LIMITED` is an admissibility verdict, not an explanation.** eqFP611-red
   spent 80.27 of 81 and found **0 hits against ridge's 15.35**. A 0.9 % budget shortfall
   does not explain a 15.35-hit gap, and the status **MUST NOT** be narrated as though it
   did. What it does mean is that the arm's mean **must not be compared directly against a
   fully budgeted method** (§5.8); the comparison is redone at a matched truncated budget.
2. **HER2's status is a property of the (landscape, endpoint) pair, not of any arm.** No
   arm run there is a registered test, and **no budget choice changes that** — the
   off-by-one in `generality_sweep.sh` (113 requested against §3C's registered 112) is real
   and should be corrected, and correcting it does **not** open the gate: verified at 112,
   it still FAILS. §3A's *"permanently untestable at this endpoint"* stands. The decision to
   express this as `INVALID_PROTOCOL` + `endpoint_not_informative` rather than as a seventh
   status is a naming choice, flagged in §10; the **discriminator** is frozen either way.
3. **ALDE's statuses carry a standing condition.** ALDE on the seven incompatible
   landscapes is `INVALID_NUMERICAL` (§5.8). Any numerical fix becomes a separately named
   **`ALDE_stabilized`** that **MUST** be validated against the authors' behaviour on
   compatible datasets: *a stabilised implementation is not automatically the published
   method.* And availability is not compatibility (rule 2.6) — `availability_probe` and
   `compatibility` stay separate manifest keys and are never collapsed.

### 5.4 Relationship to the existing run state

`manifest.run.state` is **run-level** with the vocabulary
`["SUCCESS", "PARTIAL", "FAILED", "ABORTED"]` (`:5053`). It is **not** an arm validity
status and **MUST NOT** be read as one: a run holding one degenerate arm and five clean
ones is a single `PARTIAL`, which is why a 26-of-125 arm could sit inside a run whose
headline state was read as success. Phase 0 adds the arm-level axis **beside** the run
state, and both are recorded.

### 5.5 What a conforming run records

**RECORDS** under `results.arms[<arm_id>].validity`:

```
validity.status              enum    the six of §5.1
validity.invalid_reason      enum    required unless VALID_COMPLETE
validity.decided_by          str     the rule number of §5.2 that fired, first-match
validity.evidence            {n_seeds, n_seeds_full_budget, min_distinct_measurements,
                              mean_distinct_measurements, declared_budget,
                              short_batch_rounds, degenerate_rounds,
                              rejection_histogram, adaptive_rounds_executed}
validity.budget_classification  enum {DEGENERATE, RECOVERABLE, UNCLASSIFIED, null}
validity.enters_primary_comparison  bool   true iff VALID_COMPLETE and §7.2 holds
```

**VERIFIED BY:** recomputing the whole decision from `round_metrics.csv` +
`selection_trace.csv` + the gate report, and asserting it reproduces
`validity.decided_by`. The classifier **MUST** reuse
`veritas.preflight.check_budget` / `budget_observations_from_run` rather than
reimplementing the shortfall split, which already prefers the run's own
`backfill_added` / `distinct_proposed` / `degenerate_proposals` columns and says
`UNCLASSIFIED` rather than guessing.

---

## 6. Manifest honesty and the configuration fingerprint

Rule 2.12: *"when a behaviour is gated in two places, a record that reads one of them is
worse than no record, because it will be trusted."* `component_states()` read only
`--ablate-components`, so `dpp` and `dkl` were recorded **enabled** in every run where
they were off, and a dpp-on and a dpp-off run of one library **hashed identically**. That
cost a real attribution: the HER2 null control moved by t = −8.29 and the cause could not
be read from the artefacts.

Requirements for everything Phase 0 introduces:

1. **Every new switch is recorded as its effective state**, read from the place that
   actually gates the behaviour — not from the CLI. `component_states()` (`:2518`) is the
   pattern: it reads `math_worker.DPP_DIVERSITY_ENABLED` and
   `math_worker.TWO_STAGE_GP_ENABLED` as well as the flag, and keeps `requested_spec`
   separately so "off by flag" and "off by default" stay distinguishable.
2. **The `configuration_fingerprint` MUST change** when any of `candidate_source.id`,
   `candidate_source.params`, `phase_schedule.id`, or `pools.assay_available_pool_rule`
   changes. **VERIFIED BY** extending the existing self-check — a 32-arm sweep yields 32
   distinct fingerprints with no outcome leakage (rule 2.9) — to cover the new axes.
3. **Fingerprints do not compare across `e083333`.** Quote the commit, not the
   fingerprint (rule 2.12).
4. **The pre-flight gate runs before any scored run** (rule 2.7), and its report is
   recorded under `manifest.preflight`. An arm in a run with no gate report is
   `PROVENANCE_UNRESOLVED` (§5.2 rule 9) — the gate's own four static checks would have
   refused the ALDE arm on seven of nine landscapes before a single measurement was spent.

---

## 7. Comparison validity — the rules this specification itself must obey

Each is in [`HANDOFF.md`](HANDOFF.md) §2 and each exists because it was broken here.

### 7.1 A verdict needs a paired test over seeds (rule 2.11)

**The seed is the unit of independence.** Before any arm is called better or worse:
compute the per-seed contrast `d_s = Y_s(arm) − Y_s(baseline)` on shared initial
observations, run the paired test over `{d_s}`, report **seeds better / worse / tied**,
state the interval, and resample **whole seeds** (B = 10,000) for the bootstrap. Never pool
within-seed sub-units.

`eqFP611-blue: 9.50 vs 8.79` was published as a **WIN**; the paired test on the same
numbers is **t = 0.77, better on 25 of 48 seeds** — a coin flip, and one of only two
generalization data points the project had. *"A difference of means with no test is not a
verdict."* The same rule killed the redundancy argument in `docs/545.md` §4 and inverted
two of five ablation signs. **It is the most repeated error in this project.**

**MUST:** every two numbers placed side by side in a Phase 1–4 artefact carry the test next
to them. **VERIFIED BY** `scripts/analysis/paired_stats.py` over the per-seed rows.

### 7.2 Budget equality before any yield comparison

Engine arms once measured **120–124** against baselines' **125**, and that biases *against*
our arms. A primary matched-budget comparison requires, for every arm in it:

- `validity.status == VALID_COMPLETE`;
- identical `declared_budget`, identical `initial_size`, identical `batch_size`,
  identical seed block;
- `min_distinct_measurements == mean_distinct_measurements == declared_budget` on every
  seed.

Any arm failing this is reported **outside** the primary comparison with its status
attached (§5), or re-compared at a common truncated budget with the truncation stated.

### 7.3 `random_forest` is required in any comparison that claims a win

Per (F2). **MUST:**

- `random_forest` is present, at the same budget and seeds, in any table that asserts an
  improvement;
- a **win claim** requires beating the per-landscape **better of {`ridge`,
  `random_forest`}** on the primary endpoint, under §7.1 and §7.2. Taking the maximum of
  two baselines makes the bar harder, not easier, so it needs no multiplicity correction;
- *"an improvement over random alone is not a result"* (prereg §12.2), and an improvement
  that the published fixed depth schedule A8 also achieves is not a result for the engine.

Whether the **pre-registered primary contrast** (H1's effect size, defined versus ridge)
also changes is a pre-registration amendment and is in §10.

### 7.4 Genotype canonicalisation (rule 2.1)

Mutation strings omit wild-type positions: `V183A:F184A:V227A` is a **four-site** genotype
whose fourth site is WT. Always expand to full library length, fill WT at omitted
positions, and **assert equal length**. Never zip, never truncate — this bug *silently
inverted a conclusion about hit diversity.*

**MUST:** reuse `veritas/calibration/genotypes.py` (`canonical_combo`, `hamming`, `depth`).
There are already five implementations and **a sixth MUST NOT be added** (§5.7).
Consolidation happens in Phase 6, after the scientific freeze, in a branch, with old and
new run against the full reference corpus — not during the frozen study unless a
correctness defect demands it.

### 7.5 The pre-flight gate runs before any scored run (rule 2.7)

```
.venv\Scripts\python.exe -X utf8 -m veritas.preflight --library <key> [--run-dir <dir>]
```

Non-zero exit means do not trust the run. `NOT-APPLICABLE` and `SKIPPED` are recorded
results with a stated reason, never silent passes, and `bool(CheckResult)` raises on
purpose. The eight non-CR9114 derived tables are acknowledged explicitly via
`--acknowledge-derived KEY=WHY` citing `docs/corpus-audit.md` — acknowledged, **not
silenced**.

### 7.6 Path resolution (rule 2.8)

Gitignored results, the dataset directory and external checkouts resolve against the
**real repository root**, not a worktree. Run `python -m veritas.paths` before concluding
anything is missing, use `veritas/paths.py`, and make every failure **name every path it
checked** — a bare "not found" is what let our strongest baseline get silently dropped
from a 48-seed run.

---

## 8. What a conforming run records — consolidated

One table, so a reader can check a manifest without reading the rest of this document.

| manifest path | §  | non-conformance if |
|---|---|---|
| `configuration.pools.*` | 1.5 | absent, or `assay_available_pool_size != oracle_evaluable_pool_size` under the retrospective rule |
| `configuration.candidate_source.*` | 2.8 | absent, or `params` disagrees with §2.3 without an amendment row |
| `configuration.phase_schedule.*` | 3.8 | absent, or `transition_inputs != ["measurements_at_round_start"]` for `phase_schedule_v2` |
| `configuration.components.states` | 6 | derived from the CLI alone (rule 2.12) |
| `configuration.configuration_fingerprint` | 6 | unchanged across a change to any §6.2 axis |
| `configuration.components.qualify_minimum_set_floor.enabled` | 3.4 | `true` under `phase_schedule_v2` |
| `manifest.preflight` | 6 | absent for a scored run |
| `results.arms[*].validity.*` | 5.5 | absent, or `decided_by` not reproducible from the CSVs |
| `results.retrospective.qualification_quantile` | 3.5 | recorded as an online parameter instead |
| `round_metrics.csv` new columns | 2.8, 3.8 | absent, or phase columns vary by seed for a positional rule |
| `candidate_trace.csv` `source` | 2.8 | disagrees with the per-round source split |

---

## 9. How a reader verifies conformance after the fact

Artefact-only; no re-run required except where stated.

1. **Resolve paths first.** `python -m veritas.paths`. Results live in the gitignored
   `benchmark_comparison_results/`, invisible from a worktree (rule 2.8).
2. **Read the run you care about by its manifest, never by "the newest directory"** —
   during contention that was another sweep's, and TrpB was reported as CreiLOV for an
   hour. Run IDs are UTC; the shell is local.
3. **Gate report.** `manifest.preflight` present, `exit_code == 0`, and every
   `NOT-APPLICABLE` carrying its reason.
4. **Provenance.** `repository_state_before.source_tree_sha256 ==
   repository_state_after.source_tree_sha256`. Run directory not in a recorded void window.
   Commit quoted, not a fingerprint compared across `e083333`.
5. **Pools.** Recompute `design_space_size` from `LibrarySpec` and match §1.2. Check
   `pools.identical_across_arms`.
6. **Candidate source.** Recompute `locality_share` from the distance histogram. Re-run
   `build_frontier_pool_v2` from the recorded inputs and assert an identical shortlist.
7. **Phase schedule.** Every seed of the arm transitions at the same round. No
   `qualified before round N` note.
8. **Validity.** Recompute every arm's status from §5.2 and match `validity.decided_by`.
9. **Comparison.** For each claimed contrast: both arms `VALID_COMPLETE`; budgets equal on
   every seed; `random_forest` present; paired statistic, seeds better/worse/tied and the
   interval all reported.
10. **Anchors.** Reproduce the rule 2.2 shell figures before trusting any depth number.

---

## 10. Still open

Listed rather than answered. This project's recorded failure mode is **confident
wrongness, not indecision.**

1. **Whether the pre-registered primary contrast moves from ridge to the better of
   {ridge, `random_forest`}.** §7.3 requires RF in any *win claim*, which this document can
   fix. But H1's effect size is *defined* as "versus ridge" (prereg §6.4), and prereg O-10
   records "ridge is the reference baseline" as a deliberate choice. Changing it is a
   pre-registration amendment that alters a frozen hypothesis, and §5 is silent on
   `random_forest` entirely. **Needs a human.**
2. **Whether the secondary depth covariate's amendment is admissible at all**, given that
   the licence §5.4 relies on ("no outcomes have been inspected") is false as of
   2026-11-01 (§4.5). The options are: record it as a post-outcome amendment and label the
   six already-read-out landscapes exploratory; or reserve a new landscape under prereg
   §4.4. **Needs a human** — it is the same class of call as the publish/raise question in
   [`HANDOFF.md`](HANDOFF.md) §6.
3. **Whether `f` should be constant or phase-dependent.** §2.3 fixes `f = 0.20` for both
   phases because a constant is the only version that does not confound the reserve with
   the schedule. Whether phase 2 should carry a smaller reserve is a real question and no
   evidence here addresses it.
4. **Whether landscape-gate failure deserves its own status name.** §5.8's six are
   per-arm; HER2's failure is per-(landscape, endpoint). §5.3 expresses it as
   `INVALID_PROTOCOL` + `endpoint_not_informative`, which is distinguishable and
   checkable. Whether the vocabulary should grow a seventh member is a naming question for
   the reviewer, not a semantic one.
5. **The observed-progress switching margin and its forced checkpoint** (§3.6). §5.2 defers
   the whole arm until the fixed schedule is benchmarked, so fixing a margin now would be
   inventing a parameter for an arm that does not yet have a reason to exist.
6. **What `assay_available_pool` is for HER2 prospectively.** `design_space` is
   1.82×10¹¹ and 11,276 genotypes carry labels. §1.3 answers the retrospective case by
   disclosure. A prospective pool rule for a library that is 6×10⁻⁸ of its own design space
   has no answer in §5 and none is invented here.
7. **`C` and `L` are engineering starting points, not optima** — §5.1 says so for `C`
   explicitly. They **MUST NOT** be tuned on the final-test landscapes (prereg §4.4,
   [`HANDOFF.md`](HANDOFF.md) §7 Phase 4). How they *should* be chosen is unanswered.
8. **Whether HER2 can be given a usable endpoint at all.** §5.3 (§3A) stands:
   "permanently untestable at this endpoint." Whether a different endpoint —
   measurement-efficiency curves, or yield normalised against the eligible hit-set size,
   which §5.3 names as the fix for budget saturation — rescues it is unstudied, and
   **early stopping MUST NOT be used to paper over a badly defined metric** (§5.3).

---

## 11. Order of implementation, and what each item blocks

Candidate-space **correctness** comes before any performance comparison
([`HANDOFF.md`](HANDOFF.md) §7). Within Phase 0:

| # | item | §  | blocks |
|---|---|---|---|
| **0.1** | the three pool semantics, materialised and recorded | 1 | **everything.** `frontier_pool_v2` is defined over `assay_available_pool`; the depth covariates over `oracle_evaluable_pool`; the validity statuses read the rejection histogram that only exists because the two differ |
| **0.2** | arm validity statuses — §5.8 says *immediately* | 5 | re-representing (not deleting) the 2026-10-31 void results; Phase 4's primary comparison; any re-quote of a pre-Phase-1 generative number |
| **0.3** | policy ID registry and the legacy freeze | 2.7 | `frontier_pool_v2` and `phase_schedule_v2` landing without silently redefining an existing arm, which has already happened once |
| **0.4** | manifest + fingerprint extension, effective states only | 6, 8 | attribution in Phases 1, 2 and 3. Without it a two-component change hashes identically and the next t = −8.29 is again unattributable |
| **0.5** | `frontier_pool_v2` specification → Phase 1 implementation | 2 | Phase 1's required tests: saturated 4-site; sparse 10-site HER2-like; a **disconnected Hamming-1 graph**; mixed per-position alphabets; a small library with few unmeasured candidates left; duplicate lead neighbourhoods; exhaustion; deterministic tie-breaking; exact budget fulfilment when feasible; no hidden-outcome access. *Acceptance is correctness, not performance* |
| **0.6** | the fixed phase schedule → Phase 2 | 3 | Phase 2; and the retrospective relocation of the quantile. Independent of 0.5 in principle, but **MUST NOT** be benchmarked simultaneously with it — §5.3 forbids changing generation and qualification at once |
| **0.7** | depth covariates, with the §4.5 amendment record | 4 | the H1/H2 read-out in Phase 4, and the §10.1 absolute-versus-proportional determination that prereg §6.7.4 P2 is blocked on. Does **not** block Phase 1 |
| **0.8** | the conformance checker (`scripts/analysis/phase0_conformance.py`, to be written — `scripts/analysis/` is another agent's directory today) | 8, 9 | any Phase 4 claim. A rule nobody can check from the artefacts is not a rule, which is the whole argument of `docs/preflight.md` |

**Commit all of it before viewing any outcome** (§7). After Phase 0: Phase 1
candidate-space correctness, Phase 2 the schedule, Phase 3 frozen-state diagnostics,
Phase 4 the independent benchmark with **ridge, random forest, random**, other compatible
baselines, reliability, `frontier_legacy_v1` where applicable and `frontier_pool_v2`;
Phase 5 stopping as a separate pre-registered study; Phase 6 genotype consolidation.
**Do not tune the new policy on the final independent landscapes.**

**Explicitly out of scope** (§7): a new GP architecture, a foundation model, a new
interaction-probing algorithm, an automated policy-learning controller. And per (F1), out
of scope for a *different* reason: a representation change is no longer the indicated next
step — the featurisation is not the binding constraint, the training distribution the
engine selects for itself is, which is what §2 addresses.

---

## 12. Amendments

Appended, dated, with a reason. The text above is not edited in place. **Any amendment
written after an outcome of the arm it governs has been examined must say so in its first
sentence.**

| date | section | change | reason |
|---|---|---|---|
| 2026-11-01 | — | initial freeze: pool semantics, `frontier_pool_v2`, the fixed phase schedule, the depth covariates, the arm validity statuses, the policy ID registry and the conformance record | discharges Phase 0 of `HANDOFF.md` §7, outstanding and blocking everything downstream |

**No `frontier_pool_v2`, `frontier_pool_v2_stratified` or `phase_schedule_v2` arm had been
run, and no outcome of any of them existed, when this document was frozen.** Outcomes of
the **legacy** arms had been examined and are recorded in `HANDOFF.md` §1 — which is why
§4.5 exists and says so.
