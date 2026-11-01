# Day 1 — Candidate reachability: is the bottleneck generation or ranking?

> Status: dataset-side analysis COMPLETE. Run-side confirmation (CR9114 48 seeds)
> IN FLIGHT at time of writing; §6 is the only section awaiting it.
>
> Scope note: this document answers a reviewer's challenge to an earlier claim of ours.
> Two of our own claims are retracted in §1. Nothing here is a product claim.

## 0. The question

> Is VeritasBio failing because it cannot **propose** the right sequences, or because it
> proposes them and **ranks** them poorly?

Everything below exists to answer that, and only that.

## 1. Two claims of ours that were wrong

**Retracted: "the top-1% hits on CR9114 start at depth 9."** They do not. The first
top-1% hit is at **depth 3**, and the **global optimum is at depth 7**
(`S31N:T57S:T71S:K74I:S75F:T76S:S77N`, score 9.8345). The claim came from reading
aggregated depth *shells* (`≤8` vs `9–16`) instead of the per-depth table. The
per-depth table is §3.

Independent cross-check: `benchmark_panoptes_isolated.py:516` already records that same
variant as the best observed in the prepared table. Our own spec disagreed with our
claim and we did not look.

**Retracted: "the generator's reachable depth climbs ~1 per round to about 6."** It does
not climb. It is a **hard ceiling at depth 4**, independent of budget, rounds, and
library size. Evidence in §2.

**Also softened: "search-space construction helps, ranking hurts."** Too neat. The
structural prior and the Pareto ranking are themselves *ranking* mechanisms, and they
are the two most load-bearing components we have measured. The real distinction may be
between biologically-informed ranking and poorly-calibrated exploration; we have not
established which, and the phrase should not be used until we have.

## 2. The depth ceiling is hard, and it is measured not inferred

`depth_returns.csv` records, per seed per round per depth, how many candidates were
`proposed`, how many `selected`, and how many were `top1_hits`. No new instrumentation
was required.

| run | library | positions | adaptive rounds | max depth ever proposed | top-1% hits found |
|---|---|---|---|---|---|
| `run_20261031_164431_b6c27377` | CreiLOV | 15 | **18** | **4** | 103 over 48 seeds |
| `run_20261031_165647_2804608c` | eqFP611-red | 13 | **11** | **4** | **0** over 48 seeds |
| `run_20261031_180556_27465f43` | TrpB | 4 | 4 | 4 (= saturation) | 49.83 top-1% mean |

CreiLOV ran **18 adaptive rounds** on a 15-position library and never once proposed a
depth-5 candidate. That is the whole argument for the ceiling being structural: it is
not a climb rate that ran out of rounds.

**Mechanism, from the code.** The ceiling is a **schema validator**, not a property of
the construct set. `controller/prompts.py:100`:

```python
max_mutation_depth: int = Field(default=2, ge=1, le=4,
                               description="Max simultaneous mutations (e.g. 2 for double mutants)")
```

Combined with `max_depth = max(strategy.max_mutation_depth, 3)` at
`workers/math_worker.py:1178`, the effective depth is **always in {3, 4}** — for every
library, every budget, every round count. Two places attempt to raise it and cannot:
`benchmark_panoptes_isolated.py:2369` requests `min(len(targetable), ...)` which is
already pre-clipped to 4, and `veritas/datasets/schema.py:290` computes
`effective_max_mutation_depth = n_sites` (16 for CR9114) which `controller/prompts.py:400`
then clamps with `max(1, min(strategy.max_mutation_depth, measurable))` — a `min`, so it
only ever clips down. The policy-gradient path carries the same cap independently at
`workers/math_worker.py:2024`.

An earlier draft of this document attributed the ceiling to the quadruplet crossover at
`workers/math_worker.py:941`. That was wrong, and worse: **that block is dead code.**
Lines 901–918 already accept `2 <= len(combined) <= max_mutations` over the unbounded
lead list and share the same `seen` set, so every notation lines 941–961 could emit has
already been emitted. When `max_depth == 3` the block cannot run at all (it requires
`max_mutations >= 4`); when `max_depth == 4` lines 901–918 get there first. It should be
deleted.

Separately, the stacking path at `workers/math_worker.py:1197` begins:

```python
for lead in lead_variants[:8]:
    if len(lead.mutations) >= max_depth:
        continue
```

A deep lead is not merely unable to grow — it is skipped entirely. So a *measured*
depth-9 variant contributes nothing to candidate generation.

**This answers the reviewer's strongest objection.** It is correct that the initial
25-variant random sample almost always contains deep variants — for CR9114,
P(at least one initial variant at depth ≥ 9) = **100.000%**, and P(at least one initial
variant is itself top-1%) = **22.2%**. Deep variants are always measured. They cannot
become productive leads, because of the two lines above. The ceiling lives in the
operator, not in the starting state, and CreiLOV's 18 rounds demonstrate it empirically.

A second axis — truncation in *position* space (`available_pos[:6]`,
`combinations(available_pos[:4], 2)`, `lead_variants[:16]`) and ordering bias from the
`max_candidates=2500` early return — is under separate audit and is **not** included in
the numbers below. The depth-4 figures here are therefore an *upper bound* on
reachability; the true reachable set may be smaller.

## 3. CR9114-H1, complete depth table

65,093 variants; top-1% = 651 variants at threshold 9.6088; 1,675 variants sit at the
assay floor (exactly 7.0, i.e. no detectable binding).

| depth | variants | top-1% hits | hit rate | enrichment | max score |
|---|---|---|---|---|---|
| 0 | — | — | — | — | — |
| 1 | 16 | 0 | 0.000 % | 0.00x | 9.3277 |
| 2 | 118 | 0 | 0.000 % | 0.00x | 9.4925 |
| 3 | 552 | 1 | 0.181 % | 0.18x | 9.6159 |
| **4** | 1,796 | 3 | 0.167 % | 0.17x | **9.6331 ← generator ceiling** |
| 5 | 4,302 | 10 | 0.232 % | 0.23x | 9.6804 |
| 6 | 7,915 | 39 | 0.493 % | 0.49x | 9.7512 |
| **7** | 11,339 | 60 | 0.529 % | 0.53x | **9.8345 ← global optimum** |
| 8 | 12,797 | 114 | 0.891 % | 0.89x | 9.7936 |
| 9 | 11,394 | 139 | 1.220 % | 1.22x | 9.7556 |
| 10 | 7,984 | 110 | 1.378 % | 1.38x | 9.7370 |
| 11 | 4,363 | 89 | 2.040 % | 2.04x | 9.6917 |
| 12 | 1,820 | 57 | 3.132 % | 3.13x | 9.7035 |
| 13 | 560 | 20 | 3.571 % | 3.57x | 9.6380 |
| 14 | 120 | 8 | 6.667 % | 6.67x | 9.6586 |
| 15 | 16 | 1 | 6.250 % | 6.25x | 9.6097 |
| 16 | 1 | 0 | 0.000 % | — | 9.5923 |

Enrichment rises **monotonically** with depth. The optimum is at depth 7 — *mid-library*,
not at either extreme.

### 3a. Why this kills reverse ablation as the primary fix

The depth-16 variant is the fully-mutated affinity-matured reference. It scores
**9.5923**, rank **1179 / 65,093 = top 1.81 %** — it is **not** in the top 1 %, and it
sits **9 mutations away** from the depth-7 optimum.

So a reverse ladder seeded from the matured antibody starts *outside* the hit set and
needs nine sequential single-mutation removals to reach the optimum. With four adaptive
rounds it reaches depth 12 at best. That is the mirror image of the forward ladder's
failure, exactly as the reviewer predicted, and worse than they assumed — they supposed
the matured reference might be too close to the answer. It is not close to it.

**Conclusion: direction is a red herring. The defect is sequential locality.** Any
operator that moves one mutation per round fails here regardless of which way it points.
Depth-stratified global proposal is the intervention that follows from the data.

## 4. Reachability across all nine landscapes

Top-1 % computed per landscape; **HER2 excludes the 25,114 rows whose `DMS_score` is
exactly 0.0** per rule 2.3 (absent labels, not measurements of zero) — this gives a
113-variant top-1 % tier, which matches `docs/dataset-inventory.md`. An earlier draft of
this analysis wrongly included them and reported 364.

| landscape | parent | positions | N | top-1 % | hits at depth ≤4 | % of hits reachable |
|---|---|---|---|---|---|---|
| **TrpB** | functional | **4** | 159,999 | 1,600 | 1,600 | **100 %** |
| **GB1** | functional | **4** | 149,361 | 1,494 | 1,494 | **100 %** |
| eqFP611-blue | functional | 13 | 8,191 | 82 | 51 | 62.2 % |
| CreiLOV | functional | 15 | 165,427 | 1,654 | 220 | 13.3 % |
| CR6261-H1 | non-functional | 11 | 1,886 | 19 | 2 | 10.5 % |
| HER2 | functional | 10 | 11,277 | 113 | 9 | 7.96 % |
| eqFP611-red | non-functional | 13 | 8,191 | 82 | 2 | 2.4 % |
| CR9114-H1 | non-functional | 16 | 65,093 | 651 | 4 | **0.61 %** |
| CR6261-H9 | non-functional | 11 | 1,841 | 18 | 0 | **0 %** |

GB1 is 100 % reachable whether or not its own 29,477 zero-score rows are excluded (no
rule mandates excluding them, and they do not reach the top 1 % either way), so that
ambiguity does not affect any conclusion here.

### 4a. The finding that reframes the project

**TrpB and GB1 are 4-site libraries.** Depth 4 is full saturation. The ceiling is
structurally invisible on exactly the two landscapes the engine was built, tuned,
ablated and published on.

Seven of nine landscapes are reachability-limited. Every documented "failure to
generalize" is a landscape where the generator cannot reach most of the answer.

The honest framing is therefore not "we overfit to TrpB." It is: **we built a search
operator with a radius of 4 and validated it only on libraries 4 wide.**

### 4b. The comparison that matters

On CR9114, confined to depth ≤4, the generator can reach 2,482 of 65,093 variants
(3.81 %) containing **4 of 651 top-1 % hits**. So the engine's *maximum achievable*
top-1 % recall on CR9114 is 4 hits — with a perfect surrogate, perfect acquisition and
unlimited rounds.

Uniform random sampling at 125 measurements expects `125 x 651/65093 = 1.25` hits, with
no ceiling. The generator is **handicapped relative to uniform random sampling before
any ranking happens.** That is the baseline comparison that matters and §6 tests it
directly.

## 5. Methodological positions adopted

Taken from reviewer critique and adopted:

- Test **operators** (forward / reverse / depth-stratified), not directions.
- **No automatic mode selection** from parent functionality. It is an explicit
  configuration with a documented rationale until independently validated. The rule was
  derived from these same nine landscapes and has no out-of-sample test.
- **No ESM pseudo-likelihood** as a functionality proxy. Sequence plausibility is not
  antigen-specific affinity; a germline antibody is a perfectly plausible sequence that
  does not bind H1.
- **No claim that GP uncertainty is intrinsically harmful.** The four-rung ladder
  (49.83 / 50.94 / 50.56 / 61.10) rules out mis-tuning via the exploration coefficient
  and nothing more. Calibration remains a separate, unrun experiment.
- **Baseline fairness.** No arm receives the matured reference unless every arm does.
  Handing one arm a near-optimal starting sequence is the same class of error as the
  answer key removed from the prior in Phase 2.
- **CR9114 is reclassified as development data** under rule 2.4, because we are designing
  an intervention against its observed failure. Any improvement there is mechanism
  evidence, not generalization evidence. A new landscape must be reserved.

## 6. Run-side confirmation — COMPLETE. Generation is the bottleneck.

`run_20261031_182558_8eb87d23`: 48 paired seeds on `cr9114_h1`, methods
`random,ridge,panoptes_default`, 4 rounds x 25. Run state PARTIAL.

### 6a. Verdict

**Candidate generation is the confirmed bottleneck. Acquisition is exonerated on this
landscape, because it never received an eligible candidate to rank.**

| method | depths proposed | top-1 % hits (48 seeds) | status |
|---|---|---|---|
| `panoptes_default` | **1–4 only** | **0** | **DEGENERATE — null arm, 26 of 125 measurements** |
| `random` | 1–15 | 42 | valid |
| `ridge` | 3–15 | **88** | valid |

`panoptes_default` built 757,829 candidates across the run, **none above depth 4**. Of
what its ranked lists offered to the eligibility judge, **25,712 were rejected and 32
accepted — 100 % of rejections were `not_in_eligible_universe`**, no other reason
appearing even once. 192 rounds committed short batches.

Prediction 1 (no proposal above depth 4) — **confirmed**, 757,829 proposals, zero.

Prediction 2 (`panoptes_default` underperforms `random`) — **cannot be evaluated as a
performance result.** The arm reached 26 of 125 measurements and the harness flags it
DEGENERATE. It is a null arm: it did not search the landscape. Reporting 0-vs-42 as a
performance comparison would be the exact error the degeneracy check exists to prevent.
The correct statement is that the arm **could not execute the protocol at all**.

### 6b. The dominant defect is the ALPHABET, not the depth ceiling

`generate_candidate_variants` enumerates a hardcoded 19-residue alphabet per position
(`workers/math_worker.py:1118`). CR9114 is a **binary** library: each of its 16 positions
admits exactly 2 residues, the germline and the matured one. So at depth 1 the generator
proposes ~19 x 16 = 304 single mutants where only **16** exist — 94.7 % nonexistent — and
the ineligibility compounds multiplicatively with depth. That is why 99.9 % of proposals
are not library members and why the arm starves.

`docs/dataset-inventory.md` already states this plainly: *"The new libraries are binary,
not saturated."* Four of the six newer landscapes carry exactly 2 residues per site.

| landscape | sites | alphabet per site |
|---|---|---|
| **TrpB** | **4** | **20 uniform (saturated)** |
| **GB1** | **4** | **20 uniform (saturated)** |
| HER2 | 10 | 20,20,19,19,18,13,12,10,9,5 |
| eqFP611 blue / red | 13 | **2 (binary)** |
| CreiLOV | 15 | 2 (x12), 3 (x2), 5 (x1) |
| CR9114-H1 | 16 | **2 (binary)** |
| CR6261-H1 / H9 | 11 | **2 (binary)** |

### 6c. Generation defects, ranked

Established by code audit plus a read-only in-process probe of the pure function on a
synthetic 16-position / 2-cluster strategy.

| # | defect | location | kind | effect |
|---|---|---|---|---|
| 1 | **Alphabet mismatch** — hardcoded 19-residue alphabet against 2–5 residues per site | `math_worker.py:1118` | parameter exists, unused | ~99.9 % ineligibility on binary libraries; **causes the degeneracy** |
| 2 | **Depth capped at 4** by schema validator | `prompts.py:100` (`le=4`) | hard | depths 5–16 unreachable; caps CR9114 recall at 4 of 651 |
| 3 | **Budget starvation** — stages demand 17,954 slots against `max_candidates=2500` | `math_worker.py:1106`, early returns 1191–1286 | default, never overridden | intra- and inter-cluster doubles get **zero** budget whenever any lead exists |
| 4 | **Only 8 leads** ever extracted, by inline literal | `math_worker.py:1808` | hard | depth > 2 exists only via leads, so 8 single leads yield **zero** depth-4 proposals |
| 5 | **Position slices** take the numerically lowest sites, not the best | `:1200` `[:6]`, `:1219` `[:4]`, `:1248` `[:4]`, `:1267-8` `[:3]` | hard | at depth ≥ 2 coverage collapses to **8 of 16** positions |
| 6 | **Alphabet prefixes are alphabetical** — `target_aas[:10]`, `target_aas[:8]` | `:1247`, `:1275`, `:1278` | hard | `M N Q R S T V W Y` can never appear in an intra-cluster double; 11 of 19 excluded from inter-cluster |
| 7 | **No WT reversion** — stacking only appends | `:1204-1207` and every `mut_aa != wt_aa` guard | hard | from a depth-3 lead, its depth-2 sub-genotypes are unreachable |

Measured from the probe (16 positions, `max_mutation_depth=4`, cap 2500):

| lead state | total emitted | depth histogram |
|---|---|---|
| no leads (round 1) | 1,975 | `{1: 288, 2: 1687}` — **depth ≥ 3 unreachable without leads** |
| 8 single leads | 2,500 | `{1: 288, 2: 218, 3: 1994}` — **zero depth 4** |
| 8 double leads | 2,500 | `{1: 288, 3: 216, 4: 1996}` — **zero depth 2** |

At the real default the function emits exactly 2,500 variants: 12.7 % of its own uncapped
reach and **3.8 % of the 65,536 genotypes of the library**. On a binary library roughly
**1,360 of those 2,500 are orderable at all**, and the orderable majority are
recombinations of already-measured leads rather than exploration.

Defect 7 is notable because the correct implementation already exists:
`generate_hamming1_lead_neighbors` (`math_worker.py:964-1042`) handles reversion properly
via `current_aa = lead_mut_map.get(pos, wt_aa)`. `generate_candidate_variants` never calls
it.

### 6d. The shape the generator assumes

Both of the generator's built-in assumptions — **few sites** and **full 20-residue
alphabet** — hold on TrpB and GB1 and on no other landscape in the corpus. The two
development landscapes are the only two of the generator's native shape.

This is the same pattern as §4a, now on a second axis, and it is the finding of Day 1:
the engine is not overfitted to TrpB's *biology*, it is hardcoded to TrpB's *library
geometry*.

### 6e. Baseline to beat

`ridge` found 88 top-1 % hits over 48 seeds (~1.83/seed) against `random`'s 42
(~0.875/seed), both searching depths 3–15 freely. Ridge is the target, not random.

## 7. Proposed intervention and its control

**Library-aware candidate generation.** `generate_candidate_variants` must take the
library's actual geometry as input instead of assuming TrpB's:

1. **Per-site allowed residues** from the `LibrarySpec`, replacing the hardcoded
   19-residue alphabet. This alone should lift CR9114 from a null arm to a running one,
   and it is a prerequisite for measuring anything else.
2. **Depth-stratified proposal** across the library's real depth range, replacing the
   radius-4 expansion from leads.
3. **Truncation limits** raised or removed as the §6c item 3 audit directs.

Surrogate, prior, Pareto ranking, acquisition, budget, initial samples and seeds all held
fixed. Items 1 and 2 should be separately switchable, so their effects can be attributed
rather than confounded — the Phase 2 ablation table was once invalidated by exactly that
mistake.

This is testable on **seven** landscapes with **two negative controls**:

> **Pre-registration.** Depth-stratified generation improves top-1 % recall on the seven
> reachability-limited landscapes and changes **nothing** on TrpB and GB1, which are
> already 100 % reachable and where a stratified sampler has nothing to add.

Three outcomes, all informative:

| outcome | conclusion |
|---|---|
| helps on the seven, flat on TrpB/GB1 | reachability confirmed as the mechanism, with controls |
| helps everywhere, TrpB/GB1 included | reachability is not the mechanism; the cause is elsewhere |
| helps nowhere | generation is not the bottleneck; return to acquisition |

### 7a. Pre-registered generality sweep

Recorded before any result returned. 48 paired seeds per landscape, arms
`random,ridge,panoptes_default`, pre-flight gate run per landscape (rule 2.7). Budgets are
per-landscape from §4a, not a default: eqFP611's whole top-1 % tier is 82 variants and
CR6261's is 19, both under a 125-measurement budget, so running them at 125 would measure
endpoint saturation rather than selection.

| # | landscape | sites | alphabet | budget | parent | prediction |
|---|---|---|---|---|---|---|
| 1 | CreiLOV | 15 | **mixed** 2x12, 3x2, 5x1 | 125 | functional | **improves** — only mixed-alphabet library, tests the per-site path rather than a uniform binary one |
| 2 | eqFP611-red | 13 | binary | 81 | **non**-functional | **improves from zero** — prior 48-seed run found 0 hits on every seed |
| 3 | eqFP611-blue | 13 | binary | 81 | functional | **improves** — same geometry as red, so the parent axis varies with geometry fixed |
| 4 | **GB1** | **4** | **saturated** | 125 | functional | **NO CHANGE** — negative control |
| 5 | CR6261-H1 | 11 | binary | 19 | **non**-functional | direction only; **low power by construction** |
| 6 | CR6261-H9 | 11 | binary | 19 | **non**-functional | same protein as H1, so **not independent of it** |
| 7 | **HER2** | 10 | mixed | 113 | functional | **NO CHANGE** — null control: it is `pool_restricted`, so `candidate_pool` is supplied and `generate_candidate_variants` is never called. Nothing in this change set can reach it. |

Two controls, and they fail differently. **GB1 moving** means these are not geometry fixes
and the diagnosis in §4a/§6d is wrong. **HER2 moving** means something changed that should
have been unreachable, i.e. an unintended side effect, since its generator is never invoked.

Interpretation limits fixed in advance:

- **CR9114 is development data** (rule 2.4): the intervention was designed against its
  measured failure. Its improvement is mechanism evidence only.
- **CR6261-H1 and H9 are one protein against two antigens.** They are not independent of
  each other and neither is independent evidence of generality.
- **Raw hit counts are not comparable across landscapes.** The d>=2 ceiling varies 65-fold
  (10 to 647, §4a), so any cross-landscape statement must normalise by the ceiling.
- **Hit counts must be compared on the same basis.** The smoke-test figures in §9b are
  "including the initial sample"; the `ridge` and `random` figures in §6a are sums of
  `top1_hits` over adaptive rounds only. Those two are not the same quantity and were
  briefly conflated while reporting the smoke test.

## 8. Prior art: this is a bug, not a field limitation

Literature review, sources in §8c.

### 8a. For enumerable libraries the field scores the whole library, every round

- **MLDE** (Wu et al., *PNAS* 2019): "the models simulate the fitnesses of **all possible
  sequences** and rank the sequences by fitness." GB1: "149,361 out of a total 20^4 =
  160,000 variants."
- **ALDE** (Yang et al., *Nat Commun* 2025): "an acquisition function is applied to the
  trained model to **rank all sequences in the design space**"; UCB "evaluated across all
  discrete x in the design space."
- **CR9114 specifically** is already an ML benchmark (AbBiBench, arXiv:2506.04235), and
  the one active-learning paper run on it uses a pool-based loop that "scores remaining
  unlabelled variants by acquisition" — no locality restriction.

**CR9114 is 2^16 = 65,536 variants — smaller than the GB1 landscape (160,000) that this
entire literature scores exhaustively.** A depth-4 cap reaching 3.81 % of an enumerable
library is not a design choice any published method makes.

### 8b. Even the methods that *are* local do not cap depth

- **AntBO** (Khan et al., *Cell Reports Methods* 2023) uses a Hamming trust region on an
  11-residue CDRH3 — with **d_max = 30**, a maximum radius exceeding the sequence length,
  plus an ablation (AntBO NT) removing the constraint entirely.
- **Casmopolitan** (Wan et al., ICML 2021): trust-region size ranges over "0 and the
  **dimensionality of the problem**."
- **FolDE** (arXiv:2510.24053, 2025) expands from single mutants of measured hits — our
  architecture — and documents that real campaigns "routinely explore sequence space **six
  or more mutations away**, with some venturing as far as **29**."

The field's restriction is on **which positions** are mutable (3–5 sites, then enumerate
all 20^k, so depths 0…k are all candidates). Ours is the opposite: 16 positions mutable,
depth ≤ 4. No published method imposes our restriction.

### 8c. Consequence for the fix: enumerate, do not stratify

Eight of our nine landscapes are enumerable:

| landscape | genotype space | enumerable? |
|---|---|---|
| CR6261-H1 / H9 | 2,048 | yes |
| eqFP611 blue / red | 8,192 | yes |
| CR9114-H1 | 65,536 | yes |
| TrpB | 160,000 | yes |
| GB1 | 160,000 | yes |
| CreiLOV | 184,320 | yes |
| **HER2** | **1.8e11** | **no** |

The benchmark harness **already holds the eligible universe** as `view.index_of`.
`panoptes_default` ignores it, generates candidates from biological heuristics, and has
99.9 % of them rejected as `not_in_eligible_universe`. The fix for eight of nine
landscapes is therefore not a stratified proposer — it is to **rank the eligible universe
directly**, which is both the field standard and far simpler than anything proposed in §7.

A proposal mechanism is needed only for HER2, which is independently degenerate (113-variant
top-1 % tier against a 125 budget).

### 8d. What survives as a research finding

"Depth-stratified candidate proposal" has **no established name** in the literature and no
located implementation. More usefully:

- **Kirjner et al.** (ICLR 2024, arXiv:2307.00494) define the **mutational gap** — minimum
  mutations from the training set to the optimum — and state the premise directly: "prior
  methods often **constrain optimization to a small mutational radius**, but this
  **drastically limits the design space**." They use the gap to grade task difficulty.
- **Nobody cross-tabulates proposal mechanism against the depth of the optimum.** The
  ingredients are all published separately; the combination is not.
- **No search-strategy benchmark exists on the Phillips CR9114/CR6261 landscapes.** They
  appear in ML work only as *scoring* benchmarks, never as *search* benchmarks.

So the contribution is not "reverse ablation works on antibodies." It is: *the depth of a
library's optimum relative to a proposer's reach determines whether local-expansion search
works at all, measured across nine real landscapes* — a gap Kirjner et al. defined and did
not cross. The bug was preventing us from measuring it cleanly.

## 9. What was changed, and what it did

Five defects fixed. Each is separately switchable or is a no-op on saturated libraries,
so attribution is preserved. 671 tests pass.

| # | fix | location |
|---|---|---|
| 1 | `allowed_residues` parameter: per-site allele restriction honoured by every construct | `math_worker.generate_candidate_variants` |
| 2 | depth bound `le=4` → `le=64`; the real limit is applied downstream per library | `controller/prompts.py:100` |
| 3 | position slices scale with alphabet size, so a stage's CANDIDATE budget stays fixed | `math_worker._position_budget` |
| 4 | **per-site restriction in the policy-gradient sampler**, the majority candidate source | `math_worker.sample_policy_gradient_candidates` |
| 5 | planner no longer pre-clips depth to the literal 4 | `benchmark_panoptes_isolated.py` `_AnnotationDrivenPlanner` |

Plus: `allowed_aa_per_position` declared for the six non-saturated libraries, and a
`--production-candidate-source {generate,enumerate}` switch (§9c).

### 9a. Fix 4 was the one that mattered, and the audit missed it

The code audit was scoped to `generate_candidate_variants`, so it found the combinatorial
generator's defects and not the sampler's. The run's own `candidate_trace.csv` settled it:

| candidate source | share of the engine's real ballot |
|---|---|
| `gflownet` (policy-gradient sampler) | **6,506 of 7,415** |
| `combinatorial` (`generate_candidate_variants`) | 909 of 7,415 |

After fixes 1–3 and 5, **84.85 % of the mutations in the engine's ballot still used a
residue the library does not contain** — because the majority source was untouched. The
combinatorial generator in isolation was by then 99.4 % eligible (693 of 697 proposals in
the universe; the 4 misses are genotypes genuinely absent from the 65,093-of-65,536 table).

Lesson for the next audit: scope it to the *measured* dominant path, not the path that
looks like it should be dominant. `candidate_trace.csv` already recorded which was which.

### 9b. Effect on CR9114, 4-seed smoke test

| measure | before | after |
|---|---|---|
| unique measurements of 125 | **26** (DEGENERATE) | **125, 125, 125, 125** |
| rounds committing a short batch | 192 of 192 | **0** |
| top-1 % hits per seed (incl. start) | **0** | **7, 4, 5, 8** |
| best fitness reached | 9.53 | **9.83451 = the global optimum** |

Seed 3 reached `S31N:T57S:T71S:K74I:S75F:T76S:S77N` — the depth-7 global optimum of §3,
which the depth-4 ceiling made unreachable by construction.

Generator pool on CR9114 geometry, and position coverage at depth >= 2:

| lead state | pool before | pool after | positions at depth >= 2 |
|---|---|---|---|
| no leads | 37 | **136** (= 16 singles + all `C(16,2)` = 120 doubles) | 8/16 → **16/16** |
| 1 depth-6 lead | 49 | **191** (now reaches depth 7 and 8) | 8/16 → **16/16** |
| 8 single leads | 115 | **640** | 8/16 → **16/16** |

### 9c. The task-symmetry switch

`--production-candidate-source enumerate` hands the production arm the same untested
eligible pool that `RandomSelector` and `RidgeSelector` always receive, since both call
`view.untested(revealed)` unconditionally. The default stays `generate`, so no existing
number changes meaning. The arm records a note when the switch is used, because ranking an
enumerated pool is strictly easier than generating under constraints and a result obtained
that way does not validate generation.

`pool_restricted` was deliberately NOT flipped in the library specs: that would silently
redefine what the benchmark measures and invalidate comparison with every prior run.

### 9d. A test was pinning the defect in place

`tests/test_library_registry.py` asserted `one_hot_width == n_sites * 20` — it *required*
every library to be saturated, which is the assumption that caused the bug. Replaced with
the real invariant (width follows the declared alphabets) plus a new check that each site's
parent residue appears in its own declared allele set. Same pattern as the
`assertIn("GLY B", ...)` test that pinned the placeholder bug in Phase 1.

### 9e. Still outstanding

- **`max(strategy.max_mutation_depth, 3)`** (`math_worker.py:1178`) silently raises the
  schema default of 2 to 3. It only ever *raises* depth so it does not block reachability,
  but it means a caller asking for doubles gets triples. Left alone deliberately: changing
  it would move the TrpB and GB1 numbers and confound this measurement.
- **Dead code** at `math_worker.py:941-961`, the quadruplet crossover. Lines 901–918
  already emit every notation it could produce, over the unbounded lead list and sharing
  the same `seen` set; when `max_depth == 3` it cannot run at all. Should be deleted.
- **Budget starvation** (§6c item 3): stages demand ~17,954 slots against
  `max_candidates=2500`, and the early return fires inside pair-stacking on the *first*
  lead, so the intra- and inter-cluster stages get zero budget whenever any lead exists.
  Not yet addressed.
- **No WT reversion** (§6c item 7): `generate_hamming1_lead_neighbors` implements it
  correctly at `math_worker.py:964-1042` and is never called from the generator.
- **48-seed CR9114 comparison against `random` and `ridge`, plus a TrpB regression
  check**, both in flight at time of writing. The TrpB arm should be bit-identical:
  it declares no per-site alphabet, so `allowed_residues` is `None`, `_position_budget`
  returns its base, and its four annotated sites make the planner change a no-op.

## 10. Data-provenance review: are the prepared tables safe for a DEPTH claim?

Running the pre-flight gate (rule 2.7) surfaced a **pre-existing** failure that no prior
run in this session had triggered, because no prior run invoked the gate.

`corpus.derived_table` FAILED with **7 landscapes flagged** — every landscape in the
generality sweep except HER2. CR9114 is examined and *not* flagged, because it was already
re-prepared from the raw eLife table.

| landscape | derived table read | raw alternative at the pin |
|---|---|---|
| CR6261-H1 / H9 | `Data/DMS/Clean/bnAbs_CR6261_H*_clean.xlsx` | `Data/DMS/Raw/bnAbs_CR6261_raw.csv` |
| CreiLOV | `Data/DMS/Clean/CreiLOV_clean.xlsx` | `Data/DMS/Raw/CreiLOV_raw.xlsx` |
| eqFP611 blue / red / combined | `Data/DMS/Clean/eqFP611_clean.xlsx` | `Data/DMS/Raw/eqFP611_raw.xlsx` |
| gb1_combingym | `Data/DMS/Clean/GB1_clean.xlsx` | `Data/DMS/Raw/GB1_raw.xlsx` |

**Why this is the right failure to care about here.** The motivating defect
(`docs/cr9114-reprep.md`) is that CombinGym's clean CR9114 workbook dropped 16,384
genotypes with **depth-correlated attrition that reached 100 % at depth 16**, removing the
matured antibody the landscape is about. Every claim in this document is about *depth*
reachability. If coverage falls with depth, an engine proposing deep variants is penalised
by the TABLE and the result is an artifact.

### 10a. Measured: observed vs expected genotype coverage, per depth

Expected counts come from the design's generating function — a site with `k` alleles
contributes `(1 + (k-1)x)`, so the `x^d` coefficient is the number of depth-`d` genotypes
in the complete design. Depth 0 is excluded: it is the wild type, whose mutation string
would be empty, and every one of these tables omits that row by construction.

| landscape | shallow-half coverage | deep-half coverage | worst single depth | verdict |
|---|---|---|---|---|
| eqFP611 red | **100 %** | **100 %** | 100 % | complete |
| eqFP611 blue | **100 %** | **100 %** | 100 % | complete |
| TrpB | **100 %** | **100 %** | 100 % | complete |
| CR9114-H1 | 99.1 % | 99.7 % | 98.3 % | flat |
| GB1 | 96.7 % | 93.3 % | 93.0 % | flat |
| CR6261-H1 | 92.4 % | 91.9 % | 81.8 % | flat |
| CR6261-H9 | 89.6 % | **90.2 %** | 87.3 % | inverse — deep half higher |
| **CreiLOV** | 92.6 % | 88.1 % | 62.5 % | **4.5-point gradient — caveat** |
| HER2 | — | — | — | sparse sample of 1.8e11; coverage not meaningful |

**No depth is at zero on any landscape.** None shows the CR9114 pattern. CreiLOV has the
only real gradient, 4.5 points, which sits just inside the 5-point threshold the check
script uses — it is carried as a caveat rather than called clean, and any CreiLOV-specific
depth claim should state it.

On that basis the seven flagged landscapes are acknowledged via
`--acknowledge-derived <key>=<why>`, which records the reason in the check's own result.
The gate then returns **4 passed, 0 FAILED, exit 0**.

### 10b. A sweep was discarded rather than reported

The first generality sweep was killed after one landscape. It had started CreiLOV *before*
`_one_hot_matrix` was corrected, so landscape 1 would have carried the 20-wide-per-site
encoding and landscapes 2–9 the declared-alphabet one. A comparison set assembled from two
code states is not a comparison set. The sweep was relaunched with all nine landscapes on
one code state, CR9114 included so its `ridge` arm gets the corrected 32-column encoding
instead of 320 columns of which 288 were identically zero.

## 11. Void runs: 2026-10-31, 16:30–18:55 local

**Every run directory created in this window is void and must not be quoted.** Ten of
them: `run_20261031_203023`, `_203648`, `_204209`, `_204522`, `_204603`, `_205330`,
`_210252`, `_210417`, `_211552`, `_225000`.

Cause: **two generality sweeps executed concurrently.** The first was stopped with the
harness's task-stop, which ended the *task* but left the underlying shell alive; it
continued launching landscapes for another three hours while a second sweep ran. Three
separate defects follow from that, each sufficient on its own to void the results:

1. **CPU contention.** Two 48-seed sweeps shared the machine, so every `round_seconds` and
   `cumulative_seconds` figure is inflated by an unrecorded amount.
2. **Mixed code state.** The stale sweep began before `_one_hot_matrix` was corrected, and
   each landscape is a fresh process, so earlier landscapes carry the 20-wide-per-site
   encoding and later ones the declared-alphabet encoding. One sweep, two encodings.
3. **Source tree changed mid-run.** `repository_state_before` vs `..._after` differ on
   `source_tree_sha256` for at least `_204522`, because this document was being edited
   while the run executed. The harness flags this itself and exits non-zero. It is right
   to: afterwards nobody can say which code produced the numbers.

The stale sweep also ran without the `--acknowledge-derived` arguments, so its pre-flight
gate failed on `corpus.derived_table` (§10) and its landscapes were never gated.

### 11a. Two reporting errors this produced

- Progress was read from "the newest run directory", which during contention was the
  *other* sweep's. `run_20261031_211552` was reported as CreiLOV at seed 2; it is
  **TrpB** — the diagnostics line naming 159,999 eligible variants and "depth 4 of 4"
  is TrpB's universe, not CreiLOV's 165,427. The "~2.6 min per seed, ~2 hours for
  CreiLOV" estimate derived from it is withdrawn.
- Run IDs are stamped in **UTC** while the shell reports local time (UTC-4 here), so a
  directory named `run_20261031_211552` was created at 17:15 local. Comparing the two
  made runs appear to start before runs that preceded them.

### 11b. Process rules adopted

- **One sweep at a time, enforced in the script**, which now refuses to start if any
  benchmark process is already running rather than waiting for one.
- **Stopping a sweep requires a tree kill** (`taskkill /F /T` on the shell PID, then
  verifying zero `benchmark_panoptes` processes remain). Stopping the task is not enough.
- **No repository edits while a sweep runs** — including to this file. The source hash
  covers the tree, so a markdown edit voids every landscape still to run.

## 12. Reproducing this document

```
python -m veritas.paths                      # resolve dataset and results dirs first
python -I scripts/analysis/hit_depth_reachability.py          # §4 reachability table
python -I scripts/analysis/cr9114_depth_facts.py       # §3 depth table, optimum, §4b arithmetic
python -I scripts/analysis/depth_profile.py \
    benchmark_comparison_results/run_20261031_164431_b6c27377=CreiLOV \
    benchmark_comparison_results/run_20261031_165647_2804608c=eqFP611-red
```

Depth is the count of substitutions in the mutation string. That is valid for these
tables *only* because wild-type positions are omitted (rule 2.1); `hit_depth.py` asserts
this by rejecting any synonymous `X###X` entry and by checking that no variant's depth
exceeds the library's position count. If either assertion fires, the depth definition is
wrong and every number here is void.
