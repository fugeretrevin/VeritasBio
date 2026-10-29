# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Three audiences, and the product is expected to serve all three rather than
specialising to one:

- **Academic lab scientists** — protein engineers and postdocs choosing which variants
  to physically make next under a tight materials and reagent budget. Least
  software-fluent of the three, and therefore the **design tiebreaker**: when a decision
  trades expertise for legibility, resolve it for this user, which keeps the console
  usable for the other two.
- **Biotech / pharma discovery teams** — multi-seat, with compliance obligations, LIMS
  integration and internal review.
- **Computational biologists** — able to drive the CLI directly, but wanting campaign
  bookkeeping, budget enforcement and the audit trail handled for them.

The shared situation is the bench-to-screen loop: a combinatorial mutant library is far
too large to screen, and the user must decide which small set of wells is worth making.

## Product Purpose

PanoptesBio ranks a combinatorial mutant library computationally before anything reaches
the bench, and returns the small micro-batch worth physically making. It exists because
foundation-model inference across combinatorial mutation libraries produces runaway
cloud bills (a 238-residue protein yields >1.1×10⁷ double and triple mutants), and
because open-ended agents have no cost predictability.

Success is a user who screens a fraction of a percent of a library, finds variants well
above chance, and can account for every dollar of GPU spend that got them there.

## Positioning

A **deterministic, fail-closed six-step active-learning cycle** in place of an
open-ended agent. The mechanism a neighbouring product could not truthfully copy is the
separation enforced in code:

- An LLM (Meta Llama 3.1, strict JSON mode only) **only** proposes biophysical search
  boundaries. It never does arithmetic, executes code, or selects tools.
- A Gaussian process (5-factor physicochemical featuriser, Matérn 5/2) selects the
  micro-batch by Expected Improvement per dollar. It never touches external APIs.
- A fail-closed Python gatekeeper plus an immutable SQLite ledger (integer micro-USD,
  immutability triggers, ACID locks) authorises every paid call. The Modal worker never
  decides spend or batch size.

## Operating Context

Work is organised into **campaigns**: one protein, one objective, one budget. A campaign
runs in rounds; each round is the same fixed six steps (ingest → propose → select →
authorise → score → record).

- Compute runs on Modal (NVIDIA A10G) with ESM-2 650M/150M, batched inside the container
  keep-alive window so spend goes to tensors rather than cold boots.
- The human reviews and overrules the planner's proposed residue clusters before
  anything runs; banned positions are refused by the validator for the whole campaign.
- Results leave the system as physical wet-lab artifacts: Twist/IDT synthesis order CSVs,
  Opentrons OT-2 pipetting scripts, plate maps, and SHA-256 audit certificates.
  Measured assay values are uploaded back and the surrogate re-fits around them.
- Benchling LIMS synchronisation and an FDA 21 CFR Part 11 Merkle hash chain exist in the
  codebase for the regulated audience.

## Capabilities and Constraints

**Confirmed, enforced in code (non-negotiable engineering invariants):**

- Fail-closed on parse error, unknown models, missing budgets, or non-canonical mutations.
- No paid execution without `authorize_spend(...) -> True`. The ledger is never bypassed.
- Refusals are written to the ledger with the same weight as approvals.
- Budget is stored as integer micro-USD; a call that would exceed quota is refused
  outright rather than trimmed.
- Stop conditions: budget exhaustion, round cap, adaptive regret stopping (max EI < 0.005
  for 2 rounds), fitness plateau (no >0.001 gain over 3 rounds after round 3).

**Terminology** (use these words, they are the product's own): campaign, round,
micro-batch, variant, residue cluster, ban list, authorisation, ledger, gatekeeper,
planner, surrogate, wet-lab handoff, quota.

**Undecided — do not invent:**

- Pricing, plan tiers, and what a quota costs. The hosted beta runs on PanoptesBio's own
  Modal and model credits with a per-account cap; no commercial model is set.
- Team / multi-seat behaviour in the hosted console (the compliance schema supports
  multi-tenancy; the console does not yet express it).
- Whether the Benchling and 21 CFR Part 11 layers are exposed in the hosted console or
  remain codebase capabilities.

## Brand Commitments

- Name is **PanoptesBio** (renamed from VeritasBio, 28 Oct 2026). Apache 2.0. The rename is
  applied in `web/prototype.html`, this file and DESIGN.md; the engine, README, deck and
  scripts still carry the old name and are pending.
- Built for the Meta Global AI Developer Hackathon; Llama 3.1 and ESM-2 are Meta
  foundation models, Modal Labs provides the GPU infrastructure. These attributions are
  factual and stay.
- **Binding visual constraint, volunteered:** the interface is deliberately
  **single-theme dark** — Outfit display / Inter body / Fira Code mono. The dark-only
  commitment is intentional and is not to be "fixed" with a light theme. The accent moved
  from cyan `#00f0ff` to blue `#4d9fff` on 28 Oct 2026 (chosen from six verified
  colourways); DESIGN.md is normative for the palette.
- Voice is plainspoken and specific, and states limits in its own copy rather than
  burying them (see Product Principles).

## Evidence on Hand

Real, reproducible, in-repo — these are the proof assets and they must not be
embellished:

- **avGFP DMS** (Sarkisyan et al. 2016, *Nature* / ProteinGym). `scripts/` +
  `scripts/benchmark_results.png`. In-sample ρ +0.9718, cost $0.0002 vs $0.0184 brute
  force. The README explicitly flags these correlations as **fit, not predictive skill**.
- **TrpB active site** (Johnston et al. 2024), the headline benchmark. Combinatorially
  complete 20⁴ = 160,000 library at residues 183/184/227/228 of Tm9D8*. 125 of 159,999
  screened (0.078%), 36 top-1% hits vs 1.17 random = **30.9× enrichment**, best fitness
  0.8005 against parent 0.4081, cost $0.00175. **Held-out Spearman is −0.0365** and is
  published alongside the in-sample +0.9944. `scripts/download_trpb_benchmark.py`,
  `scripts/run_trpb_benchmark.py`, `tests/test_trpb_benchmark.py`.
- **GB1 binding domain** (Olson et al. 2014), 149,361 variants, 4-site epistasis.
- **Trastuzumab CDR-H3** (Mason et al. 2021, *Nat. Biomed. Eng.*), 36,391 variants, with
  developability liability scoring.
- 148/148 tests passing. Data provenance is asserted in tests: the TrpB downloader
  refuses to emit a benchmark whose recovered saturated positions disagree with the
  published active site.

**Absences future work must not fabricate:** there are no customers, no testimonials, no
case studies, no press, no pricing, no uptime or scale claims, and no deployed production
history. The console prototype at `web/prototype.html` runs on example data and calls
nothing.

## Product Principles

1. **Fail closed, and show the refusal.** A refused call is product evidence, not an
   error to hide. The ledger records refusals with the same weight as approvals, and the
   interface shows them.
2. **The model proposes; it never decides.** What to spend, what to run and when to stop
   are separate components. The user can see each one's decision, and the separation is a
   promise to the user, not an implementation detail.
3. **The user's constraints outrank the planner's.** Positions the user bans are refused
   by the validator for the whole campaign, and the planner cannot remove them.
4. **Report the weak number next to the strong one.** A lab deciding what to physically
   make deserves to know which claim is load-bearing. The TrpB held-out ρ of −0.037 ships
   next to the 30.9× enrichment, in the README and in the interface, deliberately.
5. **Spend is always visible and bounded.** Every paid call is priced before it runs and
   charged against a hard cap. Nothing queues, nothing overdraws.

## Accessibility & Inclusion

**No formal standard set.** This is explicitly undecided, not an oversight — future work
should not assert a conformance target that nobody committed to.

Recorded as current state rather than as a requirement: as of this writing the console
prototype clears WCAG 2.1 AA contrast across all measured text styles, carries no text
below 11px, is keyboard-operable (campaign rows, tabs with roving tabindex and arrow
keys), and honours `prefers-reduced-motion`. If a conformance target is adopted later,
AA is the bar the existing work already meets.
