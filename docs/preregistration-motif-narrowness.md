# Pre-registration — tail-motif narrowness predicts the sign of the diversity term

**Registered 2026-10-31, before any arm is run on any landscape other than TrpB.**

This document fixes a prediction, its measurement, its falsification conditions and its
analysis **before** the arms that test it. It exists because three of this project's
conclusions have already been re-narrated after the fact, and because the sibling
pre-registration (`docs/preregistration-cross-protein.md`) exists for the same reason.

---

## 1. The observation this explains

On TrpB, three independent facts point the same way:

| fact | source |
| --- | --- |
| All five variants ≥0.90 carry the same pair `K227 + G228`, with tolerant flanks | `docs/PIPELINE_FOR_BIO_REVIEW.md` §5.2 |
| Ridge reaches the ≥0.90 tier in ~17% of campaigns; the engine reaches it in **0 of 48** | §6.1 |
| Ablating k-DPP batch diversity **gains** +2.55 top-1% hits (48.15 vs 45.60, n=48 paired) | this round |

The candidate explanation is that **a diversity objective is wrong when the tail is one
narrow motif**: budget spent spreading is budget spent leaving the only productive
region. Ridge concentrates implicitly, which is why it reaches a tail the engine does
not.

That explanation is not yet a result. It was derived on one landscape, and the five
interventions in `docs/ENGINE_IMPROVEMENT_QUESTIONS.md` §2 are a standing reminder that
a mechanism derived on one landscape usually fails on the next.

---

## 2. The prediction

> **k-DPP batch diversity hurts top-1% yield where a landscape's tail is a narrow
> motif, and helps or is neutral where the tail is multi-modal.**

The quantity is measured by `scripts/analysis/tail_motif_narrowness.py`:
`concentration` is the mean, over library positions, of the most-common-residue
frequency among the tail variants, scored as a **z against a within-landscape
permutation null** that holds the genotype set — and therefore every marginal residue
frequency — fixed. The null is what makes a 4-site and a 16-site library comparable;
raw concentration is not.

### Assignment, fixed now

Two quantities matter, not one, and the quantile tail rule in §3 separated them:

* **z** — is the tail more concentrated than chance?
* **conserved positions** — is that concentration located at specific residues?

A tail can be concentrated without being a motif. That distinction did not appear
under a fixed tail of 5 and is the reason the rule in §3 was fixed before registering
the prediction rather than after.

| band | criterion | prediction for ablating k-DPP | landscapes |
| --- | --- | --- | --- |
| **narrow motif** | z ≥ 3.0 **and** ≥1 conserved position | **top-1% increases** | trpb, creilov, eqfp611_red, cr6261_h9, cr6261_h1 |
| concentrated, no motif | z ≥ 3.0 **and** 0 conserved positions | **no prediction** | gb1 |
| ambiguous | 1.0 < z < 3.0 | no prediction | eqfp611_blue, cr9114_h1 |
| **diffuse** | z ≤ 1.0 | **top-1% does not increase** | her2 |

Measured at registration, with the §3 rule (top 0.01%, floor 5, cap 50):

| landscape | n | sites | tail | bound | z | p | conserved | motif |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| trpb | 159,999 | 4 | 15 | quantile | **20.72** | .0005 | 2 | `227K:228G` |
| gb1 | 149,361 | 4 | 14 | quantile | **14.98** | .0005 | **0** | (none) |
| creilov | 165,427 | 15 | 16 | quantile | **11.26** | .0005 | 9 | 9 positions |
| eqfp611_red | 8,191 | 13 | 5 | floor | **5.13** | .0005 | 10 | 10 positions |
| cr6261_h9 | 1,841 | 11 | 5 | floor | **4.06** | .001 | 8 | 8 positions |
| cr6261_h1 | 1,886 | 11 | 5 | floor | **3.49** | .003 | 8 | 8 positions |
| eqfp611_blue | 8,191 | 13 | 5 | floor | 2.86 | .008 | 7 | 7 positions |
| cr9114_h1 | 65,093 | 16 | 6 | quantile | 2.59 | .010 | 8 | 8 positions |
| **her2** | 36,391 | 10 | 5 | floor | **−1.31** | .946 | **0** | **(none)** |

**Band membership is unchanged from the earlier fixed-count-of-5 measurement** —
the same six landscapes exceed z = 3, the same two are ambiguous, and HER2 remains the
only diffuse case. Only the magnitudes moved. The assignment is therefore robust to
the tail definition, which is the property it needs to have for the bands to mean
anything.

**HER2 is the dissociation case.** Its tail concentrates *below* chance and shares no
conserved position. A prediction that only said "diversity hurts" would be
unfalsifiable by a uniform result; this one requires the sign to **flip** on HER2. If
diversity hurts everywhere including HER2, the mechanism is wrong even though the
direction looks right, and the correct conclusion is that the diversity term is simply
miscalibrated rather than mismatched to landscape structure.

**GB1 is the reason the third band exists.** It is strongly concentrated (z = 14.98)
with **no** position conserved at 0.80 among its 14 tail variants: the tail is
non-random but not located at a shared motif. The mechanism under test is specifically
about a *narrow motif*, so GB1 carries no prediction. It is worth noting that GB1 is
independently reported at **parity at every budget** (`docs/HANDOFF.md` §1), which is
what a landscape carrying no prediction should look like — but that is a consistency
observation recorded now, not a hypothesis being tested.
## 3. Measurement decisions, fixed before running

- **Tail definition.** The tail is the top **0.01% of the eligible set, floored at 5 and
  capped at 50** variants. An absolute count was rejected: `--qualify-threshold 0.70`
  was already found to sit below the global minimum on four landscapes
  (`docs/HANDOFF.md` §3), and a fixed count has the same defect in the other direction —
  5 of 1,886 is a different object from 5 of 165,427. Where the floor binds, it is
  recorded in the result.
- **Null.** 2,000 permutations, seed 0, sampling the same tail size from the full
  genotype set.
- **Bands.** z ≥ 3.0 and z ≤ 1.0, as tabled above. Set now, not after seeing arm results.
- **Conserved position.** Most-common-residue frequency ≥ 0.80 among tail variants.
- **Genotypes.** `veritas.calibration.genotypes.canonical_combo`, per handoff rule 2.1.
  Mutation strings omit wild-type positions; combos are expanded to full library length
  and asserted equal-length. No sixth implementation.
- **Arms.** `panoptes_default` with and without `--ablate-components dpp`, identical
  seeds, 48 paired seeds, 25 initial + 4×25.
- **Statistic.** Paired t on top-1% hits per seed, |t| ≥ 2.01. Primary endpoint is
  top-1% yield; best fitness is secondary and reported either way.

---

## 3A. AMENDMENT 2026-10-31 — the falsification test is not currently runnable

Recorded the same day this document was registered, before any confirmatory arm was
accepted, because the defect is in the registration rather than in a result.

**Two constraints were missed at registration.**

**1. HER2's primary endpoint is degenerate.** `veritas.preflight --library her2`
reports:

> the d≥2 ceiling is 113 against a raw count of 113 (0.00% removed, floor 5%):
> DEGENERATE; a budget of 125 measurements is not smaller than the 113-variant hit
> tier: an arm can reach the whole tier, so the endpoint cannot discriminate

An arm can measure the entire hit tier, so top-1% yield cannot separate two methods
on HER2 at this budget. HER2 is the **only** diffuse landscape of the nine — every
other one scores z ≥ 2.09 — so §4's requirement that the sign **flip** on a
multi-modal landscape has no available test. Per handoff rule 2.7 a failed gate means
the run is not to be trusted, so no HER2 arm is run.

**2. Only three of the nine landscapes are runnable at all.** The benchmark's
`LIBRARIES` registry holds `trpb`, `gb1` and `her2`; the other six are registered in
`veritas.calibration.landscapes.REGISTRY` for diagnostics but are not wired into the
runner. The z values in §2 were measurable on all nine because the instrument reads
only a prepared table; the *arms* are not.

**What this leaves testable now:** `gb1` alone, which §2 assigns to the
concentrated-no-motif band and which therefore carries **no prediction**. Running it
tests the band assignment — a landscape carrying no prediction should show no strong
effect — but it cannot confirm or falsify the mechanism.

**The registration is therefore one-sided until one of three things happens**, and it
must not be reported as confirmed before then:

1. A non-degenerate endpoint is defined for HER2 and registered *before* its arm runs.
   Best-fitness is the obvious candidate and is already a secondary endpoint, but
   substituting it post hoc for the primary is exactly the move this document exists to
   prevent, so it requires its own amendment.
2. The six diagnostic-only landscapes are wired into the benchmark. Blockers: a
   `LibrarySpec` per landscape, and four structure files that are identified in
   `docs/structures.md` but not present locally (3M24, 8QI9, 5WL2, 5WKZ).
3. A further diffuse landscape is found. None of the nine qualifies, so this means a
   new dataset.

**No confirmation count is claimed in this document until a sign-flip test exists.**
The narrow-band evidence on its own is consistent with the diversity term simply being
miscalibrated everywhere, which §4 already names as the alternative the flip is there
to exclude.

---

## 3B. RESULT 2026-10-31 — GB1, and the band refinement it disconfirms

First arm run after registration. GB1, 48 paired seeds, `panoptes_default` with and
without `--ablate-components dpp`:

| metric | baseline | dpp off | Δ | t | verdict |
| --- | --- | --- | --- | --- | --- |
| best fitness | 5.5104 | 5.5224 | +0.012 | 0.09 | not separable |
| **top-1% hits** | 14.65 | **20.90** | **+6.25** | **3.92** | **significant** |

Better on 36 of 48 seeds. For comparison, TrpB gives +2.54 (t = 2.37), so **the effect
on GB1 is more than twice as large.**

### This disconfirms the three-band refinement, not the core prediction

§2 assigned GB1 to *concentrated, no motif* (z = 14.98, **zero** positions conserved at
0.80) and therefore gave it **no prediction**. The stated reason was that the mechanism
under test is about a *narrow motif* specifically. GB1 then behaved exactly like a
narrow landscape, and more strongly than the landscape the mechanism was derived on.

So the refinement is wrong as stated. The variable that tracks the effect is
**concentration (z)**, not whether concentration happens to sit at conserved residues:

| landscape | z | conserved positions | dpp off, Δ top-1% |
| --- | --- | --- | --- |
| trpb | 20.72 | 2 | +2.54 (t = 2.37) |
| gb1 | 14.98 | **0** | **+6.25 (t = 3.92)** |

The three-band scheme was introduced in this document *after* the quantile tail rule
exposed GB1's zero conserved positions — a post-hoc refinement to a pre-registered
prediction, and it is the part that failed. **The two-band version registered first
(narrow vs diffuse, on z alone) survives.** Recorded here rather than silently
reverted, because quietly dropping a band after seeing a result is the failure mode
this document exists to prevent.

Revised assignment, for the arms not yet run: **z ≥ 3.0 predicts that removing k-DPP
increases top-1% yield**, with no sub-condition on conserved positions. The diffuse
prediction (z ≤ 1.0) is unchanged and still has no available test (§3A).

### Two things this does not establish

- **It is not generalisation evidence.** GB1 is development data under handoff rule
  2.4, as are TrpB and HER2. Two development landscapes agreeing is mechanism evidence.
- **It does not distinguish the two explanations §4 cares about.** Both landscapes are
  concentrated, so "diversity is mismatched to tail structure" and "diversity is simply
  miscalibrated" both predict this result. Only a diffuse landscape separates them, and
  §3A records that none is currently runnable.

### One consequence worth flagging

`docs/HANDOFF.md` §1 reports GB1 at **parity at every budget** — the engine's weakest
headline result. Removing the diversity term there is worth +6.25 top-1% hits. That
parity may therefore be partly an artifact of a mis-signed component rather than a
property of the landscape.

---

## 4. What falsifies this

The prediction is wrong if **any** of the following holds:

1. k-DPP ablation does **not** increase top-1% yield on the narrow band — specifically if
   fewer than 4 of the 6 narrow landscapes move in the predicted direction.
2. k-DPP ablation **increases** top-1% yield on HER2, the multi-modal case. That would
   mean diversity is miscalibrated everywhere rather than mismatched to tail structure,
   which is a different and weaker claim.
3. The effect correlates with library width or dataset size rather than with z. Both are
   recorded so this can be checked, and width is already known to confound here —
   `docs/PIPELINE_FOR_BIO_REVIEW.md` withdrew a width claim once.

**TrpB does not count as evidence.** The prediction was derived from it, so it is
excluded from the confirmation count and reported separately. Of the nine landscapes,
**eight are untested** and all eight are blind at registration.

---

## 5. What it buys if it holds

A design rule rather than a leaderboard row:

> Measure tail-motif narrowness before a campaign — it needs **no labels beyond the
> landscape itself, no structure, and no baseline** — then set the diversity weight
> from it.

That is testable on every landscape with a prepared table, including the seven where
ALDE cannot function and the four final-test landscapes with no resolved structure. It
is also the first claim in this project that would generalise across proteins, as
distinct from the within-protein performance result, which `docs/HANDOFF.md` §8 predicts
will regress.

**It does not claim the engine is better than ridge.** It claims the engine's diversity
term has a sign that depends on a measurable property of the landscape, and that the
property can be measured in advance.

---

## 6. Scope and limits

- Nine landscapes, of which TrpB and GB1 are development data and HER2 was reclassified
  as development (handoff rule 2.4). So a confirmation here is **mechanism evidence, not
  generalisation evidence**, and must not be reported as the latter.
- `eqfp611_blue` and `eqfp611_red` share a source table and differ only in endpoint
  column, as do the CR6261 H1/H9 pair. They are not independent landscapes and are not
  counted as such.
- Five further landscapes in the corpus freeze (`hiv1_protease`, `rhla_rha`,
  `rhla_rls`, `rhla_pct_rha`, `sacas9`, `spcas9`, `eqfp611_combined`) have no prepared
  table yet and are therefore outside this registration. Adding them later requires
  recording their z **before** their arm runs.
- **The floor binds on four of nine landscapes.** At top 0.01%, a 1,886-row library
  yields 0.19 variants, so cr6261_h1, cr6261_h9, eqfp611_blue, eqfp611_red and her2
  all fall back to the floor of 5 — a tail of 0.27% rather than 0.01%. No single
  quantile spans an 88× range in library size without binding a bound, and the bound
  that applied is recorded per landscape rather than hidden. z remains internally
  valid for each landscape because the permutation null is drawn at that landscape's
  own tail size; what the floor costs is the claim that "narrow" means precisely the
  same fraction everywhere.
- Tail narrowness is a property of the *measured* table. On any landscape whose labels
  are partly imputed, the tail may be an artifact of imputation rather than of biology;
  TrpB is clear (`scripts/check_trpb_imputation.py`), the rest are unchecked.
