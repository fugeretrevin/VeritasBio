# Next Steps

State as of the commits following `docs/545.md`. Every figure here is reproducible from
this repository; regenerate rather than cite. Where something is unverified, it says so.

**Baseline:** `benchmark_panoptes_isolated.py`, 48 paired seeds, identical 25-variant
start per seed, 4 rounds × 25, 125 unique measurements, paired t-tests (critical
|t| = 2.01). Default arm: **0.7450 best fitness / 44.94 top-1%**.

---

## 1. Closed since 545

| item | resolution |
|---|---|
| **Live end-to-end "never measured"** (545 §8) | It *was* measured — 48 campaigns, live Modal A10G `esm2_650m` (410.71 GPU-s), 384 live Muse calls (2,276,858 provider-reported tokens), audit chain valid 48/48. The fix had been applied in `.venv`, which is gitignored, so it never entered the repo. Now `controller/tls_trust.py`, opt-in via `VERITAS_USE_OS_TRUSTSTORE=1`, verified both ways against the live endpoint |
| **`spread=8.0` beating raw** (545 §8) | Artifact. At 48 seeds spread-8 gives 0.7484/45.92 against the corrected default's 0.7450/44.94 — the same arm within noise. The n=24 "win" compared spread-8 against a "raw" baseline that was itself silently calibrated at 0.25. No tuned-bound optimum exists |
| **Burned Part 11 signing key** (545 §8) | Revoked in code: `part11_secret_problem()` refuses that specific value, naming git history and giving the rotation command, and also refuses unset or <32-char secrets. Signing fails closed. **Does not rotate the deployed secret — that is still an operational action** |
| **Default arm ran the reverted prior config** | New. `ablation not in (...)` with `ablation=None` evaluated True, forcing calibration ON in every default run. The default command reproduced the arm 545 §2 measured as worse — 33.83 top-1% against the documented 45.60 |

---

## 2. Needs validation

### 2.1 The baseline figure itself
545 and the README quote **45.60** top-1%; a clean default run now gives **44.94**
(SD 5.27, so ~0.9 SE). Probably noise, but the two were measured on different trees and
one of those trees had the calibration defect. **Re-measure once and fix the quoted
number in both places.** ~2 minutes.

### 2.2 Artifacts that exist only locally
`benchmark_e2e_results/` and `benchmark_comparison_results/` are both gitignored, so the
48-campaign live run and every ablation folder vanish on a clean checkout. The numbers
are in commit messages and docs; the evidence is not. **Decide whether a summary
manifest per run should be committed** (the full folders are large, the manifests are
small).

### 2.3 The imputation argument
545 §1 rules the ≥0.90 tail measured rather than imputed, via a clean argument: a KNN
mean cannot exceed the maximum of its neighbours. Four of five are excluded outright;
`CIKG` survives only because a 76-way equidistant tie-break would have had to pick the
five highest neighbours (1 in 18,474,840). This is sound inference but it is **not the
upstream flags**. If the authors' supplementary data becomes available, confirm.

### 2.4 qualify-then-optimize is measured but unshipped
Built as `panoptes_qualify_then_optimize` and measured at 48 seeds:

| | baseline | qualify-then-optimize |
|---|---|---|
| P(≥0.50), P(≥0.70) at every budget | — | **identical** |
| seeds reaching ≥0.80 | 3/48 | **13/48** |
| seeds reaching ≥0.90 | **0/48** | **8/48** |
| best variant found | 0.8326 | **1.0000** |
| top-1% hits | 45.23 | 39.90 (−5.33, t = −3.20) |
| best fitness | 0.7457 | 0.7570 (t = 0.94, **not** significant) |

Run the baseline policy unchanged until a measured variant reaches 0.70, then switch the
remaining budget to greedy Hamming-1 search around the incumbent. Because every
trajectory is identical to the baseline until the protected event has occurred, the
reliability curves are preserved *by construction* — and they are, byte-identical at
every budget.

**It closes the tail gap to parity with ridge for a known, significant cost in
breadth.** Two caveats before shipping: the 0.70 threshold is favourable here
specifically because 88% of campaigns qualify by measurement 50, leaving three batches
for the frontier phase (on a slower landscape there would be no budget to switch into);
and the figures above were measured on the pre-fix tree, so they need one re-run.

**Recommendation: ship it as a second mode rather than a replacement**, and describe it
as *threshold attainment at the top improves, average peak does not* — the mean
best-fitness gain is not significant and loses 25 of 38 decided seeds.

---

## 3. Open problems, hardest last

### 3.1 Attributing the advantage — needs a joint design
Single-component ablations cannot do it. The advantage over random is 44.5 hits; the
five single removals cost **68.4 between them, ~154% of the total**, so the components
are redundant rather than additive. All five removals are significant and worse, which
establishes only the negative: the advantage is not concentrated in any one stage.
A 2^5 joint design is the honest instrument. 32 arms × 48 seeds ≈ 1.5k campaigns, which
is minutes at 2.3 s/seed for the non-LLM arms.

### 3.2 Candidate generation — "enumerate, don't rank"
The pool is singles and pairs; the ≥0.70 tier is 97% triples and quadruples. Verified
by attrition trace:

| leads supplied | pool depths | ≥0.70 in pool | ≥0.90 in pool |
|---|---|---|---|
| none | {1, 2} | 1/64 | 0/5 |
| top 8 of a random start | {1, 2, 4} | 1/64 | 0/5 |
| **+ the double `V227K:S228G`** | {1, 2, 3, 4} | **26/64** | **5/5** |

So the machinery is sound and lead-gated: once that one depth-2 variant is a lead,
lead-stacking generates all five tail variants. Its fitness is 0.6098 against a parent
of 0.408, it is in the candidate pool every round, and acquisition declines it because
its K227 component has a negative *marginal* effect.

**Six interventions have failed here.** Do not attempt a seventh without the §3.1
instrument. The pattern across all six is that acquisition re-imposes its marginal
preference on any slot that is reserved but still ranked; the one thing that worked
(§2.4) avoided the conflict rather than tuning through it. Specifically avoid:
widening the chemistry set; exempting active sites from the prior's charge penalties;
depth-stratified sampling (regressed via a base-rate error — P(good | depth) *peaks* at
depth 2–3); identity-coverage floors; and conditional scans that rank rather than
enumerate.

### 3.3 The prior expresses low mean where it should express high variance
`K227` has mean Δ −0.0206 over 8,000 matched backgrounds and is better than parent in
only 38.9% — so "K227 is bad" is correct on average. But its effect is conditional on
`G228` (53.2% positive across the 400 G228 backgrounds against 9.0% across the 399 S228
ones) and it has the largest upside of any residue at that position (max Δ +0.819).
The honest statement is "conditionally beneficial, condition unknown," which is
epistemic uncertainty over an interaction, not a low mean and not observation noise.
**No implementation is proposed.** Note that bounding the prior's magnitude is already
closed in the negative (545 §2) and deleting it costs 12.8 top-1% hits (545 §4).

### 3.4 Generalisation
One target, one four-site library. avGFP, GB1 and trastuzumab datasets are held but the
controlled protocol has never been run on them. **Nothing currently claimed generalises
beyond TrpB**, including the K227/G228 finding and the two-mode result.

### 3.5 Gatekeeper routing
545 §5 fixed self-priced spend on `/api/cloudlab/dispatch` and
`/api/campaign/ingest-curves`. 545 §8 still lists cloud-lab and curve-ingest as
spending outside the authorisation gate — confirm whether that line is stale or whether
a path remains.

### 3.6 Event-loop blocking
`ingest_plate_endpoint` and `cloudlab_webhook_endpoint` genuinely await request parsing
and so stay `async`, but their post-await SQLite work still runs on the event loop.
Explicitly **not fixed** in 545 §5.

---

## 4. Recommended order

1. **Re-measure and reconcile the baseline figure** (§2.1) — 2 minutes, and two
   documents currently quote a number from a tree with a known defect.
2. **Ship qualify-then-optimize as a second mode** (§2.4) — built, measured, and the
   only route to the tail we have. Re-run on the fixed tree first.
3. **Rotate the deployed Part 11 secret** — the code refuses the burned value now, but
   rotation is operational and still outstanding.
4. **Decide on committing run manifests** (§2.2) — otherwise the live-run evidence is
   one `git clean` from gone, which is how 545 came to record it as never measured.
5. **Build the 2^5 joint ablation** (§3.1) — the instrument every remaining engine
   question depends on.
6. Only then revisit candidate generation (§3.2).

Items 1–4 are hours. Item 5 is the one that unblocks the rest. Item 6 has a 0-for-6
record and should not be attempted before 5.
