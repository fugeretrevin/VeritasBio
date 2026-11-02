# Demo script — 3 minutes, recorded

Written for: whoever records and edits the submission video (you), not for a judge. The
spoken lines are in quotes and are meant to be read aloud or re-voiced; everything else
is direction.

**Deadline:** Nov 4, 9:00 PM ET. Judging Nov 4–8, winners Nov 9.
**Theme:** Push the Frontier.
**Judged on** (the only criteria the organiser has actually published): indie track —
*"creativity, technical ambition, and whether it works"*; startup track — *"technical
execution and whether it holds up in production."* No rubric weights are published.
Secondary write-ups circulate a 30/25/25/20 split (Impact / Meta tooling / Demo /
Technical depth); treat it as unconfirmed and do not design around it. What all versions
share is **"whether it works"**, so every claim in this video is a command running live.

**Target: 2:50.** Under three minutes with room to breathe. ~430 spoken words.

---

## The one decision that shapes the whole video

Do not claim the engine beats the baselines. It does not, on most landscapes, and the
measurements are in this repository: against one-hot ridge at 48 paired seeds it loses on
CreiLOV (−50.29, 0 of 48 seeds), loses on GB1 (−12.90), scores zero on eqFP611-red, ties
on eqFP611-blue (+0.73, t=0.77), and wins on CR9114 (+3.17) — which is development data
under rule 2.4. A judge who reads the repo finds that in ten minutes, and a video that
claimed otherwise would be over at that moment.

So the video claims the thing that is true and is *harder*: **we built the instrument
that tells you a campaign cannot succeed before you pay for it, and we pointed it at
ourselves.**

Be precise about how much that instrument explained, because we pre-registered a
prediction here and it was **falsified**. The hypothesis was "the spatial viability veto
is the cause of the three losses." Measured at 48 paired seeds with the veto ablated:

* **eqFP611-red** — the veto fully explains the **zero** (ceiling 0 of 82; the arm scored
  exactly 0.0000 on all 48 seeds, which is arithmetic, not a measurement). Ablated it
  scores 4.21, better on 28 seeds and worse on **none**. But ridge gets 15.35, so
  **73% of that loss is not the veto.**
* **CreiLOV** — not the cause. Ablating the veto *costs* 2.04 hits (t=−2.33) and widens
  the gap to −52.33.
* **GB1** — not the cause. Ablating costs 4.27 hits (t=−3.22).

On two of three landscapes **the veto was helping**, by raising the candidate pool's hit
density 2.8× and 1.8×. That is the opposite of what we expected.

That is the story: a gate that refuses plans, turned on its own authors, killing our
favourite explanation for our own bad results. It is stronger than a win, and it cannot
be falsified by the repo, because it *is* the repo.

### The strongest available beat, pending one confirmation run

**Do not record this until the 48-seed confirmation lands.** Measured 2026-11-01 at **12
paired seeds** (branch `worktree-agent-a83565bcde59a81f3`, `bb925c9`), the CreiLOV loss is
root-caused, and the cause is a **fabricated training row**. `select_next_micro_batch`
injects one row per round for the wild type with the literal label `y = 0.0` whatever the
assay's scale, at variance `1e-4` against a real measurement's `1e-2`. CreiLOV's real
labels occupy **2.793–4.196**, so that row sits **2.793 raw units below the minimum — 11.62
label SDs, or 1.99× the entire label range — and 15.68 label SDs from the centre** of the
distribution, and every round's surrogate is fitted through it. (An earlier draft said
"11.6 units below the minimum, 8.3× the label range". Both were unit-mixing; the
15.68-SDs-from-centre figure is the correct one and is what the argument rests on.)

| ablation (12 paired seeds, ridge constant at 59.167) | engine | Δ vs default | t |
|---|---|---|---|
| default | 9.083 | — | — |
| `wt_anchor` | **38.833** | **+29.75** | **12.17**, 12/12 seeds |
| `wt_anchor,prior` | 43.083 | +34.00 | 12.44 |

> **The 50.167 three-way row was here and is now VOID — do not say it on camera.** Its
> third component, `physchem_features`, was deleted from the engine at the integration
> merge on 2026-11-01 in favour of the narrower `biophysical_features`, so the figure
> cannot be reproduced on `main`. The two rows above are unaffected: neither uses it.

**Confirmed at 48 paired seeds** on two independent branches with two different drivers:
**+31.771 raw** (t=+35.88, 48/0/0) and **+16.04 at d≥2** (t=+24.29, 48/0/0), accounting for
63.2% of CreiLOV's raw deficit and 86.1% of its diversity-adjusted one. The engine still
loses CreiLOV in every mode.

> **WITHDRAWN 2026-11-01: "TrpB and GB1 are the only two landscapes where this defect is
> invisible."** That was my sentence and the explanation behind it is falsified. The
> label-scale theory predicted that distance from zero controls the damage. It does not:
> **CR9114 violates the scale by 2.47 label-ranges — worse than CreiLOV's 1.99 — and
> measures a tight null** (+0.042, t=+0.09, against an MDE of 0.924). A second reframe, on
> how wrong the zero is about the point it pins, was tested and **also falsified**: CR9114
> ranks 2nd of nine on that metric and is null, while GB1 ranks 7th and carries the only
> separable gain from *repairing* the anchor. Do not narrate a scale story on camera.
>
> What survived is narrower and still worth saying: **the row is the dominant cause of the
> CreiLOV loss specifically**, and the two landscapes every quotable engine figure came
> from are the two where it happens not to bite. That is a statement about our luck in
> which landscapes we reported, not a theory of the defect.

**The repair beats the deletion, and ablation could not have shown it.** Anchoring at the
*measured* parent label rather than removing the row is worth **+6.50 raw on GB1**
(t=+3.70, 33/15/0) where removal is worth **nothing** (+0.542, p=0.687). So part of this is
information in the true label, not just removal of a false one.

Why it is still not in the shot list: a `measured`-anchor figure means the arm receives the
parent's true label **outside its measurement budget**, and on TrpB that label is *itself a
top-1% hit*. Until it is charged (125 → 1 + 124) the number cannot be quoted as a clean
win. Three of the nine landscapes also have **censored** parent labels at an assay floor,
where the repair would anchor at a known-wrong value. The honest on-camera version is the
ablation result above, which is confirmed, paired, and reproduced twice.

---

## Pre-flight, before you hit record

```powershell
# 1. Confirm the live planner is reachable (one real API call, ~5k tokens).
.venv\Scripts\python.exe -X utf8 -m veritas.check --library trpb --from-model --skip-viability

# 2. Warm the structure cache so the viability check is fast on camera.
.venv\Scripts\python.exe -X utf8 -m veritas.check --library eqfp611_red --depth 5
.venv\Scripts\python.exe -X utf8 -m veritas.check --library creilov --depth 3

# 3. Confirm the suite is green, and show the number on screen later.
.venv\Scripts\python.exe -X utf8 -m pytest tests/ -q
```

Terminal at ~100 columns, large font. If a live call fails on camera, **stop and re-run**
— do not cut to a cached result and narrate it as live. The whole thesis is that an
unverified claim is worthless.

**Append `2>$null` to every command you record.** Importing the Modal client prints two
`DeprecationWarning` lines about `asyncio.set_event_loop_policy` to stderr before the
tool prints anything, and they are not suppressible from inside the process —
`filterwarnings`, `catch_warnings` + `simplefilter`, and `python -W` were each tried and
each failed. Redirecting stderr also hides the tool's own error messages, so run each
command once *without* the redirect first to confirm it succeeded.

---

## Shot list

### 0:00 – 0:22 · The problem, in one concrete sentence

*On screen:* the plain terminal.

> "An AI agent can design a protein engineering campaign. It will give you a confident,
> biologically sensible plan: mutate these residues, go this deep. What it cannot tell
> you is whether that plan is capable of working — and neither can any metric anyone
> reports. We built the thing that can."

### 0:22 – 0:55 · The model proposes, live

*Run:*
```powershell
.venv\Scripts\python.exe -X utf8 -m veritas.check --library trpb --from-model
```

*Let the token count and `live:` provenance line land on screen. Point at them.*

> "That is a live call to Meta's model API. It read the target and proposed a plan. First
> check: can the assay even read the positions it asked for? Yes — all four. And that
> check passes because an earlier fix removed its cause. Before that fix, the model
> proposed distal allosteric residues and whole rounds were thrown away as not in the
> eligible universe."

*Pause on the `REACH` block.*

### 0:55 – 1:35 · The catch the model could not have avoided

> "Second check is the one that bites. The model asked for a mutation depth. Where do the
> good variants actually live?"

*Run:*
```powershell
.venv\Scripts\python.exe -X utf8 -m veritas.check --library creilov --depth 3 --all-landscapes
```

*Hold on the `ONE PLAN, EVERY LANDSCAPE` table.*

> "On CreiLOV, depth three reaches seventy-five of one thousand six hundred and
> fifty-five hits. Four and a half percent. Any proposer that honours that depth finds
> seventy-five at most, whatever its model and whatever its budget. One plan, nine
> measured landscapes, sound on none of them.
>
> And the model is not being careless. Depth three is a reasonable exploration budget for
> a protein. How deep a library's answers sit is not in the sequence, not in the
> structure, not in the annotation — it is a property of the measured landscape. Only the
> gate knows it."

> **CORRECTED 2026-11-01, and the correction is why this beat is now the WEAKEST of the
> three.** An earlier draft said "a perfect model with an unlimited budget still finds
> seventy-five, because the answers are deeper than the plan." That is true of a proposer
> that *honours* its declared depth and **false of our own engine**, which does not:
> `recombine_lead_variants` stacks measured leads and carries the candidate pool past the
> declared depth from round 1. Instrumented over the real selection path on all nine
> landscapes, the pre-veto ballot spans depths **1–12** on CreiLOV and reaches **99.9%**
> of its hit set; the depth ceiling costs **1 hit of 1,655** there and **0** on the other
> eight. The tool now prints that distinction itself.
>
> So REACH is a correct statement about a *plan* and the wrong thing to narrate as a
> prediction about our engine. Keep the beat only in the depth-honouring framing above —
> "any proposer that honours that depth" — or cut it for the stronger beat below.

### 1:35 – 2:15 · We pointed it at ourselves, and it found our own bug

*Run:*
```powershell
.venv\Scripts\python.exe -X utf8 -m veritas.check --library eqfp611_red --depth 5
```

*Hold on `VIABILITY  REFUSE` and the composed ceiling line.*

> "Third check is our own engine. It has a spatial viability filter that rejects
> co-mutations too far apart in 3D. Reasonable rule. On eqFP611-red it removes every
> single one of the eighty-two hits. Zero are attainable — before any ranking happens.
>
> Our engine scored exactly zero on that landscape, on all forty-eight seeds. We had it
> written down as a performance loss. So we turned the filter off and measured it. The
> zero was the filter — it went to four point two one, and not one seed got worse. But
> the baseline gets fifteen point three five, so it is still a loss; the filter was about
> a quarter of it. Our first guess was that the filter explained the whole thing. It
> didn't, and we only know that because we checked."

*Then, the integrity beat — scroll to a `CAVEAT` or `UNTRUSTED` structure line:*

> "And the gate audits its own input. It will tell you when the structure it is reasoning
> over is not the molecule on the plate — a wrong chain, a homolog that does not align —
> because a check that cannot see its own bad input is the exact failure it exists to
> catch."

### 2:15 – 2:40 · What skipping it cost, measured

*On screen: the receipt, as plain text.*

> "This is what one campaign cost us before these checks existed. Seven hundred and
> fifty-seven thousand candidate proposals. Twenty-five thousand rejections, every single
> one for the same reason. Twenty-six of one hundred and twenty-five measurements spent.
> Zero of six hundred and fifty-one hits found.
>
> Here is the part that should worry you: mean yield and best fitness looked unremarkable.
> Not broken — unremarkable. Those two numbers are identical in shape for a search that
> lost and a search that was never allowed to start."

### 2:40 – 2:50 · Close on what is real

*On screen: `pytest -q` green, then `evidence/`.*

> "Live Meta model API. ESM-2 on Modal GPUs. Every number in the repository has a
> committed run behind it, including the ones that contradict what we believed a week
> ago. The agent plans. The gate decides whether the plan is worth paying for."

---

## Lines to avoid

- "Our engine outperforms…" — it does not, on three of six landscapes. Measured.
- **"eqFP611-red was not a loss, it was the filter."** Half true and the half that is
  false is the important half. The *zero* was the filter; the *loss* was not — ablated,
  the arm still trails ridge by 11.15 at t=−5.76. An earlier draft of this script said
  "it was not a loss" and that was wrong.
- "The spatial filter was holding the engine back." On CreiLOV and GB1 it was **helping**,
  measurably. Do not generalise eqFP611-red to the other landscapes.
- "State of the art", "revolutionary", "breakthrough" — nothing here needs them, and one
  unearned superlative makes a judge re-read everything else sceptically.
- Any claim about the deployed ESM-2 path that is not additive-verified. The deployed
  Modal app was stale for days and silently returned scores from superseded code; it was
  redeployed and verified on 2026-11-01 (additivity residual exactly 0.000000).
- Antibody-target planner calls. A prompt carrying affinity evidence against influenza
  hemagglutinin trips the provider's content policy. That is correct behaviour on Meta's
  part and not something to design around or show being circumvented.

## If you have 30 seconds spare

Cut from 1:35 rather than from the close. The eqFP611-red beat and the self-audit beat
are the two that cannot be got anywhere else; the live planner call at 0:22 is the one a
judge will most expect to see, so it stays.
