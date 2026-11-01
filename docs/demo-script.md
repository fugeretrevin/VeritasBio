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
that tells you a campaign cannot succeed before you pay for it, and we found that out by
pointing it at ourselves.** Three of those five results are explained by ceilings this
tool measures. That is a stronger story than a win, and it cannot be falsified by the
repo, because it *is* the repo.

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
> fifty-five hits. Four and a half percent. A perfect model with an unlimited budget
> still finds seventy-five, because the answers are deeper than the plan. One plan, nine
> measured landscapes, sound on none of them.
>
> And the model is not being careless. Depth three is a reasonable exploration budget for
> a protein. How deep a library's answers sit is not in the sequence, not in the
> structure, not in the annotation — it is a property of the measured landscape. Only the
> gate knows it."

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
> Our engine scored zero on that landscape. We had it written down as a loss. It was not
> a loss. There was nothing it was permitted to select. Until last week no flag could
> even turn that filter off."

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
