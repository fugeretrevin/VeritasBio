---
name: PanoptesBio
description: Single-theme dark console for a budget-gated protein engineering platform.
colors:
  bg: "#070a12"
  bg-panel: "#0e1424"
  bg-card: "#151d33"
  bg-card-hi: "#1e2a48"
  border: "rgba(255,255,255,0.08)"
  border-hi: "rgba(255,255,255,0.16)"
  fg: "#f7f9fc"
  fg-muted: "#9aa8bf"
  fg-dim: "#8795ae"
  accent: "#4d9fff"
  accent-dim: "rgba(77,159,255,0.14)"
  emerald: "#34d399"
  emerald-dim: "rgba(52,211,153,0.12)"
  amber: "#fbbf24"
  amber-dim: "rgba(251,191,36,0.12)"
  rose: "#f43f5e"
  rose-hi: "#fb7185"
  rose-dim: "rgba(244,63,94,0.12)"
  violet: "#c084fc"
  accent-hi: "#7ab4ff"
  accent-ink: "#05101f"
  accent-edge: "rgba(77,159,255,0.32)"
  accent-ring: "rgba(77,159,255,0.10)"
  emerald-edge: "rgba(52,211,153,0.3)"
  amber-edge: "rgba(251,191,36,0.3)"
  rose-edge: "rgba(244,63,94,0.3)"
  bg-raised-hi: "#27365c"
  bg-sunk: "#04070e"
  thumb-hi: "#2e3e66"
  wash-faint: "rgba(255,255,255,0.02)"
  wash-soft: "rgba(255,255,255,0.04)"
  wash: "rgba(255,255,255,0.05)"
  track: "rgba(255,255,255,0.07)"
typography:
  display:
    fontFamily: "Outfit, ui-sans-serif, system-ui, sans-serif"
    fontSize: "clamp(2.375rem, 7vw, 4.125rem)"
    fontWeight: 800
    lineHeight: 1.2
    letterSpacing: "-0.03em"
  headline:
    fontFamily: "Outfit, ui-sans-serif, system-ui, sans-serif"
    fontSize: "clamp(1.5rem, 4vw, 2.125rem)"
    fontWeight: 700
    lineHeight: 1.2
    letterSpacing: "-0.02em"
  title:
    fontFamily: "Outfit, ui-sans-serif, system-ui, sans-serif"
    fontSize: "1.5rem"
    fontWeight: 700
    lineHeight: 1.2
    letterSpacing: "-0.02em"
  body:
    fontFamily: "Inter, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.9375rem"
    fontWeight: 400
    lineHeight: 1.6
    letterSpacing: "0em"
  label:
    fontFamily: "Outfit, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.6875rem"
    fontWeight: 700
    lineHeight: 1.4
    letterSpacing: "0.08em"
  code:
    fontFamily: "Fira Code, ui-monospace, SF Mono, Menlo, monospace"
    fontSize: "0.8125rem"
    fontWeight: 400
    lineHeight: 1.75
    letterSpacing: "0em"
rounded:
  xs: "3px"
  chip: "4px"
  sm: "6px"
  md: "9px"
  lg: "10px"
  pill: "999px"
  circle: "50%"
spacing:
  xs: "6px"
  sm: "10px"
  md: "14px"
  lg: "18px"
  xl: "22px"
  2xl: "26px"
  section: "56px"
  section-top: "76px"
components:
  button-primary:
    backgroundColor: "{colors.accent}"
    textColor: "{colors.accent-ink}"
    typography: "{typography.label}"
    rounded: "{rounded.sm}"
    padding: "10px 18px"
  button-primary-hover:
    backgroundColor: "{colors.accent-hi}"
    textColor: "{colors.accent-ink}"
  button-default:
    backgroundColor: "{colors.bg-card-hi}"
    textColor: "{colors.fg}"
    rounded: "{rounded.sm}"
    padding: "10px 18px"
  button-default-hover:
    backgroundColor: "{colors.bg-raised-hi}"
    textColor: "{colors.fg}"
  button-ghost:
    backgroundColor: "transparent"
    textColor: "{colors.fg-muted}"
    rounded: "{rounded.sm}"
    padding: "10px 18px"
  button-ghost-hover:
    backgroundColor: "{colors.accent-dim}"
    textColor: "{colors.accent}"
  card:
    backgroundColor: "{colors.bg-card}"
    textColor: "{colors.fg}"
    rounded: "{rounded.lg}"
    padding: "18px"
  input:
    backgroundColor: "{colors.bg}"
    textColor: "{colors.fg}"
    rounded: "{rounded.sm}"
    padding: "9px 12px"
  pill-running:
    backgroundColor: "{colors.accent-dim}"
    textColor: "{colors.accent}"
    typography: "{typography.label}"
    rounded: "{rounded.pill}"
    padding: "3px 10px"
  pill-complete:
    backgroundColor: "{colors.emerald-dim}"
    textColor: "{colors.emerald}"
    typography: "{typography.label}"
    rounded: "{rounded.pill}"
    padding: "3px 10px"
  pill-refused:
    backgroundColor: "{colors.rose-dim}"
    textColor: "{colors.rose-hi}"
    typography: "{typography.label}"
    rounded: "{rounded.pill}"
    padding: "3px 10px"
  nav-item-active:
    backgroundColor: "{colors.accent-dim}"
    textColor: "{colors.accent}"
    rounded: "{rounded.sm}"
    padding: "8px 10px"
---

# Design System: PanoptesBio

## Overview

**Creative North Star: "The Ledger in Plain Sight"**

The obvious thing is obvious. Whatever the user came to do — read a campaign's state, confirm which residues the planner may touch, decide what goes to the bench — sits in front, uncluttered, at full contrast. Nothing competes with it for attention. The surface is quiet enough that a single blue element reads as *the* live thing on the screen.

Underneath that clean surface, every decision is accounted for and nothing is hidden. The ledger shows refusals next to approvals. The benchmark shows the held-out correlation next to the enrichment. A round that was priced too high and rejected stays on screen as a card of its own. This is not decoration or humility theatre: it is the product's actual claim, so the visual system has to make the record inspectable rather than tuck it behind a disclosure.

The register is precise, sober and unembellished. Density is earned, never performed — tables run dense because a scientist is comparing rows, while prose stays at a readable measure. Accent colour is scarce and always means something. The world is deliberately **single-theme dark** and will not be given a light mode.

**Key Characteristics:**

- Near-black blue ground (`#070a12`) with a four-step tonal surface ramp instead of shadows
- One scarce blue accent reserved for live state, current selection and primary action
- Outfit for display and labels, Inter for body, Fira Code only for code, data and measurement
- Semantic colour vocabulary that is stable across every surface: blue running, emerald complete, amber waiting, rose refused, neutral idle
- Evidence stays visible — refusals, limits and negative results are designed for, not hidden

## Colors

A near-black blue ground carrying cool slate neutrals and 400-level semantic accents under a single blue signal. Chosen on 28 Oct 2026 from six candidate colourways, each verified against WCAG AA before the choice was made.

### Primary

- **Signal Blue** (`#4d9fff`): The only true brand colour. Reserved for live state (a running round, a pulsing status dot), current selection (active rail item, current tab), primary action, and in-scope residues. It is never used for decoration, never for a heading, and never as a gradient. Its token is `--accent`, deliberately hue-neutral: the accent has moved once already and the name should survive the next move.
- **Blue Wash** (`rgba(77,159,255,0.14)`): The tint that carries Signal Blue onto a surface — active nav items, running pills, in-scope residue chips. Always paired with full-strength blue text, never used as a standalone fill. Blue sits darker than the cyan it replaced, so the tint runs two points stronger to hold the same presence.

### Secondary

The semantic states. These are a vocabulary, not a palette: each one means exactly one thing everywhere it appears.

- **Verified Emerald** (`#34d399`): Complete, accepted, approved, improved. Also the colour of a best-fitness figure that beat its parent.
- **Waiting Amber** (`#fbbf24`): Blocked on something outside the system — a campaign waiting on bench data — and the prototype banner. Also the warning threshold on a budget meter.
- **Refused Rose** (`#f43f5e`): A ledger refusal, a hard-banned position, a rejected import row. Used at full strength for borders and marks.
- **Legible Rose** (`#fb7185`): Rose text on a rose tint. `#f43f5e` on `rgba(244,63,94,0.12)` fails AA for the small tracked caps that use it; this lighter step measures 5.57:1 over the card. Text gets this value; borders keep Refused Rose.

### Tertiary

- **Analytic Violet** (`#c084fc`): A third figure colour, used only to separate a third statistic in a row of metrics so two adjacent numbers never share a hue. It carries no state meaning and must not acquire one. It is deliberately the thinnest-used token in the palette — two figures — and that is the whole job: a third hue exists so a row of three numbers never forces an accent or a status colour into decorative duty. If a future row needs no third figure, retire the token rather than finding it work.

### Neutral

The surface ramp is the depth system. Each step is a distinct plane, not a shade.

- **Ground** (`#070a12`): The page itself, and the inset of a field or a log box. The deepest plane.
- **Panel** (`#0e1424`): The console rail and the honest-limits callout — a second neutral layer for chrome that frames content rather than being content.
- **Card** (`#151d33`): The default content surface. Campaigns, rounds, clusters, exports.
- **Raised Card** (`#1e2a48`): Hover and the default button fill — one step up from Card, the only "lift" in the system.
- **Hairline** (`rgba(255,255,255,0.08)`) and **Hairline Strong** (`rgba(255,255,255,0.16)`): Every boundary in the system. The strong value is for interactive edges (inputs, ghost buttons, drop zones).
- **Primary Text** (`#f7f9fc`), **Muted Text** (`#9aa8bf`), **Dim Text** (`#8795ae`): Three steps only. Dim Text is set to clear AA on every surface including Raised Card, its tightest at 4.69:1; do not darken it.

### Derived values

Eleven further tokens carry values derived from the five above rather than new hues:
`--accent-hi` (primary hover), `--accent-ink` (text on an accent fill), `--accent-edge`
and `--accent-ring` (every accent border and the live ring), `--emerald-edge`,
`--amber-edge` and `--rose-edge`, `--bg-raised-hi` (default button hover), `--bg-sunk`
(the log box, inset below the ground), `--thumb-hi`, and the four white washes
`--wash-faint` / `--wash-soft` / `--wash` / `--track`. These were literals until the
accent tints had drifted to seven different alphas; they are tokens now so the next
palette move is one edit.

One literal is deliberate and stays: the prototype bar is pure `#000`, outside the
system on purpose so the scaffolding never reads as product chrome. It does not ship.

### Named Rules

**The One Live Thing Rule.** Signal Blue marks what is happening *now*. If a screen shows two different blue elements claiming liveness, one of them is wrong. Current selection and primary action may also carry it, but a screen with blue on more than roughly a tenth of its area has lost the signal.

**The Stable Meaning Rule.** A colour means the same thing on every surface. Emerald is never "a nice green for a heading"; rose is never "an attention grabber". If a new state needs a colour and none of the five fits, the state is probably neutral.

**The Tinted Text Rule.** Any accent text sitting on its own tint must be contrast-checked against the *composited* background, not against the tint value. Two of the five accents needed a lighter step to clear AA; assume the next one will too.

## Typography

**Display Font:** Outfit (fallback `ui-sans-serif, system-ui, sans-serif`)
**Body Font:** Inter (fallback `ui-sans-serif, system-ui, sans-serif`)
**Label/Mono Font:** Fira Code (fallback `ui-monospace, SF Mono, Menlo, monospace`)

**Character:** Outfit's geometric warmth keeps large figures and headings from reading as cold engineering output, while Inter does the unglamorous work of dense product copy. The pairing is conventional on purpose — in a tool that asks a scientist to trust a number, typography that draws attention to itself is a liability.

### The scale

Nine steps, in `rem` against the reader's own base size. Every size on the page is one
of these tokens; there are no literal font sizes left.

| Token | rem | at 16px | Used for |
|---|---|---|---|
| `--t-xs` | 0.6875 | 11px | Labels, chips, pills, table headers |
| `--t-sm` | 0.75 | 12px | Hints, meta, field labels |
| `--t-md` | 0.8125 | 13px | Dense copy, the event log |
| `--t-base` | 0.875 | 14px | Default UI text, buttons, nav |
| `--t-lg` | 0.9375 | 15px | Body copy, card headings |
| `--t-xl` | 1.1875 | 19px | Subheads, wordmark, quota amount |
| `--t-2xl` | 1.3125 | 21px | Metric values |
| `--t-3xl` | 1.5 | 24px | Pane and auth headings |
| `--t-4xl` | 1.875 | 30px | Stat values |

The three fluid steps (display, headline, and the hero claim figures) interpolate
*between* scale steps, so their `clamp()` endpoints are deliberately not scale tokens —
`1rem`, `2.75rem` and `4.125rem` exist only as the ends of a range. Fluid type is
confined to the landing page; nothing in the console uses it.

### Hierarchy

- **Display** (Outfit 800, `clamp(2.375rem, 7vw, 4.125rem)`, 1.2, `-0.03em`): The landing hero only. One per page.
- **Headline** (Outfit 700, `clamp(1.5rem, 4vw, 2.125rem)`, 1.2, `-0.02em`): Landing section headings. Fluid sizing is correct here and only here.
- **Title** (Outfit 700, `--t-3xl`, 1.2, `-0.02em`): Console pane headings. Fixed, not fluid — a console heading that shrinks inside a narrow work area looks broken, not responsive.
- **Subtitle** (Outfit 600, `--t-lg` to `--t-xl`): Card and cluster headings, round names, export titles.
- **Body** (Inter 400, `--t-lg`, 1.6): Default copy. Prose holds a 65–75ch measure; the landing lede caps at 60ch, section copy at 62ch, cluster rationale at 72ch.
- **Label** (Outfit 700, `--t-xs`, `0.08em`, uppercase): Table headers, metric keys, status pills, stage chips, rail group headings.
- **Code** (Fira Code 400, `--t-md`, ligatures off): Mutation notation, accessions, job IDs, worker names, the event log.

### Named Rules

**The Eleven Pixel Floor.** No step below `--t-xs`, and tracked uppercase especially. Labels were lifted off 9.5px and 10px because tracked caps at that size stop being readable and start being texture.

**The No Literal Sizes Rule.** Every `font-size` is a scale token. A new size is a decision about the scale, made here first — the page once carried 21 distinct sizes for nine roles (11.5, 12.5, 13.5, 14.5, 16, 17, 18, 26px among them), which is how a scale rots.

**The Reader's Base Size Rule.** The scale is in `rem` and nothing sets `html { font-size }`, so a reader who raises their browser default gets a bigger page. Layout must survive it: every row that can run out of space wraps. Verified at 200% text on a 400px viewport with no horizontal overflow.

**The Mono Means Data Rule.** Fira Code marks code, identifiers, measurements and tabular figures — `T203F:E222G`, `P42212`, `7f3a9c`, `$0.000512`. It is never a costume for "technical". Prose in a monospace face is a bug.

**The Tabular Figures Rule.** Every number that sits in a column or changes in place carries `font-variant-numeric: tabular-nums`. Costs, fitness values, counts and timestamps must not jitter as they update.

## Layout

The two surfaces use different spatial models on purpose.

**Landing** is a single centred scroll column: `max-width: 1080px`, 16px side gutters that respect `env(safe-area-inset-*)`. Sections are divided by a 1px hairline and carry asymmetric vertical padding — `76px` above, `56px` below — so a heading always has more air above it than beneath it.

**Console** is a fixed 236px rail plus a scrolling work area, collapsing to a single column at **860px** where the rail becomes a horizontal nav grid. The work area is a `grid` with a 22px gap. Panes that are read-and-decide surfaces (new campaign, quota) cap at `880px`; the campaigns pane stays full width because it is a table.

Rhythm runs on six tokens — `--sp-xs` 6px, `--sp-sm` 10px, `--sp-md` 14px, `--sp-lg` 18px, `--sp-xl` 22px, `--sp-2xl` 26px. Related content groups tightly (`xs`–`sm`); distinct groups separate generously (`lg`–`2xl`). The markup carried ten different inline margin values before these existed, which is how a grid ends up fighting its own children.

**Breakpoints:** `620px` (hero claim row stacks, the prototype note is dropped), `860px` (console rail collapses), `900px` (run detail's two-column grid collapses).

### Named Rules

**The Shrinkable Item Rule.** Every grid or flex item that can contain a table carries `min-width: 0`. A grid item defaults to `min-width: auto`, which sizes it to the min-content of its widest `nowrap` table and drags the whole page into horizontal scroll. An `overflow-x: auto` wrapper cannot save a parent that refuses to shrink.

**The Sticky Clearance Rule.** The prototype bar is sticky, so `html` carries `scroll-padding-top: 64px`. Any in-page scroll target — an anchor, a pane switch calling `scrollIntoView` — must clear fixed chrome, or it lands underneath it.

**The No Horizontal Scroll Rule.** The page holds at 360px with no horizontal overflow. Wide tables scroll inside their own wrapper; the document never does.

## Elevation & Depth

**This system is flat.** Depth comes from the four-step tonal surface ramp (`#060910` → `#0d1321` → `#141d30` → `#1c2944`) plus 1px hairline borders — not from shadows. A card is "above" the page because it is lighter and outlined, not because it casts anything.

The one exception is a **state ring**, not an elevation: a live round carries `box-shadow: 0 0 0 1px rgba(77,159,255,0.10)` as a second hairline reinforcing its blue border. Zero offset, zero blur. It is a ring, and it is the only `box-shadow` in the console.

The presentation deck (`presentation.html`) does use a glow vocabulary — `0 0 24px rgba(0,240,255,0.35)` halos, `0 10px 30px rgba(0,0,0,0.5)` drops. **That vocabulary is scoped to presentation and marketing surfaces only.** The deck has not yet been brought onto this palette, so its halos are still cyan; when it is updated, the halo hue follows `--accent`. It is a legitimate authored effect where the job is to persuade, and it does not transfer into product UI, where a coloured halo would be decoration competing with the one live signal.

### Named Rules

**The Tonal-Not-Shadow Rule.** Product UI conveys depth by surface step and hairline. If a console element needs a shadow to read as separate, the real problem is that it is on the wrong surface step.

**The Deck-Only Glow Rule.** Accent halos belong to the deck. A zero-offset coloured glow in the console is decoration, and the console has no decoration budget.

## Shapes

Softly rounded rectangles throughout, with a deliberately narrow vocabulary: **small (`--r-sm`, 6px)** for controls — buttons, inputs, residue chips, nav items; **medium (`--r-md`, 9px)** for containers — cards, clusters, rounds, exports; **large (`--r-lg`, 10px)** for the outermost card and drop zone. Stage chips take `--r-chip` (4px), tighter than a control because they are read as a strip rather than pressed. Status pills are fully rounded (`--r-pill`); step numbers and the live-state dot are circles (`50%`); the budget meter uses `--r-xs` (3px).

Borders do the work that shadows would elsewhere: a 1px hairline on nearly every surface, stepping to `rgba(255,255,255,0.16)` on interactive edges and to a 30–40% accent tint when an element carries state. A **dashed** border marks something provisional or additive — an excluded cluster, an "add position" control, the measurements drop zone.

### Named Rules

**The No Accent Rail Rule.** No coloured `border-left` or `border-right` above 1px on cards, list items or callouts. Accent rails were removed from stats and the honest-limits callout deliberately; state belongs in the figure's own colour, not in a stripe beside it.

**The Dashed Means Provisional Rule.** A dashed border is reserved for the not-yet and the set-aside. A solid border is a committed thing.

## Components

Buttons, cards and inputs are **refined and restrained** — quiet at rest, with state carried by colour and border rather than motion or depth. The control recedes until it matters. This holds hardest on the dense technical surfaces, where a row of eager controls would compete with the data.

### Buttons

- **Shape:** Small radius (`6px`); `8px` on the larger default size.
- **Primary:** Signal Blue fill with near-black ink (`#05101f`), `10px 18px`. Exactly one per view.
- **Default:** Raised Card fill (`#1c2944`), primary text. Hover lifts to `#253457`.
- **Ghost:** Transparent with a strong hairline and muted text; hover shifts border and text to blue over a blue wash.
- **Hover / Active:** `background` and `border-color` transition at `150ms`; `:active` nudges `translateY(1px)`. Both are suppressed under `prefers-reduced-motion`.
- **Disabled:** `opacity: 0.4`, `cursor: not-allowed`, pointer events off.
- **Small variant:** `12px` type, `6px 12px` padding.

### Chips and Pills

- **Pills** (status): Fully rounded, Outfit 700 at `11px` / `0.08em` uppercase, accent text on its own 12% tint with a 30% accent border. Five states: running (blue, with a pulsing dot), complete (emerald), waiting (amber), refused (rose), idle (neutral).
- **Stage chips** (the six loop steps): `4px` radius, same label type, tinted by state — done (emerald), now (blue, bordered), failed (rose), pending (neutral).
- **Residue chips:** Fira Code at `12.5px`, `6px` radius. In-scope is blue on blue wash; excluded is neutral with a line-through; banned is Legible Rose on rose wash.

### Cards / Containers

- **Corner:** Medium (`9px`); `10px` on the base `.card`.
- **Background:** Card (`#141d30`), or Panel (`#0d1321`) for framing chrome and set-aside content.
- **Border:** 1px hairline; accent-tinted at 30–40% when the card carries state.
- **Shadow:** None. See Elevation & Depth.
- **Padding:** `--sp-lg` default, `--sp-md` for denser containers. A card wrapping a full-bleed table takes the `.flush` variant (`padding: 0`) and lets the cells supply their own inset.

### Inputs / Fields

- **Style:** Ground fill (`#060910`) inset into the card — inputs are *recessed*, not raised — with a strong hairline and `6px` radius.
- **Focus:** Border shifts to Signal Blue. Caret is blue; `accent-color` is blue.
- **Label:** `12px` Inter 600 in muted text, above the field. Hints sit below in dim text at `12px`.

### Navigation

- **Rail:** Panel background, 236px. Items are left-aligned buttons at `14px` with `6px` radius; hover is a 4% white wash, current is blue on blue wash at weight 600. Group headings are `11px` tracked caps in dim text. Sign out is separated by a hairline and set in dim text — it is not a destination.
- **Tabs:** Underline style with a 2px bottom border; current is blue text and a blue rule. Implemented as a real `role="tablist"` with `aria-selected`, roving `tabindex` and arrow-key navigation. The campaign detail runs eight of them — Rounds, Planner, Search space, Variants, Ledger, Wet-lab handoff, Measurements, Economics — which is past the five a tab bar comfortably holds. They are kept because each is a distinct output of the engine rather than a slice of one; the bar wraps rather than scrolling, so every tab stays visible at any width.
- **Mobile:** Below 860px the rail becomes a horizontal flex row. Group headings and the divider are hidden — they organise a vertical list and only add noise to a three-item row — and Sign out is pushed to the far end with `margin-left: auto` so it never reads as a fourth destination. Buttons keep their natural width; an auto-fit grid would stretch them to equal columns and make Sign out look like one of the places you can go.

### Signature Component: The Budget Meter

The product's central promise rendered as a 6px bar. It **fills as money is spent** — the direction every budget meter moves — over a 7% white track, and it is placed directly beneath the sentence containing the word "spent" so proximity binds the bar to its meaning. Colour is the warning and nothing else: blue by default, amber past the warning threshold, rose near exhaustion. A non-zero fill carries `min-width: 3px` so a 1% spend reads as "barely any" rather than "none". It always carries `role="img"` and an `aria-label` stating the percentage used.

### Charts

Two inline SVGs: a latent-space scatter and a two-bar comparison. Both cap at **480px wide**,
because a `viewBox` scales its own text — left uncapped in a wide column, a 10px label renders at
22px and the chart reads as a poster. Chart text takes `--fg-dim` from the palette rather than a
browser default, series colour follows the round it belongs to (neutral, blue, emerald), and every
axis label names a value the chart actually reaches.

### The Residue Track

238 cells, one per position, built in script rather than markup. Catalytic positions are marked by
**height** rather than hue — they are a structural category, not a state, and the five status
colours are spoken for. In-scope positions take the accent, banned positions take rose, matching the
residue chips in the scope step so the two readings of the same fact agree.

### Signature Component: The Authorisation Projection

The budget step's second card, and the clearest statement of what the product is for. The user sets
six numbers — authorisation, per-round ceiling, model, batch size, round cap, no-improvement rule —
and the card answers in the same breath: the per-variant rate, the cost of one round, the cost if
every round runs, and the share of remaining quota. The last line names **which limit actually
bites**, in amber, because that is the number the user is really setting. Set a batch the ceiling
cannot cover and it says the ledger would refuse the first round; pick the mock scorer and it says
nothing binds. A budget field alone would make spend a number to type; this makes it a consequence
to read.

### Signature Component: The Stage Strip

Six chips naming the fixed loop steps — Ingest, Propose, Select, Authorise, Score, Record. It appears **only where the stage is the news**: on a running round and on a refused one (where the failing stage is the whole story). A completed round hides it, because six green chips only repeat what the "Complete" pill already said.

## Do's and Don'ts

### Do:

- **Do** keep Signal Blue scarce. One live thing per screen; primary action and current selection may share it, but roughly a tenth of the screen is the ceiling.
- **Do** convey depth with the surface ramp and 1px hairlines, never with a shadow, in product UI.
- **Do** check accent text against the **composited** background. `#f43f5e` on `rgba(244,63,94,0.12)` over the card is 4.11:1, not the ratio the raw values suggest — which is why text uses Legible Rose.
- **Do** give every grid or flex item that may contain a table `min-width: 0`.
- **Do** use `tabular-nums` on every figure in a column or one that updates in place.
- **Do** set `scroll-padding-top` to clear the sticky bar on any in-page scroll target.
- **Do** show refusals, limits and negative results at full contrast. The ledger's refused rows and the benchmark's held-out correlation are designed for, not disclosed away.
- **Do** keep tracked uppercase labels at 11px or above.
- **Do** theme the surfaces the browser draws — selection, caret, `accent-color`, scrollbars. Default light scrollbars on a near-black page are the giveaway that nobody looked.

### Don't:

- **Don't** add a light theme or a theme toggle. The dark-only commitment is intentional and binding.
- **Don't** reach for gradient text, neon halos, glass or blur as decoration, or animated mesh backgrounds. Confirmed anti-reference.
- **Don't** soften the product into consumer-SaaS friendliness — pastel rounding, illustrations, mascots, emoji as iconography, chirpy microcopy. Confirmed anti-reference.
- **Don't** wear monospace as a costume. Fira Code is for code, identifiers, data and measurement; no CRT scanlines, no ASCII chrome. Confirmed anti-reference.
- **Don't** put a kicker or eyebrow above a heading. A standalone label that *is* the heading of a card is fine; a tracked-caps line stacked on top of an `h2` is not, and these were removed deliberately.
- **Don't** give a card a coloured accent rail above 1px.
- **Don't** de-emphasise with `opacity` on anything containing text. An excluded cluster was at `0.45` and became unreadable — move it to the Panel surface and keep its text legible instead.
- **Don't** use fluid `clamp()` type in the console. Fluid sizing belongs to the landing page; a console heading that shrinks in a narrow work area looks broken.
- **Don't** let a status colour drift from its meaning, and don't give Analytic Violet a state meaning — it exists only to separate a third figure.
- **Don't** add motion that doesn't carry state. The pulsing dot on a running pill is the one authored moment, and it honours `prefers-reduced-motion`.

---

*Known divergence (1):* the presentation deck (`presentation.html`) now carries the PanoptesBio name but is still on the cyan accent and its glow vocabulary. The console is normative for colour; the deck's palette is pending.

*Known divergence (2):* the Gradio dashboard (`ui/dashboard.py`) implements a different palette — slate-900/800 surfaces (`#0f172a`, `#1e293b`), sky-400 (`#38bdf8`) in place of Signal Blue, violet-600, and JetBrains Mono in place of Fira Code. It shares only the Outfit/Inter pairing and the semantic accents. It predates this system and is **not** a second sanctioned theme; reconcile it to the tokens above rather than treating it as precedent.
