---
name: PanoptesBio Design System
description: Institutional-grade biosecurity, deterministic FinOps cloud economics, and high-fidelity molecular biophysics.
colors:
  primary: "#0284c7"
  primary-bright: "#38bdf8"
  emerald-primary: "#059669"
  emerald-bright: "#34d399"
  amber-primary: "#d97706"
  amber-bright: "#fbbf24"
  neutral-bg: "#070a13"
  neutral-card: "#0d1322"
  neutral-elevated: "#0f172a"
  text-main: "#f8fafc"
  text-muted: "#cbd5e1"
  text-dim: "#94a3b8"
  border-subtle: "rgba(56, 189, 248, 0.16)"
  border-active: "rgba(56, 189, 248, 0.45)"
  border-emerald: "rgba(52, 211, 153, 0.38)"
typography:
  display:
    fontFamily: "'Outfit', sans-serif"
    fontSize: "2.75rem"
    fontWeight: 700
    lineHeight: 1.15
    letterSpacing: "-0.03em"
  headline:
    fontFamily: "'Outfit', sans-serif"
    fontSize: "1.5rem"
    fontWeight: 700
    lineHeight: 1.25
    letterSpacing: "-0.02em"
  body:
    fontFamily: "'Outfit', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.55
    letterSpacing: "0em"
  mono:
    fontFamily: "'JetBrains Mono', monospace"
    fontSize: "0.8125rem"
    fontWeight: 500
    lineHeight: 1.45
    letterSpacing: "0em"
rounded:
  sm: "6px"
  md: "10px"
  lg: "16px"
  full: "9999px"
spacing:
  xs: "4px"
  sm: "8px"
  md: "16px"
  lg: "24px"
  xl: "36px"
components:
  button-primary:
    backgroundColor: "{colors.primary}"
    textColor: "#ffffff"
    rounded: "{rounded.sm}"
    padding: "10px 20px"
  button-primary-hover:
    backgroundColor: "{colors.primary-bright}"
  button-secondary:
    backgroundColor: "rgba(15, 23, 42, 0.8)"
    textColor: "{colors.text-main}"
    rounded: "{rounded.sm}"
    padding: "10px 18px"
---

# PanoptesBio Design System

## Overview

PanoptesBio is an institutional-grade cloud-economics engine and auditable active-learning platform for computational biologics. The visual design conveys high-stakes mathematical rigor, cryptographic immutability (FDA 21 CFR Part 11), and biophysical truth. It refuses the generic, ungrounded dark-mode AI tropes (fuzzy neon halos, floating unanchored spheres) in favor of deep obsidian density, high-contrast typography, razor-sharp tabular numbers, and precision laboratory metadata.

## Colors

- **Ground (Obsidian Depth):** `#070a13` serves as the primary canvas, layered with `#0d1322` for glassmorphic cards and `#0f172a` for elevated control panels.
- **Biophysical Cyan (Active State & Navigation):** Primary `#0284c7`, Bright `#38bdf8`. Used for interactive states, primary action triggers, and molecular focus highlights.
- **Cryptographic Emerald (Audit & Validation):** Primary `#059669`, Bright `#34d399`. Used for verified FDA Part 11 Merkle hash blocks, live server connectivity dots, and positive fitness discoveries.
- **FinOps Amber (Spend & Constraints):** Primary `#d97706`, Bright `#fbbf24`. Used for budget burn alerts, Pareto frontier rank-0 non-dominated points, and developability risk tags.
- **Alert Red (Fail-Closed Gates):** `#ef4444`. Used for rejected budget attempts, steric clash warnings, and sequence liability alerts.

### 12 Laboratory & Biophysical Color Themes
The header includes an interactive theme picker allowing users to toggle between 12 calibrated laboratory palettes:
1. **Panoptes Obsidian (Default):** Deep void (`#070a13`) with cyan (`#0284c7`) and sentinel emerald (`#10b981`).
2. **Cryo-EM Glacial:** Vitrified ice (`#040914`) with electron density cyan (`#06b6d4`) and azure (`#38bdf8`).
3. **avGFP Bioluminescence:** Aequorea victoria green (`#040f09`) with 509nm fluorophore green (`#10b981`) and lime (`#84cc16`).
4. **Active-Site Catalysis:** PLP cofactor core (`#0f0a05`) with catalytic amber (`#d97706`) and gold (`#fbbf24`).
5. **AlphaFold Topology:** Deep invariant point attention (`#090616`) with structural indigo (`#6366f1`) and violet (`#a855f7`).
6. **Rhodopsin Retinal:** Opsin proton pump (`#100508`) with photochemical crimson (`#e11d48`) and rose (`#fb7185`).
7. **Microfluidic Flow-Cell:** Droplet sorting dark (`#040d10`) with laminar stream teal (`#0d9488`) and turquoise (`#2dd4bf`).
8. **Cy5 Far-Red Laser:** Single-molecule FRET (`#0d040e`) with optical laser ruby (`#f43f5e`) and hot magenta (`#ec4899`).
9. **HPLC Chromatography:** Size-exclusion UV trace (`#0d0a04`) with analytical bronze (`#ca8a04`) and yellow gold (`#eab308`).
10. **Cleanroom Titanium:** ISO 5 pharmaceutical cleanroom (`#090b10`) with passivated steel (`#64748b`) and slate (`#94a3b8`).
11. **Superoxide Cadmium:** Reactive radical photobleaching resistance (`#100704`) with cadmium orange (`#ea580c`).
12. **Cobalt Oligo Synthesis:** Magnetic bead purification (`#040714`) with phosphoramidite cobalt (`#2563eb`) and blue (`#60a5fa`).

## Typography

- **Display & Headlines:** `Outfit` (700 weight, `-0.02em` to `-0.03em` tracking). Delivers sharp institutional authority without looking like generic tech monoculture.
- **Body & UI Text:** `Outfit` / System Sans (`#f8fafc` text main, `#cbd5e1` muted secondary text, `#94a3b8` metadata). Maintains WCAG AA compliant contrast ratio ($\ge 4.5:1$).
- **Data & Cryptography:** `JetBrains Mono` (`font-variant-numeric: tabular-nums`). Exclusively used for micro-USD budget figures, PDB coordinates, mutation strings (`S65T:Y66H`), and HMAC-SHA256 hashes.

## Layout

- Max container width constrained to `1440px` with responsive padding (`24px 20px 48px`).
- Generous separation between macro-sections (`48px` gap) and tight grouping within functional modules (`8px` to `16px`).
- Responsive multi-column grid adapting from 1 column on mobile to 3–4 columns on desktop displays.

## Elevation & Depth

- Cards utilize subtle, multi-layered dark elevation: `box-shadow: 0 8px 30px rgba(0, 0, 0, 0.45), 0 0 0 1px rgba(255, 255, 255, 0.05)`.
- Subtle glassmorphism (`backdrop-filter: blur(16px)`) with `1px` translucent borders (`rgba(56, 189, 248, 0.16)`).
- Zero-blur block shadows and bright neon halos are strictly banned.

## Shapes

- Radii scale:
  - Small controls and inputs: `6px`
  - Cards and modules: `10px` to `16px`
  - Status badges and identity pills: `9999px` (full pill)

## Components

- **Header:** Features the PanoptesBio identity with live Firebase domain pill, active operator identity badge, mode switcher (`Overview & Proof` vs `Scientific Studio`), and sync controls.
- **KPI Bar:** 5-metric quantitative ribbon tracking remaining budget, cumulative spend, top discovered fitness score, capital savings percentage, and efficiency multiple.
- **Landing Showcase:** Executive persuade surface featuring real-world 4-target benchmark proof carousel, interactive cloud economics calculator, and the 6-step loop architecture.
- **Scientific Studio:** 5 focused workbenches (Active Learning Campaign Studio, De Novo Diffusion, Pareto Frontier Studio, Physical Biology & LIMS, and FDA Part 11 Audit Trail).

## Do's and Don'ts

- **DO** use `font-variant-numeric: tabular-nums` for all financial figures, scores, and residue positions.
- **DO** verify contrast ratios $\ge 4.5:1$ against dark surfaces.
- **DO** present real measured numbers from the production DMS datasets.
- **DON'T** use gradient text or colored halo glows.
- **DON'T** put kickers or eyebrows above main section titles.
- **DON'T** introduce generic SVG doodles or ungrounded decorative illustrations.

