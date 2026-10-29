---
name: dev-ui-viewer
description: >-
  Use this skill when building or updating user interfaces, Gradio web dashboards,
  3Dmol WebGL protein structure visualizers, or interactive fitness landscape plots.
---

# UI & 3D Visualization Development Skill

This skill governs the presentation layer (`ui/dashboard.py`, `ui/structure_viewer.py`, `ui/landscape_plot.py`).

## Core Invariants
1. **Never Block on Cloud Compute:**
   - Gradio UI handlers must NEVER synchronously wait for Modal GPU executions.
   - Run computations asynchronously or in background tasks; poll the SQLite ledger using generator `yield` patterns to stream status updates.
2. **Aesthetic Excellence:**
   - Curated dark-theme palettes (e.g. Slate/Emerald/Indigo), modern typography (Inter/Outfit), smooth CSS transitions, and glassmorphic cards.
3. **Authentic 3Dmol WebGL Grounding:**
   - Always load authentic crystallographic PDB coordinates (e.g. `1EMA.pdb`, `2YBR.pdb`, `8EH0.pdb`) into the 3Dmol viewer.
   - Highlight active sites, CDR loops, or binding interfaces with distinct cartoon/stick representations (e.g., amber gold for paratope residues).
4. **Human Approval Actions:**
   - Provide dedicated interactive buttons for human sign-offs (e.g., "Approve Lead & Generate Twist Synthesis Order", "Export Opentrons Robot Protocol").

## Key Files
- [ui/dashboard.py](file:///c:/Users/stoof/source/repos/PanoptesBio/ui/dashboard.py): Main Gradio blocks application with responsive tabs.
- [ui/structure_viewer.py](file:///c:/Users/stoof/source/repos/PanoptesBio/ui/structure_viewer.py): 3Dmol WebGL HTML renderer with custom selection styles.
- [ui/landscape_plot.py](file:///c:/Users/stoof/source/repos/PanoptesBio/ui/landscape_plot.py): Plotly interactive 2D/3D fitness landscape and Pareto frontier plots.

## Verification Checklist
```powershell
python -m unittest tests/test_dashboard.py
```
Ensure:
- Gradio application instantiates cleanly without schema errors.
- 3Dmol viewer outputs valid HTML/JavaScript strings.
