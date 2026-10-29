"""
ui/structure_viewer.py — Interactive 3D Molecular Structure Viewer (3Dmol.js).

Purpose:
    Generates an enterprise-grade WebGL 3D molecular visualization component
    for the Gradio dashboard using 3Dmol.js.
    Renders the authentic crystal structure of avGFP (PDB ID: 1EMA),
    with live highlights for the catalytic chromophore and selected active-learning mutations.
    Uses an iframe srcdoc wrapper so browser executes 3Dmol scripts cleanly without Gradio warnings.
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import List, Optional, Tuple

PDB_FILE_PATH = Path(__file__).resolve().parent.parent / "workers" / "data" / "1EMA.pdb"


def resolve_target_pdb_context(
    target_accession: str = "P42212",
    pdb_path: Optional[Path] = None,
) -> Tuple[Path, str, str, List[int], str, str]:
    """
    Dynamically resolves the 3D crystal structure, structural domain annotations,
    catalytic sites, and color scheme for any protein target.

    Returns:
        Tuple of (file_path, title, scaffold_label, catalytic_residues, catalytic_label, scaffold_color).
    """
    acc_lower = str(target_accession).lower()
    data_dir = Path(__file__).resolve().parent.parent / "workers" / "data"

    if pdb_path and pdb_path.is_file():
        stem = pdb_path.stem.upper()
        return (
            pdb_path,
            f"3D Structure (PDB: {stem})",
            "Engineered Scaffold",
            [],
            "Active Site",
            "#0284c7",
        )

    if any(k in acc_lower for k in ("gb1", "spg1", "2qmt")):
        gb1_pdb = data_dir / "2QMT.pdb"
        return (
            gb1_pdb,
            "3D Crystal Structure (PDB: 2QMT · 1.10 Å Resolution)",
            "Immunoglobulin G-Binding Protein G (GB1)",
            [39, 40, 41, 54],
            "IgG Fc Binding Interface (39, 40, 41, 54)",
            "#8b5cf6",  # Royal violet
        )
    elif any(k in acc_lower for k in ("trastuzumab", "her2", "herceptin", "1n8z")):
        her2_pdb = data_dir / "1N8Z.pdb"
        return (
            her2_pdb,
            "3D Crystal Structure (PDB: 1N8Z · 2.50 Å Complex)",
            "Herceptin Fab (VH Domain) : HER2 Complex",
            [99, 100, 101, 102, 103, 104, 105, 106, 107, 108],
            "CDR-H3 Hypervariable Loop (Residues 99-108)",
            "#ec4899",  # Rose pink
        )
    elif any(k in acc_lower for k in ("scorpion", "cn2", "9004g", "2ybr")):
        scorp_pdb = data_dir / "2YBR.pdb"
        return (
            scorp_pdb,
            "3D Crystal Structure (PDB: 2YBR · 2.55 Å Complex)",
            "Neutralizing scFv 9004G : Centruroides noxius Cn2 Scorpion Toxin Complex",
            [99, 100, 101, 102, 103, 104, 105, 106],
            "CDR-H3 Paratope Loop (Residues 99-106 GGVGSFDT · Blocks NaV Channel Interface)",
            "#d97706",  # Amber gold
        )
    elif any(k in acc_lower for k in ("trpb", "tm9d8", "tryptophan", "p50909")):
        trpb_pdb = data_dir / "8EH0.pdb"
        return (
            trpb_pdb,
            "3D Crystal Structure (PDB: 8EH0 · 1.70 Å Resolution)",
            "TIM-Barrel (α/β)8 Fold",
            [83],
            "PLP Catalytic Lys83",
            "#2563eb",  # Sapphire blue
        )
    else:
        gfp_pdb = data_dir / "1EMA.pdb"
        return (
            gfp_pdb,
            "3D Crystal Structure (PDB: 1EMA · 1.90 Å Resolution)",
            "11-Stranded β-Barrel",
            [65, 66, 67],
            "SYG Chromophore Core",
            "#059669",  # Emerald green
        )


def load_pdb_content(pdb_path: Optional[Path] = None) -> str:
    """Loads and sanitizes PDB file content (defaults to 1EMA)."""
    target_path = pdb_path if pdb_path is not None else PDB_FILE_PATH
    if not target_path.is_file():
        return ""
    try:
        with open(target_path, "r", encoding="utf-8") as f:
            lines = [line for line in f if line.startswith(("ATOM", "HETATM", "TER", "END"))]
            return "".join(lines)
    except Exception:
        return ""


def generate_3d_viewer_html(
    active_mutations: Optional[List[int]] = None,
    viewer_id: str = "veritasbio-mol-viewer",
    spin: bool = True,
    target_accession: str = "P42212",
    pdb_path: Optional[Path] = None,
) -> str:
    """
    Builds the interactive HTML/JS block for 3Dmol.js wrapped in an isolated iframe srcdoc.
    Supports dynamic PDB structures (1EMA for avGFP, 8EH0 for TrpB) and catalytic labeling.
    """
    (
        resolved_pdb_path,
        title_str,
        scaffold_label,
        catalytic_residues,
        catalytic_label,
        scaffold_color,
    ) = resolve_target_pdb_context(target_accession=target_accession, pdb_path=pdb_path)

    pdb_data = load_pdb_content(resolved_pdb_path)
    escaped_pdb = json.dumps(pdb_data)

    if active_mutations is not None:
        effective_mutations = active_mutations
    else:
        effective_mutations = [183, 184, 227, 228] if "8EH0" in str(resolved_pdb_path) else [65, 66, 163, 203]

    mutations_json = json.dumps(effective_mutations)
    catalytic_json = json.dumps(catalytic_residues)
    cat_label_escaped = json.dumps(catalytic_label)

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap">
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0;
            padding: 0;
            background: #0b1329;
            font-family: 'Inter', sans-serif;
            overflow: hidden;
            user-select: none;
        }}
        .top-toolbar {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 16px;
            background: #0f172a;
            border-bottom: 1px solid #1e293b;
            color: #f8fafc;
            font-size: 13px;
        }}
        .btn {{
            background: #1e293b;
            color: #cbd5e1;
            border: 1px solid #334155;
            padding: 4px 10px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 11px;
            font-weight: 500;
            transition: all 0.15s ease;
        }}
        .btn:hover {{
            background: #334155;
            color: #ffffff;
        }}
        .footer-legend {{
            display: flex;
            justify-content: space-around;
            align-items: center;
            padding: 8px 16px;
            background: #0b1329;
            border-top: 1px solid #1e293b;
            font-family: 'JetBrains Mono', monospace;
            font-size: 11px;
            color: #94a3b8;
        }}
    </style>
</head>
<body>
    <div style="width: 100%; height: 530px; display: flex; flex-direction: column;">
        <!-- Top Toolbar -->
        <div class="top-toolbar">
            <div style="display: flex; align-items: center; gap: 8px;">
                <span style="display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: #10b981; box-shadow: 0 0 8px #10b981;"></span>
                <span style="font-weight: 600; letter-spacing: 0.3px;">{title_str}</span>
                <span style="background: #1e293b; color: #94a3b8; padding: 2px 8px; border-radius: 9999px; font-size: 11px;">{scaffold_label}</span>
            </div>
            <div style="display: flex; gap: 8px;">
                <button onclick="toggleMolSpin()" id="btn-spin" class="btn" style="color: #38bdf8;">Pause Rotation</button>
                <button onclick="toggleMolSurface()" id="btn-surface" class="btn" style="color: #a78bfa;">Toggle Surface</button>
                <button onclick="resetMolCamera()" class="btn">Reset View</button>
            </div>
        </div>

        <!-- 3Dmol WebGL Container -->
        <div id="{viewer_id}" style="flex: 1; width: 100%; position: relative;"></div>

        <!-- Legend Footer -->
        <div class="footer-legend">
            <div style="display: flex; align-items: center; gap: 6px;">
                <span style="display: inline-block; width: 12px; height: 12px; border-radius: 2px; background: {scaffold_color};"></span>
                <span>{scaffold_label}</span>
            </div>
            <div style="display: flex; align-items: center; gap: 6px;">
                <span style="display: inline-block; width: 12px; height: 12px; border-radius: 2px; background: #f59e0b;"></span>
                <span>{catalytic_label}</span>
            </div>
            <div style="display: flex; align-items: center; gap: 6px;">
                <span style="display: inline-block; width: 12px; height: 12px; border-radius: 2px; background: #f43f5e;"></span>
                <span>Cofactor / Ligand (HETATM)</span>
            </div>
            <div style="display: flex; align-items: center; gap: 6px;">
                <span style="display: inline-block; width: 12px; height: 12px; border-radius: 50%; background: #06b6d4; box-shadow: 0 0 6px #06b6d4;"></span>
                <span>Active Candidate Mutations</span>
            </div>
        </div>
    </div>

    <!-- 3Dmol.js CDN Injection -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
    <script>
    (function() {{
        let viewerContainer = document.getElementById("{viewer_id}");
        if (!viewerContainer) return;

        let glviewer = null;
        let isSpinning = {"true" if spin else "false"};
        let hasSurface = false;
        let surfObj = null;

        function initViewer() {{
            if (typeof $3Dmol === 'undefined') {{
                setTimeout(initViewer, 100);
                return;
            }}

            viewerContainer.innerHTML = "";
            glviewer = $3Dmol.createViewer(viewerContainer, {{
                backgroundColor: "#0b1329",
                id: "{viewer_id}"
            }});

            let pdbString = {escaped_pdb};
            if (!pdbString || pdbString.length < 50) {{
                viewerContainer.innerHTML = "<div style='color:#ef4444; padding:20px; font-family:sans-serif;'>PDB file not found or empty.</div>";
                return;
            }}

            glviewer.addModel(pdbString, "pdb");

            // 1. Base scaffold cartoon ribbon style
            glviewer.setStyle({{}}, {{
                cartoon: {{
                    color: "{scaffold_color}",
                    opacity: 0.85,
                    thickness: 0.4
                }}
            }});

            // 2. Catalytic invariant residues in golden sticks
            let catalyticResidues = {catalytic_json};
            let catLabel = {cat_label_escaped};
            if (catalyticResidues && catalyticResidues.length > 0) {{
                glviewer.setStyle({{resi: catalyticResidues}}, {{
                    stick: {{
                        color: "#f59e0b",
                        radius: 0.35
                    }},
                    cartoon: {{
                        color: "#f59e0b",
                        opacity: 1.0
                    }}
                }});
                let centerResi = catalyticResidues[Math.floor(catalyticResidues.length / 2)];
                glviewer.addLabel(catLabel, {{
                    position: {{resi: centerResi}},
                    backgroundColor: "rgba(15, 23, 42, 0.85)",
                    fontColor: "#fbbf24",
                    fontSize: 11,
                    borderThickness: 1,
                    borderColor: "#f59e0b"
                }});
            }}

            // 2b. Catalytic cofactors & bound ligands (HETATM: PLP, 0JO, HEM, etc.) in rose crimson sticks & spheres
            glviewer.addStyle({{hetflag: true}}, {{
                stick: {{
                    color: "#f43f5e",
                    radius: 0.32
                }},
                sphere: {{
                    color: "#fb7185",
                    radius: 0.45,
                    opacity: 0.95
                }}
            }});

            // 3. Highlight active candidate mutations in glowing cyan spheres
            let activeMutPositions = {mutations_json};
            for (let i = 0; i < activeMutPositions.length; i++) {{
                let pos = activeMutPositions[i];
                if (catalyticResidues && catalyticResidues.indexOf(pos) !== -1) continue;
                
                glviewer.addStyle({{resi: pos}}, {{
                    sphere: {{
                        color: "#06b6d4",
                        radius: 0.85,
                        opacity: 0.95
                    }},
                    stick: {{
                        color: "#38bdf8",
                        radius: 0.25
                    }}
                }});
                
                glviewer.addLabel("Pos " + pos, {{
                    position: {{resi: pos}},
                    backgroundColor: "rgba(15, 23, 42, 0.85)",
                    fontColor: "#38bdf8",
                    fontSize: 10,
                    borderThickness: 1,
                    borderColor: "#0284c7"
                }});
            }}

            glviewer.zoomTo();
            glviewer.render();

            if (isSpinning) {{
                glviewer.spin("y", 0.6);
            }}
        }}

        window.toggleMolSpin = function() {{
            if (!glviewer) return;
            isSpinning = !isSpinning;
            let btn = document.getElementById("btn-spin");
            if (isSpinning) {{
                glviewer.spin("y", 0.6);
                if (btn) btn.innerText = "Pause Rotation";
            }} else {{
                glviewer.spin(false);
                if (btn) btn.innerText = "Resume Rotation";
            }}
        }};

        window.toggleMolSurface = function() {{
            if (!glviewer) return;
            let btn = document.getElementById("btn-surface");
            if (!hasSurface) {{
                surfObj = glviewer.addSurface($3Dmol.SurfaceType.VDW, {{
                    opacity: 0.45,
                    color: "#0284c7"
                }}, {{}});
                hasSurface = true;
                if (btn) btn.innerText = "Hide Surface";
            }} else {{
                if (surfObj) glviewer.removeSurface(surfObj);
                hasSurface = false;
                if (btn) btn.innerText = "Toggle Surface";
            }}
            glviewer.render();
        }};

        window.resetMolCamera = function() {{
            if (!glviewer) return;
            glviewer.zoomTo();
            glviewer.render();
        }};

        initViewer();
    }})();
    </script>
</body>
</html>"""

    escaped_srcdoc = html.escape(inner_html, quote=True)
    return f'<iframe srcdoc="{escaped_srcdoc}" style="width: 100%; height: 540px; border: none; border-radius: 12px; background: #0b1329; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.3);"></iframe>'


def generate_dual_comparison_viewer_html(
    target_accession: str = "TRPB_TM9D8",
    mutant_positions: Optional[List[int]] = None,
    lead_name: str = "Epistatic Lead Candidate",
    viewer_id: str = "veritasbio-dual-viewer",
    spin: bool = True,
    pdb_path: Optional[Path] = None,
) -> str:
    """
    Builds a synchronized dual-viewport 3D molecular viewer (WT Parent vs Optimized Lead)
    with explicit HETATM cofactor catalytic coordination shell rendering.
    """
    (
        resolved_pdb_path,
        title_str,
        scaffold_label,
        catalytic_residues,
        catalytic_label,
        scaffold_color,
    ) = resolve_target_pdb_context(target_accession=target_accession, pdb_path=pdb_path)

    pdb_data = load_pdb_content(resolved_pdb_path)
    escaped_pdb = json.dumps(pdb_data)

    effective_mutations = mutant_positions or ([183, 184, 227, 228] if "8EH0" in str(resolved_pdb_path) else [65, 66, 163, 203])
    mutations_json = json.dumps(effective_mutations)
    catalytic_json = json.dumps(catalytic_residues)
    spin_js = "true" if spin else "false"

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0; padding: 0; background: #0b1329; font-family: 'Inter', sans-serif;
            overflow: hidden; user-select: none; color: #f8fafc;
        }}
        .header {{
            display: flex; justify-content: space-between; align-items: center;
            padding: 8px 16px; background: #0f172a; border-bottom: 1px solid #1e293b;
        }}
        .title {{ font-size: 13px; font-weight: 600; color: #38bdf8; }}
        .badge {{
            background: rgba(56, 189, 248, 0.15); border: 1px solid rgba(56, 189, 248, 0.3);
            color: #38bdf8; padding: 2px 8px; border-radius: 4px; font-size: 11px;
            font-family: 'JetBrains Mono', monospace;
        }}
        .dual-container {{
            display: grid; grid-template-columns: 1fr 1fr; width: 100%; height: 450px;
            position: relative; background: #070a13;
        }}
        .viewport-box {{
            width: 100%; height: 100%; position: relative;
            border-right: 1px solid rgba(255,255,255,0.06);
        }}
        .viewport-tag {{
            position: absolute; top: 10px; left: 12px; z-index: 10;
            background: rgba(15, 23, 42, 0.85); border: 1px solid rgba(51, 65, 85, 0.6);
            padding: 4px 10px; border-radius: 6px; font-size: 11px; font-weight: 600;
        }}
        .toolbar {{
            position: absolute; bottom: 10px; left: 50%; transform: translateX(-50%);
            z-index: 20; display: flex; gap: 8px;
        }}
        .btn {{
            background: #1e293b; color: #cbd5e1; border: 1px solid #334155;
            padding: 5px 12px; border-radius: 6px; cursor: pointer; font-size: 11px;
            font-weight: 500; transition: all 0.15s ease;
        }}
        .btn:hover {{ background: #334155; color: #ffffff; }}
        .legend {{
            display: flex; justify-content: space-around; align-items: center;
            padding: 8px 16px; background: #0b1329; border-top: 1px solid #1e293b;
            font-family: 'JetBrains Mono', monospace; font-size: 11px; color: #94a3b8;
        }}
        .dot {{ display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; }}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <span class="title">Dual-Structure Comparison · WT Scaffold vs Mutated Lead</span>
            <span style="font-size: 11px; color: #94a3b8; margin-left: 8px;">Target: {html.escape(scaffold_label)}</span>
        </div>
        <span class="badge">{html.escape(lead_name)}</span>
    </div>

    <div class="dual-container">
        <!-- Left Viewport: Wild-Type Parent -->
        <div class="viewport-box" id="viewport-wt">
            <div class="viewport-tag" style="color: #60a5fa;">
                ● Wild-Type Parent (PDB Baseline)
            </div>
        </div>

        <!-- Right Viewport: Mutated Lead Candidate -->
        <div class="viewport-box" id="viewport-mutant">
            <div class="viewport-tag" style="color: #34d399;">
                ★ Engineered Lead ({len(effective_mutations)} Mutations)
            </div>
        </div>

        <!-- Float Toolbar -->
        <div class="toolbar">
            <button class="btn" id="btn-spin" onclick="toggleSpin()">Pause Spin</button>
            <button class="btn" onclick="resetViews()">Center Both Views</button>
        </div>
    </div>

    <div class="legend">
        <div><span class="dot" style="background: #38bdf8;"></span>Scaffold Cartoon</div>
        <div><span class="dot" style="background: #fbbf24;"></span>Active Site / Mutated Hotspots</div>
        <div><span class="dot" style="background: #f43f5e;"></span>HETATM Cofactor (PLP/0JO Catalytic Core)</div>
        <div><span class="dot" style="background: #a855f7;"></span>4.5 Å Catalytic Proximity Shell</div>
    </div>

    <script>
    (function() {{
        let viewerWT = null;
        let viewerMutant = null;
        let isSpinning = {spin_js};
        const pdbData = {escaped_pdb};
        const mutations = {mutations_json};
        const catalytic = {catalytic_json};

        function init() {{
            const elWT = document.getElementById("viewport-wt");
            const elMut = document.getElementById("viewport-mutant");

            viewerWT = $3Dmol.createViewer(elWT, {{ backgroundColor: "#070a13" }});
            viewerMutant = $3Dmol.createViewer(elMut, {{ backgroundColor: "#070a13" }});

            // Model 1: Wild-Type
            viewerWT.addModel(pdbData, "pdb");
            viewerWT.setStyle({{}}, {{ cartoon: {{ color: "#2563eb", opacity: 0.82 }} }});
            if (catalytic && catalytic.length > 0) {{
                viewerWT.addStyle({{ resi: catalytic }}, {{ stick: {{ color: "#fbbf24", radius: 0.3 }} }});
            }}
            viewerWT.addStyle({{ hetflag: true }}, {{
                stick: {{ color: "#f43f5e", radius: 0.34 }},
                sphere: {{ color: "#fb7185", radius: 0.45, opacity: 0.95 }}
            }});
            viewerWT.zoomTo();
            viewerWT.render();

            // Model 2: Mutated Lead
            viewerMutant.addModel(pdbData, "pdb");
            viewerMutant.setStyle({{}}, {{ cartoon: {{ color: "#06b6d4", opacity: 0.88 }} }});
            if (mutations && mutations.length > 0) {{
                viewerMutant.addStyle({{ resi: mutations }}, {{ stick: {{ color: "#f59e0b", radius: 0.32 }} }});
            }}
            // Render Cofactor + 4.5 Å proximity shell
            viewerMutant.addStyle({{ hetflag: true }}, {{
                stick: {{ color: "#f43f5e", radius: 0.38 }},
                sphere: {{ color: "#fb7185", radius: 0.50, opacity: 0.95 }}
            }});
            viewerMutant.addStyle({{ within: {{ distance: 4.5, sel: {{ hetflag: true }} }} }}, {{
                stick: {{ color: "#a855f7", radius: 0.22 }}
            }});

            viewerMutant.zoomTo();
            viewerMutant.render();

            if (isSpinning) {{
                viewerWT.spin("y", 0.4);
                viewerMutant.spin("y", 0.4);
            }}
        }}

        window.toggleSpin = function() {{
            isSpinning = !isSpinning;
            const btn = document.getElementById("btn-spin");
            if (isSpinning) {{
                if (viewerWT) viewerWT.spin("y", 0.4);
                if (viewerMutant) viewerMutant.spin("y", 0.4);
                if (btn) btn.innerText = "Pause Spin";
            }} else {{
                if (viewerWT) viewerWT.spin(false);
                if (viewerMutant) viewerMutant.spin(false);
                if (btn) btn.innerText = "Resume Spin";
            }}
        }};

        window.resetViews = function() {{
            if (viewerWT) {{ viewerWT.zoomTo(); viewerWT.render(); }}
            if (viewerMutant) {{ viewerMutant.zoomTo(); viewerMutant.render(); }}
        }};

        init();
    }})();
    </script>
</body>
</html>"""

    escaped_srcdoc = html.escape(inner_html, quote=True)
    return f'<iframe srcdoc="{escaped_srcdoc}" style="width: 100%; height: 550px; border: none; border-radius: 12px; background: #0b1329; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.3);"></iframe>'



def generate_de_novo_complex_viewer_html(
    target_accession: str = "SCFV_9004G_CN2",
    binder_pdb_content: Optional[str] = None,
    hotspot_residues: Optional[List[int]] = None,
    binder_name: str = "VB-DIFF-001 (De Novo 3-Helix Bundle)",
    spin: bool = True,
) -> str:
    """
    Renders an interactive 3D WebGL complex of the target protein epitope docked
    with the de novo designed binder backbone (RFdiffusion + ProteinMPNN).
    """
    from workers.diffusion_worker import GeometricBackboneSampler

    resolved_pdb_path, title_str, scaffold_label, default_hotspots, cat_label, _ = resolve_target_pdb_context(
        target_accession=target_accession
    )
    target_pdb_str = load_pdb_content(resolved_pdb_path)

    # If no binder provided, generate an initial geometric 3-helix bundle for live preview
    if not binder_pdb_content:
        sampler = GeometricBackboneSampler(random_seed=42)
        pdb_b, _, _ = sampler.sample_3helix_bundle(length=65, design_idx=1)
        binder_pdb_content = pdb_b

    # Combine into a single complex PDB: Target as Chain A, Binder as Chain B
    # Filter target lines and ensure Chain A
    target_lines = []
    for line in target_pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if len(line) >= 22:
                line = line[:21] + "A" + line[22:]
            target_lines.append(line)

    binder_lines = []
    for line in binder_pdb_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if len(line) >= 22:
                line = line[:21] + "B" + line[22:]
            binder_lines.append(line)

    complex_pdb = "\n".join(target_lines) + "\nTER\n" + "\n".join(binder_lines) + "\nEND\n"
    escaped_pdb = json.dumps(complex_pdb)

    effective_hotspots = hotspot_residues if hotspot_residues is not None else default_hotspots
    hotspots_json = json.dumps(effective_hotspots)
    spin_js = "true" if spin else "false"

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: #0b1329;
            color: #f8fafc;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            width: 100vw;
            height: 100vh;
            overflow: hidden;
            display: flex;
            flex-direction: column;
        }}
        .header {{
            padding: 10px 16px;
            background: rgba(15, 23, 42, 0.85);
            border-bottom: 1px solid rgba(51, 65, 85, 0.6);
            display: flex;
            justify-content: space-between;
            align-items: center;
            z-index: 10;
        }}
        .title {{ font-size: 13px; font-weight: 600; color: #38bdf8; }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 11px;
            font-weight: 500;
            background: rgba(6, 182, 212, 0.2);
            color: #22d3ee;
            border: 1px solid rgba(6, 182, 212, 0.4);
        }}
        .viewport {{ flex: 1; width: 100%; position: relative; }}
        .controls {{
            position: absolute;
            bottom: 12px;
            left: 12px;
            display: flex;
            gap: 8px;
            z-index: 20;
        }}
        .btn {{
            background: rgba(30, 41, 59, 0.85);
            color: #e2e8f0;
            border: 1px solid #475569;
            padding: 5px 12px;
            border-radius: 6px;
            font-size: 11px;
            cursor: pointer;
            transition: all 0.15s ease;
        }}
        .btn:hover {{ background: #334155; color: #ffffff; }}
        .legend {{
            position: absolute;
            top: 12px;
            right: 12px;
            background: rgba(15, 23, 42, 0.85);
            border: 1px solid rgba(51, 65, 85, 0.6);
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 11px;
            z-index: 20;
            line-height: 1.6;
        }}
        .dot {{ display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; }}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <span class="title">De Novo Docked Complex · RFdiffusion + ProteinMPNN</span>
            <span style="font-size: 11px; color: #94a3b8; margin-left: 8px;">Target: {html.escape(scaffold_label)}</span>
        </div>
        <span class="badge">{html.escape(binder_name)}</span>
    </div>
    <div class="viewport" id="viewport-box">
        <div class="legend">
            <div><span class="dot" style="background: #94a3b8;"></span>Target Epitope (Chain A)</div>
            <div><span class="dot" style="background: #f59e0b;"></span>Target Hotspots (Contacts)</div>
            <div><span class="dot" style="background: #06b6d4;"></span>Designed De Novo Binder (Chain B)</div>
        </div>
        <div class="controls">
            <button class="btn" id="btn-spin" onclick="toggleSpin()">Pause Spin</button>
            <button class="btn" onclick="resetCamera()">Center Camera</button>
        </div>
    </div>
    <script>
    (function() {{
        let glviewer = null;
        let isSpinning = {spin_js};
        const pdbData = {escaped_pdb};
        const hotspots = {hotspots_json};

        function init() {{
            const el = document.getElementById("viewport-box");
            glviewer = $3Dmol.createViewer(el, {{ backgroundColor: "#0b1329" }});
            glviewer.addModel(pdbData, "pdb");

            // 1. Style Target Protein (Chain A): Translucent slate gray
            glviewer.setStyle({{ chain: "A" }}, {{
                cartoon: {{ color: "#94a3b8", opacity: 0.70 }}
            }});

            // 2. Highlight Hotspots on Chain A: Golden yellow sticks
            if (hotspots && hotspots.length > 0) {{
                glviewer.addStyle({{ chain: "A", resi: hotspots }}, {{
                    stick: {{ color: "#f59e0b", radius: 0.28 }}
                }});
            }}

            // 3. Style Designed De Novo Binder (Chain B): Bright electric cyan cartoon
            glviewer.setStyle({{ chain: "B" }}, {{
                cartoon: {{ color: "#06b6d4", opacity: 1.0, thickness: 0.6 }}
            }});

            // Show interface sidechains on Binder
            glviewer.addStyle({{ chain: "B", within: {{ distance: 8.5, sel: {{ chain: "A" }} }} }}, {{
                stick: {{ color: "#22d3ee", radius: 0.20 }}
            }});

            glviewer.zoomTo();
            glviewer.render();
            if (isSpinning) glviewer.spin("y", 0.5);
        }}

        window.toggleSpin = function() {{
            if (!glviewer) return;
            isSpinning = !isSpinning;
            const btn = document.getElementById("btn-spin");
            if (isSpinning) {{
                glviewer.spin("y", 0.5);
                if (btn) btn.innerText = "Pause Spin";
            }} else {{
                glviewer.spin(false);
                if (btn) btn.innerText = "Resume Spin";
            }}
        }};

        window.resetCamera = function() {{
            if (!glviewer) return;
            glviewer.zoomTo();
            glviewer.render();
        }};

        init();
    }})();
    </script>
</body>
</html>"""

    escaped_srcdoc = html.escape(inner_html, quote=True)
    return f'<iframe srcdoc="{escaped_srcdoc}" style="width: 100%; height: 540px; border: none; border-radius: 12px; background: #0b1329; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.3);"></iframe>'


def generate_dual_comparison_viewer_html(
    target_accession: str = "TRPB_TM9D8",
    mutant_positions: Optional[List[int]] = None,
    lead_name: str = "Epistatic Lead Candidate",
    viewer_id: str = "veritasbio-dual-viewer",
    spin: bool = True,
    pdb_path: Optional[Path] = None,
) -> str:
    """
    Builds a synchronized dual-viewport 3D molecular viewer (WT Parent vs Optimized Lead)
    with explicit HETATM cofactor catalytic coordination shell rendering.
    """
    (
        resolved_pdb_path,
        title_str,
        scaffold_label,
        catalytic_residues,
        catalytic_label,
        scaffold_color,
    ) = resolve_target_pdb_context(target_accession=target_accession, pdb_path=pdb_path)

    pdb_data = load_pdb_content(resolved_pdb_path)
    escaped_pdb = json.dumps(pdb_data)

    effective_mutations = mutant_positions or ([183, 184, 227, 228] if "8EH0" in str(resolved_pdb_path) else [65, 66, 163, 203])
    mutations_json = json.dumps(effective_mutations)
    catalytic_json = json.dumps(catalytic_residues)
    spin_js = "true" if spin else "false"

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0; padding: 0; background: #0b1329; font-family: 'Inter', sans-serif;
            overflow: hidden; user-select: none; color: #f8fafc;
        }}
        .header {{
            display: flex; justify-content: space-between; align-items: center;
            padding: 8px 16px; background: #0f172a; border-bottom: 1px solid #1e293b;
        }}
        .title {{ font-size: 13px; font-weight: 600; color: #38bdf8; }}
        .badge {{
            background: rgba(56, 189, 248, 0.15); border: 1px solid rgba(56, 189, 248, 0.3);
            color: #38bdf8; padding: 2px 8px; border-radius: 4px; font-size: 11px;
            font-family: 'JetBrains Mono', monospace;
        }}
        .dual-container {{
            display: grid; grid-template-columns: 1fr 1fr; width: 100%; height: 450px;
            position: relative; background: #070a13;
        }}
        .viewport-box {{
            width: 100%; height: 100%; position: relative;
            border-right: 1px solid rgba(255,255,255,0.06);
        }}
        .viewport-tag {{
            position: absolute; top: 10px; left: 12px; z-index: 10;
            background: rgba(15, 23, 42, 0.85); border: 1px solid rgba(51, 65, 85, 0.6);
            padding: 4px 10px; border-radius: 6px; font-size: 11px; font-weight: 600;
        }}
        .toolbar {{
            position: absolute; bottom: 10px; left: 50%; transform: translateX(-50%);
            z-index: 20; display: flex; gap: 8px;
        }}
        .btn {{
            background: #1e293b; color: #cbd5e1; border: 1px solid #334155;
            padding: 5px 12px; border-radius: 6px; cursor: pointer; font-size: 11px;
            font-weight: 500; transition: all 0.15s ease;
        }}
        .btn:hover {{ background: #334155; color: #ffffff; }}
        .legend {{
            display: flex; justify-content: space-around; align-items: center;
            padding: 8px 16px; background: #0b1329; border-top: 1px solid #1e293b;
            font-family: 'JetBrains Mono', monospace; font-size: 11px; color: #94a3b8;
        }}
        .dot {{ display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; }}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <span class="title">Dual-Structure Comparison · WT Scaffold vs Mutated Lead</span>
            <span style="font-size: 11px; color: #94a3b8; margin-left: 8px;">Target: {html.escape(scaffold_label)}</span>
        </div>
        <span class="badge">{html.escape(lead_name)}</span>
    </div>

    <div class="dual-container">
        <!-- Left Viewport: Wild-Type Parent -->
        <div class="viewport-box" id="viewport-wt">
            <div class="viewport-tag" style="color: #60a5fa;">
                ● Wild-Type Parent (PDB Baseline)
            </div>
        </div>

        <!-- Right Viewport: Mutated Lead Candidate -->
        <div class="viewport-box" id="viewport-mutant">
            <div class="viewport-tag" style="color: #34d399;">
                ★ Engineered Lead ({len(effective_mutations)} Mutations)
            </div>
        </div>

        <!-- Float Toolbar -->
        <div class="toolbar">
            <button class="btn" id="btn-spin" onclick="toggleSpin()">Pause Spin</button>
            <button class="btn" onclick="resetViews()">Center Both Views</button>
        </div>
    </div>

    <div class="legend">
        <div><span class="dot" style="background: #38bdf8;"></span>Scaffold Cartoon</div>
        <div><span class="dot" style="background: #fbbf24;"></span>Active Site / Mutated Hotspots</div>
        <div><span class="dot" style="background: #f43f5e;"></span>HETATM Cofactor (PLP/0JO Catalytic Core)</div>
        <div><span class="dot" style="background: #a855f7;"></span>4.5 Å Catalytic Proximity Shell</div>
    </div>

    <script>
    (function() {{
        let viewerWT = null;
        let viewerMutant = null;
        let isSpinning = {spin_js};
        const pdbData = {escaped_pdb};
        const mutations = {mutations_json};
        const catalytic = {catalytic_json};

        function init() {{
            const elWT = document.getElementById("viewport-wt");
            const elMut = document.getElementById("viewport-mutant");

            viewerWT = $3Dmol.createViewer(elWT, {{ backgroundColor: "#070a13" }});
            viewerMutant = $3Dmol.createViewer(elMut, {{ backgroundColor: "#070a13" }});

            // Model 1: Wild-Type
            viewerWT.addModel(pdbData, "pdb");
            viewerWT.setStyle({{}}, {{ cartoon: {{ color: "#2563eb", opacity: 0.82 }} }});
            if (catalytic && catalytic.length > 0) {{
                viewerWT.addStyle({{ resi: catalytic }}, {{ stick: {{ color: "#fbbf24", radius: 0.3 }} }});
            }}
            viewerWT.addStyle({{ hetflag: true }}, {{
                stick: {{ color: "#f43f5e", radius: 0.34 }},
                sphere: {{ color: "#fb7185", radius: 0.45, opacity: 0.95 }}
            }});
            viewerWT.zoomTo();
            viewerWT.render();

            // Model 2: Mutated Lead
            viewerMutant.addModel(pdbData, "pdb");
            viewerMutant.setStyle({{}}, {{ cartoon: {{ color: "#06b6d4", opacity: 0.88 }} }});
            if (mutations && mutations.length > 0) {{
                viewerMutant.addStyle({{ resi: mutations }}, {{ stick: {{ color: "#f59e0b", radius: 0.32 }} }});
            }}
            // Render Cofactor + 4.5 Å proximity shell
            viewerMutant.addStyle({{ hetflag: true }}, {{
                stick: {{ color: "#f43f5e", radius: 0.38 }},
                sphere: {{ color: "#fb7185", radius: 0.50, opacity: 0.95 }}
            }});
            viewerMutant.addStyle({{ within: {{ distance: 4.5, sel: {{ hetflag: true }} }} }}, {{
                stick: {{ color: "#a855f7", radius: 0.22 }}
            }});

            viewerMutant.zoomTo();
            viewerMutant.render();

            if (isSpinning) {{
                viewerWT.spin("y", 0.4);
                viewerMutant.spin("y", 0.4);
            }}
        }}

        window.toggleSpin = function() {{
            isSpinning = !isSpinning;
            const btn = document.getElementById("btn-spin");
            if (isSpinning) {{
                if (viewerWT) viewerWT.spin("y", 0.4);
                if (viewerMutant) viewerMutant.spin("y", 0.4);
                if (btn) btn.innerText = "Pause Spin";
            }} else {{
                if (viewerWT) viewerWT.spin(false);
                if (viewerMutant) viewerMutant.spin(false);
                if (btn) btn.innerText = "Resume Spin";
            }}
        }};

        window.resetViews = function() {{
            if (viewerWT) {{ viewerWT.zoomTo(); viewerWT.render(); }}
            if (viewerMutant) {{ viewerMutant.zoomTo(); viewerMutant.render(); }}
        }};

        init();
    }})();
    </script>
</body>
</html>"""

    escaped_srcdoc = html.escape(inner_html, quote=True)
    return f'<iframe srcdoc="{escaped_srcdoc}" style="width: 100%; height: 550px; border: none; border-radius: 12px; background: #0b1329; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.3);"></iframe>'


