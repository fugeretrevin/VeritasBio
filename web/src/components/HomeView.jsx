import React, { useState, useEffect, useRef } from 'react';
import {
  ArrowRight,
  ChevronLeft,
  ChevronRight,
  Play,
  Check,
  Sparkles,
  ShieldCheck,
  Layers,
  Cpu,
  FlaskConical,
} from 'lucide-react';

// =============================================================================
// 1. Curated Multi-Target Showcase for Automated Hero 3D Molecular Choreography
// =============================================================================
const SHOWCASE_PROTEINS = [
  {
    id: 'TRPB_TM9D8',
    pdbId: '8EH0',
    shortLabel: 'TrpB Synthase',
    name: 'Tryptophan Synthase β-Subunit (Tm9D8*)',
    category: '4-Site Combinatorial Enzyme · 159,999 Variants',
    hotspots: [183, 184, 227, 228],
    wtResidues: 'V183 · F184 · V227 · S228',
    leadMutation: 'V183A : F184I : V227K : S228G',
    wtMetric: '0.402 WT Baseline',
    leadMetric: '1.0000 Global Peak (98% Seed Hit Rate)',
    summary:
      'Captures the K227 × G228 sign-epistasis flip across a 20⁴ active-site library while protecting the PLP catalytic cofactor.',
    backboneColor: '#64748b',
    optimizedBackboneColor: '#1e293b',
    hotspotColor: '#0d9488',
    cofactorColor: '#0284c7',
  },
  {
    id: 'SCFV_9004G_CN2',
    pdbId: '2YBR',
    shortLabel: 'Antivenom scFv',
    name: 'Neutralizing Antibody 9004G · Toxin Cn2 Complex',
    category: 'Immunotherapy Biologic · Interface Maturation',
    hotspots: [98, 99, 100, 101, 102],
    wtResidues: 'CDR-H3 Paratope Loop (Res 98–102)',
    leadMutation: 'Y100W : S101A : G102V',
    wtMetric: 'Kd = 1.8 nM (Parent scFv)',
    leadMetric: 'Kd = 95 pM (19x Affinity Gain)',
    summary:
      'Repacks the heavy-chain CDR-H3 loop against Centruroides noxius neurotoxin Cn2 while eliminating deamidation liabilities.',
    backboneColor: '#64748b',
    optimizedBackboneColor: '#0f172a',
    hotspotColor: '#0d9488',
    cofactorColor: '#0369a1',
  },
  {
    id: 'P42212',
    pdbId: '1EMA',
    shortLabel: 'avGFP Barrel',
    name: 'Aequorea victoria Green Fluorescent Protein',
    category: 'Optical Biosensor · 54,025 DMS Landscape',
    hotspots: [64, 65, 68, 203],
    wtResidues: 'F64 · S65 · V68 · T203',
    leadMutation: 'F64L : S65T : V68L : T203Y',
    wtMetric: '1.0x Relative Brightness',
    leadMetric: '3.2x Quantum Yield (ρ = 0.85)',
    summary:
      'Stabilizes the central α-helix chromophore hydrogen-bond network inside the 11-stranded β-barrel.',
    backboneColor: '#64748b',
    optimizedBackboneColor: '#1e293b',
    hotspotColor: '#059669',
    cofactorColor: '#0284c7',
  },
  {
    id: 'HER2_TRASTUZUMAB',
    pdbId: '1N8Z',
    shortLabel: 'Trastuzumab Fab',
    name: 'Trastuzumab Fab · HER2 Extracellular Domain',
    category: 'Oncology mAb · Developability Co-Optimization',
    hotspots: [99, 100, 101, 102],
    wtResidues: 'W99 · G100 · D101 · G102',
    leadMutation: 'W99F : G100A : D101E : G102Y',
    wtMetric: 'DG101 Isomerization Liability',
    leadMetric: '+1.45 ΔΔG · Zero CMC Liabilities',
    summary:
      'Simultaneously optimizes HER2 binding energy and removes the clinical D101 isomerization hotspot.',
    backboneColor: '#64748b',
    optimizedBackboneColor: '#0f172a',
    hotspotColor: '#0d9488',
    cofactorColor: '#0369a1',
  },
];

// =============================================================================
// 2. Continuous Capability Rail Items (Clean SVG Diagrams + Tonal Cards)
// =============================================================================
const CAPABILITY_CARDS = [
  {
    title: 'Calibrated Ensemble Surrogate',
    metric: '98% Seed Reliability',
    desc: 'Combines amino-acid physicochemical descriptors with pair-interaction priors so epistasis is learned from the very first 96-well plate.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <path d="M 12 50 Q 60 46, 96 28 T 174 16 T 228 10" fill="none" stroke="#0d9488" strokeWidth="2.2" />
        <path d="M 12 50 Q 60 46, 96 28 T 174 16 T 228 10 L 228 56 L 12 56 Z" fill="rgba(13, 148, 136, 0.1)" />
        <circle cx="52" cy="43" r="3.5" fill="#141413" />
        <circle cx="116" cy="23" r="3.5" fill="#141413" />
        <circle cx="178" cy="15" r="4" fill="#0d9488" />
        <circle cx="220" cy="11" r="5" fill="#0d9488" />
      </svg>
    ),
  },
  {
    title: 'Qualify-Then-Optimize Frontier',
    metric: '38% Global Peak Hit',
    desc: 'Reserves 3 exploratory slots per batch to test unseen mutations at epistatic positions, unlocking peaks like K227 × G228 without hurting reliability.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <rect x="14" y="14" width="90" height="36" rx="6" fill="#efeae1" />
        <text x="59" y="36" textAnchor="middle" fontSize="10" fontWeight="600" fill="#141413" fontFamily="monospace">21 EXPLOIT</text>
        <path d="M 108 32 L 132 32" stroke="#0d9488" strokeWidth="1.8" />
        <rect x="132" y="14" width="94" height="36" rx="6" fill="rgba(13, 148, 136, 0.14)" />
        <text x="179" y="36" textAnchor="middle" fontSize="10" fontWeight="700" fill="#0f766e" fontFamily="monospace">+3 FRONTIER</text>
      </svg>
    ),
  },
  {
    title: 'CPU Sidechain Clash Rescue',
    metric: '< 50 ms L-BFGS',
    desc: 'Relaxes mutated rotamers in local 3D space before scoring, rescuing high-affinity bulky substitutions that rigid distance filters falsely reject.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <circle cx="56" cy="32" r="16" fill="rgba(220, 38, 38, 0.12)" />
        <circle cx="72" cy="32" r="16" fill="rgba(220, 38, 38, 0.12)" />
        <path d="M 102 32 L 138 32" stroke="#141413" strokeWidth="1.5" strokeDasharray="3 3" />
        <circle cx="168" cy="32" r="16" fill="rgba(13, 148, 136, 0.15)" />
        <circle cx="204" cy="32" r="16" fill="rgba(13, 148, 136, 0.15)" />
      </svg>
    ),
  },
  {
    title: '4D Biologics Pareto Surface',
    metric: 'NSGA-II Rank 0',
    desc: 'Co-optimizes binding fitness, CMC developability (deamidation, isomerization, aggregation), codon synthesizability, and fold stability.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <path d="M 28 50 Q 80 18, 210 14" fill="none" stroke="#0d9488" strokeWidth="2" strokeDasharray="4 3" />
        <circle cx="48" cy="40" r="4.5" fill="#0d9488" />
        <circle cx="102" cy="24" r="4.5" fill="#0d9488" />
        <circle cx="168" cy="16" r="4.5" fill="#0d9488" />
        <circle cx="95" cy="44" r="3.5" fill="#9c9589" />
        <circle cx="148" cy="38" r="3.5" fill="#9c9589" />
      </svg>
    ),
  },
  {
    title: 'Fail-Closed Micro-USD Ledger',
    metric: '1 USD = 10⁶ µUSD',
    desc: 'Every paid Modal A10G tensor pass must clear an integer micro-USD SQLite gatekeeper stamped with an HMAC-SHA256 Merkle audit block.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <rect x="18" y="16" width="56" height="32" rx="5" fill="#efeae1" />
        <rect x="92" y="16" width="56" height="32" rx="5" fill="#efeae1" />
        <rect x="166" y="16" width="56" height="32" rx="5" fill="#0d9488" />
        <path d="M 74 32 L 92 32 M 148 32 L 166 32" stroke="#141413" strokeWidth="1.5" />
      </svg>
    ),
  },
  {
    title: 'Zero-Touch CRO Handshake',
    metric: 'Twist + Opentrons',
    desc: 'Compiles codon-optimized Golden Gate gene fragments, perimeter-hydrated 96-well Opentrons OT-2 scripts, and 4PL plate reader ingestion.',
    renderGraphic: () => (
      <svg viewBox="0 0 240 64" style={{ width: '100%', height: '56px', display: 'block' }}>
        <g fill="#0d9488">
          <circle cx="44" cy="20" r="3.5" /><circle cx="64" cy="20" r="3.5" /><circle cx="84" cy="20" r="3.5" /><circle cx="104" cy="20" r="3.5" />
          <circle cx="44" cy="34" r="3.5" /><circle cx="64" cy="34" r="3.5" /><circle cx="84" cy="34" r="3.5" /><circle cx="104" cy="34" r="3.5" />
          <circle cx="44" cy="48" r="3.5" /><circle cx="64" cy="48" r="3.5" /><circle cx="84" cy="48" r="3.5" /><circle cx="104" cy="48" r="3.5" />
        </g>
        <path d="M 124 34 L 154 34" stroke="#141413" strokeWidth="1.5" />
        <rect x="158" y="16" width="62" height="36" rx="6" fill="#141413" />
        <text x="189" y="38" textAnchor="middle" fontSize="9.5" fontWeight="600" fill="#fcfbf8" fontFamily="monospace">CSV + .PY</text>
      </svg>
    ),
  },
];

// =============================================================================
// 3. Sign-Epistasis Matrix Data (Real TrpB 227 × 228 Data from 48-Seed Study)
// =============================================================================
const EPISTASIS_CELLS = [
  {
    id: 'M227_S228',
    pos227: 'M227 (Local Trap)',
    pos228: 'S228 (Wild-Type)',
    genotype: 'V183I : F184V : V227M : S228S',
    meanFitness: '0.284',
    peakFitness: '0.512',
    role: 'Permissive Single-Mutant Step',
    tagColor: '#52504a',
    bgTone: 'var(--bg-dark)',
    explanation:
      'In the wild-type S228 background, M227 is well-tolerated (+0.11 over WT). Greedy and standard GP acquisition functions heavily sample M227 early because it looks safe.',
  },
  {
    id: 'M227_G228',
    pos227: 'M227 (Local Trap)',
    pos228: 'G228 (Glycine Shift)',
    genotype: 'V183I : F184V : V227M : S228G',
    meanFitness: '0.391',
    peakFitness: '0.619',
    role: 'Deceptive Local Optimum (#96)',
    tagColor: '#b45309',
    bgTone: 'rgba(180, 83, 9, 0.08)',
    explanation:
      'Pairing M227 with G228 reaches fitness 0.619 (rank #96 out of 159,999). Standard active-learning engines converge here and stop, missing the true global peak.',
  },
  {
    id: 'K227_S228',
    pos227: 'K227 (Lysine)',
    pos228: 'S228 (Wild-Type)',
    genotype: 'V183A : F184I : V227K : S228S',
    meanFitness: '0.031',
    peakFitness: '0.094',
    role: 'Deleterious in WT Background (-0.50)',
    tagColor: '#dc2626',
    bgTone: 'rgba(220, 38, 38, 0.06)',
    explanation:
      'Without G228, introducing the positively charged Lysine K227 sterically clashes with S228 (-0.50 fitness penalty). Additive models permanently blacklist K227.',
  },
  {
    id: 'K227_G228',
    pos227: 'K227 (Lysine)',
    pos228: 'G228 (Glycine Shift)',
    genotype: 'V183A : F184I : V227K : S228G',
    meanFitness: '0.842',
    peakFitness: '1.0000',
    role: 'True Global Optimum (#1 / 159,999)',
    tagColor: '#0d9488',
    bgTone: 'rgba(13, 148, 136, 0.12)',
    explanation:
      'Once G228 clears steric volume in the catalytic pocket, K227 flips from -0.50 to +0.82! PanoptesBio’s qualify-then-optimize rule probes K227 inside the G228 background and hits 1.0000.',
  },
];

// =============================================================================
// 4. Multi-Dataset Empirical Summaries (Streamlined, No Intermediate Clutter)
// =============================================================================
const DATASET_TARGETS = [
  {
    id: 'trpb',
    name: 'TrpB Synthase (Tm9D8*)',
    citation: 'Johnston et al. 2024 · PDB 8EH0 · 159,999 Combinatorial Variants',
    wtLabel: 'Wild-Type Parent (V183 · F184 · V227 · S228)',
    wtScore: '0.402',
    leadLabel: 'V183A : F184I : V227K : S228G',
    leadScore: '1.0000 (Global Peak #1)',
    hitRate: '98% of seeds hit ≥ 0.70 (vs 58% ALDE)',
    halfBudget: '88% hit rate in 50 assays (vs 4% ALDE)',
    spend: '$0.0094',
    takeaway:
      'Across 48 paired random starts (24 initial + 3 × 24 active learning = 96 measurements out of 159,999), PanoptesBio reaches the ≥ 0.70 top tier on 47/48 seeds and cuts median assays to hit by 2x.',
  },
  {
    id: 'gfp',
    name: 'avGFP Fluorophore',
    citation: 'Sarkisyan et al. 2016 · PDB 1EMA · 54,025 Empirical Variants',
    wtLabel: 'Wild-Type Chromophore (F64 · S65 · V68 · T203)',
    wtScore: '1.00x',
    leadLabel: 'F64L : S65T : V68L : T203Y',
    leadScore: '3.20x Brightness (ρ = 0.85)',
    hitRate: 'Top 0.05% superfolder variant in 72 assays',
    halfBudget: 'Zero chromophore-quenching false positives',
    spend: '$0.0125',
    takeaway:
      'Maintains out-of-sample Spearman ρ = 0.85 across multi-mutant epistasis while protecting the buried CRO66 hydrogen-bond network.',
  },
  {
    id: 'gb1',
    name: 'Streptococcal GB1 Domain',
    citation: 'Wu et al. 2016 · PDB 2QMT · 149,361 Combinatorial Variants',
    wtLabel: 'Wild-Type IgG Fc Interface (V39 · D40 · G41 · V54)',
    wtScore: '1.00',
    leadLabel: 'V39W : D40E : G41A : V54I',
    leadScore: '5.84 Relative Binding Fitness',
    hitRate: '9.2x enrichment over random screening',
    halfBudget: 'Rescued V39W aromatic packing via CPU L-BFGS',
    spend: '$0.0062',
    takeaway:
      'Demonstrates why local sidechain relaxation matters: rigid steric filters reject V39W, whereas <50ms CPU relaxation resolves the rotamer cone and captures the top Fc binder.',
  },
  {
    id: 'her2',
    name: 'Trastuzumab CDR-H3 (HER2)',
    citation: 'PDB 1N8Z · Therapeutic Antibody Paratope Maturation',
    wtLabel: 'Clinical Parent (Contains D101-G102 Isomerization Motif)',
    wtScore: '0.42',
    leadLabel: 'W99F : G100A : D101E : G102Y',
    leadScore: '+1.45 ΔΔG · 100% CMC Clean',
    hitRate: 'Eliminates DG isomerization & NG deamidation',
    halfBudget: 'Rank 0 on 4-objective Pareto frontier',
    spend: '$0.0064',
    takeaway:
      'Simultaneously improves predicted HER2 binding affinity and removes the clinical D101 chemical degradation liability before ordering a single cDNA fragment.',
  },
];

// =============================================================================
// 5. Right-Side Hover Dot Navigation Sections
// =============================================================================
const NAV_SECTIONS = [
  { id: 'mission', label: 'Overview' },
  { id: 'capabilities', label: 'Capabilities' },
  { id: 'benchmarks', label: 'Economics' },
  { id: 'datasets', label: '48-Seed Proof' },
  { id: 'how-it-works', label: 'How It Works' },
  { id: 'how-to-use', label: 'Quick Start' },
  { id: 'next-steps', label: 'Launch Studio' },
];

export default function HomeView({ onNavigate, onOpenAuth }) {
  // Hero 3D Showcase State
  const [proteinIdx, setProteinIdx] = useState(0);
  const [animPhase, setAnimPhase] = useState(0); // 0: WT, 1: Active Site, 2: Lead
  const [pdbText, setPdbText] = useState(null);
  const viewerContainerRef = useRef(null);
  const viewerInstanceRef = useRef(null);

  // Capabilities Marquee Slow-on-Hover Ref
  const marqueeTrackRef = useRef(null);

  // Proof Section Mode ('reliability' | 'epistasis' | 'datasets')
  const [proofTab, setProofTab] = useState('reliability');
  const [selectedEpistasisId, setSelectedEpistasisId] = useState('K227_G228');
  const [datasetIdx, setDatasetIdx] = useState(0);

  // Interactive Pipeline Stage (0, 1, 2)
  const [activePipelineStage, setActivePipelineStage] = useState(1);

  // Right Dot Rail Active Section & Hover State
  const [activeSection, setActiveSection] = useState('mission');
  const [hoveredDot, setHoveredDot] = useState(null);

  const activeProtein = SHOWCASE_PROTEINS[proteinIdx];
  const activeDataset = DATASET_TARGETS[datasetIdx];
  const activeEpistasis =
    EPISTASIS_CELLS.find((c) => c.id === selectedEpistasisId) || EPISTASIS_CELLS[3];

  // Track scroll position for right-side navigation dots
  useEffect(() => {
    const handleScroll = () => {
      const scrollY = window.scrollY + 260;
      let current = 'mission';
      for (const sec of NAV_SECTIONS) {
        const el = document.getElementById(sec.id);
        if (el && el.offsetTop <= scrollY) {
          current = sec.id;
        }
      }
      setActiveSection(current);
    };
    window.addEventListener('scroll', handleScroll, { passive: true });
    handleScroll();
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  // Automatically cycle 3D choreography phases (WT -> Pocket -> Optimized Lead -> Next Protein)
  useEffect(() => {
    const phaseTimer = setInterval(() => {
      setAnimPhase((prev) => {
        if (prev < 2) return prev + 1;
        setProteinIdx((pIdx) => (pIdx + 1) % SHOWCASE_PROTEINS.length);
        return 0;
      });
    }, 4600);
    return () => clearInterval(phaseTimer);
  }, []);

  // Load PDB structure when proteinIdx changes
  useEffect(() => {
    let canceled = false;
    async function fetchPdb() {
      try {
        const res = await fetch(`/api/structure/${activeProtein.id}`);
        if (res.ok) {
          const data = await res.json();
          if (!canceled && data.pdb_content) {
            setPdbText(data.pdb_content);
            return;
          }
        }
      } catch (err) {
        // Fallback below
      }
      try {
        const cdnRes = await fetch(`https://files.rcsb.org/download/${activeProtein.pdbId}.pdb`);
        if (cdnRes.ok && !canceled) {
          setPdbText(await cdnRes.text());
        }
      } catch (err) {
        console.warn('PDB fetch fallback failed:', err);
      }
    }
    fetchPdb();
    return () => {
      canceled = true;
    };
  }, [activeProtein.id, activeProtein.pdbId]);

  // Render & Choreograph 3Dmol Viewer (Non-interactive automated camera)
  useEffect(() => {
    if (!viewerContainerRef.current || !pdbText || typeof window.$3Dmol === 'undefined') return;

    try {
      viewerContainerRef.current.innerHTML = '';
      const viewer = window.$3Dmol.createViewer(viewerContainerRef.current, {
        backgroundColor: '#fcfbf8',
        antialias: true,
      });
      viewerInstanceRef.current = viewer;
      viewer.addModel(pdbText, 'pdb');

      if (animPhase === 0) {
        viewer.setStyle({}, { cartoon: { color: activeProtein.backboneColor, opacity: 0.88 } });
        viewer.addStyle(
          { hetflag: true },
          { stick: { color: activeProtein.cofactorColor, radius: 0.32 } }
        );
        viewer.zoomTo();
      } else if (animPhase === 1) {
        viewer.setStyle({}, { cartoon: { color: activeProtein.backboneColor, opacity: 0.45 } });
        viewer.addStyle(
          { resi: activeProtein.hotspots },
          {
            stick: { color: activeProtein.cofactorColor, radius: 0.36 },
            sphere: { color: activeProtein.cofactorColor, radius: 0.42, opacity: 0.85 },
          }
        );
        viewer.addStyle(
          { hetflag: true },
          {
            stick: { color: '#141413', radius: 0.35 },
            sphere: { color: activeProtein.cofactorColor, radius: 0.48, opacity: 0.9 },
          }
        );
        viewer.zoomTo({ resi: activeProtein.hotspots });
        viewer.zoom(0.82);
      } else {
        viewer.setStyle(
          {},
          { cartoon: { color: activeProtein.optimizedBackboneColor, opacity: 0.92 } }
        );
        viewer.addStyle(
          { resi: activeProtein.hotspots },
          {
            stick: { color: activeProtein.hotspotColor, radius: 0.4 },
            sphere: { color: '#10b981', radius: 0.46, opacity: 0.92 },
          }
        );
        viewer.addStyle(
          { within: { distance: 4.5, sel: { resi: activeProtein.hotspots } } },
          { stick: { color: '#64748b', radius: 0.2 } }
        );
        viewer.addStyle(
          { hetflag: true },
          {
            stick: { color: activeProtein.cofactorColor, radius: 0.36 },
            sphere: { color: activeProtein.cofactorColor, radius: 0.48, opacity: 0.94 },
          }
        );
        viewer.zoomTo();
      }

      viewer.render();
      viewer.spin('y', 0.45);
    } catch (err) {
      console.error('Hero 3Dmol render error:', err);
    }

    return () => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.spin(false);
      }
    };
  }, [pdbText, animPhase, activeProtein]);

  const scrollToSection = (id) => {
    const el = document.getElementById(id);
    if (el) {
      el.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  };

  const handleMarqueeSpeed = (rate) => {
    if (!marqueeTrackRef.current) return;
    const animations = marqueeTrackRef.current.getAnimations();
    animations.forEach((anim) => {
      anim.playbackRate = rate;
    });
  };

  const pipelineStages = [
    {
      num: 'Stage I',
      title: 'Structural & Prior Grounding',
      subtitle: 'PDB Coordinate Masking + Meta Llama 3.1 Boundary JSON',
      inputLabel: 'Target PDB (e.g., 8EH0) + Cofactor Coordinates (PLP / Heme)',
      operationLabel:
        'Builds an 8.0 Å Cα contact graph, locks catalytic residues within 4.2 Å of the cofactor, and prompts Llama 3.1 in strict JSON mode to propose mutable residue shells.',
      outputLabel: 'Validated Combinatorial Search Space (159,999 Candidates · 0 Active-Site Clashes)',
    },
    {
      num: 'Stage II',
      title: 'Calibrated Ensemble Surrogate & Epistasis Screen',
      subtitle: 'Pair-Interaction GP + Qualify-Then-Optimize Frontier + CPU L-BFGS',
      inputLabel: '159,999 Candidate Genotypes + Prior Wet-Lab Plate Readings',
      operationLabel:
        'Fits a calibrated 5-model ridge ensemble with explicit AA×AA pair features, reserves 3 frontier slots per batch to test unseen mutations at epistatic sites, and relaxes sidechain clashes on CPU in <50ms.',
      outputLabel: 'Pareto Rank-0 Micro-Batch (24 High-Information Variants · $0.000 GPU Screening Cost)',
    },
    {
      num: 'Stage III',
      title: 'Micro-USD Spend Gate & Zero-Touch CRO Export',
      subtitle: 'SQLite Immutability Gate + Modal A10G Tensor Pass + Wet-Lab Files',
      inputLabel: '24 Selected Variants + Estimated Compute Cost ($0.0031 / 3,100 µUSD)',
      operationLabel:
        'Atomically verifies remaining budget in SQLite, scores batch on Modal A10G ESM-2 650M, stamps an HMAC-SHA256 Merkle block, and compiles Twist Golden Gate CSV + Opentrons OT-2 scripts.',
      outputLabel: '96-Well Synthesis Order + Robot Protocol + 21 CFR Part 11 Audit Proof',
    },
  ];

  const currentStage = pipelineStages[activePipelineStage];

  return (
    <div style={{ position: 'relative', display: 'flex', flexDirection: 'column', gap: '136px' }}>
      {/* =====================================================================
          RIGHT-SIDE FLOATING DOT RAIL (Dots only; label pops out on hover)
          ===================================================================== */}
      <nav
        aria-label="Section navigation"
        style={{
          position: 'fixed',
          right: '22px',
          top: '50%',
          transform: 'translateY(-50%)',
          zIndex: 900,
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'flex-end',
          gap: '12px',
        }}
      >
        {NAV_SECTIONS.map((sec) => {
          const isActive = activeSection === sec.id;
          const isHovered = hoveredDot === sec.id;
          return (
            <button
              key={sec.id}
              type="button"
              onClick={() => scrollToSection(sec.id)}
              onMouseEnter={() => setHoveredDot(sec.id)}
              onMouseLeave={() => setHoveredDot(null)}
              aria-label={sec.label}
              style={{
                position: 'relative',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'flex-end',
                background: 'none',
                border: 'none',
                padding: '4px',
                cursor: 'pointer',
              }}
            >
              {isHovered && (
                <span
                  style={{
                    position: 'absolute',
                    right: '24px',
                    whiteSpace: 'nowrap',
                    padding: '4px 10px',
                    borderRadius: '6px',
                    background: '#141413',
                    color: '#fcfbf8',
                    fontSize: '12px',
                    fontWeight: 500,
                    pointerEvents: 'none',
                  }}
                >
                  {sec.label}
                </span>
              )}
              <span
                style={{
                  width: isActive ? '10px' : '7px',
                  height: isActive ? '10px' : '7px',
                  borderRadius: '50%',
                  background: isActive ? '#0d9488' : '#b8b2a6',
                  transition: 'all 0.2s ease',
                  display: 'block',
                }}
              />
            </button>
          );
        })}
      </nav>

      {/* =====================================================================
          SECTION 1: HERO (#mission) — Open Editorial Split + Automated 3D Showcase
          ===================================================================== */}
      <section
        id="mission"
        style={{
          scrollMarginTop: '96px',
          paddingTop: '16px',
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(360px, 1fr))',
          gap: '56px',
          alignItems: 'center',
        }}
      >
        {/* Left: Editorial Statement + Actions + Unboxed Metric Strip */}
        <div>
          <h1
            className="title-display"
            style={{
              fontSize: 'var(--step-hero)',
              color: 'var(--text-main)',
              marginBottom: '22px',
            }}
          >
            Direct protein evolution,{' '}
            <span style={{ color: '#0d9488' }}>
              without burning
            </span>{' '}
            wet-lab capital.
          </h1>

          <p
            style={{
              fontSize: 'var(--step-lead)',
              color: 'var(--text-muted)',
              lineHeight: 1.55,
              maxWidth: '58ch',
              marginBottom: '32px',
            }}
          >
            PanoptesBio pairs 3D structural active learning with a fail-closed micro-USD spend ledger, compressing 160,000-variant combinatorial landscapes into a single 96-well synthesis plate.
          </p>

          <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap', marginBottom: '44px' }}>
            <button
              type="button"
              onClick={() => onNavigate && onNavigate('studio')}
              className="btn btn-primary"
              style={{ padding: '13px 24px', fontSize: '15px', fontWeight: 600 }}
            >
              <Play size={15} />
              <span>Launch Studio</span>
            </button>

            <button
              type="button"
              onClick={() => scrollToSection('datasets')}
              className="btn btn-secondary"
              style={{ padding: '13px 22px', fontSize: '15px' }}
            >
              <span>Inspect 48-Seed Proof</span>
              <ArrowRight size={15} />
            </button>
          </div>

          {/* Unboxed Editorial Baseline Strip (Separated by vertical rules, no cards) */}
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(3, 1fr)',
              gap: '24px',
              paddingTop: '24px',
              borderTop: '1px solid var(--border-subtle)',
            }}
          >
            <div>
              <div
                className="mono tabular"
                style={{ fontSize: 'clamp(1.75rem, 2.6vw, 2.35rem)', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.1 }}
              >
                98%
              </div>
              <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '6px' }}>
                Seed reliability reaching ≥ 0.70 (vs 58% ALDE)
              </div>
            </div>

            <div>
              <div
                className="mono tabular"
                style={{ fontSize: 'clamp(1.75rem, 2.6vw, 2.35rem)', fontWeight: 600, color: '#0d9488', lineHeight: 1.1 }}
              >
                50 assays
              </div>
              <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '6px' }}>
                Median budget to hit top tier (2x faster)
              </div>
            </div>

            <div>
              <div
                className="mono tabular"
                style={{ fontSize: 'clamp(1.75rem, 2.6vw, 2.35rem)', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.1 }}
              >
                $0.0094
              </div>
              <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '6px' }}>
                Modal A10G GPU spend per 4-round campaign
              </div>
            </div>
          </div>
        </div>

        {/* Right: Automated Non-Interactive 3D Choreography Surface */}
        <div
          style={{
            background: 'var(--bg-elevated)',
            borderRadius: '18px',
            padding: '22px',
            display: 'flex',
            flexDirection: 'column',
            gap: '16px',
          }}
        >
          {/* Top Target Switcher Row */}
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '8px' }}>
            <div>
              <div style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)' }}>
                {activeProtein.name}
              </div>
              <div style={{ fontSize: '12.5px', color: 'var(--text-dim)' }}>
                {activeProtein.category}
              </div>
            </div>

            <div style={{ display: 'flex', gap: '4px', background: 'var(--bg-card)', padding: '3px', borderRadius: '8px' }}>
              {SHOWCASE_PROTEINS.map((prot, idx) => (
                <button
                  key={prot.id}
                  type="button"
                  onClick={() => {
                    setProteinIdx(idx);
                    setAnimPhase(0);
                  }}
                  className="mono"
                  style={{
                    padding: '5px 9px',
                    fontSize: '11.5px',
                    fontWeight: proteinIdx === idx ? 600 : 500,
                    borderRadius: '6px',
                    border: 'none',
                    cursor: 'pointer',
                    background: proteinIdx === idx ? 'var(--bg-elevated)' : 'transparent',
                    color: proteinIdx === idx ? 'var(--text-main)' : 'var(--text-muted)',
                  }}
                >
                  {prot.pdbId}
                </button>
              ))}
            </div>
          </div>

          {/* Non-Interactive Automated WebGL Viewport */}
          <div
            style={{
              position: 'relative',
              width: '100%',
              height: '340px',
              borderRadius: '12px',
              overflow: 'hidden',
              background: '#fcfbf8',
              pointerEvents: 'none',
            }}
          >
            <div ref={viewerContainerRef} style={{ width: '100%', height: '100%' }} />

            {/* Subtle Phase Pill Overlay */}
            <div
              style={{
                position: 'absolute',
                top: '14px',
                left: '14px',
                zIndex: 10,
                display: 'flex',
                gap: '6px',
              }}
            >
              {[
                { idx: 0, label: '1. Wild-Type Scaffold' },
                { idx: 1, label: '2. Epistatic Pocket' },
                { idx: 2, label: '3. Optimized Lead' },
              ].map((step) => {
                const active = animPhase === step.idx;
                return (
                  <span
                    key={step.idx}
                    style={{
                      padding: '4px 10px',
                      borderRadius: '6px',
                      fontSize: '11.5px',
                      fontWeight: active ? 600 : 500,
                      background: active ? '#141413' : 'rgba(239, 234, 225, 0.85)',
                      color: active ? '#fcfbf8' : 'var(--text-muted)',
                      transition: 'all 0.25s ease',
                    }}
                  >
                    {step.label}
                  </span>
                );
              })}
            </div>
          </div>

          {/* Bottom Genotype & Gain Readout */}
          <div
            style={{
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
              flexWrap: 'wrap',
              gap: '12px',
              padding: '14px 16px',
              borderRadius: '10px',
              background: 'var(--bg-card)',
            }}
          >
            <div>
              <div style={{ fontSize: '11.5px', color: 'var(--text-dim)', marginBottom: '2px' }}>
                {animPhase === 2 ? 'Discovered Lead Genotype' : 'Wild-Type Pocket Residues'}
              </div>
              <div className="mono" style={{ fontSize: '13.5px', fontWeight: 600, color: 'var(--text-main)' }}>
                {animPhase === 2 ? activeProtein.leadMutation : activeProtein.wtResidues}
              </div>
            </div>

            <div style={{ textAlign: 'right' }}>
              <div style={{ fontSize: '11.5px', color: 'var(--text-dim)', marginBottom: '2px' }}>
                Empirical Activity
              </div>
              <div className="mono tabular" style={{ fontSize: '13.5px', fontWeight: 600, color: '#0d9488' }}>
                {animPhase === 2 ? activeProtein.leadMetric : activeProtein.wtMetric}
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* =====================================================================
          SECTION 2: CORE CAPABILITIES (#capabilities) — Open-Canvas Continuous Rail
          ===================================================================== */}
      <section id="capabilities" style={{ scrollMarginTop: '96px' }}>
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'flex-end',
            flexWrap: 'wrap',
            gap: '20px',
            marginBottom: '36px',
          }}
        >
          <div>
            <h2
              className="title-display"
              style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '10px' }}
            >
              Six engines, one{' '}
              <span style={{ color: '#0d9488' }}>
                deterministic
              </span>{' '}
              loop.
            </h2>
            <p style={{ fontSize: 'var(--step-body)', color: 'var(--text-muted)', maxWidth: '62ch' }}>
              Every module is isolated by strict contracts: the LLM proposes boundaries, calibrated math selects variants, and the SQLite ledger gates every micro-dollar of GPU spend.
            </p>
          </div>
        </div>

        {/* Continuous Marquee Rail (Slows gently on hover, no outer box wrapper) */}
        <div
          className="marquee-container"
          onMouseEnter={() => handleMarqueeSpeed(0.32)}
          onMouseLeave={() => handleMarqueeSpeed(1.0)}
        >
          <div className="marquee-track" ref={marqueeTrackRef}>
            {[...CAPABILITY_CARDS, ...CAPABILITY_CARDS].map((card, i) => (
              <div
                key={`${card.title}-${i}`}
                style={{
                  width: '340px',
                  flexShrink: 0,
                  background: 'var(--bg-elevated)',
                  borderRadius: '14px',
                  padding: '24px',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-between',
                  gap: '18px',
                }}
              >
                <div
                  style={{
                    background: 'var(--bg-dark)',
                    borderRadius: '10px',
                    padding: '12px 14px',
                  }}
                >
                  {card.renderGraphic()}
                </div>

                <div>
                  <div
                    style={{
                      display: 'flex',
                      justifyContent: 'space-between',
                      alignItems: 'baseline',
                      gap: '8px',
                      marginBottom: '8px',
                    }}
                  >
                    <h3 style={{ fontSize: '17px', fontWeight: 600, color: 'var(--text-main)' }}>
                      {card.title}
                    </h3>
                  </div>
                  <p style={{ fontSize: '13.5px', color: 'var(--text-muted)', lineHeight: 1.55, marginBottom: '14px' }}>
                    {card.desc}
                  </p>
                  <div
                    className="mono tabular"
                    style={{
                      fontSize: '12px',
                      fontWeight: 600,
                      color: '#0d9488',
                      paddingTop: '12px',
                      borderTop: '1px solid var(--border-subtle)',
                    }}
                  >
                    {card.metric}
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* =====================================================================
          SECTION 3: EMPIRICAL COST BENCHMARKS (#benchmarks) — Asymmetric Editorial Split
          ===================================================================== */}
      <section
        id="benchmarks"
        style={{
          scrollMarginTop: '96px',
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))',
          gap: '48px',
          alignItems: 'start',
        }}
      >
        {/* Left: Editorial Savings Anchor */}
        <div>
          <h2
            className="title-display"
            style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '14px' }}
          >
            99.1% lower campaign cost,{' '}
            <span style={{ color: '#0d9488' }}>
              pinned first.
            </span>
          </h2>

          <p style={{ fontSize: 'var(--step-body)', color: 'var(--text-muted)', lineHeight: 1.6, marginBottom: '28px' }}>
            Wet-lab DNA synthesis ($45–$65 per variant) dominates protein engineering budgets. By resolving epistatic dead-ends in memory and ordering only one 96-well plate of Rank-0 leads, PanoptesBio eliminates $138,700 in wasted synthesis per campaign.
          </p>

          {/* Single Tonal Savings Callout Block */}
          <div
            style={{
              background: 'var(--bg-card)',
              borderRadius: '14px',
              padding: '28px',
            }}
          >
            <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginBottom: '6px' }}>
              Net Capital Saved per Target Campaign
            </div>
            <div
              className="mono tabular"
              style={{
                fontSize: 'clamp(2.5rem, 4.2vw, 3.5rem)',
                fontWeight: 600,
                color: '#0d9488',
                lineHeight: 1.05,
                marginBottom: '18px',
              }}
            >
              $138,700
            </div>

            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr',
                gap: '16px',
                paddingTop: '16px',
                borderTop: '1px solid var(--border-subtle)',
              }}
            >
              <div>
                <div className="mono tabular" style={{ fontSize: '18px', fontWeight: 600, color: 'var(--text-main)' }}>
                  50 – 96 wells
                </div>
                <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginTop: '2px' }}>
                  Total DNA variants synthesized (vs 3,000+ brute-force)
                </div>
              </div>
              <div>
                <div className="mono tabular" style={{ fontSize: '18px', fontWeight: 600, color: 'var(--text-main)' }}>
                  10 days
                </div>
                <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginTop: '2px' }}>
                  Target-to-validated lead turnaround (vs 20 weeks)
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Right: Clean Stacked Horizontal Comparison Ledger (PanoptesBio Pinned at Top) */}
        <div
          style={{
            background: 'var(--bg-elevated)',
            borderRadius: '16px',
            padding: '32px',
            display: 'flex',
            flexDirection: 'column',
            gap: '28px',
          }}
        >
          {[
            {
              name: 'PanoptesBio Active Learning + 96-Well Multiplexing',
              tag: 'Recommended · Pinned Benchmark',
              totalCost: '$1,290',
              breakdown: '$0.01 Modal A10G GPU + $1,290 single 96-well Twist synthesis plate',
              turnaround: '10 Days',
              barPct: 3.5,
              accent: '#0d9488',
              highlight: true,
            },
            {
              name: 'Standard Cloud Bio-Wrapper Platforms',
              tag: 'Ungated GPU Inference + 6-Plate Screen',
              totalCost: '$31,900',
              breakdown: '$3,100 unbudgeted cloud GPU compute + $28,800 across 6 synthesis plates',
              turnaround: '6 Weeks',
              barPct: 38,
              accent: '#78736a',
              highlight: false,
            },
            {
              name: 'Traditional Wet-Lab Directed Evolution',
              tag: 'Random Mutagenesis · 3,000+ Variants',
              totalCost: '$140,000',
              breakdown: '32+ 96-well plates, manual colony picking, and multi-month assay cycles',
              turnaround: '20 Weeks',
              barPct: 100,
              accent: '#a8a29e',
              highlight: false,
            },
          ].map((row, idx) => (
            <div
              key={row.name}
              style={{
                paddingBottom: idx < 2 ? '24px' : 0,
                borderBottom: idx < 2 ? '1px solid var(--border-subtle)' : 'none',
              }}
            >
              <div
                style={{
                  display: 'flex',
                  justifyContent: 'space-between',
                  alignItems: 'baseline',
                  gap: '12px',
                  marginBottom: '6px',
                }}
              >
                <div>
                  <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
                    <span
                      style={{
                        fontSize: '16px',
                        fontWeight: 600,
                        color: row.highlight ? '#0d9488' : 'var(--text-main)',
                      }}
                    >
                      {row.name}
                    </span>
                  </div>
                  <div style={{ fontSize: '12.5px', color: 'var(--text-dim)', marginTop: '2px' }}>
                    {row.breakdown}
                  </div>
                </div>

                <div style={{ textAlign: 'right', flexShrink: 0 }}>
                  <div
                    className="mono tabular"
                    style={{
                      fontSize: row.highlight ? '24px' : '20px',
                      fontWeight: 600,
                      color: row.highlight ? '#0d9488' : 'var(--text-main)',
                    }}
                  >
                    {row.totalCost}
                  </div>
                  <div className="mono" style={{ fontSize: '11.5px', color: 'var(--text-dim)' }}>
                    {row.turnaround}
                  </div>
                </div>
              </div>

              {/* Clean Tonal Bar */}
              <div
                style={{
                  width: '100%',
                  height: row.highlight ? '12px' : '8px',
                  borderRadius: '999px',
                  background: 'var(--bg-card)',
                  overflow: 'hidden',
                  marginTop: '10px',
                }}
              >
                <div
                  style={{
                    width: `${row.barPct}%`,
                    height: '100%',
                    borderRadius: '999px',
                    background: row.accent,
                  }}
                />
              </div>
            </div>
          ))}
        </div>
      </section>

      {/* =====================================================================
          SECTION 4: EMPIRICAL DATASET PROOF & LIVE SIGN-EPISTASIS ARTEFACT (#datasets)
          ===================================================================== */}
      <section
        id="datasets"
        style={{
          scrollMarginTop: '96px',
          background: 'var(--bg-card)',
          borderRadius: '20px',
          padding: 'clamp(28px, 4vw, 48px)',
        }}
      >
        {/* Header + Interactive Proof Mode Switcher */}
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'flex-end',
            flexWrap: 'wrap',
            gap: '20px',
            marginBottom: '32px',
          }}
        >
          <div>
            <h2
              className="title-display"
              style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '8px' }}
            >
              Verified on 48 paired seeds &{' '}
              <span style={{ color: '#0d9488' }}>
                sign-epistatic
              </span>{' '}
              landscapes.
            </h2>
            <p style={{ fontSize: 'var(--step-body)', color: 'var(--text-muted)', maxWidth: '64ch' }}>
              Tested head-to-head against published active-learning baselines on the complete 159,999-variant TrpB empirical fitness landscape (Johnston et al. 2024) and three additional benchmark families.
            </p>
          </div>

          {/* 3-Mode Interactive Switcher */}
          <div
            style={{
              display: 'inline-flex',
              gap: '4px',
              padding: '4px',
              borderRadius: '10px',
              background: 'var(--bg-elevated)',
            }}
          >
            {[
              { id: 'reliability', label: '48-Seed Reliability' },
              { id: 'epistasis', label: 'K227 × G228 Sign-Epistasis' },
              { id: 'datasets', label: '4-Target Benchmarks' },
            ].map((tab) => (
              <button
                key={tab.id}
                type="button"
                onClick={() => setProofTab(tab.id)}
                style={{
                  padding: '8px 14px',
                  fontSize: '13px',
                  fontWeight: proofTab === tab.id ? 600 : 500,
                  borderRadius: '7px',
                  border: 'none',
                  cursor: 'pointer',
                  background: proofTab === tab.id ? '#141413' : 'transparent',
                  color: proofTab === tab.id ? '#fcfbf8' : 'var(--text-muted)',
                  transition: 'all 0.15s ease',
                }}
              >
                {tab.label}
              </button>
            ))}
          </div>
        </div>

        {/* MODE 1: 48-Seed Reliability vs ALDE / Ridge / BoTorch */}
        {proofTab === 'reliability' && (
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))',
              gap: '28px',
              background: 'var(--bg-elevated)',
              borderRadius: '14px',
              padding: '32px',
            }}
          >
            <div style={{ display: 'flex', flexDirection: 'column', justifyContent: 'space-between', gap: '20px' }}>
              <div>
                <h3 style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '10px' }}>
                  Why reliability matters more than lucky seeds
                </h3>
                <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6, marginBottom: '18px' }}>
                  A wet-lab team runs a campaign once, not 50 times hoping for a lucky initialization. Across 48 paired random starts (96 total assays out of 159,999 variants), PanoptesBio reaches the ≥ 0.70 top tier on <strong>47 of 48 seeds (98%)</strong> and hits <strong>88% success at half budget (50 assays)</strong>.
                </p>
              </div>

              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1fr',
                  gap: '16px',
                  padding: '18px',
                  borderRadius: '10px',
                  background: 'var(--bg-card)',
                }}
              >
                <div>
                  <div style={{ fontSize: '12px', color: 'var(--text-dim)' }}>Reliability Mode (Default)</div>
                  <div className="mono tabular" style={{ fontSize: '22px', fontWeight: 600, color: '#0d9488', marginTop: '2px' }}>
                    98% Hit Rate
                  </div>
                  <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginTop: '2px' }}>
                    47/48 seeds ≥ 0.70 · 0.619 floor
                  </div>
                </div>
                <div>
                  <div style={{ fontSize: '12px', color: 'var(--text-dim)' }}>Frontier Mode (Peak Hunter)</div>
                  <div className="mono tabular" style={{ fontSize: '22px', fontWeight: 600, color: 'var(--text-main)', marginTop: '2px' }}>
                    38% Exact 1.000
                  </div>
                  <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginTop: '2px' }}>
                    Finds #1/159,999 global peak (92% ≥ 0.70)
                  </div>
                </div>
              </div>
            </div>

            {/* Right: Paired 48-Seed Bar Comparison */}
            <div style={{ display: 'flex', flexDirection: 'column', gap: '18px', justifyContent: 'center' }}>
              {[
                {
                  label: 'PanoptesBio (Reliability Mode)',
                  fullPct: 98,
                  halfPct: 88,
                  medianAssays: '50 assays',
                  color: '#0d9488',
                  bold: true,
                },
                {
                  label: 'PanoptesBio (Qualify-Then-Optimize Frontier)',
                  fullPct: 92,
                  halfPct: 60,
                  medianAssays: '49 assays (38% hit 1.000)',
                  color: '#0f766e',
                  bold: true,
                },
                {
                  label: 'One-Hot Ridge Surrogate Baseline',
                  fullPct: 67,
                  halfPct: 6,
                  medianAssays: '96 assays',
                  color: '#78736a',
                  bold: false,
                },
                {
                  label: 'ALDE (Published DNN Ensemble Baseline)',
                  fullPct: 58,
                  halfPct: 4,
                  medianAssays: '96 assays',
                  color: '#a8a29e',
                  bold: false,
                },
              ].map((item) => (
                <div key={item.label}>
                  <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '13px', marginBottom: '6px' }}>
                    <span style={{ fontWeight: item.bold ? 600 : 500, color: 'var(--text-main)' }}>
                      {item.label}
                    </span>
                    <span className="mono tabular" style={{ fontWeight: 600, color: item.bold ? '#0d9488' : 'var(--text-muted)' }}>
                      {item.fullPct}% full · {item.halfPct}% @ 50 assays
                    </span>
                  </div>
                  <div
                    style={{
                      width: '100%',
                      height: '10px',
                      borderRadius: '999px',
                      background: 'var(--bg-card)',
                      overflow: 'hidden',
                      position: 'relative',
                    }}
                  >
                    <div
                      style={{
                        width: `${item.fullPct}%`,
                        height: '100%',
                        background: item.color,
                        borderRadius: '999px',
                      }}
                    />
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* MODE 2: Interactive K227 × G228 Sign-Epistasis Matrix */}
        {proofTab === 'epistasis' && (
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))',
              gap: '28px',
              background: 'var(--bg-elevated)',
              borderRadius: '14px',
              padding: '32px',
              alignItems: 'center',
            }}
          >
            {/* Left: Interactive 2x2 Epistasis Grid */}
            <div>
              <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-muted)', marginBottom: '12px' }}>
                Click any residue pair at Positions 227 × 228 to inspect sign epistasis:
              </div>

              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '12px' }}>
                {EPISTASIS_CELLS.map((cell) => {
                  const isSelected = cell.id === activeEpistasis.id;
                  return (
                    <button
                      key={cell.id}
                      type="button"
                      onClick={() => setSelectedEpistasisId(cell.id)}
                      style={{
                        textAlign: 'left',
                        padding: '18px',
                        borderRadius: '12px',
                        background: cell.bgTone,
                        outline: isSelected ? `2px solid ${cell.tagColor}` : 'none',
                        border: 'none',
                        cursor: 'pointer',
                        transition: 'all 0.15s ease',
                      }}
                    >
                      <div className="mono" style={{ fontSize: '12px', fontWeight: 600, color: cell.tagColor, marginBottom: '4px' }}>
                        {cell.pos227} × {cell.pos228}
                      </div>
                      <div className="mono tabular" style={{ fontSize: '24px', fontWeight: 600, color: 'var(--text-main)' }}>
                        {cell.peakFitness}
                      </div>
                      <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginTop: '4px' }}>
                        {cell.role}
                      </div>
                    </button>
                  );
                })}
              </div>
            </div>

            {/* Right: Detailed Epistatic Explanation */}
            <div
              style={{
                background: 'var(--bg-card)',
                borderRadius: '12px',
                padding: '26px',
                display: 'flex',
                flexDirection: 'column',
                gap: '14px',
              }}
            >
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '8px' }}>
                <span className="mono" style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)' }}>
                  {activeEpistasis.genotype}
                </span>
                <span
                  className="mono tabular"
                  style={{
                    fontSize: '13px',
                    fontWeight: 600,
                    color: activeEpistasis.tagColor,
                  }}
                >
                  Peak: {activeEpistasis.peakFitness} (Mean: {activeEpistasis.meanFitness})
                </span>
              </div>

              <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
                {activeEpistasis.explanation}
              </p>

              <div
                style={{
                  paddingTop: '12px',
                  borderTop: '1px solid var(--border-subtle)',
                  fontSize: '13px',
                  color: 'var(--text-main)',
                }}
              >
                <strong>How PanoptesBio solves this:</strong> Once the surrogate qualifies a strong multi-site scaffold (predicted ≥ 0.45), 3 batch slots systematically probe under-sampled amino acids at unlocked positions, discovering <strong>K227 inside G228</strong> without wasting plates on random noise.
              </div>
            </div>
          </div>
        )}

        {/* MODE 3: Streamlined Single-Card Dataset Browser */}
        {proofTab === 'datasets' && (
          <div
            style={{
              background: 'var(--bg-elevated)',
              borderRadius: '14px',
              padding: '32px',
            }}
          >
            {/* Target Selector Pills + Prev/Next */}
            <div
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                flexWrap: 'wrap',
                gap: '12px',
                marginBottom: '24px',
                paddingBottom: '18px',
                borderBottom: '1px solid var(--border-subtle)',
              }}
            >
              <div style={{ display: 'flex', gap: '8px', flexWrap: 'wrap' }}>
                {DATASET_TARGETS.map((t, idx) => (
                  <button
                    key={t.id}
                    type="button"
                    onClick={() => setDatasetIdx(idx)}
                    style={{
                      padding: '7px 14px',
                      borderRadius: '8px',
                      border: 'none',
                      cursor: 'pointer',
                      fontSize: '13px',
                      fontWeight: datasetIdx === idx ? 600 : 500,
                      background: datasetIdx === idx ? '#141413' : 'var(--bg-card)',
                      color: datasetIdx === idx ? '#fcfbf8' : 'var(--text-muted)',
                    }}
                  >
                    {t.name}
                  </button>
                ))}
              </div>

              <div style={{ display: 'flex', gap: '6px' }}>
                <button
                  type="button"
                  onClick={() =>
                    setDatasetIdx((prev) => (prev - 1 + DATASET_TARGETS.length) % DATASET_TARGETS.length)
                  }
                  className="btn btn-secondary"
                  style={{ padding: '7px 10px' }}
                  aria-label="Previous dataset"
                >
                  <ChevronLeft size={15} />
                </button>
                <button
                  type="button"
                  onClick={() => setDatasetIdx((prev) => (prev + 1) % DATASET_TARGETS.length)}
                  className="btn btn-secondary"
                  style={{ padding: '7px 10px' }}
                  aria-label="Next dataset"
                >
                  <ChevronRight size={15} />
                </button>
              </div>
            </div>

            {/* Active Target Baseline -> Final Lead Summary */}
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))',
                gap: '28px',
                alignItems: 'center',
              }}
            >
              <div>
                <div className="mono" style={{ fontSize: '12px', color: '#0d9488', fontWeight: 600, marginBottom: '4px' }}>
                  {activeDataset.citation}
                </div>
                <h3 style={{ fontSize: '24px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '12px' }}>
                  {activeDataset.name}
                </h3>
                <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
                  {activeDataset.takeaway}
                </p>
              </div>

              <div
                style={{
                  background: 'var(--bg-card)',
                  borderRadius: '12px',
                  padding: '22px',
                  display: 'flex',
                  flexDirection: 'column',
                  gap: '14px',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline' }}>
                  <span style={{ fontSize: '13px', color: 'var(--text-muted)' }}>Starting Baseline</span>
                  <span className="mono tabular" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-main)' }}>
                    {activeDataset.wtScore} ({activeDataset.wtLabel})
                  </span>
                </div>

                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline' }}>
                  <span style={{ fontSize: '13px', color: 'var(--text-muted)' }}>Discovered Lead</span>
                  <span className="mono tabular" style={{ fontSize: '15px', fontWeight: 600, color: '#0d9488' }}>
                    {activeDataset.leadLabel} → {activeDataset.leadScore}
                  </span>
                </div>

                <div
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    paddingTop: '12px',
                    borderTop: '1px solid var(--border-subtle)',
                    fontSize: '12.5px',
                    color: 'var(--text-muted)',
                  }}
                >
                  <span>{activeDataset.hitRate}</span>
                  <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
                    GPU Spend: {activeDataset.spend}
                  </span>
                </div>
              </div>
            </div>
          </div>
        )}
      </section>

      {/* =====================================================================
          SECTION 5A: HOW IT WORKS (#how-it-works) — Single Unified Visual Flow
          ===================================================================== */}
      <section id="how-it-works" style={{ scrollMarginTop: '96px' }}>
        <div style={{ marginBottom: '36px' }}>
          <h2
            className="title-display"
            style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '10px' }}
          >
            How raw sequence data becomes a{' '}
            <span style={{ color: '#0d9488' }}>
              validated plate.
            </span>
          </h2>
          <p style={{ fontSize: 'var(--step-body)', color: 'var(--text-muted)', maxWidth: '64ch' }}>
            Click any stage below to follow how a wild-type PDB structure flows through our 3-stage deterministic pipeline into synthesis-ready cDNA orders.
          </p>
        </div>

        {/* Unified 3-Stage Selector Bar */}
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))',
            gap: '16px',
            marginBottom: '20px',
          }}
        >
          {pipelineStages.map((st, idx) => {
            const isSelected = activePipelineStage === idx;
            return (
              <button
                key={st.num}
                type="button"
                onClick={() => setActivePipelineStage(idx)}
                style={{
                  textAlign: 'left',
                  padding: '24px',
                  borderRadius: '14px',
                  border: 'none',
                  cursor: 'pointer',
                  background: isSelected ? '#141413' : 'var(--bg-elevated)',
                  color: isSelected ? '#fcfbf8' : 'var(--text-main)',
                  transition: 'all 0.18s ease',
                  display: 'flex',
                  flexDirection: 'column',
                  gap: '10px',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <span
                    className="mono"
                    style={{
                      fontSize: '12px',
                      fontWeight: 600,
                      color: isSelected ? '#5eead4' : '#0d9488',
                    }}
                  >
                    {st.num}
                  </span>
                  <ArrowRight size={16} style={{ opacity: isSelected ? 0.9 : 0.35 }} />
                </div>
                <div style={{ fontSize: '17px', fontWeight: 600 }}>{st.title}</div>
                <div style={{ fontSize: '13px', opacity: isSelected ? 0.8 : 0.68, lineHeight: 1.45 }}>
                  {st.subtitle}
                </div>
              </button>
            );
          })}
        </div>

        {/* Unified Visual Transformation Detail Panel (No nested 1px sub-boxes) */}
        <div
          style={{
            background: 'var(--bg-elevated)',
            borderRadius: '16px',
            padding: '32px',
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(240px, 1fr))',
            gap: '28px',
            alignItems: 'start',
          }}
        >
          <div>
            <div style={{ fontSize: '12.5px', fontWeight: 600, color: 'var(--text-dim)', marginBottom: '6px' }}>
              Incoming Data
            </div>
            <div className="mono" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.5 }}>
              {currentStage.inputLabel}
            </div>
          </div>

          <div>
            <div style={{ fontSize: '12.5px', fontWeight: 600, color: 'var(--text-dim)', marginBottom: '6px' }}>
              Biophysical Transformation
            </div>
            <div style={{ fontSize: '14px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
              {currentStage.operationLabel}
            </div>
          </div>

          <div>
            <div style={{ fontSize: '12.5px', fontWeight: 600, color: '#0d9488', marginBottom: '6px' }}>
              Stage Output
            </div>
            <div className="mono" style={{ fontSize: '14px', fontWeight: 600, color: '#0d9488', lineHeight: 1.5 }}>
              {currentStage.outputLabel}
            </div>
          </div>
        </div>
      </section>

      {/* =====================================================================
          SECTION 5B: HOW TO USE (#how-to-use) — Open 3-Column Operator Guide
          ===================================================================== */}
      <section id="how-to-use" style={{ scrollMarginTop: '96px' }}>
        <div style={{ marginBottom: '36px' }}>
          <h2
            className="title-display"
            style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '10px' }}
          >
            From accession ID to synthesis order in{' '}
            <span style={{ color: '#0d9488' }}>
              three steps.
            </span>
          </h2>
        </div>

        <div
          style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))',
            gap: '32px',
            paddingTop: '24px',
            borderTop: '1px solid var(--border-subtle)',
          }}
        >
          {[
            {
              step: 'First',
              title: 'Select your target & budget cap',
              body: 'Pick a benchmark or upload your PDB scaffold in the Studio, then set a hard micro-USD compute cap (e.g., $1.00) and batch size.',
              actionLabel: 'Open Active Learning Studio',
            },
            {
              step: 'Second',
              title: 'Run the loop & inspect Pareto leads',
              body: 'Watch real-time SSE telemetry as the ensemble surrogate ranks 160,000+ variants, relaxes steric clashes on CPU, and scores top leads.',
              actionLabel: 'Explore Pareto Frontier',
            },
            {
              step: 'Third',
              title: 'Download Twist CSV & robot scripts',
              body: 'Export codon-optimized Twist Golden Gate fragments, Opentrons OT-2 pipetting protocols, and 21 CFR Part 11 signed Merkle proofs.',
              actionLabel: 'Open Wet-Lab Export Studio',
            },
          ].map((item) => (
            <div key={item.step} style={{ display: 'flex', flexDirection: 'column', justifyContent: 'space-between', gap: '16px' }}>
              <div>
                <div
                  className="mono"
                  style={{ fontSize: '12px', fontWeight: 600, color: '#0d9488', marginBottom: '10px' }}
                >
                  {item.step}
                </div>
                <h3 style={{ fontSize: '19px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '8px' }}>
                  {item.title}
                </h3>
                <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
                  {item.body}
                </p>
              </div>

              <div>
                <button
                  type="button"
                  onClick={() => onNavigate && onNavigate('studio')}
                  style={{
                    background: 'none',
                    border: 'none',
                    padding: 0,
                    fontSize: '13.5px',
                    fontWeight: 600,
                    color: 'var(--text-main)',
                    cursor: 'pointer',
                    display: 'inline-flex',
                    alignItems: 'center',
                    gap: '6px',
                  }}
                >
                  <span>{item.actionLabel}</span>
                  <ArrowRight size={14} color="#0d9488" />
                </button>
              </div>
            </div>
          ))}
        </div>
      </section>

      {/* =====================================================================
          SECTION 6: FINAL CALL TO ACTION (#next-steps) — High-Contrast Dark Ink Panel
          ===================================================================== */}
      <section
        id="next-steps"
        style={{
          scrollMarginTop: '96px',
          background: '#141413',
          color: '#fcfbf8',
          borderRadius: '20px',
          padding: 'clamp(36px, 5vw, 64px)',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '32px',
        }}
      >
        <div style={{ maxWidth: '600px' }}>
          <h2
            className="title-display"
            style={{
              fontSize: 'clamp(2rem, 3.4vw, 3rem)',
              color: '#fcfbf8',
              marginBottom: '12px',
            }}
          >
            Ready to evolve your next{' '}
            <span style={{ color: '#5eead4' }}>
              therapeutic lead?
            </span>
          </h2>
          <p style={{ fontSize: '15.5px', color: '#b8b2a6', lineHeight: 1.6 }}>
            Launch an active-learning or de novo binder campaign in seconds, inspect the 4-objective Pareto frontier, and export synthesis-ready files under a cryptographic spend gate.
          </p>
        </div>

        <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap', alignItems: 'center' }}>
          <button
            type="button"
            onClick={() => onNavigate && onNavigate('studio')}
            className="btn"
            style={{
              background: '#0d9488',
              color: '#fcfbf8',
              padding: '14px 28px',
              fontSize: '15px',
              fontWeight: 600,
            }}
          >
            <Play size={16} />
            <span>Launch Studio</span>
          </button>

          <button
            type="button"
            onClick={() => onOpenAuth && onOpenAuth()}
            className="btn"
            style={{
              background: 'rgba(252, 251, 248, 0.1)',
              color: '#fcfbf8',
              padding: '14px 22px',
              fontSize: '14.5px',
            }}
          >
            <span>Create Operator Account</span>
            <ArrowRight size={15} />
          </button>
        </div>
      </section>
    </div>
  );
}
