import React, { useState, useEffect, useRef } from 'react';
import { ArrowRight, ChevronLeft, ChevronRight } from 'lucide-react';

export default function HomeView({ onNavigate }) {
  const [activeSection, setActiveSection] = useState('hero');
  const [hoveredDot, setHoveredDot] = useState(null);
  const [isRailHovered, setIsRailHovered] = useState(false);

  // Automated 3D Showcase State (0 = Base Wild-Type, 1 = Hotspot Mapping, 2 = Optimized Lead)
  const [pdbContent, setPdbContent] = useState('');
  const [autoPhase, setAutoPhase] = useState(0);
  const container3DRef = useRef(null);
  const viewer3DRef = useRef(null);

  // Continuous Marquee Ref for smooth slow-down on hover (never pauses)
  const marqueeTrackRef = useRef(null);

  // One-at-a-time Dataset Proof Stepper State
  const [datasetIdx, setDatasetIdx] = useState(0);

  // Interactive Benchmark Library Slider
  const [candidateCount, setCandidateCount] = useState(160000);

  // Interactive Visual Pipeline Node State
  const [activePipelineNode, setActivePipelineNode] = useState(1);

  const sections = [
    { id: 'hero', label: 'Our Mission' },
    { id: 'capabilities', label: 'Core Capabilities' },
    { id: 'benchmarks', label: 'Cost Benchmarks' },
    { id: 'datasets', label: 'Dataset Proof' },
    { id: 'how-it-works', label: 'How It Works' },
    { id: 'how-to-use', label: 'How To Use' },
    { id: 'next-steps', label: 'Launch Studio' },
  ];

  useEffect(() => {
    const handleScroll = () => {
      const scrollPos = window.scrollY + 260;
      for (const section of sections) {
        const el = document.getElementById(section.id);
        if (el) {
          const top = el.offsetTop;
          const height = el.offsetHeight;
          if (scrollPos >= top && scrollPos < top + height) {
            setActiveSection(section.id);
            break;
          }
        }
      }
    };
    window.addEventListener('scroll', handleScroll, { passive: true });
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  const scrollToSection = (id) => {
    const el = document.getElementById(id);
    if (el) {
      el.scrollIntoView({ behavior: 'smooth' });
    }
  };

  // Load crystallographic PDB coordinates for automated hero showcase
  useEffect(() => {
    let isCancelled = false;
    const fetchPdb = async () => {
      try {
        const res = await fetch('/api/structure/8EH0');
        if (res.ok) {
          const data = await res.json();
          if (!isCancelled && data.pdb_content) {
            setPdbContent(data.pdb_content);
          }
        }
      } catch (err) {
        console.warn('PDB fetch error:', err);
      }
    };
    fetchPdb();
    return () => {
      isCancelled = true;
    };
  }, []);

  // Automated 3-Phase Timer (cycles Base -> Hotspots -> Optimized Lead every 4.5s)
  useEffect(() => {
    const timer = setInterval(() => {
      setAutoPhase((prev) => (prev + 1) % 3);
    }, 4500);
    return () => clearInterval(timer);
  }, []);

  // Initialize 3Dmol Viewer Once PDB Loads (non-interactive, automated motion on warm off-white surface)
  useEffect(() => {
    if (!container3DRef.current || !pdbContent) return;

    let retryTimer = null;
    const initViewer = () => {
      if (typeof window.$3Dmol === 'undefined') {
        retryTimer = setTimeout(initViewer, 200);
        return;
      }

      try {
        container3DRef.current.innerHTML = '';
        const viewer = window.$3Dmol.createViewer(container3DRef.current, {
          backgroundColor: '#f2efe9',
          antialias: true,
        });
        viewer3DRef.current = viewer;

        viewer.addModel(pdbContent, 'pdb');
        viewer.setStyle({}, { cartoon: { color: '#334155', opacity: 0.92 } });
        viewer.zoomTo();
        viewer.render();
        viewer.spin('y', 0.55);
      } catch (err) {
        console.error('Error initializing automated 3Dmol showcase:', err);
      }
    };

    initViewer();
    return () => {
      if (retryTimer) clearTimeout(retryTimer);
      if (viewer3DRef.current) {
        try {
          viewer3DRef.current.spin(false);
        } catch (e) {
          // ignore cleanup errors
        }
      }
    };
  }, [pdbContent]);

  // Apply Smooth Automated Visual Transitions When autoPhase Changes
  useEffect(() => {
    const viewer = viewer3DRef.current;
    if (!viewer || typeof window.$3Dmol === 'undefined') return;

    try {
      viewer.removeAllSurfaces();
      const hotspots = [183, 184, 227, 228];

      if (autoPhase === 0) {
        viewer.setStyle({}, { cartoon: { color: '#334155', opacity: 0.92 } });
        viewer.addStyle({ hetflag: true }, { stick: { color: '#64748b', radius: 0.25 } });
        viewer.zoomTo({}, 900);
        viewer.render();
      } else if (autoPhase === 1) {
        viewer.setStyle({}, { cartoon: { color: '#94a3b8', opacity: 0.55 } });
        viewer.addStyle(
          { resi: hotspots },
          {
            stick: { color: '#0369a1', radius: 0.38 },
            sphere: { color: '#0369a1', radius: 0.42, opacity: 0.85 },
          }
        );
        viewer.addStyle(
          { hetflag: true },
          {
            stick: { color: '#141413', radius: 0.34 },
            sphere: { color: '#0369a1', radius: 0.48, opacity: 0.9 },
          }
        );
        viewer.zoomTo({ resi: [83, 183, 184, 227, 228] }, 1100);
        viewer.render();
      } else if (autoPhase === 2) {
        viewer.setStyle({}, { cartoon: { color: '#141413', opacity: 0.88 } });
        viewer.addStyle(
          { resi: hotspots },
          {
            stick: { color: '#059669', radius: 0.42 },
            sphere: { color: '#10b981', radius: 0.48, opacity: 0.92 },
          }
        );
        viewer.addStyle(
          { within: { distance: 4.5, sel: { resi: hotspots } } },
          { stick: { color: '#0369a1', radius: 0.22 } }
        );
        viewer.addSurface(
          window.$3Dmol.SurfaceType.VDW,
          { opacity: 0.28, color: '#10b981' },
          { resi: [83, 183, 184, 227, 228] }
        );
        viewer.zoomTo({}, 1000);
        viewer.render();
      }
    } catch (err) {
      console.warn('3D phase transition warning:', err);
    }
  }, [autoPhase]);

  const phaseDescriptions = [
    {
      title: 'Wild-Type Base Scaffold',
      detail: 'Native Tm9D8* Tryptophan Synthase backbone (1.00x baseline turnover)',
      accent: '#334155',
    },
    {
      title: 'Catalytic & Allosteric Mapping',
      detail: 'Locking the 4.2 Å PLP cofactor shell & targeting residues 183, 184, 227, 228',
      accent: '#0369a1',
    },
    {
      title: 'Optimized Lead (V183L:F184I:V227C:S228G)',
      detail: 'L-BFGS relaxed active-site pocket · 7.16x catalytic turnover (100% library peak)',
      accent: '#059669',
    },
  ];

  // Smooth slow-down on hover for the continuous capabilities rail
  const handleMarqueeSpeed = (rate) => {
    if (!marqueeTrackRef.current) return;
    try {
      const animations = marqueeTrackRef.current.getAnimations();
      animations.forEach((anim) => {
        if (typeof anim.updatePlaybackRate === 'function') {
          anim.updatePlaybackRate(rate);
        } else {
          anim.playbackRate = rate;
        }
      });
    } catch (e) {
      // Ignore if Web Animations API not supported
    }
  };

  const capabilities = [
    {
      title: 'Micro-USD Spend Gatekeeper',
      category: 'Deterministic FinOps',
      subtitle: 'Every cloud GPU pass requires atomic integer micro-USD pre-authorization before execution.',
      metric: '$0.000061',
      visual: 'ledger',
    },
    {
      title: '3D Deep Kernel Gaussian Process',
      category: 'Surrogate Modeling',
      subtitle: 'Couples 1280-dim ESM-2 representations with C-alpha contact graphs and heteroscedastic assay noise.',
      metric: 'ρ = 0.82',
      visual: 'kernel',
    },
    {
      title: 'k-DPP Epistatic Diversity',
      category: 'Batch Acquisition',
      subtitle: 'Determinantal point processes force orthogonal exploration across multi-mutant fitness peaks.',
      metric: '13.0x',
      visual: 'diversity',
    },
    {
      title: 'CPU Elastic Relaxation',
      category: 'Structural Biophysics',
      subtitle: 'Two-stage L-BFGS sidechain and backbone minimizer resolves steric clashes in under 50ms.',
      metric: '42 ms',
      visual: 'relax',
    },
    {
      title: '4D Pareto Trade-Off Engine',
      category: 'Clinical Developability',
      subtitle: 'Simultaneously optimizes affinity, thermal stability, synthesizability, and zero chemical liabilities.',
      metric: 'Rank 0',
      visual: 'pareto',
    },
    {
      title: 'Automated Wet-Lab Translation',
      category: 'Physical Execution',
      subtitle: 'Compiles codon-optimized Twist DNA orders, Opentrons OT-2 scripts, and 4PL dose-response curves.',
      metric: '96-Well',
      visual: 'wetlab',
    },
  ];

  const marqueeItems = [...capabilities, ...capabilities];

  const renderCapabilityVisual = (type) => {
    switch (type) {
      case 'ledger':
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            <rect x="0" y="10" width="120" height="8" rx="4" fill="#e4e0d7" />
            <rect x="0" y="10" width="28" height="8" rx="4" fill="#059669" />
            <circle cx="28" cy="14" r="6" fill="#141413" />
            <text x="40" y="32" fontSize="9" fill="#78746d" fontWeight="500">
              Pre-authorized ceiling locked
            </text>
          </svg>
        );
      case 'kernel':
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            <path d="M 4 28 Q 30 6, 58 20 T 116 8" fill="none" stroke="#0369a1" strokeWidth="2.2" />
            <path d="M 4 28 Q 30 2, 58 14 T 116 4 L 116 14 Q 88 26, 58 26 T 4 32 Z" fill="rgba(3, 105, 161, 0.12)" />
            <circle cx="58" cy="20" r="3.5" fill="#141413" />
            <circle cx="96" cy="11" r="3.5" fill="#059669" />
          </svg>
        );
      case 'diversity':
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            <circle cx="20" cy="18" r="4" fill="#141413" />
            <circle cx="52" cy="10" r="4" fill="#0369a1" />
            <circle cx="84" cy="24" r="4" fill="#059669" />
            <circle cx="106" cy="12" r="4" fill="#141413" />
            <line x1="20" y1="18" x2="52" y2="10" stroke="#a8a29e" strokeDasharray="2 2" />
            <line x1="52" y1="10" x2="84" y2="24" stroke="#a8a29e" strokeDasharray="2 2" />
            <line x1="84" y1="24" x2="106" y2="12" stroke="#a8a29e" strokeDasharray="2 2" />
          </svg>
        );
      case 'relax':
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            <path d="M 6 6 C 24 6, 28 28, 56 28 L 114 28" fill="none" stroke="#059669" strokeWidth="2.2" />
            <circle cx="6" cy="6" r="3.5" fill="#78746d" />
            <circle cx="56" cy="28" r="4" fill="#059669" />
            <text x="64" y="22" fontSize="9" fill="#047857" fontWeight="600">
              -18.4 kcal/mol
            </text>
          </svg>
        );
      case 'pareto':
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            <path d="M 14 30 Q 32 8, 106 8" fill="none" stroke="#141413" strokeWidth="2" strokeDasharray="3 2" />
            <circle cx="24" cy="20" r="3.5" fill="#059669" />
            <circle cx="54" cy="11" r="3.5" fill="#059669" />
            <circle cx="90" cy="8" r="3.5" fill="#059669" />
            <circle cx="65" cy="26" r="2.5" fill="#cfc9be" />
            <circle cx="42" cy="28" r="2.5" fill="#cfc9be" />
          </svg>
        );
      default:
        return (
          <svg viewBox="0 0 120 36" style={{ width: '100%', height: '36px' }}>
            {[0, 1, 2, 3, 4, 5, 6, 7].map((c) =>
              [0, 1].map((r) => (
                <circle
                  key={`${c}-${r}`}
                  cx={12 + c * 13}
                  cy={12 + r * 12}
                  r="3.5"
                  fill={(c + r) % 3 === 0 ? '#059669' : '#cfc9be'}
                />
              ))
            )}
          </svg>
        );
    }
  };

  // Empirical Datasets (Direct Baseline vs Discovered Peak — no intermediate round clutter)
  const datasetWalkthroughs = [
    {
      id: 'trpb',
      name: 'TrpB Tryptophan Synthase',
      tag: 'Industrial Biocatalysis · 159,999 Variants',
      whatItIs:
        'A complete 20⁴ combinatorial wet-lab dataset (Johnston et al. 2024) measuring four active-site entrance residues (183, 184, 227, 228) on Tryptophan Synthase.',
      whyItMatters:
        'Engineered TrpB synthesizes chiral pharmaceutical building blocks in water in a single step, replacing multi-step toxic heavy-metal organic reactions.',
      howWeDid:
        'Evaluated just 40 variants out of 159,999 possibilities and discovered the exact #1 global empirical optimum at $0.0125 total compute spend.',
      spend: '$0.0125',
      enrichment: '13.0x vs. Random',
      testedFraction: '40 / 159,999 (0.025%)',
      wtScore: '1.00x',
      wtBarPct: 14,
      peakScore: '7.16x (100% Global Optimum)',
      peakBarPct: 100,
      improvementBadge: '+616% Catalytic Activity',
      leadMutation: 'V183L : F184I : V227C : S228G',
    },
    {
      id: 'her2',
      name: 'Trastuzumab Antibody (HER2)',
      tag: 'Oncology Biologics · 36,391 Variants',
      whatItIs:
        'Deep mutational scanning dataset of the CDR-H3 binding loop of Trastuzumab (Herceptin) against the HER2 breast cancer receptor (Mason et al. 2021).',
      whyItMatters:
        'Monoclonal antibodies often gain binding affinity at the cost of aggregation or chemical instability. Optimizing affinity while keeping zero clinical liabilities is essential for drug formulation.',
      howWeDid:
        'Identified a double-mutant CDR-H3 lead with +45% binding improvement while maintaining a 1.00/1.00 clinical developability score and zero unpaired cysteines.',
      spend: '$0.0064',
      enrichment: '7.1x vs. Random',
      testedFraction: '40 / 36,391 (0.11%)',
      wtScore: '1.00x',
      wtBarPct: 55,
      peakScore: '1.45x Affinity + 1.00 CMC Score',
      peakBarPct: 100,
      improvementBadge: '+45% Binding · 0 Liabilities',
      leadMutation: 'A106G : M107L (CDR-H3)',
    },
    {
      id: 'gfp',
      name: 'avGFP Fluorescent Reporter',
      tag: 'Biosensor Folding · 51,714 Variants',
      whatItIs:
        'The gold-standard ProteinGym benchmark (Sarkisyan et al. 2016, Nature) measuring fluorescence intensity across 51,714 multi-mutant Aequorea victoria GFP sequences.',
      whyItMatters:
        'Fluorescence requires a tightly packed 11-stranded beta-barrel; a single steric clash quenches the internal chromophore completely.',
      howWeDid:
        'By combining 3D contact graphs with CPU sidechain relaxation, our engine avoided barrel-destabilizing mutations and recovered top-decile superfolder variants under $0.02.',
      spend: '$0.0125',
      enrichment: '8.4x vs. Random',
      testedFraction: '40 / 51,714 (0.077%)',
      wtScore: '1.000 Log Brightness',
      wtBarPct: 58,
      peakScore: '1.029 Superfolder Peak',
      peakBarPct: 100,
      improvementBadge: 'Top-1% Superfolder Brightness',
      leadMutation: 'T97A : N146D : M153T : N212D',
    },
    {
      id: 'gb1',
      name: 'Protein G (GB1 IgG Domain)',
      tag: 'Epistatic Landscape · 149,361 Variants',
      whatItIs:
        'An exhaustive 4-site combinatorial landscape (Wu et al. 2016) known for extreme sign epistasis where 98% of variants are non-functional.',
      whyItMatters:
        'Standard greedy optimization gets trapped in local dead valleys when individual mutations look deleterious on their own but are synergistic together.',
      howWeDid:
        'Using k-DPP determinantal diversity and Kriging-believer batch conditioning, PanoptesBio jumped across the fitness valley to capture the synergistic quadruple mutant.',
      spend: '$0.0064',
      enrichment: '9.2x vs. Random',
      testedFraction: '40 / 149,361 (0.027%)',
      wtScore: '0.000 (98% Dead Valley)',
      wtBarPct: 12,
      peakScore: '+0.0442 Epistatic Peak',
      peakBarPct: 100,
      improvementBadge: 'Escaped 98% Non-Functional Valley',
      leadMutation: 'V39F : D40A : G41A : V54A',
    },
  ];

  const currentDataset = datasetWalkthroughs[datasetIdx];

  // Visual System Pipeline Nodes ("How It Works")
  const pipelineNodes = [
    {
      id: 0,
      stage: 'Base Input',
      title: 'Sequence, 3D Structure & Budget Cap',
      shortDesc: 'Raw PDB coordinates + FASTA + hard USD ceiling',
      inputLabel: 'Operator Specification',
      inputValue: 'PDB 8EH0 (TrpB) + $1.00 Max Ceiling',
      processLabel: 'Spatial & Allosteric Pre-Filter',
      processValue: 'Locks 4.2 Å PLP cofactor shell & maps Kirchhoff normal-mode hinges',
      outputLabel: 'Feasible Search Space',
      outputValue: 'Residues [183, 184, 227, 228] unlocked; catalytic triad protected',
    },
    {
      id: 1,
      stage: 'Prior & Surrogate',
      title: 'Llama 3.1 + 3D Kernel Gaussian Process',
      shortDesc: 'Hypothesis generation & Bayesian uncertainty ranking',
      inputLabel: 'Feasible Search Space',
      inputValue: '4 coupled active-site entrance positions (160,000 combinations)',
      processLabel: 'GFlowNet + 3D Contact DKL + k-DPP',
      processValue: 'Projects 1280D embeddings + C-alpha graph; selects orthogonal peaks',
      outputLabel: 'Proposed Micro-Batch',
      outputValue: 'Top 10 diverse multi-mutants ranked by Expected Improvement',
    },
    {
      id: 2,
      stage: 'Spend Gate & Physics',
      title: 'Gatekeeper + CPU Relaxation + GPU Scoring',
      shortDesc: 'Micro-USD authorization, 42ms L-BFGS & ESM-2 scoring',
      inputLabel: 'Proposed Micro-Batch',
      inputValue: '10 candidate variants (Est. cost: $0.000061 / 61 µUSD)',
      processLabel: 'Fail-Closed Ledger + 2-Stage L-BFGS + Modal A10G',
      processValue: 'Pre-authorizes spend, relaxes steric clashes on CPU, scores on ESM-2 650M',
      outputLabel: 'Scored & Relaxed Leads',
      outputValue: 'V183L:F184I:V227C:S228G (-18.4 kcal/mol, 7.16x predicted turnover)',
    },
    {
      id: 3,
      stage: 'Physical Output',
      title: 'Pareto Filter & Turnkey Wet-Lab Package',
      shortDesc: 'Rank 0 leads exported to Twist DNA & Opentrons',
      inputLabel: 'Scored & Relaxed Leads',
      inputValue: '40 evaluated variants across 4 active-learning rounds',
      processLabel: 'NSGA-II 4D Pareto + Codon & Plate Compiler',
      processValue: 'Verifies 0 unpaired cysteines, optimizes codons, builds 36-well moat layout',
      outputLabel: 'Turnkey Lab Deliverables',
      outputValue: 'Twist 96-Well CSV + Opentrons OT-2 Script + Cryptographic Audit Proof',
    },
  ];

  const activeNodeData = pipelineNodes[activePipelineNode];

  const bruteforceCloudCost = Math.max(15, candidateCount * 0.02 * 0.000306 + 6.0);
  const wetlabCost = candidateCount * 0.35 + 5000;
  const panoptesCost = 40 * 0.02 * 0.000306 + 0.006242;
  const multipleCheaper = (bruteforceCloudCost / panoptesCost).toFixed(0);
  const wetlabSaved = Math.max(0, wetlabCost - 420);

  return (
    <div style={{ position: 'relative' }}>
      {/* Right-Side Dot Rail: Dots only at rest; labels pop out on hover so nothing overlaps */}
      <aside
        onMouseEnter={() => setIsRailHovered(true)}
        onMouseLeave={() => {
          setIsRailHovered(false);
          setHoveredDot(null);
        }}
        style={{
          position: 'fixed',
          right: '18px',
          top: '50%',
          transform: 'translateY(-50%)',
          zIndex: 900,
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'flex-end',
          gap: '10px',
          padding: '10px 6px',
        }}
        aria-label="Section navigation"
      >
        {sections.map((sec) => {
          const isActive = activeSection === sec.id;
          const showLabel = isRailHovered || hoveredDot === sec.id;
          return (
            <button
              key={sec.id}
              onClick={() => scrollToSection(sec.id)}
              onMouseEnter={() => setHoveredDot(sec.id)}
              onFocus={() => setHoveredDot(sec.id)}
              onBlur={() => setHoveredDot(null)}
              style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'flex-end',
                gap: '10px',
                background: 'none',
                border: 'none',
                padding: '2px 4px',
                cursor: 'pointer',
                position: 'relative',
              }}
              aria-label={sec.label}
            >
              {/* Pop-up Label on Hover */}
              <span
                style={{
                  fontSize: '12px',
                  fontWeight: isActive ? 600 : 500,
                  color: isActive ? 'var(--text-main)' : 'var(--text-muted)',
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  padding: '4px 10px',
                  borderRadius: '6px',
                  whiteSpace: 'nowrap',
                  opacity: showLabel ? 1 : 0,
                  transform: showLabel ? 'translateX(0)' : 'translateX(8px)',
                  pointerEvents: showLabel ? 'auto' : 'none',
                  transition: 'opacity 0.16s ease, transform 0.16s ease',
                }}
              >
                {sec.label}
              </span>

              {/* Minimal Dot Indicator */}
              <span
                style={{
                  width: '8px',
                  height: isActive ? '22px' : '8px',
                  borderRadius: '9999px',
                  background: isActive ? '#059669' : '#c7c1b6',
                  transition: 'all 0.2s ease',
                  flexShrink: 0,
                }}
              />
            </button>
          );
        })}
      </aside>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '120px', paddingBottom: '60px' }}>
        {/* =========================================================================
            SECTION 1: HERO - EDITORIAL STATEMENT & AUTOMATED 3D TRANSITION STAGE
            ========================================================================= */}
        <section id="hero" style={{ paddingTop: '16px' }}>
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(360px, 1fr))',
              gap: '56px',
              alignItems: 'center',
            }}
          >
            {/* Left Column: Mission & Clean Typographic Figures */}
            <div>
              <h1
                className="title-display"
                style={{
                  fontSize: 'var(--step-hero)',
                  fontWeight: 600,
                  lineHeight: 0.98,
                  letterSpacing: '-0.035em',
                  color: 'var(--text-main)',
                  marginBottom: '24px',
                }}
              >
                Protein discovery without the brute-force bill.
              </h1>

              <p
                style={{
                  fontSize: 'var(--step-lead)',
                  lineHeight: 1.5,
                  color: 'var(--text-muted)',
                  marginBottom: '36px',
                  maxWidth: '50ch',
                }}
              >
                Unguided generative pipelines burn thousands of dollars screening dead sequences.
                PanoptesBio pairs structural biophysics, 3D Deep Kernel active learning, and a
                fail-closed micro-USD ledger to reach empirical optima under $0.04.
              </p>

              {/* Primary Action */}
              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '20px',
                  flexWrap: 'wrap',
                  marginBottom: '48px',
                }}
              >
                <button
                  onClick={() => onNavigate('studio')}
                  className="btn btn-primary"
                  style={{
                    padding: '15px 30px',
                    fontSize: '16px',
                    fontWeight: 500,
                  }}
                >
                  <span>Launch Studio</span>
                  <ArrowRight size={17} />
                </button>

                <button
                  onClick={() => onNavigate('experiments')}
                  className="btn btn-secondary"
                  style={{
                    padding: '15px 24px',
                    fontSize: '15px',
                  }}
                >
                  <span>Explore Novel Leads</span>
                </button>
              </div>

              {/* 3 Focused Scientific Metrics */}
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
                    className="tabular"
                    style={{
                      fontSize: '28px',
                      fontWeight: 600,
                      letterSpacing: '-0.03em',
                      color: 'var(--text-main)',
                    }}
                  >
                    $0.0376
                  </div>
                  <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '2px' }}>
                    Total 4-target GPU spend
                  </div>
                </div>
                <div>
                  <div
                    className="tabular"
                    style={{
                      fontSize: '28px',
                      fontWeight: 600,
                      letterSpacing: '-0.03em',
                      color: '#059669',
                    }}
                  >
                    7.16x
                  </div>
                  <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '2px' }}>
                    TrpB global peak reached
                  </div>
                </div>
                <div>
                  <div
                    className="tabular"
                    style={{
                      fontSize: '28px',
                      fontWeight: 600,
                      letterSpacing: '-0.03em',
                      color: 'var(--text-main)',
                    }}
                  >
                    13.0x
                  </div>
                  <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginTop: '2px' }}>
                    Top-1% hit enrichment
                  </div>
                </div>
              </div>
            </div>

            {/* Right Column: Automated Non-Interactive 3D Structural Evolution Showcase */}
            <div
              style={{
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
              }}
            >
              <div
                style={{
                  width: '100%',
                  maxWidth: '520px',
                  background: 'var(--bg-card)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '20px',
                  padding: '20px',
                  position: 'relative',
                  overflow: 'hidden',
                }}
              >
                <div
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    marginBottom: '12px',
                    fontSize: '12px',
                  }}
                >
                  <span style={{ fontWeight: 600, color: 'var(--text-main)' }}>
                    PDB 8EH0 · Tryptophan Synthase Evolution
                  </span>
                  <span
                    style={{
                      padding: '3px 10px',
                      borderRadius: '9999px',
                      background: 'var(--bg-elevated)',
                      border: '1px solid var(--border-subtle)',
                      color: phaseDescriptions[autoPhase].accent,
                      fontWeight: 600,
                    }}
                  >
                    Stage {autoPhase + 1} of 3
                  </span>
                </div>

                {/* Non-interactive 3Dmol Canvas */}
                <div
                  style={{
                    width: '100%',
                    height: '360px',
                    borderRadius: '12px',
                    overflow: 'hidden',
                    position: 'relative',
                    pointerEvents: 'none',
                    userSelect: 'none',
                  }}
                >
                  <div ref={container3DRef} style={{ width: '100%', height: '100%' }} />
                </div>

                {/* Automated Stage Progression Footer */}
                <div
                  style={{
                    marginTop: '14px',
                    background: 'var(--bg-elevated)',
                    border: '1px solid var(--border-subtle)',
                    borderRadius: '10px',
                    padding: '14px 16px',
                  }}
                >
                  <div
                    style={{
                      fontSize: '14px',
                      fontWeight: 600,
                      color: 'var(--text-main)',
                      marginBottom: '4px',
                    }}
                  >
                    {phaseDescriptions[autoPhase].title}
                  </div>
                  <div style={{ fontSize: '13px', color: 'var(--text-muted)', marginBottom: '12px' }}>
                    {phaseDescriptions[autoPhase].detail}
                  </div>

                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '8px' }}>
                    {['Base Scaffold', 'Hotspot Shell', 'Optimized Lead'].map((lbl, idx) => {
                      const active = autoPhase === idx;
                      return (
                        <div key={lbl}>
                          <div
                            style={{
                              height: '4px',
                              borderRadius: '2px',
                              background: active ? '#059669' : 'var(--border-subtle)',
                              transition: 'background-color 0.3s ease',
                              marginBottom: '4px',
                            }}
                          />
                          <div
                            style={{
                              fontSize: '11px',
                              fontWeight: active ? 600 : 500,
                              color: active ? 'var(--text-main)' : 'var(--text-dim)',
                            }}
                          >
                            {lbl}
                          </div>
                        </div>
                      );
                    })}
                  </div>
                </div>
              </div>
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 2: CORE CAPABILITIES - CONTINUOUS RAIL (SLOWS ON HOVER)
            ========================================================================= */}
        <section
          id="capabilities"
          style={{
            background: 'var(--bg-card)',
            border: '1px solid var(--border-subtle)',
            padding: '48px 32px',
            borderRadius: '16px',
          }}
        >
          <div style={{ marginBottom: '32px' }}>
            <h2
              className="title-display"
              style={{
                fontSize: 'var(--step-head)',
                color: 'var(--text-main)',
                marginBottom: '6px',
              }}
            >
              Core engine capabilities.
            </h2>
            <p style={{ fontSize: '16px', color: 'var(--text-muted)' }}>
              Six coupled subsystems linking structural physics to wet-lab synthesis.
            </p>
          </div>

          <div
            className="marquee-container"
            onMouseEnter={() => handleMarqueeSpeed(0.32)}
            onMouseLeave={() => handleMarqueeSpeed(1.0)}
          >
            <div ref={marqueeTrackRef} className="marquee-track" style={{ gap: '20px' }}>
              {marqueeItems.map((cap, idx) => (
                <div
                  key={idx}
                  style={{
                    width: '320px',
                    flexShrink: 0,
                    background: 'var(--bg-elevated)',
                    border: '1px solid var(--border-subtle)',
                    borderRadius: '12px',
                    padding: '24px',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'space-between',
                    gap: '16px',
                  }}
                >
                  <div>
                    <div
                      style={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        marginBottom: '12px',
                      }}
                    >
                      <span
                        style={{
                          fontSize: '12px',
                          fontWeight: 600,
                          color: 'var(--cyan-bright)',
                        }}
                      >
                        {cap.category}
                      </span>
                      <span
                        className="tabular badge badge-emerald"
                        style={{
                          fontSize: '12px',
                          fontWeight: 600,
                        }}
                      >
                        {cap.metric}
                      </span>
                    </div>

                    <h3
                      style={{
                        fontSize: '18px',
                        fontWeight: 600,
                        color: 'var(--text-main)',
                        marginBottom: '8px',
                      }}
                    >
                      {cap.title}
                    </h3>
                    <p style={{ fontSize: '14px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                      {cap.subtitle}
                    </p>
                  </div>

                  <div
                    style={{
                      paddingTop: '12px',
                      borderTop: '1px solid var(--border-subtle)',
                    }}
                  >
                    {renderCapabilityVisual(cap.visual)}
                  </div>
                </div>
              ))}
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 3: COST BENCHMARKS - PANOPTESBIO AT THE TOP & SAVINGS EMPHASIZED
            ========================================================================= */}
        <section id="benchmarks">
          <div
            style={{
              background: 'var(--bg-card)',
              border: '1px solid var(--border-subtle)',
              borderRadius: '16px',
              padding: '48px 40px',
            }}
          >
            {/* Top Digestible Header + Giant Savings Callout */}
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
                gap: '32px',
                alignItems: 'end',
                marginBottom: '36px',
              }}
            >
              <div>
                <h2
                  className="title-display"
                  style={{
                    fontSize: 'var(--step-head)',
                    color: 'var(--text-main)',
                    marginBottom: '12px',
                  }}
                >
                  Save 99.96% on discovery compute.
                </h2>
                <p style={{ fontSize: '16px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                  By testing 40 mathematically selected variants instead of brute-forcing the entire
                  library, PanoptesBio eliminates thousands in wasted GPU and wet-lab spend.
                </p>
              </div>

              {/* Emphasized Savings Highlight Box */}
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid rgba(5, 150, 105, 0.35)',
                  borderRadius: '12px',
                  padding: '20px 24px',
                  display: 'flex',
                  justifyContent: 'space-between',
                  alignItems: 'center',
                  gap: '16px',
                }}
              >
                <div>
                  <div style={{ fontSize: '12px', fontWeight: 600, color: '#047857', marginBottom: '4px' }}>
                    NET CAPITAL SAVED PER CAMPAIGN
                  </div>
                  <div
                    className="tabular"
                    style={{ fontSize: '32px', fontWeight: 600, color: '#059669', letterSpacing: '-0.03em' }}
                  >
                    ${wetlabSaved.toLocaleString(undefined, { maximumFractionDigits: 0 })}+
                  </div>
                  <div style={{ fontSize: '13px', color: 'var(--text-muted)' }}>
                    {multipleCheaper}x cheaper compute than unguided cloud runs
                  </div>
                </div>
              </div>
            </div>

            {/* Stacked Bars: PanoptesBio First at the Top! */}
            <div style={{ display: 'flex', flexDirection: 'column', gap: '18px', marginBottom: '32px' }}>
              {/* Bar 1 (TOP): PanoptesBio Active Learning Engine */}
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '2px solid #059669',
                  borderRadius: '12px',
                  padding: '20px 24px',
                }}
              >
                <div
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'baseline',
                    flexWrap: 'wrap',
                    gap: '12px',
                    marginBottom: '10px',
                  }}
                >
                  <div style={{ display: 'flex', alignItems: 'center', gap: '10px', flexWrap: 'wrap' }}>
                    <span style={{ fontSize: '17px', fontWeight: 600, color: 'var(--text-main)' }}>
                      PanoptesBio Active Learning Engine
                    </span>
                    <span className="badge badge-emerald">
                      40 Guided Evaluations · {multipleCheaper}x Lower Cost
                    </span>
                  </div>
                  <span className="tabular" style={{ fontSize: '26px', fontWeight: 600, color: '#059669' }}>
                    ${panoptesCost.toFixed(4)}
                  </span>
                </div>
                <div
                  style={{
                    width: '100%',
                    height: '14px',
                    background: 'var(--bg-card)',
                    borderRadius: '7px',
                    overflow: 'hidden',
                  }}
                >
                  <div
                    style={{
                      width: '2.5%',
                      minWidth: '18px',
                      height: '100%',
                      background: '#059669',
                      borderRadius: '7px',
                    }}
                  />
                </div>
              </div>

              {/* Bar 2: Standard Cloud GPU Competitors */}
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '12px',
                  padding: '18px 24px',
                }}
              >
                <div
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'baseline',
                    marginBottom: '8px',
                  }}
                >
                  <div>
                    <span style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)' }}>
                      Standard Cloud GPU Competitors
                    </span>
                    <span style={{ fontSize: '13px', color: 'var(--text-dim)', marginLeft: '10px' }}>
                      Exhaustive forward passes across all {candidateCount.toLocaleString()} variants
                    </span>
                  </div>
                  <span className="tabular" style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-muted)' }}>
                    ${bruteforceCloudCost.toFixed(2)}
                  </span>
                </div>
                <div
                  style={{
                    width: '100%',
                    height: '12px',
                    background: 'var(--bg-card)',
                    borderRadius: '6px',
                    overflow: 'hidden',
                  }}
                >
                  <div
                    style={{
                      width: `${Math.min(82, Math.max(14, (candidateCount / 250000) * 75))}%`,
                      height: '100%',
                      background: '#64748b',
                      borderRadius: '6px',
                      transition: 'width 0.2s ease',
                    }}
                  />
                </div>
              </div>

              {/* Bar 3: Physical Wet-Lab Screening */}
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '12px',
                  padding: '18px 24px',
                }}
              >
                <div
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'baseline',
                    marginBottom: '8px',
                  }}
                >
                  <div>
                    <span style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)' }}>
                      Conventional Wet-Lab Library Synthesis & Assay
                    </span>
                    <span style={{ fontSize: '13px', color: 'var(--text-dim)', marginLeft: '10px' }}>
                      Unfiltered clonal gene synthesis, cloning, and plate consumables
                    </span>
                  </div>
                  <span className="tabular" style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-dim)' }}>
                    ${wetlabCost.toLocaleString(undefined, { maximumFractionDigits: 0 })}
                  </span>
                </div>
                <div
                  style={{
                    width: '100%',
                    height: '12px',
                    background: 'var(--bg-card)',
                    borderRadius: '6px',
                    overflow: 'hidden',
                  }}
                >
                  <div
                    style={{
                      width: `${Math.min(100, Math.max(28, (candidateCount / 250000) * 100))}%`,
                      height: '100%',
                      background: '#a8a29e',
                      borderRadius: '6px',
                      transition: 'width 0.2s ease',
                    }}
                  />
                </div>
              </div>
            </div>

            {/* Interactive Library Size Slider */}
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                gap: '24px',
                flexWrap: 'wrap',
                paddingTop: '12px',
              }}
            >
              <span style={{ fontSize: '14px', color: 'var(--text-muted)' }}>
                Test with different combinatorial library sizes:
              </span>
              <div style={{ flex: 1, minWidth: '240px', display: 'flex', alignItems: 'center', gap: '16px' }}>
                <input
                  type="range"
                  min="5000"
                  max="250000"
                  step="5000"
                  value={candidateCount}
                  onChange={(e) => setCandidateCount(Number(e.target.value))}
                  aria-label="Combinatorial library size"
                  style={{ flex: 1, accentColor: '#059669', cursor: 'pointer' }}
                />
                <span className="tabular" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-main)', minWidth: '130px', textAlign: 'right' }}>
                  {candidateCount.toLocaleString()} variants
                </span>
              </div>
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 4: DATASET PROOF - DIRECT BASELINE VS OPTIMUM (NO ROUND CLUTTER)
            ========================================================================= */}
        <section id="datasets">
          <div
            style={{
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'flex-end',
              flexWrap: 'wrap',
              gap: '20px',
              marginBottom: '24px',
            }}
          >
            <div>
              <h2
                className="title-display"
                style={{
                  fontSize: 'var(--step-head)',
                  color: 'var(--text-main)',
                  marginBottom: '10px',
                }}
              >
                Validated on real laboratory datasets.
              </h2>
              <p style={{ fontSize: '16px', color: 'var(--text-muted)', maxWidth: '58ch' }}>
                Step through each published wet-lab benchmark to see what the target protein does
                and the final lead discovered by PanoptesBio.
              </p>
            </div>

            {/* Previous / Next Controls */}
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
              <span style={{ fontSize: '13px', fontWeight: 500, color: 'var(--text-muted)', marginRight: '4px' }}>
                {datasetIdx + 1} of {datasetWalkthroughs.length}
              </span>
              <button
                onClick={() =>
                  setDatasetIdx((prev) => (prev - 1 + datasetWalkthroughs.length) % datasetWalkthroughs.length)
                }
                className="btn btn-secondary"
                style={{ padding: '8px 14px' }}
                aria-label="Previous dataset"
              >
                <ChevronLeft size={16} />
                <span>Previous</span>
              </button>
              <button
                onClick={() => setDatasetIdx((prev) => (prev + 1) % datasetWalkthroughs.length)}
                className="btn btn-primary"
                style={{ padding: '8px 16px' }}
                aria-label="Next dataset"
              >
                <span>Next Dataset</span>
                <ChevronRight size={16} />
              </button>
            </div>
          </div>

          {/* Dataset Selector Pills */}
          <div style={{ display: 'flex', gap: '8px', flexWrap: 'wrap', marginBottom: '20px' }}>
            {datasetWalkthroughs.map((ds, idx) => {
              const active = idx === datasetIdx;
              return (
                <button
                  key={ds.id}
                  onClick={() => setDatasetIdx(idx)}
                  className={`btn ${active ? 'btn-primary' : 'btn-secondary'}`}
                  style={{ padding: '8px 16px', fontSize: '13px' }}
                >
                  {ds.name}
                </button>
              );
            })}
          </div>

          {/* Single Focused Dataset Card — Direct Baseline vs Discovered Optimum */}
          <div
            style={{
              background: 'var(--bg-card)',
              border: '1px solid var(--border-subtle)',
              borderRadius: '16px',
              padding: '40px',
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))',
              gap: '44px',
              alignItems: 'center',
            }}
          >
            {/* Left: Plain-English Explanation */}
            <div>
              <div style={{ marginBottom: '12px' }}>
                <span className="badge badge-cyan">{currentDataset.tag}</span>
              </div>

              <h3
                className="title-display"
                style={{
                  fontSize: '30px',
                  color: 'var(--text-main)',
                  marginBottom: '20px',
                }}
              >
                {currentDataset.name}
              </h3>

              <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
                <div>
                  <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '3px' }}>
                    What This Dataset Is
                  </div>
                  <p style={{ fontSize: '15px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                    {currentDataset.whatItIs}
                  </p>
                </div>

                <div>
                  <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '3px' }}>
                    Why It Matters
                  </div>
                  <p style={{ fontSize: '15px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                    {currentDataset.whyItMatters}
                  </p>
                </div>

                <div>
                  <div style={{ fontSize: '13px', fontWeight: 600, color: '#059669', marginBottom: '3px' }}>
                    Result
                  </div>
                  <p style={{ fontSize: '15px', color: 'var(--text-main)', lineHeight: 1.55 }}>
                    {currentDataset.howWeDid}
                  </p>
                </div>
              </div>
            </div>

            {/* Right: Clean Direct Outcome Comparison (Wild-Type Baseline vs. Discovered Lead) */}
            <div
              style={{
                background: 'var(--bg-elevated)',
                border: '1px solid var(--border-subtle)',
                borderRadius: '14px',
                padding: '32px',
                display: 'flex',
                flexDirection: 'column',
                gap: '24px',
              }}
            >
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '8px' }}>
                <span style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)' }}>
                  Empirical Benchmark Outcome
                </span>
                <span className="badge badge-emerald">{currentDataset.improvementBadge}</span>
              </div>

              {/* Direct 2-Bar Comparison: Wild-Type vs Discovered Lead */}
              <div style={{ display: 'flex', flexDirection: 'column', gap: '18px' }}>
                <div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '13px', marginBottom: '6px' }}>
                    <span style={{ color: 'var(--text-muted)' }}>Wild-Type Parent Baseline</span>
                    <span className="tabular" style={{ fontWeight: 600, color: 'var(--text-muted)' }}>
                      {currentDataset.wtScore}
                    </span>
                  </div>
                  <div style={{ width: '100%', height: '12px', background: 'var(--bg-card)', borderRadius: '6px', overflow: 'hidden' }}>
                    <div
                      style={{
                        width: `${currentDataset.wtBarPct}%`,
                        height: '100%',
                        background: '#a8a29e',
                        borderRadius: '6px',
                      }}
                    />
                  </div>
                </div>

                <div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '14px', marginBottom: '6px' }}>
                    <span style={{ fontWeight: 600, color: 'var(--text-main)' }}>PanoptesBio Discovered Lead</span>
                    <span className="tabular" style={{ fontWeight: 600, color: '#059669' }}>
                      {currentDataset.peakScore}
                    </span>
                  </div>
                  <div style={{ width: '100%', height: '14px', background: 'var(--bg-card)', borderRadius: '7px', overflow: 'hidden' }}>
                    <div
                      style={{
                        width: `${currentDataset.peakBarPct}%`,
                        height: '100%',
                        background: '#059669',
                        borderRadius: '7px',
                      }}
                    />
                  </div>
                </div>
              </div>

              {/* 3 Key Telemetry Metrics */}
              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: 'repeat(3, 1fr)',
                  gap: '12px',
                  padding: '14px',
                  background: 'var(--bg-card)',
                  borderRadius: '10px',
                }}
              >
                <div>
                  <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Evaluated</div>
                  <div className="tabular" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-main)' }}>
                    {currentDataset.testedFraction}
                  </div>
                </div>
                <div>
                  <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Compute Spend</div>
                  <div className="tabular" style={{ fontSize: '14px', fontWeight: 600, color: '#059669' }}>
                    {currentDataset.spend}
                  </div>
                </div>
                <div>
                  <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Hit Enrichment</div>
                  <div className="tabular" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--cyan-bright)' }}>
                    {currentDataset.enrichment}
                  </div>
                </div>
              </div>

              <div
                style={{
                  paddingTop: '14px',
                  borderTop: '1px solid var(--border-subtle)',
                  display: 'flex',
                  justifyContent: 'space-between',
                  alignItems: 'center',
                  fontSize: '13px',
                }}
              >
                <span style={{ color: 'var(--text-muted)' }}>Discovered Lead Genotype</span>
                <span
                  className="mono"
                  style={{
                    fontWeight: 600,
                    color: 'var(--text-main)',
                    background: 'var(--bg-card)',
                    padding: '4px 10px',
                    borderRadius: '6px',
                    border: '1px solid var(--border-subtle)',
                  }}
                >
                  {currentDataset.leadMutation}
                </span>
              </div>
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 5A: HOW IT WORKS - VISUAL DATA TRANSFORMATION PIPELINE
            ========================================================================= */}
        <section
          id="how-it-works"
          style={{
            background: 'var(--bg-card)',
            border: '1px solid var(--border-subtle)',
            borderRadius: '16px',
            padding: '48px 40px',
          }}
        >
          <div style={{ maxWidth: '680px', marginBottom: '36px' }}>
            <h2
              className="title-display"
              style={{
                fontSize: 'var(--step-head)',
                color: 'var(--text-main)',
                marginBottom: '10px',
              }}
            >
              How data flows through the engine.
            </h2>
            <p style={{ fontSize: '16px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
              Follow a protein from raw crystal coordinates and a $1.00 budget ceiling through each
              mathematical filter to a ready-to-synthesize 96-well plate. Click any stage to inspect
              the live data transformation.
            </p>
          </div>

          {/* Visual Horizontal Pipeline Schematic */}
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))',
              gap: '14px',
              marginBottom: '28px',
            }}
          >
            {pipelineNodes.map((node, idx) => {
              const isSelected = activePipelineNode === idx;
              return (
                <button
                  key={node.id}
                  onClick={() => setActivePipelineNode(idx)}
                  style={{
                    textAlign: 'left',
                    padding: '20px',
                    borderRadius: '12px',
                    background: isSelected ? '#141413' : 'var(--bg-elevated)',
                    color: isSelected ? '#fdfcfa' : 'var(--text-main)',
                    border: isSelected ? '1px solid #141413' : '1px solid var(--border-subtle)',
                    cursor: 'pointer',
                    position: 'relative',
                    transition: 'all 0.18s ease',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'space-between',
                    gap: '12px',
                  }}
                >
                  <div>
                    <div
                      style={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        marginBottom: '8px',
                        fontSize: '12px',
                        fontWeight: 600,
                        color: isSelected ? '#a7f3d0' : '#059669',
                      }}
                    >
                      <span>{node.stage}</span>
                      <span>→</span>
                    </div>
                    <div style={{ fontSize: '16px', fontWeight: 600, marginBottom: '6px', lineHeight: 1.3 }}>
                      {node.title}
                    </div>
                  </div>
                  <div
                    style={{
                      fontSize: '13px',
                      color: isSelected ? 'rgba(253, 252, 250, 0.78)' : 'var(--text-muted)',
                      lineHeight: 1.45,
                    }}
                  >
                    {node.shortDesc}
                  </div>
                </button>
              );
            })}
          </div>

          {/* Visual Input -> Transformation -> Output Inspector */}
          <div
            style={{
              background: 'var(--bg-elevated)',
              border: '1px solid var(--border-subtle)',
              borderRadius: '14px',
              padding: '28px 32px',
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(240px, 1fr))',
              gap: '24px',
              alignItems: 'stretch',
            }}
          >
            {/* Column 1: Stage Input */}
            <div
              style={{
                padding: '18px 20px',
                borderRadius: '10px',
                background: 'var(--bg-card)',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
              }}
            >
              <div>
                <div style={{ fontSize: '12px', fontWeight: 600, color: 'var(--text-dim)', marginBottom: '6px' }}>
                  STAGE INPUT · {activeNodeData.inputLabel}
                </div>
                <div style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.45 }}>
                  {activeNodeData.inputValue}
                </div>
              </div>
              <div style={{ marginTop: '14px', fontSize: '12px', color: 'var(--text-muted)' }}>
                Feeds into subsystem →
              </div>
            </div>

            {/* Column 2: System Transformation */}
            <div
              style={{
                padding: '18px 20px',
                borderRadius: '10px',
                background: 'rgba(3, 105, 161, 0.06)',
                border: '1px solid rgba(3, 105, 161, 0.2)',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
              }}
            >
              <div>
                <div style={{ fontSize: '12px', fontWeight: 600, color: '#0369a1', marginBottom: '6px' }}>
                  ENGINE OPERATION · {activeNodeData.processLabel}
                </div>
                <div style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.45 }}>
                  {activeNodeData.processValue}
                </div>
              </div>
              <div style={{ marginTop: '14px', fontSize: '12px', color: '#0369a1', fontWeight: 500 }}>
                Deterministic verification pass →
              </div>
            </div>

            {/* Column 3: Stage Output */}
            <div
              style={{
                padding: '18px 20px',
                borderRadius: '10px',
                background: 'rgba(5, 150, 105, 0.08)',
                border: '1px solid rgba(5, 150, 105, 0.25)',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
              }}
            >
              <div>
                <div style={{ fontSize: '12px', fontWeight: 600, color: '#047857', marginBottom: '6px' }}>
                  VERIFIED OUTPUT · {activeNodeData.outputLabel}
                </div>
                <div style={{ fontSize: '15px', fontWeight: 600, color: 'var(--text-main)', lineHeight: 1.45 }}>
                  {activeNodeData.outputValue}
                </div>
              </div>
              <div style={{ marginTop: '14px', fontSize: '12px', color: '#047857', fontWeight: 600 }}>
                ✓ Committed to SQLite Merkle Ledger
              </div>
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 5B: HOW TO USE - BROKEN DOWN AS SIMPLY AS POSSIBLE
            ========================================================================= */}
        <section id="how-to-use">
          <div style={{ maxWidth: '640px', marginBottom: '36px' }}>
            <h2
              className="title-display"
              style={{
                fontSize: 'var(--step-head)',
                color: 'var(--text-main)',
                marginBottom: '10px',
              }}
            >
              How to use PanoptesBio in three steps.
            </h2>
            <p style={{ fontSize: '16px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
              No command-line setup or manual primer math required. Everything runs directly from the
              browser Studio.
            </p>
          </div>

          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
              gap: '24px',
            }}
          >
            {/* Step 1 */}
            <div
              style={{
                background: 'var(--bg-card)',
                border: '1px solid var(--border-subtle)',
                borderRadius: '14px',
                padding: '32px 28px',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                gap: '20px',
              }}
            >
              <div>
                <div
                  style={{
                    display: 'inline-block',
                    padding: '4px 10px',
                    borderRadius: '6px',
                    background: 'var(--bg-elevated)',
                    border: '1px solid var(--border-subtle)',
                    fontSize: '12px',
                    fontWeight: 600,
                    color: 'var(--text-main)',
                    marginBottom: '16px',
                  }}
                >
                  Step 1
                </div>
                <h3 style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '10px' }}>
                  Pick a protein & lock your budget
                </h3>
                <p style={{ fontSize: '15px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                  Select a target scaffold (enzyme, antibody, or antivenom) and set your dollar cap
                  (for example, $1.00). The fail-closed gatekeeper guarantees the engine cannot
                  overspend.
                </p>
              </div>
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '8px',
                  padding: '12px 14px',
                  fontSize: '13px',
                  display: 'flex',
                  justifyContent: 'space-between',
                }}
              >
                <span style={{ color: 'var(--text-muted)' }}>Default Ceiling</span>
                <strong className="tabular" style={{ color: 'var(--text-main)' }}>$1.00 USD Max</strong>
              </div>
            </div>

            {/* Step 2 */}
            <div
              style={{
                background: 'var(--bg-card)',
                border: '1px solid var(--border-subtle)',
                borderRadius: '14px',
                padding: '32px 28px',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                gap: '20px',
              }}
            >
              <div>
                <div
                  style={{
                    display: 'inline-block',
                    padding: '4px 10px',
                    borderRadius: '6px',
                    background: 'var(--bg-elevated)',
                    border: '1px solid var(--border-subtle)',
                    fontSize: '12px',
                    fontWeight: 600,
                    color: 'var(--text-main)',
                    marginBottom: '16px',
                  }}
                >
                  Step 2
                </div>
                <h3 style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '10px' }}>
                  Stream the active learning loop
                </h3>
                <p style={{ fontSize: '15px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                  Click Launch. Watch in real time as the engine proposes mutations, relaxes sidechain
                  clashes on CPU, scores top candidates, and plots the Pareto trade-off frontier.
                </p>
              </div>
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '8px',
                  padding: '12px 14px',
                  fontSize: '13px',
                  display: 'flex',
                  justifyContent: 'space-between',
                }}
              >
                <span style={{ color: 'var(--text-muted)' }}>Live Telemetry</span>
                <strong style={{ color: '#059669' }}>3D Structure + Pareto Rank 0</strong>
              </div>
            </div>

            {/* Step 3 */}
            <div
              style={{
                background: 'var(--bg-card)',
                border: '1px solid var(--border-subtle)',
                borderRadius: '14px',
                padding: '32px 28px',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                gap: '20px',
              }}
            >
              <div>
                <div
                  style={{
                    display: 'inline-block',
                    padding: '4px 10px',
                    borderRadius: '6px',
                    background: 'var(--bg-elevated)',
                    border: '1px solid var(--border-subtle)',
                    fontSize: '12px',
                    fontWeight: 600,
                    color: 'var(--text-main)',
                    marginBottom: '16px',
                  }}
                >
                  Step 3
                </div>
                <h3 style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '10px' }}>
                  Export lab orders or upload plate data
                </h3>
                <p style={{ fontSize: '15px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                  Download ready-to-order Twist DNA spreadsheets and Opentrons robot scripts—or drop
                  your wet-lab plate reader CSV into the Studio to calibrate the next round.
                </p>
              </div>
              <div
                style={{
                  background: 'var(--bg-elevated)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: '8px',
                  padding: '12px 14px',
                  fontSize: '13px',
                  display: 'flex',
                  justifyContent: 'space-between',
                }}
              >
                <span style={{ color: 'var(--text-muted)' }}>One-Click Exports</span>
                <strong style={{ color: 'var(--cyan-bright)' }}>Twist CSV · Opentrons .py</strong>
              </div>
            </div>
          </div>
        </section>

        {/* =========================================================================
            SECTION 6: SLEEK CENTERED CALL TO ACTION - LAUNCH STUDIO FRONT AND CENTER
            ========================================================================= */}
        <section id="next-steps">
          <div
            style={{
              background: 'var(--bg-card)',
              border: '1px solid var(--border-subtle)',
              borderRadius: '20px',
              padding: '72px 32px',
              textAlign: 'center',
              display: 'flex',
              flexDirection: 'column',
              alignItems: 'center',
            }}
          >
            <h2
              className="title-display"
              style={{
                fontSize: 'clamp(2.2rem, 3.8vw, 3.5rem)',
                color: 'var(--text-main)',
                marginBottom: '16px',
                maxWidth: '20ch',
              }}
            >
              Ready to run a budgeted discovery campaign?
            </h2>
            <p
              style={{
                fontSize: '17px',
                color: 'var(--text-muted)',
                maxWidth: '52ch',
                marginBottom: '36px',
                lineHeight: 1.55,
              }}
            >
              Launch the interactive Studio to stream active learning rounds, optimize multi-objective
              Pareto frontiers, and export turnkey DNA synthesis manifests.
            </p>

            {/* Primary Front-and-Center CTA */}
            <button
              onClick={() => onNavigate('studio')}
              className="btn btn-primary"
              style={{
                padding: '18px 44px',
                fontSize: '17px',
                fontWeight: 600,
                borderRadius: '10px',
                marginBottom: '20px',
              }}
            >
              <span>Launch Studio</span>
              <ArrowRight size={19} />
            </button>

            <div style={{ display: 'flex', gap: '24px', fontSize: '14px', color: 'var(--text-muted)' }}>
              <button
                onClick={() => onNavigate('experiments')}
                style={{
                  background: 'none',
                  border: 'none',
                  color: 'var(--text-muted)',
                  cursor: 'pointer',
                  fontWeight: 500,
                }}
              >
                Browse Novel Experiments →
              </button>
              <button
                onClick={() => onNavigate('hackathon')}
                style={{
                  background: 'none',
                  border: 'none',
                  color: 'var(--text-muted)',
                  cursor: 'pointer',
                  fontWeight: 500,
                }}
              >
                View Hackathon Architecture →
              </button>
            </div>
          </div>
        </section>
      </div>
    </div>
  );
}
