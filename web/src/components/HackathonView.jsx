import React, { useState } from 'react';
import {
  Award,
  ShieldCheck,
  Zap,
  ExternalLink,
  GitBranch,
  ArrowRight,
  CheckCircle2,
  Database,
  Cpu,
  Lock,
  Activity,
  FileSpreadsheet,
} from 'lucide-react';

export default function HackathonView() {
  const [selectedNode, setSelectedNode] = useState(2);

  const pipelineNodes = [
    {
      title: 'UniProt & PDB Anchor',
      role: 'Structural Grounding',
      module: 'workers/rag_worker.py',
      summary: 'Extracts verified catalytic triads, cofactor coordinates (PLP/Chromophore), and Cα contact graphs before any model runs.',
      input: 'Target Accession (e.g., 8EH0, 1N8Z)',
      output: 'Protected Residue Mask + 8.0 Å Adjacency Graph',
      guarantee: 'Zero active-site obliteration',
      icon: Database,
    },
    {
      title: 'Meta Llama 3.1 Prior',
      role: 'Search Boundary Proposal',
      module: 'controller/prompts.py',
      summary: 'Operates strictly in deterministic JSON schema mode to propose mutable residue boundaries. Never executes math or calls paid cloud APIs.',
      input: 'Structural Context + Round Residuals',
      output: 'Strict JSON Combinatorial Search Space',
      guarantee: 'Fail-closed Pydantic schema gate',
      icon: Activity,
    },
    {
      title: '3D Kernel GP & CPU Relaxation',
      role: 'Biophysical Surrogate & Filter',
      module: 'workers/math_worker.py + relaxation_worker.py',
      summary: 'Evaluates 160,000+ epistatic combinations in memory via 5D Z-scales, k-DPP diversity, and <50ms CPU L-BFGS steric clash rescue.',
      input: 'Combinatorial Space (20⁴ Landscape)',
      output: 'Top-K Pareto Micro-Batch (8–24 Leads)',
      guarantee: '$0.000 GPU spend during screening',
      icon: Cpu,
    },
    {
      title: 'SQLite Micro-USD Gatekeeper',
      role: 'Deterministic FinOps Ledger',
      module: 'ledger/database.py + controller/gatekeeper.py',
      summary: 'Verifies remaining micro-USD allocation ($1 = 1,000,000 µUSD) with SQLite immutability triggers before unlocking serverless GPU containers.',
      input: 'Micro-Batch Cost Estimate ($0.0031)',
      output: 'Cryptographic Spend Authorization Token',
      guarantee: 'Zero unbudgeted cloud execution',
      icon: Lock,
    },
    {
      title: 'Modal A10G & CRO Export',
      role: 'Tensor Scoring & Wet-Lab Handoff',
      module: 'workers/esm_worker.py + wetlab_worker.py',
      summary: 'Executes ESM-2 650M pseudo-log-likelihood scoring on warm A10G GPUs and compiles Twist Golden Gate fragments + Opentrons OT-2 scripts.',
      input: 'Authorized Micro-Batch Sequences',
      output: 'Empirical Fitness + Twist CSV + Part 11 Merkle Hash',
      guarantee: 'HMAC-SHA256 auditable provenance',
      icon: FileSpreadsheet,
    },
  ];

  const targetBenchmarks = [
    {
      name: 'TrpB Synthase (Tm9D8*)',
      pdb: 'PDB: 8EH0',
      space: '159,999 combinatorial variants',
      wtScore: 1.0,
      leadScore: 7.16,
      maxScale: 8.0,
      enrichment: '13.0x',
      spend: '$0.0125',
      highlight: '100% of global empirical peak discovered in 4 rounds',
    },
    {
      name: 'Trastuzumab Fab (HER2)',
      pdb: 'PDB: 1N8Z',
      space: 'CDR-H3 paratope maturation',
      wtScore: 0.42,
      leadScore: 1.45,
      maxScale: 1.8,
      enrichment: '7.1x',
      spend: '$0.0064',
      highlight: 'Zero PTM deamidation or aggregation liabilities',
    },
    {
      name: 'Green Fluorescent Protein',
      pdb: 'PDB: 1EMA',
      space: 'Sarkisyan DMS landscape',
      wtScore: 0.35,
      leadScore: 1.03,
      maxScale: 1.3,
      enrichment: '8.4x',
      spend: '$0.0125',
      highlight: 'Superfolder epistatic triple-mutant recovered',
    },
    {
      name: 'Streptococcal GB1 Domain',
      pdb: 'PDB: 2QMT',
      space: 'IgG Fc-binding interface',
      wtScore: 0.01,
      leadScore: 0.044,
      maxScale: 0.06,
      enrichment: '9.2x',
      spend: '$0.0062',
      highlight: 'Rescued core packing clash via CPU L-BFGS relaxation',
    },
  ];

  const sprintVisuals = [
    {
      day: 'Hours 01 – 24',
      title: 'FinOps Ledger & Fail-Closed Gatekeeper',
      metric: '46 / 46 Ledger Checks',
      detail: 'Engineered integer micro-USD SQLite ledger with monotonic spend triggers and cold/warm Modal A10G cost equations.',
      svg: (
        <svg viewBox="0 0 220 72" style={{ width: '100%', height: '68px', display: 'block' }}>
          <rect x="8" y="16" width="58" height="40" rx="6" fill="#f2efe9" stroke="#141413" strokeWidth="1.5" />
          <text x="37" y="34" textAnchor="middle" fontSize="9" fontWeight="700" fill="#141413" fontFamily="monospace">PROPOSAL</text>
          <text x="37" y="46" textAnchor="middle" fontSize="8" fill="#575653" fontFamily="monospace">24 vars</text>
          <path d="M 66 36 L 86 36" stroke="#141413" strokeWidth="1.5" />
          <rect x="86" y="12" width="52" height="48" rx="6" fill="#0d9488" />
          <text x="112" y="33" textAnchor="middle" fontSize="8.5" fontWeight="700" fill="#fdfcfa" fontFamily="monospace">µUSD GATE</text>
          <text x="112" y="46" textAnchor="middle" fontSize="8" fill="#fdfcfa" fontFamily="monospace">$0.0031</text>
          <path d="M 138 36 L 158 36" stroke="#0d9488" strokeWidth="1.5" />
          <rect x="158" y="16" width="54" height="40" rx="6" fill="#fdfcfa" stroke="#0d9488" strokeWidth="1.5" />
          <text x="185" y="34" textAnchor="middle" fontSize="8.5" fontWeight="700" fill="#0d9488" fontFamily="monospace">APPROVED</text>
          <text x="185" y="46" textAnchor="middle" fontSize="7.5" fill="#575653" fontFamily="monospace">SHA-256</text>
        </svg>
      ),
    },
    {
      day: 'Hours 25 – 48',
      title: '3D Contact DKL Surrogate & Modal A10G Loop',
      metric: '13.0x Hit Enrichment',
      detail: 'Connected Meta Llama 3.1 boundary proposals to 3D Cα contact kernels, k-DPP batch diversity, and live Modal A10G ESM-2 650M scoring.',
      svg: (
        <svg viewBox="0 0 220 72" style={{ width: '100%', height: '68px', display: 'block' }}>
          <path d="M 12 58 Q 48 54, 74 36 T 136 22 T 206 14" fill="none" stroke="#0d9488" strokeWidth="2.2" />
          <path d="M 12 58 Q 48 54, 74 36 T 136 22 T 206 14 L 206 62 L 12 62 Z" fill="rgba(13, 148, 136, 0.12)" />
          <circle cx="45" cy="50" r="3.5" fill="#141413" />
          <circle cx="95" cy="29" r="3.5" fill="#141413" />
          <circle cx="152" cy="19" r="4" fill="#0d9488" />
          <circle cx="198" cy="15" r="5" fill="#0d9488" stroke="#fdfcfa" strokeWidth="1.5" />
          <text x="16" y="20" fontSize="8.5" fontWeight="700" fill="#141413" fontFamily="monospace">GP Posterior Climb</text>
        </svg>
      ),
    },
    {
      day: 'Hours 49 – 72',
      title: '4D Pareto, CPU L-BFGS & Wet-Lab CRO Studio',
      metric: '188 / 188 Tests Passing',
      detail: 'Added <50ms CPU sidechain clash relaxation, NSGA-II 4-objective Pareto sorting, and automated Twist Golden Gate + Opentrons OT-2 exports.',
      svg: (
        <svg viewBox="0 0 220 72" style={{ width: '100%', height: '68px', display: 'block' }}>
          <rect x="14" y="14" width="84" height="44" rx="6" fill="#f2efe9" stroke="#141413" strokeWidth="1.2" />
          <circle cx="34" cy="44" r="3" fill="#8c8a84" />
          <circle cx="52" cy="30" r="4" fill="#0d9488" />
          <circle cx="74" cy="24" r="4" fill="#0d9488" />
          <path d="M 28 48 Q 48 26, 82 22" fill="none" stroke="#0d9488" strokeWidth="1.5" strokeDasharray="3 2" />
          <path d="M 98 36 L 122 36" stroke="#141413" strokeWidth="1.5" />
          <rect x="122" y="14" width="84" height="44" rx="6" fill="#fdfcfa" stroke="#0d9488" strokeWidth="1.5" />
          <g fill="#0d9488">
            <circle cx="138" cy="26" r="2.5" /><circle cx="150" cy="26" r="2.5" /><circle cx="162" cy="26" r="2.5" /><circle cx="174" cy="26" r="2.5" /><circle cx="186" cy="26" r="2.5" />
            <circle cx="138" cy="36" r="2.5" /><circle cx="150" cy="36" r="2.5" /><circle cx="162" cy="36" r="2.5" /><circle cx="174" cy="36" r="2.5" /><circle cx="186" cy="36" r="2.5" />
            <circle cx="138" cy="46" r="2.5" /><circle cx="150" cy="46" r="2.5" /><circle cx="162" cy="46" r="2.5" /><circle cx="174" cy="46" r="2.5" /><circle cx="186" cy="46" r="2.5" />
          </g>
        </svg>
      ),
    },
  ];

  const activeNode = pipelineNodes[selectedNode];
  const ActiveIcon = activeNode.icon;

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '28px' }}>
      {/* Hero Banner + Live Summary Metrics */}
      <div className="glass-panel" style={{ padding: '36px', background: 'var(--bg-elevated)' }}>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))', gap: '32px', alignItems: 'center' }}>
          <div>
            <div style={{ display: 'inline-flex', alignItems: 'center', gap: '8px', padding: '5px 12px', borderRadius: '9999px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)', color: 'var(--text-main)', fontSize: '12px', fontWeight: 600, marginBottom: '16px' }}>
              <Award size={15} color="#0d9488" />
              <span>Meta Global AI Developer Hackathon</span>
            </div>

            <h1 style={{ fontSize: 'clamp(26px, 3vw, 34px)', fontWeight: 700, color: 'var(--text-main)', letterSpacing: '-0.025em', lineHeight: 1.15, marginBottom: '14px' }}>
              Deterministic Cloud-Economics for Protein Engineering
            </h1>

            <p style={{ color: 'var(--text-muted)', fontSize: '15px', lineHeight: 1.65, marginBottom: '24px' }}>
              Unconstrained AI agents waste thousands of dollars scoring non-functional proteins. PanoptesBio separates boundary reasoning (Meta Llama 3.1) from biophysical math (3D Kernel GP) and enforces a fail-closed SQLite micro-USD spend gate before touching Modal A10G GPUs.
            </p>

            <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap', alignItems: 'center' }}>
              <a
                href="https://github.com/fugeretrevin/VeritasBio"
                target="_blank"
                rel="noopener noreferrer"
                className="btn btn-primary"
                style={{ padding: '10px 18px', fontSize: '13px', textDecoration: 'none' }}
              >
                <GitBranch size={15} />
                <span>View GitHub Repository</span>
                <ExternalLink size={13} />
              </a>

              <div style={{ display: 'inline-flex', alignItems: 'center', gap: '6px', padding: '8px 14px', borderRadius: '8px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)', fontSize: '12px', fontWeight: 600, color: 'var(--emerald-primary)' }}>
                <ShieldCheck size={15} />
                <span>188 / 188 Automated Tests Passing</span>
              </div>
            </div>
          </div>

          {/* Right Visual Telemetry Matrix */}
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '12px' }}>
            <div style={{ padding: '20px', borderRadius: '10px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '6px' }}>4-Target Live GPU Spend</div>
              <div className="mono" style={{ fontSize: '28px', fontWeight: 700, color: '#0d9488', fontVariantNumeric: 'tabular-nums' }}>$0.0376</div>
              <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>160 ESM-2 650M passes on A10G</div>
            </div>

            <div style={{ padding: '20px', borderRadius: '10px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '6px' }}>TrpB Global Optimum</div>
              <div className="mono" style={{ fontSize: '28px', fontWeight: 700, color: 'var(--text-main)', fontVariantNumeric: 'tabular-nums' }}>7.1597</div>
              <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>100% peak hit in 159,999 library</div>
            </div>

            <div style={{ padding: '20px', borderRadius: '10px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '6px' }}>Discovery Enrichment</div>
              <div className="mono" style={{ fontSize: '28px', fontWeight: 700, color: 'var(--text-main)', fontVariantNumeric: 'tabular-nums' }}>13.0x</div>
              <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>Top-1% hit rate vs random screen</div>
            </div>

            <div style={{ padding: '20px', borderRadius: '10px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '6px' }}>Compute Cost Reduction</div>
              <div className="mono" style={{ fontSize: '28px', fontWeight: 700, color: '#0d9488', fontVariantNumeric: 'tabular-nums' }}>2,463x</div>
              <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>99.96% saved vs brute-force cloud</div>
            </div>
          </div>
        </div>
      </div>

      {/* Interactive Visual Architecture Pipeline */}
      <div className="glass-panel" style={{ padding: '36px', background: 'var(--bg-elevated)' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-end', flexWrap: 'wrap', gap: '12px', marginBottom: '24px' }}>
          <div>
            <h2 style={{ fontSize: '22px', fontWeight: 700, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
              Interactive System Blueprint
            </h2>
            <p style={{ color: 'var(--text-muted)', fontSize: '14px', margin: 0 }}>
              Click any node in the unidirectional 5-stage pipeline to inspect its exact input, output, and fail-closed guarantee.
            </p>
          </div>
          <span className="mono" style={{ fontSize: '12px', color: 'var(--text-muted)', padding: '6px 10px', borderRadius: '6px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)' }}>
            Strict Separation of Concerns (ARCHITECTURE.md §3)
          </span>
        </div>

        {/* Visual Node Diagram */}
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(175px, 1fr))', gap: '12px', marginBottom: '24px' }}>
          {pipelineNodes.map((node, idx) => {
            const IconComponent = node.icon;
            const isSelected = selectedNode === idx;
            return (
              <button
                key={node.title}
                type="button"
                onClick={() => setSelectedNode(idx)}
                style={{
                  textAlign: 'left',
                  padding: '18px',
                  borderRadius: '10px',
                  cursor: 'pointer',
                  background: isSelected ? 'var(--text-main)' : 'var(--bg-card)',
                  color: isSelected ? 'var(--bg-elevated)' : 'var(--text-main)',
                  border: isSelected ? '1px solid var(--text-main)' : '1px solid var(--border-subtle)',
                  transition: 'all 0.18s ease',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-between',
                  gap: '14px',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div style={{ width: '34px', height: '34px', borderRadius: '8px', display: 'flex', alignItems: 'center', justifyContent: 'center', background: isSelected ? 'rgba(255,255,255,0.14)' : 'var(--bg-elevated)', color: isSelected ? '#5eead4' : '#0d9488' }}>
                    <IconComponent size={18} />
                  </div>
                  <ArrowRight size={15} style={{ opacity: isSelected ? 0.9 : 0.35 }} />
                </div>
                <div>
                  <div style={{ fontSize: '11px', fontWeight: 600, opacity: isSelected ? 0.8 : 0.65, marginBottom: '4px' }}>
                    {node.role}
                  </div>
                  <div style={{ fontSize: '15px', fontWeight: 700, lineHeight: 1.25 }}>
                    {node.title}
                  </div>
                </div>
              </button>
            );
          })}
        </div>

        {/* Active Node Inspection Panel */}
        <div style={{ padding: '24px', borderRadius: '10px', background: 'var(--bg-card)', border: '1px solid var(--border-subtle)', display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))', gap: '24px', alignItems: 'center' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '10px' }}>
              <ActiveIcon size={18} color="#0d9488" />
              <span style={{ fontSize: '18px', fontWeight: 700, color: 'var(--text-main)' }}>{activeNode.title}</span>
              <span className="mono" style={{ fontSize: '11px', padding: '3px 8px', borderRadius: '4px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)', color: 'var(--text-muted)' }}>
                {activeNode.module}
              </span>
            </div>
            <p style={{ fontSize: '14px', color: 'var(--text-muted)', lineHeight: 1.6, margin: 0 }}>
              {activeNode.summary}
            </p>
          </div>

          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(160px, 1fr))', gap: '12px' }}>
            <div style={{ padding: '14px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>Stage Input</div>
              <div className="mono" style={{ fontSize: '12px', fontWeight: 600, color: 'var(--text-main)' }}>{activeNode.input}</div>
            </div>
            <div style={{ padding: '14px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>Stage Output</div>
              <div className="mono" style={{ fontSize: '12px', fontWeight: 600, color: '#0d9488' }}>{activeNode.output}</div>
            </div>
            <div style={{ padding: '14px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>Invariant Guarantee</div>
              <div style={{ fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', display: 'flex', alignItems: 'center', gap: '6px' }}>
                <CheckCircle2 size={14} color="#0d9488" />
                <span>{activeNode.guarantee}</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Visual 4-Target Live Production Benchmark Matrix */}
      <div className="glass-panel" style={{ padding: '36px', background: 'var(--bg-elevated)' }}>
        <div style={{ marginBottom: '24px' }}>
          <h2 style={{ fontSize: '22px', fontWeight: 700, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
            4-Target Live Production Benchmark Matrix
          </h2>
          <p style={{ color: 'var(--text-muted)', fontSize: '14px', margin: 0 }}>
            End-to-end validation across enzymes, therapeutic antibodies, and fluorescent reporters using live Meta Llama 3.1 + Modal A10G ESM-2 650M.
          </p>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))', gap: '16px' }}>
          {targetBenchmarks.map((item) => {
            const wtPct = Math.max(8, Math.round((item.wtScore / item.maxScale) * 100));
            const leadPct = Math.max(15, Math.round((item.leadScore / item.maxScale) * 100));
            return (
              <div
                key={item.name}
                style={{
                  padding: '22px',
                  borderRadius: '10px',
                  background: 'var(--bg-card)',
                  border: '1px solid var(--border-subtle)',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-between',
                  gap: '16px',
                }}
              >
                <div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '6px' }}>
                    <span className="mono" style={{ fontSize: '11px', fontWeight: 600, color: '#0d9488' }}>{item.pdb}</span>
                    <span className="mono" style={{ fontSize: '11px', color: 'var(--text-muted)' }}>Spend: {item.spend}</span>
                  </div>
                  <h3 style={{ fontSize: '16px', fontWeight: 700, color: 'var(--text-main)', marginBottom: '2px' }}>{item.name}</h3>
                  <div style={{ fontSize: '12px', color: 'var(--text-dim)' }}>{item.space}</div>
                </div>

                {/* Visual Baseline vs Lead Comparison Bars */}
                <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
                  <div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '4px' }}>
                      <span>Wild-Type Baseline</span>
                      <span className="mono" style={{ fontVariantNumeric: 'tabular-nums' }}>{item.wtScore}</span>
                    </div>
                    <div style={{ height: '8px', borderRadius: '999px', background: 'var(--bg-elevated)', overflow: 'hidden', border: '1px solid var(--border-subtle)' }}>
                      <div style={{ width: `${wtPct}%`, height: '100%', background: '#8c8a84', borderRadius: '999px' }} />
                    </div>
                  </div>

                  <div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '11px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '4px' }}>
                      <span>Discovered Lead Peak</span>
                      <span className="mono" style={{ color: '#0d9488', fontVariantNumeric: 'tabular-nums' }}>{item.leadScore} ({item.enrichment})</span>
                    </div>
                    <div style={{ height: '10px', borderRadius: '999px', background: 'var(--bg-elevated)', overflow: 'hidden', border: '1px solid var(--border-subtle)' }}>
                      <div style={{ width: `${leadPct}%`, height: '100%', background: '#0d9488', borderRadius: '999px' }} />
                    </div>
                  </div>
                </div>

                <div style={{ fontSize: '12px', color: 'var(--text-muted)', paddingTop: '10px', borderTop: '1px solid var(--border-subtle)' }}>
                  {item.highlight}
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Visual 72-Hour Engineering Sprint Progression */}
      <div className="glass-panel" style={{ padding: '36px', background: 'var(--bg-elevated)' }}>
        <div style={{ marginBottom: '24px' }}>
          <h2 style={{ fontSize: '22px', fontWeight: 700, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
            72-Hour Engineering Progression
          </h2>
          <p style={{ color: 'var(--text-muted)', fontSize: '14px', margin: 0 }}>
            How PanoptesBio evolved from an immutable micro-USD ledger to a closed-loop wet-lab CRO platform in three 24-hour phases.
          </p>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))', gap: '18px' }}>
          {sprintVisuals.map((phase) => (
            <div
              key={phase.day}
              style={{
                padding: '22px',
                borderRadius: '10px',
                background: 'var(--bg-card)',
                border: '1px solid var(--border-subtle)',
                display: 'flex',
                flexDirection: 'column',
                gap: '14px',
              }}
            >
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                <span className="mono" style={{ fontSize: '12px', fontWeight: 700, color: '#0d9488' }}>{phase.day}</span>
                <span className="mono" style={{ fontSize: '11px', padding: '3px 8px', borderRadius: '999px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)', color: 'var(--text-main)' }}>
                  {phase.metric}
                </span>
              </div>

              <div style={{ padding: '12px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
                {phase.svg}
              </div>

              <div>
                <h3 style={{ fontSize: '16px', fontWeight: 700, color: 'var(--text-main)', marginBottom: '6px' }}>
                  {phase.title}
                </h3>
                <p style={{ fontSize: '13px', color: 'var(--text-muted)', lineHeight: 1.55, margin: 0 }}>
                  {phase.detail}
                </p>
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
