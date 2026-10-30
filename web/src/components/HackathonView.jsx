import React, { useState } from 'react';
import {
  ShieldCheck,
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
      summary:
        'Extracts verified catalytic triads, cofactor coordinates (PLP/Chromophore), and Cα contact graphs before any model runs.',
      input: 'Target Accession (e.g., 8EH0, 1N8Z)',
      output: 'Protected Residue Mask + 8.0 Å Adjacency Graph',
      guarantee: 'Zero active-site obliteration',
      icon: Database,
    },
    {
      title: 'Meta Llama 3.1 Prior',
      role: 'Boundary Proposal',
      module: 'controller/prompts.py',
      summary:
        'Operates strictly in deterministic JSON schema mode to propose mutable residue boundaries. Never executes math or calls paid cloud APIs.',
      input: 'Structural Context + Round Residuals',
      output: 'Strict JSON Combinatorial Search Space',
      guarantee: 'Fail-closed Pydantic schema gate',
      icon: Activity,
    },
    {
      title: 'Calibrated Ensemble & L-BFGS',
      role: 'Epistasis Surrogate',
      module: 'workers/math_worker.py + relaxation_worker.py',
      summary:
        'Evaluates 160,000+ epistatic combinations via pair-interaction features, qualify-then-optimize frontier slots, and <50ms CPU L-BFGS clash rescue.',
      input: 'Combinatorial Space (20⁴ Landscape)',
      output: 'Top-K Pareto Micro-Batch (24 Leads)',
      guarantee: '$0.000 GPU spend during screening',
      icon: Cpu,
    },
    {
      title: 'SQLite Micro-USD Gatekeeper',
      role: 'Deterministic FinOps',
      module: 'ledger/database.py + controller/gatekeeper.py',
      summary:
        'Verifies remaining micro-USD allocation ($1 = 1,000,000 µUSD) with SQLite immutability triggers before unlocking serverless GPU containers.',
      input: 'Micro-Batch Cost Estimate ($0.0031)',
      output: 'Cryptographic Spend Authorization Token',
      guarantee: 'Zero unbudgeted cloud execution',
      icon: Lock,
    },
    {
      title: 'Modal A10G & CRO Export',
      role: 'Tensor & Wet-Lab',
      module: 'workers/esm_worker.py + wetlab_worker.py',
      summary:
        'Executes ESM-2 650M pseudo-log-likelihood scoring on warm A10G GPUs and compiles Twist Golden Gate fragments + Opentrons OT-2 scripts.',
      input: 'Authorized Micro-Batch Sequences',
      output: 'Empirical Fitness + Twist CSV + Part 11 Merkle Hash',
      guarantee: 'HMAC-SHA256 auditable provenance',
      icon: FileSpreadsheet,
    },
  ];

  const targetBenchmarks = [
    {
      name: 'TrpB Synthase (Tm9D8*)',
      pdb: 'PDB: 8EH0 · 48 Paired Seeds',
      space: '159,999 combinatorial variants',
      wtScore: 0.402,
      leadScore: 1.0,
      maxScale: 1.0,
      enrichment: '98% Hit Rate',
      spend: '$0.0094',
      highlight: '47/48 seeds reach ≥ 0.70 (vs 58% ALDE); discovers exact global peak V183A:F184I:V227K:S228G',
    },
    {
      name: 'Trastuzumab Fab (HER2)',
      pdb: 'PDB: 1N8Z',
      space: 'CDR-H3 paratope maturation',
      wtScore: 0.42,
      leadScore: 1.45,
      maxScale: 1.6,
      enrichment: '7.1x',
      spend: '$0.0064',
      highlight: 'Simultaneously optimizes HER2 binding and removes D101 isomerization liability',
    },
    {
      name: 'Green Fluorescent Protein',
      pdb: 'PDB: 1EMA',
      space: 'Sarkisyan 54K DMS landscape',
      wtScore: 0.35,
      leadScore: 1.03,
      maxScale: 1.2,
      enrichment: '8.4x',
      spend: '$0.0125',
      highlight: 'Superfolder epistatic triple-mutant recovered with ρ = 0.85 out-of-sample correlation',
    },
    {
      name: 'Streptococcal GB1 Domain',
      pdb: 'PDB: 2QMT',
      space: '149,361 IgG Fc interface variants',
      wtScore: 0.18,
      leadScore: 0.92,
      maxScale: 1.0,
      enrichment: '9.2x',
      spend: '$0.0062',
      highlight: 'Rescued aromatic core packing clash via <50ms CPU L-BFGS sidechain relaxation',
    },
  ];

  const sprintVisuals = [
    {
      day: 'Stage I',
      title: 'FinOps Ledger & Fail-Closed Gatekeeper',
      metric: '46 / 46 Ledger Checks',
      detail:
        'Engineered integer micro-USD SQLite ledger with monotonic spend triggers, cold/warm Modal A10G cost equations, and HMAC-SHA256 Merkle audit chains.',
      svg: (
        <svg viewBox="0 0 220 68" style={{ width: '100%', height: '64px', display: 'block' }}>
          <rect x="8" y="14" width="58" height="40" rx="6" fill="#efeae1" />
          <text x="37" y="32" textAnchor="middle" fontSize="9" fontWeight="700" fill="#141413" fontFamily="monospace">PROPOSAL</text>
          <text x="37" y="44" textAnchor="middle" fontSize="8" fill="#52504a" fontFamily="monospace">24 vars</text>
          <path d="M 66 34 L 86 34" stroke="#141413" strokeWidth="1.5" />
          <rect x="86" y="10" width="52" height="48" rx="6" fill="#0d9488" />
          <text x="112" y="31" textAnchor="middle" fontSize="8.5" fontWeight="700" fill="#fcfbf8" fontFamily="monospace">µUSD GATE</text>
          <text x="112" y="44" textAnchor="middle" fontSize="8" fill="#fcfbf8" fontFamily="monospace">$0.0031</text>
          <path d="M 138 34 L 158 34" stroke="#0d9488" strokeWidth="1.5" />
          <rect x="158" y="14" width="54" height="40" rx="6" fill="#141413" />
          <text x="185" y="32" textAnchor="middle" fontSize="8.5" fontWeight="700" fill="#5eead4" fontFamily="monospace">APPROVED</text>
          <text x="185" y="44" textAnchor="middle" fontSize="7.5" fill="#fcfbf8" fontFamily="monospace">SHA-256</text>
        </svg>
      ),
    },
    {
      day: 'Stage II',
      title: 'Calibrated Pair-Interaction Ensemble & Frontier Rule',
      metric: '98% Seed Hit Rate',
      detail:
        'Combined amino-acid physicochemical descriptors with AA×AA pair interactions and qualify-then-optimize frontier slots to solve K227 × G228 sign epistasis.',
      svg: (
        <svg viewBox="0 0 220 68" style={{ width: '100%', height: '64px', display: 'block' }}>
          <path d="M 12 54 Q 48 50, 74 32 T 136 18 T 206 10" fill="none" stroke="#0d9488" strokeWidth="2.2" />
          <path d="M 12 54 Q 48 50, 74 32 T 136 18 T 206 10 L 206 58 L 12 58 Z" fill="rgba(13, 148, 136, 0.12)" />
          <circle cx="45" cy="46" r="3.5" fill="#141413" />
          <circle cx="95" cy="25" r="3.5" fill="#141413" />
          <circle cx="152" cy="15" r="4" fill="#0d9488" />
          <circle cx="198" cy="11" r="5" fill="#0d9488" />
          <text x="16" y="18" fontSize="8.5" fontWeight="700" fill="#141413" fontFamily="monospace">48-Seed Climb</text>
        </svg>
      ),
    },
    {
      day: 'Stage III',
      title: '4D Pareto, CPU L-BFGS & Wet-Lab CRO Studio',
      metric: '192+ Tests Passing',
      detail:
        'Added <50ms CPU sidechain clash relaxation, NSGA-II 4-objective Pareto sorting, and automated Twist Golden Gate + Opentrons OT-2 exports.',
      svg: (
        <svg viewBox="0 0 220 68" style={{ width: '100%', height: '64px', display: 'block' }}>
          <rect x="14" y="12" width="84" height="44" rx="6" fill="#efeae1" />
          <circle cx="34" cy="42" r="3" fill="#8c8a84" />
          <circle cx="52" cy="28" r="4" fill="#0d9488" />
          <circle cx="74" cy="22" r="4" fill="#0d9488" />
          <path d="M 28 46 Q 48 24, 82 20" fill="none" stroke="#0d9488" strokeWidth="1.5" strokeDasharray="3 2" />
          <path d="M 98 34 L 122 34" stroke="#141413" strokeWidth="1.5" />
          <rect x="122" y="12" width="84" height="44" rx="6" fill="rgba(13, 148, 136, 0.12)" />
          <g fill="#0d9488">
            <circle cx="138" cy="24" r="2.5" /><circle cx="150" cy="24" r="2.5" /><circle cx="162" cy="24" r="2.5" /><circle cx="174" cy="24" r="2.5" /><circle cx="186" cy="24" r="2.5" />
            <circle cx="138" cy="34" r="2.5" /><circle cx="150" cy="34" r="2.5" /><circle cx="162" cy="34" r="2.5" /><circle cx="174" cy="34" r="2.5" /><circle cx="186" cy="34" r="2.5" />
            <circle cx="138" cy="44" r="2.5" /><circle cx="150" cy="44" r="2.5" /><circle cx="162" cy="44" r="2.5" /><circle cx="174" cy="44" r="2.5" /><circle cx="186" cy="44" r="2.5" />
          </g>
        </svg>
      ),
    },
  ];

  const activeNode = pipelineNodes[selectedNode];
  const ActiveIcon = activeNode.icon;

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '56px' }}>
      {/* Open Editorial Hero + Telemetry Matrix (No pill badge above H1) */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))',
          gap: '40px',
          alignItems: 'center',
          paddingBottom: '36px',
          borderBottom: '1px solid var(--border-subtle)',
        }}
      >
        <div>
          <h1
            className="title-display"
            style={{
              fontSize: 'var(--step-head)',
              color: 'var(--text-main)',
              marginBottom: '14px',
            }}
          >
            Deterministic cloud-economics for{' '}
            <span style={{ color: '#0d9488' }}>
              protein engineering.
            </span>
          </h1>

          <p style={{ color: 'var(--text-muted)', fontSize: '15.5px', lineHeight: 1.6, marginBottom: '24px' }}>
            Built for the Meta Global AI Developer Hackathon. PanoptesBio separates boundary reasoning (Meta Llama 3.1) from calibrated biophysical math (pair-interaction GP ensemble) and enforces a fail-closed SQLite micro-USD spend gate before touching Modal A10G GPUs.
          </p>

          <div style={{ display: 'flex', gap: '12px', flexWrap: 'wrap', alignItems: 'center' }}>
            <a
              href="https://github.com/fugeretrevin/VeritasBio"
              target="_blank"
              rel="noopener noreferrer"
              className="btn btn-primary"
              style={{ padding: '11px 20px', fontSize: '13.5px', textDecoration: 'none' }}
            >
              <GitBranch size={15} />
              <span>View GitHub Repository</span>
              <ExternalLink size={13} />
            </a>

            <div
              style={{
                display: 'inline-flex',
                alignItems: 'center',
                gap: '6px',
                padding: '9px 14px',
                borderRadius: '8px',
                background: 'var(--bg-card)',
                fontSize: '12.5px',
                fontWeight: 600,
                color: '#0d9488',
              }}
            >
              <ShieldCheck size={15} />
              <span>192+ Automated Tests Passing</span>
            </div>
          </div>
        </div>

        {/* Right Telemetry Matrix (Tonal Surface Cards, No 1px Borders) */}
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '14px' }}>
          <div style={{ padding: '22px', borderRadius: '14px', background: 'var(--bg-elevated)' }}>
            <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginBottom: '4px' }}>
              48-Seed Reliability
            </div>
            <div className="mono tabular" style={{ fontSize: '30px', fontWeight: 600, color: '#0d9488' }}>
              98%
            </div>
            <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>
              47/48 seeds hit ≥ 0.70 (vs 58% ALDE)
            </div>
          </div>

          <div style={{ padding: '22px', borderRadius: '14px', background: 'var(--bg-elevated)' }}>
            <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginBottom: '4px' }}>
              Half-Budget Success
            </div>
            <div className="mono tabular" style={{ fontSize: '30px', fontWeight: 600, color: 'var(--text-main)' }}>
              88%
            </div>
            <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>
              In 50 assays (vs 4% ALDE / 6% Ridge)
            </div>
          </div>

          <div style={{ padding: '22px', borderRadius: '14px', background: 'var(--bg-elevated)' }}>
            <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginBottom: '4px' }}>
              4-Target Live GPU Spend
            </div>
            <div className="mono tabular" style={{ fontSize: '30px', fontWeight: 600, color: 'var(--text-main)' }}>
              $0.0376
            </div>
            <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>
              160 ESM-2 650M passes on Modal A10G
            </div>
          </div>

          <div style={{ padding: '22px', borderRadius: '14px', background: 'var(--bg-elevated)' }}>
            <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', marginBottom: '4px' }}>
              Global Peak Discovery
            </div>
            <div className="mono tabular" style={{ fontSize: '30px', fontWeight: 600, color: '#0d9488' }}>
              1.0000
            </div>
            <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '4px' }}>
              Rank #1 / 159,999 (V183A:F184I:V227K:S228G)
            </div>
          </div>
        </div>
      </div>

      {/* Interactive System Blueprint */}
      <div>
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'flex-end',
            flexWrap: 'wrap',
            gap: '12px',
            marginBottom: '24px',
          }}
        >
          <div>
            <h2 style={{ fontSize: '24px', fontWeight: 600, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
              Interactive System Blueprint
            </h2>
            <p style={{ color: 'var(--text-muted)', fontSize: '14.5px', margin: 0 }}>
              Click any node in the unidirectional 5-stage pipeline to inspect its exact input, output, and fail-closed guarantee.
            </p>
          </div>
        </div>

        {/* Visual Node Diagram */}
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(185px, 1fr))',
            gap: '12px',
            marginBottom: '18px',
          }}
        >
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
                  padding: '20px',
                  borderRadius: '14px',
                  cursor: 'pointer',
                  background: isSelected ? '#141413' : 'var(--bg-elevated)',
                  color: isSelected ? '#fcfbf8' : 'var(--text-main)',
                  border: 'none',
                  transition: 'all 0.18s ease',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-between',
                  gap: '16px',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div
                    style={{
                      width: '36px',
                      height: '36px',
                      borderRadius: '8px',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      background: isSelected ? 'rgba(255,255,255,0.12)' : 'var(--bg-card)',
                      color: isSelected ? '#5eead4' : '#0d9488',
                    }}
                  >
                    <IconComponent size={18} />
                  </div>
                  <ArrowRight size={15} style={{ opacity: isSelected ? 0.9 : 0.35 }} />
                </div>
                <div>
                  <div className="mono" style={{ fontSize: '11px', opacity: isSelected ? 0.8 : 0.65, marginBottom: '4px' }}>
                    {node.role}
                  </div>
                  <div style={{ fontSize: '15.5px', fontWeight: 600, lineHeight: 1.25 }}>
                    {node.title}
                  </div>
                </div>
              </button>
            );
          })}
        </div>

        {/* Active Node Inspection Surface */}
        <div
          style={{
            padding: '28px',
            borderRadius: '16px',
            background: 'var(--bg-card)',
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))',
            gap: '24px',
            alignItems: 'center',
          }}
        >
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '10px', flexWrap: 'wrap' }}>
              <ActiveIcon size={18} color="#0d9488" />
              <span style={{ fontSize: '18px', fontWeight: 600, color: 'var(--text-main)' }}>
                {activeNode.title}
              </span>
              <span
                className="mono"
                style={{
                  fontSize: '11.5px',
                  padding: '3px 9px',
                  borderRadius: '6px',
                  background: 'var(--bg-elevated)',
                  color: 'var(--text-muted)',
                }}
              >
                {activeNode.module}
              </span>
            </div>
            <p style={{ fontSize: '14px', color: 'var(--text-muted)', lineHeight: 1.6, margin: 0 }}>
              {activeNode.summary}
            </p>
          </div>

          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(160px, 1fr))', gap: '12px' }}>
            <div style={{ padding: '16px', borderRadius: '10px', background: 'var(--bg-elevated)' }}>
              <div style={{ fontSize: '11.5px', color: 'var(--text-dim)', marginBottom: '4px' }}>Incoming Data</div>
              <div className="mono" style={{ fontSize: '12.5px', fontWeight: 600, color: 'var(--text-main)' }}>
                {activeNode.input}
              </div>
            </div>
            <div style={{ padding: '16px', borderRadius: '10px', background: 'var(--bg-elevated)' }}>
              <div style={{ fontSize: '11.5px', color: 'var(--text-dim)', marginBottom: '4px' }}>Verified Output</div>
              <div className="mono" style={{ fontSize: '12.5px', fontWeight: 600, color: '#0d9488' }}>
                {activeNode.output}
              </div>
            </div>
            <div style={{ padding: '16px', borderRadius: '10px', background: 'var(--bg-elevated)' }}>
              <div style={{ fontSize: '11.5px', color: 'var(--text-dim)', marginBottom: '4px' }}>System Guarantee</div>
              <div style={{ fontSize: '12.5px', fontWeight: 600, color: 'var(--text-main)', display: 'flex', alignItems: 'center', gap: '6px' }}>
                <CheckCircle2 size={14} color="#0d9488" />
                <span>{activeNode.guarantee}</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* 4-Target Benchmark Matrix */}
      <div>
        <div style={{ marginBottom: '24px' }}>
          <h2 style={{ fontSize: '24px', fontWeight: 600, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
            4-Target Live Production Benchmark Matrix
          </h2>
          <p style={{ color: 'var(--text-muted)', fontSize: '14.5px', margin: 0 }}>
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
                  padding: '24px',
                  borderRadius: '14px',
                  background: 'var(--bg-elevated)',
                  display: 'flex',
                  flexDirection: 'column',
                  justifyContent: 'space-between',
                  gap: '16px',
                }}
              >
                <div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '6px' }}>
                    <span className="mono" style={{ fontSize: '11.5px', fontWeight: 600, color: '#0d9488' }}>
                      {item.pdb}
                    </span>
                    <span className="mono" style={{ fontSize: '11.5px', color: 'var(--text-muted)' }}>
                      Spend: {item.spend}
                    </span>
                  </div>
                  <h3 style={{ fontSize: '17px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '2px' }}>
                    {item.name}
                  </h3>
                  <div style={{ fontSize: '12.5px', color: 'var(--text-dim)' }}>{item.space}</div>
                </div>

                {/* Visual Baseline vs Lead Comparison Bars */}
                <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
                  <div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '11.5px', color: 'var(--text-muted)', marginBottom: '4px' }}>
                      <span>Wild-Type Baseline</span>
                      <span className="mono tabular">{item.wtScore}</span>
                    </div>
                    <div style={{ height: '8px', borderRadius: '999px', background: 'var(--bg-card)', overflow: 'hidden' }}>
                      <div style={{ width: `${wtPct}%`, height: '100%', background: '#8c8a84', borderRadius: '999px' }} />
                    </div>
                  </div>

                  <div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '11.5px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '4px' }}>
                      <span>Discovered Lead Peak</span>
                      <span className="mono tabular" style={{ color: '#0d9488' }}>
                        {item.leadScore} ({item.enrichment})
                      </span>
                    </div>
                    <div style={{ height: '10px', borderRadius: '999px', background: 'var(--bg-card)', overflow: 'hidden' }}>
                      <div style={{ width: `${leadPct}%`, height: '100%', background: '#0d9488', borderRadius: '999px' }} />
                    </div>
                  </div>
                </div>

                <div style={{ fontSize: '12.5px', color: 'var(--text-muted)', paddingTop: '12px', borderTop: '1px solid var(--border-subtle)' }}>
                  {item.highlight}
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Engineering Progression */}
      <div>
        <div style={{ marginBottom: '24px' }}>
          <h2 style={{ fontSize: '24px', fontWeight: 600, color: 'var(--text-main)', letterSpacing: '-0.02em', marginBottom: '6px' }}>
            Engineering Progression
          </h2>
          <p style={{ color: 'var(--text-muted)', fontSize: '14.5px', margin: 0 }}>
            How PanoptesBio evolved from an immutable micro-USD ledger to a 48-seed calibrated epistatic discovery platform.
          </p>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))', gap: '18px' }}>
          {sprintVisuals.map((phase) => (
            <div
              key={phase.day}
              style={{
                padding: '24px',
                borderRadius: '14px',
                background: 'var(--bg-elevated)',
                display: 'flex',
                flexDirection: 'column',
                gap: '14px',
              }}
            >
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                <span className="mono" style={{ fontSize: '12px', fontWeight: 600, color: '#0d9488' }}>
                  {phase.day}
                </span>
                <span
                  className="mono"
                  style={{
                    fontSize: '11.5px',
                    padding: '3px 10px',
                    borderRadius: '6px',
                    background: 'var(--bg-card)',
                    color: 'var(--text-main)',
                    fontWeight: 600,
                  }}
                >
                  {phase.metric}
                </span>
              </div>

              <div style={{ padding: '12px', borderRadius: '10px', background: 'var(--bg-dark)' }}>
                {phase.svg}
              </div>

              <div>
                <h3 style={{ fontSize: '16.5px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '6px' }}>
                  {phase.title}
                </h3>
                <p style={{ fontSize: '13.5px', color: 'var(--text-muted)', lineHeight: 1.55, margin: 0 }}>
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
