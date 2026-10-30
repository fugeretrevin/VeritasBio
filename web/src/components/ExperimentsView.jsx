import React, { useState } from 'react';
import {
  FlaskConical,
  CheckCircle2,
  ArrowRight,
  FileSpreadsheet,
  Bot,
  Activity,
} from 'lucide-react';

export default function ExperimentsView() {
  const [selectedCampaign, setSelectedCampaign] = useState(0);
  const [selectedLeadIdx, setSelectedLeadIdx] = useState(0);

  const campaigns = [
    {
      id: 'SCFV_9004G_CN2',
      pdbId: '2YBR',
      title: 'Scorpion Antivenom Neutralizing scFv (9004G · Toxin Cn2)',
      subtitle: 'Immunotherapy Biologic · Interface Affinity & Thermal Maturation',
      summary:
        'Optimizes CDR-H3 paratope contact density against Centruroides noxius sodium-channel neurotoxin Cn2 while eliminating clinical deamidation and aggregation liabilities.',
      wtMetric: 'Kd = 1.8 nM (Tm = 62.4 °C)',
      targetMetric: 'Kd ≤ 95 pM (Tm = 74.1 °C)',
      barWt: 22,
      barLead: 92,
      leads: [
        {
          mutation: 'Y100W : S101A : G102V',
          predictedGain: '19x Affinity Gain (Kd ≈ 95 pM)',
          developability: '100% CMC Clean · Zero Unpaired Cys',
          rationale:
            'W100 deepens π-cation packing against Cn2 Lys35 while A101/V102 rigidify the CDR-H3 turn, reducing entropic binding penalty.',
        },
        {
          mutation: 'R98K : Y100F : S101T : G102A',
          predictedGain: '11x Affinity Gain · +8.2 °C ΔTm',
          developability: 'Zero NG/NS Deamidation Motifs',
          rationale:
            'Balances sub-nanomolar toxin neutralization with improved resistance to thermal aggregation for ambient-temperature field storage.',
        },
        {
          mutation: 'Y100W : D101E : G102V',
          predictedGain: '15x Affinity Gain · Isomerization-Free',
          developability: 'Removes DG101 Isomerization Liability',
          rationale:
            'Substitutes the Asp101-Gly102 chemical degradation hotspot with Glu101 while preserving salt-bridge geometry.',
        },
      ],
      protocolSteps: [
        {
          step: 'Phase I',
          title: 'Twist Golden Gate cDNA Assembly',
          detail: 'Order 380 bp BsaI-flanked variable fragments codon-optimized for E. coli BL21(DE3) periplasmic expression.',
        },
        {
          step: 'Phase II',
          title: 'Opentrons OT-2 96-Well Expression',
          detail: 'Automated colony inoculation and IPTG induction with 36 perimeter H₂O hydration barrier wells to prevent edge evaporation.',
        },
        {
          step: 'Phase III',
          title: 'BLI / SPR Off-Rate Screening',
          detail: 'Measure kon / koff kinetics against biotinylated recombinant Cn2 toxin and ingest CSV curves back into PanoptesBio.',
        },
      ],
    },
    {
      id: 'TRPB_TM9D8',
      pdbId: '8EH0',
      title: 'Tryptophan Synthase β-Subunit (Tm9D8* Active-Site Pocket)',
      subtitle: 'Biocatalysis · 4-Site Sign-Epistatic Landscape (159,999 Variants)',
      summary:
        'Validates our qualify-then-optimize frontier acquisition rule on the complete 20⁴ combinatorial pocket (residues 183, 184, 227, 228) surrounding the PLP cofactor.',
      wtMetric: '0.402 Relative Catalytic Activity',
      targetMetric: '1.0000 Global Peak (#1 / 159,999)',
      barWt: 40,
      barLead: 100,
      leads: [
        {
          mutation: 'V183A : F184I : V227K : S228G',
          predictedGain: '1.0000 Global Optimum (2.49x over WT)',
          developability: 'Zero PLP Cofactor Clashes (> 4.2 Å)',
          rationale:
            'Unlocks the K227 × G228 sign-epistasis flip: G228 clears active-site steric volume so K227 flips from -0.50 deleterious to +0.82 catalytic gain.',
        },
        {
          mutation: 'V183I : F184L : V227K : S228G',
          predictedGain: '0.9410 Activity (Rank #2 Peak)',
          developability: 'High Thermal Stability (Tm > 78 °C)',
          rationale:
            'Hydrophobic core packing at 183I/184L stabilizes the β-subunit dimer interface alongside the K227:G228 catalytic dyad.',
        },
        {
          mutation: 'V183I : F184V : V227M : S228G',
          predictedGain: '0.6192 Activity (98% Seed Floor)',
          developability: 'Guaranteed Round-2 Convergence',
          rationale:
            'Serves as the ultra-reliable baseline scaffold reached by 47/48 random starts before frontier exploration probes position 227.',
        },
      ],
      protocolSteps: [
        {
          step: 'Phase I',
          title: 'Site-Directed Golden Gate Library',
          detail: 'Synthesize top 24 Pareto-optimal TrpB cassettes with C-terminal 6xHis purification tags.',
        },
        {
          step: 'Phase II',
          title: 'PLP-Supplemented Lysis & Reaction',
          detail: 'Opentrons OT-2 dispenses indole + L-serine substrate mix with 100 µM pyridoxal phosphate cofactor.',
        },
        {
          step: 'Phase III',
          title: 'UV-Vis Spectrophotometric Readout',
          detail: 'Quantify L-tryptophan formation at 290 nm with automated Grubbs/Huber outlier filtering and WT B02 normalization.',
        },
      ],
    },
    {
      id: 'HER2_TRASTUZUMAB',
      pdbId: '1N8Z',
      title: 'Trastuzumab Fab CDR-H3 (HER2 Extracellular Domain)',
      subtitle: 'Oncology Biologic · Affinity & CMC Liability Co-Optimization',
      summary:
        'Eliminates the clinical Asp101-Gly102 isomerization liability in Trastuzumab CDR-H3 while improving electrostatic complementarity against HER2 Domain IV.',
      wtMetric: '+0.42 Baseline (Contains DG101 Liability)',
      targetMetric: '+1.45 ΔΔG · Zero Chemical Liabilities',
      barWt: 29,
      barLead: 94,
      leads: [
        {
          mutation: 'W99F : G100A : D101E : G102Y',
          predictedGain: '+1.45 ΔΔG · Rank 0 Pareto Optimal',
          developability: '100% CMC Clean · pI = 8.4',
          rationale:
            'Replaces the isomerization-prone D101-G102 dipeptide with E101-Y102 and relaxes the Y102 aromatic ring via CPU L-BFGS.',
        },
        {
          mutation: 'W99Y : G100S : D101E : G102V',
          predictedGain: '+1.28 ΔΔG · Low Viscosity Score',
          developability: 'Reduced Hydrophobic Patch Index',
          rationale:
            'Minimizes surface-exposed aromaticity to support high-concentration (>150 mg/mL) subcutaneous formulation.',
        },
      ],
      protocolSteps: [
        {
          step: 'Phase I',
          title: 'CHO-K1 Heavy-Chain Variable Synthesis',
          detail: 'Mammalian codon-optimized VH fragments cloned into human IgG1 expression vectors.',
        },
        {
          step: 'Phase II',
          title: 'Automated Protein A Micro-Purification',
          detail: '96-well magnetic bead capture and pH 4.0 thermal stress challenge to verify zero isomerization.',
        },
        {
          step: 'Phase III',
          title: '4PL ELISA &Octet HER2 Titration',
          detail: 'Ingest 8-point serial dilution curves into PanoptesBio to fit IC50 / Kd and update the GP posterior.',
        },
      ],
    },
  ];

  const activeCampaign = campaigns[selectedCampaign];
  const activeLead = activeCampaign.leads[selectedLeadIdx] || activeCampaign.leads[0];

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '48px' }}>
      {/* Editorial Header + Target Campaign Switcher */}
      <div
        style={{
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'flex-end',
          flexWrap: 'wrap',
          gap: '24px',
          paddingBottom: '24px',
          borderBottom: '1px solid var(--border-subtle)',
        }}
      >
        <div>
          <h1
            className="title-display"
            style={{ fontSize: 'var(--step-head)', color: 'var(--text-main)', marginBottom: '10px' }}
          >
            Prospective wet-lab{' '}
            <span style={{ color: '#0d9488' }}>
              validation dossiers.
            </span>
          </h1>
          <p style={{ fontSize: 'var(--step-body)', color: 'var(--text-muted)', maxWidth: '64ch' }}>
            Synthesis-ready lead candidates generated by PanoptesBio across neutralizing antivenoms, industrial biocatalysts, and therapeutic antibodies, complete with Twist cDNA designs and Opentrons assay protocols.
          </p>
        </div>

        {/* Target Switcher Pills */}
        <div
          style={{
            display: 'inline-flex',
            gap: '6px',
            padding: '5px',
            borderRadius: '10px',
            background: 'var(--bg-card)',
          }}
        >
          {campaigns.map((camp, idx) => (
            <button
              key={camp.id}
              type="button"
              onClick={() => {
                setSelectedCampaign(idx);
                setSelectedLeadIdx(0);
              }}
              style={{
                padding: '8px 14px',
                borderRadius: '7px',
                border: 'none',
                cursor: 'pointer',
                fontSize: '13px',
                fontWeight: selectedCampaign === idx ? 600 : 500,
                background: selectedCampaign === idx ? '#141413' : 'transparent',
                color: selectedCampaign === idx ? '#fcfbf8' : 'var(--text-muted)',
              }}
            >
              {camp.pdbId} · {camp.id.split('_')[0]}
            </button>
          ))}
        </div>
      </div>

      {/* Main Campaign Dossier Split: Left Overview + Visual Gain, Right Interactive Lead Inspector */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))',
          gap: '32px',
          alignItems: 'start',
        }}
      >
        {/* Left: Target Summary & Visual Baseline-to-Lead Comparison */}
        <div
          style={{
            background: 'var(--bg-elevated)',
            borderRadius: '16px',
            padding: '32px',
            display: 'flex',
            flexDirection: 'column',
            gap: '24px',
          }}
        >
          <div>
            <div className="mono" style={{ fontSize: '12px', fontWeight: 600, color: '#0d9488', marginBottom: '6px' }}>
              PDB: {activeCampaign.pdbId} · {activeCampaign.subtitle}
            </div>
            <h2 style={{ fontSize: '23px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '12px', lineHeight: 1.25 }}>
              {activeCampaign.title}
            </h2>
            <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
              {activeCampaign.summary}
            </p>
          </div>

          {/* Visual Baseline vs Lead Bars */}
          <div
            style={{
              background: 'var(--bg-card)',
              borderRadius: '12px',
              padding: '22px',
              display: 'flex',
              flexDirection: 'column',
              gap: '16px',
            }}
          >
            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12.5px', marginBottom: '6px' }}>
                <span style={{ color: 'var(--text-muted)' }}>Wild-Type Parent Baseline</span>
                <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
                  {activeCampaign.wtMetric}
                </span>
              </div>
              <div style={{ height: '8px', borderRadius: '999px', background: 'var(--bg-elevated)', overflow: 'hidden' }}>
                <div style={{ width: `${activeCampaign.barWt}%`, height: '100%', background: '#8c8a84', borderRadius: '999px' }} />
              </div>
            </div>

            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12.5px', marginBottom: '6px' }}>
                <span style={{ fontWeight: 600, color: 'var(--text-main)' }}>Rank-0 Engineered Lead</span>
                <span className="mono tabular" style={{ fontWeight: 600, color: '#0d9488' }}>
                  {activeCampaign.targetMetric}
                </span>
              </div>
              <div style={{ height: '10px', borderRadius: '999px', background: 'var(--bg-elevated)', overflow: 'hidden' }}>
                <div style={{ width: `${activeCampaign.barLead}%`, height: '100%', background: '#0d9488', borderRadius: '999px' }} />
              </div>
            </div>
          </div>
        </div>

        {/* Right: Interactive Lead Genotype Selector */}
        <div
          style={{
            background: 'var(--bg-elevated)',
            borderRadius: '16px',
            padding: '32px',
            display: 'flex',
            flexDirection: 'column',
            gap: '20px',
          }}
        >
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '10px' }}>
            <h3 style={{ fontSize: '17px', fontWeight: 600, color: 'var(--text-main)' }}>
              Nominated Synthesis Candidates
            </h3>
            <div style={{ display: 'flex', gap: '6px' }}>
              {activeCampaign.leads.map((ld, idx) => (
                <button
                  key={ld.mutation}
                  type="button"
                  onClick={() => setSelectedLeadIdx(idx)}
                  className="mono"
                  style={{
                    padding: '6px 12px',
                    borderRadius: '6px',
                    border: 'none',
                    cursor: 'pointer',
                    fontSize: '12px',
                    fontWeight: selectedLeadIdx === idx ? 600 : 500,
                    background: selectedLeadIdx === idx ? '#141413' : 'var(--bg-card)',
                    color: selectedLeadIdx === idx ? '#fcfbf8' : 'var(--text-muted)',
                  }}
                >
                  Lead #{idx + 1}
                </button>
              ))}
            </div>
          </div>

          {/* Active Lead Card */}
          <div
            style={{
              background: 'var(--bg-card)',
              borderRadius: '12px',
              padding: '24px',
              display: 'flex',
              flexDirection: 'column',
              gap: '14px',
            }}
          >
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', flexWrap: 'wrap', gap: '8px' }}>
              <span className="mono" style={{ fontSize: '18px', fontWeight: 600, color: 'var(--text-main)' }}>
                {activeLead.mutation}
              </span>
              <span className="mono tabular" style={{ fontSize: '13.5px', fontWeight: 600, color: '#0d9488' }}>
                {activeLead.predictedGain}
              </span>
            </div>

            <p style={{ fontSize: '14.5px', color: 'var(--text-muted)', lineHeight: 1.6 }}>
              {activeLead.rationale}
            </p>

            <div
              style={{
                paddingTop: '12px',
                borderTop: '1px solid var(--border-subtle)',
                display: 'flex',
                alignItems: 'center',
                gap: '8px',
                fontSize: '13px',
                fontWeight: 600,
                color: '#0d9488',
              }}
            >
              <CheckCircle2 size={15} />
              <span>{activeLead.developability}</span>
            </div>
          </div>
        </div>
      </div>

      {/* Bottom 3-Column Wet-Lab Execution Timeline (Open Tonal Surface) */}
      <div
        style={{
          background: 'var(--bg-card)',
          borderRadius: '18px',
          padding: '32px',
        }}
      >
        <div style={{ marginBottom: '24px' }}>
          <h3 style={{ fontSize: '20px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '4px' }}>
            Standardized Wet-Lab Validation Protocol
          </h3>
          <p style={{ fontSize: '14px', color: 'var(--text-muted)' }}>
            Every nominated candidate exports directly to automated CRO synthesis and liquid-handling robotics.
          </p>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))', gap: '20px' }}>
          {activeCampaign.protocolSteps.map((st) => (
            <div
              key={st.step}
              style={{
                background: 'var(--bg-elevated)',
                borderRadius: '12px',
                padding: '22px',
                display: 'flex',
                flexDirection: 'column',
                gap: '8px',
              }}
            >
              <span className="mono" style={{ fontSize: '12px', fontWeight: 600, color: '#0d9488' }}>
                {st.step}
              </span>
              <div style={{ fontSize: '16px', fontWeight: 600, color: 'var(--text-main)' }}>
                {st.title}
              </div>
              <p style={{ fontSize: '13.5px', color: 'var(--text-muted)', lineHeight: 1.55 }}>
                {st.detail}
              </p>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
