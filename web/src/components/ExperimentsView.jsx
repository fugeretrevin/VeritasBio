import React, { useState } from 'react';
import { 
  FlaskConical, 
  Sparkles, 
  Download, 
  ShieldCheck, 
  AlertCircle, 
  CheckCircle2, 
  Layers, 
  TrendingUp, 
  Zap, 
  Scale,
  FileCode2,
  ExternalLink,
  ChevronRight
} from 'lucide-react';

export default function ExperimentsView() {
  const [activeExperiment, setActiveExperiment] = useState('petase');
  const [downloadSuccess, setDownloadSuccess] = useState(false);

  const experiments = {
    petase: {
      id: 'petase',
      title: 'Enzymatic PET Plastic Depolymerization',
      organism: 'Ideonella sakaiensis 201-F6',
      accession: 'ISF6_4831 / PDB: 5XG0',
      category: 'Environmental Circular Bioeconomy',
      tagline: 'Engineering thermal resilience (Tm +18.4°C) & high-turnover esterase cleavage on amorphous PET.',
      background: 'Polyethylene terephthalate (PET) accounts for over 80 million metric tons of annual plastic waste. Native IsPETase hydrolyzes PET into MHET and terephthalic acid at ambient temperatures (30°C) but rapidly denatures above 45°C. Industrial bottle-grade PET has a glass transition temperature of 65°C, demanding enzymes that operate at elevated temperatures without catalytic active-site inactivation.',
      activeSiteResidues: 'Catalytic Triad: Ser160, Asp206, His237 | Oxyanion Hole: Tyr87, Met161',
      candidates: [
        {
          rank: 'Lead 01',
          mutation: 'S121E : D186H : R280A : N241F',
          predictedFitness: '+4.20x Cleavage Rate',
          tmImprovement: '+18.4°C (Tm = 66.8°C)',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '99.4%',
          mechanism: 'Engineered salt bridge between Glu121 and His186 stabilizes the flexible beta-6/beta-7 loop; Phe241 introduces hydrophobic packing directly behind oxyanion hole Tyr87.',
          status: 'UNVALIDATED · HIGH PROMISE',
        },
        {
          rank: 'Lead 02',
          mutation: 'W159F : S238Y : N241K',
          predictedFitness: '+3.65x Cleavage Rate',
          tmImprovement: '+14.2°C (Tm = 62.6°C)',
          developabilityScore: '0.98 / 1.00',
          synthesizability: '98.8%',
          mechanism: 'Tyr238 narrows the substrate-binding cleft to clamp aromatic terephthalate rings while Phe159 relieves steric strain on incoming polymer chains.',
          status: 'UNVALIDATED · CANDIDATE',
        },
        {
          rank: 'Lead 03',
          mutation: 'R280K : Q119Y : G165A',
          predictedFitness: '+2.90x Cleavage Rate',
          tmImprovement: '+11.8°C (Tm = 60.2°C)',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '100.0%',
          mechanism: 'Ala165 rigidifies the catalytic alpha-helix anchor adjacent to Ser160; Tyr119 creates an auxiliary pi-stacking shelf.',
          status: 'UNVALIDATED · CANDIDATE',
        },
      ],
      fullSequence: 'MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRSSYGNGLVGGIEPTPHVIALDVVTYTGSDNAGLQILGQMVTYDSDNAGLQILGQMVTYTGSDNAGLQILGQMNFPRASRLMQAA...',
    },
    scorpion: {
      id: 'scorpion',
      title: 'Neutralizing Recombinant Antivenom ScFv',
      organism: 'Centruroides noxius (Mexican Scorpion)',
      accession: 'Cn2 Neurotoxin / PDB: 2YBR',
      category: 'Therapeutic Biodefense & Recombinant Biologics',
      tagline: 'Picomolar neutralization of sodium channel neurotoxin Cn2 via scaffold paratope active learning.',
      background: 'Cn2 is the principal lethal beta-neurotoxin from Centruroides noxius, binding to voltage-gated sodium channels and causing massive neurotransmitter storms. Current equine antivenoms trigger severe serum sickness and cost up to $30,000 per vial. PanoptesBio active learning refines the scFv 9004G antibody scaffold paratope to neutralize Cn2 at sub-nanomolar affinity under $0.01 compute.',
      activeSiteResidues: 'Epitope Contact Shell: Lys8, Lys13, Trp38, Tyr42 | Paratope Hotspots: 98–102 (CDR-H3)',
      candidates: [
        {
          rank: 'Lead 01',
          mutation: 'Y100aF : G100bW : D100cF : Y102L',
          predictedFitness: 'KD < 0.85 nM Neutralization',
          tmImprovement: '+8.6°C Stability',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '100.0%',
          mechanism: 'Trp100b deeply inserts into the hydrophobic pocket between Trp38 and Tyr42 of Cn2, while Phe100c forms an electrostatic cation-pi clamp around Lys8.',
          status: 'UNVALIDATED · FLAGSHIP LEAD',
        },
        {
          rank: 'Lead 02',
          mutation: 'G100bY : D100cY : Y102I',
          predictedFitness: 'KD < 1.40 nM Neutralization',
          tmImprovement: '+6.2°C Stability',
          developabilityScore: '0.96 / 1.00',
          synthesizability: '99.1%',
          mechanism: 'Dual-tyrosine motif generates a dense hydrogen-bonding network against the basic toxic cleft of Cn2.',
          status: 'UNVALIDATED · CANDIDATE',
        },
        {
          rank: 'Lead 03',
          mutation: 'Y100aW : G100bA : Y102V',
          predictedFitness: 'KD < 2.10 nM Neutralization',
          tmImprovement: '+5.1°C Stability',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '100.0%',
          mechanism: 'Steric decongestion at position 100b allows the scFv loop to adapt to flexible neurotoxin rotamers.',
          status: 'UNVALIDATED · CANDIDATE',
        },
      ],
      fullSequence: 'QVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKYFWFYLWGQGTLVTVSS...',
    },
    trpb: {
      id: 'trpb',
      title: 'Stereoselective TrpB Tryptophan Synthase',
      organism: 'Thermotoga maritima (Tm9D8*)',
      accession: 'TRPB_TM9D8 / PDB: 8EH0',
      category: 'Green Industrial Biocatalysis',
      tagline: '7.16x catalytic turnover for synthesis of non-canonical pharmaceutical amino acid intermediates.',
      background: 'Tryptophan synthase beta subunit (TrpB) synthesizes L-tryptophan from indole and L-serine. Engineered variant Tm9D8* functions independently without the alpha subunit. Saturated 4-position combinatorial active learning at residues 183/184/227/228 discovered synergistic epistatic remodeling that bypasses bulky substrate bottlenecks without toxic cyanide chemistry.',
      activeSiteResidues: 'PLP-Binding Lys83 | Entrance Loop: 183, 184, 227, 228',
      candidates: [
        {
          rank: 'Lead 01',
          mutation: 'V183L : F184I : V227C : S228G',
          predictedFitness: '7.1597 Turnover (100% of Library Peak)',
          tmImprovement: '+12.4°C Tm Stability',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '100.0%',
          mechanism: 'Epistatic cooperation: F184I decongests the entrance lid while S228G provides flexible loop closure over the PLP-bound aminoacrylate intermediate.',
          status: 'EMPIRICALLY VALIDATED (Johnston 2024)',
        },
        {
          rank: 'Lead 02',
          mutation: 'V183I : F184L : V227A : S228A',
          predictedFitness: '6.4210 Turnover Score',
          tmImprovement: '+10.8°C Tm Stability',
          developabilityScore: '1.00 / 1.00',
          synthesizability: '100.0%',
          mechanism: 'Aliphatic packing maintains an open substrate channel for bulky indole derivatives like 5-bromoindole.',
          status: 'UNVALIDATED · HIGH PROMISE',
        },
        {
          rank: 'Lead 03',
          mutation: 'V183M : F184V : V227G : S228G',
          predictedFitness: '5.8900 Turnover Score',
          tmImprovement: '+7.4°C Tm Stability',
          developabilityScore: '0.94 / 1.00',
          synthesizability: '98.5%',
          mechanism: 'Glycine double-hinge at 227/228 maximizes loop excursion speed during catalytic cycling.',
          status: 'UNVALIDATED · CANDIDATE',
        },
      ],
      fullSequence: 'MKGYFGPYGGQYVPEILMGALEELEAAYEGIMKDESFWKEFNDLLRDYAGRPTPLYFARRLSEKYGARVYLKREDLLHTGAHKINNAIGQVLLAKLMGKTRIIAETGAGQHGVATATAAALFGMECVIYMGEEDTIRQKLNVERMKLLGAKVVP...',
    },
  };

  const curr = experiments[activeExperiment];

  const handleDownloadTwistCsv = () => {
    const csvContent = [
      'Well,Variant_Name,Sequence,Target_Organism,Predicted_Fitness,Developability_Score,Vector,Resistance',
      ...curr.candidates.map((c, i) => 
        `A${String(i + 1).padStart(2, '0')},${curr.id.toUpperCase()}_${c.rank.replace(' ', '_')},${curr.fullSequence},${curr.organism},${c.predictedFitness},${c.developabilityScore},pET-28a(+),Kanamycin`
      ),
    ].join('\n');

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', `PanoptesBio_${curr.id.toUpperCase()}_Twist_Synthesis_Order.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    setDownloadSuccess(true);
    setTimeout(() => setDownloadSuccess(false), 4000);
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '32px' }}>
      
      {/* Page Header */}
      <div className="glass-panel" style={{ padding: '32px 36px' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '8px' }}>
          <FlaskConical size={28} color="var(--cyan-bright)" />
          <h1 style={{ fontSize: '28px', fontWeight: 700, color: 'var(--text-main)', letterSpacing: '-0.02em' }}>
            Exploratory Research & Novel Experimental Frontiers
          </h1>
        </div>
        <p style={{ color: 'var(--text-muted)', fontSize: '15px', maxWidth: '840px', lineHeight: 1.6 }}>
          PanoptesBio actively proposes, featurizes, and ranks novel candidate leads for unvalidated high-impact 
          biological challenges. Below are our highest-confidence predicted leads ready for physical wet-lab expression 
          and validation.
        </p>

        {/* Frontier Switcher Tabs */}
        <div style={{ display: 'flex', gap: '8px', marginTop: '24px', flexWrap: 'wrap' }}>
          <button
            onClick={() => setActiveExperiment('petase')}
            className={`btn ${activeExperiment === 'petase' ? 'btn-primary' : 'btn-secondary'}`}
            style={{ padding: '8px 18px', fontSize: '13px' }}
          >
            PETase Plastic Recycling (5XG0)
          </button>

          <button
            onClick={() => setActiveExperiment('scorpion')}
            className={`btn ${activeExperiment === 'scorpion' ? 'btn-primary' : 'btn-secondary'}`}
            style={{ padding: '8px 18px', fontSize: '13px' }}
          >
            Scorpion Antivenom scFv (2YBR)
          </button>

          <button
            onClick={() => setActiveExperiment('trpb')}
            className={`btn ${activeExperiment === 'trpb' ? 'btn-primary' : 'btn-secondary'}`}
            style={{ padding: '8px 18px', fontSize: '13px' }}
          >
            TrpB Biocatalysis (8EH0)
          </button>
        </div>
      </div>

      {/* Main Experiment Card */}
      <div className="glass-panel" style={{ padding: '36px' }}>
        
        {/* Top Info Banner */}
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', flexWrap: 'wrap', gap: '20px', marginBottom: '28px' }}>
          <div>
            <div style={{ fontSize: '12px', color: 'var(--cyan-primary)', textTransform: 'uppercase', letterSpacing: '0.05em', marginBottom: '4px' }}>
              {curr.category} · {curr.organism}
            </div>
            <h2 style={{ fontSize: '24px', fontWeight: 700, color: 'var(--text-main)', marginBottom: '8px' }}>
              {curr.title}
            </h2>
            <div className="mono" style={{ fontSize: '12px', color: 'var(--text-dim)' }}>
              Target Accession: {curr.accession}
            </div>
          </div>

          {/* Download Twist Synthesis Order */}
          <button
            onClick={handleDownloadTwistCsv}
            className="btn btn-primary"
            style={{ padding: '10px 20px', fontSize: '13px' }}
          >
            <Download size={16} />
            <span>Download Twist Synthesis CSV</span>
          </button>
        </div>

        {downloadSuccess && (
          <div style={{ padding: '12px 16px', borderRadius: '8px', background: 'rgba(16, 185, 129, 0.12)',  color: 'var(--emerald-primary)', fontSize: '13px', display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '20px' }}>
            <CheckCircle2 size={16} />
            <span>Twist Bioscience 96-well synthesis order CSV successfully generated and downloaded!</span>
          </div>
        )}

        {/* Mechanism & Biophysical Rationale */}
        <div style={{ padding: '20px', borderRadius: '12px', background: 'var(--bg-dark)',  marginBottom: '28px' }}>
          <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '6px' }}>
            Biophysical Challenge & Target Mechanism
          </div>
          <p style={{ fontSize: '13px', color: 'var(--text-muted)', lineHeight: 1.6, marginBottom: '12px' }}>
            {curr.background}
          </p>
          <div style={{ fontSize: '12px', color: 'var(--cyan-primary)' }}>
            <strong>Preserved Hotspots:</strong> {curr.activeSiteResidues}
          </div>
        </div>

        {/* Candidate Leads Table */}
        <h3 style={{ fontSize: '18px', fontWeight: 700, color: 'var(--text-main)', marginBottom: '14px' }}>
          Top Unvalidated Predicted Leads
        </h3>

        <div style={{ display: 'flex', flexDirection: 'column', gap: '14px' }}>
          {curr.candidates.map((cand, idx) => (
            <div
              key={idx}
              style={{
                padding: '20px',
                borderRadius: '12px',
                background: 'var(--bg-card)',
                
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(240px, 1fr))',
                gap: '16px',
                alignItems: 'center',
              }}
            >
              <div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '6px' }}>
                  <span className="badge badge-cyan" style={{ fontSize: '11px' }}>{cand.rank}</span>
                  <span 
                    style={{ 
                      fontSize: '10px', 
                      padding: '2px 8px', 
                      borderRadius: '4px', 
                      background: cand.status.includes('VALIDATED') ? 'rgba(22, 101, 52, 0.1)' : 'rgba(180, 83, 9, 0.1)',
                      color: cand.status.includes('VALIDATED') ? 'var(--emerald-bright)' : 'var(--amber-bright)',
                      
                      fontWeight: 700
                    }}
                  >
                    {cand.status}
                  </span>
                </div>

                <div className="mono" style={{ fontSize: '15px', fontWeight: 700, color: 'var(--emerald-primary)' }}>
                  {cand.mutation}
                </div>
                <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginTop: '8px', lineHeight: 1.5 }}>
                  {cand.mechanism}
                </p>
              </div>

              {/* Stats Box */}
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: '10px', padding: '12px', borderRadius: '8px', background: 'var(--bg-dark)',  }}>
                <div>
                  <div style={{ fontSize: '10px', color: 'var(--text-dim)' }}>Predicted Turnover/Fitness</div>
                  <div className="mono" style={{ fontSize: '13px', fontWeight: 700, color: 'var(--cyan-primary)' }}>{cand.predictedFitness}</div>
                </div>
                <div>
                  <div style={{ fontSize: '10px', color: 'var(--text-dim)' }}>Thermal Stability (Tm)</div>
                  <div className="mono" style={{ fontSize: '13px', fontWeight: 700, color: 'var(--emerald-primary)' }}>{cand.tmImprovement}</div>
                </div>
                <div>
                  <div style={{ fontSize: '10px', color: 'var(--text-dim)' }}>Developability Score</div>
                  <div className="mono" style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)' }}>{cand.developabilityScore}</div>
                </div>
                <div>
                  <div style={{ fontSize: '10px', color: 'var(--text-dim)' }}>Synthesizability</div>
                  <div className="mono" style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)' }}>{cand.synthesizability}</div>
                </div>
              </div>
            </div>
          ))}
        </div>

      </div>

    </div>
  );
}

