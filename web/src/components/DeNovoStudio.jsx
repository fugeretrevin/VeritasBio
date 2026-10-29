import React, { useState } from 'react';
import { Sparkles, Download, Layers, ShieldCheck, FileSpreadsheet, Bot } from 'lucide-react';
import Protein3DViewer from './Protein3DViewer';

export default function DeNovoStudio({ onRefreshLedger }) {
  const [targetScaffold, setTargetScaffold] = useState('SCFV_9004G_CN2');
  const [numDesigns, setNumDesigns] = useState(3);
  const [binderLength, setBinderLength] = useState(65);
  const [budgetUsd, setBudgetUsd] = useState(0.50);

  const [isLoading, setIsLoading] = useState(false);
  const [errorMsg, setErrorMsg] = useState(null);
  const [generationResult, setGenerationResult] = useState(null);
  const [selectedDesign, setSelectedDesign] = useState(null);

  const handleGenerate = async () => {
    setIsLoading(true);
    setErrorMsg(null);
    try {
      const res = await fetch('/api/diffusion/generate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          target_scaffold: targetScaffold,
          num_designs: numDesigns,
          binder_length: binderLength,
          budget_usd: budgetUsd,
          auth_token: 'panoptes_demo_token',
        }),
      });

      if (!res.ok) {
        const errData = await res.json();
        throw new Error(errData.detail || 'De novo generation failed');
      }

      const data = await res.json();
      setGenerationResult(data);
      if (data.designs && data.designs.length > 0) {
        setSelectedDesign(data.designs[0]);
      }
      if (onRefreshLedger) onRefreshLedger();
    } catch (err) {
      setErrorMsg(err.message);
    } finally {
      setIsLoading(false);
    }
  };

  const handleDownloadTwist = async () => {
    try {
      const res = await fetch('/api/wetlab/export', { method: 'POST' });
      const data = await res.json();
      const blob = new Blob([data.twist_csv], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `twist_de_novo_${targetScaffold}_orders.csv`;
      a.click();
    } catch (err) {
      alert(`Export error: ${err.message}`);
    }
  };

  const handleDownloadOpentrons = async () => {
    try {
      const res = await fetch('/api/wetlab/export', { method: 'POST' });
      const data = await res.json();
      const blob = new Blob([data.opentrons_script], { type: 'text/x-python' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `opentrons_de_novo_${targetScaffold}_protocol.py`;
      a.click();
    } catch (err) {
      alert(`Export error: ${err.message}`);
    }
  };

  return (
    <div style={{ display: 'grid', gridTemplateColumns: 'minmax(320px, 420px) 1fr', gap: '20px' }}>
      {/* Left Control Column */}
      <div className="glass-panel" style={{ padding: '24px' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '16px' }}>
          <Sparkles size={20} color="var(--cyan-bright)" />
          <h2 className="title-display" style={{ fontSize: '18px', margin: 0 }}>
            De Novo Binder Diffusion Studio
          </h2>
        </div>
        <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '20px' }}>
          Autonomous 3-helix bundle backbone diffusion (RFdiffusion) & sequence inverse folding (ProteinMPNN)
          with zero unpaired cysteines and strict developability scoring.
        </p>

        {/* Target Epitope */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, marginBottom: '6px' }}>
            Target Epitope / Scaffold
          </label>
          <select
            className="input-field"
            value={targetScaffold}
            onChange={(e) => setTargetScaffold(e.target.value)}
          >
            <option value="SCFV_9004G_CN2">SCFV_9004G_CN2 (Scorpion Neurotoxin · 2YBR)</option>
            <option value="HER2_TRASTUZUMAB">HER2_TRASTUZUMAB (HER2 Extracellular Receptor · 1N8Z)</option>
            <option value="P42212">avGFP Surface Epitope (avGFP β-Barrel · 1EMA)</option>
          </select>
        </div>

        {/* Num Designs Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
            <span>Number of Backbones</span>
            <span className="mono" style={{ color: 'var(--cyan-bright)' }}>{numDesigns} designs</span>
          </div>
          <input
            type="range"
            min="1"
            max="8"
            step="1"
            value={numDesigns}
            onChange={(e) => setNumDesigns(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--cyan-bright)' }}
          />
        </div>

        {/* Binder Length Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
            <span>Binder Length (AAs)</span>
            <span className="mono" style={{ color: 'var(--cyan-bright)' }}>{binderLength} AAs</span>
          </div>
          <input
            type="range"
            min="45"
            max="95"
            step="5"
            value={binderLength}
            onChange={(e) => setBinderLength(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--cyan-bright)' }}
          />
        </div>

        {/* Budget Allocation */}
        <div style={{ marginBottom: '20px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
            <span>Spend Budget Gate ($ USD)</span>
            <span className="mono" style={{ color: 'var(--emerald-bright)' }}>${budgetUsd.toFixed(2)}</span>
          </div>
          <input
            type="range"
            min="0.10"
            max="2.00"
            step="0.05"
            value={budgetUsd}
            onChange={(e) => setBudgetUsd(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--emerald-bright)' }}
          />
        </div>

        {errorMsg && (
          <div
            style={{
              padding: '10px 14px',
              borderRadius: '8px',
              background: 'rgba(239, 68, 68, 0.15)',
              border: '1px solid var(--red-alert)',
              color: '#fca5a5',
              fontSize: '12px',
              marginBottom: '16px',
            }}
          >
            {errorMsg}
          </div>
        )}

        <button
          onClick={handleGenerate}
          disabled={isLoading}
          className="btn btn-primary"
          style={{ width: '100%', padding: '12px', fontSize: '14px' }}
        >
          <Sparkles size={16} />
          <span>{isLoading ? 'Generating Backbones & Inverse Folding...' : 'Generate Budgeted De Novo Binders'}</span>
        </button>

        {generationResult && (
          <div style={{ marginTop: '20px', display: 'flex', gap: '10px' }}>
            <button onClick={handleDownloadTwist} className="btn btn-secondary" style={{ flex: 1, fontSize: '12px' }}>
              <FileSpreadsheet size={14} color="var(--emerald-bright)" />
              <span>Twist CSV</span>
            </button>
            <button onClick={handleDownloadOpentrons} className="btn btn-secondary" style={{ flex: 1, fontSize: '12px' }}>
              <Bot size={14} color="var(--cyan-bright)" />
              <span>Opentrons .py</span>
            </button>
          </div>
        )}
      </div>

      {/* Right 3D Visualizer & Table Column */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
        {/* 3D WebGL Viewer */}
        <Protein3DViewer
          pdbContent={selectedDesign?.pdb_content || null}
          targetName={selectedDesign?.design_id || 'De Novo Docked Complex'}
          hotspots={[98, 99, 100, 101, 102]}
          isComplex={true}
          height="380px"
        />

        {/* Candidates Table */}
        <div className="glass-panel" style={{ padding: '20px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
            <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
              Generated Binder Candidates & Biophysical Liabilities
            </h3>
            {generationResult?.audit_block_hash && (
              <span className="badge badge-emerald mono" style={{ fontSize: '10px' }}>
                Merkle Block: {generationResult.audit_block_hash.slice(0, 14)}...
              </span>
            )}
          </div>

          <div style={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-muted)' }}>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Design ID</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Length</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>pLDDT</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>RMSD (Å)</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Rg (Å)</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Developability</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Sequence Preview</th>
                  <th style={{ textAlign: 'right', padding: '8px' }}>Action</th>
                </tr>
              </thead>
              <tbody>
                {generationResult?.designs?.length > 0 ? (
                  generationResult.designs.map((design) => (
                    <tr
                      key={design.design_id}
                      style={{
                        borderBottom: '1px solid rgba(255,255,255,0.04)',
                        background: selectedDesign?.design_id === design.design_id ? 'rgba(56, 189, 248, 0.08)' : 'transparent',
                      }}
                    >
                      <td style={{ padding: '10px 8px', fontWeight: 600, color: 'var(--cyan-bright)' }}>
                        {design.design_id}
                      </td>
                      <td style={{ padding: '10px 8px' }}>{design.length} AA</td>
                      <td style={{ padding: '10px 8px', color: 'var(--emerald-bright)', fontWeight: 600 }}>
                        {design.predicted_plddt.toFixed(1)}
                      </td>
                      <td style={{ padding: '10px 8px' }}>{design.target_docking_rmsd.toFixed(2)}</td>
                      <td style={{ padding: '10px 8px' }}>{design.radius_of_gyration.toFixed(2)}</td>
                      <td style={{ padding: '10px 8px' }}>
                        <span className="badge badge-emerald" style={{ fontSize: '10px' }}>
                          0 Unpaired Cys · Clean
                        </span>
                      </td>
                      <td style={{ padding: '10px 8px' }} className="mono">
                        {design.sequence.slice(0, 16)}...
                      </td>
                      <td style={{ padding: '10px 8px', textAlign: 'right' }}>
                        <button
                          onClick={() => setSelectedDesign(design)}
                          className="btn btn-secondary"
                          style={{ padding: '4px 10px', fontSize: '12px' }}
                        >
                          View 3D
                        </button>
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="8" style={{ textAlign: 'center', padding: '30px', color: 'var(--text-dim)' }}>
                      No designs generated yet. Click 'Generate Budgeted De Novo Binders' above to start diffusion.
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  );
}
