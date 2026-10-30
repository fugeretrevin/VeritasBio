import React, { useState, useEffect } from 'react';
import { Sparkles, FileSpreadsheet, Bot, ArrowRight } from 'lucide-react';
import Protein3DViewer from './Protein3DViewer';

export default function DeNovoStudio({
  targetAccession = 'SCFV_9004G_CN2',
  onTargetChange,
  currentUser,
  onRefreshLedger,
  onRunComplete,
  onNavigateTab,
}) {
  const [targetScaffold, setTargetScaffold] = useState(
    ['SCFV_9004G_CN2', 'HER2_TRASTUZUMAB', 'P42212'].includes(targetAccession)
      ? targetAccession
      : 'SCFV_9004G_CN2'
  );
  const [numDesigns, setNumDesigns] = useState(3);
  const [binderLength, setBinderLength] = useState(65);
  const [budgetUsd, setBudgetUsd] = useState(0.50);

  const [isLoading, setIsLoading] = useState(false);
  const [errorMsg, setErrorMsg] = useState(null);
  const [generationResult, setGenerationResult] = useState(null);
  const [selectedDesign, setSelectedDesign] = useState(null);

  useEffect(() => {
    if (['SCFV_9004G_CN2', 'HER2_TRASTUZUMAB', 'P42212'].includes(targetAccession)) {
      setTargetScaffold(targetAccession);
    }
  }, [targetAccession]);

  const handleGenerate = async () => {
    setIsLoading(true);
    setErrorMsg(null);
    try {
      const tokenToUse =
        currentUser?.id_token ||
        (currentUser?.uid ? `panoptes_session_${currentUser.uid}` : null);

      const res = await fetch('/api/diffusion/generate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          target_scaffold: targetScaffold,
          num_designs: numDesigns,
          binder_length: binderLength,
          budget_usd: budgetUsd,
          auth_token: tokenToUse,
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
        if (onRunComplete) {
          onRunComplete({
            target: targetScaffold,
            topScore: data.designs[0].interface_contact_pairs || 1.45,
          });
        }
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
    <div style={{ display: 'grid', gridTemplateColumns: '340px 1fr', gap: '20px', alignItems: 'start' }}>
      {/* Left Control Column */}
      <div className="glass-panel" style={{ padding: '22px' }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '18px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <Sparkles size={16} color="var(--text-main)" />
            <h2 style={{ fontSize: '15px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              De Novo Binder Parameters
            </h2>
          </div>
          <span className="mono" style={{ fontSize: '11px', color: 'var(--text-dim)' }}>
            3-Helix
          </span>
        </div>

        {/* Target Epitope */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-muted)', marginBottom: '6px' }}>
            Target Epitope
          </label>
          <select
            className="input-field"
            value={targetScaffold}
            onChange={(e) => {
              setTargetScaffold(e.target.value);
              if (onTargetChange) onTargetChange(e.target.value);
            }}
            style={{ fontSize: '12.5px', padding: '8px 10px' }}
          >
            <option value="SCFV_9004G_CN2">Scorpion Neurotoxin Cn2 (2YBR)</option>
            <option value="HER2_TRASTUZUMAB">HER2 Extracellular Domain (1N8Z)</option>
            <option value="P42212">avGFP β-Barrel Epitope (1EMA)</option>
          </select>
        </div>

        {/* Num Designs Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Backbone Candidates</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              {numDesigns} designs
            </span>
          </div>
          <input
            type="range"
            min="1"
            max="8"
            step="1"
            value={numDesigns}
            onChange={(e) => setNumDesigns(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--text-main)' }}
          />
        </div>

        {/* Binder Length Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Scaffold Length</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              {binderLength} AA
            </span>
          </div>
          <input
            type="range"
            min="45"
            max="95"
            step="5"
            value={binderLength}
            onChange={(e) => setBinderLength(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--text-main)' }}
          />
        </div>

        {/* Budget Allocation */}
        <div style={{ marginBottom: '18px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Spend Gate Cap</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              ${budgetUsd.toFixed(2)}
            </span>
          </div>
          <input
            type="range"
            min="0.10"
            max="2.00"
            step="0.05"
            value={budgetUsd}
            onChange={(e) => setBudgetUsd(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--text-main)' }}
          />
        </div>

        {errorMsg && (
          <div
            style={{
              padding: '10px 12px',
              borderRadius: '8px',
              background: '#fef2f2',
              border: '1px solid #fecaca',
              color: '#dc2626',
              fontSize: '12px',
              marginBottom: '14px',
            }}
          >
            {errorMsg}
          </div>
        )}

        <button
          onClick={handleGenerate}
          disabled={isLoading}
          className="btn btn-primary"
          style={{ width: '100%', padding: '11px', fontSize: '13.5px' }}
        >
          <Sparkles size={15} />
          <span>{isLoading ? 'Generating Binders...' : 'Generate De Novo Binders'}</span>
        </button>

        {generationResult && (
          <div style={{ marginTop: '14px', display: 'flex', gap: '8px' }}>
            <button onClick={handleDownloadTwist} className="btn btn-secondary" style={{ flex: 1, fontSize: '12px', padding: '8px' }}>
              <FileSpreadsheet size={13} color="var(--emerald-primary)" />
              <span>Twist CSV</span>
            </button>
            <button onClick={handleDownloadOpentrons} className="btn btn-secondary" style={{ flex: 1, fontSize: '12px', padding: '8px' }}>
              <Bot size={13} color="var(--text-main)" />
              <span>Opentrons .py</span>
            </button>
          </div>
        )}
      </div>

      {/* Right 3D Visualizer & Table Column */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
        <Protein3DViewer
          pdbContent={selectedDesign?.pdb_content || null}
          targetName={selectedDesign?.design_id || `${targetScaffold} Complex`}
          hotspots={[98, 99, 100, 101, 102]}
          isComplex={true}
          height="320px"
        />

        {/* Candidates Table */}
        <div className="glass-panel" style={{ padding: '18px 20px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Designed Binder Candidates
            </h3>
            {generationResult?.audit_block_hash && (
              <span className="mono" style={{ fontSize: '11px', color: 'var(--text-dim)' }}>
                Block #{generationResult.audit_block_hash.slice(0, 10)}
              </span>
            )}
          </div>

          <div style={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-dim)' }}>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Design ID</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Length</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Contacts</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Min Dist</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Rg</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Sequence</th>
                  <th style={{ textAlign: 'right', padding: '7px 8px' }}>Inspect</th>
                </tr>
              </thead>
              <tbody>
                {generationResult?.designs?.length > 0 ? (
                  generationResult.designs.map((design) => (
                    <tr
                      key={design.design_id}
                      style={{
                        borderBottom: '1px solid var(--border-subtle)',
                        background:
                          selectedDesign?.design_id === design.design_id
                            ? 'var(--bg-elevated)'
                            : 'transparent',
                      }}
                    >
                      <td className="mono" style={{ padding: '8px', fontWeight: 600, color: 'var(--text-main)' }}>
                        {design.design_id}
                      </td>
                      <td style={{ padding: '8px' }}>{design.length} AA</td>
                      <td className="mono tabular" style={{ padding: '8px', color: 'var(--emerald-primary)', fontWeight: 600 }}>
                        {design.interface_contact_pairs ?? 0}
                      </td>
                      <td className="mono tabular" style={{ padding: '8px' }}>
                        {design.min_interface_distance_angstrom != null
                          ? `${design.min_interface_distance_angstrom.toFixed(2)} Å`
                          : '—'}
                      </td>
                      <td className="mono tabular" style={{ padding: '8px' }}>
                        {design.radius_of_gyration.toFixed(2)} Å
                      </td>
                      <td style={{ padding: '8px' }} className="mono">
                        {design.sequence.slice(0, 18)}...
                      </td>
                      <td style={{ padding: '8px', textAlign: 'right' }}>
                        <button
                          onClick={() => setSelectedDesign(design)}
                          className="btn btn-secondary"
                          style={{ padding: '4px 10px', fontSize: '11.5px' }}
                        >
                          View 3D
                        </button>
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="7" style={{ textAlign: 'center', padding: '24px', color: 'var(--text-dim)' }}>
                      Click "Generate De Novo Binders" on the left to sample docked 3-helix bundles.
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
