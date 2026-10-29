import React, { useState, useEffect, useRef } from 'react';
import {
  Play,
  Terminal,
  Award,
  Compass,
  RefreshCw,
  UploadCloud,
  FileSpreadsheet,
  CheckCircle2,
  AlertCircle,
} from 'lucide-react';
import Protein3DViewer from './Protein3DViewer';

export default function CampaignStudio({ onRefreshLedger }) {
  const [targetAccession, setTargetAccession] = useState('TRPB_TM9D8');
  const [budgetUsd, setBudgetUsd] = useState(1.00);
  const [maxIterations, setMaxIterations] = useState(4);
  const [batchSize, setBatchSize] = useState(15);
  const [enablePocketHopping, setEnablePocketHopping] = useState(true);

  const [isRunning, setIsRunning] = useState(false);
  const [terminalLogs, setTerminalLogs] = useState([
    '● PanoptesBio active learning kernel initialized.',
    '● Standing by. Select target benchmark and click Launch.',
  ]);

  const [targetPdb, setTargetPdb] = useState(null);
  const [topVariants, setTopVariants] = useState([]);
  const activeStreamRef = useRef(null);
  const terminalEndRef = useRef(null);

  // Microplate Reader Ingestion State
  const [plateUploading, setPlateUploading] = useState(false);
  const [plateResult, setPlateResult] = useState(null);
  const [plateError, setPlateError] = useState(null);
  const [isDragging, setIsDragging] = useState(false);
  const [normalizeToWt, setNormalizeToWt] = useState(true);
  const fileInputRef = useRef(null);

  useEffect(() => {
    if (terminalEndRef.current) {
      terminalEndRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }, [terminalLogs]);

  useEffect(() => {
    return () => {
      if (activeStreamRef.current) {
        activeStreamRef.current.close();
      }
    };
  }, []);

  // Fetch PDB structure when target changes
  useEffect(() => {
    async function loadStructure() {
      try {
        const res = await fetch(`/api/structure/${targetAccession}`);
        if (res.ok) {
          const data = await res.json();
          setTargetPdb(data.pdb_content);
        }
      } catch (err) {
        console.error('Failed to load target PDB:', err);
      }
    }
    loadStructure();
  }, [targetAccession]);

  const handleStop = () => {
    if (activeStreamRef.current) {
      activeStreamRef.current.close();
      activeStreamRef.current = null;
    }
    setIsRunning(false);
    setTerminalLogs((prev) => [...prev, '■ Campaign stream stopped by operator.']);
    if (onRefreshLedger) onRefreshLedger();
  };

  const handleLaunch = () => {
    if (isRunning) return;
    setIsRunning(true);
    setTerminalLogs((prev) => [
      ...prev,
      `--- Starting Real-Time Stream: ${targetAccession} (Budget: $${budgetUsd.toFixed(2)}, Cap: ${maxIterations} rounds) ---`,
    ]);

    const params = new URLSearchParams({
      target_accession: targetAccession,
      total_budget_usd: budgetUsd.toString(),
      max_iterations: maxIterations.toString(),
      batch_size: batchSize.toString(),
      target_model: 'mock_esm',
      enable_pocket_hopping: enablePocketHopping.toString(),
      auth_token: 'panoptes_demo_token',
    });

    const url = `/api/campaign/stream?${params.toString()}`;
    const es = new EventSource(url);
    activeStreamRef.current = es;

    es.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        if (data.step === 'complete') {
          setTerminalLogs((prev) => [
            ...prev,
            '✓ ' + (data.message || 'Optimization campaign completed successfully.'),
          ]);
          if (data.top_variants) {
            setTopVariants(data.top_variants);
          }
          if (onRefreshLedger) onRefreshLedger();
          es.close();
          activeStreamRef.current = null;
          setIsRunning(false);
        } else if (data.step === 'error') {
          setTerminalLogs((prev) => [...prev, `[ERROR] ${data.message}`]);
          es.close();
          activeStreamRef.current = null;
          setIsRunning(false);
        } else {
          const remStr = data.remaining_budget_usd !== undefined ? ` (Rem: $${data.remaining_budget_usd.toFixed(4)})` : '';
          let tag = `[Round ${data.iteration} · ${data.step}]`;
          if (data.step === 'reflection_complete') tag = `[Llama 3.1 Reflection · R${data.iteration}]`;
          if (data.step === 'authorization') tag = `[Ledger Gate · R${data.iteration}]`;
          if (data.step === 'strategy') tag = `[Strategic Prior · R${data.iteration}]`;
          if (data.step === 'selection') tag = `[Quant DKL · R${data.iteration}]`;
          if (data.step === 'execution') tag = `[A10G Muscle · R${data.iteration}]`;

          setTerminalLogs((prev) => [...prev, `${tag} ${data.message}${remStr}`]);

          if (data.top_variants && data.top_variants.length > 0) {
            setTopVariants(data.top_variants);
          }
        }
      } catch (err) {
        console.error('Failed to parse SSE event payload:', err);
      }
    };

    es.onerror = (err) => {
      console.warn('SSE EventSource disconnected or finished:', err);
      es.close();
      activeStreamRef.current = null;
      setIsRunning(false);
      if (onRefreshLedger) onRefreshLedger();
    };
  };

  // Microplate Ingestion Handlers
  const handlePlateIngestion = async (fileOrText, isText = false) => {
    setPlateUploading(true);
    setPlateError(null);
    setPlateResult(null);

    try {
      let res;
      if (isText) {
        res = await fetch('/api/campaign/ingest-plate', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            csv_content: fileOrText,
            target_accession: targetAccession,
            normalize_to_wt: normalizeToWt,
            auth_token: 'panoptes_demo_token',
          }),
        });
      } else {
        const formData = new FormData();
        formData.append('file', fileOrText);
        formData.append('target_accession', targetAccession);
        formData.append('normalize_to_wt', normalizeToWt ? 'true' : 'false');
        formData.append('auth_token', 'panoptes_demo_token');

        res = await fetch('/api/campaign/ingest-plate', {
          method: 'POST',
          body: formData,
        });
      }

      if (!res.ok) {
        const errData = await res.json().catch(() => ({ detail: 'Failed to ingest microplate CSV' }));
        throw new Error(errData.detail || 'Microplate ingestion failed');
      }

      const data = await res.json();
      setPlateResult(data);

      setTerminalLogs((prev) => [
        ...prev,
        `✓ [WET-LAB INGESTION] Microplate processed: ${data.variants_ingested} variants (${data.barrier_wells_ignored} barrier wells verified).`,
        `  ↳ Champion: ${data.top_variant} (Fitness: +${Number(data.top_fitness_discovered).toFixed(4)}) · Merkle Block #${data.block_index ?? 'N/A'}`,
      ]);

      if (data.top_variants_ranked && data.top_variants_ranked.length > 0) {
        setTopVariants(
          data.top_variants_ranked.map((v) => ({
            variant_sequence: v.variant_sequence || v.mutations,
            computed_score: v.normalized_score,
            job_id: data.job_id,
          }))
        );
      }

      if (onRefreshLedger) onRefreshLedger();
    } catch (err) {
      setPlateError(err.message);
      setTerminalLogs((prev) => [...prev, `[ERROR] Plate Ingestion: ${err.message}`]);
    } finally {
      setPlateUploading(false);
    }
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setIsDragging(false);
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      handlePlateIngestion(e.dataTransfer.files[0], false);
    }
  };

  const handleFileChange = (e) => {
    if (e.target.files && e.target.files.length > 0) {
      handlePlateIngestion(e.target.files[0], false);
    }
  };

  const handleLoadSamplePlate = () => {
    let sampleCsv = '';
    if (targetAccession === 'TRPB_TM9D8') {
      sampleCsv = [
        'Well,Sample_Type,Sample_Name,Reading',
        'A01,evaporation_barrier,H2O,0.01',
        'A02,evaporation_barrier,H2O,0.01',
        'A12,evaporation_barrier,H2O,0.01',
        'B02,control_wt,WT,1.00',
        'B03,variant,V183A:F184A:V227A:S228A,0.74',
        'B04,variant,F184V:S228A,1.82',
        'B05,variant,V183I:F184V:V227I:S228A,2.15',
        'H12,evaporation_barrier,H2O,0.01',
      ].join('\n');
    } else {
      sampleCsv = [
        'Well,Sample_Type,Sample_Name,Reading',
        'A01,evaporation_barrier,H2O,0.01',
        'A02,evaporation_barrier,H2O,0.01',
        'A12,evaporation_barrier,H2O,0.01',
        'B02,control_wt,WT,1.00',
        'B03,variant,S65T:V68L,1.85',
        'B04,variant,T203Y,1.40',
        'B05,variant,F64L:S65T:T203Y,2.45',
        'H12,evaporation_barrier,H2O,0.01',
      ].join('\n');
    }
    handlePlateIngestion(sampleCsv, true);
  };

  const getHotspots = (acc) => {
    if (acc === 'TRPB_TM9D8') return [183, 184, 227, 228];
    if (acc === 'SCFV_9004G_CN2') return [98, 99, 100, 101, 102];
    if (acc === 'HER2_TRASTUZUMAB') return [99, 100, 101, 102];
    if (acc === 'SPG1_STRSG') return [39, 40, 41, 54];
    return [64, 65, 68, 203];
  };

  return (
    <div style={{ display: 'grid', gridTemplateColumns: 'minmax(340px, 440px) 1fr', gap: '20px' }}>
      {/* Left Column: Search Controls & Wet-Lab Plate Ingestion */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
        {/* Card 1: Autonomous Active Learning Controls */}
        <div className="glass-panel" style={{ padding: '24px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '16px' }}>
            <Compass size={20} color="var(--cyan-bright)" />
            <h2 className="title-display" style={{ fontSize: '18px', margin: 0 }}>
              Autonomous Search Controls
            </h2>
          </div>
          <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '20px' }}>
            Biophysical Active Learning loop governed by Meta Llama 3.1 hypothesis planning, Deep Kernel Learning,
            and fail-closed gatekeeping.
          </p>

          {/* Target Benchmark Dropdown */}
          <div style={{ marginBottom: '16px' }}>
            <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, marginBottom: '6px' }}>
              Target Protein Benchmark
            </label>
            <select
              className="input-field"
              value={targetAccession}
              onChange={(e) => setTargetAccession(e.target.value)}
            >
              <option value="TRPB_TM9D8">TRPB_TM9D8 (TrpB Active Site · 160K DMS · 8EH0)</option>
              <option value="SCFV_9004G_CN2">SCFV_9004G_CN2 (Scorpion Antivenom scFv · 2YBR)</option>
              <option value="HER2_TRASTUZUMAB">HER2_TRASTUZUMAB (Trastuzumab CDR-H3 · 1N8Z)</option>
              <option value="SPG1_STRSG">SPG1_STRSG (Protein G GB1 · 150K DMS · 2QMT)</option>
              <option value="P42212">P42212 (avGFP Fluorescent Protein · 54K DMS · 1EMA)</option>
            </select>
          </div>

          {/* Budget Allocation Slider */}
          <div style={{ marginBottom: '16px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
              <span>Compute Budget Allocation ($ USD)</span>
              <span className="mono" style={{ color: 'var(--cyan-bright)' }}>${budgetUsd.toFixed(2)}</span>
            </div>
            <input
              type="range"
              min="0.25"
              max="10.00"
              step="0.25"
              value={budgetUsd}
              onChange={(e) => setBudgetUsd(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--cyan-bright)' }}
            />
          </div>

          {/* Max Iterations Slider */}
          <div style={{ marginBottom: '16px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
              <span>Maximum Active Learning Iterations</span>
              <span className="mono" style={{ color: 'var(--text-main)' }}>{maxIterations} rounds</span>
            </div>
            <input
              type="range"
              min="1"
              max="10"
              step="1"
              value={maxIterations}
              onChange={(e) => setMaxIterations(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--text-main)' }}
            />
          </div>

          {/* Micro-Batch Size */}
          <div style={{ marginBottom: '20px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '6px' }}>
              <span>Micro-Batch Size</span>
              <span className="mono" style={{ color: 'var(--emerald-bright)' }}>{batchSize} variants/round</span>
            </div>
            <input
              type="range"
              min="5"
              max="50"
              step="5"
              value={batchSize}
              onChange={(e) => setBatchSize(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--emerald-bright)' }}
            />
          </div>

          {/* Hierarchical Pocket Hopping Checkbox */}
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '10px',
              padding: '10px 14px',
              borderRadius: '8px',
              background: 'rgba(15, 23, 42, 0.6)',
              border: '1px solid var(--border-subtle)',
              marginBottom: '20px',
            }}
          >
            <input
              type="checkbox"
              id="hop_check"
              checked={enablePocketHopping}
              onChange={(e) => setEnablePocketHopping(e.target.checked)}
              style={{ accentColor: 'var(--cyan-bright)', width: '16px', height: '16px' }}
            />
            <label htmlFor="hop_check" style={{ fontSize: '12px', cursor: 'pointer' }}>
              <span style={{ fontWeight: 600, color: 'var(--text-main)', display: 'block' }}>
                Hierarchical Pocket-Hopping
              </span>
              <span style={{ color: 'var(--text-dim)', fontSize: '12px' }}>
                Lock lead mutations on regret saturation & sweep coupled secondary loop
              </span>
            </label>
          </div>

          <div style={{ display: 'flex', gap: '8px' }}>
            <button
              onClick={handleLaunch}
              disabled={isRunning}
              className="btn btn-primary"
              style={{ flex: 1, padding: '12px', fontSize: '14px' }}
            >
              <Play size={16} />
              <span>{isRunning ? 'Streaming Optimization...' : 'Launch Active Learning Loop'}</span>
            </button>
            {isRunning && (
              <button
                onClick={handleStop}
                className="btn btn-secondary"
                style={{ padding: '12px 16px', fontSize: '14px', borderColor: '#ef4444', color: '#ef4444' }}
              >
                Stop
              </button>
            )}
          </div>
        </div>

        {/* Card 2: Microplate Reader Data Ingestion (Wet-Lab Handoff) */}
        <div className="glass-panel" style={{ padding: '24px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '12px' }}>
            <FileSpreadsheet size={20} color="var(--emerald-bright)" />
            <h2 className="title-display" style={{ fontSize: '17px', margin: 0 }}>
              Physical Microplate Ingestion
            </h2>
          </div>
          <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '16px' }}>
            Ingest 96/384-well microplate reader experimental data (Tecan Spark, BioTek Synergy).
            Automatically verifies anti-evaporation hydration barriers and writes to the Merkle audit ledger.
          </p>

          {/* Drag and Drop Zone */}
          <div
            onDragOver={(e) => {
              e.preventDefault();
              setIsDragging(true);
            }}
            onDragLeave={() => setIsDragging(false)}
            onDrop={handleDrop}
            onClick={() => fileInputRef.current?.click()}
            style={{
              border: `2px dashed ${isDragging ? 'var(--cyan-bright)' : 'rgba(255,255,255,0.15)'}`,
              borderRadius: '10px',
              padding: '24px 16px',
              textAlign: 'center',
              cursor: 'pointer',
              background: isDragging ? 'rgba(56, 189, 248, 0.08)' : 'rgba(15, 23, 42, 0.4)',
              transition: 'all 0.2s ease',
              marginBottom: '16px',
            }}
          >
            <input
              ref={fileInputRef}
              type="file"
              accept=".csv,.txt"
              onChange={handleFileChange}
              style={{ display: 'none' }}
            />
            <UploadCloud size={32} color={isDragging ? 'var(--cyan-bright)' : 'var(--text-dim)'} style={{ margin: '0 auto 8px' }} />
            <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '4px' }}>
              {plateUploading ? 'Ingesting microplate data...' : 'Drop plate reader CSV here'}
            </div>
            <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>
              Supports Tecan, BioTek, and VeritasBio robotic CSV exports
            </div>
          </div>

          {/* Controls & Quick Actions */}
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
            <label style={{ display: 'flex', alignItems: 'center', gap: '8px', fontSize: '12px', cursor: 'pointer' }}>
              <input
                type="checkbox"
                checked={normalizeToWt}
                onChange={(e) => setNormalizeToWt(e.target.checked)}
                style={{ accentColor: 'var(--emerald-bright)' }}
              />
              <span style={{ color: 'var(--text-muted)' }}>Normalize to WT Parent (Well B2)</span>
            </label>

            <button
              onClick={handleLoadSamplePlate}
              disabled={plateUploading}
              className="btn btn-secondary"
              style={{ padding: '6px 12px', fontSize: '11px', gap: '6px' }}
            >
              <span>⚡ Load Sample 96-Well Plate</span>
            </button>
          </div>

          {/* Plate Error Banner */}
          {plateError && (
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: '8px',
                padding: '10px 14px',
                borderRadius: '8px',
                background: 'rgba(239, 68, 68, 0.1)',
                border: '1px solid rgba(239, 68, 68, 0.3)',
                color: '#f87171',
                fontSize: '12px',
                marginBottom: '14px',
              }}
            >
              <AlertCircle size={16} />
              <span>{plateError}</span>
            </div>
          )}

          {/* Plate Ingestion Success Panel */}
          {plateResult && (
            <div
              style={{
                padding: '14px',
                borderRadius: '8px',
                background: 'rgba(16, 185, 129, 0.08)',
                border: '1px solid rgba(16, 185, 129, 0.25)',
                fontSize: '12px',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px', color: 'var(--emerald-bright)', fontWeight: 600, marginBottom: '8px' }}>
                <CheckCircle2 size={16} />
                <span>Microplate Ingested · Ledger Updated</span>
              </div>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '8px', color: 'var(--text-muted)' }}>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Variants Recorded</span>
                  <strong style={{ color: 'var(--text-main)' }}>{plateResult.variants_ingested}</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Perimeter Barriers Verified</span>
                  <strong style={{ color: 'var(--text-main)' }}>{plateResult.barrier_wells_ignored} wells</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>WT Parent Baseline</span>
                  <strong style={{ color: 'var(--text-main)' }}>{plateResult.wt_baseline_score ?? 1.00}</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Top Hit Discovered</span>
                  <strong style={{ color: 'var(--emerald-bright)' }}>{plateResult.top_variant} (+{Number(plateResult.top_fitness_discovered).toFixed(4)})</strong>
                </div>
              </div>
              {plateResult.block_hash && (
                <div style={{ marginTop: '10px', paddingTop: '8px', borderTop: '1px solid rgba(255,255,255,0.06)', fontSize: '11px', color: 'var(--text-dim)' }}>
                  Merkle Block #{plateResult.block_index}: <span className="mono" style={{ color: 'var(--cyan-bright)' }}>{plateResult.block_hash.slice(0, 18)}...</span>
                </div>
              )}
            </div>
          )}
        </div>
      </div>

      {/* Right Column: 3D Visualization + Real-Time Terminal Log */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
        {/* 3D WebGL Protein Viewer */}
        <Protein3DViewer
          pdbContent={targetPdb}
          targetName={targetAccession}
          hotspots={getHotspots(targetAccession)}
          isComplex={false}
          height="340px"
        />

        {/* Terminal Live Stream */}
        <div className="terminal-window">
          <div className="terminal-header">
            <span className="terminal-circle" style={{ background: '#ef4444' }} />
            <span className="terminal-circle" style={{ background: '#f59e0b' }} />
            <span className="terminal-circle" style={{ background: '#10b981' }} />
            <span style={{ fontSize: '12px', color: '#94a3b8', marginLeft: '6px' }}>
              panoptesbio-kernel :: live-event-stream.log
            </span>
          </div>
          <div className="terminal-body">
            {terminalLogs.map((log, idx) => (
              <div key={idx} style={{ marginBottom: '3px' }}>
                {log}
              </div>
            ))}
            <div ref={terminalEndRef} />
          </div>
        </div>

        {/* Top Discovered Variants Table */}
        {topVariants.length > 0 && (
          <div className="glass-panel" style={{ padding: '20px' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '12px' }}>
              <Award size={18} color="var(--emerald-bright)" />
              <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
                Top Discovered Epistatic Variants
              </h3>
            </div>
            <div style={{ overflowX: 'auto' }}>
              <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
                <thead>
                  <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-muted)' }}>
                    <th style={{ textAlign: 'left', padding: '8px' }}>Rank</th>
                    <th style={{ textAlign: 'left', padding: '8px' }}>Mutations / Sequence</th>
                    <th style={{ textAlign: 'left', padding: '8px' }}>Fitness Score</th>
                    <th style={{ textAlign: 'left', padding: '8px' }}>Job ID</th>
                  </tr>
                </thead>
                <tbody>
                  {topVariants.map((varItem, idx) => (
                    <tr key={idx} style={{ borderBottom: '1px solid rgba(255,255,255,0.04)' }}>
                      <td style={{ padding: '8px', fontWeight: 600, color: 'var(--cyan-bright)' }}>
                        #{idx + 1}
                      </td>
                      <td style={{ padding: '8px' }} className="mono">
                        {varItem.variant_sequence?.slice(0, 24)}...
                      </td>
                      <td style={{ padding: '8px', color: 'var(--emerald-bright)', fontWeight: 600 }}>
                        +{Number(varItem.computed_score).toFixed(4)}
                      </td>
                      <td style={{ padding: '8px', color: 'var(--text-dim)' }} className="mono">
                        {varItem.job_id}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
