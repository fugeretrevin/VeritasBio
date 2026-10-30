import React, { useState, useEffect, useRef } from 'react';
import {
  Play,
  Award,
  ArrowRight,
  Sliders,
} from 'lucide-react';
import Protein3DViewer from './Protein3DViewer';

export default function CampaignStudio({
  targetAccession = 'TRPB_TM9D8',
  onTargetChange,
  currentUser,
  onRefreshLedger,
  onRunComplete,
  onNavigateTab,
}) {
  const [budgetUsd, setBudgetUsd] = useState(1.00);
  const [maxIterations, setMaxIterations] = useState(4);
  const [batchSize, setBatchSize] = useState(15);
  const [enablePocketHopping, setEnablePocketHopping] = useState(true);

  const [isRunning, setIsRunning] = useState(false);
  const [terminalLogs, setTerminalLogs] = useState([
    'Ready. Configure parameters on the left and click Run Campaign.',
  ]);

  const [targetPdb, setTargetPdb] = useState(null);
  const [topVariants, setTopVariants] = useState([]);
  const activeStreamRef = useRef(null);
  const terminalEndRef = useRef(null);

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
    setTerminalLogs((prev) => [...prev, '■ Campaign stopped by operator.']);
    if (onRefreshLedger) onRefreshLedger();
  };

  const handleLaunch = async () => {
    if (isRunning) return;
    setIsRunning(true);
    setTerminalLogs([
      `▶ Starting ${targetAccession} (${maxIterations} rounds · $${budgetUsd.toFixed(2)} cap · ${batchSize} variants/batch)`,
    ]);

    const tokenToUse =
      currentUser?.id_token ||
      (currentUser?.uid ? `panoptes_session_${currentUser.uid}` : null);

    let streamUrl;
    try {
      const prepareRes = await fetch('/api/campaign/prepare-stream', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          target_accession: targetAccession,
          total_budget_usd: budgetUsd,
          max_iterations: maxIterations,
          batch_size: batchSize,
          target_model: 'mock_esm',
          enable_pocket_hopping: enablePocketHopping,
          auth_token: tokenToUse,
        }),
      });
      if (!prepareRes.ok) {
        const detail = await prepareRes.text();
        throw new Error(`Could not authorize campaign (${prepareRes.status}): ${detail}`);
      }
      const prepared = await prepareRes.json();
      streamUrl = prepared.stream_url;
    } catch (err) {
      setTerminalLogs((prev) => [...prev, `[Error] ${err.message}`]);
      setIsRunning(false);
      return;
    }

    const es = new EventSource(streamUrl);
    activeStreamRef.current = es;

    es.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        if (data.step === 'complete') {
          setTerminalLogs((prev) => [
            ...prev,
            '[Done] ' + (data.message || 'Optimization campaign completed.'),
          ]);
          if (data.top_variants) {
            setTopVariants(data.top_variants);
            if (onRunComplete && data.top_variants.length > 0) {
              onRunComplete({
                target: targetAccession,
                topScore: data.top_variants[0].computed_score,
              });
            }
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
          const remStr =
            data.remaining_budget_usd !== undefined
              ? ` (Rem: $${data.remaining_budget_usd.toFixed(4)})`
              : '';
          let tag = `[R${data.iteration} · ${data.step}]`;
          if (data.step === 'reflection_complete') tag = `[R${data.iteration} · Reflection]`;
          if (data.step === 'authorization') tag = `[R${data.iteration} · Spend Gate]`;
          if (data.step === 'strategy') tag = `[R${data.iteration} · Llama Prior]`;
          if (data.step === 'selection') tag = `[R${data.iteration} · Surrogate]`;
          if (data.step === 'execution') tag = `[R${data.iteration} · A10G Scoring]`;

          setTerminalLogs((prev) => [...prev, `${tag} ${data.message}${remStr}`]);

          if (data.top_variants && data.top_variants.length > 0) {
            setTopVariants(data.top_variants);
          }
        }
      } catch (err) {
        console.error('Failed to parse SSE event payload:', err);
      }
    };

    es.onerror = () => {
      es.close();
      activeStreamRef.current = null;
      setIsRunning(false);
      if (onRefreshLedger) onRefreshLedger();
    };
  };

  const getHotspots = (acc) => {
    if (acc === 'TRPB_TM9D8') return [183, 184, 227, 228];
    if (acc === 'SCFV_9004G_CN2') return [98, 99, 100, 101, 102];
    if (acc === 'HER2_TRASTUZUMAB') return [99, 100, 101, 102];
    if (acc === 'SPG1_STRSG') return [39, 40, 41, 54];
    return [64, 65, 68, 203];
  };

  return (
    <div style={{ display: 'grid', gridTemplateColumns: '340px 1fr', gap: '20px', alignItems: 'start' }}>
      {/* Left Column: Focused Search Parameters Card */}
      <div className="glass-panel" style={{ padding: '24px' }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '18px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <Sliders size={16} color="var(--text-main)" />
            <h2 style={{ fontSize: '15px', fontWeight: 600, margin: 0, color: 'var(--text-main)' }}>
              Campaign Parameters
            </h2>
          </div>
          <span className="mono" style={{ fontSize: '11px', color: 'var(--text-dim)' }}>
            {targetAccession}
          </span>
        </div>

        {/* Target Benchmark Dropdown */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-muted)', marginBottom: '6px' }}>
            Target Protein
          </label>
          <select
            className="input-field"
            value={targetAccession}
            onChange={(e) => onTargetChange && onTargetChange(e.target.value)}
            style={{ fontSize: '12.5px', padding: '8px 10px' }}
          >
            <option value="TRPB_TM9D8">TrpB Synthase (8EH0 · 160K DMS)</option>
            <option value="SCFV_9004G_CN2">Scorpion Antivenom scFv (2YBR)</option>
            <option value="HER2_TRASTUZUMAB">Trastuzumab CDR-H3 (1N8Z)</option>
            <option value="SPG1_STRSG">Protein G GB1 (2QMT · 150K DMS)</option>
            <option value="P42212">avGFP Fluorophore (1EMA · 54K DMS)</option>
          </select>
        </div>

        {/* Budget Allocation Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Compute Spend Cap</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              ${budgetUsd.toFixed(2)}
            </span>
          </div>
          <input
            type="range"
            min="0.25"
            max="10.00"
            step="0.25"
            value={budgetUsd}
            onChange={(e) => setBudgetUsd(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--text-main)' }}
          />
        </div>

        {/* Max Iterations Slider */}
        <div style={{ marginBottom: '16px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Active Learning Rounds</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              {maxIterations} rounds
            </span>
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
        <div style={{ marginBottom: '18px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '5px' }}>
            <span style={{ color: 'var(--text-muted)' }}>Batch Size per Round</span>
            <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              {batchSize} variants
            </span>
          </div>
          <input
            type="range"
            min="5"
            max="50"
            step="5"
            value={batchSize}
            onChange={(e) => setBatchSize(Number(e.target.value))}
            style={{ width: '100%', accentColor: 'var(--text-main)' }}
          />
        </div>

        {/* Hierarchical Pocket Hopping Checkbox */}
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: '10px',
            padding: '12px',
            borderRadius: '8px',
            background: 'var(--bg-card)',
            marginBottom: '18px',
          }}
        >
          <input
            type="checkbox"
            id="hop_check"
            checked={enablePocketHopping}
            onChange={(e) => setEnablePocketHopping(e.target.checked)}
            style={{ accentColor: '#0d9488', width: '15px', height: '15px', cursor: 'pointer' }}
          />
          <label htmlFor="hop_check" style={{ fontSize: '12px', cursor: 'pointer', lineHeight: 1.35 }}>
            <span style={{ fontWeight: 600, color: 'var(--text-main)', display: 'block' }}>
              Qualify-Then-Optimize Frontier
            </span>
            <span style={{ color: 'var(--text-dim)', fontSize: '11px' }}>
              Probe unseen epistatic residues once scaffold qualifies
            </span>
          </label>
        </div>

        <div style={{ display: 'flex', gap: '8px' }}>
          <button
            onClick={handleLaunch}
            disabled={isRunning}
            className="btn btn-primary"
            style={{ flex: 1, padding: '11px', fontSize: '13.5px' }}
          >
            <Play size={15} />
            <span>{isRunning ? 'Running Loop...' : 'Run Campaign'}</span>
          </button>
          {isRunning && (
            <button
              onClick={handleStop}
              className="btn btn-secondary"
              style={{ padding: '11px 14px', fontSize: '13px', color: '#dc2626' }}
            >
              Stop
            </button>
          )}
        </div>
      </div>

      {/* Right Column: 3D Structure + Clean Execution Log + Discovered Leads */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
        <Protein3DViewer
          pdbContent={targetPdb}
          targetName={targetAccession}
          hotspots={getHotspots(targetAccession)}
          isComplex={false}
          height="320px"
        />

        {/* Clean Execution Stream */}
        <div className="terminal-window">
          <div className="terminal-header" style={{ justifyContent: 'space-between' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <span
                style={{
                  width: '7px',
                  height: '7px',
                  borderRadius: '50%',
                  background: isRunning ? '#0d9488' : '#8c8a84',
                }}
              />
              <span className="mono" style={{ fontSize: '11.5px', fontWeight: 600, color: 'var(--text-muted)' }}>
                {isRunning ? 'LIVE STREAM · ACTIVE LEARNING' : 'EXECUTION LOG'}
              </span>
            </div>
            <button
              type="button"
              onClick={() => setTerminalLogs(['Log cleared.'])}
              style={{
                background: 'none',
                border: 'none',
                fontSize: '11px',
                color: 'var(--text-dim)',
                cursor: 'pointer',
              }}
            >
              Clear
            </button>
          </div>
          <div className="terminal-body" style={{ height: '165px' }}>
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
          <div className="glass-panel" style={{ padding: '18px 20px' }}>
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '12px' }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                <Award size={16} color="var(--emerald-primary)" />
                <h3 style={{ fontSize: '14px', fontWeight: 600, margin: 0, color: 'var(--text-main)' }}>
                  Discovered Lead Variants
                </h3>
              </div>
              {onNavigateTab && (
                <button
                  onClick={() => onNavigateTab('wetlab')}
                  className="btn btn-secondary"
                  style={{ padding: '5px 10px', fontSize: '11.5px' }}
                >
                  <span>Prepare Synthesis Order</span>
                  <ArrowRight size={12} />
                </button>
              )}
            </div>
            <div style={{ overflowX: 'auto' }}>
              <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
                <thead>
                  <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-dim)' }}>
                    <th style={{ textAlign: 'left', padding: '7px 8px' }}>Rank</th>
                    <th style={{ textAlign: 'left', padding: '7px 8px' }}>Sequence Preview</th>
                    <th style={{ textAlign: 'left', padding: '7px 8px' }}>Fitness Score</th>
                    <th style={{ textAlign: 'left', padding: '7px 8px' }}>Authorized Job</th>
                  </tr>
                </thead>
                <tbody>
                  {topVariants.map((varItem, idx) => (
                    <tr key={idx} style={{ borderBottom: '1px solid var(--border-subtle)' }}>
                      <td className="mono" style={{ padding: '7px 8px', fontWeight: 600, color: 'var(--text-main)' }}>
                        #{idx + 1}
                      </td>
                      <td style={{ padding: '7px 8px' }} className="mono">
                        {varItem.variant_sequence?.slice(0, 26)}...
                      </td>
                      <td className="mono tabular" style={{ padding: '7px 8px', color: 'var(--emerald-primary)', fontWeight: 600 }}>
                        +{Number(varItem.computed_score).toFixed(4)}
                      </td>
                      <td style={{ padding: '7px 8px', color: 'var(--text-dim)' }} className="mono">
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
