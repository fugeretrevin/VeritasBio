import React, { useState, useEffect, useMemo } from 'react';
import { Sliders, Award, AlertTriangle, ShieldCheck, CheckCircle2, RefreshCw, Download, ArrowUpRight, Sparkles } from 'lucide-react';

export default function ParetoFrontierStudio({ onRefreshLedger }) {
  const [selectedTarget, setSelectedTarget] = useState('SCFV_9004G_CN2');
  const [loading, setLoading] = useState(false);
  const [data, setData] = useState(null);
  const [selectedPoint, setSelectedPoint] = useState(null);

  // 4 Objective Weights (0 to 100)
  const [wFitness, setWFitness] = useState(40);
  const [wDevelopability, setWDevelopability] = useState(30);
  const [wSynthesizability, setWSynthesizability] = useState(15);
  const [wStability, setWStability] = useState(15);

  const targets = [
    { id: 'SCFV_9004G_CN2', label: 'Scorpion Antivenom scFv (2YBR)', desc: 'Paratope optimization & disulfide preservation' },
    { id: 'HER2_TRASTUZUMAB', label: 'Trastuzumab : HER2 (1N8Z)', desc: 'CDR-H3 binding vs CMC viscosity trade-offs' },
    { id: 'TRPB_TM9D8', label: 'Tryptophan Synthase (8EH0)', desc: 'Catalytic turnover vs PLP core clashes' },
    { id: 'P42212', label: 'avGFP Chromophore (1EMA)', desc: 'Fluorophore quantum yield vs folding stability' },
    { id: 'SPG1_STRSG', label: 'Protein G (GB1) (2QMT)', desc: 'Core packing vs epistatic fitness' },
  ];

  const fetchParetoFrontier = async () => {
    setLoading(true);
    try {
      const payload = {
        target_protein: selectedTarget,
        weights: {
          fitness: wFitness / 100,
          developability: wDevelopability / 100,
          synthesizability: wSynthesizability / 100,
          stability: wStability / 100,
        },
      };

      const res = await fetch('/api/campaign/pareto-frontier', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });

      if (res.ok) {
        const json = await res.json();
        setData(json);
        if (json.frontier_points && json.frontier_points.length > 0) {
          setSelectedPoint(json.frontier_points[0]);
        }
      } else {
        console.error('Failed to load Pareto frontier:', await res.text());
      }
    } catch (err) {
      console.error('Network error fetching Pareto frontier:', err);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchParetoFrontier();
  }, [selectedTarget]);

  // Apply Presets
  const applyPreset = (presetName) => {
    if (presetName === 'therapeutic') {
      setWFitness(30);
      setWDevelopability(50);
      setWSynthesizability(10);
      setWStability(10);
    } else if (presetName === 'affinity') {
      setWFitness(70);
      setWDevelopability(15);
      setWSynthesizability(10);
      setWStability(5);
    } else if (presetName === 'balanced') {
      setWFitness(40);
      setWDevelopability(30);
      setWSynthesizability(15);
      setWStability(15);
    } else if (presetName === 'stability') {
      setWFitness(25);
      setWDevelopability(15);
      setWSynthesizability(10);
      setWStability(50);
    }
  };

  // Re-fetch on slider release
  const handleRecalculate = () => {
    fetchParetoFrontier();
  };

  // Scatter chart coordinates calculation
  const chartPoints = useMemo(() => {
    if (!data?.chart_data || data.chart_data.length === 0) return [];
    const xs = data.chart_data.map(p => p.x_fitness);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const spanX = maxX - minX > 0.001 ? maxX - minX : 1.0;

    return data.chart_data.map(p => {
      // Map to 0-100 range for SVG viewbox (500 wide, 300 high with 40px margins)
      const normX = (p.x_fitness - minX) / spanX;
      const normY = p.y_developability; // naturally 0 to 1

      const svgX = 50 + normX * 410;
      const svgY = 250 - normY * 200; // Inverted Y for SVG

      return {
        ...p,
        svgX,
        svgY,
      };
    });
  }, [data]);

  // Frontier points sorted by X for drawing the trade-off curve
  const frontierLinePoints = useMemo(() => {
    const front = chartPoints.filter(p => p.is_frontier);
    return front.sort((a, b) => a.svgX - b.svgX);
  }, [chartPoints]);

  const frontierPathD = useMemo(() => {
    if (frontierLinePoints.length < 2) return '';
    return frontierLinePoints.reduce((acc, curr, idx) => {
      return idx === 0 ? `M ${curr.svgX} ${curr.svgY}` : `${acc} L ${curr.svgX} ${curr.svgY}`;
    }, '');
  }, [frontierLinePoints]);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
      {/* Header Panel */}
      <div className="glass-panel" style={{ padding: '20px 24px' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '16px' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
              <h2 style={{ fontSize: '20px', fontWeight: 700, margin: 0, color: 'var(--text-primary)' }}>
                Multi-Objective Pareto Frontier Studio
              </h2>
              <span className="badge badge-success" style={{ display: 'inline-flex', alignItems: 'center', gap: '4px' }}>
                <Sparkles size={12} />
                <span>NSGA-II Kung/Deb Engine</span>
              </span>
            </div>
            <p style={{ margin: '6px 0 0', fontSize: '13px', color: 'var(--text-secondary)' }}>
              Trade-off optimization between Predicted Fitness, CMC Biologics Developability, Synthesizability, and Stability.
            </p>
          </div>

          <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
            <select
              value={selectedTarget}
              onChange={(e) => setSelectedTarget(e.target.value)}
              className="select"
              style={{ minWidth: '240px', padding: '8px 12px', fontSize: '13px' }}
            >
              {targets.map(t => (
                <option key={t.id} value={t.id}>{t.label}</option>
              ))}
            </select>

            <button
              onClick={fetchParetoFrontier}
              disabled={loading}
              className="btn btn-secondary"
              style={{ padding: '8px 14px' }}
            >
              <RefreshCw size={14} className={loading ? 'animate-spin' : ''} />
              <span>Refresh</span>
            </button>
          </div>
        </div>
      </div>

      {/* KPI Stats Row */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))', gap: '16px' }}>
        <div className="glass-panel" style={{ padding: '16px 20px' }}>
          <div style={{ fontSize: '12px', color: 'var(--text-secondary)', fontWeight: 600 }}>Frontier Leads (Rank 0)</div>
          <div style={{ fontSize: '26px', fontWeight: 700, color: '#f59e0b', marginTop: '4px' }}>
            {data?.rank_0_count ?? 0}
            <span style={{ fontSize: '13px', color: 'var(--text-secondary)', fontWeight: 400, marginLeft: '6px' }}>
              / {data?.points?.length ?? 0} candidates
            </span>
          </div>
          <div style={{ fontSize: '11px', color: '#10b981', marginTop: '4px' }}>
            Non-dominated optimal trade-offs
          </div>
        </div>

        <div className="glass-panel" style={{ padding: '16px 20px' }}>
          <div style={{ fontSize: '12px', color: 'var(--text-secondary)', fontWeight: 600 }}>2D Hypervolume Indicator</div>
          <div style={{ fontSize: '26px', fontWeight: 700, color: '#38bdf8', marginTop: '4px' }}>
            {data?.hypervolume ? data.hypervolume.toFixed(4) : '0.0000'}
          </div>
          <div style={{ fontSize: '11px', color: 'var(--text-secondary)', marginTop: '4px' }}>
            Frontier coverage metric relative to nadir
          </div>
        </div>

        <div className="glass-panel" style={{ padding: '16px 20px' }}>
          <div style={{ fontSize: '12px', color: 'var(--text-secondary)', fontWeight: 600 }}>Dominated Rejected Leads</div>
          <div style={{ fontSize: '26px', fontWeight: 700, color: '#ec4899', marginTop: '4px' }}>
            {data?.dominated_count ?? 0}
          </div>
          <div style={{ fontSize: '11px', color: 'var(--text-secondary)', marginTop: '4px' }}>
            Filtered out by superior alternatives
          </div>
        </div>

        <div className="glass-panel" style={{ padding: '16px 20px' }}>
          <div style={{ fontSize: '12px', color: 'var(--text-secondary)', fontWeight: 600 }}>Top Weighted Utility</div>
          <div style={{ fontSize: '26px', fontWeight: 700, color: '#10b981', marginTop: '4px' }}>
            {data?.points?.[0] ? data.points[0].weighted_utility.toFixed(3) : '1.000'}
          </div>
          <div style={{ fontSize: '11px', color: 'var(--text-secondary)', marginTop: '4px' }}>
            Balanced composite lead score
          </div>
        </div>
      </div>

      {/* Main Grid: Controls & 2D Frontier Scatter Plot */}
      <div style={{ display: 'grid', gridTemplateColumns: '320px 1fr', gap: '20px' }}>
        {/* Left: Dynamic Multi-Objective Weight Sliders */}
        <div className="glass-panel" style={{ padding: '20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
            <Sliders size={16} color="#38bdf8" />
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0 }}>Objective Weighting</h3>
          </div>

          {/* Quick Presets */}
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '6px', marginBottom: '18px' }}>
            <button onClick={() => applyPreset('therapeutic')} className="btn btn-secondary" style={{ fontSize: '11px', padding: '6px' }}>
              Therapeutic CMC
            </button>
            <button onClick={() => applyPreset('affinity')} className="btn btn-secondary" style={{ fontSize: '11px', padding: '6px' }}>
              Max Affinity
            </button>
            <button onClick={() => applyPreset('balanced')} className="btn btn-secondary" style={{ fontSize: '11px', padding: '6px' }}>
              Balanced 40/30
            </button>
            <button onClick={() => applyPreset('stability')} className="btn btn-secondary" style={{ fontSize: '11px', padding: '6px' }}>
              High Stability
            </button>
          </div>

          <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
                <span style={{ color: 'var(--text-primary)', fontWeight: 600 }}>Binding Affinity / Fitness</span>
                <span style={{ color: '#38bdf8', fontWeight: 700 }}>{wFitness}%</span>
              </div>
              <input
                type="range"
                min="5"
                max="90"
                value={wFitness}
                onChange={(e) => setWFitness(Number(e.target.value))}
                style={{ width: '100%', accentColor: '#38bdf8' }}
              />
            </div>

            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
                <span style={{ color: 'var(--text-primary)', fontWeight: 600 }}>Biologics Developability</span>
                <span style={{ color: '#10b981', fontWeight: 700 }}>{wDevelopability}%</span>
              </div>
              <input
                type="range"
                min="5"
                max="90"
                value={wDevelopability}
                onChange={(e) => setWDevelopability(Number(e.target.value))}
                style={{ width: '100%', accentColor: '#10b981' }}
              />
              <div style={{ fontSize: '10px', color: 'var(--text-secondary)', marginTop: '2px' }}>
                Penalizes NG/NS deamidation, DG isomerization, odd cysteines
              </div>
            </div>

            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
                <span style={{ color: 'var(--text-primary)', fontWeight: 600 }}>Synthesizability & Expression</span>
                <span style={{ color: '#f59e0b', fontWeight: 700 }}>{wSynthesizability}%</span>
              </div>
              <input
                type="range"
                min="5"
                max="90"
                value={wSynthesizability}
                onChange={(e) => setWSynthesizability(Number(e.target.value))}
                style={{ width: '100%', accentColor: '#f59e0b' }}
              />
            </div>

            <div>
              <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
                <span style={{ color: 'var(--text-primary)', fontWeight: 600 }}>Conformational Stability</span>
                <span style={{ color: '#a855f7', fontWeight: 700 }}>{wStability}%</span>
              </div>
              <input
                type="range"
                min="5"
                max="90"
                value={wStability}
                onChange={(e) => setWStability(Number(e.target.value))}
                style={{ width: '100%', accentColor: '#a855f7' }}
              />
            </div>

            <button
              onClick={handleRecalculate}
              className="btn btn-primary"
              style={{ marginTop: '10px', width: '100%', justifyContent: 'center' }}
            >
              <RefreshCw size={14} />
              <span>Apply Weights & Re-rank</span>
            </button>
          </div>
        </div>

        {/* Right: 2D Pareto Scatter Trade-off Viewport */}
        <div className="glass-panel" style={{ padding: '20px', display: 'flex', flexDirection: 'column' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0 }}>
              Pareto Frontier: Predicted Fitness vs. CMC Developability
            </h3>
            <div style={{ display: 'flex', gap: '14px', fontSize: '11px', color: 'var(--text-secondary)' }}>
              <span style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                <span style={{ width: '10px', height: '10px', borderRadius: '50%', background: '#f59e0b', boxShadow: '0 0 6px #f59e0b' }}></span>
                Rank 0 (Optimal Frontier)
              </span>
              <span style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                <span style={{ width: '8px', height: '8px', borderRadius: '50%', background: '#64748b' }}></span>
                Dominated / Sub-optimal
              </span>
            </div>
          </div>

          {/* Interactive SVG Chart */}
          <div style={{ position: 'relative', width: '100%', height: '320px', background: 'rgba(15, 23, 42, 0.4)', borderRadius: '8px', border: '1px solid var(--border-color)', overflow: 'hidden' }}>
            <svg viewBox="0 0 500 300" style={{ width: '100%', height: '100%' }}>
              {/* Grid Lines */}
              <line x1="50" y1="50" x2="470" y2="50" stroke="rgba(255,255,255,0.06)" strokeDasharray="4 4" />
              <line x1="50" y1="150" x2="470" y2="150" stroke="rgba(255,255,255,0.06)" strokeDasharray="4 4" />
              <line x1="50" y1="250" x2="470" y2="250" stroke="rgba(255,255,255,0.15)" />
              <line x1="50" y1="50" x2="50" y2="250" stroke="rgba(255,255,255,0.15)" />

              {/* Axis Labels */}
              <text x="260" y="285" fill="#94a3b8" fontSize="11" textAnchor="middle">Predicted Fitness / Binding Affinity →</text>
              <text x="20" y="150" fill="#94a3b8" fontSize="11" textAnchor="middle" transform="rotate(-90 20 150)">Developability Score (0 to 1.0) →</text>

              {/* Trade-off Frontier Curve Line */}
              {frontierPathD && (
                <path
                  d={frontierPathD}
                  fill="none"
                  stroke="#f59e0b"
                  strokeWidth="2"
                  strokeDasharray="3 3"
                  opacity="0.8"
                />
              )}

              {/* Points */}
              {chartPoints.map((pt, i) => {
                const isSelected = selectedPoint?.notation === pt.label;
                const isFrontier = pt.is_frontier;
                const fillColor = isFrontier ? '#f59e0b' : '#64748b';
                const r = isSelected ? 8 : (isFrontier ? 6 : 4);

                return (
                  <g key={i} onClick={() => setSelectedPoint(data?.points?.find(p => p.notation === pt.label))} style={{ cursor: 'pointer' }}>
                    {isFrontier && (
                      <circle cx={pt.svgX} cy={pt.svgY} r={r + 3} fill="rgba(245, 158, 11, 0.25)" />
                    )}
                    <circle
                      cx={pt.svgX}
                      cy={pt.svgY}
                      r={r}
                      fill={fillColor}
                      stroke={isSelected ? '#ffffff' : (isFrontier ? '#fef3c7' : '#334155')}
                      strokeWidth={isSelected ? 2 : 1}
                    />
                    {isFrontier && (
                      <text x={pt.svgX} y={pt.svgY - 10} fill="#fef3c7" fontSize="9" textAnchor="middle" fontWeight="bold">
                        {pt.label}
                      </text>
                    )}
                  </g>
                );
              })}
            </svg>
          </div>

          {/* Selected Point Inspector */}
          {selectedPoint && (
            <div style={{ marginTop: '14px', padding: '12px 16px', background: 'rgba(30, 41, 59, 0.5)', borderRadius: '8px', border: '1px solid var(--border-color)', display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '12px' }}>
              <div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                  <span style={{ fontSize: '15px', fontWeight: 700, color: 'var(--text-primary)' }}>{selectedPoint.notation}</span>
                  <span className={`badge ${selectedPoint.rank === 0 ? 'badge-warning' : 'badge-secondary'}`}>
                    {selectedPoint.rank === 0 ? 'Rank 0 (Optimal)' : `Rank ${selectedPoint.rank}`}
                  </span>
                  <span style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>
                    Weighted Utility: <strong style={{ color: '#10b981' }}>{selectedPoint.weighted_utility.toFixed(3)}</strong>
                  </span>
                </div>
                <div style={{ display: 'flex', gap: '16px', fontSize: '11px', color: 'var(--text-secondary)', marginTop: '4px' }}>
                  <span>Fitness: <strong>{selectedPoint.fitness.toFixed(3)}</strong></span>
                  <span>Developability: <strong style={{ color: selectedPoint.developability_score > 0.85 ? '#10b981' : '#f59e0b' }}>{(selectedPoint.developability_score * 100).toFixed(0)}%</strong></span>
                  <span>Synthesizability: <strong>{(selectedPoint.synthesizability_score * 100).toFixed(0)}%</strong></span>
                  <span>Stability: <strong>{selectedPoint.stability_score.toFixed(2)}</strong></span>
                </div>
              </div>

              <div>
                {selectedPoint.liabilities?.clinical_liabilities?.length > 0 ? (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '6px', color: '#f59e0b', fontSize: '11px' }}>
                    <AlertTriangle size={14} />
                    <span>{selectedPoint.liabilities.clinical_liabilities.join(' • ')}</span>
                  </div>
                ) : (
                  <div style={{ display: 'flex', alignItems: 'center', gap: '6px', color: '#10b981', fontSize: '11px' }}>
                    <CheckCircle2 size={14} />
                    <span>Zero Clinical Liabilities (Formulation Ready)</span>
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Rank 0 Non-Dominated Leads Table */}
      <div className="glass-panel" style={{ padding: '20px' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '14px' }}>
          <div>
            <h3 style={{ fontSize: '15px', fontWeight: 700, margin: 0 }}>
              Rank 0 Non-Dominated Leads (Ready for Physical Synthesis)
            </h3>
            <p style={{ margin: '4px 0 0', fontSize: '12px', color: 'var(--text-secondary)' }}>
              These variants represent the optimal Pareto surface. No variant exists with both higher affinity and better developability.
            </p>
          </div>

          <button
            onClick={() => {
              if (onRefreshLedger) onRefreshLedger();
              alert(`Added ${data?.frontier_points?.length ?? 0} Pareto-optimal leads to synthesis queue!`);
            }}
            className="btn btn-primary"
            style={{ padding: '8px 16px' }}
          >
            <Download size={14} />
            <span>Export Rank 0 to Twist cDNA</span>
          </button>
        </div>

        <div style={{ overflowX: 'auto' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '13px', textAlign: 'left' }}>
            <thead>
              <tr style={{ borderBottom: '1px solid var(--border-color)', color: 'var(--text-secondary)' }}>
                <th style={{ padding: '10px 12px' }}>Variant Notation</th>
                <th style={{ padding: '10px 12px' }}>Pareto Rank</th>
                <th style={{ padding: '10px 12px' }}>Predicted Fitness</th>
                <th style={{ padding: '10px 12px' }}>Developability Score</th>
                <th style={{ padding: '10px 12px' }}>Synthesizability</th>
                <th style={{ padding: '10px 12px' }}>Weighted Utility</th>
                <th style={{ padding: '10px 12px' }}>Clinical Liabilities Status</th>
              </tr>
            </thead>
            <tbody>
              {data?.frontier_points?.map((p, idx) => (
                <tr
                  key={idx}
                  style={{
                    borderBottom: '1px solid rgba(255,255,255,0.04)',
                    background: selectedPoint?.notation === p.notation ? 'rgba(56, 189, 248, 0.08)' : 'transparent',
                    cursor: 'pointer',
                  }}
                  onClick={() => setSelectedPoint(p)}
                >
                  <td style={{ padding: '10px 12px', fontWeight: 600, color: '#f59e0b' }}>
                    {p.notation}
                  </td>
                  <td style={{ padding: '10px 12px' }}>
                    <span className="badge badge-warning" style={{ fontSize: '10px', padding: '2px 6px' }}>
                      Rank 0 (Optimal)
                    </span>
                  </td>
                  <td style={{ padding: '10px 12px', color: '#38bdf8', fontWeight: 600 }}>
                    {p.fitness.toFixed(3)}
                  </td>
                  <td style={{ padding: '10px 12px' }}>
                    <span style={{ color: p.developability_score > 0.85 ? '#10b981' : '#f59e0b', fontWeight: 600 }}>
                      {(p.developability_score * 100).toFixed(0)}%
                    </span>
                  </td>
                  <td style={{ padding: '10px 12px' }}>
                    {(p.synthesizability_score * 100).toFixed(0)}%
                  </td>
                  <td style={{ padding: '10px 12px', color: '#10b981', fontWeight: 700 }}>
                    {p.weighted_utility.toFixed(3)}
                  </td>
                  <td style={{ padding: '10px 12px' }}>
                    {p.liabilities?.clinical_liabilities?.length > 0 ? (
                      <span style={{ color: '#f59e0b', fontSize: '11px', display: 'flex', alignItems: 'center', gap: '4px' }}>
                        <AlertTriangle size={12} />
                        {p.liabilities.clinical_liabilities[0]}
                      </span>
                    ) : (
                      <span style={{ color: '#10b981', fontSize: '11px', display: 'flex', alignItems: 'center', gap: '4px' }}>
                        <CheckCircle2 size={12} />
                        Zero Liabilities
                      </span>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

