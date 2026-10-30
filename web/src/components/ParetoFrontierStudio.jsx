import React, { useState, useEffect, useMemo } from 'react';
import { Sliders, AlertTriangle, CheckCircle2, RefreshCw, ArrowRight } from 'lucide-react';

export default function ParetoFrontierStudio({
  targetAccession = 'TRPB_TM9D8',
  onRefreshLedger,
  onNavigateTab,
}) {
  const [loading, setLoading] = useState(false);
  const [data, setData] = useState(null);
  const [selectedPoint, setSelectedPoint] = useState(null);

  // 4 Objective Weights (0 to 100)
  const [wFitness, setWFitness] = useState(40);
  const [wDevelopability, setWDevelopability] = useState(30);
  const [wSynthesizability, setWSynthesizability] = useState(15);
  const [wStability, setWStability] = useState(15);

  const fetchParetoFrontier = async () => {
    setLoading(true);
    try {
      const payload = {
        target_protein: targetAccession,
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
      }
    } catch (err) {
      console.error('Network error fetching Pareto frontier:', err);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchParetoFrontier();
  }, [targetAccession]);

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

  const chartPoints = useMemo(() => {
    if (!data?.chart_data || data.chart_data.length === 0) return [];
    const xs = data.chart_data.map((p) => p.x_fitness);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const spanX = maxX - minX > 0.001 ? maxX - minX : 1.0;

    return data.chart_data.map((p) => {
      const normX = (p.x_fitness - minX) / spanX;
      const normY = p.y_developability;

      const svgX = 55 + normX * 405;
      const svgY = 240 - normY * 195;

      return {
        ...p,
        svgX,
        svgY,
      };
    });
  }, [data]);

  const frontierLinePoints = useMemo(() => {
    const front = chartPoints.filter((p) => p.is_frontier);
    return front.sort((a, b) => a.svgX - b.svgX);
  }, [chartPoints]);

  const frontierPathD = useMemo(() => {
    if (frontierLinePoints.length < 2) return '';
    return frontierLinePoints.reduce((acc, curr, idx) => {
      return idx === 0 ? `M ${curr.svgX} ${curr.svgY}` : `${acc} L ${curr.svgX} ${curr.svgY}`;
    }, '');
  }, [frontierLinePoints]);

  return (
    <div style={{ display: 'grid', gridTemplateColumns: '340px 1fr', gap: '20px', alignItems: 'start' }}>
      {/* Left Column: Multi-Objective Weights + Compact Summary */}
      <div className="glass-panel" style={{ padding: '22px' }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '16px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <Sliders size={16} color="var(--text-main)" />
            <h2 style={{ fontSize: '15px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Objective Weights
            </h2>
          </div>
          <span className="mono" style={{ fontSize: '11px', color: 'var(--text-dim)' }}>
            {targetAccession}
          </span>
        </div>

        {/* Presets */}
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '6px', marginBottom: '18px' }}>
          <button onClick={() => applyPreset('therapeutic')} className="btn btn-secondary" style={{ fontSize: '11.5px', padding: '6px' }}>
            Therapeutic CMC
          </button>
          <button onClick={() => applyPreset('affinity')} className="btn btn-secondary" style={{ fontSize: '11.5px', padding: '6px' }}>
            Max Affinity
          </button>
          <button onClick={() => applyPreset('balanced')} className="btn btn-secondary" style={{ fontSize: '11.5px', padding: '6px' }}>
            Balanced
          </button>
          <button onClick={() => applyPreset('stability')} className="btn btn-secondary" style={{ fontSize: '11.5px', padding: '6px' }}>
            High Stability
          </button>
        </div>

        {/* Sliders */}
        <div style={{ display: 'flex', flexDirection: 'column', gap: '15px' }}>
          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
              <span style={{ color: 'var(--text-muted)' }}>Binding Affinity / Fitness</span>
              <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>{wFitness}%</span>
            </div>
            <input
              type="range"
              min="5"
              max="90"
              value={wFitness}
              onChange={(e) => setWFitness(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--text-main)' }}
            />
          </div>

          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
              <span style={{ color: 'var(--text-muted)' }}>CMC Developability</span>
              <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>{wDevelopability}%</span>
            </div>
            <input
              type="range"
              min="5"
              max="90"
              value={wDevelopability}
              onChange={(e) => setWDevelopability(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--text-main)' }}
            />
          </div>

          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
              <span style={{ color: 'var(--text-muted)' }}>Synthesizability</span>
              <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>{wSynthesizability}%</span>
            </div>
            <input
              type="range"
              min="5"
              max="90"
              value={wSynthesizability}
              onChange={(e) => setWSynthesizability(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--text-main)' }}
            />
          </div>

          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '12px', marginBottom: '4px' }}>
              <span style={{ color: 'var(--text-muted)' }}>Conformational Stability</span>
              <span className="mono tabular" style={{ fontWeight: 600, color: 'var(--text-main)' }}>{wStability}%</span>
            </div>
            <input
              type="range"
              min="5"
              max="90"
              value={wStability}
              onChange={(e) => setWStability(Number(e.target.value))}
              style={{ width: '100%', accentColor: 'var(--text-main)' }}
            />
          </div>

          <button
            onClick={fetchParetoFrontier}
            disabled={loading}
            className="btn btn-primary"
            style={{ marginTop: '4px', width: '100%', padding: '10px', fontSize: '13px' }}
          >
            <RefreshCw size={14} />
            <span>{loading ? 'Re-ranking...' : 'Re-rank Frontier'}</span>
          </button>
        </div>

        {/* Compact Inline Frontier Metrics */}
        <div
          style={{
            marginTop: '18px',
            paddingTop: '16px',
            borderTop: '1px solid var(--border-subtle)',
            display: 'grid',
            gridTemplateColumns: '1fr 1fr',
            gap: '10px',
          }}
        >
          <div style={{ padding: '10px 12px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
            <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Rank 0 Leads</div>
            <div className="mono tabular" style={{ fontSize: '17px', fontWeight: 700, color: 'var(--text-main)' }}>
              {data?.rank_0_count ?? 0} <span style={{ fontSize: '11px', color: 'var(--text-dim)', fontWeight: 400 }}>/ {data?.points?.length ?? 0}</span>
            </div>
          </div>
          <div style={{ padding: '10px 12px', borderRadius: '8px', background: 'var(--bg-elevated)', border: '1px solid var(--border-subtle)' }}>
            <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Hypervolume</div>
            <div className="mono tabular" style={{ fontSize: '17px', fontWeight: 700, color: '#0d9488' }}>
              {data?.hypervolume ? data.hypervolume.toFixed(3) : '0.000'}
            </div>
          </div>
        </div>
      </div>

      {/* Right Column: Interactive Trade-off Plot + Rank 0 Table */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
        <div className="glass-panel" style={{ padding: '20px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px', flexWrap: 'wrap', gap: '8px' }}>
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Fitness vs. Biologics Developability Frontier
            </h3>
            <div style={{ display: 'flex', gap: '14px', fontSize: '11.5px', color: 'var(--text-muted)' }}>
              <span style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                <span style={{ width: '9px', height: '9px', borderRadius: '50%', background: '#0d9488' }} />
                Rank 0 Optimal
              </span>
              <span style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                <span style={{ width: '8px', height: '8px', borderRadius: '50%', background: '#8c8a84' }} />
                Dominated
              </span>
            </div>
          </div>

          {/* Interactive SVG Chart */}
          <div
            style={{
              position: 'relative',
              width: '100%',
              height: '285px',
              background: 'var(--bg-elevated)',
              borderRadius: '8px',
              border: '1px solid var(--border-subtle)',
              overflow: 'hidden',
            }}
          >
            <svg viewBox="0 0 500 280" style={{ width: '100%', height: '100%' }}>
              <line x1="55" y1="45" x2="465" y2="45" stroke="#e4e0d7" strokeDasharray="4 4" />
              <line x1="55" y1="142" x2="465" y2="142" stroke="#e4e0d7" strokeDasharray="4 4" />
              <line x1="55" y1="240" x2="465" y2="240" stroke="#c8c3b8" />
              <line x1="55" y1="45" x2="55" y2="240" stroke="#c8c3b8" />

              <text x="260" y="268" fill="#575653" fontSize="10.5" textAnchor="middle">
                Predicted Fitness / Binding Affinity →
              </text>
              <text x="20" y="142" fill="#575653" fontSize="10.5" textAnchor="middle" transform="rotate(-90 20 142)">
                Developability Score →
              </text>

              {frontierPathD && (
                <path
                  d={frontierPathD}
                  fill="none"
                  stroke="#0d9488"
                  strokeWidth="1.8"
                  strokeDasharray="4 3"
                />
              )}

              {chartPoints.map((pt, i) => {
                const isSelected = selectedPoint?.notation === pt.label;
                const isFrontier = pt.is_frontier;
                const fillColor = isFrontier ? '#0d9488' : '#8c8a84';
                const r = isSelected ? 7.5 : isFrontier ? 5.5 : 4;

                return (
                  <g
                    key={i}
                    onClick={() => setSelectedPoint(data?.points?.find((p) => p.notation === pt.label))}
                    style={{ cursor: 'pointer' }}
                  >
                    {isFrontier && (
                      <circle cx={pt.svgX} cy={pt.svgY} r={r + 3} fill="rgba(13, 148, 136, 0.15)" />
                    )}
                    <circle
                      cx={pt.svgX}
                      cy={pt.svgY}
                      r={r}
                      fill={fillColor}
                      stroke={isSelected ? '#141413' : '#fdfcfa'}
                      strokeWidth={isSelected ? 2 : 1.2}
                    />
                    {isFrontier && (
                      <text
                        x={pt.svgX}
                        y={pt.svgY - 10}
                        fill="#141413"
                        fontSize="9"
                        textAnchor="middle"
                        fontWeight="600"
                        fontFamily="monospace"
                      >
                        {pt.label}
                      </text>
                    )}
                  </g>
                );
              })}
            </svg>
          </div>

          {/* Selected Lead Inspector */}
          {selectedPoint && (
            <div
              style={{
                marginTop: '12px',
                padding: '12px 16px',
                background: 'var(--bg-elevated)',
                borderRadius: '8px',
                border: '1px solid var(--border-subtle)',
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                flexWrap: 'wrap',
                gap: '12px',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: '16px', flexWrap: 'wrap' }}>
                <span className="mono" style={{ fontSize: '14px', fontWeight: 700, color: 'var(--text-main)' }}>
                  {selectedPoint.notation}
                </span>
                <span style={{ fontSize: '12px', color: 'var(--text-muted)' }}>
                  Fitness: <strong className="mono tabular" style={{ color: 'var(--text-main)' }}>{selectedPoint.fitness.toFixed(3)}</strong>
                </span>
                <span style={{ fontSize: '12px', color: 'var(--text-muted)' }}>
                  Developability: <strong className="mono tabular" style={{ color: '#0d9488' }}>{(selectedPoint.developability_score * 100).toFixed(0)}%</strong>
                </span>
                <span style={{ fontSize: '12px', color: 'var(--text-muted)' }}>
                  Utility: <strong className="mono tabular" style={{ color: 'var(--text-main)' }}>{selectedPoint.weighted_utility.toFixed(3)}</strong>
                </span>
              </div>

              <div style={{ fontSize: '11.5px' }}>
                {selectedPoint.liabilities?.clinical_liabilities?.length > 0 ? (
                  <span style={{ display: 'flex', alignItems: 'center', gap: '5px', color: '#b45309' }}>
                    <AlertTriangle size={13} />
                    <span>{selectedPoint.liabilities.clinical_liabilities[0]}</span>
                  </span>
                ) : (
                  <span style={{ display: 'flex', alignItems: 'center', gap: '5px', color: '#0d9488', fontWeight: 600 }}>
                    <CheckCircle2 size={13} />
                    <span>Zero CMC Liabilities</span>
                  </span>
                )}
              </div>
            </div>
          )}
        </div>

        {/* Compact Rank 0 Table */}
        <div className="glass-panel" style={{ padding: '18px 20px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Rank 0 Pareto-Optimal Candidates
            </h3>
            {onNavigateTab && (
              <button
                onClick={() => onNavigateTab('wetlab')}
                className="btn btn-primary"
                style={{ padding: '6px 12px', fontSize: '12px' }}
              >
                <span>Order Rank 0 Leads</span>
                <ArrowRight size={13} />
              </button>
            )}
          </div>

          <div style={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px', textAlign: 'left' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-dim)' }}>
                  <th style={{ padding: '8px' }}>Variant</th>
                  <th style={{ padding: '8px' }}>Fitness</th>
                  <th style={{ padding: '8px' }}>Developability</th>
                  <th style={{ padding: '8px' }}>Synthesizability</th>
                  <th style={{ padding: '8px' }}>Utility</th>
                  <th style={{ padding: '8px' }}>Liability Check</th>
                </tr>
              </thead>
              <tbody>
                {data?.frontier_points?.map((p, idx) => (
                  <tr
                    key={idx}
                    style={{
                      borderBottom: '1px solid var(--border-subtle)',
                      background:
                        selectedPoint?.notation === p.notation ? 'var(--bg-elevated)' : 'transparent',
                      cursor: 'pointer',
                    }}
                    onClick={() => setSelectedPoint(p)}
                  >
                    <td className="mono" style={{ padding: '8px', fontWeight: 600, color: 'var(--text-main)' }}>
                      {p.notation}
                    </td>
                    <td className="mono tabular" style={{ padding: '8px', fontWeight: 600 }}>
                      {p.fitness.toFixed(3)}
                    </td>
                    <td className="mono tabular" style={{ padding: '8px', color: '#0d9488', fontWeight: 600 }}>
                      {(p.developability_score * 100).toFixed(0)}%
                    </td>
                    <td className="mono tabular" style={{ padding: '8px' }}>
                      {(p.synthesizability_score * 100).toFixed(0)}%
                    </td>
                    <td className="mono tabular" style={{ padding: '8px', fontWeight: 700, color: 'var(--text-main)' }}>
                      {p.weighted_utility.toFixed(3)}
                    </td>
                    <td style={{ padding: '8px' }}>
                      {p.liabilities?.clinical_liabilities?.length > 0 ? (
                        <span style={{ color: '#b45309', fontSize: '11px' }}>
                          {p.liabilities.clinical_liabilities[0]}
                        </span>
                      ) : (
                        <span style={{ color: '#0d9488', fontSize: '11px', fontWeight: 600 }}>
                          Clean
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
    </div>
  );
}
