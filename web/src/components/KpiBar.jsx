import React from 'react';
import { DollarSign, Flame, Award, TrendingUp } from 'lucide-react';

export default function KpiBar({
  remainingBudget = 5.0,
  cumulativeSpent = 0.0,
  topFitness = 0.0,
  savedPct = '99.5%',
  speedup = '184.3x Cheaper',
}) {
  return (
    <div
      style={{
        display: 'grid',
        gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))',
        gap: '16px',
        marginBottom: '20px',
      }}
    >
      {/* Metric 1: Remaining Budget */}
      <div className="glass-panel" style={{ padding: '16px 20px', position: 'relative', overflow: 'hidden' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <div style={{ fontSize: '12px', color: 'var(--text-muted)', fontWeight: 500, marginBottom: '4px' }}>
              REMAINING BUDGET
            </div>
            <div className="title-display mono" style={{ fontSize: '26px', color: 'var(--cyan-bright)' }}>
              ${Number(remainingBudget).toFixed(4)}
            </div>
          </div>
          <div
            style={{
              padding: '8px',
              borderRadius: '8px',
              background: 'rgba(2, 132, 199, 0.15)',
              color: 'var(--cyan-bright)',
            }}
          >
            <DollarSign size={20} />
          </div>
        </div>
        <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '6px' }}>
          Exact micro-USD integer precision
        </div>
      </div>

      {/* Metric 2: Cumulative Spend */}
      <div className="glass-panel" style={{ padding: '16px 20px', position: 'relative', overflow: 'hidden' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <div style={{ fontSize: '12px', color: 'var(--text-muted)', fontWeight: 500, marginBottom: '4px' }}>
              CUMULATIVE CLOUD SPEND
            </div>
            <div className="title-display mono" style={{ fontSize: '26px', color: 'var(--text-main)' }}>
              ${Number(cumulativeSpent).toFixed(4)}
            </div>
          </div>
          <div
            style={{
              padding: '8px',
              borderRadius: '8px',
              background: 'rgba(245, 158, 11, 0.15)',
              color: 'var(--amber-bright)',
            }}
          >
            <Flame size={20} />
          </div>
        </div>
        <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '6px' }}>
          Modal A10G cold/warm amortized
        </div>
      </div>

      {/* Metric 3: Discovered Top Fitness */}
      <div className="glass-panel" style={{ padding: '16px 20px', position: 'relative', overflow: 'hidden' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <div style={{ fontSize: '12px', color: 'var(--text-muted)', fontWeight: 500, marginBottom: '4px' }}>
              TOP DISCOVERED FITNESS
            </div>
            <div className="title-display mono" style={{ fontSize: '26px', color: 'var(--emerald-bright)' }}>
              +{Number(topFitness).toFixed(4)}
            </div>
          </div>
          <div
            style={{
              padding: '8px',
              borderRadius: '8px',
              background: 'rgba(16, 185, 129, 0.15)',
              color: 'var(--emerald-bright)',
            }}
          >
            <Award size={20} />
          </div>
        </div>
        <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '6px' }}>
          Lead-locked epistatic optimum
        </div>
      </div>

      {/* Metric 4: Cloud Valuation Saved */}
      <div className="glass-panel" style={{ padding: '16px 20px', position: 'relative', overflow: 'hidden' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <div style={{ fontSize: '12px', color: 'var(--text-muted)', fontWeight: 500, marginBottom: '4px' }}>
              CLOUD CAPITAL SAVINGS
            </div>
            <div className="title-display" style={{ fontSize: '24px', color: '#0284c7' }}>
              {savedPct}{' '}
              <span style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-muted)' }}>
                ({speedup})
              </span>
            </div>
          </div>
          <div
            style={{
              padding: '8px',
              borderRadius: '8px',
              background: 'rgba(2, 132, 199, 0.1)',
              color: '#0284c7',
            }}
          >
            <TrendingUp size={20} />
          </div>
        </div>
        <div style={{ fontSize: '12px', color: 'var(--text-dim)', marginTop: '6px' }}>
          Compared to unguided brute force
        </div>
      </div>
    </div>
  );
}
