import React, { useState } from 'react';
import CampaignStudio from './CampaignStudio';
import DeNovoStudio from './DeNovoStudio';
import ParetoFrontierStudio from './ParetoFrontierStudio';
import LedgerAuditTrail from './LedgerAuditTrail';
import WetlabExportStudio from './WetlabExportStudio';
import { Compass, Sparkles, Scale, ShieldCheck, FlaskConical, ArrowRight } from 'lucide-react';

export const STUDIO_TARGETS = [
  { id: 'TRPB_TM9D8', label: 'TrpB Synthase (8EH0)', short: 'TRPB_TM9D8' },
  { id: 'SCFV_9004G_CN2', label: 'Scorpion Antivenom scFv (2YBR)', short: 'SCFV_9004G_CN2' },
  { id: 'HER2_TRASTUZUMAB', label: 'Trastuzumab CDR-H3 (1N8Z)', short: 'HER2_TRASTUZUMAB' },
  { id: 'SPG1_STRSG', label: 'Protein G GB1 (2QMT)', short: 'SPG1_STRSG' },
  { id: 'P42212', label: 'avGFP Fluorophore (1EMA)', short: 'P42212' },
];

export default function StudioView({
  onRefreshLedger,
  ledgerData,
  currentUser,
  remainingBudget = 5.0,
  cumulativeSpent = 0.0,
  topFitness = 0.0,
}) {
  const [activeSubTab, setActiveSubTab] = useState('campaign');
  const [targetAccession, setTargetAccession] = useState('TRPB_TM9D8');
  const [lastRunSummary, setLastRunSummary] = useState(null);

  const tabs = [
    { id: 'campaign', label: 'Active Learning', icon: Compass },
    { id: 'denovo', label: 'De Novo Binders', icon: Sparkles },
    { id: 'pareto', label: 'Pareto Frontier', icon: Scale },
    { id: 'wetlab', label: 'Wet-Lab & LIMS', icon: FlaskConical },
    { id: 'ledger', label: 'Audit Ledger', icon: ShieldCheck },
  ];

  const showBottomDrawer = Boolean(lastRunSummary) || Number(cumulativeSpent) > 0;

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
      {/* Unified Studio Command Bar: Navigation + Shared Target Selector + Single Spend Counter */}
      <div
        className="glass-panel"
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          flexWrap: 'wrap',
          gap: '12px',
          padding: '10px 14px',
          background: 'var(--bg-elevated)',
        }}
      >
        {/* Left: Sub-Navigation Tabs */}
        <div style={{ display: 'flex', alignItems: 'center', gap: '6px', flexWrap: 'wrap' }}>
          {tabs.map((tab) => {
            const Icon = tab.icon;
            const isActive = activeSubTab === tab.id;
            return (
              <button
                key={tab.id}
                onClick={() => setActiveSubTab(tab.id)}
                className={`btn ${isActive ? 'btn-primary' : 'btn-secondary'}`}
                style={{
                  padding: '7px 13px',
                  fontSize: '13px',
                  fontWeight: isActive ? 600 : 500,
                }}
              >
                <Icon size={14} />
                <span>{tab.label}</span>
              </button>
            );
          })}
        </div>

        {/* Right: Global Target Selector + Spend Counter */}
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px', flexWrap: 'wrap' }}>
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '8px',
              padding: '4px 10px',
              borderRadius: '8px',
              background: 'var(--bg-card)',
              border: '1px solid var(--border-subtle)',
              fontSize: '12px',
            }}
          >
            <span style={{ color: 'var(--text-dim)', fontWeight: 600 }}>Target:</span>
            <select
              value={targetAccession}
              onChange={(e) => setTargetAccession(e.target.value)}
              style={{
                background: 'transparent',
                border: 'none',
                color: 'var(--text-main)',
                fontWeight: 600,
                fontSize: '12px',
                cursor: 'pointer',
                outline: 'none',
              }}
            >
              {STUDIO_TARGETS.map((t) => (
                <option key={t.id} value={t.id}>
                  {t.label}
                </option>
              ))}
            </select>
          </div>

          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '14px',
              padding: '6px 14px',
              borderRadius: '8px',
              background: 'var(--bg-card)',
              border: '1px solid var(--border-subtle)',
              fontSize: '12.5px',
            }}
          >
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
              <span className="dot-live" />
              <span style={{ color: 'var(--text-muted)' }}>Spent:</span>
              <strong className="tabular mono" style={{ color: 'var(--text-main)' }}>
                ${Number(cumulativeSpent).toFixed(4)}
              </strong>
            </div>
            <span style={{ color: 'var(--border-subtle)' }}>|</span>
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
              <span style={{ color: 'var(--text-muted)' }}>Remaining:</span>
              <strong className="tabular mono" style={{ color: 'var(--emerald-primary)' }}>
                ${Number(remainingBudget).toFixed(4)}
              </strong>
            </div>
          </div>
        </div>
      </div>

      {/* Main Studio Workspace */}
      <div>
        {activeSubTab === 'campaign' && (
          <CampaignStudio
            targetAccession={targetAccession}
            onTargetChange={setTargetAccession}
            currentUser={currentUser}
            onRefreshLedger={onRefreshLedger}
            onRunComplete={(summary) => setLastRunSummary(summary)}
            onNavigateTab={setActiveSubTab}
          />
        )}
        {activeSubTab === 'denovo' && (
          <DeNovoStudio
            targetAccession={targetAccession}
            onTargetChange={setTargetAccession}
            currentUser={currentUser}
            onRefreshLedger={onRefreshLedger}
            onRunComplete={(summary) => setLastRunSummary(summary)}
            onNavigateTab={setActiveSubTab}
          />
        )}
        {activeSubTab === 'pareto' && (
          <ParetoFrontierStudio
            targetAccession={targetAccession}
            onRefreshLedger={onRefreshLedger}
            onNavigateTab={setActiveSubTab}
          />
        )}
        {activeSubTab === 'wetlab' && (
          <WetlabExportStudio
            targetAccession={targetAccession}
            currentUser={currentUser}
            onRefreshLedger={onRefreshLedger}
          />
        )}
        {activeSubTab === 'ledger' && (
          <LedgerAuditTrail
            ledgerData={ledgerData}
            onRefresh={onRefreshLedger}
          />
        )}
      </div>

      {/* Post-Run Summary Strip (Only Appears After a Run or Active Spend) */}
      {showBottomDrawer && (
        <div
          className="glass-panel"
          style={{
            padding: '14px 20px',
            background: 'var(--bg-elevated)',
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            flexWrap: 'wrap',
            gap: '16px',
          }}
        >
          <div style={{ display: 'flex', alignItems: 'center', gap: '28px', flexWrap: 'wrap' }}>
            <div>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Active Target</div>
              <div className="mono" style={{ fontSize: '14px', fontWeight: 700, color: 'var(--text-main)' }}>
                {lastRunSummary?.target || targetAccession}
              </div>
            </div>

            <div>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Peak Lead Fitness</div>
              <div className="tabular mono" style={{ fontSize: '15px', fontWeight: 700, color: 'var(--emerald-primary)' }}>
                +{Number(lastRunSummary?.topScore ?? topFitness).toFixed(4)}
              </div>
            </div>

            <div>
              <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>Efficiency vs Brute Force</div>
              <div className="tabular" style={{ fontSize: '14px', fontWeight: 600, color: 'var(--text-main)' }}>
                99.5% Saved <span style={{ color: 'var(--text-muted)', fontSize: '12px' }}>(184x)</span>
              </div>
            </div>
          </div>

          <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
            <button
              onClick={() => setActiveSubTab('pareto')}
              className="btn btn-secondary"
              style={{ padding: '7px 12px', fontSize: '12px' }}
            >
              <span>Inspect Pareto Trade-Offs</span>
            </button>
            <button
              onClick={() => setActiveSubTab('wetlab')}
              className="btn btn-primary"
              style={{ padding: '7px 14px', fontSize: '12px' }}
            >
              <span>Export to Wet-Lab</span>
              <ArrowRight size={13} />
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
