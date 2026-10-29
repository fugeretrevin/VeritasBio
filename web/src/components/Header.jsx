import React from 'react';
import { ShieldCheck, UserCheck, RefreshCw, KeyRound, Cpu } from 'lucide-react';

export default function Header({
  currentUser,
  onOpenAuth,
  onRefresh,
  isRefreshing,
}) {
  return (
    <header
      className="glass-panel"
      style={{
        padding: '16px 24px',
        marginBottom: '20px',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        flexWrap: 'wrap',
        gap: '16px',
      }}
    >
      <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
        <div
          style={{
            width: '42px',
            height: '42px',
            borderRadius: '12px',
            background: 'linear-gradient(135deg, #0284c7 0%, #10b981 100%)',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            boxShadow: '0 4px 14px rgba(0, 0, 0, 0.4), inset 0 1px 0 rgba(255, 255, 255, 0.2)',
          }}
        >
          <Cpu size={24} color="#ffffff" />
        </div>

        <div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
            <h1 className="title-display" style={{ fontSize: '24px', margin: 0, letterSpacing: '-0.02em' }}>
              PANOPTES<span style={{ color: 'var(--cyan-bright)' }}>BIO</span>
            </h1>
            <span className="badge badge-cyan">v3.2 Cloud Engine</span>
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }} className="badge badge-emerald">
              <span className="dot-live" />
              <span>panoptesbio.firebaseapp.com</span>
            </div>
          </div>
          <p style={{ margin: 0, fontSize: '12px', color: 'var(--text-muted)' }}>
            Deterministic Cloud Economics & Biologics Evolution Engine · FDA 21 CFR Part 11 Cryptographic Merkle Chain
          </p>
        </div>
      </div>

      <div style={{ display: 'flex', alignItems: 'center', gap: '12px', flexWrap: 'wrap' }}>
        {/* Active Operator Identity Pill */}
        <div
          style={{
            padding: '6px 14px',
            display: 'flex',
            alignItems: 'center',
            gap: '8px',
            borderRadius: 'var(--radius-full)',
            background: 'rgba(15, 23, 42, 0.7)',
            border: '1px solid var(--border-subtle)',
          }}
        >
          <UserCheck size={16} color="var(--emerald-bright)" />
          <div style={{ fontSize: '12px' }}>
            <span style={{ fontWeight: 600, color: 'var(--text-main)' }}>
              {currentUser?.display_name || 'Dr. Lead Biologist'}
            </span>
            <span style={{ color: 'var(--text-dim)', marginLeft: '6px' }}>
              ({currentUser?.role || 'lead_biologist'})
            </span>
          </div>
          <ShieldCheck size={14} color="var(--cyan-bright)" title="FDA 21 CFR Part 11 Certified" />
        </div>

        <button
          onClick={onOpenAuth}
          className="btn btn-secondary"
          title="Switch Operator Identity or Generate 21 CFR §11.200 Electronic Signature"
        >
          <KeyRound size={15} />
          <span>Identity & Sign</span>
        </button>

        <button
          onClick={onRefresh}
          className="btn btn-secondary"
          disabled={isRefreshing}
          title="Sync Ledger State from SQLite"
        >
          <RefreshCw size={15} className={isRefreshing ? 'spin' : ''} />
          <span>Sync State</span>
        </button>
      </div>
    </header>
  );
}
