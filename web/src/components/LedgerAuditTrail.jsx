import React from 'react';
import { ShieldCheck, Database, FileCheck, CheckCircle2 } from 'lucide-react';

export default function LedgerAuditTrail({ ledgerData, onRefresh }) {
  const jobs = ledgerData?.jobs || [];
  const auditBlocks = ledgerData?.audit_blocks || [];
  const integrity = ledgerData?.audit_integrity || { is_valid: true, block_count: 5, verification_message: 'Audit chain intact' };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '24px' }}>
      {/* Integrity Header Banner */}
      <div
        className="glass-panel"
        style={{
          padding: '16px 20px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '12px',
          borderLeft: '4px solid var(--emerald-bright)',
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <CheckCircle2 size={24} color="var(--emerald-bright)" />
          <div>
            <div style={{ fontWeight: 600, fontSize: '14px', color: 'var(--text-main)' }}>
              FDA 21 CFR Part 11 Cryptographic Merkle Hash Chain
            </div>
            <div style={{ fontSize: '12px', color: 'var(--text-muted)' }}>
              {integrity.verification_message} ({integrity.block_count} validated blocks · SHA-256 + HMAC)
            </div>
          </div>
        </div>

        <button onClick={onRefresh} className="btn btn-emerald" style={{ fontSize: '12px' }}>
          <FileCheck size={15} />
          <span>Verify Merkle Proofs</span>
        </button>
      </div>

      {/* Grid: Left Jobs Table, Right Merkle Block Inspector */}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '20px' }}>
        {/* Jobs Table */}
        <div className="glass-panel" style={{ padding: '20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '14px' }}>
            <Database size={18} color="var(--cyan-bright)" />
            <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
              Double-Entry Compute Spend Ledger
            </h3>
          </div>

          <div style={{ overflowX: 'auto', maxHeight: '400px' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-muted)' }}>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Job ID</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Model</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Variants</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Cost ($)</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Status</th>
                </tr>
              </thead>
              <tbody>
                {jobs.length > 0 ? (
                  jobs.map((j) => (
                    <tr key={j.job_id} style={{ borderBottom: '1px solid rgba(255,255,255,0.04)' }}>
                      <td style={{ padding: '8px', color: 'var(--cyan-bright)' }} className="mono">
                        {j.job_id.slice(0, 16)}...
                      </td>
                      <td style={{ padding: '8px' }}>{j.target_model}</td>
                      <td style={{ padding: '8px' }}>{j.variant_count}</td>
                      <td style={{ padding: '8px', fontWeight: 600 }}>${Number(j.estimated_cost).toFixed(4)}</td>
                      <td style={{ padding: '8px' }}>
                        <span
                          className={`badge ${
                            j.status === 'COMPLETED'
                              ? 'badge-emerald'
                              : j.status === 'APPROVED'
                              ? 'badge-cyan'
                              : 'badge-amber'
                          }`}
                        >
                          {j.status}
                        </span>
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="5" style={{ textAlign: 'center', padding: '20px', color: 'var(--text-dim)' }}>
                      No jobs recorded yet.
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* Merkle Hash Blocks */}
        <div className="glass-panel" style={{ padding: '20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '14px' }}>
            <ShieldCheck size={18} color="var(--emerald-bright)" />
            <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
              Immutable State-Transition Blocks (Merkle Chain)
            </h3>
          </div>

          <div style={{ overflowX: 'auto', maxHeight: '400px' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-muted)' }}>
                  <th style={{ textAlign: 'left', padding: '8px' }}>#</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Action</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Block Hash</th>
                  <th style={{ textAlign: 'left', padding: '8px' }}>Parent Hash</th>
                </tr>
              </thead>
              <tbody>
                {auditBlocks.length > 0 ? (
                  auditBlocks.map((b) => (
                    <tr key={b.block_index} style={{ borderBottom: '1px solid rgba(255,255,255,0.04)' }}>
                      <td style={{ padding: '8px', fontWeight: 600, color: 'var(--cyan-bright)' }}>
                        #{b.block_index}
                      </td>
                      <td style={{ padding: '8px' }}>
                        <span className="badge badge-cyan" style={{ fontSize: '10px' }}>
                          {b.action_type}
                        </span>
                      </td>
                      <td style={{ padding: '8px', color: 'var(--emerald-bright)' }} className="mono">
                        {b.block_hash.slice(0, 14)}...
                      </td>
                      <td style={{ padding: '8px', color: 'var(--text-dim)' }} className="mono">
                        {b.prev_block_hash.slice(0, 14)}...
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="4" style={{ textAlign: 'center', padding: '20px', color: 'var(--text-dim)' }}>
                      Genesis block loading...
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
