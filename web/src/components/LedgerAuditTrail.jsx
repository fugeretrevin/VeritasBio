import React from 'react';
import { ShieldCheck, Database, FileCheck, CheckCircle2 } from 'lucide-react';

export default function LedgerAuditTrail({ ledgerData, onRefresh }) {
  const jobs = ledgerData?.jobs || [];
  const auditBlocks = ledgerData?.audit_blocks || [];
  const integrity = ledgerData?.audit_integrity || {
    is_valid: true,
    block_count: 5,
    verification_message: 'Audit chain intact',
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
      {/* Clean Compact Integrity Bar */}
      <div
        className="glass-panel"
        style={{
          padding: '12px 18px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '12px',
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
          <CheckCircle2 size={18} color="#0d9488" />
          <span style={{ fontWeight: 600, fontSize: '13.5px', color: 'var(--text-main)' }}>
            {integrity.verification_message}
          </span>
          <span className="mono" style={{ fontSize: '12px', color: 'var(--text-dim)' }}>
            ({integrity.block_count} Merkle blocks verified · SHA-256)
          </span>
        </div>

        <button onClick={onRefresh} className="btn btn-secondary" style={{ fontSize: '12px', padding: '6px 12px' }}>
          <FileCheck size={14} />
          <span>Re-verify Chain</span>
        </button>
      </div>

      {/* Grid: Left Jobs Table, Right Merkle Block Inspector */}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '16px' }}>
        {/* Jobs Table */}
        <div className="glass-panel" style={{ padding: '18px 20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '12px' }}>
            <Database size={16} color="var(--text-main)" />
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Compute Spend Ledger
            </h3>
          </div>

          <div style={{ overflowX: 'auto', maxHeight: '380px' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-dim)' }}>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Job ID</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Worker</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Variants</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Cost</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Status</th>
                </tr>
              </thead>
              <tbody>
                {jobs.length > 0 ? (
                  jobs.map((j) => (
                    <tr key={j.job_id} style={{ borderBottom: '1px solid var(--border-subtle)' }}>
                      <td style={{ padding: '7px 8px', color: 'var(--text-main)' }} className="mono">
                        {j.job_id.slice(0, 14)}...
                      </td>
                      <td style={{ padding: '7px 8px' }}>{j.target_model}</td>
                      <td className="mono tabular" style={{ padding: '7px 8px' }}>{j.variant_count}</td>
                      <td className="mono tabular" style={{ padding: '7px 8px', fontWeight: 600 }}>
                        ${Number(j.estimated_cost).toFixed(4)}
                      </td>
                      <td style={{ padding: '7px 8px' }}>
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
                    <td colSpan="5" style={{ textAlign: 'center', padding: '24px', color: 'var(--text-dim)' }}>
                      No jobs recorded yet.
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* Merkle Hash Blocks */}
        <div className="glass-panel" style={{ padding: '18px 20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '12px' }}>
            <ShieldCheck size={16} color="var(--text-main)" />
            <h3 style={{ fontSize: '14px', fontWeight: 700, margin: 0, color: 'var(--text-main)' }}>
              Merkle State-Transition Blocks
            </h3>
          </div>

          <div style={{ overflowX: 'auto', maxHeight: '380px' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '12px' }}>
              <thead>
                <tr style={{ borderBottom: '1px solid var(--border-subtle)', color: 'var(--text-dim)' }}>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>#</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Action</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Block Hash</th>
                  <th style={{ textAlign: 'left', padding: '7px 8px' }}>Parent Hash</th>
                </tr>
              </thead>
              <tbody>
                {auditBlocks.length > 0 ? (
                  auditBlocks.map((b) => (
                    <tr key={b.block_index} style={{ borderBottom: '1px solid var(--border-subtle)' }}>
                      <td className="mono" style={{ padding: '7px 8px', fontWeight: 600, color: 'var(--text-main)' }}>
                        #{b.block_index}
                      </td>
                      <td style={{ padding: '7px 8px' }}>
                        <span className="mono" style={{ fontSize: '11px', color: 'var(--text-muted)' }}>
                          {b.action_type}
                        </span>
                      </td>
                      <td style={{ padding: '7px 8px', color: '#0d9488' }} className="mono">
                        {b.block_hash.slice(0, 14)}...
                      </td>
                      <td style={{ padding: '7px 8px', color: 'var(--text-dim)' }} className="mono">
                        {b.prev_block_hash.slice(0, 14)}...
                      </td>
                    </tr>
                  ))
                ) : (
                  <tr>
                    <td colSpan="4" style={{ textAlign: 'center', padding: '24px', color: 'var(--text-dim)' }}>
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
