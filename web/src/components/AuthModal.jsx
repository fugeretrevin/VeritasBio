import React, { useState } from 'react';
import { X, ShieldCheck, KeyRound, Check, Lock } from 'lucide-react';

export default function AuthModal({ isOpen, onClose, currentUser, onSelectUser }) {
  const [email, setEmail] = useState('lead.scientist@panoptesbio.com');
  const [role, setRole] = useState('lead_biologist');
  const [actionMeaning, setActionMeaning] = useState('Approval of Modal A10G Compute Spend');
  const [generatedSig, setGeneratedSig] = useState(null);
  const [isSigning, setIsSigning] = useState(false);

  if (!isOpen) return null;

  const handleApplyIdentity = () => {
    const cleanUser = email.split('@')[0] || 'scientist';
    const user = {
      uid: `usr_${cleanUser}`,
      email: email,
      display_name: `Dr. ${cleanUser.replace('.', ' ').toUpperCase()}`,
      organization_id: 'panoptesbio',
      role: role,
      is_part11_certified: true,
    };
    onSelectUser(user);
  };

  const handleGenerateSignature = async () => {
    setIsSigning(true);
    try {
      const res = await fetch('/api/auth/sign', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          id_token: `mock_token_${role}_${email.split('@')[0]}`,
          action: actionMeaning,
          meaning: 'Authenticated operator electronic signature binding action under 21 CFR Part 11 §11.200',
        }),
      });
      const data = await res.json();
      setGeneratedSig(data.signature);
    } catch (err) {
      alert(`Signature failed: ${err.message}`);
    } finally {
      setIsSigning(false);
    }
  };

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        width: '100vw',
        height: '100vh',
        backgroundColor: 'rgba(7, 10, 19, 0.82)',
        backdropFilter: 'blur(8px)',
        zIndex: 1000,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '20px',
      }}
    >
      <div
        className="glass-panel"
        style={{
          width: '100%',
          maxWidth: '560px',
          padding: '28px',
          position: 'relative',
          border: '1px solid var(--border-active)',
          boxShadow: '0 20px 60px rgba(0, 0, 0, 0.6)',
        }}
      >
        <button
          onClick={onClose}
          style={{
            position: 'absolute',
            top: '16px',
            right: '16px',
            background: 'none',
            border: 'none',
            color: 'var(--text-muted)',
            cursor: 'pointer',
          }}
        >
          <X size={20} />
        </button>

        <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
          <ShieldCheck size={24} color="var(--cyan-bright)" />
          <h2 className="title-display" style={{ fontSize: '20px', margin: 0 }}>
            PanoptesBio Cloud Identity & FDA Part 11
          </h2>
        </div>
        <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '20px' }}>
          Multi-tenant Firebase Authentication (`panoptesbio.firebaseapp.com`) & 21 CFR §11.200
          dual-identifier electronic signatures.
        </p>

        {/* Corporate Email */}
        <div style={{ marginBottom: '14px' }}>
          <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, marginBottom: '6px' }}>
            Operator Corporate Email
          </label>
          <input
            type="email"
            className="input-field"
            value={email}
            onChange={(e) => setEmail(e.target.value)}
          />
        </div>

        {/* FDA Part 11 Role */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, marginBottom: '6px' }}>
            FDA Part 11 Role & Privileges
          </label>
          <select
            className="input-field"
            value={role}
            onChange={(e) => setRole(e.target.value)}
          >
            <option value="lead_biologist">Lead Computational Biologist (Full Optimization Scope)</option>
            <option value="lab_director">Laboratory Director (Budget Approval & CRO Sign-Off)</option>
            <option value="compliance_officer">Compliance Officer (Merkle Audit Trail Verification)</option>
          </select>
        </div>

        <div style={{ display: 'flex', gap: '10px', marginBottom: '24px' }}>
          <button onClick={handleApplyIdentity} className="btn btn-primary" style={{ flex: 1 }}>
            <Check size={16} />
            <span>Apply Operator Identity</span>
          </button>
        </div>

        {/* Signature Generator */}
        <div
          style={{
            borderTop: '1px solid var(--border-subtle)',
            paddingTop: '16px',
          }}
        >
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '10px' }}>
            <Lock size={16} color="var(--amber-bright)" />
            <h4 style={{ fontSize: '13px', margin: 0, fontWeight: 600 }}>
              Generate 21 CFR §11.200 Electronic Signature
            </h4>
          </div>

          <div style={{ marginBottom: '12px' }}>
            <select
              className="input-field"
              value={actionMeaning}
              onChange={(e) => setActionMeaning(e.target.value)}
            >
              <option value="Approval of Modal A10G Compute Spend">Approval of Modal A10G Compute Spend</option>
              <option value="Authorization of Twist cDNA Synthesis Manifest">Authorization of Twist cDNA Synthesis Manifest</option>
              <option value="Approval of Opentrons OT-2 Neutralization Protocol">Approval of Opentrons OT-2 Neutralization Protocol</option>
              <option value="Certification of FDA Part 11 Audit Trail Integrity">Certification of FDA Part 11 Audit Trail Integrity</option>
            </select>
          </div>

          <button
            onClick={handleGenerateSignature}
            disabled={isSigning}
            className="btn btn-secondary"
            style={{ width: '100%', marginBottom: '12px' }}
          >
            <KeyRound size={14} color="var(--cyan-bright)" />
            <span>{isSigning ? 'Computing HMAC-SHA256...' : 'Sign Electronic Proof'}</span>
          </button>

          {generatedSig && (
            <div
              className="mono"
              style={{
                fontSize: '11px',
                padding: '12px',
                borderRadius: '8px',
                background: '#060911',
                border: '1px solid var(--border-subtle)',
                color: 'var(--emerald-bright)',
                wordBreak: 'break-all',
              }}
            >
              <div>Operator: {generatedSig.operator_email} ({generatedSig.role})</div>
              <div>Timestamp: {generatedSig.signed_timestamp}</div>
              <div style={{ marginTop: '4px', color: 'var(--cyan-bright)' }}>
                HMAC-SHA256: {generatedSig.electronic_signature_hash}
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
