import React, { useState, useEffect } from 'react';
import { X, ShieldCheck, KeyRound, Lock, UserPlus, LogIn, LogOut, AlertCircle, CheckCircle2 } from 'lucide-react';

export default function AuthModal({
  isOpen,
  initialMode = 'signup',
  onClose,
  currentUser,
  onSelectUser,
  onSignOut,
}) {
  const [mode, setMode] = useState(initialMode); // 'signup' | 'signin' | 'profile'
  const [displayName, setDisplayName] = useState('');
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [organization, setOrganization] = useState('PanoptesBio Lab');
  const [role, setRole] = useState('lead_biologist');

  const [errorMsg, setErrorMsg] = useState(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  const [actionMeaning, setActionMeaning] = useState('Approval of Modal A10G Compute Spend');
  const [generatedSig, setGeneratedSig] = useState(null);
  const [isSigning, setIsSigning] = useState(false);

  useEffect(() => {
    if (isOpen) {
      setErrorMsg(null);
      if (currentUser?.uid) {
        setMode('profile');
      } else {
        setMode(initialMode === 'profile' ? 'signup' : initialMode);
      }
    }
  }, [isOpen, initialMode, currentUser]);

  if (!isOpen) return null;

  const handleSubmitAuth = async (e) => {
    e.preventDefault();
    setErrorMsg(null);
    setIsSubmitting(true);

    try {
      const endpoint = mode === 'signup' ? '/api/auth/register' : '/api/auth/login';
      const payload =
        mode === 'signup'
          ? {
              email,
              password,
              display_name: displayName,
              organization_id: organization,
              role,
            }
          : {
              email,
              password,
            };

      const res = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });

      const data = await res.json();
      if (!res.ok) {
        throw new Error(data.detail || 'Authentication failed');
      }

      const authenticatedUser = {
        ...data.user,
        id_token: data.id_token,
      };
      onSelectUser(authenticatedUser);
    } catch (err) {
      setErrorMsg(err.message);
    } finally {
      setIsSubmitting(false);
    }
  };

  const handleGenerateSignature = async () => {
    setIsSigning(true);
    setErrorMsg(null);
    try {
      const tokenToUse =
        currentUser?.id_token ||
        (currentUser?.uid ? `panoptes_session_${currentUser.uid}` : null);

      const res = await fetch('/api/auth/sign', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          id_token: tokenToUse,
          action: actionMeaning,
          meaning: 'Authenticated operator electronic signature binding action under 21 CFR Part 11 §11.200',
        }),
      });
      const data = await res.json();
      if (!res.ok) {
        throw new Error(data.detail || 'Failed to generate signature');
      }
      setGeneratedSig(data.signature);
    } catch (err) {
      setErrorMsg(err.message);
    } finally {
      setIsSigning(false);
    }
  };

  return (
    <div
      style={{
        position: 'fixed',
        inset: 0,
        backgroundColor: 'rgba(20, 20, 19, 0.42)',
        backdropFilter: 'blur(6px)',
        zIndex: 1000,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '20px',
      }}
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div
        style={{
          width: '100%',
          maxWidth: '460px',
          background: 'var(--bg-elevated)',
          borderRadius: '16px',
          padding: '32px',
          position: 'relative',
        }}
      >
        <button
          onClick={onClose}
          style={{
            position: 'absolute',
            top: '20px',
            right: '20px',
            background: 'none',
            border: 'none',
            color: 'var(--text-muted)',
            cursor: 'pointer',
          }}
          aria-label="Close modal"
        >
          <X size={20} />
        </button>

        {mode === 'profile' && currentUser?.uid ? (
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '6px' }}>
              <ShieldCheck size={22} color="#0d9488" />
              <h2 className="title-display" style={{ fontSize: '21px', margin: 0, color: 'var(--text-main)' }}>
                Operator Account
              </h2>
            </div>
            <p style={{ fontSize: '13px', color: 'var(--text-muted)', marginBottom: '20px' }}>
              Signed in to PanoptesBio Cloud Identity (panoptesbio.firebaseapp.com).
            </p>

            <div
              style={{
                background: 'var(--bg-card)',
                borderRadius: '10px',
                padding: '16px 18px',
                marginBottom: '20px',
                display: 'flex',
                flexDirection: 'column',
                gap: '8px',
                fontSize: '13px',
              }}
            >
              <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                <span style={{ color: 'var(--text-muted)' }}>Name</span>
                <strong style={{ color: 'var(--text-main)' }}>{currentUser.display_name}</strong>
              </div>
              <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                <span style={{ color: 'var(--text-muted)' }}>Email</span>
                <span style={{ color: 'var(--text-main)' }}>{currentUser.email}</span>
              </div>
              <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                <span style={{ color: 'var(--text-muted)' }}>Role</span>
                <span className="badge badge-emerald" style={{ fontSize: '11px' }}>
                  {currentUser.role}
                </span>
              </div>
            </div>

            {/* Electronic Signature Section */}
            <div
              style={{
                borderTop: '1px solid var(--border-subtle)',
                paddingTop: '18px',
                marginBottom: '20px',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '10px' }}>
                <Lock size={15} color="var(--text-main)" />
                <h4 style={{ fontSize: '13px', margin: 0, fontWeight: 600, color: 'var(--text-main)' }}>
                  21 CFR §11.200 Electronic Signature
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
                  <option value="Approval of Opentrons OT-2 Protocol">Approval of Opentrons OT-2 Protocol</option>
                  <option value="Certification of Audit Trail Integrity">Certification of Audit Trail Integrity</option>
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
                    fontSize: '12px',
                    padding: '12px 14px',
                    borderRadius: '8px',
                    background: 'var(--bg-card)',
                    color: 'var(--text-main)',
                    wordBreak: 'break-all',
                  }}
                >
                  <div style={{ color: '#0d9488', fontWeight: 600, marginBottom: '4px', display: 'flex', alignItems: 'center', gap: '6px' }}>
                    <CheckCircle2 size={14} />
                    <span>Signed by {generatedSig.operator_email}</span>
                  </div>
                  <div style={{ color: 'var(--text-muted)', fontSize: '11px' }}>
                    Timestamp: {generatedSig.signed_timestamp}
                  </div>
                  <div style={{ marginTop: '6px', color: 'var(--cyan-bright)', fontSize: '11px' }}>
                    HMAC-SHA256: {generatedSig.electronic_signature_hash}
                  </div>
                </div>
              )}
            </div>

            <button
              onClick={onSignOut}
              className="btn btn-secondary"
              style={{ width: '100%', color: '#dc2626' }}
            >
              <LogOut size={15} />
              <span>Sign Out</span>
            </button>
          </div>
        ) : (
          <div>
            <h2 className="title-display" style={{ fontSize: '24px', margin: '0 0 6px', color: 'var(--text-main)' }}>
              {mode === 'signup' ? 'Create your account' : 'Welcome back'}
            </h2>
            <p style={{ fontSize: '14px', color: 'var(--text-muted)', marginBottom: '22px' }}>
              {mode === 'signup'
                ? 'Register an operator identity to launch campaigns and sign synthesis manifests.'
                : 'Sign in with your registered PanoptesBio credentials.'}
            </p>

            {/* Mode Switcher Tabs */}
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1fr',
                gap: '6px',
                background: 'var(--bg-card)',
                padding: '4px',
                borderRadius: '8px',
                marginBottom: '20px',
              }}
            >
              <button
                type="button"
                onClick={() => {
                  setMode('signup');
                  setErrorMsg(null);
                }}
                style={{
                  padding: '8px',
                  fontSize: '13px',
                  fontWeight: mode === 'signup' ? 600 : 500,
                  background: mode === 'signup' ? 'var(--bg-elevated)' : 'transparent',
                  color: mode === 'signup' ? 'var(--text-main)' : 'var(--text-muted)',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                }}
              >
                Create Account
              </button>
              <button
                type="button"
                onClick={() => {
                  setMode('signin');
                  setErrorMsg(null);
                }}
                style={{
                  padding: '8px',
                  fontSize: '13px',
                  fontWeight: mode === 'signin' ? 600 : 500,
                  background: mode === 'signin' ? 'var(--bg-elevated)' : 'transparent',
                  color: mode === 'signin' ? 'var(--text-main)' : 'var(--text-muted)',
                  border: 'none',
                  borderRadius: '6px',
                  cursor: 'pointer',
                }}
              >
                Sign In
              </button>
            </div>

            {errorMsg && (
              <div
                style={{
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                  padding: '10px 12px',
                  borderRadius: '8px',
                  background: '#fef2f2',
                  color: '#dc2626',
                  fontSize: '13px',
                  marginBottom: '16px',
                }}
              >
                <AlertCircle size={16} />
                <span>{errorMsg}</span>
              </div>
            )}

            <form onSubmit={handleSubmitAuth} style={{ display: 'flex', flexDirection: 'column', gap: '14px' }}>
              {mode === 'signup' && (
                <div>
                  <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '5px' }}>
                    Full Name
                  </label>
                  <input
                    type="text"
                    required
                    placeholder="e.g. Dr. Elena Rostova"
                    className="input-field"
                    value={displayName}
                    onChange={(e) => setDisplayName(e.target.value)}
                  />
                </div>
              )}

              <div>
                <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '5px' }}>
                  Email Address
                </label>
                <input
                  type="email"
                  required
                  placeholder="you@laboratory.org"
                  className="input-field"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                />
              </div>

              <div>
                <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '5px' }}>
                  Password
                </label>
                <input
                  type="password"
                  required
                  minLength={6}
                  placeholder="At least 6 characters"
                  className="input-field"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
              </div>

              {mode === 'signup' && (
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '12px' }}>
                  <div>
                    <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '5px' }}>
                      Organization
                    </label>
                    <input
                      type="text"
                      className="input-field"
                      value={organization}
                      onChange={(e) => setOrganization(e.target.value)}
                    />
                  </div>
                  <div>
                    <label style={{ display: 'block', fontSize: '12px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '5px' }}>
                      Laboratory Role
                    </label>
                    <select
                      className="input-field"
                      value={role}
                      onChange={(e) => setRole(e.target.value)}
                    >
                      <option value="lead_biologist">Lead Biologist</option>
                      <option value="lab_director">Lab Director</option>
                      <option value="compliance_officer">Compliance Officer</option>
                    </select>
                  </div>
                </div>
              )}

              <button
                type="submit"
                disabled={isSubmitting}
                className="btn btn-primary"
                style={{ width: '100%', padding: '12px', fontSize: '14px', marginTop: '6px' }}
              >
                {mode === 'signup' ? <UserPlus size={16} /> : <LogIn size={16} />}
                <span>
                  {isSubmitting
                    ? 'Authenticating...'
                    : mode === 'signup'
                    ? 'Create Account'
                    : 'Sign In'}
                </span>
              </button>
            </form>
          </div>
        )}
      </div>
    </div>
  );
}
