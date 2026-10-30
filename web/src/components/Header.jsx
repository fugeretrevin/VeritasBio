import React from 'react';
import {
  UserCheck,
  RefreshCw,
  LogOut,
  UserPlus,
} from 'lucide-react';

export default function Header({
  currentUser,
  onOpenAuth,
  onSignOut,
  onRefresh,
  isRefreshing,
  activePage = 'home',
  onNavigate,
}) {
  const navTabs = [
    { id: 'home', label: 'Home' },
    { id: 'studio', label: 'Studio' },
    { id: 'experiments', label: 'Experiments' },
    { id: 'hackathon', label: 'Hackathon' },
  ];

  return (
    <header
      style={{
        position: 'sticky',
        top: 0,
        left: 0,
        right: 0,
        width: '100%',
        zIndex: 950,
        background: 'rgba(249, 248, 245, 0.94)',
        backdropFilter: 'blur(12px)',
        borderBottom: '1px solid var(--border-subtle)',
      }}
    >
      <div
        style={{
          maxWidth: '1400px',
          margin: '0 auto',
          padding: '0 28px',
          height: '66px',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          gap: '20px',
        }}
      >
        {/* Brand Identity */}
        <div
          style={{ display: 'flex', alignItems: 'center', gap: '11px', cursor: 'pointer' }}
          onClick={() => onNavigate && onNavigate('home')}
        >
          <div
            style={{
              width: '28px',
              height: '28px',
              borderRadius: '6px',
              background: '#0f172a',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
            }}
            aria-hidden="true"
          >
            <span
              style={{
                width: '10px',
                height: '10px',
                borderRadius: '50%',
                background: '#10b981',
                display: 'block',
              }}
            />
          </div>

          <span
            className="title-display"
            style={{
              fontSize: '19px',
              fontWeight: 600,
              letterSpacing: '-0.03em',
              color: 'var(--text-main)',
            }}
          >
            PanoptesBio
          </span>
        </div>

        {/* Navigation Links attached directly inside top bar */}
        <nav
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: '4px',
            height: '100%',
          }}
          aria-label="Main navigation"
        >
          {navTabs.map((tab) => {
            const isActive = activePage === tab.id;
            return (
              <button
                key={tab.id}
                onClick={() => onNavigate && onNavigate(tab.id)}
                style={{
                  height: '100%',
                  padding: '0 16px',
                  fontSize: '14px',
                  fontFamily: 'var(--font-sans)',
                  fontWeight: isActive ? 600 : 500,
                  color: isActive ? '#0f172a' : 'var(--text-muted)',
                  background: 'transparent',
                  border: 'none',
                  borderBottom: isActive ? '2px solid #0f172a' : '2px solid transparent',
                  cursor: 'pointer',
                  transition: 'color 0.15s ease, border-color 0.15s ease',
                }}
              >
                {tab.label}
              </button>
            );
          })}
        </nav>

        {/* Right Controls: Sync & Real Account Authentication */}
        <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
          {activePage === 'studio' && (
            <button
              onClick={onRefresh}
              className="btn btn-secondary"
              disabled={isRefreshing}
              title="Sync Ledger State"
              style={{ padding: '7px 12px', fontSize: '13px' }}
            >
              <RefreshCw size={13} className={isRefreshing ? 'spin' : ''} />
              <span>Sync</span>
            </button>
          )}

          {currentUser?.uid ? (
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <button
                onClick={() => onOpenAuth && onOpenAuth('profile')}
                className="btn btn-secondary"
                style={{
                  padding: '7px 12px',
                  fontSize: '13px',
                  display: 'flex',
                  alignItems: 'center',
                  gap: '7px',
                }}
                title="View account profile and electronic signatures"
              >
                <UserCheck size={14} color="#059669" />
                <span style={{ fontWeight: 500, color: 'var(--text-main)' }}>
                  {currentUser.display_name || currentUser.email}
                </span>
              </button>
              <button
                onClick={onSignOut}
                className="btn btn-secondary"
                style={{ padding: '7px 10px', fontSize: '13px', color: 'var(--text-muted)' }}
                title="Sign out"
              >
                <LogOut size={14} />
              </button>
            </div>
          ) : (
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <button
                onClick={() => onOpenAuth && onOpenAuth('signin')}
                className="btn btn-secondary"
                style={{
                  padding: '7px 14px',
                  fontSize: '13px',
                  fontWeight: 500,
                }}
              >
                <span>Sign In</span>
              </button>
              <button
                onClick={() => onOpenAuth && onOpenAuth('signup')}
                className="btn btn-primary"
                style={{
                  padding: '7px 15px',
                  fontSize: '13px',
                  fontWeight: 500,
                }}
              >
                <UserPlus size={14} />
                <span>Create Account</span>
              </button>
            </div>
          )}
        </div>
      </div>
    </header>
  );
}
