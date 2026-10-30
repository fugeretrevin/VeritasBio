import React, { useState, useEffect } from 'react';
import Header from './components/Header';
import HomeView from './components/HomeView';
import StudioView from './components/StudioView';
import ExperimentsView from './components/ExperimentsView';
import HackathonView from './components/HackathonView';
import AuthModal from './components/AuthModal';

export default function App() {
  const [activePage, setActivePage] = useState('home'); // 'home' | 'studio' | 'experiments' | 'hackathon'
  const [isAuthOpen, setIsAuthOpen] = useState(false);
  const [authInitialMode, setAuthInitialMode] = useState('signup'); // 'signin' | 'signup'
  const [isRefreshing, setIsRefreshing] = useState(false);

  // Real user state: starts null unless user has created an account or signed in
  const [currentUser, setCurrentUser] = useState(() => {
    try {
      const saved = localStorage.getItem('panoptes_auth_session');
      if (saved) {
        const parsed = JSON.parse(saved);
        if (parsed && parsed.uid && parsed.uid !== 'usr_panoptes_001') {
          return parsed;
        }
      }
    } catch (err) {
      // Ignore storage read errors
    }
    return null;
  });

  const [healthData, setHealthData] = useState(null);
  const [ledgerData, setLedgerData] = useState(null);

  const refreshState = async () => {
    setIsRefreshing(true);
    try {
      const [hRes, lRes] = await Promise.all([
        fetch('/api/health'),
        fetch('/api/ledger'),
      ]);
      if (hRes.ok) setHealthData(await hRes.json());
      if (lRes.ok) setLedgerData(await lRes.json());
    } catch (err) {
      console.error('Failed to sync backend state:', err);
    } finally {
      setIsRefreshing(false);
    }
  };

  useEffect(() => {
    refreshState();
    const interval = setInterval(refreshState, 12000);
    return () => clearInterval(interval);
  }, []);

  const handleOpenAuth = (mode = 'signup') => {
    setAuthInitialMode(mode);
    setIsAuthOpen(true);
  };

  const handleAuthenticatedUser = (userObj) => {
    setCurrentUser(userObj);
    try {
      if (userObj) {
        localStorage.setItem('panoptes_auth_session', JSON.stringify(userObj));
      } else {
        localStorage.removeItem('panoptes_auth_session');
      }
    } catch (err) {
      // Ignore storage errors
    }
  };

  const remBudget = healthData?.budget?.remaining_usd ?? 5.0;
  const cumSpent = healthData?.budget?.cumulative_spent_usd ?? 0.0;
  const topScore = ledgerData?.top_variants?.[0]?.computed_score ?? 1.482;

  return (
    <div style={{ minHeight: '100svh', display: 'flex', flexDirection: 'column', background: 'var(--bg-dark)' }}>
      {/* Top Header Attached Directly to Viewport Top Edge */}
      <Header
        currentUser={currentUser}
        onOpenAuth={handleOpenAuth}
        onSignOut={() => handleAuthenticatedUser(null)}
        onRefresh={refreshState}
        isRefreshing={isRefreshing}
        activePage={activePage}
        onNavigate={(page) => {
          setActivePage(page);
          window.scrollTo({ top: 0, behavior: 'smooth' });
        }}
      />

      {/* Main Content Container */}
      <div style={{ maxWidth: '1360px', width: '100%', margin: '0 auto', padding: '36px 28px 64px', flex: 1 }}>
        {/* Main Page Viewport */}
        <main>
          {activePage === 'home' && (
            <HomeView
              onNavigate={(page) => {
                setActivePage(page);
                window.scrollTo({ top: 0, behavior: 'smooth' });
              }}
              onOpenAuth={() => handleOpenAuth('signup')}
            />
          )}

          {activePage === 'studio' && (
            <StudioView
              onRefreshLedger={refreshState}
              ledgerData={ledgerData}
              currentUser={currentUser}
              remainingBudget={remBudget}
              cumulativeSpent={cumSpent}
              topFitness={topScore}
            />
          )}

          {activePage === 'experiments' && (
            <ExperimentsView />
          )}

          {activePage === 'hackathon' && (
            <HackathonView />
          )}
        </main>

        {/* Footer */}
        <footer
          style={{
            marginTop: '96px',
            paddingTop: '24px',
            borderTop: '1px solid var(--border-subtle)',
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            flexWrap: 'wrap',
            gap: '16px',
            fontSize: '13px',
            color: 'var(--text-dim)',
          }}
        >
          <div>
            PanoptesBio Platform · Deterministic Cloud-Economics Engine for Biologics
          </div>
          <div style={{ display: 'flex', gap: '20px', flexWrap: 'wrap', alignItems: 'center' }}>
            <button
              onClick={() => {
                setActivePage('home');
                window.scrollTo({ top: 0, behavior: 'smooth' });
              }}
              style={{ background: 'none', border: 'none', color: 'inherit', cursor: 'pointer', fontSize: 'inherit' }}
            >
              Home
            </button>
            <button
              onClick={() => {
                setActivePage('studio');
                window.scrollTo({ top: 0, behavior: 'smooth' });
              }}
              style={{ background: 'none', border: 'none', color: 'inherit', cursor: 'pointer', fontSize: 'inherit' }}
            >
              Studio
            </button>
            <button
              onClick={() => {
                setActivePage('experiments');
                window.scrollTo({ top: 0, behavior: 'smooth' });
              }}
              style={{ background: 'none', border: 'none', color: 'inherit', cursor: 'pointer', fontSize: 'inherit' }}
            >
              Experiments
            </button>
            <button
              onClick={() => {
                setActivePage('hackathon');
                window.scrollTo({ top: 0, behavior: 'smooth' });
              }}
              style={{ background: 'none', border: 'none', color: 'inherit', cursor: 'pointer', fontSize: 'inherit' }}
            >
              Hackathon
            </button>
          </div>
        </footer>
      </div>

      {/* Real Account Registration, Sign-In & FDA Part 11 Signature Modal */}
      <AuthModal
        isOpen={isAuthOpen}
        initialMode={authInitialMode}
        onClose={() => setIsAuthOpen(false)}
        currentUser={currentUser}
        onSelectUser={(u) => {
          handleAuthenticatedUser(u);
          setIsAuthOpen(false);
        }}
        onSignOut={() => {
          handleAuthenticatedUser(null);
          setIsAuthOpen(false);
        }}
      />
    </div>
  );
}
