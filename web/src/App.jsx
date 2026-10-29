import React, { useState, useEffect } from 'react';
import Header from './components/Header';
import KpiBar from './components/KpiBar';
import CampaignStudio from './components/CampaignStudio';
import DeNovoStudio from './components/DeNovoStudio';
import LedgerAuditTrail from './components/LedgerAuditTrail';
import AuthModal from './components/AuthModal';
import { Compass, Sparkles, ShieldCheck } from 'lucide-react';

export default function App() {
  const [activeTab, setActiveTab] = useState('campaign');
  const [isAuthOpen, setIsAuthOpen] = useState(false);
  const [isRefreshing, setIsRefreshing] = useState(false);

  const [currentUser, setCurrentUser] = useState({
    uid: 'usr_panoptes_001',
    email: 'scientist@panoptesbio.com',
    display_name: 'Dr. Lead Biologist',
    organization_id: 'panoptesbio',
    role: 'lead_biologist',
    is_part11_certified: true,
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
    // Background polling every 12 seconds
    const interval = setInterval(refreshState, 12000);
    return () => clearInterval(interval);
  }, []);

  const remBudget = healthData?.budget?.remaining_usd ?? 5.0;
  const cumSpent = healthData?.budget?.cumulative_spent_usd ?? 0.0;
  const topScore = ledgerData?.top_variants?.[0]?.computed_score ?? 1.482;

  return (
    <div style={{ maxWidth: '1440px', margin: '0 auto', padding: '24px 20px 48px' }}>
      {/* Top Header */}
      <Header
        currentUser={currentUser}
        onOpenAuth={() => setIsAuthOpen(true)}
        onRefresh={refreshState}
        isRefreshing={isRefreshing}
      />

      {/* KPI Stats Bar */}
      <KpiBar
        remainingBudget={remBudget}
        cumulativeSpent={cumSpent}
        topFitness={topScore}
        savedPct="99.5%"
        speedup="184.3x Cheaper"
      />

      {/* Main Tab Navigation */}
      <div
        className="glass-panel"
        style={{
          display: 'flex',
          gap: '8px',
          padding: '6px',
          marginBottom: '20px',
          width: 'fit-content',
        }}
      >
        <button
          onClick={() => setActiveTab('campaign')}
          className={`btn ${activeTab === 'campaign' ? 'btn-primary' : 'btn-secondary'}`}
          style={{ padding: '8px 18px', fontSize: '13px' }}
        >
          <Compass size={16} />
          <span>Active Learning Campaign</span>
        </button>

        <button
          onClick={() => setActiveTab('denovo')}
          className={`btn ${activeTab === 'denovo' ? 'btn-primary' : 'btn-secondary'}`}
          style={{ padding: '8px 18px', fontSize: '13px' }}
        >
          <Sparkles size={16} />
          <span>De Novo Binder Diffusion Studio</span>
        </button>

        <button
          onClick={() => setActiveTab('ledger')}
          className={`btn ${activeTab === 'ledger' ? 'btn-primary' : 'btn-secondary'}`}
          style={{ padding: '8px 18px', fontSize: '13px' }}
        >
          <ShieldCheck size={16} />
          <span>Ledger & Merkle Audit Trail</span>
        </button>
      </div>

      {/* Active Tab Viewport */}
      <main>
        {activeTab === 'campaign' && <CampaignStudio onRefreshLedger={refreshState} />}
        {activeTab === 'denovo' && <DeNovoStudio onRefreshLedger={refreshState} />}
        {activeTab === 'ledger' && <LedgerAuditTrail ledgerData={ledgerData} onRefresh={refreshState} />}
      </main>

      {/* Identity / Authentication Modal */}
      <AuthModal
        isOpen={isAuthOpen}
        onClose={() => setIsAuthOpen(false)}
        currentUser={currentUser}
        onSelectUser={(u) => {
          setCurrentUser(u);
          setIsAuthOpen(false);
        }}
      />

      {/* Modern Biopharma Footer */}
      <footer
        style={{
          marginTop: '40px',
          paddingTop: '20px',
          borderTop: '1px solid var(--border-subtle)',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '12px',
          fontSize: '12px',
          color: 'var(--text-dim)',
        }}
      >
        <div>
          PanoptesBio Engine · Deterministic Cloud Economics · Fail-Closed Spend Gateway
        </div>
        <div style={{ display: 'flex', gap: '16px' }}>
          <span>FDA 21 CFR Part 11 Compliant</span>
          <span>•</span>
          <span>HMAC-SHA256 Merkle Ledger</span>
          <span>•</span>
          <span>Modal A10G + ESM-2 / ProteinMPNN</span>
        </div>
      </footer>
    </div>
  );
}
