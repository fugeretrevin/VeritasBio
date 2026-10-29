import React, { useState, useEffect } from 'react';
import {
  FlaskConical,
  Download,
  Copy,
  Check,
  FileSpreadsheet,
  Bot,
  ShieldCheck,
  Layers,
  Sparkles,
  ExternalLink,
} from 'lucide-react';

export default function WetlabExportStudio({ onRefreshLedger }) {
  const [targetName, setTargetName] = useState('TRPB_TM9D8_Winners');
  const [hostOrganism, setHostOrganism] = useState('Escherichia coli (BL21-DE3)');
  const [plateFormat, setPlateFormat] = useState('96_well');
  const [isExporting, setIsExporting] = useState(false);
  const [exportData, setExportData] = useState(null);
  const [activePreviewTab, setActivePreviewTab] = useState('twist');
  const [copied, setCopied] = useState(false);
  const [errorMessage, setErrorMessage] = useState(null);

  const fetchExports = async () => {
    setIsExporting(true);
    setErrorMessage(null);
    try {
      const params = new URLSearchParams({
        target_name: targetName,
        plate_format: plateFormat,
      });
      const res = await fetch(`/api/wetlab/export?${params.toString()}`, {
        method: 'POST',
      });
      if (!res.ok) {
        throw new Error('Failed to generate wet-lab translation manifests');
      }
      const data = await res.json();
      setExportData(data);
      if (onRefreshLedger) onRefreshLedger();
    } catch (err) {
      setErrorMessage(err.message);
    } finally {
      setIsExporting(false);
    }
  };

  useEffect(() => {
    fetchExports();
  }, [plateFormat]);

  const handleDownloadFile = (content, filename, mimeType) => {
    if (!content) return;
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  };

  const handleDownloadAll = () => {
    if (!exportData) return;
    handleDownloadFile(
      exportData.twist_csv,
      `twist_bioscience_${targetName}_orders.csv`,
      'text/csv'
    );
    setTimeout(() => {
      handleDownloadFile(
        exportData.opentrons_script,
        `opentrons_${targetName}_protocol.py`,
        'text/x-python'
      );
    }, 250);
    setTimeout(() => {
      handleDownloadFile(
        exportData.fda_certificate_json,
        `fda_part11_${targetName}_certificate.json`,
        'application/json'
      );
    }, 500);
    setTimeout(() => {
      handleDownloadFile(
        exportData.benchling_manifest,
        `benchling_lims_${targetName}_manifest.csv`,
        'text/csv'
      );
    }, 750);
  };

  const handleCopy = () => {
    const textToCopy = getActiveContent();
    if (!textToCopy) return;
    navigator.clipboard.writeText(textToCopy);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const getActiveContent = () => {
    if (!exportData) return '';
    switch (activePreviewTab) {
      case 'twist':
        return exportData.twist_csv || '';
      case 'opentrons':
        return exportData.opentrons_script || '';
      case 'fda':
        return exportData.fda_certificate_json || '';
      case 'benchling':
        return exportData.benchling_manifest || '';
      case 'plate':
        return exportData.plate_csv || '';
      default:
        return '';
    }
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '24px' }}>
      {/* Top Header Card */}
      <div
        className="glass-panel"
        style={{
          padding: '24px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '16px',
          borderLeft: '4px solid var(--cyan-bright)',
        }}
      >
        <div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '6px' }}>
            <FlaskConical size={22} color="var(--cyan-bright)" />
            <h2 className="title-display" style={{ fontSize: '20px', margin: 0 }}>
              Physical Wet-Lab & CRO Translation Studio
            </h2>
          </div>
          <p style={{ fontSize: '13px', color: 'var(--text-muted)', margin: 0, maxWidth: '800px' }}>
            Bridge computational leads directly into live molecular biology: automated Twist Bioscience cDNA synthesis orders,
            Opentrons OT-2 liquid-handling protocols, Benchling LIMS manifests, and FDA 21 CFR Part 11 cryptographic audit certificates.
          </p>
        </div>

        <div style={{ display: 'flex', gap: '10px' }}>
          <button
            onClick={fetchExports}
            disabled={isExporting}
            className="btn btn-secondary"
            style={{ fontSize: '13px', padding: '10px 16px' }}
          >
            <Sparkles size={16} />
            <span>{isExporting ? 'Compiling...' : 'Regenerate Manifests'}</span>
          </button>
          <button
            onClick={handleDownloadAll}
            disabled={!exportData || isExporting}
            className="btn btn-primary"
            style={{ fontSize: '13px', padding: '10px 18px' }}
          >
            <Download size={16} />
            <span>Download All Artifacts</span>
          </button>
        </div>
      </div>

      {/* Configuration Controls Bar */}
      <div
        className="glass-panel"
        style={{
          padding: '16px 20px',
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))',
          gap: '16px',
          alignItems: 'center',
        }}
      >
        <div>
          <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>
            Target Campaign Name
          </label>
          <input
            type="text"
            className="input-field"
            value={targetName}
            onChange={(e) => setTargetName(e.target.value)}
            style={{ fontSize: '12px', padding: '8px 12px' }}
          />
        </div>

        <div>
          <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>
            Synthesis Expression Host
          </label>
          <select
            className="input-field"
            value={hostOrganism}
            onChange={(e) => setHostOrganism(e.target.value)}
            style={{ fontSize: '12px', padding: '8px 12px' }}
          >
            <option value="Escherichia coli (BL21-DE3)">Escherichia coli (BL21-DE3 / K12)</option>
            <option value="Chinese Hamster Ovary (CHO-K1)">Chinese Hamster Ovary (CHO-K1)</option>
            <option value="Saccharomyces cerevisiae">Saccharomyces cerevisiae (Baker's Yeast)</option>
            <option value="Pichia pastoris">Pichia pastoris (Komagataella phaffii)</option>
          </select>
        </div>

        <div>
          <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-dim)', marginBottom: '4px' }}>
            Microplate Screening Format
          </label>
          <select
            className="input-field"
            value={plateFormat}
            onChange={(e) => setPlateFormat(e.target.value)}
            style={{ fontSize: '12px', padding: '8px 12px' }}
          >
            <option value="96_well">96-Well Plate (56 leads + 36 hydration barriers)</option>
            <option value="384_well">384-Well Plate (308 leads + 76 hydration barriers)</option>
          </select>
        </div>
      </div>

      {/* 4 Cards Grid */}
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(280px, 1fr))', gap: '16px' }}>
        {/* Card 1: Twist Bioscience */}
        <div className="glass-panel" style={{ padding: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
              <FileSpreadsheet size={20} color="var(--cyan-bright)" />
              <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
                Twist Bioscience cDNA Orders
              </h3>
            </div>
            <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '16px' }}>
              Codon-optimized gene fragments with Gibson assembly overhangs, 5' flanking primers, and non-clashing synthesis adapters.
            </p>
          </div>
          <button
            onClick={() =>
              handleDownloadFile(
                exportData?.twist_csv,
                `twist_bioscience_${targetName}_orders.csv`,
                'text/csv'
              )
            }
            disabled={!exportData}
            className="btn btn-secondary"
            style={{ width: '100%', fontSize: '12px', gap: '8px', padding: '10px' }}
          >
            <Download size={14} />
            <span>Download Twist CSV</span>
          </button>
        </div>

        {/* Card 2: Opentrons OT-2 Robot */}
        <div className="glass-panel" style={{ padding: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
              <Bot size={20} color="var(--emerald-bright)" />
              <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
                Opentrons OT-2 Protocol
              </h3>
            </div>
            <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '16px' }}>
              Python liquid-handling robot script with perimeter anti-evaporation hydration barriers, 8-channel multi-pipetting, and blow-out logic.
            </p>
          </div>
          <button
            onClick={() =>
              handleDownloadFile(
                exportData?.opentrons_script,
                `opentrons_${targetName}_protocol.py`,
                'text/x-python'
              )
            }
            disabled={!exportData}
            className="btn btn-secondary"
            style={{ width: '100%', fontSize: '12px', gap: '8px', padding: '10px' }}
          >
            <Download size={14} />
            <span>Download Opentrons Script</span>
          </button>
        </div>

        {/* Card 3: Benchling LIMS Manifest */}
        <div className="glass-panel" style={{ padding: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
              <Layers size={20} color="#a855f7" />
              <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
                Benchling LIMS Registry
              </h3>
            </div>
            <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '16px' }}>
              Bulk import manifest with well barcodes, coordinates, concentrations, and developability liability flags for enterprise inventory.
            </p>
          </div>
          <button
            onClick={() =>
              handleDownloadFile(
                exportData?.benchling_manifest,
                `benchling_lims_${targetName}_manifest.csv`,
                'text/csv'
              )
            }
            disabled={!exportData}
            className="btn btn-secondary"
            style={{ width: '100%', fontSize: '12px', gap: '8px', padding: '10px' }}
          >
            <Download size={14} />
            <span>Download Benchling CSV</span>
          </button>
        </div>

        {/* Card 4: FDA 21 CFR Part 11 Audit Certificate */}
        <div className="glass-panel" style={{ padding: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
          <div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
              <ShieldCheck size={20} color="#f59e0b" />
              <h3 className="title-display" style={{ fontSize: '15px', margin: 0 }}>
                FDA 21 CFR Part 11 Certificate
              </h3>
            </div>
            <p style={{ fontSize: '12px', color: 'var(--text-muted)', marginBottom: '16px' }}>
              Cryptographically signed JSON certificate containing dual-identifier electronic signatures, HMAC-SHA256 Merkle proofs, and timestamps.
            </p>
          </div>
          <button
            onClick={() =>
              handleDownloadFile(
                exportData?.fda_certificate_json,
                `fda_part11_${targetName}_certificate.json`,
                'application/json'
              )
            }
            disabled={!exportData}
            className="btn btn-secondary"
            style={{ width: '100%', fontSize: '12px', gap: '8px', padding: '10px' }}
          >
            <Download size={14} />
            <span>Download FDA Certificate</span>
          </button>
        </div>
      </div>

      {/* Live Preview Window */}
      <div className="glass-panel" style={{ padding: '24px' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px', flexWrap: 'wrap', gap: '10px' }}>
          {/* Sub Tabs */}
          <div style={{ display: 'flex', gap: '8px', background: 'rgba(15, 23, 42, 0.6)', padding: '4px', borderRadius: '8px' }}>
            <button
              onClick={() => setActivePreviewTab('twist')}
              className={`btn ${activePreviewTab === 'twist' ? 'btn-primary' : 'btn-secondary'}`}
              style={{ padding: '6px 14px', fontSize: '12px' }}
            >
              Twist CSV
            </button>
            <button
              onClick={() => setActivePreviewTab('opentrons')}
              className={`btn ${activePreviewTab === 'opentrons' ? 'btn-primary' : 'btn-secondary'}`}
              style={{ padding: '6px 14px', fontSize: '12px' }}
            >
              Opentrons Python
            </button>
            <button
              onClick={() => setActivePreviewTab('benchling')}
              className={`btn ${activePreviewTab === 'benchling' ? 'btn-primary' : 'btn-secondary'}`}
              style={{ padding: '6px 14px', fontSize: '12px' }}
            >
              Benchling LIMS
            </button>
            <button
              onClick={() => setActivePreviewTab('fda')}
              className={`btn ${activePreviewTab === 'fda' ? 'btn-primary' : 'btn-secondary'}`}
              style={{ padding: '6px 14px', fontSize: '12px' }}
            >
              FDA Certificate JSON
            </button>
            <button
              onClick={() => setActivePreviewTab('plate')}
              className={`btn ${activePreviewTab === 'plate' ? 'btn-primary' : 'btn-secondary'}`}
              style={{ padding: '6px 14px', fontSize: '12px' }}
            >
              Plate Map CSV
            </button>
          </div>

          <button
            onClick={handleCopy}
            disabled={!exportData}
            className="btn btn-secondary"
            style={{ fontSize: '12px', padding: '6px 14px', gap: '6px' }}
          >
            {copied ? <Check size={14} color="var(--emerald-bright)" /> : <Copy size={14} />}
            <span>{copied ? 'Copied to Clipboard' : 'Copy File Content'}</span>
          </button>
        </div>

        {/* Text Viewport */}
        <div
          style={{
            background: '#070a13',
            border: '1px solid rgba(255,255,255,0.08)',
            borderRadius: '8px',
            padding: '16px',
            maxHeight: '440px',
            overflowY: 'auto',
            fontFamily: 'var(--font-mono)',
            fontSize: '12px',
            lineHeight: '1.6',
            color: '#e2e8f0',
            whiteSpace: 'pre-wrap',
            wordBreak: 'break-all',
          }}
        >
          {isExporting ? (
            <div style={{ color: 'var(--text-dim)', textAlign: 'center', padding: '40px' }}>
              Compiling physical synthesis sequences & robotic trajectories...
            </div>
          ) : getActiveContent() ? (
            getActiveContent()
          ) : (
            <div style={{ color: 'var(--text-dim)', textAlign: 'center', padding: '40px' }}>
              No manifest data available. Click "Regenerate Manifests" above.
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

