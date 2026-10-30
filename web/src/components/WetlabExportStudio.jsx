import React, { useState, useEffect, useRef } from 'react';
import {
  Download,
  Copy,
  Check,
  FileSpreadsheet,
  Bot,
  ShieldCheck,
  Layers,
  UploadCloud,
  CheckCircle2,
  AlertCircle,
} from 'lucide-react';

export default function WetlabExportStudio({
  targetAccession = 'TRPB_TM9D8',
  currentUser,
  onRefreshLedger,
}) {
  const [workflowMode, setWorkflowMode] = useState('outbound'); // 'outbound' | 'inbound'
  const [hostOrganism, setHostOrganism] = useState('Escherichia coli (BL21-DE3)');
  const [plateFormat, setPlateFormat] = useState('96_well');
  const [isExporting, setIsExporting] = useState(false);
  const [exportData, setExportData] = useState(null);
  const [activePreviewTab, setActivePreviewTab] = useState('twist');
  const [copied, setCopied] = useState(false);

  // Inbound Plate Ingestion State
  const [plateUploading, setPlateUploading] = useState(false);
  const [plateResult, setPlateResult] = useState(null);
  const [plateError, setPlateError] = useState(null);
  const [isDragging, setIsDragging] = useState(false);
  const [normalizeToWt, setNormalizeToWt] = useState(true);
  const fileInputRef = useRef(null);

  const targetName = `${targetAccession}_Leads`;

  const fetchExports = async () => {
    setIsExporting(true);
    try {
      const params = new URLSearchParams({
        target_name: targetName,
        plate_format: plateFormat,
      });
      const res = await fetch(`/api/wetlab/export?${params.toString()}`, {
        method: 'POST',
      });
      if (res.ok) {
        const data = await res.json();
        setExportData(data);
        if (onRefreshLedger) onRefreshLedger();
      }
    } catch (err) {
      console.error('Export error:', err);
    } finally {
      setIsExporting(false);
    }
  };

  useEffect(() => {
    fetchExports();
  }, [targetAccession, plateFormat]);

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
    handleDownloadFile(exportData.twist_csv, `twist_${targetName}.csv`, 'text/csv');
    setTimeout(() => {
      handleDownloadFile(exportData.opentrons_script, `opentrons_${targetName}.py`, 'text/x-python');
    }, 200);
    setTimeout(() => {
      handleDownloadFile(exportData.fda_certificate_json, `fda_part11_${targetName}.json`, 'application/json');
    }, 400);
    setTimeout(() => {
      handleDownloadFile(exportData.benchling_manifest, `benchling_${targetName}.csv`, 'text/csv');
    }, 600);
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

  const handleCopy = () => {
    const textToCopy = getActiveContent();
    if (!textToCopy) return;
    navigator.clipboard.writeText(textToCopy);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handlePlateIngestion = async (fileOrText, isText = false) => {
    setPlateUploading(true);
    setPlateError(null);
    setPlateResult(null);

    const tokenToUse =
      currentUser?.id_token ||
      (currentUser?.uid ? `panoptes_session_${currentUser.uid}` : null);

    try {
      let res;
      if (isText) {
        res = await fetch('/api/campaign/ingest-plate', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            csv_content: fileOrText,
            target_accession: targetAccession,
            normalize_to_wt: normalizeToWt,
            auth_token: tokenToUse,
          }),
        });
      } else {
        const formData = new FormData();
        formData.append('file', fileOrText);
        formData.append('target_accession', targetAccession);
        formData.append('normalize_to_wt', normalizeToWt ? 'true' : 'false');
        formData.append('auth_token', tokenToUse);

        res = await fetch('/api/campaign/ingest-plate', {
          method: 'POST',
          body: formData,
        });
      }

      const data = await res.json();
      if (!res.ok) {
        throw new Error(data.detail || 'Plate ingestion failed');
      }

      setPlateResult(data);
      if (onRefreshLedger) onRefreshLedger();
    } catch (err) {
      setPlateError(err.message);
    } finally {
      setPlateUploading(false);
    }
  };

  const handleLoadSamplePlate = () => {
    const sampleCsv =
      targetAccession === 'TRPB_TM9D8'
        ? [
            'Well,Sample_Type,Sample_Name,Reading',
            'A01,evaporation_barrier,H2O,0.01',
            'A02,evaporation_barrier,H2O,0.01',
            'A12,evaporation_barrier,H2O,0.01',
            'B02,control_wt,WT,1.00',
            'B03,variant,V183A:F184A:V227A:S228A,0.74',
            'B04,variant,F184V:S228A,1.82',
            'B05,variant,V183I:F184V:V227I:S228A,2.15',
            'H12,evaporation_barrier,H2O,0.01',
          ].join('\n')
        : [
            'Well,Sample_Type,Sample_Name,Reading',
            'A01,evaporation_barrier,H2O,0.01',
            'A02,evaporation_barrier,H2O,0.01',
            'A12,evaporation_barrier,H2O,0.01',
            'B02,control_wt,WT,1.00',
            'B03,variant,S65T:V68L,1.85',
            'B04,variant,T203Y,1.40',
            'B05,variant,F64L:S65T:T203Y,2.45',
            'H12,evaporation_barrier,H2O,0.01',
          ].join('\n');
    handlePlateIngestion(sampleCsv, true);
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
      {/* Clean Mode Switcher & Configuration Bar */}
      <div
        className="glass-panel"
        style={{
          padding: '14px 18px',
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '12px',
        }}
      >
        <div
          style={{
            display: 'inline-flex',
            gap: '4px',
            padding: '4px',
            borderRadius: '8px',
            background: 'var(--bg-elevated)',
            border: '1px solid var(--border-subtle)',
          }}
        >
          <button
            type="button"
            onClick={() => setWorkflowMode('outbound')}
            className={`btn ${workflowMode === 'outbound' ? 'btn-primary' : 'btn-secondary'}`}
            style={{ padding: '6px 14px', fontSize: '12.5px', border: 'none' }}
          >
            Outbound Synthesis & Robot Orders
          </button>
          <button
            type="button"
            onClick={() => setWorkflowMode('inbound')}
            className={`btn ${workflowMode === 'inbound' ? 'btn-primary' : 'btn-secondary'}`}
            style={{ padding: '6px 14px', fontSize: '12.5px', border: 'none' }}
          >
            Inbound Plate Reader CSV
          </button>
        </div>

        {workflowMode === 'outbound' && (
          <div style={{ display: 'flex', alignItems: 'center', gap: '10px', flexWrap: 'wrap' }}>
            <select
              className="input-field"
              value={hostOrganism}
              onChange={(e) => setHostOrganism(e.target.value)}
              style={{ width: 'auto', fontSize: '12px', padding: '6px 10px' }}
            >
              <option value="Escherichia coli (BL21-DE3)">Host: E. coli (BL21-DE3)</option>
              <option value="Chinese Hamster Ovary (CHO-K1)">Host: CHO-K1 Mammalian</option>
              <option value="Saccharomyces cerevisiae">Host: S. cerevisiae</option>
            </select>

            <select
              className="input-field"
              value={plateFormat}
              onChange={(e) => setPlateFormat(e.target.value)}
              style={{ width: 'auto', fontSize: '12px', padding: '6px 10px' }}
            >
              <option value="96_well">96-Well Format</option>
              <option value="384_well">384-Well Format</option>
            </select>

            <button
              onClick={handleDownloadAll}
              disabled={!exportData || isExporting}
              className="btn btn-primary"
              style={{ fontSize: '12.5px', padding: '7px 14px' }}
            >
              <Download size={14} />
              <span>Download All 4 Files</span>
            </button>
          </div>
        )}
      </div>

      {workflowMode === 'outbound' ? (
        <div style={{ display: 'grid', gridTemplateColumns: '320px 1fr', gap: '16px', alignItems: 'start' }}>
          {/* Left: Compact Artifact Selector List */}
          <div className="glass-panel" style={{ padding: '16px', display: 'flex', flexDirection: 'column', gap: '10px' }}>
            {[
              {
                id: 'twist',
                title: 'Twist Bioscience cDNA CSV',
                sub: 'Codon-optimized fragments + adapters',
                icon: FileSpreadsheet,
                file: `twist_${targetName}.csv`,
                mime: 'text/csv',
                content: exportData?.twist_csv,
              },
              {
                id: 'opentrons',
                title: 'Opentrons OT-2 Python',
                sub: 'Perimeter hydration barrier + pipetting',
                icon: Bot,
                file: `opentrons_${targetName}.py`,
                mime: 'text/x-python',
                content: exportData?.opentrons_script,
              },
              {
                id: 'benchling',
                title: 'Benchling LIMS Registry',
                sub: 'Well barcodes & liability flags',
                icon: Layers,
                file: `benchling_${targetName}.csv`,
                mime: 'text/csv',
                content: exportData?.benchling_manifest,
              },
              {
                id: 'fda',
                title: 'FDA 21 CFR Part 11 Proof',
                sub: 'HMAC-SHA256 Merkle audit JSON',
                icon: ShieldCheck,
                file: `fda_part11_${targetName}.json`,
                mime: 'application/json',
                content: exportData?.fda_certificate_json,
              },
            ].map((item) => {
              const Icon = item.icon;
              const isSelected = activePreviewTab === item.id;
              return (
                <div
                  key={item.id}
                  onClick={() => setActivePreviewTab(item.id)}
                  style={{
                    padding: '12px 14px',
                    borderRadius: '8px',
                    cursor: 'pointer',
                    background: isSelected ? 'var(--bg-elevated)' : 'transparent',
                    border: isSelected ? '1px solid var(--text-main)' : '1px solid var(--border-subtle)',
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    gap: '10px',
                  }}
                >
                  <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                    <Icon size={16} color={isSelected ? '#0d9488' : 'var(--text-muted)'} />
                    <div>
                      <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)' }}>{item.title}</div>
                      <div style={{ fontSize: '11px', color: 'var(--text-dim)' }}>{item.sub}</div>
                    </div>
                  </div>
                  <button
                    type="button"
                    onClick={(e) => {
                      e.stopPropagation();
                      handleDownloadFile(item.content, item.file, item.mime);
                    }}
                    className="btn btn-secondary"
                    style={{ padding: '5px 8px', fontSize: '11px' }}
                    title={`Download ${item.file}`}
                  >
                    <Download size={13} />
                  </button>
                </div>
              );
            })}
          </div>

          {/* Right: Live Code / CSV Inspector */}
          <div className="glass-panel" style={{ padding: '18px 20px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
              <span className="mono" style={{ fontSize: '12px', fontWeight: 600, color: 'var(--text-main)' }}>
                Preview · {activePreviewTab.toUpperCase()} ({targetName})
              </span>
              <button
                onClick={handleCopy}
                disabled={!exportData}
                className="btn btn-secondary"
                style={{ fontSize: '11.5px', padding: '5px 12px', gap: '6px' }}
              >
                {copied ? <Check size={13} color="#0d9488" /> : <Copy size={13} />}
                <span>{copied ? 'Copied' : 'Copy Content'}</span>
              </button>
            </div>

            <div
              style={{
                background: 'var(--bg-elevated)',
                border: '1px solid var(--border-subtle)',
                borderRadius: '8px',
                padding: '16px',
                height: '330px',
                overflowY: 'auto',
                fontFamily: 'var(--font-mono)',
                fontSize: '12px',
                lineHeight: '1.6',
                color: 'var(--text-main)',
                whiteSpace: 'pre-wrap',
                wordBreak: 'break-all',
              }}
            >
              {isExporting ? 'Compiling manifests...' : getActiveContent() || 'No manifest loaded.'}
            </div>
          </div>
        </div>
      ) : (
        /* Inbound Microplate Reader CSV Ingestion View */
        <div className="glass-panel" style={{ padding: '28px', maxWidth: '760px' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px', flexWrap: 'wrap', gap: '10px' }}>
            <div>
              <h3 style={{ fontSize: '16px', fontWeight: 700, margin: '0 0 4px', color: 'var(--text-main)' }}>
                Ingest 96 / 384-Well Microplate Reader Data
              </h3>
              <p style={{ fontSize: '12.5px', color: 'var(--text-muted)', margin: 0 }}>
                Upload experimental optical density or fluorescence CSV returns to condition the next active learning round.
              </p>
            </div>

            <button
              onClick={handleLoadSamplePlate}
              disabled={plateUploading}
              className="btn btn-secondary"
              style={{ padding: '7px 12px', fontSize: '12px' }}
            >
              <span>Load Sample 96-Well Plate</span>
            </button>
          </div>

          <div
            onDragOver={(e) => {
              e.preventDefault();
              setIsDragging(true);
            }}
            onDragLeave={() => setIsDragging(false)}
            onDrop={(e) => {
              e.preventDefault();
              setIsDragging(false);
              if (e.dataTransfer.files?.length > 0) {
                handlePlateIngestion(e.dataTransfer.files[0], false);
              }
            }}
            onClick={() => fileInputRef.current?.click()}
            style={{
              border: `1.5px dashed ${isDragging ? '#0d9488' : 'var(--border-subtle)'}`,
              borderRadius: '10px',
              padding: '32px 20px',
              textAlign: 'center',
              cursor: 'pointer',
              background: 'var(--bg-elevated)',
              marginBottom: '16px',
            }}
          >
            <input
              ref={fileInputRef}
              type="file"
              accept=".csv,.txt"
              onChange={(e) => {
                if (e.target.files?.length > 0) {
                  handlePlateIngestion(e.target.files[0], false);
                }
              }}
              style={{ display: 'none' }}
            />
            <UploadCloud size={28} color="#0d9488" style={{ margin: '0 auto 8px' }} />
            <div style={{ fontSize: '13.5px', fontWeight: 600, color: 'var(--text-main)', marginBottom: '4px' }}>
              {plateUploading ? 'Ingesting microplate data...' : 'Drop plate reader CSV here or click to browse'}
            </div>
            <div style={{ fontSize: '11.5px', color: 'var(--text-dim)' }}>
              Automatically strips perimeter hydration barrier wells (Row A/H, Col 1/12) and normalizes to Wild-Type control
            </div>
          </div>

          <label style={{ display: 'inline-flex', alignItems: 'center', gap: '8px', fontSize: '12px', cursor: 'pointer', marginBottom: '12px' }}>
            <input
              type="checkbox"
              checked={normalizeToWt}
              onChange={(e) => setNormalizeToWt(e.target.checked)}
              style={{ accentColor: '#0d9488' }}
            />
            <span style={{ color: 'var(--text-muted)' }}>Normalize readings to Wild-Type control well (B02)</span>
          </label>

          {plateError && (
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: '8px',
                padding: '10px 14px',
                borderRadius: '8px',
                background: '#fef2f2',
                border: '1px solid #fecaca',
                color: '#dc2626',
                fontSize: '12px',
                marginTop: '8px',
              }}
            >
              <AlertCircle size={15} />
              <span>{plateError}</span>
            </div>
          )}

          {plateResult && (
            <div
              style={{
                marginTop: '12px',
                padding: '16px',
                borderRadius: '8px',
                background: 'var(--bg-elevated)',
                border: '1px solid var(--border-subtle)',
                fontSize: '12.5px',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px', color: '#0d9488', fontWeight: 600, marginBottom: '10px' }}>
                <CheckCircle2 size={16} />
                <span>Microplate Ingested · Evidence Committed to SQLite Ledger</span>
              </div>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '12px' }}>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Variants Recorded</span>
                  <strong className="mono tabular">{plateResult.variants_ingested}</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Barriers Ignored</span>
                  <strong className="mono tabular">{plateResult.barrier_wells_ignored} wells</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>WT Baseline</span>
                  <strong className="mono tabular">{plateResult.wt_baseline_score ?? 1.0}</strong>
                </div>
                <div>
                  <span style={{ color: 'var(--text-dim)', display: 'block', fontSize: '11px' }}>Top Hit</span>
                  <strong className="mono tabular" style={{ color: '#0d9488' }}>
                    {plateResult.top_variant} (+{Number(plateResult.top_fitness_discovered).toFixed(3)})
                  </strong>
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
