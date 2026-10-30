import React, { useEffect, useRef, useState } from 'react';
import { Columns, Square, Orbit, Sparkles, RefreshCw } from 'lucide-react';

/**
 * Protein3DViewer — Interactive WebGL Molecular Structure Component (3Dmol.js).
 * Renders target proteins, active-site hotspots, de novo dual-chain complexes,
 * and synchronized dual-structure comparison (WT Parent vs Mutated Lead) on a cohesive light canvas.
 */
export default function Protein3DViewer({
  pdbContent,
  targetName = 'Target Protein',
  hotspots = [],
  isComplex = false,
  height = '420px',
}) {
  const containerRef = useRef(null);
  const containerWTRef = useRef(null);
  const containerMutRef = useRef(null);
  const viewerRef = useRef(null);
  const viewerWTRef = useRef(null);
  const viewerMutRef = useRef(null);

  const [isDualMode, setIsDualMode] = useState(false);
  const [showCofactorShell, setShowCofactorShell] = useState(true);
  const [isSpinning, setIsSpinning] = useState(true);

  // Initialize Single Viewport
  useEffect(() => {
    if (isDualMode || !containerRef.current || !pdbContent) return;
    if (typeof window.$3Dmol === 'undefined') {
      console.warn('3Dmol.js is not yet available on window');
      return;
    }

    try {
      containerRef.current.innerHTML = '';
      const config = { backgroundColor: '#fdfcfa', antialias: true };
      const viewer = window.$3Dmol.createViewer(containerRef.current, config);
      viewerRef.current = viewer;

      viewer.addModel(pdbContent, 'pdb');

      if (isComplex) {
        viewer.setStyle({ chain: 'A' }, { cartoon: { color: '#64748b', opacity: 0.82 } });
        viewer.setStyle({ chain: 'B' }, { cartoon: { color: '#0284c7', opacity: 0.95 } });

        if (hotspots && hotspots.length > 0) {
          viewer.addStyle(
            { chain: 'A', resi: hotspots },
            { stick: { color: '#059669', radius: 0.3 } }
          );
        }
      } else {
        viewer.setStyle({}, { cartoon: { color: '#334155', opacity: 0.9 } });

        if (hotspots && hotspots.length > 0) {
          viewer.addStyle(
            { resi: hotspots },
            {
              stick: { color: '#059669', radius: 0.34 },
              sphere: { color: '#10b981', radius: 0.38, opacity: 0.88 },
            }
          );
        }

        viewer.addStyle(
          { hetflag: true },
          {
            stick: { color: '#0284c7', radius: 0.34 },
            sphere: { color: '#0284c7', radius: 0.45, opacity: 0.92 },
          }
        );

        if (showCofactorShell) {
          viewer.addStyle(
            { within: { distance: 4.5, sel: { hetflag: true } } },
            { stick: { color: '#0f172a', radius: 0.22 } }
          );
        }
      }

      viewer.zoomTo();
      viewer.render();
      if (isSpinning) viewer.spin('y', 0.4);
    } catch (err) {
      console.error('Failed to initialize 3Dmol viewer:', err);
    }

    return () => {
      if (viewerRef.current) viewerRef.current.spin(false);
    };
  }, [pdbContent, hotspots, isComplex, isDualMode, showCofactorShell]);

  // Initialize Dual Viewports
  useEffect(() => {
    if (!isDualMode || !containerWTRef.current || !containerMutRef.current || !pdbContent) return;
    if (typeof window.$3Dmol === 'undefined') return;

    try {
      containerWTRef.current.innerHTML = '';
      containerMutRef.current.innerHTML = '';

      const viewerWT = window.$3Dmol.createViewer(containerWTRef.current, { backgroundColor: '#fdfcfa', antialias: true });
      const viewerMut = window.$3Dmol.createViewer(containerMutRef.current, { backgroundColor: '#fdfcfa', antialias: true });

      viewerWTRef.current = viewerWT;
      viewerMutRef.current = viewerMut;

      // Model 1: Wild-Type Parent
      viewerWT.addModel(pdbContent, 'pdb');
      viewerWT.setStyle({}, { cartoon: { color: '#64748b', opacity: 0.85 } });
      if (hotspots && hotspots.length > 0) {
        viewerWT.addStyle({ resi: [hotspots[0]] }, { stick: { color: '#0284c7', radius: 0.3 } });
      }
      viewerWT.addStyle(
        { hetflag: true },
        {
          stick: { color: '#0f172a', radius: 0.32 },
          sphere: { color: '#0284c7', radius: 0.45, opacity: 0.9 },
        }
      );
      viewerWT.zoomTo();
      viewerWT.render();

      // Model 2: Mutated Lead
      viewerMut.addModel(pdbContent, 'pdb');
      viewerMut.setStyle({}, { cartoon: { color: '#0f172a', opacity: 0.9 } });
      if (hotspots && hotspots.length > 0) {
        viewerMut.addStyle(
          { resi: hotspots },
          {
            stick: { color: '#059669', radius: 0.35 },
            sphere: { color: '#10b981', radius: 0.4, opacity: 0.9 },
          }
        );
      }
      viewerMut.addStyle(
        { hetflag: true },
        {
          stick: { color: '#0284c7', radius: 0.38 },
          sphere: { color: '#0284c7', radius: 0.48, opacity: 0.92 },
        }
      );
      if (showCofactorShell) {
        viewerMut.addStyle(
          { within: { distance: 4.5, sel: { hetflag: true } } },
          { stick: { color: '#475569', radius: 0.22 } }
        );
      }
      viewerMut.zoomTo();
      viewerMut.render();

      if (isSpinning) {
        viewerWT.spin('y', 0.4);
        viewerMut.spin('y', 0.4);
      }
    } catch (err) {
      console.error('Failed to initialize dual 3Dmol viewers:', err);
    }

    return () => {
      if (viewerWTRef.current) viewerWTRef.current.spin(false);
      if (viewerMutRef.current) viewerMutRef.current.spin(false);
    };
  }, [pdbContent, hotspots, isDualMode, showCofactorShell]);

  const toggleSpin = () => {
    const nextSpin = !isSpinning;
    setIsSpinning(nextSpin);
    if (viewerRef.current) viewerRef.current.spin(nextSpin ? 'y' : false, 0.4);
    if (viewerWTRef.current) viewerWTRef.current.spin(nextSpin ? 'y' : false, 0.4);
    if (viewerMutRef.current) viewerMutRef.current.spin(nextSpin ? 'y' : false, 0.4);
  };

  const centerViews = () => {
    if (viewerRef.current) {
      viewerRef.current.zoomTo();
      viewerRef.current.render();
    }
    if (viewerWTRef.current) {
      viewerWTRef.current.zoomTo();
      viewerWTRef.current.render();
    }
    if (viewerMutRef.current) {
      viewerMutRef.current.zoomTo();
      viewerMutRef.current.render();
    }
  };

  return (
    <div
      className="glass-panel"
      style={{
        overflow: 'hidden',
        position: 'relative',
        background: '#ffffff',
      }}
    >
      {/* Top Floating Control Bar */}
      <div
        style={{
          position: 'absolute',
          top: '12px',
          left: '14px',
          right: '14px',
          zIndex: 10,
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          flexWrap: 'wrap',
          gap: '8px',
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
          <span className="badge badge-cyan">{targetName}</span>
          {isComplex ? (
            <span className="badge badge-emerald">Docked De Novo Complex (Chains A:B)</span>
          ) : (
            <span className="badge badge-secondary">Hotspots: [{hotspots.join(', ')}]</span>
          )}
        </div>

        {/* Viewport Control Buttons */}
        <div style={{ display: 'flex', gap: '6px' }}>
          {!isComplex && (
            <button
              onClick={() => setIsDualMode(!isDualMode)}
              className="btn btn-secondary"
              style={{ padding: '5px 10px', fontSize: '12px', gap: '5px' }}
            >
              {isDualMode ? <Square size={13} /> : <Columns size={13} />}
              <span>{isDualMode ? 'Single View' : 'Dual WT vs Lead'}</span>
            </button>
          )}

          <button
            onClick={() => setShowCofactorShell(!showCofactorShell)}
            className="btn btn-secondary"
            style={{
              padding: '5px 10px',
              fontSize: '12px',
              gap: '5px',
              color: showCofactorShell ? '#0284c7' : 'var(--text-dim)',
              borderColor: showCofactorShell ? '#bae6fd' : 'var(--border-subtle)',
            }}
          >
            <Sparkles size={13} />
            <span>4.5 Å Shell</span>
          </button>

          <button
            onClick={toggleSpin}
            className="btn btn-secondary"
            style={{ padding: '5px 10px', fontSize: '12px', gap: '5px' }}
          >
            <Orbit size={13} />
            <span>{isSpinning ? 'Pause' : 'Spin'}</span>
          </button>

          <button
            onClick={centerViews}
            className="btn btn-secondary"
            style={{ padding: '5px 8px', fontSize: '12px' }}
            title="Center Viewport"
          >
            <RefreshCw size={13} />
          </button>
        </div>
      </div>

      {/* Main Viewport Container */}
      {!isDualMode ? (
        <div
          ref={containerRef}
          style={{
            width: '100%',
            height: height,
            position: 'relative',
          }}
        />
      ) : (
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', width: '100%', height: height }}>
          {/* Left: WT Parent */}
          <div style={{ position: 'relative', borderRight: '1px solid var(--border-subtle)' }}>
            <div
              style={{
                position: 'absolute',
                bottom: '10px',
                left: '12px',
                zIndex: 10,
                fontSize: '11px',
                color: '#475569',
                fontWeight: 600,
                background: 'rgba(255, 255, 255, 0.92)',
                border: '1px solid var(--border-subtle)',
                padding: '3px 9px',
                borderRadius: '6px',
              }}
            >
              ● Wild-Type Parent (Baseline)
            </div>
            <div ref={containerWTRef} style={{ width: '100%', height: '100%' }} />
          </div>

          {/* Right: Mutated Lead */}
          <div style={{ position: 'relative' }}>
            <div
              style={{
                position: 'absolute',
                bottom: '10px',
                left: '12px',
                zIndex: 10,
                fontSize: '11px',
                color: '#059669',
                fontWeight: 600,
                background: 'rgba(255, 255, 255, 0.92)',
                border: '1px solid var(--border-subtle)',
                padding: '3px 9px',
                borderRadius: '6px',
              }}
            >
              ★ Engineered Lead Candidate
            </div>
            <div ref={containerMutRef} style={{ width: '100%', height: '100%' }} />
          </div>
        </div>
      )}

      {/* Bottom Hint */}
      <div
        style={{
          position: 'absolute',
          bottom: '8px',
          right: '12px',
          fontSize: '11px',
          color: 'var(--text-dim)',
          zIndex: 5,
        }}
      >
        WebGL · 3Dmol.js · Drag to rotate · Scroll to zoom
      </div>
    </div>
  );
}
