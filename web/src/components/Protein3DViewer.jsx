import React, { useEffect, useRef, useState } from 'react';
import { Columns, Square, Orbit, Sparkles, RefreshCw } from 'lucide-react';

/**
 * Protein3DViewer — Interactive WebGL Molecular Structure Component (3Dmol.js).
 * Renders target proteins, active-site hotspots, de novo dual-chain complexes,
 * and synchronized dual-structure comparison (WT Parent vs Mutated Lead).
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
      const config = { backgroundColor: '#070a13' };
      const viewer = window.$3Dmol.createViewer(containerRef.current, config);
      viewerRef.current = viewer;

      viewer.addModel(pdbContent, 'pdb');

      if (isComplex) {
        // Dual-chain de novo complex styling
        viewer.setStyle({ chain: 'A' }, { cartoon: { color: '#64748b', opacity: 0.75 } });
        viewer.setStyle({ chain: 'B' }, { cartoon: { color: '#38bdf8', opacity: 0.95 } });

        if (hotspots && hotspots.length > 0) {
          viewer.addStyle(
            { chain: 'A', resi: hotspots },
            { stick: { color: '#fbbf24', radius: 0.28 } }
          );
        }
      } else {
        // Single scaffold styling: Electric cyan cartoon
        viewer.setStyle({}, { cartoon: { color: '#0284c7', opacity: 0.88 } });

        // Highlight catalytic / active-site residues in Amber Gold sticks
        if (hotspots && hotspots.length > 0) {
          viewer.addStyle(
            { resi: hotspots },
            { stick: { color: '#fbbf24', radius: 0.3 } }
          );
        }

        // Render cofactors / bound ligands (HETATM: PLP, 0JO, HEM) in bright rose crimson
        viewer.addStyle(
          { hetflag: true },
          {
            stick: { color: '#f43f5e', radius: 0.34 },
            sphere: { color: '#fb7185', radius: 0.45, opacity: 0.95 },
          }
        );

        // 4.5 Å Catalytic proximity coordination shell
        if (showCofactorShell) {
          viewer.addStyle(
            { within: { distance: 4.5, sel: { hetflag: true } } },
            { stick: { color: '#a855f7', radius: 0.22 } }
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

      const viewerWT = window.$3Dmol.createViewer(containerWTRef.current, { backgroundColor: '#070a13' });
      const viewerMut = window.$3Dmol.createViewer(containerMutRef.current, { backgroundColor: '#070a13' });

      viewerWTRef.current = viewerWT;
      viewerMutRef.current = viewerMut;

      // Model 1: Wild-Type Parent
      viewerWT.addModel(pdbContent, 'pdb');
      viewerWT.setStyle({}, { cartoon: { color: '#2563eb', opacity: 0.82 } });
      if (hotspots && hotspots.length > 0) {
        viewerWT.addStyle({ resi: [hotspots[0]] }, { stick: { color: '#fbbf24', radius: 0.3 } });
      }
      viewerWT.addStyle(
        { hetflag: true },
        {
          stick: { color: '#f43f5e', radius: 0.32 },
          sphere: { color: '#fb7185', radius: 0.45, opacity: 0.95 },
        }
      );
      viewerWT.zoomTo();
      viewerWT.render();

      // Model 2: Mutated Lead
      viewerMut.addModel(pdbContent, 'pdb');
      viewerMut.setStyle({}, { cartoon: { color: '#06b6d4', opacity: 0.9 } });
      if (hotspots && hotspots.length > 0) {
        viewerMut.addStyle({ resi: hotspots }, { stick: { color: '#f59e0b', radius: 0.32 } });
      }
      viewerMut.addStyle(
        { hetflag: true },
        {
          stick: { color: '#f43f5e', radius: 0.38 },
          sphere: { color: '#fb7185', radius: 0.50, opacity: 0.95 },
        }
      );
      if (showCofactorShell) {
        viewerMut.addStyle(
          { within: { distance: 4.5, sel: { hetflag: true } } },
          { stick: { color: '#a855f7', radius: 0.22 } }
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
    <div className="glass-panel" style={{ overflow: 'hidden', position: 'relative' }}>
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
            <span className="badge badge-amber">Hotspots: [{hotspots.join(', ')}]</span>
          )}
        </div>

        {/* Viewport Control Buttons */}
        <div style={{ display: 'flex', gap: '6px' }}>
          {!isComplex && (
            <button
              onClick={() => setIsDualMode(!isDualMode)}
              className="btn btn-secondary"
              style={{ padding: '4px 10px', fontSize: '11px', gap: '4px' }}
            >
              {isDualMode ? <Square size={13} /> : <Columns size={13} />}
              <span>{isDualMode ? 'Single View' : 'Dual WT vs Lead'}</span>
            </button>
          )}

          <button
            onClick={() => setShowCofactorShell(!showCofactorShell)}
            className="btn btn-secondary"
            style={{
              padding: '4px 10px',
              fontSize: '11px',
              gap: '4px',
              color: showCofactorShell ? '#c084fc' : 'var(--text-dim)',
              borderColor: showCofactorShell ? 'rgba(192, 132, 252, 0.4)' : 'var(--border-subtle)',
            }}
          >
            <Sparkles size={13} />
            <span>4.5 Å Shell</span>
          </button>

          <button
            onClick={toggleSpin}
            className="btn btn-secondary"
            style={{ padding: '4px 10px', fontSize: '11px', gap: '4px' }}
          >
            <Orbit size={13} />
            <span>{isSpinning ? 'Pause' : 'Spin'}</span>
          </button>

          <button
            onClick={centerViews}
            className="btn btn-secondary"
            style={{ padding: '4px 8px', fontSize: '11px' }}
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
          <div style={{ position: 'relative', borderRight: '1px solid rgba(255,255,255,0.08)' }}>
            <div
              style={{
                position: 'absolute',
                bottom: '10px',
                left: '12px',
                zIndex: 10,
                fontSize: '11px',
                color: '#60a5fa',
                fontWeight: 600,
                background: 'rgba(15, 23, 42, 0.8)',
                padding: '2px 8px',
                borderRadius: '4px',
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
                color: '#34d399',
                fontWeight: 600,
                background: 'rgba(15, 23, 42, 0.8)',
                padding: '2px 8px',
                borderRadius: '4px',
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
          fontFamily: 'var(--font-mono)',
          zIndex: 5,
        }}
      >
        WebGL · 3Dmol.js · Drag to rotate · Scroll to zoom
      </div>
    </div>
  );
}
