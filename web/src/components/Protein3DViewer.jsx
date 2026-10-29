import React, { useEffect, useRef } from 'react';

/**
 * Protein3DViewer — Interactive WebGL Molecular Structure Component (3Dmol.js).
 * Renders target proteins, active-site hotspots, and de novo dual-chain complexes.
 */
export default function Protein3DViewer({
  pdbContent,
  targetName = 'Target Protein',
  hotspots = [],
  isComplex = false,
  height = '420px',
}) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!containerRef.current || !pdbContent) return;

    // Check if 3Dmol is loaded on window
    if (typeof window.$3Dmol === 'undefined') {
      console.warn('3Dmol.js is not yet available on window');
      return;
    }

    try {
      // Clear previous viewer instance
      containerRef.current.innerHTML = '';

      const config = { backgroundColor: '#070a13' };
      const viewer = window.$3Dmol.createViewer(containerRef.current, config);
      viewerRef.current = viewer;

      viewer.addModel(pdbContent, 'pdb');

      if (isComplex) {
        // Dual-chain de novo complex styling
        // Chain A (Target Receptor): Translucent slate cartoon
        viewer.setStyle({ chain: 'A' }, { cartoon: { color: '#64748b', opacity: 0.75 } });
        // Chain B (Designed De Novo Binder): Electric Cyan cartoon
        viewer.setStyle({ chain: 'B' }, { cartoon: { color: '#38bdf8', opacity: 0.95 } });

        // Target Hotspots: Amber Gold licorice sticks
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
      }

      viewer.zoomTo();
      viewer.render();
      viewer.spin('y', 0.4); // Subtle gentle rotation
    } catch (err) {
      console.error('Failed to initialize 3Dmol viewer:', err);
    }

    return () => {
      if (viewerRef.current) {
        viewerRef.current.spin(false);
      }
    };
  }, [pdbContent, hotspots, isComplex]);

  return (
    <div className="glass-panel" style={{ overflow: 'hidden', position: 'relative' }}>
      <div
        style={{
          position: 'absolute',
          top: '12px',
          left: '14px',
          zIndex: 10,
          display: 'flex',
          alignItems: 'center',
          gap: '8px',
        }}
      >
        <span className="badge badge-cyan">{targetName}</span>
        {isComplex ? (
          <span className="badge badge-emerald">Docked De Novo Complex (Chains A:B)</span>
        ) : (
          <span className="badge badge-amber">Hotspots: [{hotspots.join(', ')}]</span>
        )}
      </div>

      <div
        ref={containerRef}
        style={{
          width: '100%',
          height: height,
          position: 'relative',
        }}
      />

      <div
        style={{
          position: 'absolute',
          bottom: '8px',
          right: '12px',
          fontSize: '12px',
          color: 'var(--text-dim)',
          fontFamily: 'var(--font-mono)',
        }}
      >
        WebGL · 3Dmol.js · Drag to rotate · Scroll to zoom
      </div>
    </div>
  );
}
