import React, { useState, useEffect, useRef } from 'react';
import { Palette, Check, ChevronDown } from 'lucide-react';
import { COLOR_THEMES, applyTheme, getSavedThemeId } from '../themes';

export default function ThemePicker() {
  const [isOpen, setIsOpen] = useState(false);
  const [activeThemeId, setActiveThemeId] = useState(getSavedThemeId);
  const dropdownRef = useRef(null);

  useEffect(() => {
    // Apply saved theme on initial mount
    const saved = getSavedThemeId();
    applyTheme(saved);
    setActiveThemeId(saved);
  }, []);

  useEffect(() => {
    // Close on click outside or Escape
    function handleClickOutside(event) {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
        setIsOpen(false);
      }
    }
    function handleKeyDown(event) {
      if (event.key === 'Escape') {
        setIsOpen(false);
      }
    }

    if (isOpen) {
      document.addEventListener('mousedown', handleClickOutside);
      document.addEventListener('keydown', handleKeyDown);
    }
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [isOpen]);

  const handleSelectTheme = (themeId) => {
    applyTheme(themeId);
    setActiveThemeId(themeId);
    setIsOpen(false);
  };

  const currentTheme = COLOR_THEMES.find((t) => t.id === activeThemeId) || COLOR_THEMES[0];

  return (
    <div style={{ position: 'relative' }} ref={dropdownRef}>
      {/* Trigger Button in Header */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="btn btn-secondary"
        style={{
          padding: '6px 12px',
          display: 'flex',
          alignItems: 'center',
          gap: '8px',
          fontSize: '12px',
          borderRadius: 'var(--radius-full)',
        }}
        title="Select Biophysical Theme (12 Themes Available)"
        aria-expanded={isOpen}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '3px' }}>
          <span
            style={{
              width: '9px',
              height: '9px',
              borderRadius: '50%',
              backgroundColor: currentTheme.primary,
              boxShadow: `0 0 6px ${currentTheme.primary}`,
            }}
          />
          <span
            style={{
              width: '9px',
              height: '9px',
              borderRadius: '50%',
              backgroundColor: currentTheme.accent,
            }}
          />
        </div>
        <Palette size={14} color="var(--cyan-bright)" />
        <span style={{ fontWeight: 500 }}>{currentTheme.name}</span>
        <ChevronDown size={13} style={{ transform: isOpen ? 'rotate(180deg)' : 'none', transition: 'transform 0.15s' }} />
      </button>

      {/* Popover Menu */}
      {isOpen && (
        <div
          className="glass-panel"
          style={{
            position: 'absolute',
            top: 'calc(100% + 8px)',
            right: 0,
            width: '320px',
            maxHeight: '440px',
            overflowY: 'auto',
            padding: '12px',
            borderRadius: '12px',
            background: 'var(--bg-card)',
            border: '1px solid var(--border-active)',
            boxShadow: '0 10px 25px rgba(35, 30, 20, 0.12), 0 1px 3px rgba(35, 30, 20, 0.05)',
            zIndex: 9999,
          }}
        >
          <div style={{ padding: '4px 8px 10px', borderBottom: '1px solid var(--border-subtle)', marginBottom: '8px' }}>
            <div style={{ fontSize: '13px', fontWeight: 600, color: 'var(--text-main)' }}>
              Biophysical Color Themes
            </div>
            <div style={{ fontSize: '11px', color: 'var(--text-dim)', marginTop: '2px' }}>
              12 laboratory-calibrated palettes (WCAG AA compliant)
            </div>
          </div>

          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            {COLOR_THEMES.map((theme) => {
              const isSelected = theme.id === activeThemeId;
              return (
                <button
                  key={theme.id}
                  onClick={() => handleSelectTheme(theme.id)}
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    padding: '8px 10px',
                    borderRadius: '8px',
                    background: isSelected ? 'var(--bg-elevated)' : 'transparent',
                    border: isSelected ? '1px solid var(--border-active)' : '1px solid transparent',
                    cursor: 'pointer',
                    textAlign: 'left',
                    transition: 'all 0.15s ease',
                  }}
                  onMouseEnter={(e) => {
                    if (!isSelected) e.currentTarget.style.background = 'var(--bg-card-hover)';
                  }}
                  onMouseLeave={(e) => {
                    if (!isSelected) e.currentTarget.style.background = 'transparent';
                  }}
                >
                  <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                    {/* Swatch Preview Pill */}
                    <div
                      style={{
                        width: '26px',
                        height: '26px',
                        borderRadius: '6px',
                        backgroundColor: theme.bgPreview,
                        border: '1px solid var(--border-subtle)',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        gap: '2px',
                        padding: '2px',
                      }}
                    >
                      <span
                        style={{
                          width: '7px',
                          height: '7px',
                          borderRadius: '50%',
                          backgroundColor: theme.primary,
                        }}
                      />
                      <span
                        style={{
                          width: '7px',
                          height: '7px',
                          borderRadius: '50%',
                          backgroundColor: theme.accent,
                        }}
                      />
                    </div>

                    <div>
                      <div style={{ fontSize: '12px', fontWeight: isSelected ? 600 : 500, color: 'var(--text-main)' }}>
                        {theme.name}
                      </div>
                      <div style={{ fontSize: '10px', color: 'var(--text-dim)', marginTop: '1px' }}>
                        {theme.tagline}
                      </div>
                    </div>
                  </div>

                  {isSelected && <Check size={15} color="var(--cyan-primary)" />}
                </button>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}

