// =============================================================================
// PanoptesBio — Warm Alabaster & Laboratory Color Themes
// =============================================================================

export const COLOR_THEMES = [
  {
    id: 'panoptes-white',
    name: 'Panoptes Alabaster',
    tagline: 'Warm Off-White & Soft Shading · Editorial Swiss Lab',
    primary: '#141413',
    accent: '#0d9488',
    bgPreview: '#f7f5f0',
    vars: {
      '--bg-dark': '#f7f5f0',
      '--bg-card': '#efeae1',
      '--bg-card-hover': '#e6dfd3',
      '--bg-elevated': '#fcfbf8',
      '--bg-input': '#fcfbf8',
      '--border-subtle': '#e2dcd0',
      '--border-active': '#9c9589',
      '--border-color': '#e2dcd0',
      '--border-emerald': '#a7f3d0',
      '--text-main': '#141413',
      '--text-primary': '#141413',
      '--text-muted': '#52504a',
      '--text-secondary': '#52504a',
      '--text-dim': '#78736a',
      '--cyan-primary': '#141413',
      '--cyan-bright': '#0369a1',
      '--cyan-glow': 'rgba(3, 105, 161, 0.08)',
      '--emerald-primary': '#0d9488',
      '--emerald-bright': '#0f766e',
      '--emerald-glow': 'rgba(13, 148, 136, 0.08)',
      '--amber-primary': '#141413',
      '--amber-bright': '#b45309',
      '--shadow-card': 'none',
      '--shadow-glow': 'none',
    },
  },
];

export function applyTheme(themeId) {
  const theme = COLOR_THEMES.find((t) => t.id === themeId) || COLOR_THEMES[0];
  const root = document.documentElement;

  Object.entries(theme.vars).forEach(([key, value]) => {
    root.style.setProperty(key, value);
  });

  try {
    localStorage.setItem('panoptes_theme', theme.id);
  } catch (err) {
    // Ignore localStorage errors
  }

  return theme;
}

export function getSavedThemeId() {
  return 'panoptes-white';
}
