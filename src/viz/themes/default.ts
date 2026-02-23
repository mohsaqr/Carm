/**
 * Default (light) theme for Carm visualizations.
 * Modern palette inspired by Tableau / Observable Plot design language.
 * Clean typography + generous whitespace + beautiful colors.
 */

/**
 * Carm modern palette — 8 perceptually distinct, visually beautiful colors.
 * Based on Tableau 10 (gold standard qualitative palette) with the pale yellow
 * replaced by dusty mauve for better legibility on white backgrounds.
 */
export const CARM_PALETTE = [
  '#4e79a7',  // cornflower blue
  '#f28e2b',  // amber orange
  '#e15759',  // soft crimson
  '#76b7b2',  // dusty teal
  '#59a14f',  // forest green
  '#af7aa1',  // dusty mauve
  '#ff9da7',  // rose pink
  '#9c755f',  // warm umber
] as const

/** @deprecated Use CARM_PALETTE. Kept for backwards compatibility. */
export const OKABE_ITO = CARM_PALETTE

export type ThemeName = 'light' | 'dark'

export interface CarmTheme {
  readonly background: string
  readonly surface: string
  readonly text: string
  readonly textMuted: string
  readonly textAnnotation: string
  readonly gridLine: string
  readonly axisLine: string
  readonly colors: readonly string[]
  readonly fontFamily: string
  readonly fontFamilyMono: string
  readonly fontSize: number
  readonly fontSizeSmall: number
  readonly fontSizeTitle: number
  readonly marginTop: number
  readonly marginRight: number
  readonly marginBottom: number
  readonly marginLeft: number
  readonly pointOpacity: number
  readonly violinOpacity: number
  readonly ciOpacity: number
}

export const DEFAULT_THEME: CarmTheme = {
  background: '#ffffff',
  surface: '#f8f9fa',
  text: '#1a1a2e',
  textMuted: '#6c757d',
  textAnnotation: '#495057',
  gridLine: '#eaeef3',
  axisLine: '#c4cdd6',
  colors: CARM_PALETTE,
  fontFamily: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
  fontFamilyMono: "'JetBrains Mono', 'Fira Code', 'Cascadia Code', Menlo, Consolas, monospace",
  fontSize: 12,
  fontSizeSmall: 11,
  fontSizeTitle: 16,
  marginTop: 58,
  marginRight: 32,
  marginBottom: 84,
  marginLeft: 64,
  pointOpacity: 0.55,
  violinOpacity: 0.72,
  ciOpacity: 0.15,
}

export const DARK_THEME: CarmTheme = {
  ...DEFAULT_THEME,
  background: '#14142b',
  surface: '#1e1e3f',
  text: '#e9ecef',
  textMuted: '#9aa5b1',
  textAnnotation: '#ced4da',
  gridLine: '#252548',
  axisLine: '#3d3d6b',
}

/** Inject CSS custom properties into a container element. */
export function applyTheme(container: HTMLElement, theme: CarmTheme = DEFAULT_THEME): void {
  const vars: Record<string, string> = {
    '--js-bg': theme.background,
    '--js-surface': theme.surface,
    '--js-text': theme.text,
    '--js-text-muted': theme.textMuted,
    '--js-text-annotation': theme.textAnnotation,
    '--js-grid': theme.gridLine,
    '--js-axis': theme.axisLine,
    '--js-font': theme.fontFamily,
    '--js-font-mono': theme.fontFamilyMono,
    '--js-font-size': `${theme.fontSize}px`,
    '--js-font-size-sm': `${theme.fontSizeSmall}px`,
    '--js-font-size-title': `${theme.fontSizeTitle}px`,
  }
  theme.colors.forEach((c, i) => { vars[`--js-color-${i}`] = c })
  for (const [k, v] of Object.entries(vars)) {
    container.style.setProperty(k, v)
  }
  container.style.background = theme.background
}

/** Get a color from the palette by index (wraps around). */
export function getColor(index: number, theme: CarmTheme = DEFAULT_THEME): string {
  return theme.colors[index % theme.colors.length]!
}

/** D3-compatible color scale from theme palette. */
export function themeColorScale(theme: CarmTheme = DEFAULT_THEME): (i: number) => string {
  return (i: number) => getColor(i, theme)
}
