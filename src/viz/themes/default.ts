/**
 * Default (light) theme for Carm visualizations.
 * Okabe-Ito colorblind-safe palette + clean typography + generous whitespace.
 */

/** Okabe-Ito colorblind-safe palette (8 colors). */
export const OKABE_ITO = [
  '#E69F00',  // orange
  '#56B4E9',  // sky blue
  '#009E73',  // bluish green
  '#F0E442',  // yellow
  '#0072B2',  // blue
  '#D55E00',  // vermillion
  '#CC79A7',  // reddish purple
  '#000000',  // black
] as const

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
  text: '#212529',
  textMuted: '#6c757d',
  textAnnotation: '#495057',
  gridLine: '#e9ecef',
  axisLine: '#adb5bd',
  colors: OKABE_ITO,
  fontFamily: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
  fontFamilyMono: "'JetBrains Mono', 'Fira Code', 'Cascadia Code', Menlo, Consolas, monospace",
  fontSize: 12,
  fontSizeSmall: 10,
  fontSizeTitle: 15,
  marginTop: 48,
  marginRight: 32,
  marginBottom: 72,
  marginLeft: 64,
  pointOpacity: 0.5,
  violinOpacity: 0.7,
  ciOpacity: 0.15,
}

export const DARK_THEME: CarmTheme = {
  ...DEFAULT_THEME,
  background: '#1a1a2e',
  surface: '#16213e',
  text: '#e9ecef',
  textMuted: '#adb5bd',
  textAnnotation: '#ced4da',
  gridLine: '#2d3748',
  axisLine: '#4a5568',
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
