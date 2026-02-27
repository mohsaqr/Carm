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
export declare const CARM_PALETTE: readonly ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#af7aa1", "#ff9da7", "#9c755f"];
/** @deprecated Use CARM_PALETTE. Kept for backwards compatibility. */
export declare const OKABE_ITO: readonly ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#af7aa1", "#ff9da7", "#9c755f"];
export type ThemeName = 'light' | 'dark';
export interface CarmTheme {
    readonly background: string;
    readonly surface: string;
    readonly text: string;
    readonly textMuted: string;
    readonly textAnnotation: string;
    readonly gridLine: string;
    readonly axisLine: string;
    readonly colors: readonly string[];
    readonly fontFamily: string;
    readonly fontFamilyMono: string;
    readonly fontSize: number;
    readonly fontSizeSmall: number;
    readonly fontSizeTitle: number;
    readonly marginTop: number;
    readonly marginRight: number;
    readonly marginBottom: number;
    readonly marginLeft: number;
    readonly pointOpacity: number;
    readonly violinOpacity: number;
    readonly ciOpacity: number;
}
export declare const DEFAULT_THEME: CarmTheme;
export declare const DARK_THEME: CarmTheme;
/** Inject CSS custom properties into a container element. */
export declare function applyTheme(container: HTMLElement, theme?: CarmTheme): void;
/** Get a color from the palette by index (wraps around). */
export declare function getColor(index: number, theme?: CarmTheme): string;
/** D3-compatible color scale from theme palette. */
export declare function themeColorScale(theme?: CarmTheme): (i: number) => string;
//# sourceMappingURL=default.d.ts.map