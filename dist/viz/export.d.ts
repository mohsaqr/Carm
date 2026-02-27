/**
 * SVG + PNG export for Carm plots.
 * Exports at publication quality (300 DPI for PNG).
 */
/** Download SVG source of a plot container. */
export declare function exportSVG(container: HTMLElement, filename?: string): void;
/** Export plot container as PNG at specified DPI (default 300). */
export declare function exportPNG(container: HTMLElement, filename?: string, dpi?: number): Promise<void>;
//# sourceMappingURL=export.d.ts.map