/**
 * Statistical annotation helpers for D3 plots.
 * Renders APA text, regression equations, and stat result boxes.
 */
import type * as d3 from 'd3';
import type { CarmTheme } from '../themes/default.js';
type SVGSelection = d3.Selection<SVGSVGElement, unknown, null, undefined>;
type GSelection = d3.Selection<SVGGElement, unknown, null, undefined>;
/**
 * Add a title + statistical subtitle at the top of an SVG.
 *
 * Layout (y positions, px from SVG top):
 *   20 — decorative accent bar (2×20 px, theme color)
 *   22 — title baseline  (16 px bold, near-black)
 *   43 — subtitle baseline (11 px italic sans-serif, slate)
 *
 * Theme marginTop should be ≥ 58 so the plot area starts below the subtitle.
 */
export declare function addSubtitle(svg: SVGSelection, title: string, subtitle: string, _width: number, theme?: CarmTheme): void;
/**
 * Add an italic caption at the bottom-left of the SVG.
 * Used for data source, method notes, sample size.
 */
export declare function addCaption(svg: SVGSelection, text: string, _width: number, height: number, theme?: CarmTheme): void;
/**
 * Add a regression equation text annotation on the plot area.
 * Uses monospace font so numbers align cleanly.
 */
export declare function addRegressionEquation(g: GSelection, intercept: number, slope: number, r2: number, x: number, y: number, theme?: CarmTheme): void;
/**
 * Add an n= label below a group's x-position.
 */
export declare function addNLabel(g: GSelection, n: number, x: number, y: number, theme?: CarmTheme): void;
/**
 * Add a stat annotation box (pill) directly on the plot area.
 * E.g. AUC = 0.82, r = .91, p < .001
 */
export declare function addStatBadge(g: GSelection, lines: string[], x: number, y: number, theme?: CarmTheme): void;
export {};
//# sourceMappingURL=annotations.d.ts.map