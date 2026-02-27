/**
 * Mosaic plot — column widths proportional to column marginals, row heights
 * proportional to within-column conditional proportions.
 *
 * Cells coloured by discretised Pearson residual (RdBu ColorBrewer bins):
 *   > +4   dark blue     #2166ac
 *  +2..+4  medium blue   #74add1
 *   0..+2  light blue    #d1e5f0
 *  -2.. 0  light peach   #fddbc7
 *  -4..-2  salmon        #d6604d
 *   < -4   dark red      #b2182b
 *
 * Borders:
 *   |residual| ≥ 4  → thick solid (2.5 px)
 *   |residual| ≥ 2  → dashed (stroke-dasharray 4 2)
 *   otherwise       → thin solid (0.5 px, background colour)
 */
import type { CarmTheme } from '../themes/default.js';
import type { StatResult } from '../../core/types.js';
export interface MosaicPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
export interface MosaicPlotData {
    /** table[i][j] = count for row i, column j */
    readonly table: readonly (readonly number[])[];
    readonly rowLabels: readonly string[];
    readonly colLabels: readonly string[];
    readonly testResult?: StatResult;
}
export declare function renderMosaicPlot(container: HTMLElement, data: MosaicPlotData, config?: MosaicPlotConfig): void;
//# sourceMappingURL=mosaic-plot.d.ts.map