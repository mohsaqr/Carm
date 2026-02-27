/**
 * Factor Analysis visualizations (6 plot types).
 * Scree with parallel analysis, loadings heatmap, CFA/SEM path diagram,
 * communality bar chart, factor correlation matrix, fit indices dashboard.
 *
 * All plots follow Carm viz conventions: dynamic D3 import, CarmTheme,
 * addSubtitle/addCaption annotations, tooltip interactions.
 */
import type { CarmTheme } from '../themes/default.js';
import type { FAResult, CFAResult, FADiagnostics } from '../../core/types.js';
export type FAPlotType = 'scree' | 'loadings' | 'path' | 'communality' | 'factor-correlation' | 'fit-indices';
/**
 * Obsessive control over every visual element of the SEM path diagram.
 * Every field has a sensible default — override only what you need.
 */
export interface FAPathStyle {
    readonly factorRx?: number;
    readonly factorRy?: number;
    readonly factorStroke?: number;
    readonly factorFontSize?: number;
    readonly factorFontWeight?: string;
    readonly factorGlow?: boolean;
    readonly factorGlowRadius?: number;
    readonly factorHighlight?: boolean;
    readonly itemWidth?: number;
    readonly itemHeight?: number;
    readonly itemRadius?: number;
    readonly itemStroke?: number;
    readonly itemFontSize?: number;
    readonly itemFontWeight?: string;
    readonly itemAccentWidth?: number;
    readonly itemGradient?: boolean;
    readonly errorRx?: number;
    readonly errorRy?: number;
    readonly errorStroke?: number;
    readonly errorFontSize?: number;
    readonly errorSelfLoop?: boolean;
    readonly errorSelfLoopSize?: number;
    readonly arrowMinWidth?: number;
    readonly arrowMaxWidth?: number;
    readonly arrowMinOpacity?: number;
    readonly arrowMaxOpacity?: number;
    readonly arrowMarkerSize?: number;
    readonly arrowCurvature?: number;
    readonly loadingFontSize?: number;
    readonly loadingFontWeight?: string;
    readonly loadingPillRadius?: number;
    readonly loadingLabelOffset?: number;
    readonly loadingPosition?: number;
    readonly halo?: boolean;
    readonly haloColor?: string;
    readonly haloWidth?: number;
    readonly covColor?: string;
    readonly covMinWidth?: number;
    readonly covMaxWidth?: number;
    readonly covMarkerSize?: number;
    readonly covLabelFontSize?: number;
    readonly covNestSpacing?: number;
    readonly errorArrowWidth?: number;
    readonly errorArrowMarkerSize?: number;
    readonly itemGap?: number;
    readonly factorGroupGap?: number;
    readonly arrowSpan?: number;
    readonly errorSpan?: number;
    readonly covArcReserve?: number;
    readonly topPadding?: number;
    readonly bottomPadding?: number;
    readonly rightPadding?: number;
    readonly factorColors?: readonly string[];
    readonly fitCardWidth?: number;
    readonly fitCardHeight?: number;
    readonly crossLoadingThreshold?: number;
    readonly crossLoadingDash?: string;
    readonly crossLoadingOpacity?: number;
    readonly showGroupBrackets?: boolean;
    readonly groupBracketOpacity?: number;
    readonly showShadows?: boolean;
}
export interface FAPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly type?: FAPlotType;
    readonly variableLabels?: readonly string[];
    readonly factorLabels?: readonly string[];
    readonly loadingThreshold?: number;
    readonly showParallelAnalysis?: boolean;
    readonly showErrorTerms?: boolean;
    readonly showFitBox?: boolean;
    readonly diagnostics?: FADiagnostics;
    readonly pathStyle?: FAPathStyle;
}
export declare function renderFAPlot(container: HTMLElement, data: FAResult | CFAResult, config?: FAPlotConfig): void;
//# sourceMappingURL=fa-plot.d.ts.map