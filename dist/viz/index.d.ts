import * as d3 from 'd3';
import { j as PairwiseResult, G as GroupData, S as StatResult, k as RegressionResult, b as DescriptiveResult, f as FrequencyRow, R as RegressionCoef, L as LMMResult, i as PCAResult } from '../types-DC8rlZlK.js';
import { C as CorrelationMatrix } from '../correlation-CicTScAu.js';

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
declare const CARM_PALETTE: readonly ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#af7aa1", "#ff9da7", "#9c755f"];
/** @deprecated Use CARM_PALETTE. Kept for backwards compatibility. */
declare const OKABE_ITO: readonly ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#af7aa1", "#ff9da7", "#9c755f"];
type ThemeName = 'light' | 'dark';
interface CarmTheme {
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
declare const DEFAULT_THEME: CarmTheme;
declare const DARK_THEME: CarmTheme;
/** Inject CSS custom properties into a container element. */
declare function applyTheme(container: HTMLElement, theme?: CarmTheme): void;
/** Get a color from the palette by index (wraps around). */
declare function getColor(index: number, theme?: CarmTheme): string;
/** D3-compatible color scale from theme palette. */
declare function themeColorScale(theme?: CarmTheme): (i: number) => string;

/**
 * Hover tooltip component.
 * Shows a floating tooltip near the cursor with formatted content.
 */

declare function showTooltip(event: MouseEvent, content: string, theme?: CarmTheme): void;
declare function hideTooltip(): void;
declare function formatTooltipRow(label: string, value: string | number): string;

/**
 * Statistical annotation helpers for D3 plots.
 * Renders APA text, regression equations, and stat result boxes.
 */

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
declare function addSubtitle(svg: SVGSelection, title: string, subtitle: string, _width: number, theme?: CarmTheme): void;
/**
 * Add an italic caption at the bottom-left of the SVG.
 * Used for data source, method notes, sample size.
 */
declare function addCaption(svg: SVGSelection, text: string, _width: number, height: number, theme?: CarmTheme): void;
/**
 * Add a regression equation text annotation on the plot area.
 * Uses monospace font so numbers align cleanly.
 */
declare function addRegressionEquation(g: GSelection, intercept: number, slope: number, r2: number, x: number, y: number, theme?: CarmTheme): void;
/**
 * Add an n= label below a group's x-position.
 */
declare function addNLabel(g: GSelection, n: number, x: number, y: number, theme?: CarmTheme): void;
/**
 * Add a stat annotation box (pill) directly on the plot area.
 * E.g. AUC = 0.82, r = .91, p < .001
 */
declare function addStatBadge(g: GSelection, lines: string[], x: number, y: number, theme?: CarmTheme): void;

/**
 * Significance bracket layout for group comparison plots.
 * Renders stacked brackets above groups with p-value labels.
 * Non-overlapping bracket stacking algorithm.
 */

interface BracketConfig {
    readonly groupPositions: ReadonlyMap<string, number>;
    readonly yBase: number;
    readonly bracketHeight: number;
    readonly significantOnly: boolean;
}
/**
 * Render significance brackets onto an SVG group.
 * Uses greedy level-assignment to avoid bracket overlap.
 */
declare function renderBrackets(g: d3.Selection<SVGGElement, unknown, null, undefined>, comparisons: readonly PairwiseResult[], config: BracketConfig, theme?: CarmTheme): void;
/** Total bracket height needed (to set top margin). */
declare function totalBracketHeight(comparisons: readonly PairwiseResult[], significantOnly: boolean, bracketHeight: number): number;

/**
 * Smart axis rendering for Carm plots.
 * Auto tick formatting, minimal grid, clean style.
 */

type AnyScale = d3.ScaleLinear<number, number> | d3.ScaleBand<string>;
/** Render a styled x-axis with label. */
declare function renderXAxis(g: d3.Selection<SVGGElement, unknown, null, undefined>, _scale: AnyScale, height: number, label: string, width: number, theme?: CarmTheme): void;
/** Render a styled y-axis with label. */
declare function renderYAxis(g: d3.Selection<SVGGElement, unknown, null, undefined>, _scale: d3.ScaleLinear<number, number>, height: number, label: string, theme?: CarmTheme): void;
/** Render horizontal grid lines. */
declare function renderGridLines(g: d3.Selection<SVGGElement, unknown, null, undefined>, scale: d3.ScaleLinear<number, number>, width: number, theme?: CarmTheme): void;

/**
 * Violin + Box + Jitter combo plot (ggbetweenstats style).
 * Shows distribution shape (violin), summary (box), individual points (jitter),
 * and significance brackets between groups.
 */

interface ViolinBoxConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showJitter?: boolean;
    readonly showBrackets?: boolean;
    readonly significantBracketsOnly?: boolean;
    readonly jitterWidth?: number;
    readonly violinBandwidth?: number;
}
interface ViolinBoxData {
    readonly groups: readonly GroupData[];
    readonly testResult?: StatResult;
    readonly pairwise?: readonly PairwiseResult[];
}
/**
 * Render a violin + box + jitter plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - group data + optional stat results
 * @param config - visual configuration
 */
declare function renderViolinBox(container: HTMLElement, data: ViolinBoxData, config?: ViolinBoxConfig): void;

/**
 * Scatter plot with regression line, CI band, and marginal distributions.
 * ggscatterstats style: stat result in subtitle, regression equation on plot.
 */

interface ScatterStatsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showMarginals?: boolean;
    readonly showCI?: boolean;
    readonly pointSize?: number;
}
interface ScatterStatsData {
    readonly x: readonly number[];
    readonly y: readonly number[];
    readonly labels?: readonly string[];
    readonly correlationResult?: StatResult;
    readonly regressionResult?: RegressionResult;
}
declare function renderScatterStats(container: HTMLElement, data: ScatterStatsData, config?: ScatterStatsConfig): void;

/**
 * Histogram with density overlay and optional normality curve.
 */

interface HistogramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly bins?: number;
    readonly showDensity?: boolean;
    readonly showNormalCurve?: boolean;
    readonly color?: string;
}
interface HistogramData {
    readonly values: readonly number[];
    readonly descriptives?: DescriptiveResult;
}
declare function renderHistogram(container: HTMLElement, data: HistogramData, config?: HistogramConfig): void;

/**
 * Annotated bar chart with counts, percentages, and significance brackets.
 */

interface BarStatsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showPercentages?: boolean;
    readonly showCounts?: boolean;
    readonly orientation?: 'vertical' | 'horizontal';
}
interface BarStatsData {
    readonly rows: readonly FrequencyRow[];
    readonly testResult?: StatResult;
}
declare function renderBarStats(container: HTMLElement, data: BarStatsData, config?: BarStatsConfig): void;

/**
 * Correlation heatmap (correlogram) with significance stars and clustering.
 */

interface CorrelogramConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showValues?: boolean;
    readonly showSignificance?: boolean;
}
declare function renderCorrelogram(container: HTMLElement, data: CorrelationMatrix, config?: CorrelogramConfig): void;

/**
 * Coefficient dot-and-whisker plot (ggcoefstats style).
 * Shows regression/LMM coefficients with CI bars and p-value annotations.
 */

interface CoefPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showZeroLine?: boolean;
    readonly excludeIntercept?: boolean;
}
declare function renderCoefPlot(container: HTMLElement, coefficients: readonly RegressionCoef[], config?: CoefPlotConfig): void;

/**
 * QQ plot (quantile-quantile plot) with confidence band.
 * Tests normality assumption visually.
 */

interface QQPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly showCI?: boolean;
}
declare function renderQQPlot(container: HTMLElement, values: readonly number[], config?: QQPlotConfig): void;

/**
 * 4-panel regression diagnostics plot:
 * 1. Residuals vs Fitted
 * 2. QQ of residuals
 * 3. Scale-Location (sqrt|residuals| vs fitted)
 * 4. Residuals vs Leverage
 */

interface ResidualPanelConfig {
    readonly title?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
declare function renderResidualPanel(container: HTMLElement, result: RegressionResult, leverage: readonly number[], config?: ResidualPanelConfig): void;

/**
 * Raincloud plot: half-violin + jitter + box.
 * Premium visualization for showing full distribution shape.
 */

interface RaincloudConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface RaincloudData {
    readonly groups: readonly GroupData[];
    readonly testResult?: StatResult;
}
declare function renderRaincloud(container: HTMLElement, data: RaincloudData, config?: RaincloudConfig): void;

/**
 * LMM visualization: caterpillar plot of BLUPs + variance component chart.
 */

interface MixedPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
declare function renderMixedPlot(container: HTMLElement, result: LMMResult, blups: ReadonlyArray<{
    group: string | number;
    blup: number;
}>, config?: MixedPlotConfig): void;

/**
 * PCA visualizations: biplot, scree plot, loadings heatmap.
 */

interface PCAPlotConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly variableLabels?: readonly string[];
    readonly observationLabels?: readonly string[];
    readonly type?: 'biplot' | 'scree' | 'loadings';
}
declare function renderPCAPlot(container: HTMLElement, pca: PCAResult, config?: PCAPlotConfig): void;

/**
 * Interactive distribution explorer.
 * Shows PDF/CDF for normal, t, chi-square, F, binomial, Poisson distributions.
 * Uses jStat for accurate distribution functions.
 */

type DistributionName = 'normal' | 't' | 'chi-square' | 'F' | 'uniform' | 'exponential';
interface DistributionParams {
    readonly distribution: DistributionName;
    readonly params: Readonly<Record<string, number>>;
    readonly showPDF?: boolean;
    readonly showCDF?: boolean;
    readonly highlightX?: number;
}
interface DistributionConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
declare function renderDistribution(container: HTMLElement, params: DistributionParams, config?: DistributionConfig): void;

/**
 * Density plot: KDE curves for one or more series with optional rug ticks.
 * Uses Silverman's rule-of-thumb bandwidth: 1.06 * σ * n^(-1/5)
 */

interface DensityConfig {
    readonly bandwidth?: number;
    readonly showRug?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface DensitySeries {
    readonly label: string;
    readonly values: readonly number[];
}
interface DensityData {
    readonly series: readonly DensitySeries[];
    readonly testResult?: StatResult;
}
declare function renderDensity(container: HTMLElement, data: DensityData, config?: DensityConfig): void;

/**
 * Standalone box plot: Q1/Q3/median/whiskers (1.5 IQR)/outliers per group.
 * No violin. Uses quantile() from core/math.
 */

interface BoxplotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface BoxplotGroup {
    readonly label: string;
    readonly values: readonly number[];
}
interface BoxplotData {
    readonly groups: readonly BoxplotGroup[];
    readonly testResult?: StatResult;
}
declare function renderBoxplot(container: HTMLElement, data: BoxplotData, config?: BoxplotConfig): void;

/**
 * Lollipop chart: stem + circle for each category.
 * Horizontal layout (categories on y-axis, values on x-axis).
 * Optional descending sort.
 */

interface LollipopConfig {
    readonly sorted?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface LollipopData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
declare function renderLollipop(container: HTMLElement, data: LollipopData, config?: LollipopConfig): void;

/**
 * Cleveland dot plot: horizontal layout with labels on y-axis.
 * Supports single dots or paired dots connected by a line (group comparison).
 */

interface DotPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface DotPlotData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly group2?: readonly number[];
    readonly group1Label?: string;
    readonly group2Label?: string;
    readonly testResult?: StatResult;
}
declare function renderDotPlot(container: HTMLElement, data: DotPlotData, config?: DotPlotConfig): void;

/**
 * Grouped or stacked bar chart: multiple series per category.
 * type='grouped' renders bars side-by-side; type='stacked' stacks them.
 */

interface GroupedBarConfig {
    readonly type?: 'grouped' | 'stacked';
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface GroupedBarSeries {
    readonly label: string;
    readonly values: readonly number[];
}
interface GroupedBarData {
    readonly categories: readonly string[];
    readonly series: readonly GroupedBarSeries[];
    readonly testResult?: StatResult;
}
declare function renderGroupedBar(container: HTMLElement, data: GroupedBarData, config?: GroupedBarConfig): void;

/**
 * Multi-line chart with optional area fill.
 * Each series has its own line + color. Points at each data location.
 */

interface LineChartConfig {
    readonly showArea?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface LineChartSeries {
    readonly label: string;
    readonly x: readonly number[];
    readonly y: readonly number[];
}
interface LineChartData {
    readonly series: readonly LineChartSeries[];
    readonly testResult?: StatResult;
}
declare function renderLineChart(container: HTMLElement, data: LineChartData, config?: LineChartConfig): void;

/**
 * Bubble chart: scatter plot where bubble size is proportional to r value.
 * Color by group if group field is present. Tooltip on hover.
 */

interface BubbleChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface BubblePoint {
    readonly x: number;
    readonly y: number;
    readonly r: number;
    readonly label?: string;
    readonly group?: string;
}
interface BubbleChartData {
    readonly points: readonly BubblePoint[];
    readonly testResult?: StatResult;
}
declare function renderBubbleChart(container: HTMLElement, data: BubbleChartData, config?: BubbleChartConfig): void;

/**
 * Pareto chart: descending bars + cumulative % line on dual axis.
 * Includes 80% threshold line (the "vital few" cutoff).
 */

interface ParetoConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface ParetoData {
    readonly labels: readonly string[];
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
declare function renderPareto(container: HTMLElement, data: ParetoData, config?: ParetoConfig): void;

/**
 * Funnel chart: trapezoidal segments stacked vertically.
 * Width proportional to value. Labels on left, value + % drop on right.
 */

interface FunnelConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface FunnelStage {
    readonly label: string;
    readonly value: number;
}
interface FunnelData {
    readonly stages: readonly FunnelStage[];
    readonly testResult?: StatResult;
}
declare function renderFunnel(container: HTMLElement, data: FunnelData, config?: FunnelConfig): void;

/**
 * Pie / donut chart with polyline labels and optional test result subtitle.
 * Donut mode uses innerRadius = outerRadius * 0.55.
 */

interface PieChartConfig {
    readonly donut?: boolean;
    readonly showPercentages?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface PieSlice {
    readonly label: string;
    readonly value: number;
}
interface PieChartData {
    readonly slices: readonly PieSlice[];
    readonly testResult?: StatResult;
}
declare function renderPieChart(container: HTMLElement, data: PieChartData, config?: PieChartConfig): void;

/**
 * Area chart with optional stacking.
 * Non-stacked: overlapping semi-transparent areas + lines per series.
 * Stacked: cumulative areas per series.
 */

interface AreaChartConfig {
    readonly stacked?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface AreaSeries {
    readonly label: string;
    readonly x: readonly number[];
    readonly y: readonly number[];
}
interface AreaChartData {
    readonly series: readonly AreaSeries[];
    readonly testResult?: StatResult;
}
declare function renderAreaChart(container: HTMLElement, data: AreaChartData, config?: AreaChartConfig): void;

/**
 * Forest plot: per-study CI lines with square points + pooled diamond.
 * Vertical reference line at 0. Numeric estimates shown on the right.
 */

interface ForestPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface ForestStudy {
    readonly label: string;
    readonly estimate: number;
    readonly ciLow: number;
    readonly ciHigh: number;
    readonly weight?: number;
}
interface ForestPooled {
    readonly estimate: number;
    readonly ciLow: number;
    readonly ciHigh: number;
}
interface ForestPlotData {
    readonly studies: readonly ForestStudy[];
    readonly pooled?: ForestPooled;
    readonly testResult?: StatResult;
}
declare function renderForestPlot(container: HTMLElement, data: ForestPlotData, config?: ForestPlotConfig): void;

/**
 * ROC curve: TPR vs FPR with AUC annotation.
 * Includes diagonal reference line, filled area under curve.
 */

interface ROCCurveConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface ROCCurveData {
    readonly fpr: readonly number[];
    readonly tpr: readonly number[];
    readonly auc: number;
    readonly testResult?: StatResult;
}
declare function renderROCCurve(container: HTMLElement, data: ROCCurveData, config?: ROCCurveConfig): void;

/**
 * Strip plot: pure jittered dot plot per group.
 * Shows individual data points with a horizontal mean line per group.
 * No violin, no box — raw distribution only.
 */

interface StripPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface StripGroup {
    readonly label: string;
    readonly values: readonly number[];
}
interface StripPlotData {
    readonly groups: readonly StripGroup[];
    readonly testResult?: StatResult;
}
declare function renderStripPlot(container: HTMLElement, data: StripPlotData, config?: StripPlotConfig): void;

/**
 * Beeswarm plot — individual points per group arranged to avoid overlap.
 * Points are sorted by value and offset horizontally so they spread
 * left/right of the group centre axis without colliding.
 */

interface SwarmPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly pointRadius?: number;
    readonly bandwidth?: number;
}
interface SwarmPlotData {
    readonly groups: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly testResult?: StatResult;
}
/**
 * Render a beeswarm plot with per-group colour and collision-avoidance layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - groups of numeric values + optional stat result
 * @param config - visual configuration
 */
declare function renderSwarmPlot(container: HTMLElement, data: SwarmPlotData, config?: SwarmPlotConfig): void;

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

interface MosaicPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface MosaicPlotData {
    /** table[i][j] = count for row i, column j */
    readonly table: readonly (readonly number[])[];
    readonly rowLabels: readonly string[];
    readonly colLabels: readonly string[];
    readonly testResult?: StatResult;
}
declare function renderMosaicPlot(container: HTMLElement, data: MosaicPlotData, config?: MosaicPlotConfig): void;

/**
 * Scatter matrix (pairs plot) — n×n grid of mini-plots.
 * Off-diagonal cells: scatter plot of variable i vs variable j.
 * Diagonal cells: histogram of variable i.
 * Each cell is a self-contained SVG <g> fitting in the available space.
 */

interface PairPlotConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface PairPlotData {
    /** data[varIndex][obsIndex] — each inner array is one variable's observations */
    readonly data: readonly (readonly number[])[];
    readonly labels: readonly string[];
    readonly testResult?: StatResult;
}
/**
 * Render a scatter matrix (pairs plot).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - variable arrays + labels + optional stat result
 * @param config - visual configuration
 */
declare function renderPairPlot(container: HTMLElement, data: PairPlotData, config?: PairPlotConfig): void;

/**
 * Radar / spider chart — one polygon per series, one spoke per axis.
 * Values are normalized per axis to [0,1] using min/max across all series.
 * Spokes are uniformly spaced angularly. Filled polygons with opacity.
 * Axis labels at spoke tips. Legend for series.
 */

interface RadarChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly levels?: number;
}
interface RadarChartData {
    readonly series: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly axes: readonly string[];
    readonly testResult?: StatResult;
}
/**
 * Render a radar / spider chart.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - series (each with values per axis) + axis names + optional stat result
 * @param config - visual configuration
 */
declare function renderRadarChart(container: HTMLElement, data: RadarChartData, config?: RadarChartConfig): void;

/**
 * Parallel coordinates plot — one vertical axis per variable, evenly spaced
 * horizontally. Each data row is a polyline through its per-axis values.
 * Axes are individually min-max scaled. Rows are optionally coloured by group.
 * Lines are semi-transparent to show density.
 */

interface ParallelCoordsConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface ParallelCoordsData {
    /** rows[obsIndex][varIndex] */
    readonly rows: readonly (readonly number[])[];
    readonly axes: readonly string[];
    /** optional group index per row (0-based) for colour coding */
    readonly groups?: readonly number[];
    readonly testResult?: StatResult;
}
/**
 * Render a parallel coordinates plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - observations (rows × vars), axis names, optional groups + stat result
 * @param config - visual configuration
 */
declare function renderParallelCoords(container: HTMLElement, data: ParallelCoordsData, config?: ParallelCoordsConfig): void;

/**
 * Treemap — rectangle area proportional to value, laid out using D3's
 * squarified treemap algorithm. Cells are labelled (truncated when small)
 * and coloured by group when provided, otherwise by index.
 * White borders separate cells. Tooltip shows label, value, and %.
 */

interface TreemapConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface TreemapChild {
    readonly label: string;
    readonly value: number;
    readonly group?: string;
}
interface TreemapData {
    readonly children: readonly TreemapChild[];
    readonly testResult?: StatResult;
}
/**
 * Render a treemap using D3's squarified layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - leaf nodes with labels, values, optional groups + stat result
 * @param config - visual configuration
 */
declare function renderTreemap(container: HTMLElement, data: TreemapData, config?: TreemapConfig): void;

/**
 * Waffle chart — 10×10 grid of small squares where each square represents 1%
 * of the total. Squares are coloured by slice proportionally (floor-rounded,
 * remainder assigned to largest slice). Legend below shows colour, label, %.
 */

interface WaffleChartConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface WaffleChartData {
    readonly slices: readonly {
        readonly label: string;
        readonly value: number;
    }[];
    readonly testResult?: StatResult;
}
/**
 * Render a waffle chart (10×10 proportional grid).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - named slices with values + optional stat result
 * @param config - visual configuration
 */
declare function renderWaffleChart(container: HTMLElement, data: WaffleChartData, config?: WaffleChartConfig): void;

/**
 * Sparkline — minimal inline trend line for dashboard embedding.
 * Tiny margins, no axes, no labels, no grid. Just the line, an optional
 * shaded area, and a highlighted end-point dot.
 */

interface SparklineConfig {
    readonly showArea?: boolean;
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
}
interface SparklineData {
    readonly values: readonly number[];
    readonly testResult?: StatResult;
}
/**
 * Render a sparkline (minimal inline trend line).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - numeric time series values + optional stat result
 * @param config - visual configuration (very minimal by design)
 */
declare function renderSparkline(container: HTMLElement, data: SparklineData, config?: SparklineConfig): void;

/**
 * Sunburst chart — radial hierarchy visualization where each ring represents
 * one level of depth. Arc angles are proportional to node values. Center hole
 * shows total. Depth-1 nodes are listed in a right-side legend.
 * Inspired by ggstatsplot philosophy: every plot tells the full statistical story.
 */

interface SunburstNode {
    readonly name: string;
    readonly value?: number;
    readonly children?: readonly SunburstNode[];
}
interface SunburstData {
    readonly root: SunburstNode;
    readonly testResult?: StatResult;
}
interface SunburstConfig {
    readonly title?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly innerRadius?: number;
    readonly maxDepth?: number;
}
declare function renderSunburst(container: HTMLElement, data: SunburstData, config?: SunburstConfig): void;

/**
 * Marimekko (Mekko) chart — a 100% stacked bar chart where column WIDTHS are
 * proportional to each column's total value, so both axes encode data.
 *
 * - Column width  = that category's total as a fraction of the grand total.
 * - Within each column, stacked segments sum to 100%, coloured by series.
 * - Cells receive percentage (or raw value) labels when large enough.
 * - A "spine" bar below the x-axis shows relative column widths at a glance.
 *
 * Usage:
 *   renderMarimekko(container, data, config)
 */

interface MarimekkoData {
    /** Labels for each column (X variable). */
    readonly categories: readonly string[];
    /** One series per stacked segment. values[j] corresponds to categories[j]. */
    readonly series: readonly {
        readonly label: string;
        readonly values: readonly number[];
    }[];
    readonly testResult?: StatResult;
}
interface MarimekkoConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    /** Show percentage labels inside cells when the cell is large enough. Default: true. */
    readonly showPercentLabels?: boolean;
    /** Show raw values instead of percentages in cell labels. Default: false. */
    readonly showValueLabels?: boolean;
}
declare function renderMarimekko(container: HTMLElement, data: MarimekkoData, config?: MarimekkoConfig): void;

/**
 * Chord diagram: shows flows between groups as ribbons inside a circular layout.
 * Each group occupies an arc on the outer ring; ribbons inside represent flows.
 * Based on the D3 chord layout (d3.chord / d3.ribbon / d3.arc).
 */

interface ChordDiagramData {
    readonly matrix: readonly (readonly number[])[];
    readonly labels: readonly string[];
    readonly testResult?: StatResult;
}
interface ChordDiagramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly innerRadius?: number;
    readonly padAngle?: number;
}
declare function renderChordDiagram(container: HTMLElement, data: ChordDiagramData, config?: ChordDiagramConfig): void;

/**
 * Arc diagram: nodes on a horizontal baseline connected by arcs above them.
 * Short connections produce shallow arcs; distant connections produce tall arcs.
 * Useful for visualising network topology when node order carries meaning.
 */

interface ArcDiagramNode {
    readonly id: string;
    readonly label?: string;
    readonly group?: string;
}
interface ArcDiagramEdge {
    readonly source: string;
    readonly target: string;
    readonly value?: number;
}
interface ArcDiagramData {
    readonly nodes: readonly ArcDiagramNode[];
    readonly edges: readonly ArcDiagramEdge[];
    readonly testResult?: StatResult;
}
interface ArcDiagramConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly nodeRadius?: number;
    readonly sortByGroup?: boolean;
}
declare function renderArcDiagram(container: HTMLElement, data: ArcDiagramData, config?: ArcDiagramConfig): void;

/**
 * Alluvial (Sankey-like) plot — nodes stacked by stage, flows as cubic bezier
 * ribbons. Every stage is a vertical column; node height is proportional to
 * total flow value. Ribbons connect adjacent stages and are coloured by source
 * node colour at reduced opacity so layering is visible.
 *
 * Layout:
 *   stages are evenly spaced across the plot width.
 *   within each stage nodes are sorted largest→smallest, stacked top→bottom
 *   with a configurable gap between them.
 *   each node tracks separate y-offsets for outgoing and incoming ribbons so
 *   ribbons stack neatly without crossing inside a node.
 *
 * Colour: all unique node labels are collected up-front and assigned a fixed
 * colour from the theme palette so the same category always has the same
 * colour regardless of which stage it appears in.
 */

interface AlluvialNode {
    readonly id: string;
    readonly stage: number;
    readonly label: string;
}
interface AlluvialFlow {
    readonly source: string;
    readonly target: string;
    readonly value: number;
}
interface AlluvialData {
    readonly nodes: readonly AlluvialNode[];
    readonly flows: readonly AlluvialFlow[];
    readonly stageLabels?: readonly string[];
    readonly testResult?: StatResult;
}
interface AlluvialConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly nodePadding?: number;
    readonly nodeWidth?: number;
}
/**
 * Render an alluvial / Sankey-style plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes, flows, optional stage labels and stat result
 * @param config    - visual configuration
 */
declare function renderAlluvialPlot(container: HTMLElement, data: AlluvialData, config?: AlluvialConfig): void;

/**
 * Hierarchical edge bundling (Holten 2006).
 *
 * Leaf nodes are arranged on a circle. Non-leaf nodes form the hierarchy used
 * to bundle edges — they are not drawn. Edges between leaf nodes are drawn as
 * cubic B-spline paths routed through the common ancestor path in the
 * hierarchy, which causes edges that share ancestry to visually cluster
 * ("bundle") together.
 *
 * Layout:
 *   d3.stratify builds the hierarchy from the flat node list.
 *   d3.cluster assigns angular positions (x) and radial depth (y).
 *   Only leaf nodes (no children) are shown as circles with radial labels.
 *   Each edge walks the path [source → LCA → target] through the tree and
 *   passes those anchor points to d3.lineRadial with d3.curveBundle.beta(β).
 *
 * Colour:
 *   Nodes are coloured by their `group` field if provided, otherwise by their
 *   immediate parent id. The same colour mapping is reused for edges (source
 *   colour at low opacity so overlapping edges remain readable).
 *
 * Interactivity:
 *   Hovering a node highlights all its edges (full opacity, thicker stroke)
 *   and dims all other edges. Clicking a node locks the highlight; clicking
 *   again unlocks it.
 */

interface BundleNode {
    readonly id: string;
    readonly parent: string;
    readonly label?: string;
    readonly group?: string;
}
interface BundleEdge {
    readonly source: string;
    readonly target: string;
}
interface EdgeBundlingData {
    readonly nodes: readonly BundleNode[];
    readonly edges: readonly BundleEdge[];
    readonly testResult?: StatResult;
}
interface EdgeBundlingConfig {
    readonly title?: string;
    readonly xLabel?: string;
    readonly yLabel?: string;
    readonly caption?: string;
    readonly width?: number;
    readonly height?: number;
    readonly theme?: CarmTheme;
    readonly bundlingStrength?: number;
    readonly nodeRadius?: number;
}
/**
 * Render a hierarchical edge bundling diagram.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes (with parent hierarchy), edges between leaves, optional stat result
 * @param config    - visual configuration
 */
declare function renderEdgeBundling(container: HTMLElement, data: EdgeBundlingData, config?: EdgeBundlingConfig): void;

/**
 * SVG + PNG export for Carm plots.
 * Exports at publication quality (300 DPI for PNG).
 */
/** Download SVG source of a plot container. */
declare function exportSVG(container: HTMLElement, filename?: string): void;
/** Export plot container as PNG at specified DPI (default 300). */
declare function exportPNG(container: HTMLElement, filename?: string, dpi?: number): Promise<void>;

export { type AlluvialConfig, type AlluvialData, type AlluvialFlow, type AlluvialNode, type ArcDiagramConfig, type ArcDiagramData, type ArcDiagramEdge, type ArcDiagramNode, type AreaChartConfig, type AreaChartData, type AreaSeries, type BarStatsConfig, type BarStatsData, type BoxplotConfig, type BoxplotData, type BoxplotGroup, type BracketConfig, type BubbleChartConfig, type BubbleChartData, type BubblePoint, type BundleEdge, type BundleNode, CARM_PALETTE, type CarmTheme, type ChordDiagramConfig, type ChordDiagramData, type CoefPlotConfig, type CorrelogramConfig, DARK_THEME, DEFAULT_THEME, type DensityConfig, type DensityData, type DensitySeries, type DistributionConfig, type DistributionName, type DistributionParams, type DotPlotConfig, type DotPlotData, type EdgeBundlingConfig, type EdgeBundlingData, type ForestPlotConfig, type ForestPlotData, type ForestPooled, type ForestStudy, type FunnelConfig, type FunnelData, type FunnelStage, type GroupedBarConfig, type GroupedBarData, type GroupedBarSeries, type HistogramConfig, type HistogramData, type LineChartConfig, type LineChartData, type LineChartSeries, type LollipopConfig, type LollipopData, type MarimekkoConfig, type MarimekkoData, type MixedPlotConfig, type MosaicPlotConfig, type MosaicPlotData, OKABE_ITO, type PCAPlotConfig, type PairPlotConfig, type PairPlotData, type ParallelCoordsConfig, type ParallelCoordsData, type ParetoConfig, type ParetoData, type PieChartConfig, type PieChartData, type PieSlice, type QQPlotConfig, type ROCCurveConfig, type ROCCurveData, type RadarChartConfig, type RadarChartData, type RaincloudConfig, type RaincloudData, type ResidualPanelConfig, type ScatterStatsConfig, type ScatterStatsData, type SparklineConfig, type SparklineData, type StripGroup, type StripPlotConfig, type StripPlotData, type SunburstConfig, type SunburstData, type SunburstNode, type SwarmPlotConfig, type SwarmPlotData, type ThemeName, type TreemapChild, type TreemapConfig, type TreemapData, type ViolinBoxConfig, type ViolinBoxData, type WaffleChartConfig, type WaffleChartData, addCaption, addNLabel, addRegressionEquation, addStatBadge, addSubtitle, applyTheme, exportPNG, exportSVG, formatTooltipRow, getColor, hideTooltip, renderAlluvialPlot, renderArcDiagram, renderAreaChart, renderBarStats, renderBoxplot, renderBrackets, renderBubbleChart, renderChordDiagram, renderCoefPlot, renderCorrelogram, renderDensity, renderDistribution, renderDotPlot, renderEdgeBundling, renderForestPlot, renderFunnel, renderGridLines, renderGroupedBar, renderHistogram, renderLineChart, renderLollipop, renderMarimekko, renderMixedPlot, renderMosaicPlot, renderPCAPlot, renderPairPlot, renderParallelCoords, renderPareto, renderPieChart, renderQQPlot, renderROCCurve, renderRadarChart, renderRaincloud, renderResidualPanel, renderScatterStats, renderSparkline, renderStripPlot, renderSunburst, renderSwarmPlot, renderTreemap, renderViolinBox, renderWaffleChart, renderXAxis, renderYAxis, showTooltip, themeColorScale, totalBracketHeight };
