/**
 * Factor Analysis visualizations (6 plot types).
 * Scree with parallel analysis, loadings heatmap, CFA/SEM path diagram,
 * communality bar chart, factor correlation matrix, fit indices dashboard.
 *
 * All plots follow Carm viz conventions: dynamic D3 import, CarmTheme,
 * addSubtitle/addCaption annotations, tooltip interactions.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import { formatStat, formatP } from '../../core/apa.js'
import type { FAResult, CFAResult, FADiagnostics } from '../../core/types.js'

// ─── Config ──────────────────────────────────────────────────────────────

export type FAPlotType = 'scree' | 'loadings' | 'path' | 'communality' | 'factor-correlation' | 'fit-indices'

/**
 * Obsessive control over every visual element of the SEM path diagram.
 * Every field has a sensible default — override only what you need.
 */
export interface FAPathStyle {
  // ── Factor ellipses ──
  readonly factorRx?: number              // default: 56
  readonly factorRy?: number              // default: 28
  readonly factorStroke?: number          // default: 2.5
  readonly factorFontSize?: number        // default: 14
  readonly factorFontWeight?: string      // default: '800'
  readonly factorGlow?: boolean           // default: true — outer glow on factor ellipses
  readonly factorGlowRadius?: number      // default: 6
  readonly factorHighlight?: boolean      // default: true — inner specular shine arc

  // ── Item rectangles ──
  readonly itemWidth?: number             // default: 96
  readonly itemHeight?: number            // default: 30
  readonly itemRadius?: number            // default: 5 — corner radius
  readonly itemStroke?: number            // default: 1.2
  readonly itemFontSize?: number          // default: 11
  readonly itemFontWeight?: string        // default: '600'
  readonly itemAccentWidth?: number       // default: 4 — left color bar width
  readonly itemGradient?: boolean         // default: true — subtle left→right gradient fill

  // ── Error ellipses ──
  readonly errorRx?: number               // default: 14
  readonly errorRy?: number               // default: 12
  readonly errorStroke?: number           // default: 1.2
  readonly errorFontSize?: number         // default: 9
  readonly errorSelfLoop?: boolean        // default: true — self-loop variance arc
  readonly errorSelfLoopSize?: number     // default: 12 — height of self-loop arc above error

  // ── Loading arrows (factor → item) ──
  readonly arrowMinWidth?: number         // default: 0.8
  readonly arrowMaxWidth?: number         // default: 2.5
  readonly arrowMinOpacity?: number       // default: 0.4
  readonly arrowMaxOpacity?: number       // default: 0.95
  readonly arrowMarkerSize?: number       // default: 6 — arrowhead size
  readonly arrowCurvature?: number        // default: 0.4 — Bézier control point offset (0=straight, 1=extreme)

  // ── Loading labels ──
  readonly loadingFontSize?: number       // default: 9.5
  readonly loadingFontWeight?: string     // default: '700'
  readonly loadingPillRadius?: number     // default: 7
  readonly loadingLabelOffset?: number    // default: -6 — vertical offset from arrow midpoint (negative = above)
  readonly loadingPosition?: number       // default: 0.6 — t param along Bézier (0=factor, 1=item)

  // ── Text halo (white outline behind labels for legibility) ──
  readonly halo?: boolean                 // default: true — enable white halo on all path labels
  readonly haloColor?: string             // default: '#ffffff'
  readonly haloWidth?: number             // default: 3.5 — stroke width of halo outline

  // ── Covariance arcs ──
  readonly covColor?: string              // default: '#8892a0'
  readonly covMinWidth?: number           // default: 1.0
  readonly covMaxWidth?: number           // default: 3.0
  readonly covMarkerSize?: number         // default: 6
  readonly covLabelFontSize?: number      // default: 10
  readonly covNestSpacing?: number        // default: 20 — extra X offset per nesting level

  // ── Error arrows (error → item) ──
  readonly errorArrowWidth?: number       // default: 0.8
  readonly errorArrowMarkerSize?: number  // default: 5

  // ── Layout spacing ──
  readonly itemGap?: number               // default: 6 — vertical gap between items
  readonly factorGroupGap?: number        // default: 32 — vertical gap between factor groups
  readonly arrowSpan?: number             // default: 120 — horizontal distance factor→item
  readonly errorSpan?: number             // default: 36 — horizontal distance item→error
  readonly covArcReserve?: number         // default: 70 — left space for covariance arcs
  readonly topPadding?: number            // default: 74
  readonly bottomPadding?: number         // default: 44 (140 if showFitBox)
  readonly rightPadding?: number          // default: 50

  // ── Colors (override palette per factor) ──
  readonly factorColors?: readonly string[] // default: uses theme palette

  // ── Fit card ──
  readonly fitCardWidth?: number          // default: 210
  readonly fitCardHeight?: number         // default: 118

  // ── Cross-loadings (EFA) ──
  readonly crossLoadingThreshold?: number  // default: 0.2
  readonly crossLoadingDash?: string       // default: '5,4'
  readonly crossLoadingOpacity?: number    // default: 0.3

  // ── Group bracket lines ──
  readonly showGroupBrackets?: boolean     // default: true
  readonly groupBracketOpacity?: number    // default: 0.12

  // ── Shadow ──
  readonly showShadows?: boolean           // default: true
}

export interface FAPlotConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly type?: FAPlotType               // default: 'loadings'
  readonly variableLabels?: readonly string[]
  readonly factorLabels?: readonly string[]
  readonly loadingThreshold?: number       // default: 0.3 (for loadings heatmap)
  readonly showParallelAnalysis?: boolean  // default: true (for scree)
  readonly showErrorTerms?: boolean        // default: true (for path)
  readonly showFitBox?: boolean            // default: true (for path)
  readonly diagnostics?: FADiagnostics     // for scree plot parallel analysis overlay
  readonly pathStyle?: FAPathStyle         // obsessive control over path diagram
}

// ─── Entry Point ─────────────────────────────────────────────────────────

export function renderFAPlot(
  container: HTMLElement,
  data: FAResult | CFAResult,
  config: FAPlotConfig = {}
): void {
  import('d3').then(d3 => {
    const type = config.type ?? 'loadings'
    if (type === 'scree') renderScree(d3, container, data, config)
    else if (type === 'loadings') renderLoadings(d3, container, data, config)
    else if (type === 'path') renderPath(d3, container, data, config)
    else if (type === 'communality') renderCommunality(d3, container, data, config)
    else if (type === 'factor-correlation') renderFactorCorrelation(d3, container, data, config)
    else if (type === 'fit-indices') renderFitIndices(d3, container, data, config)
  })
}

// ─── Helper: check if CFAResult ──────────────────────────────────────────

function isCFA(data: FAResult | CFAResult): data is CFAResult {
  return 'parameterEstimates' in data
}

// ─── Plot 1: Scree with Parallel Analysis ────────────────────────────────

function renderScree(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 520, H = config.height ?? 380
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)

  const diag = config.diagnostics
  const subtitleParts: string[] = []
  if (diag) {
    subtitleParts.push(`Parallel Analysis: ${diag.parallelSuggested} factors`)
    subtitleParts.push(`MAP: ${diag.mapSuggested} factors`)
  }
  addSubtitle(svg, config.title ?? 'Scree Plot with Parallel Analysis', subtitleParts.join('; '), W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const eigenvalues = data.eigenvalues
  const comps = eigenvalues.map((_, i) => i + 1)
  const maxEv = Math.max(...eigenvalues, ...(diag?.parallelSimulated ?? []), 1)

  const xScale = d3.scaleBand<number>().domain(comps).range([0, width]).padding(0.25)
  const yScale = d3.scaleLinear().domain([0, maxEv * 1.15]).range([height, 0]).nice()

  // Gridlines
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('x1', 0).attr('x2', width).attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Kaiser criterion line (eigenvalue = 1)
  g.append('line').attr('x1', 0).attr('x2', width)
    .attr('y1', yScale(1)).attr('y2', yScale(1))
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '6,3').attr('stroke-width', 1.5)
  g.append('text').attr('x', width + 4).attr('y', yScale(1) + 4)
    .attr('font-size', 9).attr('fill', theme.textMuted).text('Kaiser')

  // Eigenvalue bars
  eigenvalues.forEach((ev, i) => {
    const x = xScale(i + 1) ?? 0
    const barColor = i < data.nFactors ? getColor(0, theme) : theme.gridLine
    g.append('rect')
      .attr('x', x).attr('y', yScale(ev))
      .attr('width', xScale.bandwidth()).attr('height', height - yScale(ev))
      .attr('fill', barColor).attr('opacity', 0.8).attr('rx', 2)
      .on('mouseover', (event: MouseEvent) => {
        const pctVar = eigenvalues.length > 0
          ? (ev / eigenvalues.reduce((s, v) => s + v, 0) * 100).toFixed(1)
          : '0'
        showTooltip(event, [
          formatTooltipRow('Component', `${i + 1}`),
          formatTooltipRow('Eigenvalue', ev),
          formatTooltipRow('% Variance', `${pctVar}%`),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Eigenvalue label on bar
    g.append('text')
      .attr('x', x + xScale.bandwidth() / 2).attr('y', yScale(ev) - 4)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text).text(ev.toFixed(2))
  })

  // Parallel analysis line (simulated 95th percentile)
  if (config.showParallelAnalysis !== false && diag?.parallelSimulated) {
    const sim = diag.parallelSimulated
    const lineData = comps.map((c, i) => ({
      x: (xScale(c) ?? 0) + xScale.bandwidth() / 2,
      y: yScale(sim[i] ?? 0),
    }))

    const line = d3.line<{ x: number; y: number }>().x(d => d.x).y(d => d.y)
    g.append('path').datum(lineData)
      .attr('d', line).attr('fill', 'none')
      .attr('stroke', getColor(2, theme)).attr('stroke-width', 2)
      .attr('stroke-dasharray', '5,3')

    lineData.forEach((pt) => {
      g.append('circle').attr('cx', pt.x).attr('cy', pt.y)
        .attr('r', 3).attr('fill', getColor(2, theme))
    })

    // Legend
    const legY = 8
    g.append('line').attr('x1', width - 120).attr('x2', width - 90)
      .attr('y1', legY).attr('y2', legY)
      .attr('stroke', getColor(2, theme)).attr('stroke-width', 2).attr('stroke-dasharray', '5,3')
    g.append('text').attr('x', width - 86).attr('y', legY + 4)
      .attr('font-size', 9).attr('fill', theme.textMuted).text('Simulated (95th)')
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Factor')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Eigenvalue')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Plot 2: Loadings Heatmap ────────────────────────────────────────────

function renderLoadings(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const threshold = config.loadingThreshold ?? 0.3
  const nItems = data.loadings.length
  const nFactors = data.nFactors

  // Sort items by primary factor assignment
  const itemOrder = Array.from({ length: nItems }, (_, i) => i)
  itemOrder.sort((a, b) => {
    const rowA = data.loadings[a]!
    const rowB = data.loadings[b]!
    const maxFactorA = rowA.indexOf(Math.max(...rowA.map(Math.abs)))
    const maxFactorB = rowB.indexOf(Math.max(...rowB.map(Math.abs)))
    if (maxFactorA !== maxFactorB) return maxFactorA - maxFactorB
    return Math.max(...rowB.map(Math.abs)) - Math.max(...rowA.map(Math.abs))
  })

  const cellW = 56, cellH = 28
  const showComm = true
  const extraCols = showComm ? 1 : 0
  const W = config.width ?? ((nFactors + extraCols) * cellW + 140)
  const H = config.height ?? (nItems * cellH + 100)
  const margin = { top: 60, right: 20, bottom: 30, left: 110 }

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'Factor Loadings', data.formatted, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu)
  const commScale = d3.scaleSequential().domain([0, 1]).interpolator(d3.interpolateGreens)

  // Variable labels
  const varLabels = config.variableLabels ?? data.variableNames
  const facLabels = config.factorLabels ?? data.factorNames

  // Column headers
  facLabels.forEach((label, ci) => {
    g.append('text').attr('x', ci * cellW + cellW / 2).attr('y', -8)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall)
      .attr('font-weight', '600').attr('fill', theme.text).text(label)
  })
  if (showComm) {
    g.append('text').attr('x', nFactors * cellW + cellW / 2).attr('y', -8)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall)
      .attr('font-weight', '600').attr('fill', theme.text).text('h²')
  }

  // Cells
  itemOrder.forEach((origIdx, rowIdx) => {
    const row = data.loadings[origIdx]!
    const label = varLabels[origIdx] ?? `V${origIdx + 1}`

    // Row label
    g.append('text').attr('x', -6).attr('y', rowIdx * cellH + cellH / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(label)

    // Loading cells
    row.forEach((val, ci) => {
      const dimmed = Math.abs(val) < threshold
      g.append('rect')
        .attr('x', ci * cellW).attr('y', rowIdx * cellH)
        .attr('width', cellW - 2).attr('height', cellH - 2).attr('rx', 3)
        .attr('fill', colorScale(val)).attr('opacity', dimmed ? 0.3 : 1)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow(label, ''),
            formatTooltipRow(facLabels[ci] ?? `F${ci + 1}`, val),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)

      g.append('text')
        .attr('x', ci * cellW + cellW / 2).attr('y', rowIdx * cellH + cellH / 2 + 4)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', Math.min(10, cellW / 5))
        .attr('fill', Math.abs(val) > 0.6 ? '#fff' : theme.text)
        .attr('opacity', dimmed ? 0.5 : 1)
        .text(val.toFixed(2))
    })

    // Communality cell
    if (showComm) {
      const comm = data.communalities[origIdx] ?? 0
      g.append('rect')
        .attr('x', nFactors * cellW).attr('y', rowIdx * cellH)
        .attr('width', cellW - 2).attr('height', cellH - 2).attr('rx', 3)
        .attr('fill', commScale(comm))

      g.append('text')
        .attr('x', nFactors * cellW + cellW / 2).attr('y', rowIdx * cellH + cellH / 2 + 4)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', Math.min(10, cellW / 5))
        .attr('fill', comm > 0.6 ? '#fff' : theme.text)
        .text(comm.toFixed(2))
    }
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Plot 3: CFA/SEM Path Diagram ───────────────────────────────────────

/** Resolve FAPathStyle with defaults. */
function resolvePathStyle(s?: FAPathStyle): Required<FAPathStyle> {
  return {
    factorRx: s?.factorRx ?? 56,
    factorRy: s?.factorRy ?? 28,
    factorStroke: s?.factorStroke ?? 2.5,
    factorFontSize: s?.factorFontSize ?? 14,
    factorFontWeight: s?.factorFontWeight ?? '800',
    factorGlow: s?.factorGlow ?? true,
    factorGlowRadius: s?.factorGlowRadius ?? 6,
    factorHighlight: s?.factorHighlight ?? true,
    itemWidth: s?.itemWidth ?? 96,
    itemHeight: s?.itemHeight ?? 30,
    itemRadius: s?.itemRadius ?? 5,
    itemStroke: s?.itemStroke ?? 1.2,
    itemFontSize: s?.itemFontSize ?? 11,
    itemFontWeight: s?.itemFontWeight ?? '600',
    itemAccentWidth: s?.itemAccentWidth ?? 4,
    itemGradient: s?.itemGradient ?? true,
    errorRx: s?.errorRx ?? 14,
    errorRy: s?.errorRy ?? 12,
    errorStroke: s?.errorStroke ?? 1.2,
    errorFontSize: s?.errorFontSize ?? 9,
    errorSelfLoop: s?.errorSelfLoop ?? true,
    errorSelfLoopSize: s?.errorSelfLoopSize ?? 12,
    arrowMinWidth: s?.arrowMinWidth ?? 0.8,
    arrowMaxWidth: s?.arrowMaxWidth ?? 2.5,
    arrowMinOpacity: s?.arrowMinOpacity ?? 0.4,
    arrowMaxOpacity: s?.arrowMaxOpacity ?? 0.95,
    arrowMarkerSize: s?.arrowMarkerSize ?? 6,
    arrowCurvature: s?.arrowCurvature ?? 0.4,
    loadingFontSize: s?.loadingFontSize ?? 9.5,
    loadingFontWeight: s?.loadingFontWeight ?? '700',
    loadingPillRadius: s?.loadingPillRadius ?? 7,
    loadingLabelOffset: s?.loadingLabelOffset ?? -6,
    loadingPosition: s?.loadingPosition ?? 0.6,
    halo: s?.halo ?? true,
    haloColor: s?.haloColor ?? '#ffffff',
    haloWidth: s?.haloWidth ?? 3.5,
    covColor: s?.covColor ?? '#8892a0',
    covMinWidth: s?.covMinWidth ?? 1.0,
    covMaxWidth: s?.covMaxWidth ?? 3.0,
    covMarkerSize: s?.covMarkerSize ?? 6,
    covLabelFontSize: s?.covLabelFontSize ?? 10,
    covNestSpacing: s?.covNestSpacing ?? 20,
    errorArrowWidth: s?.errorArrowWidth ?? 0.8,
    errorArrowMarkerSize: s?.errorArrowMarkerSize ?? 5,
    itemGap: s?.itemGap ?? 6,
    factorGroupGap: s?.factorGroupGap ?? 32,
    arrowSpan: s?.arrowSpan ?? 120,
    errorSpan: s?.errorSpan ?? 36,
    covArcReserve: s?.covArcReserve ?? 70,
    topPadding: s?.topPadding ?? 74,
    bottomPadding: s?.bottomPadding ?? 44,
    rightPadding: s?.rightPadding ?? 50,
    factorColors: s?.factorColors ?? [],
    fitCardWidth: s?.fitCardWidth ?? 210,
    fitCardHeight: s?.fitCardHeight ?? 118,
    crossLoadingThreshold: s?.crossLoadingThreshold ?? 0.2,
    crossLoadingDash: s?.crossLoadingDash ?? '5,4',
    crossLoadingOpacity: s?.crossLoadingOpacity ?? 0.3,
    showGroupBrackets: s?.showGroupBrackets ?? true,
    groupBracketOpacity: s?.groupBracketOpacity ?? 0.12,
    showShadows: s?.showShadows ?? true,
  }
}

/** Get factor color: custom array → theme palette fallback. */
function fColor(fi: number, ps: Required<FAPathStyle>, theme: CarmTheme): string {
  if (ps.factorColors.length > 0) return ps.factorColors[fi % ps.factorColors.length]!
  return getColor(fi, theme)
}

/**
 * Enterprise-grade AMOS-style SEM path diagram.
 * Every visual parameter is controllable via FAPathStyle.
 */
function renderPath(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult | CFAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const ps = resolvePathStyle(config.pathStyle)
  const showErrors = config.showErrorTerms !== false
  const showFit = config.showFitBox !== false

  const nFactors = data.nFactors
  const nItems = data.loadings.length
  const varLabels = config.variableLabels ?? data.variableNames
  const facLabels = config.factorLabels ?? data.factorNames
  const hasCFA = isCFA(data)

  // ── Assign items to primary factor ─────────────────────────────
  const factorItems: number[][] = Array.from({ length: nFactors }, () => [])
  if (hasCFA) {
    const cfaModel = (data as CFAResult).model
    for (let fi = 0; fi < nFactors; fi++) {
      const key = facLabels[fi]!
      const items = cfaModel[key]
      if (items) {
        for (let j = 0; j < items.length; j++) factorItems[fi]!.push(items[j]!)
      }
    }
  } else {
    for (let i = 0; i < nItems; i++) {
      const row = data.loadings[i]!
      let maxF = 0, maxVal = 0
      for (let f = 0; f < nFactors; f++) {
        if (Math.abs(row[f]!) > maxVal) { maxVal = Math.abs(row[f]!); maxF = f }
      }
      factorItems[maxF]!.push(i)
    }
  }

  // ── Layout from pathStyle ──────────────────────────────────────
  const groupHeights = factorItems.map(items =>
    Math.max(items.length * (ps.itemHeight + ps.itemGap) - ps.itemGap, ps.factorRy * 2)
  )
  const totalGroupH = groupHeights.reduce((s, h) => s + h, 0) + (nFactors - 1) * ps.factorGroupGap

  const bottomPad = ps.bottomPadding

  // Column X positions
  const factorCx = ps.covArcReserve + ps.factorRx + 16
  const itemLeft = factorCx + ps.factorRx + ps.arrowSpan
  const itemRight = itemLeft + ps.itemWidth
  const errorCx = showErrors ? itemRight + ps.errorSpan + ps.errorRx : itemRight + 10

  const selfLoopExtra = (showErrors && ps.errorSelfLoop) ? ps.errorSelfLoopSize + 30 : 0
  const W = config.width ?? (errorCx + (showErrors ? ps.errorRx : 0) + selfLoopExtra + ps.rightPadding)
  const H = config.height ?? (ps.topPadding + totalGroupH + bottomPad)

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H)
    .style('background', theme.background)
    .attr('xmlns', 'http://www.w3.org/2000/svg')

  addSubtitle(
    svg,
    config.title ?? (hasCFA ? 'Confirmatory Factor Analysis — Path Diagram' : 'Exploratory Factor Analysis — Path Diagram'),
    data.formatted, W, theme
  )

  const defs = svg.append('defs')
  const g = svg.append('g')

  // ── SVG Defs: filters ──────────────────────────────────────────

  if (ps.factorGlow) {
    const glowFilter = defs.append('filter').attr('id', 'fa-glow')
      .attr('x', '-40%').attr('y', '-40%').attr('width', '180%').attr('height', '180%')
    glowFilter.append('feGaussianBlur').attr('in', 'SourceGraphic').attr('stdDeviation', ps.factorGlowRadius).attr('result', 'blur')
    glowFilter.append('feColorMatrix').attr('in', 'blur').attr('type', 'saturate').attr('values', '0.35').attr('result', 'desatBlur')
    const glowMerge = glowFilter.append('feMerge')
    glowMerge.append('feMergeNode').attr('in', 'desatBlur')
    glowMerge.append('feMergeNode').attr('in', 'SourceGraphic')
  }

  if (ps.showShadows) {
    const shadowMed = defs.append('filter').attr('id', 'fa-shadow-med')
      .attr('x', '-20%').attr('y', '-15%').attr('width', '140%').attr('height', '140%')
    shadowMed.append('feDropShadow').attr('dx', 0).attr('dy', 2).attr('stdDeviation', 4).attr('flood-color', 'rgba(0,0,0,0.10)')

    const shadowSm = defs.append('filter').attr('id', 'fa-shadow-sm')
      .attr('x', '-10%').attr('y', '-10%').attr('width', '120%').attr('height', '130%')
    shadowSm.append('feDropShadow').attr('dx', 0).attr('dy', 1).attr('stdDeviation', 2).attr('flood-color', 'rgba(0,0,0,0.07)')
  }

  // ── SVG Defs: gradients ────────────────────────────────────────

  for (let f = 0; f < nFactors; f++) {
    const color = fColor(f, ps, theme)
    const grad = defs.append('radialGradient').attr('id', `fa-fgrad-${f}`)
      .attr('cx', '35%').attr('cy', '28%').attr('r', '70%')
    grad.append('stop').attr('offset', '0%').attr('stop-color', '#fff').attr('stop-opacity', 0.55)
    grad.append('stop').attr('offset', '35%').attr('stop-color', color).attr('stop-opacity', 0.18)
    grad.append('stop').attr('offset', '100%').attr('stop-color', color).attr('stop-opacity', 0.06)

    if (ps.itemGradient) {
      const lg = defs.append('linearGradient').attr('id', `fa-igrad-${f}`)
        .attr('x1', '0%').attr('y1', '0%').attr('x2', '100%').attr('y2', '0%')
      lg.append('stop').attr('offset', '0%').attr('stop-color', color).attr('stop-opacity', 0.08)
      lg.append('stop').attr('offset', '100%').attr('stop-color', color).attr('stop-opacity', 0.0)
    }
  }

  // ── SVG Defs: markers ──────────────────────────────────────────

  const mw = ps.arrowMarkerSize, mh = ps.arrowMarkerSize * 0.7
  for (let f = 0; f < nFactors; f++) {
    const color = fColor(f, ps, theme)
    defs.append('marker')
      .attr('id', `fa-arr-${f}`)
      .attr('viewBox', '0 0 10 7').attr('refX', 9.5).attr('refY', 3.5)
      .attr('markerWidth', mw).attr('markerHeight', mh).attr('orient', 'auto')
      .append('path').attr('d', 'M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z').attr('fill', color)
  }

  defs.append('marker')
    .attr('id', 'fa-arr-ns')
    .attr('viewBox', '0 0 10 7').attr('refX', 9.5).attr('refY', 3.5)
    .attr('markerWidth', mw).attr('markerHeight', mh).attr('orient', 'auto')
    .append('path').attr('d', 'M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z').attr('fill', theme.textMuted)

  // Covariance markers
  const cmw = ps.covMarkerSize, cmh = ps.covMarkerSize * 0.7
  defs.append('marker')
    .attr('id', 'fa-cov-s')
    .attr('viewBox', '0 0 10 7').attr('refX', 1).attr('refY', 3.5)
    .attr('markerWidth', cmw).attr('markerHeight', cmh).attr('orient', 'auto')
    .append('path').attr('d', 'M10,0.5 L1,3.5 L10,6.5 L8,3.5 Z').attr('fill', ps.covColor)
  defs.append('marker')
    .attr('id', 'fa-cov-e')
    .attr('viewBox', '0 0 10 7').attr('refX', 9).attr('refY', 3.5)
    .attr('markerWidth', cmw).attr('markerHeight', cmh).attr('orient', 'auto')
    .append('path').attr('d', 'M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z').attr('fill', ps.covColor)

  // Error arrow
  const emw = ps.errorArrowMarkerSize, emh = ps.errorArrowMarkerSize * 0.7
  defs.append('marker')
    .attr('id', 'fa-arr-err')
    .attr('viewBox', '0 0 10 7').attr('refX', 9).attr('refY', 3.5)
    .attr('markerWidth', emw).attr('markerHeight', emh).attr('orient', 'auto')
    .append('path').attr('d', 'M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z').attr('fill', theme.axisLine)

  // ── Compute Y positions ────────────────────────────────────────
  const factorCy: number[] = []
  const itemCy: number[] = new Array(nItems).fill(0)
  let yAccum = ps.topPadding

  for (let fi = 0; fi < nFactors; fi++) {
    const items = factorItems[fi]!
    const gh = groupHeights[fi]!
    factorCy.push(yAccum + gh / 2)
    for (let j = 0; j < items.length; j++) {
      itemCy[items[j]!] = yAccum + j * (ps.itemHeight + ps.itemGap) + ps.itemHeight / 2
    }
    yAccum += gh + ps.factorGroupGap
  }

  // ── Layer 0: Group bracket lines ───────────────────────────────
  if (ps.showGroupBrackets) {
    for (let fi = 0; fi < nFactors; fi++) {
      const items = factorItems[fi]!
      if (items.length < 2) continue
      const firstY = itemCy[items[0]!]!
      const lastY = itemCy[items[items.length - 1]!]!
      g.append('line')
        .attr('x1', itemLeft - 6).attr('y1', firstY)
        .attr('x2', itemLeft - 6).attr('y2', lastY)
        .attr('stroke', fColor(fi, ps, theme)).attr('stroke-width', 1.5)
        .attr('opacity', ps.groupBracketOpacity).attr('stroke-linecap', 'round')
    }
  }

  // ── Layer 1: Factor covariance arcs (left side) ────────────────
  const Phi = data.factorCorrelations
  for (let fi = 0; fi < nFactors; fi++) {
    for (let fj = fi + 1; fj < nFactors; fj++) {
      const corr = Phi[fi]![fj]!
      if (Math.abs(corr) < 0.01) continue

      const y1 = factorCy[fi]!, y2 = factorCy[fj]!
      const span = Math.abs(y2 - y1)
      const nestLevel = fj - fi
      const arcX = factorCx - ps.factorRx - 14 - nestLevel * ps.covNestSpacing - span * 0.04
      const absCorr = Math.abs(corr)
      const strokeW = ps.covMinWidth + absCorr * (ps.covMaxWidth - ps.covMinWidth)

      const path = `M${factorCx - ps.factorRx - 2},${y1} C${arcX},${y1} ${arcX},${y2} ${factorCx - ps.factorRx - 2},${y2}`
      g.append('path')
        .attr('d', path).attr('fill', 'none')
        .attr('stroke', ps.covColor).attr('stroke-width', strokeW)
        .attr('stroke-linecap', 'round')
        .attr('marker-start', 'url(#fa-cov-s)')
        .attr('marker-end', 'url(#fa-cov-e)')

      // Pill label at Bézier t=0.5
      const pillW = 38, pillH = 17
      const t = 0.5, t1 = 0.5
      const labelX = t1 * t1 * t1 * (factorCx - ps.factorRx - 2) + 3 * t1 * t1 * t * arcX + 3 * t1 * t * t * arcX + t * t * t * (factorCx - ps.factorRx - 2)
      const labelY = t1 * t1 * t1 * y1 + 3 * t1 * t1 * t * y1 + 3 * t1 * t * t * y2 + t * t * t * y2
      g.append('rect')
        .attr('x', labelX - pillW / 2).attr('y', labelY - pillH / 2)
        .attr('width', pillW).attr('height', pillH).attr('rx', 8)
        .attr('fill', theme.background).attr('stroke', theme.gridLine).attr('stroke-width', 0.7)
      g.append('text')
        .attr('x', labelX).attr('y', labelY + 3.5)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', ps.covLabelFontSize).attr('font-weight', '600')
        .attr('fill', ps.covColor)
        .text(corr.toFixed(2))
    }
  }

  // ── Layer 2: Loading arrows (factor → item) ────────────────────
  for (let fi = 0; fi < nFactors; fi++) {
    const items = factorItems[fi]!
    const color = fColor(fi, ps, theme)

    for (let idx = 0; idx < items.length; idx++) {
      const itemIdx = items[idx]!
      if (itemIdx >= nItems) continue
      const loading = data.loadings[itemIdx]![fi]!
      if (Math.abs(loading) < 0.01) continue

      let sig = true
      let stars = ''
      if (hasCFA) {
        const pe = (data as CFAResult).parameterEstimates.loadings[fi]?.[idx]
        if (pe) {
          sig = pe.pValue < 0.05
          stars = pe.pValue < 0.001 ? '***' : pe.pValue < 0.01 ? '**' : pe.pValue < 0.05 ? '*' : ''
        }
      }

      const x1 = factorCx + ps.factorRx + 2
      const y1 = factorCy[fi]!
      const x2 = itemLeft - 2
      const y2 = itemCy[itemIdx]!
      const dy = y2 - y1

      // S-shaped Bézier controlled by arrowCurvature
      const curv = ps.arrowCurvature
      const cp1x = x1 + ps.arrowSpan * curv, cp1y = y1 + dy * 0.05
      const cp2x = x2 - ps.arrowSpan * (curv * 0.85), cp2y = y2

      const absL = Math.abs(loading)
      const strokeW = ps.arrowMinWidth + absL * (ps.arrowMaxWidth - ps.arrowMinWidth)
      const opacity = ps.arrowMinOpacity + absL * (ps.arrowMaxOpacity - ps.arrowMinOpacity)
      const arrowColor = sig ? color : theme.textMuted
      const markerId = sig ? `fa-arr-${fi}` : 'fa-arr-ns'

      g.append('path')
        .attr('d', `M${x1},${y1} C${cp1x},${cp1y} ${cp2x},${cp2y} ${x2},${y2}`)
        .attr('fill', 'none')
        .attr('stroke', arrowColor).attr('stroke-width', strokeW)
        .attr('stroke-linecap', 'round').attr('opacity', opacity)
        .attr('marker-end', `url(#${markerId})`)

      // Loading label — clean number directly on the arrow (AMOS-style)
      const mt = ps.loadingPosition, mt1 = 1 - mt
      const lx = mt1 * mt1 * mt1 * x1 + 3 * mt1 * mt1 * mt * cp1x + 3 * mt1 * mt * mt * cp2x + mt * mt * mt * x2
      const ly = mt1 * mt1 * mt1 * y1 + 3 * mt1 * mt1 * mt * cp1y + 3 * mt1 * mt * mt * cp2y + mt * mt * mt * y2

      const labelText = `${loading.toFixed(2)}${stars}`

      const lbl = g.append('text')
        .attr('x', lx).attr('y', ly + ps.loadingLabelOffset)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', ps.loadingFontSize)
        .attr('font-weight', ps.loadingFontWeight)
        .attr('fill', arrowColor)
        .text(labelText)
      if (ps.halo) lbl.attr('paint-order', 'stroke').attr('stroke', ps.haloColor).attr('stroke-width', ps.haloWidth).attr('stroke-linejoin', 'round')
    }
  }

  // ── Layer 2b: Cross-loadings (EFA only, dashed) ────────────────
  if (!hasCFA) {
    for (let i = 0; i < nItems; i++) {
      const row = data.loadings[i]!
      let primaryF = 0, maxVal = 0
      for (let f = 0; f < nFactors; f++) {
        if (Math.abs(row[f]!) > maxVal) { maxVal = Math.abs(row[f]!); primaryF = f }
      }
      for (let f = 0; f < nFactors; f++) {
        if (f === primaryF) continue
        if (Math.abs(row[f]!) < ps.crossLoadingThreshold) continue
        const x1 = factorCx + ps.factorRx + 2, y1 = factorCy[f]!
        const x2 = itemLeft - 2, y2 = itemCy[i]!
        const dy = y2 - y1
        const curv = ps.arrowCurvature
        g.append('path')
          .attr('d', `M${x1},${y1} C${x1 + ps.arrowSpan * curv},${y1 + dy * 0.05} ${x2 - ps.arrowSpan * curv * 0.85},${y2} ${x2},${y2}`)
          .attr('fill', 'none')
          .attr('stroke', theme.textMuted).attr('stroke-width', 0.9)
          .attr('stroke-dasharray', ps.crossLoadingDash).attr('opacity', ps.crossLoadingOpacity)
      }
    }
  }

  // ── Layer 3: Factor ellipses ───────────────────────────────────
  for (let fi = 0; fi < nFactors; fi++) {
    const cy = factorCy[fi]!
    const color = fColor(fi, ps, theme)

    if (ps.factorGlow) {
      g.append('ellipse').attr('cx', factorCx).attr('cy', cy)
        .attr('rx', ps.factorRx + 4).attr('ry', ps.factorRy + 3)
        .attr('fill', color).attr('opacity', 0.06)
        .attr('filter', 'url(#fa-glow)')
    }

    if (ps.showShadows) {
      g.append('ellipse').attr('cx', factorCx).attr('cy', cy)
        .attr('rx', ps.factorRx).attr('ry', ps.factorRy)
        .attr('fill', 'none').attr('filter', 'url(#fa-shadow-med)')
        .attr('stroke', 'transparent')
    }

    g.append('ellipse').attr('cx', factorCx).attr('cy', cy)
      .attr('rx', ps.factorRx).attr('ry', ps.factorRy)
      .attr('fill', `url(#fa-fgrad-${fi})`)
      .attr('stroke', color).attr('stroke-width', ps.factorStroke)

    if (ps.factorHighlight && ps.factorRx > 20 && ps.factorRy > 15) {
      g.append('ellipse').attr('cx', factorCx).attr('cy', cy - 4)
        .attr('rx', ps.factorRx - 10).attr('ry', ps.factorRy - 10)
        .attr('fill', 'none')
        .attr('stroke', '#fff').attr('stroke-width', 1).attr('opacity', 0.25)
    }

    g.append('text').attr('x', factorCx).attr('y', cy + 1)
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', ps.factorFontSize)
      .attr('font-weight', ps.factorFontWeight).attr('letter-spacing', '-0.4')
      .attr('fill', theme.text).text(facLabels[fi] ?? `F${fi + 1}`)
  }

  // ── Layer 4: Item rectangles ───────────────────────────────────
  for (let i = 0; i < nItems; i++) {
    const cy = itemCy[i]!
    let primaryFactor = 0, primVal = 0
    for (let f = 0; f < nFactors; f++) {
      if (Math.abs(data.loadings[i]![f]!) > primVal) { primVal = Math.abs(data.loadings[i]![f]!); primaryFactor = f }
    }
    const accentColor = fColor(primaryFactor, ps, theme)

    if (ps.showShadows) {
      g.append('rect')
        .attr('x', itemLeft).attr('y', cy - ps.itemHeight / 2)
        .attr('width', ps.itemWidth).attr('height', ps.itemHeight).attr('rx', ps.itemRadius)
        .attr('fill', theme.surface).attr('filter', 'url(#fa-shadow-sm)')
    }

    g.append('rect')
      .attr('x', itemLeft).attr('y', cy - ps.itemHeight / 2)
      .attr('width', ps.itemWidth).attr('height', ps.itemHeight).attr('rx', ps.itemRadius)
      .attr('fill', ps.itemGradient ? `url(#fa-igrad-${primaryFactor})` : theme.background)
      .attr('stroke', accentColor).attr('stroke-width', ps.itemStroke).attr('stroke-opacity', 0.5)

    if (ps.itemAccentWidth > 0) {
      g.append('rect')
        .attr('x', itemLeft).attr('y', cy - ps.itemHeight / 2)
        .attr('width', ps.itemAccentWidth).attr('height', ps.itemHeight)
        .attr('rx', Math.min(ps.itemAccentWidth, ps.itemRadius))
        .attr('fill', accentColor).attr('opacity', 0.85)
    }

    g.append('text')
      .attr('x', itemLeft + ps.itemWidth / 2 + ps.itemAccentWidth / 2).attr('y', cy + 1)
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', ps.itemFontSize)
      .attr('font-weight', ps.itemFontWeight).attr('fill', theme.text)
      .text(varLabels[i] ?? `V${i + 1}`)

    // Tooltip hitbox
    g.append('rect')
      .attr('x', itemLeft).attr('y', cy - ps.itemHeight / 2)
      .attr('width', ps.itemWidth).attr('height', ps.itemHeight)
      .attr('fill', 'transparent').attr('cursor', 'default')
      .on('mouseover', (event: MouseEvent) => {
        const rows = facLabels.map((fl, fi) =>
          formatTooltipRow(fl, data.loadings[i]![fi] ?? 0)
        ).join('')
        showTooltip(event, [
          formatTooltipRow(varLabels[i] ?? `V${i + 1}`, ''),
          rows,
          formatTooltipRow('h²', data.communalities[i] ?? 0),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)
  }

  // ── Layer 5: Error ellipses ────────────────────────────────────
  if (showErrors) {
    for (let i = 0; i < nItems; i++) {
      const cy = itemCy[i]!
      const uniq = data.uniqueness[i] ?? 0

      // Arrow: error → item (fixed to 1 in CFA for identification)
      const errArrX1 = errorCx - ps.errorRx - 1
      const errArrX2 = itemRight + 2
      g.append('line')
        .attr('x1', errArrX1).attr('y1', cy)
        .attr('x2', errArrX2).attr('y2', cy)
        .attr('stroke', theme.axisLine).attr('stroke-width', ps.errorArrowWidth)
        .attr('marker-end', 'url(#fa-arr-err)')

      // "1" label on the error→item arrow (fixed loading)
      const errMidX = (errArrX1 + errArrX2) / 2
      const oneLbl = g.append('text')
        .attr('x', errMidX).attr('y', cy - 5)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', ps.errorFontSize)
        .attr('font-weight', '600').attr('fill', theme.textMuted)
        .text('1')
      if (ps.halo) oneLbl.attr('paint-order', 'stroke').attr('stroke', ps.haloColor).attr('stroke-width', ps.haloWidth).attr('stroke-linejoin', 'round')

      g.append('ellipse').attr('cx', errorCx).attr('cy', cy)
        .attr('rx', ps.errorRx).attr('ry', ps.errorRy)
        .attr('fill', theme.surface).attr('stroke', theme.axisLine).attr('stroke-width', ps.errorStroke)

      // Self-loop on the outside right of the error ellipse
      if (ps.errorSelfLoop) {
        const loopRight = errorCx + ps.errorRx + ps.errorSelfLoopSize
        g.append('path')
          .attr('d', `M${errorCx + ps.errorRx},${cy - 3} C${loopRight},${cy - 14} ${loopRight},${cy + 14} ${errorCx + ps.errorRx},${cy + 3}`)
          .attr('fill', 'none')
          .attr('stroke', theme.axisLine).attr('stroke-width', ps.errorArrowWidth)
          .attr('marker-end', 'url(#fa-arr-err)')

        // Uniqueness label next to the self-loop
        const errLbl = g.append('text')
          .attr('x', loopRight + 2).attr('y', cy + ps.errorFontSize * 0.35)
          .attr('text-anchor', 'start')
          .attr('font-family', theme.fontFamilyMono).attr('font-size', ps.errorFontSize - 0.5)
          .attr('font-weight', '500').attr('fill', theme.textMuted)
          .text(uniq.toFixed(2))
        if (ps.halo) errLbl.attr('paint-order', 'stroke').attr('stroke', ps.haloColor).attr('stroke-width', ps.haloWidth).attr('stroke-linejoin', 'round')
      }

      g.append('text').attr('x', errorCx).attr('y', cy + 1)
        .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
        .attr('font-family', theme.fontFamilyMono).attr('font-size', ps.errorFontSize)
        .attr('font-weight', '600').attr('fill', theme.textMuted)
        .text(`e${i + 1}`)
    }
  }

  // ── Layer 6: Fit indices — publication table below diagram ─────
  if (showFit) {
    const fit = data.fit
    const fitColor = (good: boolean, ok: boolean) =>
      good ? '#2e7d32' : ok ? '#e65100' : '#c62828'

    const table = d3.select(container).append('table')
      .style('border-collapse', 'collapse')
      .style('margin', '12px auto 4px')
      .style('font-family', theme.fontFamilyMono)
      .style('font-size', '12px')
      .style('color', theme.text)
      .style('min-width', '420px')

    // Header
    const thead = table.append('thead')
    const headerRow = thead.append('tr')
      .style('border-bottom', `2px solid ${theme.text}`)
    const headers = ['Index', 'Value', '90% CI', 'Verdict']
    headers.forEach(h => {
      headerRow.append('th')
        .style('padding', '6px 14px')
        .style('text-align', h === 'Index' ? 'left' : 'right')
        .style('font-weight', '700')
        .style('font-size', '11px')
        .style('letter-spacing', '0.5px')
        .style('color', theme.textMuted)
        .text(h)
    })

    const tbody = table.append('tbody')

    const fitRows: Array<{
      label: string; value: string; ci: string; verdict: string; color: string
    }> = [
      {
        label: 'χ²',
        value: `${formatStat(fit.chiSq)} (df = ${fit.df})`,
        ci: '',
        verdict: `${formatP(fit.pValue)}`,
        color: theme.textAnnotation,
      },
      {
        label: 'RMSEA',
        value: formatStat(fit.rmsea, 3),
        ci: `[${formatStat(fit.rmseaCI[0], 3)}, ${formatStat(fit.rmseaCI[1], 3)}]`,
        verdict: fit.rmsea <= 0.05 ? 'Good' : fit.rmsea <= 0.08 ? 'Acceptable' : 'Poor',
        color: fitColor(fit.rmsea <= 0.05, fit.rmsea <= 0.08),
      },
      {
        label: 'CFI',
        value: formatStat(fit.cfi, 3),
        ci: '',
        verdict: fit.cfi >= 0.95 ? 'Good' : fit.cfi >= 0.90 ? 'Acceptable' : 'Poor',
        color: fitColor(fit.cfi >= 0.95, fit.cfi >= 0.90),
      },
      {
        label: 'TLI',
        value: formatStat(fit.tli, 3),
        ci: '',
        verdict: fit.tli >= 0.95 ? 'Good' : fit.tli >= 0.90 ? 'Acceptable' : 'Poor',
        color: fitColor(fit.tli >= 0.95, fit.tli >= 0.90),
      },
      {
        label: 'SRMR',
        value: formatStat(fit.srmr, 3),
        ci: '',
        verdict: fit.srmr <= 0.05 ? 'Good' : fit.srmr <= 0.08 ? 'Acceptable' : 'Poor',
        color: fitColor(fit.srmr <= 0.05, fit.srmr <= 0.08),
      },
      {
        label: 'AIC',
        value: formatStat(fit.aic, 1),
        ci: '',
        verdict: '',
        color: theme.textAnnotation,
      },
      {
        label: 'BIC',
        value: formatStat(fit.bic, 1),
        ci: '',
        verdict: '',
        color: theme.textAnnotation,
      },
    ]

    fitRows.forEach((row, i) => {
      const tr = tbody.append('tr')
        .style('border-bottom', i < fitRows.length - 1 ? `1px solid ${theme.gridLine}` : 'none')

      tr.append('td')
        .style('padding', '5px 14px')
        .style('font-weight', '600')
        .style('white-space', 'nowrap')
        .text(row.label)

      tr.append('td')
        .style('padding', '5px 14px')
        .style('text-align', 'right')
        .style('font-weight', '500')
        .text(row.value)

      tr.append('td')
        .style('padding', '5px 14px')
        .style('text-align', 'right')
        .style('color', theme.textMuted)
        .style('font-size', '11px')
        .text(row.ci)

      const verdictCell = tr.append('td')
        .style('padding', '5px 14px')
        .style('text-align', 'right')
        .style('font-weight', '700')
        .style('color', row.color)

      if (row.verdict) {
        verdictCell.append('span')
          .style('display', 'inline-block')
          .style('width', '8px')
          .style('height', '8px')
          .style('border-radius', '50%')
          .style('background', row.color)
          .style('margin-right', '6px')
          .style('vertical-align', 'middle')
        verdictCell.append('span')
          .style('vertical-align', 'middle')
          .text(row.verdict)
      }
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Plot 4: Communality Bar Chart ──────────────────────────────────────

function renderCommunality(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const nItems = data.communalities.length
  const W = config.width ?? 500
  const H = config.height ?? Math.max(nItems * 28 + 120, 250)
  const margin = { top: theme.marginTop, right: 70, bottom: theme.marginBottom, left: 110 }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'Communalities', data.formatted, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const varLabels = config.variableLabels ?? data.variableNames

  // Sort by communality descending
  const indices = Array.from({ length: nItems }, (_, i) => i)
  indices.sort((a, b) => (data.communalities[b] ?? 0) - (data.communalities[a] ?? 0))

  const names = indices.map(i => varLabels[i] ?? `V${i + 1}`)
  const values = indices.map(i => data.communalities[i] ?? 0)

  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width]).nice()
  const yScale = d3.scaleBand<string>().domain(names).range([0, height]).padding(0.25)

  // Color scale: red (0) → amber (0.4) → green (1.0)
  const commColor = (v: number) => {
    if (v < 0.3) return '#e15759'
    if (v < 0.5) return '#f28e2b'
    return '#59a14f'
  }

  // Gridlines
  g.selectAll('.grid').data(xScale.ticks(5)).join('line')
    .attr('x1', d => xScale(d)).attr('x2', d => xScale(d))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Reference line at 0.4
  g.append('line')
    .attr('x1', xScale(0.4)).attr('x2', xScale(0.4))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '6,3').attr('stroke-width', 1.5)

  // Bars
  values.forEach((val, i) => {
    const y = yScale(names[i]!) ?? 0
    g.append('rect')
      .attr('x', 0).attr('y', y)
      .attr('width', xScale(val)).attr('height', yScale.bandwidth())
      .attr('fill', commColor(val)).attr('opacity', 0.8).attr('rx', 2)
      .on('mouseover', (event: MouseEvent) => {
        const origIdx = indices[i]!
        showTooltip(event, [
          formatTooltipRow('Variable', varLabels[origIdx] ?? `V${origIdx + 1}`),
          formatTooltipRow('Communality', val),
          formatTooltipRow('Uniqueness', data.uniqueness[origIdx] ?? 0),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Value at bar end
    g.append('text')
      .attr('x', xScale(val) + 4).attr('y', y + yScale.bandwidth() / 2 + 4)
      .attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text).text(val.toFixed(2))
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Communality (h²)')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Plot 5: Factor Correlation Matrix ──────────────────────────────────

function renderFactorCorrelation(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const Phi = data.factorCorrelations
  const k = data.nFactors
  const facLabels = config.factorLabels ?? data.factorNames

  // Check if orthogonal (identity matrix)
  let isOrthogonal = true
  for (let i = 0; i < k && isOrthogonal; i++) {
    for (let j = 0; j < k && isOrthogonal; j++) {
      if (i !== j && Math.abs(Phi[i]![j]!) > 0.001) isOrthogonal = false
    }
  }

  const cellSize = 60
  const W = config.width ?? (k * cellSize + 140)
  const H = config.height ?? (k * cellSize + 120)

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'Factor Correlations (Φ)', data.rotation === 'varimax' ? 'Orthogonal rotation — Φ = I' : '', W, theme)

  if (isOrthogonal) {
    svg.append('text')
      .attr('x', W / 2).attr('y', H / 2)
      .attr('text-anchor', 'middle').attr('font-size', 14)
      .attr('fill', theme.textMuted)
      .text('Orthogonal rotation: all factor correlations = 0')
    return
  }

  const margin = { top: 60, right: 20, bottom: 20, left: 80 }
  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu)

  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) {
      const val = Phi[i]![j]!
      g.append('rect')
        .attr('x', j * cellSize).attr('y', i * cellSize)
        .attr('width', cellSize - 2).attr('height', cellSize - 2).attr('rx', 3)
        .attr('fill', i === j ? theme.surface : colorScale(val))
        .on('mouseover', (event: MouseEvent) => {
          if (i !== j) {
            showTooltip(event, [
              formatTooltipRow(`${facLabels[i]} × ${facLabels[j]}`, ''),
              formatTooltipRow('r', val),
            ].join(''), theme)
          }
        })
        .on('mouseout', hideTooltip)

      if (i !== j) {
        g.append('text')
          .attr('x', j * cellSize + cellSize / 2).attr('y', i * cellSize + cellSize / 2 + 4)
          .attr('text-anchor', 'middle')
          .attr('font-family', theme.fontFamilyMono).attr('font-size', 12)
          .attr('fill', Math.abs(val) > 0.6 ? '#fff' : theme.text)
          .text(val.toFixed(2))
      }
    }
  }

  // Labels
  facLabels.forEach((lbl, i) => {
    g.append('text').attr('x', -6).attr('y', i * cellSize + cellSize / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', theme.fontSizeSmall).attr('fill', theme.text).text(lbl)
    g.append('text').attr('x', i * cellSize + cellSize / 2).attr('y', -8)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall).attr('fill', theme.text).text(lbl)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Plot 6: Fit Indices Dashboard ──────────────────────────────────────

function renderFitIndices(
  d3: typeof D3,
  container: HTMLElement,
  data: FAResult,
  config: FAPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 500, H = config.height ?? 340
  const margin = { top: theme.marginTop, right: 40, bottom: 50, left: 80 }
  const width = W - margin.left - margin.right

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'Model Fit Indices', data.formatted, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const fit = data.fit

  // Gauge definitions: { label, value, range, goodThreshold, okThreshold, reversed }
  const gauges = [
    { label: 'RMSEA', value: fit.rmsea, range: [0, 0.15] as const, good: 0.05, ok: 0.08, reversed: true, ci: fit.rmseaCI },
    { label: 'CFI', value: fit.cfi, range: [0.8, 1.0] as const, good: 0.95, ok: 0.90, reversed: false },
    { label: 'TLI', value: fit.tli, range: [0.8, 1.0] as const, good: 0.95, ok: 0.90, reversed: false },
    { label: 'SRMR', value: fit.srmr, range: [0, 0.12] as const, good: 0.05, ok: 0.08, reversed: true },
  ]

  const barH = 24
  const gap = 50
  const startY = 10

  gauges.forEach((gauge, i) => {
    const y = startY + i * gap

    const xScale = d3.scaleLinear()
      .domain([gauge.range[0], gauge.range[1]])
      .range([0, width]).clamp(true)

    // Background zones
    if (gauge.reversed) {
      // Good (left) → OK (mid) → Poor (right)
      g.append('rect').attr('x', 0).attr('y', y).attr('width', xScale(gauge.good)).attr('height', barH)
        .attr('fill', '#59a14f').attr('opacity', 0.15).attr('rx', 4)
      g.append('rect').attr('x', xScale(gauge.good)).attr('y', y)
        .attr('width', xScale(gauge.ok) - xScale(gauge.good)).attr('height', barH)
        .attr('fill', '#f28e2b').attr('opacity', 0.15)
      g.append('rect').attr('x', xScale(gauge.ok)).attr('y', y)
        .attr('width', width - xScale(gauge.ok)).attr('height', barH)
        .attr('fill', '#e15759').attr('opacity', 0.15).attr('rx', 4)
    } else {
      // Poor (left) → OK (mid) → Good (right)
      g.append('rect').attr('x', 0).attr('y', y).attr('width', xScale(gauge.ok)).attr('height', barH)
        .attr('fill', '#e15759').attr('opacity', 0.15).attr('rx', 4)
      g.append('rect').attr('x', xScale(gauge.ok)).attr('y', y)
        .attr('width', xScale(gauge.good) - xScale(gauge.ok)).attr('height', barH)
        .attr('fill', '#f28e2b').attr('opacity', 0.15)
      g.append('rect').attr('x', xScale(gauge.good)).attr('y', y)
        .attr('width', width - xScale(gauge.good)).attr('height', barH)
        .attr('fill', '#59a14f').attr('opacity', 0.15).attr('rx', 4)
    }

    // Threshold lines
    g.append('line').attr('x1', xScale(gauge.good)).attr('x2', xScale(gauge.good))
      .attr('y1', y).attr('y2', y + barH)
      .attr('stroke', theme.axisLine).attr('stroke-dasharray', '2,2').attr('stroke-width', 1)
    g.append('line').attr('x1', xScale(gauge.ok)).attr('x2', xScale(gauge.ok))
      .attr('y1', y).attr('y2', y + barH)
      .attr('stroke', theme.axisLine).attr('stroke-dasharray', '2,2').attr('stroke-width', 1)

    // CI range for RMSEA
    if ('ci' in gauge && gauge.ci) {
      const ciLo = gauge.ci[0]
      const ciHi = gauge.ci[1]
      g.append('rect')
        .attr('x', xScale(ciLo)).attr('y', y + barH / 2 - 3)
        .attr('width', xScale(ciHi) - xScale(ciLo)).attr('height', 6)
        .attr('fill', getColor(0, theme)).attr('opacity', 0.3).attr('rx', 2)
    }

    // Value marker
    const valColor = gauge.reversed
      ? (gauge.value <= gauge.good ? '#59a14f' : gauge.value <= gauge.ok ? '#f28e2b' : '#e15759')
      : (gauge.value >= gauge.good ? '#59a14f' : gauge.value >= gauge.ok ? '#f28e2b' : '#e15759')

    g.append('circle')
      .attr('cx', xScale(gauge.value)).attr('cy', y + barH / 2)
      .attr('r', 6).attr('fill', valColor).attr('stroke', '#fff').attr('stroke-width', 2)

    // Label
    g.append('text').attr('x', -6).attr('y', y + barH / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', theme.fontSizeSmall)
      .attr('font-weight', '600').attr('fill', theme.text).text(gauge.label)

    // Value
    g.append('text').attr('x', width + 6).attr('y', y + barH / 2 + 4)
      .attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall)
      .attr('fill', valColor).attr('font-weight', '600').text(gauge.value.toFixed(3))
  })

  // Chi-square + AIC/BIC text below
  const textY = startY + gauges.length * gap + 10
  g.append('text').attr('x', 0).attr('y', textY)
    .attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textAnnotation)
    .text(`χ²(${fit.df}) = ${formatStat(fit.chiSq)}, ${formatP(fit.pValue)}`)
  g.append('text').attr('x', 0).attr('y', textY + 16)
    .attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textAnnotation)
    .text(`AIC = ${formatStat(fit.aic, 1)}, BIC = ${formatStat(fit.bic, 1)}`)

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
