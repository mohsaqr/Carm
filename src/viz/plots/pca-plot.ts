/**
 * PCA visualizations: biplot, scree plot, loadings heatmap.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import type { PCAResult } from '../../core/types.js'

export interface PCAPlotConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly variableLabels?: readonly string[]
  readonly observationLabels?: readonly string[]
  readonly type?: 'biplot' | 'scree' | 'loadings'
}

export function renderPCAPlot(
  container: HTMLElement,
  pca: PCAResult,
  config: PCAPlotConfig = {}
): void {
  import('d3').then(d3 => {
    const type = config.type ?? 'biplot'
    if (type === 'biplot') renderBiplot(d3, container, pca, config)
    else if (type === 'scree') renderScree(d3, container, pca, config)
    else renderLoadingsHeatmap(d3, container, pca, config)
  })
}

function renderBiplot(d3: typeof D3, container: HTMLElement, pca: PCAResult, config: PCAPlotConfig): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 500, H = config.height ?? 500
  const margin = { top: theme.marginTop, right: 60, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const pct1 = ((pca.varianceExplained[0] ?? 0) * 100).toFixed(1)
  const pct2 = ((pca.varianceExplained[1] ?? 0) * 100).toFixed(1)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'PCA Biplot', `PC1 (${pct1}%) × PC2 (${pct2}%)`, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const scores = pca.scores
  const xs = scores.map(s => s[0] ?? 0), ys = scores.map(s => s[1] ?? 0)
  const xExt = d3.extent(xs) as [number, number], yExt = d3.extent(ys) as [number, number]
  const xPad = (xExt[1] - xExt[0]) * 0.1, yPad = (yExt[1] - yExt[0]) * 0.1

  const xScale = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, width]).nice()
  const yScale = d3.scaleLinear().domain([yExt[0] - yPad, yExt[1] + yPad]).range([height, 0]).nice()

  // Grid
  g.selectAll('.grid-h').data(yScale.ticks(5)).join('line')
    .attr('x1', 0).attr('x2', width).attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Zero lines
  g.append('line').attr('x1', xScale(0)).attr('x2', xScale(0)).attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '4,2')
  g.append('line').attr('x1', 0).attr('x2', width).attr('y1', yScale(0)).attr('y2', yScale(0))
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '4,2')

  // Observation points
  scores.forEach((s, i) => {
    g.append('circle').attr('cx', xScale(s[0] ?? 0)).attr('cy', yScale(s[1] ?? 0))
      .attr('r', 3.5).attr('fill', getColor(0, theme)).attr('opacity', theme.pointOpacity)
    if (config.observationLabels?.[i]) {
      g.append('text').attr('x', xScale(s[0] ?? 0) + 5).attr('y', yScale(s[1] ?? 0) + 4)
        .attr('font-size', theme.fontSizeSmall - 2).attr('fill', theme.textMuted)
        .text(config.observationLabels[i]!)
    }
  })

  // Loading vectors
  const scale_ = Math.min(width, height) * 0.4
  pca.loadings.forEach((loading, vi) => {
    const lx = (loading[0] ?? 0) * scale_, ly = (loading[1] ?? 0) * scale_
    const x0 = xScale(0), y0 = yScale(0)
    const color = getColor(vi % 8 + 1, theme)
    g.append('line')
      .attr('x1', x0).attr('y1', y0)
      .attr('x2', x0 + lx).attr('y2', y0 - ly)
      .attr('stroke', color).attr('stroke-width', 1.5)
      .attr('marker-end', 'url(#arrow)')
    g.append('text')
      .attr('x', x0 + lx + 5).attr('y', y0 - ly + 4)
      .attr('font-size', theme.fontSizeSmall).attr('fill', color)
      .text(config.variableLabels?.[vi] ?? `Var${vi + 1}`)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text)
    .text(`PC1 (${pct1}%)`)
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text)
    .text(`PC2 (${pct2}%)`)

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

function renderScree(d3: typeof D3, container: HTMLElement, pca: PCAResult, config: PCAPlotConfig): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 500, H = config.height ?? 350
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'Scree Plot', 'Eigenvalue by component', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const comps = pca.eigenvalues.map((_, i) => i + 1)
  const xScale = d3.scaleBand<number>().domain(comps).range([0, width]).padding(0.3)
  const yScale = d3.scaleLinear().domain([0, Math.max(...pca.eigenvalues) * 1.1]).range([height, 0]).nice()

  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('x1', 0).attr('x2', width).attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Eigenvalue=1 reference line
  g.append('line').attr('x1', 0).attr('x2', width)
    .attr('y1', yScale(1)).attr('y2', yScale(1))
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '6,3').attr('stroke-width', 1.5)

  // Bars
  pca.eigenvalues.forEach((ev, i) => {
    const x = xScale(i + 1) ?? 0
    g.append('rect').attr('x', x).attr('y', yScale(ev))
      .attr('width', xScale.bandwidth()).attr('height', height - yScale(ev))
      .attr('fill', getColor(i, theme)).attr('opacity', 0.8).attr('rx', 2)
    g.append('text').attr('x', x + xScale.bandwidth() / 2).attr('y', yScale(ev) - 4)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text).text(ev.toFixed(2))
  })

  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Component')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Eigenvalue')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

function renderLoadingsHeatmap(d3: typeof D3, container: HTMLElement, pca: PCAResult, config: PCAPlotConfig): void {
  const theme = config.theme ?? DEFAULT_THEME
  const nVars = pca.loadings.length
  const nc = pca.nComponents
  const cellSize = 50
  const W = config.width ?? (nc * cellSize + 120), H = config.height ?? (nVars * cellSize + 100)
  const margin = { top: 60, right: 20, bottom: 40, left: 100 }

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? 'PCA Loadings', '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)
  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu)

  pca.loadings.forEach((row, vi) => {
    row.forEach((val, ci) => {
      g.append('rect').attr('x', ci * cellSize).attr('y', vi * cellSize)
        .attr('width', cellSize - 2).attr('height', cellSize - 2).attr('rx', 3)
        .attr('fill', colorScale(val))
      g.append('text').attr('x', ci * cellSize + cellSize / 2).attr('y', vi * cellSize + cellSize / 2 + 4)
        .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', Math.abs(val) > 0.6 ? '#fff' : theme.text).text(val.toFixed(2))
    })
    const label = config.variableLabels?.[vi] ?? `Var${vi + 1}`
    g.append('text').attr('x', -6).attr('y', vi * cellSize + cellSize / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', theme.fontSizeSmall).attr('fill', theme.text).text(label)
  })

  Array.from({ length: nc }, (_, ci) => {
    g.append('text').attr('x', ci * cellSize + cellSize / 2).attr('y', -8)
      .attr('text-anchor', 'middle').attr('font-size', theme.fontSizeSmall).attr('fill', theme.text).text(`PC${ci + 1}`)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
