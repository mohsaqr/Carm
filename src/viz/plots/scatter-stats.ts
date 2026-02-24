/**
 * Scatter plot with regression line, CI band, and marginal distributions.
 * ggscatterstats style: stat result in subtitle, regression equation on plot.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption, addRegressionEquation } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult, RegressionResult } from '../../core/types.js'

export interface ScatterStatsConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showMarginals?: boolean
  readonly showCI?: boolean
  readonly pointSize?: number
  readonly showEquation?: boolean
}

export interface ScatterStatsData {
  readonly x: readonly number[]
  readonly y: readonly number[]
  readonly labels?: readonly string[]
  readonly correlationResult?: StatResult
  readonly regressionResult?: RegressionResult
}

export function renderScatterStats(
  container: HTMLElement,
  data: ScatterStatsData,
  config: ScatterStatsConfig = {}
): void {
  import('d3').then(d3 => renderScatterD3(d3, container, data, config))
}

function renderScatterD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ScatterStatsData,
  config: ScatterStatsConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight + 40, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W).attr('height', H)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Scatter Plot', data.correlationResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const [xMin, xMax] = d3.extent(data.x) as [number, number]
  const [yMin, yMax] = d3.extent(data.y) as [number, number]
  const xPad = (xMax - xMin) * 0.05
  const yPad = (yMax - yMin) * 0.05

  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice()

  // Grid
  g.selectAll('.grid-h').data(yScale.ticks(6)).join('line')
    .attr('class', 'grid-h').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  // Labels
  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'X')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Y')

  // Regression line + CI band
  if (data.regressionResult) {
    const coefs = data.regressionResult.coefficients
    const b0 = coefs[0]?.estimate ?? 0
    const b1 = coefs[1]?.estimate ?? 1
    const r2 = data.regressionResult.r2

    const xRange = xScale.domain()
    const lineData: Array<[number, number]> = [
      [xRange[0]!, b0 + b1 * xRange[0]!],
      [xRange[1]!, b0 + b1 * xRange[1]!],
    ]

    if (config.showCI !== false) {
      // Approximate CI band using residual SE
      const n = data.x.length
      const se_reg = Math.sqrt(data.regressionResult.residuals.reduce((s, r) => s + r * r, 0) / (n - 2))
      const xMean = data.x.reduce((s, v) => s + v, 0) / n
      const sxx = data.x.reduce((s, v) => s + (v - xMean) ** 2, 0)

      const ciPoints = xScale.ticks(50).map(x => {
        const yHat = b0 + b1 * x
        const h = 1 / n + (x - xMean) ** 2 / sxx
        const halfCI = 1.96 * se_reg * Math.sqrt(h)
        return { x, lo: yHat - halfCI, hi: yHat + halfCI }
      })

      const area = d3.area<typeof ciPoints[0]>()
        .x(d => xScale(d.x)).y0(d => yScale(d.lo)).y1(d => yScale(d.hi))
        .curve(d3.curveBasis)

      g.append('path').datum(ciPoints).attr('d', area)
        .attr('fill', getColor(0, theme)).attr('opacity', theme.ciOpacity)
    }

    g.append('line')
      .attr('x1', xScale(lineData[0]![0]!)).attr('x2', xScale(lineData[1]![0]!))
      .attr('y1', yScale(lineData[0]![1]!)).attr('y2', yScale(lineData[1]![1]!))
      .attr('stroke', getColor(0, theme)).attr('stroke-width', 2)

    if (config.showEquation !== false) {
      addRegressionEquation(g as unknown as D3.Selection<SVGGElement, unknown, null, undefined>, b0, b1, r2, 10, 20, theme)
    }
  }

  // Data points
  const color = getColor(0, theme)
  data.x.forEach((xi, i) => {
    const yi = data.y[i] ?? 0
    g.append('circle')
      .attr('cx', xScale(xi)).attr('cy', yScale(yi))
      .attr('r', config.pointSize ?? 4)
      .attr('fill', color).attr('opacity', theme.pointOpacity)
      .attr('stroke', theme.background).attr('stroke-width', 0.5)
      .on('mouseover', (event: MouseEvent) => {
        const rows = [formatTooltipRow(config.xLabel ?? 'X', xi.toFixed(3)), formatTooltipRow(config.yLabel ?? 'Y', yi.toFixed(3))]
        if (data.labels?.[i]) rows.unshift(`<strong>${data.labels[i]}</strong>`)
        showTooltip(event, rows.join(''), theme)
      })
      .on('mouseout', hideTooltip)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
