/**
 * QQ plot (quantile-quantile plot) with confidence band.
 * Tests normality assumption visually.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { normalQuantile, sortAsc } from '../../core/math.js'

export interface QQPlotConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showCI?: boolean
}

export function renderQQPlot(
  container: HTMLElement,
  values: readonly number[],
  config: QQPlotConfig = {}
): void {
  import('d3').then(d3 => renderQQD3(d3, container, values, config))
}

function renderQQD3(
  d3: typeof D3,
  container: HTMLElement,
  values: readonly number[],
  config: QQPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 400
  const H = config.height ?? 400
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  const n = values.length
  const sorted = sortAsc(values)
  const mean_ = sorted.reduce((s, v) => s + v, 0) / n
  const sd_ = Math.sqrt(sorted.reduce((s, v) => s + (v - mean_) ** 2, 0) / (n - 1))

  // Theoretical quantiles (Blom positions)
  const points = sorted.map((y, i) => {
    const p = (i + 1 - 0.375) / (n + 0.25)
    const x = normalQuantile(Math.max(0.0001, Math.min(0.9999, p)))
    return { x, y }
  })

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'QQ Plot', 'Normal probability plot', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const xExtent = d3.extent(points.map(p => p.x)) as [number, number]
  const yExtent = d3.extent(points.map(p => p.y)) as [number, number]
  const xPad = (xExtent[1] - xExtent[0]) * 0.05
  const yPad = (yExtent[1] - yExtent[0]) * 0.05

  const xScale = d3.scaleLinear().domain([xExtent[0] - xPad, xExtent[1] + xPad]).range([0, width]).nice()
  const yScale = d3.scaleLinear().domain([yExtent[0] - yPad, yExtent[1] + yPad]).range([height, 0]).nice()

  // Grid
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Theoretical reference line: y = mean + sd * x
  const xDom = xScale.domain()
  const xd0 = xDom[0] ?? 0, xd1 = xDom[1] ?? 1
  g.append('line')
    .attr('x1', xScale(xd0)).attr('x2', xScale(xd1))
    .attr('y1', yScale(mean_ + sd_ * xd0)).attr('y2', yScale(mean_ + sd_ * xd1))
    .attr('stroke', getColor(4, theme)).attr('stroke-width', 2).attr('stroke-dasharray', '6,3')

  // Confidence band (approximate 95% simultaneous CI using Kolmogorov-Smirnov bounds)
  if (config.showCI !== false) {
    const CI_MULT = 1.36 / Math.sqrt(n)  // 95% KS bound
    const bandData = points.map(p => ({
      x: p.x,
      lo: (mean_ + sd_ * p.x) - CI_MULT * sd_ * Math.sqrt(n),
      hi: (mean_ + sd_ * p.x) + CI_MULT * sd_ * Math.sqrt(n),
    }))

    const area = d3.area<typeof bandData[0]>()
      .x(d => xScale(d.x)).y0(d => yScale(d.lo)).y1(d => yScale(d.hi))
    g.append('path').datum(bandData).attr('d', area)
      .attr('fill', getColor(4, theme)).attr('opacity', theme.ciOpacity)
  }

  // Data points
  g.selectAll('.qq-point').data(points).join('circle')
    .attr('class', 'qq-point')
    .attr('cx', p => xScale(p.x)).attr('cy', p => yScale(p.y))
    .attr('r', 3).attr('fill', getColor(0, theme))
    .attr('opacity', theme.pointOpacity)

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text('Theoretical Quantiles')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text('Sample Quantiles')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
