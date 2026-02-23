/**
 * 4-panel regression diagnostics plot:
 * 1. Residuals vs Fitted
 * 2. QQ of residuals
 * 3. Scale-Location (sqrt|residuals| vs fitted)
 * 4. Residuals vs Leverage
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import type { RegressionResult } from '../../core/types.js'
import { normalQuantile, sortAsc, mean as _mean } from '../../core/math.js'

export interface ResidualPanelConfig {
  readonly title?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export function renderResidualPanel(
  container: HTMLElement,
  result: RegressionResult,
  leverage: readonly number[],
  config: ResidualPanelConfig = {}
): void {
  import('d3').then(d3 => renderResidualD3(d3, container, result, leverage, config))
}

function renderResidualD3(
  d3: typeof D3,
  container: HTMLElement,
  result: RegressionResult,
  leverage: readonly number[],
  config: ResidualPanelConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 700
  const H = config.height ?? 600
  const panelW = W / 2, panelH = H / 2

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  // Title
  svg.append('text').attr('x', W / 2).attr('y', 18)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeTitle).attr('font-weight', '600').attr('fill', theme.text)
    .text(config.title ?? 'Regression Diagnostics')

  const color = getColor(0, theme)
  const panels = [
    { title: 'Residuals vs Fitted', row: 0, col: 0 },
    { title: 'Normal Q-Q', row: 0, col: 1 },
    { title: 'Scale-Location', row: 1, col: 0 },
    { title: 'Residuals vs Leverage', row: 1, col: 1 },
  ]

  const margin = { top: 40, right: 20, bottom: 50, left: 50 }
  const pw = panelW - margin.left - margin.right
  const ph = panelH - margin.top - margin.bottom

  panels.forEach((panel, idx) => {
    const ox = panel.col * panelW + margin.left
    const oy = panel.row * panelH + margin.top + 25
    const g = svg.append('g').attr('transform', `translate(${ox},${oy})`)

    g.append('text').attr('x', pw / 2).attr('y', -18)
      .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize).attr('font-weight', '600').attr('fill', theme.text)
      .text(panel.title)

    let xData: readonly number[]
    let yData: readonly number[]

    if (idx === 0) {
      xData = result.fitted; yData = result.residuals
    } else if (idx === 1) {
      // QQ
      const sorted = sortAsc(result.residuals)
      const n = sorted.length
      yData = sorted
      xData = sorted.map((_, i) => {
        const p = (i + 1 - 0.375) / (n + 0.25)
        return normalQuantile(Math.max(0.0001, Math.min(0.9999, p)))
      })
    } else if (idx === 2) {
      xData = result.fitted
      yData = result.residuals.map(r => Math.sqrt(Math.abs(r)))
    } else {
      xData = leverage; yData = result.residuals
    }

    const xExt = d3.extent([...xData]) as [number, number]
    const yExt = d3.extent([...yData]) as [number, number]
    const xPad = (xExt[1] - xExt[0]) * 0.05
    const yPad = (yExt[1] - yExt[0]) * 0.05

    const xs = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, pw]).nice()
    const ys = d3.scaleLinear().domain([yExt[0] - yPad, yExt[1] + yPad]).range([ph, 0]).nice()

    // Grid
    g.selectAll('.grid').data(ys.ticks(4)).join('line')
      .attr('class', 'grid').attr('x1', 0).attr('x2', pw)
      .attr('y1', d => ys(d)).attr('y2', d => ys(d))
      .attr('stroke', theme.gridLine).attr('stroke-width', 1)

    // Zero line for residual panels
    if (idx === 0 || idx === 2) {
      g.append('line').attr('x1', 0).attr('x2', pw)
        .attr('y1', ys(0)).attr('y2', ys(0))
        .attr('stroke', theme.axisLine).attr('stroke-dasharray', '4,2')
    }

    // Reference line for QQ
    if (idx === 1) {
      const mean_ = _mean(result.residuals)
      const sd_ = Math.sqrt(result.residuals.reduce((s, r) => s + (r - mean_) ** 2, 0) / (result.residuals.length - 1))
      const xd = xs.domain()
      const xd0 = xd[0] ?? 0, xd1 = xd[1] ?? 1
      g.append('line')
        .attr('x1', xs(xd0)).attr('x2', xs(xd1))
        .attr('y1', ys(mean_ + sd_ * xd0)).attr('y2', ys(mean_ + sd_ * xd1))
        .attr('stroke', getColor(4, theme)).attr('stroke-width', 1.5).attr('stroke-dasharray', '5,2')
    }

    // Points
    xData.forEach((xi, i) => {
      g.append('circle')
        .attr('cx', xs(xi)).attr('cy', ys(yData[i] ?? 0))
        .attr('r', 2.5).attr('fill', color).attr('opacity', theme.pointOpacity)
    })

    // Axes
    g.append('g').attr('transform', `translate(0,${ph})`).call(d3.axisBottom(xs).ticks(4))
      .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
    g.append('g').call(d3.axisLeft(ys).ticks(4))
      .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
  })
}
