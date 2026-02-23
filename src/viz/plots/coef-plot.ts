/**
 * Coefficient dot-and-whisker plot (ggcoefstats style).
 * Shows regression/LMM coefficients with CI bars and p-value annotations.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import type { RegressionCoef } from '../../core/types.js'

export interface CoefPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showZeroLine?: boolean
  readonly excludeIntercept?: boolean
}

export function renderCoefPlot(
  container: HTMLElement,
  coefficients: readonly RegressionCoef[],
  config: CoefPlotConfig = {}
): void {
  import('d3').then(d3 => renderCoefD3(d3, container, coefficients, config))
}

function renderCoefD3(
  d3: typeof D3,
  container: HTMLElement,
  coefficientsAll: readonly RegressionCoef[],
  config: CoefPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const coefs = config.excludeIntercept !== false
    ? coefficientsAll.filter(c => c.name !== '(Intercept)')
    : coefficientsAll

  const k = coefs.length
  const W = config.width ?? 500
  const H = config.height ?? Math.max(k * 40 + 100, 200)
  const margin = { top: theme.marginTop, right: theme.marginRight + 60, bottom: theme.marginBottom, left: 120 }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Coefficient Plot', '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allX = coefs.flatMap(c => [c.ci[0], c.ci[1], c.estimate])
  const [xMin, xMax] = d3.extent(allX) as [number, number]
  const xPad = (xMax - xMin) * 0.1
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()

  const yScale = d3.scaleBand<string>()
    .domain([...coefs].reverse().map(c => c.name))
    .range([0, height]).padding(0.3)

  // Zero line
  if (config.showZeroLine !== false) {
    const x0 = xScale(0)
    g.append('line').attr('x1', x0).attr('x2', x0)
      .attr('y1', 0).attr('y2', height)
      .attr('stroke', theme.axisLine).attr('stroke-dasharray', '4,2').attr('stroke-width', 1)
  }

  // Gridlines
  g.selectAll('.grid').data(xScale.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', d => xScale(d)).attr('x2', d => xScale(d))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // CI bars + dots
  coefs.forEach((coef, _i) => {
    const cy = (yScale(coef.name) ?? 0) + yScale.bandwidth() / 2
    const significant = coef.pValue < 0.05
    const color = significant ? getColor(0, theme) : theme.textMuted

    // CI bar
    g.append('line')
      .attr('x1', xScale(coef.ci[0])).attr('x2', xScale(coef.ci[1]))
      .attr('y1', cy).attr('y2', cy)
      .attr('stroke', color).attr('stroke-width', 2.5).attr('opacity', 0.6)

    // Dot
    g.append('circle')
      .attr('cx', xScale(coef.estimate)).attr('cy', cy)
      .attr('r', 5).attr('fill', color)

    // P-value label on right
    const pLabel = coef.pValue < 0.001 ? 'p < .001' : `p = ${coef.pValue.toFixed(3)}`
    g.append('text')
      .attr('x', width + 5).attr('y', cy + 4)
      .attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.textMuted).text(pLabel)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'Estimate (95% CI)')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
