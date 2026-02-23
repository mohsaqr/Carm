/**
 * Density plot: KDE curves for one or more series with optional rug ticks.
 * Uses Silverman's rule-of-thumb bandwidth: 1.06 * σ * n^(-1/5)
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface DensityConfig {
  readonly bandwidth?: number
  readonly showRug?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface DensitySeries {
  readonly label: string
  readonly values: readonly number[]
}

export interface DensityData {
  readonly series: readonly DensitySeries[]
  readonly testResult?: StatResult
}

export function renderDensity(
  container: HTMLElement,
  data: DensityData,
  config: DensityConfig = {}
): void {
  import('d3').then(d3 => renderDensityD3(d3, container, data, config))
}

function silvermanBW(values: readonly number[]): number {
  const n = values.length
  if (n < 2) return 1
  const m = values.reduce((s, v) => s + v, 0) / n
  const std = Math.sqrt(values.reduce((s, v) => s + (v - m) ** 2, 0) / (n - 1))
  return 1.06 * std * Math.pow(n, -0.2)
}

function kdeGaussian(values: readonly number[], xPoints: number[], bw: number): Array<[number, number]> {
  return xPoints.map(x => {
    const density = values.reduce((s, xi) =>
      s + Math.exp(-0.5 * ((x - xi) / bw) ** 2) / (Math.sqrt(2 * Math.PI) * bw), 0
    ) / values.length
    return [x, density] as [number, number]
  })
}

function renderDensityD3(
  d3: typeof D3,
  container: HTMLElement,
  data: DensityData,
  config: DensityConfig
): void {
  if (data.series.length === 0 || data.series.every(s => s.values.length === 0)) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 420
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom
  const showRug = config.showRug ?? true

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Density Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allValues = data.series.flatMap(s => [...s.values])
  const [xMin, xMax] = d3.extent(allValues) as [number, number]
  const xPad = (xMax - xMin) * 0.1
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()

  // Compute KDE for each series
  const xTicks = xScale.ticks(100)
  const seriesDensities = data.series.map(s => {
    const bw = config.bandwidth ?? silvermanBW(s.values)
    return kdeGaussian(s.values, xTicks, bw)
  })

  const maxDensity = Math.max(...seriesDensities.flatMap(pts => pts.map(p => p[1])))
  const rugHeight = showRug ? height - 12 : height
  const yScale = d3.scaleLinear().domain([0, maxDensity * 1.1]).range([rugHeight, 0]).nice()

  // Grid lines
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Draw density curves and areas
  seriesDensities.forEach((pts, i) => {
    const color = getColor(i, theme)
    const lineGen = d3.line<[number, number]>()
      .x(d => xScale(d[0])).y(d => yScale(d[1])).curve(d3.curveCatmullRom)
    const areaGen = d3.area<[number, number]>()
      .x(d => xScale(d[0])).y0(yScale(0)).y1(d => yScale(d[1])).curve(d3.curveCatmullRom)

    g.append('path').datum(pts).attr('d', areaGen)
      .attr('fill', color).attr('opacity', theme.ciOpacity)
    g.append('path').datum(pts).attr('d', lineGen)
      .attr('fill', 'none').attr('stroke', color).attr('stroke-width', 2)
  })

  // Rug ticks
  if (showRug) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme)
      const rugY = height - 10
      s.values.forEach(v => {
        g.append('line')
          .attr('x1', xScale(v)).attr('x2', xScale(v))
          .attr('y1', rugY).attr('y2', rugY + 6)
          .attr('stroke', color).attr('stroke-width', 1).attr('opacity', 0.5)
          .on('mouseover', (event: MouseEvent) => {
            showTooltip(event,
              [formatTooltipRow('Series', s.label), formatTooltipRow('Value', v.toFixed(3))].join(''),
              theme)
          })
          .on('mouseout', hideTooltip)
      })
    })
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${rugHeight})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Density')

  // Legend (if multiple series)
  if (data.series.length > 1) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme)
      const lx = width - 120
      const ly = i * 18 + 8
      g.append('line').attr('x1', lx).attr('x2', lx + 20).attr('y1', ly).attr('y2', ly)
        .attr('stroke', color).attr('stroke-width', 2)
      g.append('text').attr('x', lx + 24).attr('y', ly + 4)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text).text(s.label)
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
