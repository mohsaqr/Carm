/**
 * Histogram with density overlay and optional normality curve.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import type { DescriptiveResult } from '../../core/types.js'

export interface HistogramConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly bins?: number
  readonly showDensity?: boolean
  readonly showNormalCurve?: boolean
  readonly color?: string
}

export interface HistogramData {
  readonly values: readonly number[]
  readonly descriptives?: DescriptiveResult
}

export function renderHistogram(
  container: HTMLElement,
  data: HistogramData,
  config: HistogramConfig = {}
): void {
  import('d3').then(d3 => renderHistogramD3(d3, container, data, config))
}

function renderHistogramD3(
  d3: typeof D3,
  container: HTMLElement,
  data: HistogramData,
  config: HistogramConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 400
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(
    svg,
    config.title ?? 'Distribution',
    data.descriptives ? `M = ${data.descriptives.mean.toFixed(2)}, SD = ${data.descriptives.sd.toFixed(2)}, n = ${data.descriptives.n}` : '',
    W, theme
  )

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const [xMin, xMax] = d3.extent(data.values) as [number, number]
  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]).nice()

  const nBins = config.bins ?? Math.ceil(Math.sqrt(data.values.length))
  const histogram = d3.bin().domain(xScale.domain() as [number, number]).thresholds(nBins)
  const bins = histogram(data.values as number[])

  const maxCount = Math.max(...bins.map(b => b.length))
  const yScale = d3.scaleLinear().domain([0, maxCount * 1.1]).range([height, 0])

  // Grid
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Bars
  const barColor = config.color ?? getColor(0, theme)
  g.selectAll('.bar').data(bins).join('rect')
    .attr('class', 'bar')
    .attr('x', b => xScale(b.x0!))
    .attr('y', b => yScale(b.length))
    .attr('width', b => Math.max(0, xScale(b.x1!) - xScale(b.x0!) - 1))
    .attr('height', b => height - yScale(b.length))
    .attr('fill', barColor).attr('opacity', 0.7)
    .attr('stroke', theme.background).attr('stroke-width', 0.5)

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'Value')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text('Count')

  // Normal curve overlay
  if (config.showNormalCurve !== false && data.descriptives) {
    const m = data.descriptives.mean
    const s = data.descriptives.sd
    const n = data.descriptives.n
    const binWidth = bins[0]?.x1! - bins[0]?.x0!
    const scale_ = n * binWidth  // scale density to count

    const xTicks = d3.range(xMin, xMax, (xMax - xMin) / 100)
    const normalPoints = xTicks.map(x => ({
      x,
      y: scale_ * Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI)),
    }))

    const line = d3.line<typeof normalPoints[0]>()
      .x(d => xScale(d.x)).y(d => yScale(d.y)).curve(d3.curveBasis)

    g.append('path').datum(normalPoints).attr('d', line)
      .attr('fill', 'none').attr('stroke', getColor(4, theme))
      .attr('stroke-width', 2).attr('stroke-dasharray', '6,3')
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
