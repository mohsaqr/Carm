/**
 * Pareto chart: descending bars + cumulative % line on dual axis.
 * Includes 80% threshold line (the "vital few" cutoff).
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ParetoConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface ParetoData {
  readonly labels: readonly string[]
  readonly values: readonly number[]
  readonly testResult?: StatResult
}

export function renderPareto(
  container: HTMLElement,
  data: ParetoData,
  config: ParetoConfig = {}
): void {
  import('d3').then(d3 => renderParetoD3(d3, container, data, config))
}

interface ParetoPoint {
  label: string
  value: number
  cumPct: number
  index: number
}

function renderParetoD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ParetoData,
  config: ParetoConfig
): void {
  if (data.labels.length === 0 || data.values.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: 64, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Pareto Chart', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // Sort by descending value
  const indices = Array.from({ length: data.labels.length }, (_, i) => i)
    .sort((a, b) => (data.values[b] ?? 0) - (data.values[a] ?? 0))

  const total = indices.reduce((s, i) => s + (data.values[i] ?? 0), 0)
  let running = 0
  const points: ParetoPoint[] = indices.map((origIdx, sortPos) => {
    const val = data.values[origIdx] ?? 0
    const lbl = data.labels[origIdx] ?? ''
    running += val
    return {
      label: lbl,
      value: val,
      cumPct: (running / total) * 100,
      index: sortPos
    }
  })

  const sortedLabels = points.map(p => p.label)

  // Scales
  const xScale = d3.scaleBand<string>().domain(sortedLabels).range([0, width]).padding(0.25)
  const yScaleBar = d3.scaleLinear().domain([0, (points[0]?.value ?? 0) * 1.1]).range([height, 0]).nice()
  const yScalePct = d3.scaleLinear().domain([0, 100]).range([height, 0])

  // Grid lines
  g.selectAll<SVGLineElement, number>('.grid').data(yScaleBar.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScaleBar(d)).attr('y2', d => yScaleBar(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  const barColor = getColor(0, theme)
  const lineColor = getColor(1, theme)

  // Bars
  g.selectAll<SVGRectElement, ParetoPoint>('.bar').data(points).join('rect')
    .attr('class', 'bar')
    .attr('x', d => xScale(d.label) ?? 0)
    .attr('y', d => yScaleBar(d.value))
    .attr('width', xScale.bandwidth())
    .attr('height', d => height - yScaleBar(d.value))
    .attr('fill', barColor)
    .attr('opacity', 0.85)
    .on('mouseover', function(event: MouseEvent, d: ParetoPoint) {
      showTooltip(event, [
        formatTooltipRow('Category', d.label),
        formatTooltipRow('Value', d.value.toFixed(2)),
        formatTooltipRow('Cumulative %', d.cumPct.toFixed(1) + '%')
      ].join(''), theme)
    })
    .on('mouseout', hideTooltip)

  // Cumulative line
  const lineGen = d3.line<ParetoPoint>()
    .x(d => (xScale(d.label) ?? 0) + xScale.bandwidth() / 2)
    .y(d => yScalePct(d.cumPct))
    .curve(d3.curveMonotoneX)

  g.append('path').datum(points)
    .attr('d', lineGen)
    .attr('fill', 'none')
    .attr('stroke', lineColor)
    .attr('stroke-width', 2.5)

  // Dots on cumulative line
  g.selectAll<SVGCircleElement, ParetoPoint>('.cum-dot').data(points).join('circle')
    .attr('class', 'cum-dot')
    .attr('cx', d => (xScale(d.label) ?? 0) + xScale.bandwidth() / 2)
    .attr('cy', d => yScalePct(d.cumPct))
    .attr('r', 4)
    .attr('fill', lineColor)
    .attr('stroke', theme.background)
    .attr('stroke-width', 1.5)

  // 80% threshold line
  g.append('line')
    .attr('x1', 0).attr('x2', width)
    .attr('y1', yScalePct(80)).attr('y2', yScalePct(80))
    .attr('stroke', '#D55E00')
    .attr('stroke-width', 1.5)
    .attr('stroke-dasharray', '6,4')

  g.append('text')
    .attr('x', width + 4)
    .attr('y', yScalePct(80) + 4)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall)
    .attr('fill', '#D55E00')
    .text('80%')

  // Left axis (bar values)
  const leftAxis = g.append('g').call(d3.axisLeft(yScaleBar).ticks(6))
  leftAxis.selectAll('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  // Right axis (%)
  const rightAxis = g.append('g')
    .attr('transform', `translate(${width},0)`)
    .call(d3.axisRight(yScalePct).ticks(5).tickFormat(d => `${d}%`))
  rightAxis.selectAll('text')
    .attr('fill', lineColor)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  // Bottom axis
  const bottomAxis = g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
  bottomAxis.selectAll('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('text-anchor', 'end')
    .attr('transform', 'rotate(-35)')

  // Axis labels
  g.append('text').attr('x', width / 2).attr('y', height + 60)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')

  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Frequency')

  g.append('text')
    .attr('transform', `translate(${width + 56},${height / 2}) rotate(90)`)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', lineColor).text('Cumulative %')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
