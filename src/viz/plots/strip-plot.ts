/**
 * Strip plot: pure jittered dot plot per group.
 * Shows individual data points with a horizontal mean line per group.
 * No violin, no box — raw distribution only.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface StripPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface StripGroup {
  readonly label: string
  readonly values: readonly number[]
}

export interface StripPlotData {
  readonly groups: readonly StripGroup[]
  readonly testResult?: StatResult
}

export function renderStripPlot(
  container: HTMLElement,
  data: StripPlotData,
  config: StripPlotConfig = {}
): void {
  import('d3').then(d3 => renderStripPlotD3(d3, container, data, config))
}

function renderStripPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: StripPlotData,
  config: StripPlotConfig
): void {
  if (data.groups.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Strip Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const labels = data.groups.map(gr => gr.label)
  const allValues = data.groups.flatMap(gr => [...gr.values])
  const yMinRaw = Math.min(...allValues)
  const yMaxRaw = Math.max(...allValues)
  const yPad = (yMaxRaw - yMinRaw) * 0.1

  const xScale = d3.scaleBand<string>().domain(labels).range([0, width]).padding(0.3)
  const yScale = d3.scaleLinear().domain([yMinRaw - yPad, yMaxRaw + yPad]).range([height, 0]).nice()

  // Grid
  g.selectAll<SVGLineElement, number>('.grid').data(yScale.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  data.groups.forEach((gr, gi) => {
    if (gr.values.length === 0) return

    const color = getColor(gi, theme)
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2
    const jitterWidth = xScale.bandwidth() * 0.38
    const seed = gi * 77773

    // Mean
    const mean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length

    // Jittered points
    gr.values.forEach((v, vi) => {
      const jitter = (pseudoRnd(seed + vi * 13 + gi) - 0.5) * jitterWidth
      g.append('circle')
        .attr('cx', cx + jitter)
        .attr('cy', yScale(v))
        .attr('r', 3.5)
        .attr('fill', color)
        .attr('opacity', theme.pointOpacity)
        .on('mouseover', function(event: MouseEvent) {
          showTooltip(event, [
            formatTooltipRow('Group', gr.label),
            formatTooltipRow('Value', v.toFixed(4)),
            formatTooltipRow('Group mean', mean.toFixed(4))
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })

    // Mean line
    const lineHalfW = xScale.bandwidth() * 0.32
    g.append('line')
      .attr('x1', cx - lineHalfW).attr('x2', cx + lineHalfW)
      .attr('y1', yScale(mean)).attr('y2', yScale(mean))
      .attr('stroke', color)
      .attr('stroke-width', 2.5)
      .attr('opacity', 0.9)

    // n label below axis
    g.append('text')
      .attr('x', cx).attr('y', height + 52)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted)
      .text(`n = ${gr.values.length}`)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text')
    .attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('g').call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text')
    .attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')

  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Value')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

function pseudoRnd(seed: number): number {
  const x = Math.sin(seed) * 10000
  return x - Math.floor(x)
}
