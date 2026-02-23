/**
 * Lollipop chart: stem + circle for each category.
 * Horizontal layout (categories on y-axis, values on x-axis).
 * Optional descending sort.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface LollipopConfig {
  readonly sorted?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface LollipopData {
  readonly labels: readonly string[]
  readonly values: readonly number[]
  readonly testResult?: StatResult
}

export function renderLollipop(
  container: HTMLElement,
  data: LollipopData,
  config: LollipopConfig = {}
): void {
  import('d3').then(d3 => renderLollipopD3(d3, container, data, config))
}

function renderLollipopD3(
  d3: typeof D3,
  container: HTMLElement,
  data: LollipopData,
  config: LollipopConfig
): void {
  if (data.labels.length === 0 || data.values.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const itemCount = data.labels.length
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? Math.max(200, itemCount * 36 + theme.marginTop + theme.marginBottom)
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft + 40 }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Lollipop Chart', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // Build ordered index
  let indices = data.labels.map((_, i) => i)
  if (config.sorted) {
    indices = [...indices].sort((a, b) => (data.values[b] ?? 0) - (data.values[a] ?? 0))
  }

  const orderedLabels = indices.map(i => data.labels[i] ?? '')
  const orderedValues = indices.map(i => data.values[i] ?? 0)

  const [vMin, vMax] = d3.extent(orderedValues) as [number, number]
  const xMin = Math.min(0, vMin)
  const xMax = Math.max(0, vMax)
  const xPad = (xMax - xMin) * 0.1

  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()
  const yScale = d3.scaleBand<string>().domain(orderedLabels).range([0, height]).padding(0.35)

  // Grid lines (vertical)
  g.selectAll('.grid').data(xScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', d => xScale(d)).attr('x2', d => xScale(d))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Zero line
  g.append('line')
    .attr('x1', xScale(0)).attr('x2', xScale(0))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.axisLine).attr('stroke-width', 1.5)

  // Lollipops
  orderedLabels.forEach((label, i) => {
    const val = orderedValues[i] ?? 0
    const cy = (yScale(label) ?? 0) + yScale.bandwidth() / 2
    const color = getColor(0, theme)

    // Stem
    g.append('line')
      .attr('x1', xScale(0)).attr('x2', xScale(val))
      .attr('y1', cy).attr('y2', cy)
      .attr('stroke', color).attr('stroke-width', 2).attr('opacity', 0.7)

    // Circle
    g.append('circle')
      .attr('cx', xScale(val)).attr('cy', cy).attr('r', 6)
      .attr('fill', color).attr('stroke', theme.background).attr('stroke-width', 1.5)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('Label', label),
          formatTooltipRow('Value', val.toFixed(3)),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'Value')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -margin.left + 12)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? '')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
