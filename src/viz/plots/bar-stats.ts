/**
 * Annotated bar chart with counts, percentages, and significance brackets.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import type { StatResult, FrequencyRow } from '../../core/types.js'

export interface BarStatsConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showPercentages?: boolean
  readonly showCounts?: boolean
  readonly orientation?: 'vertical' | 'horizontal'
}

export interface BarStatsData {
  readonly rows: readonly FrequencyRow[]
  readonly testResult?: StatResult
}

export function renderBarStats(
  container: HTMLElement,
  data: BarStatsData,
  config: BarStatsConfig = {}
): void {
  import('d3').then(d3 => renderBarD3(d3, container, data, config))
}

function renderBarD3(
  d3: typeof D3,
  container: HTMLElement,
  data: BarStatsData,
  config: BarStatsConfig
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

  addSubtitle(svg, config.title ?? 'Frequency', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const categories = data.rows.map(r => String(r.value))
  const maxCount = Math.max(...data.rows.map(r => r.count))

  const xScale = d3.scaleBand<string>().domain(categories).range([0, width]).padding(0.25)
  const yScale = d3.scaleLinear().domain([0, maxCount * 1.15]).range([height, 0]).nice()

  // Grid
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Bars
  g.selectAll('.bar').data(data.rows).join('rect')
    .attr('class', 'bar')
    .attr('x', row => xScale(String(row.value)) ?? 0)
    .attr('y', row => yScale(row.count))
    .attr('width', xScale.bandwidth())
    .attr('height', row => height - yScale(row.count))
    .attr('fill', (_, i) => getColor(i, theme))
    .attr('rx', 2)

  // Count labels
  if (config.showCounts !== false) {
    g.selectAll('.label-count').data(data.rows).join('text')
      .attr('class', 'label-count')
      .attr('x', row => (xScale(String(row.value)) ?? 0) + xScale.bandwidth() / 2)
      .attr('y', row => yScale(row.count) - 4)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(row => `${row.count}`)
  }

  // Percentage labels
  if (config.showPercentages !== false) {
    g.selectAll('.label-pct').data(data.rows).join('text')
      .attr('class', 'label-pct')
      .attr('x', row => (xScale(String(row.value)) ?? 0) + xScale.bandwidth() / 2)
      .attr('y', row => yScale(row.count) - 16)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted)
      .text(row => `(${(row.relative * 100).toFixed(1)}%)`)
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'Category')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Count')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
