/**
 * Grouped or stacked bar chart: multiple series per category.
 * type='grouped' renders bars side-by-side; type='stacked' stacks them.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface GroupedBarConfig {
  readonly type?: 'grouped' | 'stacked'
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface GroupedBarSeries {
  readonly label: string
  readonly values: readonly number[]
}

export interface GroupedBarData {
  readonly categories: readonly string[]
  readonly series: readonly GroupedBarSeries[]
  readonly testResult?: StatResult
}

export function renderGroupedBar(
  container: HTMLElement,
  data: GroupedBarData,
  config: GroupedBarConfig = {}
): void {
  import('d3').then(d3 => renderGroupedBarD3(d3, container, data, config))
}

function renderGroupedBarD3(
  d3: typeof D3,
  container: HTMLElement,
  data: GroupedBarData,
  config: GroupedBarConfig
): void {
  if (data.categories.length === 0 || data.series.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const chartType = config.type ?? 'grouped'
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 440
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? (chartType === 'stacked' ? 'Stacked Bar Chart' : 'Grouped Bar Chart'),
    data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const categories = [...data.categories]
  const seriesLabels = data.series.map(s => s.label)

  // Build row data: one object per category
  const rows = categories.map((cat, ci) => {
    const obj: Record<string, number | string> = { cat }
    data.series.forEach(s => { obj[s.label] = s.values[ci] ?? 0 })
    return obj
  })

  let yMax: number
  if (chartType === 'stacked') {
    yMax = Math.max(...rows.map(row =>
      data.series.reduce((sum, s) => sum + (row[s.label] as number), 0)
    ))
  } else {
    yMax = Math.max(...data.series.flatMap(s => [...s.values].map(v => Math.abs(v))))
  }

  const xOuter = d3.scaleBand<string>().domain(categories).range([0, width]).padding(0.2)
  const xInner = d3.scaleBand<string>().domain(seriesLabels).range([0, xOuter.bandwidth()]).padding(0.05)
  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([height, 0]).nice()

  // Grid lines
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  if (chartType === 'stacked') {
    const stackGen = d3.stack<Record<string, number | string>>()
      .keys(seriesLabels)
      .value((row, key) => (row[key] as number) ?? 0)

    const stackedData = stackGen(rows)

    stackedData.forEach((layer, si) => {
      const color = getColor(si, theme)
      layer.forEach((seg, ci) => {
        const cat = categories[ci] ?? ''
        const bx = xOuter(cat) ?? 0
        const bw = xOuter.bandwidth()
        g.append('rect')
          .attr('x', bx).attr('width', bw)
          .attr('y', yScale(seg[1])).attr('height', Math.max(0, yScale(seg[0]) - yScale(seg[1])))
          .attr('fill', color).attr('opacity', theme.violinOpacity)
          .attr('stroke', theme.background).attr('stroke-width', 0.5)
          .on('mouseover', (event: MouseEvent) => {
            showTooltip(event, [
              formatTooltipRow('Category', cat),
              formatTooltipRow('Series', seriesLabels[si] ?? ''),
              formatTooltipRow('Value', ((seg[1] - seg[0])).toFixed(3)),
            ].join(''), theme)
          })
          .on('mouseout', hideTooltip)
      })
    })
  } else {
    // Grouped bars
    categories.forEach(cat => {
      const bxOuter = xOuter(cat) ?? 0
      data.series.forEach((s, si) => {
        const color = getColor(si, theme)
        const val = s.values[categories.indexOf(cat)] ?? 0
        const bx = bxOuter + (xInner(s.label) ?? 0)
        const bw = xInner.bandwidth()
        g.append('rect')
          .attr('x', bx).attr('width', bw)
          .attr('y', yScale(Math.max(0, val))).attr('height', Math.abs(yScale(0) - yScale(val)))
          .attr('fill', color).attr('opacity', theme.violinOpacity)
          .attr('stroke', color).attr('stroke-width', 0.5).attr('rx', 2)
          .on('mouseover', (event: MouseEvent) => {
            showTooltip(event, [
              formatTooltipRow('Category', cat),
              formatTooltipRow('Series', s.label),
              formatTooltipRow('Value', val.toFixed(3)),
            ].join(''), theme)
          })
          .on('mouseout', hideTooltip)
      })
    })
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xOuter))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Value')

  // Legend
  data.series.forEach((s, i) => {
    const color = getColor(i, theme)
    const lx = width - 130
    const ly = i * 18 + 8
    g.append('rect').attr('x', lx).attr('y', ly - 9).attr('width', 14).attr('height', 10)
      .attr('fill', color).attr('opacity', theme.violinOpacity).attr('rx', 2)
    g.append('text').attr('x', lx + 18).attr('y', ly)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(s.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
