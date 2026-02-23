/**
 * Cleveland dot plot: horizontal layout with labels on y-axis.
 * Supports single dots or paired dots connected by a line (group comparison).
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface DotPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface DotPlotData {
  readonly labels: readonly string[]
  readonly values: readonly number[]
  readonly group2?: readonly number[]
  readonly group1Label?: string
  readonly group2Label?: string
  readonly testResult?: StatResult
}

export function renderDotPlot(
  container: HTMLElement,
  data: DotPlotData,
  config: DotPlotConfig = {}
): void {
  import('d3').then(d3 => renderDotPlotD3(d3, container, data, config))
}

function renderDotPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: DotPlotData,
  config: DotPlotConfig
): void {
  if (data.labels.length === 0 || data.values.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const itemCount = data.labels.length
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? Math.max(200, itemCount * 36 + theme.marginTop + theme.marginBottom)
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft + 40 }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom
  const isPaired = data.group2 != null && data.group2.length === data.labels.length

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Dot Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allValues: number[] = [...data.values]
  if (isPaired) allValues.push(...data.group2!)
  const [vMin, vMax] = d3.extent(allValues) as [number, number]
  const xPad = (vMax - vMin) * 0.1

  const xScale = d3.scaleLinear().domain([vMin - xPad, vMax + xPad]).range([0, width]).nice()
  const yScale = d3.scaleBand<string>().domain([...data.labels]).range([0, height]).padding(0.35)

  // Grid lines (vertical)
  g.selectAll('.grid').data(xScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', d => xScale(d)).attr('x2', d => xScale(d))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  const color1 = getColor(0, theme)
  const color2 = getColor(1, theme)

  data.labels.forEach((label, i) => {
    const v1 = data.values[i] ?? 0
    const v2 = isPaired ? (data.group2![i] ?? 0) : null
    const cy = (yScale(label) ?? 0) + yScale.bandwidth() / 2

    // Connecting line for paired
    if (isPaired && v2 !== null) {
      g.append('line')
        .attr('x1', xScale(v1)).attr('x2', xScale(v2))
        .attr('y1', cy).attr('y2', cy)
        .attr('stroke', theme.gridLine).attr('stroke-width', 2)
    }

    // Group 1 dot
    g.append('circle')
      .attr('cx', xScale(v1)).attr('cy', cy).attr('r', 6)
      .attr('fill', color1).attr('stroke', theme.background).attr('stroke-width', 1.5)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('Label', label),
          formatTooltipRow(data.group1Label ?? 'Value', v1.toFixed(3)),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Group 2 dot (paired)
    if (isPaired && v2 !== null) {
      g.append('circle')
        .attr('cx', xScale(v2)).attr('cy', cy).attr('r', 6)
        .attr('fill', color2).attr('stroke', theme.background).attr('stroke-width', 1.5)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Label', label),
            formatTooltipRow(data.group2Label ?? 'Group 2', v2.toFixed(3)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    }
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

  // Legend for paired
  if (isPaired) {
    const legendX = width - 140
    const label1 = data.group1Label ?? 'Group 1'
    const label2 = data.group2Label ?? 'Group 2'
    g.append('circle').attr('cx', legendX + 6).attr('cy', 8).attr('r', 5).attr('fill', color1)
    g.append('text').attr('x', legendX + 16).attr('y', 12)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(label1)
    g.append('circle').attr('cx', legendX + 6).attr('cy', 26).attr('r', 5).attr('fill', color2)
    g.append('text').attr('x', legendX + 16).attr('y', 30)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(label2)
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
