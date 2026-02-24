/**
 * Standalone box plot: Q1/Q3/median/whiskers (1.5 IQR)/outliers per group.
 * No violin. Uses quantile() from core/math.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import { quantile } from '../../core/math.js'
import type { StatResult } from '../../core/types.js'

export interface BoxplotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showN?: boolean
  readonly showMean?: boolean
  readonly showOutliers?: boolean
  readonly showMedian?: boolean
}

export interface BoxplotGroup {
  readonly label: string
  readonly values: readonly number[]
}

export interface BoxplotData {
  readonly groups: readonly BoxplotGroup[]
  readonly testResult?: StatResult
}

export function renderBoxplot(
  container: HTMLElement,
  data: BoxplotData,
  config: BoxplotConfig = {}
): void {
  import('d3').then(d3 => renderBoxplotD3(d3, container, data, config))
}

function computeBoxStats(values: readonly number[]) {
  const q1 = quantile(values, 0.25)
  const med = quantile(values, 0.5)
  const q3 = quantile(values, 0.75)
  const iqr = q3 - q1
  const fence_lo = q1 - 1.5 * iqr
  const fence_hi = q3 + 1.5 * iqr
  const whisker_lo = Math.min(...values.filter(v => v >= fence_lo))
  const whisker_hi = Math.max(...values.filter(v => v <= fence_hi))
  const outliers = values.filter(v => v < fence_lo || v > fence_hi)
  return { q1, med, q3, whisker_lo, whisker_hi, outliers }
}

function renderBoxplotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: BoxplotData,
  config: BoxplotConfig
): void {
  if (data.groups.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 440
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Box Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const labels = data.groups.map(gr => gr.label)
  const allValues = data.groups.flatMap(gr => [...gr.values])
  const [yMin, yMax] = d3.extent(allValues) as [number, number]
  const yPad = (yMax - yMin) * 0.1

  const xScale = d3.scaleBand<string>().domain(labels).range([0, width]).padding(0.3)
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice()

  // Grid lines
  g.selectAll('.grid').data(yScale.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  data.groups.forEach((gr, gi) => {
    if (gr.values.length === 0) return
    const color = getColor(gi, theme)
    const bx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2
    const boxW = xScale.bandwidth() * 0.5
    const { q1, med, q3, whisker_lo, whisker_hi, outliers } = computeBoxStats(gr.values)

    // Whisker lines
    g.append('line')
      .attr('x1', bx).attr('x2', bx)
      .attr('y1', yScale(whisker_lo)).attr('y2', yScale(q1))
      .attr('stroke', color).attr('stroke-width', 1.5).attr('stroke-dasharray', '3,2')
    g.append('line')
      .attr('x1', bx).attr('x2', bx)
      .attr('y1', yScale(q3)).attr('y2', yScale(whisker_hi))
      .attr('stroke', color).attr('stroke-width', 1.5).attr('stroke-dasharray', '3,2')

    // Whisker caps
    const capW = boxW * 0.4
    g.append('line').attr('x1', bx - capW / 2).attr('x2', bx + capW / 2)
      .attr('y1', yScale(whisker_lo)).attr('y2', yScale(whisker_lo))
      .attr('stroke', color).attr('stroke-width', 1.5)
    g.append('line').attr('x1', bx - capW / 2).attr('x2', bx + capW / 2)
      .attr('y1', yScale(whisker_hi)).attr('y2', yScale(whisker_hi))
      .attr('stroke', color).attr('stroke-width', 1.5)

    // IQR box
    g.append('rect')
      .attr('x', bx - boxW / 2).attr('width', boxW)
      .attr('y', yScale(q3)).attr('height', Math.max(1, yScale(q1) - yScale(q3)))
      .attr('fill', color).attr('opacity', theme.violinOpacity)
      .attr('stroke', color).attr('stroke-width', 2).attr('rx', 2)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('Group', gr.label),
          formatTooltipRow('Median', med.toFixed(3)),
          formatTooltipRow('Q1', q1.toFixed(3)),
          formatTooltipRow('Q3', q3.toFixed(3)),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Median line + value
    if (config.showMedian !== false) {
      g.append('line')
        .attr('x1', bx - boxW / 2).attr('x2', bx + boxW / 2)
        .attr('y1', yScale(med)).attr('y2', yScale(med))
        .attr('stroke', theme.background).attr('stroke-width', 2.5)
      g.append('text')
        .attr('x', bx + boxW / 2 + 4).attr('y', yScale(med) + 3.5)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textAnnotation).text(med.toFixed(2))
    }

    // Mean diamond marker + value
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length
      const my = yScale(groupMean)
      const ds = 5
      g.append('polygon')
        .attr('points', `${bx},${my - ds} ${bx + ds},${my} ${bx},${my + ds} ${bx - ds},${my}`)
        .attr('fill', 'white').attr('stroke', color).attr('stroke-width', 1.5)
      g.append('text')
        .attr('x', bx + ds + 4).attr('y', my + 3.5)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textAnnotation).text(groupMean.toFixed(2))
    }

    // Outliers
    if (config.showOutliers !== false) {
      outliers.forEach((v) => {
        g.append('circle')
          .attr('cx', bx).attr('cy', yScale(v))
          .attr('r', 3.5).attr('fill', 'none').attr('stroke', color).attr('stroke-width', 1.5)
          .on('mouseover', (event: MouseEvent) => {
            showTooltip(event, [
              formatTooltipRow('Group', gr.label),
              formatTooltipRow('Outlier', v.toFixed(3)),
            ].join(''), theme)
          })
          .on('mouseout', hideTooltip)
      })
    }

    // n label
    if (config.showN !== false) {
      g.append('text').attr('x', bx).attr('y', height + 44)
        .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall).attr('fill', theme.textMuted)
        .text(`n = ${gr.values.length}`)
    }
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 60)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Value')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
