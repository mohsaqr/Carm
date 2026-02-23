/**
 * Pie / donut chart with polyline labels and optional test result subtitle.
 * Donut mode uses innerRadius = outerRadius * 0.55.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface PieChartConfig {
  readonly donut?: boolean
  readonly showPercentages?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface PieSlice {
  readonly label: string
  readonly value: number
}

export interface PieChartData {
  readonly slices: readonly PieSlice[]
  readonly testResult?: StatResult
}

export function renderPieChart(
  container: HTMLElement,
  data: PieChartData,
  config: PieChartConfig = {}
): void {
  import('d3').then(d3 => renderPieChartD3(d3, container, data, config))
}

function renderPieChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: PieChartData,
  config: PieChartConfig
): void {
  if (data.slices.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Pie Chart', data.testResult?.formatted ?? '', W, theme)

  const legendH = Math.ceil(data.slices.length / 3) * 22
  const plotH = H - margin.top - margin.bottom - legendH
  const cx = W / 2
  const cy = margin.top + plotH / 2

  const outerR = Math.min(W / 2, plotH / 2) * 0.72
  const innerR = (config.donut ?? false) ? outerR * 0.55 : 0
  const total = data.slices.reduce((s, sl) => s + sl.value, 0)

  const pieGen = d3.pie<PieSlice>().value(d => d.value).sort(null)
  const arcGen = d3.arc<D3.PieArcDatum<PieSlice>>()
    .innerRadius(innerR).outerRadius(outerR)
  const labelArc = d3.arc<D3.PieArcDatum<PieSlice>>()
    .innerRadius(outerR * 1.15).outerRadius(outerR * 1.15)
  const polyArc = d3.arc<D3.PieArcDatum<PieSlice>>()
    .innerRadius(outerR * 0.85).outerRadius(outerR * 0.85)

  const arcs = pieGen(data.slices as PieSlice[])
  const g = svg.append('g').attr('transform', `translate(${cx},${cy})`)

  // Slices
  arcs.forEach((arc, i) => {
    const color = getColor(i, theme)
    const pct = ((arc.data.value / total) * 100).toFixed(1)

    g.append('path')
      .attr('d', arcGen(arc) ?? '')
      .attr('fill', color)
      .attr('opacity', 0.88)
      .attr('stroke', theme.background)
      .attr('stroke-width', 2)
      .on('mouseover', function(event: MouseEvent) {
        showTooltip(event, [
          formatTooltipRow('Label', arc.data.label),
          formatTooltipRow('Value', arc.data.value.toFixed(2)),
          formatTooltipRow('Percentage', pct + '%')
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Polyline label
    const midAngle = (arc.startAngle + arc.endAngle) / 2
    const isRight = midAngle < Math.PI
    const labelPt = labelArc.centroid(arc)
    const polyPt = polyArc.centroid(arc)
    const endX = isRight ? outerR * 1.35 : -outerR * 1.35

    // Only draw label if slice is big enough
    if (arc.data.value / total > 0.03) {
      g.append('polyline')
        .attr('points', [polyPt, labelPt, [endX, labelPt[1]]].map(p => p.join(',')).join(' '))
        .attr('fill', 'none')
        .attr('stroke', theme.textMuted)
        .attr('stroke-width', 1)

      const displayStr = (config.showPercentages ?? true)
        ? `${arc.data.label} (${pct}%)`
        : arc.data.label

      g.append('text')
        .attr('x', isRight ? endX + 4 : endX - 4)
        .attr('y', labelPt[1])
        .attr('dominant-baseline', 'middle')
        .attr('text-anchor', isRight ? 'start' : 'end')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text)
        .text(displayStr)
    }
  })

  // Donut center label
  if (config.donut ?? false) {
    g.append('text').attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize + 2)
      .attr('font-weight', '600').attr('fill', theme.text)
      .text(`n = ${total.toFixed(0)}`)
  }

  // Legend
  const legendY = margin.top + plotH + 12
  const colW = W / 3
  data.slices.forEach((sl, i) => {
    const col = i % 3
    const row = Math.floor(i / 3)
    const lx = col * colW + 16
    const ly = legendY + row * 22

    svg.append('rect').attr('x', lx).attr('y', ly).attr('width', 12).attr('height', 12)
      .attr('fill', getColor(i, theme)).attr('rx', 2)
    svg.append('text').attr('x', lx + 16).attr('y', ly + 9)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(sl.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
