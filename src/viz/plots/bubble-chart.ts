/**
 * Bubble chart: scatter plot where bubble size is proportional to r value.
 * Color by group if group field is present. Tooltip on hover.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface BubbleChartConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface BubblePoint {
  readonly x: number
  readonly y: number
  readonly r: number
  readonly label?: string
  readonly group?: string
}

export interface BubbleChartData {
  readonly points: readonly BubblePoint[]
  readonly testResult?: StatResult
}

export function renderBubbleChart(
  container: HTMLElement,
  data: BubbleChartData,
  config: BubbleChartConfig = {}
): void {
  import('d3').then(d3 => renderBubbleChartD3(d3, container, data, config))
}

function renderBubbleChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: BubbleChartData,
  config: BubbleChartConfig
): void {
  if (data.points.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 460
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Bubble Chart', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const xVals = data.points.map(p => p.x)
  const yVals = data.points.map(p => p.y)
  const rVals = data.points.map(p => p.r)
  const [xMin, xMax] = d3.extent(xVals) as [number, number]
  const [yMin, yMax] = d3.extent(yVals) as [number, number]
  const [rMin, rMax] = d3.extent(rVals) as [number, number]
  const xPad = (xMax - xMin) * 0.1
  const yPad = (yMax - yMin) * 0.1

  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice()
  // Map r values to pixel radius [4, 30]
  const rScale = rMax === rMin
    ? (_: number) => 12
    : d3.scaleSqrt().domain([rMin, rMax]).range([4, 30])

  // Unique groups for color mapping
  const groups = [...new Set(data.points.map(p => p.group ?? '__default__'))]
  const groupIndex = (grp: string) => groups.indexOf(grp)
  const hasGroups = !(groups.length === 1 && groups[0] === '__default__')

  // Grid lines
  g.selectAll('.gridH').data(yScale.ticks(5)).join('line')
    .attr('class', 'gridH').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Bubbles (sorted by r desc so smaller ones render on top)
  const sorted = [...data.points].sort((a, b) => b.r - a.r)

  sorted.forEach(pt => {
    const grp = pt.group ?? '__default__'
    const color = getColor(groupIndex(grp), theme)
    const cx = xScale(pt.x)
    const cy = yScale(pt.y)
    const radius = rScale(pt.r)

    g.append('circle')
      .attr('cx', cx).attr('cy', cy).attr('r', radius)
      .attr('fill', color).attr('opacity', theme.pointOpacity)
      .attr('stroke', color).attr('stroke-width', 1.5).attr('stroke-opacity', 0.8)
      .on('mouseover', (event: MouseEvent) => {
        const rows = [
          formatTooltipRow('x', pt.x.toFixed(3)),
          formatTooltipRow('y', pt.y.toFixed(3)),
          formatTooltipRow('size', pt.r.toFixed(3)),
        ]
        if (pt.label) rows.unshift(formatTooltipRow('Label', pt.label))
        if (hasGroups && pt.group) rows.push(formatTooltipRow('Group', pt.group))
        showTooltip(event, rows.join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Optional label inside large bubbles
    if (pt.label && radius > 16) {
      g.append('text')
        .attr('x', cx).attr('y', cy + 4)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.background).attr('pointer-events', 'none')
        .text(pt.label.slice(0, 6))
    }
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'x')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'y')

  // Group legend
  if (hasGroups) {
    groups.forEach((grp, i) => {
      const color = getColor(i, theme)
      const lx = width - 130
      const ly = i * 18 + 8
      g.append('circle').attr('cx', lx + 6).attr('cy', ly).attr('r', 5)
        .attr('fill', color).attr('opacity', theme.pointOpacity)
      g.append('text').attr('x', lx + 16).attr('y', ly + 4)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text).text(grp)
    })
  }

  // Size legend (show 3 representative sizes if r values vary)
  if (rMax > rMin) {
    const sizeLegendVals = [rMin, (rMin + rMax) / 2, rMax]
    const slx = 16
    let sly = height - 70
    g.append('text').attr('x', slx).attr('y', sly - 6)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted).text('Size:')
    sizeLegendVals.forEach(rv => {
      const rad = rScale(rv)
      g.append('circle').attr('cx', slx + 16).attr('cy', sly + rad)
        .attr('r', rad).attr('fill', 'none')
        .attr('stroke', theme.textMuted).attr('stroke-width', 1)
      g.append('text').attr('x', slx + 36).attr('y', sly + rad + 4)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textMuted).text(rv.toFixed(1))
      sly += rad * 2 + 10
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
