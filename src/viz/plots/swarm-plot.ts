/**
 * Beeswarm plot — individual points per group arranged to avoid overlap.
 * Points are sorted by value and offset horizontally so they spread
 * left/right of the group centre axis without colliding.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface SwarmPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly pointRadius?: number
  readonly bandwidth?: number
  readonly showN?: boolean
  readonly showMean?: boolean
}

export interface SwarmPlotData {
  readonly groups: readonly { readonly label: string; readonly values: readonly number[] }[]
  readonly testResult?: StatResult
}

/**
 * Render a beeswarm plot with per-group colour and collision-avoidance layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - groups of numeric values + optional stat result
 * @param config - visual configuration
 */
export function renderSwarmPlot(
  container: HTMLElement,
  data: SwarmPlotData,
  config: SwarmPlotConfig = {}
): void {
  import('d3').then(d3 => renderSwarmPlotD3(d3, container, data, config))
}

function renderSwarmPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: SwarmPlotData,
  config: SwarmPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft,
  }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom
  const r = config.pointRadius ?? 4

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Beeswarm Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // Scales
  const labels = data.groups.map(gr => gr.label)
  const xScale = d3.scaleBand<string>()
    .domain(labels)
    .range([0, width])
    .padding(0.3)

  const allValues = data.groups.flatMap(gr => [...gr.values])
  const [yMin, yMax] = d3.extent(allValues) as [number, number]
  const yPad = (yMax - yMin) * 0.08 || 1
  const yScale = d3.scaleLinear()
    .domain([yMin - yPad, yMax + yPad])
    .range([height, 0])
    .nice()

  // Horizontal grid
  g.selectAll<SVGLineElement, number>('.grid')
    .data(yScale.ticks(6))
    .join('line')
    .attr('class', 'grid')
    .attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1)

  // Axes
  g.append('g')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(xScale))
    .call(ax => ax.select('.domain').attr('stroke', theme.axisLine))
    .selectAll<SVGTextElement, string>('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  g.append('g')
    .call(d3.axisLeft(yScale).ticks(6))
    .call(ax => ax.select('.domain').attr('stroke', theme.axisLine))
    .selectAll<SVGTextElement, number>('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  // Axis labels
  g.append('text')
    .attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.xLabel ?? '')

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.yLabel ?? 'Value')

  // Draw beeswarm per group
  data.groups.forEach((gr, gi) => {
    const color = getColor(gi, theme)
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2
    const maxHalfWidth = xScale.bandwidth() / 2 - r

    // Compute beeswarm positions using bucket collision-avoidance
    const positions = beeswarm([...gr.values], yScale, r, maxHalfWidth)

    positions.forEach(({ value, xOffset }) => {
      g.append('circle')
        .attr('cx', cx + xOffset)
        .attr('cy', yScale(value))
        .attr('r', r)
        .attr('fill', color)
        .attr('opacity', theme.pointOpacity)
        .attr('stroke', theme.background)
        .attr('stroke-width', 0.5)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Group', gr.label),
            formatTooltipRow('Value', value.toFixed(4)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })

    // Mean diamond marker + value
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length
      const my = yScale(groupMean)
      const ds = 5
      g.append('polygon')
        .attr('points', `${cx},${my - ds} ${cx + ds},${my} ${cx},${my + ds} ${cx - ds},${my}`)
        .attr('fill', 'white').attr('stroke', color).attr('stroke-width', 1.5)
      g.append('text')
        .attr('x', cx + ds + 4).attr('y', my + 3.5)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textAnnotation).text(groupMean.toFixed(2))
    }

    // n label below x-axis
    if (config.showN !== false) {
      svg.append('text')
        .attr('x', margin.left + cx)
        .attr('y', H - theme.marginBottom + 28)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.textMuted)
        .text(`n = ${gr.values.length}`)
    }
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Beeswarm layout ─────────────────────────────────────────────────────────

/**
 * Compute x-offsets for a set of values so that no two circles overlap.
 * Algorithm: sort by y-pixel position, then greedily place each point at the
 * smallest valid x-offset (trying 0, ±step, ±2step, …).
 */
function beeswarm(
  values: number[],
  yScale: { (v: number): number },
  r: number,
  maxHalf: number
): Array<{ value: number; xOffset: number }> {
  // Sort values ascending (bottom-up on screen when y increases downward)
  const sorted = values.slice().sort((a, b) => a - b)
  const placed: Array<{ y: number; xOffset: number }> = []
  const diameter = r * 2 + 0.5  // collision diameter with tiny gap

  const result = sorted.map(value => {
    const y = yScale(value)
    let xOffset = 0
    let placed_flag = false

    // Try offsets in order: 0, +step, -step, +2step, -2step, …
    for (let attempt = 0; attempt <= Math.ceil(maxHalf / diameter) + 1; attempt++) {
      const candidates = attempt === 0 ? [0] : [attempt * diameter, -attempt * diameter]
      for (const xOff of candidates) {
        if (Math.abs(xOff) > maxHalf) continue
        const overlaps = placed.some(p => {
          const dx = xOff - p.xOffset
          const dy = y - p.y
          return Math.sqrt(dx * dx + dy * dy) < diameter
        })
        if (!overlaps) {
          xOffset = xOff
          placed_flag = true
          break
        }
      }
      if (placed_flag) break
    }

    placed.push({ y, xOffset })
    return { value, xOffset }
  })

  return result
}
