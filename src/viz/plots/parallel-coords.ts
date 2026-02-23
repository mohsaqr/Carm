/**
 * Parallel coordinates plot — one vertical axis per variable, evenly spaced
 * horizontally. Each data row is a polyline through its per-axis values.
 * Axes are individually min-max scaled. Rows are optionally coloured by group.
 * Lines are semi-transparent to show density.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ParallelCoordsConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface ParallelCoordsData {
  /** rows[obsIndex][varIndex] */
  readonly rows: readonly (readonly number[])[]
  readonly axes: readonly string[]
  /** optional group index per row (0-based) for colour coding */
  readonly groups?: readonly number[]
  readonly testResult?: StatResult
}

/**
 * Render a parallel coordinates plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - observations (rows × vars), axis names, optional groups + stat result
 * @param config - visual configuration
 */
export function renderParallelCoords(
  container: HTMLElement,
  data: ParallelCoordsData,
  config: ParallelCoordsConfig = {}
): void {
  import('d3').then(d3 => renderParallelCoordsD3(d3, container, data, config))
}

function renderParallelCoordsD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ParallelCoordsData,
  config: ParallelCoordsConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 700, 400)
  const H = config.height ?? 480
  const margin = {
    top: theme.marginTop + 16,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft,
  }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Parallel Coordinates', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const nAxes = data.axes.length
  if (nAxes < 2 || data.rows.length === 0) return

  // X position of each axis
  const axisX = (i: number): number => (i / (nAxes - 1)) * width

  // Per-axis linear scales (min–max normalised to plot height)
  const yScales = data.axes.map((_, ai) => {
    const vals = data.rows.map(row => row[ai] ?? 0)
    const [lo, hi] = d3.extent(vals) as [number, number]
    const pad = ((hi - lo) * 0.05) || 0.5
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([height, 0]).nice()
  })

  // Draw polyline per row
  data.rows.forEach((row, ri) => {
    const groupIdx = data.groups?.[ri] ?? 0
    const color = getColor(groupIdx, theme)

    const pathData = data.axes.map((_, ai) => {
      const val = row[ai] ?? 0
      return [axisX(ai), yScales[ai]!(val)] as [number, number]
    })

    g.append('path')
      .attr('d', lineThrough(pathData))
      .attr('fill', 'none')
      .attr('stroke', color)
      .attr('stroke-width', 1.2)
      .attr('opacity', Math.max(0.1, Math.min(0.55, 30 / data.rows.length)))
      .on('mouseover', (event: MouseEvent) => {
        const tooltipRows = data.axes.map((label, ai) =>
          formatTooltipRow(label, (row[ai] ?? 0).toFixed(3))
        )
        showTooltip(event, tooltipRows.join(''), theme)
      })
      .on('mouseout', hideTooltip)
  })

  // Draw vertical axes and labels
  data.axes.forEach((label, ai) => {
    const x = axisX(ai)
    const ySc = yScales[ai]!

    // Axis line
    g.append('line')
      .attr('x1', x).attr('x2', x)
      .attr('y1', 0).attr('y2', height)
      .attr('stroke', theme.axisLine)
      .attr('stroke-width', 1.5)

    // Tick marks (5 ticks)
    const ticks = ySc.ticks(5)
    ticks.forEach(tick => {
      g.append('line')
        .attr('x1', x - 4).attr('x2', x + 4)
        .attr('y1', ySc(tick)).attr('y2', ySc(tick))
        .attr('stroke', theme.axisLine)
        .attr('stroke-width', 1)
      g.append('text')
        .attr('x', x - 6).attr('y', ySc(tick) + 4)
        .attr('text-anchor', 'end')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textMuted)
        .text(tick.toFixed(1))
    })

    // Axis label at top
    g.append('text')
      .attr('x', x).attr('y', -10)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('font-weight', '600')
      .attr('fill', theme.text)
      .text(label)
  })

  // Legend for groups (if present)
  if (data.groups) {
    const uniqueGroups = Array.from(new Set(data.groups)).sort((a, b) => a - b)
    uniqueGroups.forEach((grp, i) => {
      const lgX = 0
      const lgY = height + 32 + i * 18
      g.append('rect')
        .attr('x', lgX).attr('y', lgY)
        .attr('width', 10).attr('height', 10)
        .attr('fill', getColor(grp, theme))
        .attr('rx', 2)
      g.append('text')
        .attr('x', lgX + 14).attr('y', lgY + 9)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text)
        .text(`Group ${grp}`)
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Helper ───────────────────────────────────────────────────────────────────

/** Build an SVG path string connecting points with straight segments. */
function lineThrough(pts: ReadonlyArray<readonly [number, number]>): string {
  if (pts.length === 0) return ''
  return pts.map(([x, y], i) => `${i === 0 ? 'M' : 'L'}${x.toFixed(2)},${y.toFixed(2)}`).join(' ')
}
