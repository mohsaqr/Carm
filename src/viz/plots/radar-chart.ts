/**
 * Radar / spider chart — one polygon per series, one spoke per axis.
 * Values are normalized per axis to [0,1] using min/max across all series.
 * Spokes are uniformly spaced angularly. Filled polygons with opacity.
 * Axis labels at spoke tips. Legend for series.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface RadarChartConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly levels?: number   // number of concentric rings
}

export interface RadarChartData {
  readonly series: readonly {
    readonly label: string
    readonly values: readonly number[]
  }[]
  readonly axes: readonly string[]
  readonly testResult?: StatResult
}

/**
 * Render a radar / spider chart.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - series (each with values per axis) + axis names + optional stat result
 * @param config - visual configuration
 */
export function renderRadarChart(
  container: HTMLElement,
  data: RadarChartData,
  config: RadarChartConfig = {}
): void {
  import('d3').then(d3 => renderRadarChartD3(d3, container, data, config))
}

function renderRadarChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: RadarChartData,
  config: RadarChartConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const levels = config.levels ?? 5

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Radar Chart', data.testResult?.formatted ?? '', W, theme)

  const nAxes = data.axes.length
  if (nAxes < 3 || data.series.length === 0) return

  // Reserve space for legend on the right
  const legendWidth = 120
  const chartCentreX = (W - legendWidth) / 2
  const chartCentreY = (H - theme.marginTop - 32) / 2 + theme.marginTop + 8
  const radius = Math.min((W - legendWidth) / 2, (H - theme.marginTop - 48) / 2) * 0.72

  const g = svg.append('g')

  // Per-axis min/max (across all series) for normalisation
  const axisMin = Array.from({ length: nAxes }, (_, ai) =>
    Math.min(...data.series.map(s => s.values[ai] ?? 0))
  )
  const axisMax = Array.from({ length: nAxes }, (_, ai) =>
    Math.max(...data.series.map(s => s.values[ai] ?? 0))
  )

  const normalize = (val: number, ai: number): number => {
    const lo = axisMin[ai]!
    const hi = axisMax[ai]!
    return hi === lo ? 0.5 : (val - lo) / (hi - lo)
  }

  // Angle for each axis (start top, clockwise)
  const angleOf = (i: number): number => (i / nAxes) * 2 * Math.PI - Math.PI / 2

  // Point on radar at normalised radius t for axis ai
  const pt = (t: number, ai: number): [number, number] => {
    const a = angleOf(ai)
    return [
      chartCentreX + radius * t * Math.cos(a),
      chartCentreY + radius * t * Math.sin(a),
    ]
  }

  // Concentric level rings
  for (let lvl = 1; lvl <= levels; lvl++) {
    const t = lvl / levels
    const ringPts = Array.from({ length: nAxes }, (_, ai) => pt(t, ai))
    g.append('polygon')
      .attr('points', ringPts.map(([x, y]) => `${x},${y}`).join(' '))
      .attr('fill', 'none')
      .attr('stroke', theme.gridLine)
      .attr('stroke-width', 1)
  }

  // Spokes
  Array.from({ length: nAxes }, (_, ai) => {
    const [x2, y2] = pt(1, ai)
    g.append('line')
      .attr('x1', chartCentreX).attr('y1', chartCentreY)
      .attr('x2', x2).attr('y2', y2)
      .attr('stroke', theme.gridLine)
      .attr('stroke-width', 1)
  })

  // Axis labels at spoke tips
  data.axes.forEach((label, ai) => {
    const angle = angleOf(ai)
    const labelDist = radius * 1.15
    const lx = chartCentreX + labelDist * Math.cos(angle)
    const ly = chartCentreY + labelDist * Math.sin(angle)

    let anchor = 'middle'
    if (Math.cos(angle) > 0.1) anchor = 'start'
    else if (Math.cos(angle) < -0.1) anchor = 'end'

    g.append('text')
      .attr('x', lx).attr('y', ly + 4)
      .attr('text-anchor', anchor)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(label)
  })

  // Draw each series polygon
  data.series.forEach((series, si) => {
    const color = getColor(si, theme)
    const polyPts = data.axes.map((_, ai) => {
      const t = normalize(series.values[ai] ?? 0, ai)
      return pt(t, ai)
    })

    g.append('polygon')
      .attr('points', polyPts.map(([x, y]) => `${x},${y}`).join(' '))
      .attr('fill', color)
      .attr('fill-opacity', theme.ciOpacity + 0.05)
      .attr('stroke', color)
      .attr('stroke-width', 2)

    // Vertex dots with tooltip
    polyPts.forEach(([px, py], ai) => {
      g.append('circle')
        .attr('cx', px).attr('cy', py)
        .attr('r', 4)
        .attr('fill', color)
        .attr('stroke', theme.background)
        .attr('stroke-width', 1.5)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Series', series.label),
            formatTooltipRow('Axis', data.axes[ai] ?? `Axis ${ai}`),
            formatTooltipRow('Value', (series.values[ai] ?? 0).toFixed(3)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })
  })

  // Legend (right side)
  const lgX = W - legendWidth + 8
  data.series.forEach((series, si) => {
    const color = getColor(si, theme)
    const lgY = theme.marginTop + si * 22

    g.append('rect')
      .attr('x', lgX).attr('y', lgY)
      .attr('width', 12).attr('height', 12)
      .attr('fill', color).attr('opacity', 0.8)
      .attr('rx', 2)

    g.append('text')
      .attr('x', lgX + 16).attr('y', lgY + 10)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(series.label.length > 14 ? series.label.slice(0, 13) + '…' : series.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
