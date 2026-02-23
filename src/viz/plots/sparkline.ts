/**
 * Sparkline — minimal inline trend line for dashboard embedding.
 * Tiny margins, no axes, no labels, no grid. Just the line, an optional
 * shaded area, and a highlighted end-point dot.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import type { StatResult } from '../../core/types.js'

export interface SparklineConfig {
  readonly showArea?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface SparklineData {
  readonly values: readonly number[]
  readonly testResult?: StatResult
}

/**
 * Render a sparkline (minimal inline trend line).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - numeric time series values + optional stat result
 * @param config - visual configuration (very minimal by design)
 */
export function renderSparkline(
  container: HTMLElement,
  data: SparklineData,
  config: SparklineConfig = {}
): void {
  import('d3').then(d3 => renderSparklineD3(d3, container, data, config))
}

function renderSparklineD3(
  d3: typeof D3,
  container: HTMLElement,
  data: SparklineData,
  config: SparklineConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 200, 80)
  const H = config.height ?? 80
  const m = 8   // tiny uniform margin

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  const values = [...data.values]
  if (values.length < 2) {
    // Draw a flat dot for a single value
    if (values.length === 1) {
      const color = getColor(0, theme)
      svg.append('circle')
        .attr('cx', W / 2).attr('cy', H / 2)
        .attr('r', 3).attr('fill', color)
    }
    return
  }

  const color = getColor(0, theme)
  const [yMin, yMax] = d3.extent(values) as [number, number]
  const yPad = ((yMax - yMin) * 0.1) || 1

  const xScale = d3.scaleLinear()
    .domain([0, values.length - 1])
    .range([m, W - m])

  const yScale = d3.scaleLinear()
    .domain([yMin - yPad, yMax + yPad])
    .range([H - m, m])

  const indexed = values.map((v, i) => [i, v] as [number, number])

  // Optional area
  if (config.showArea !== false) {
    const areaFn = d3.area<[number, number]>()
      .x(d => xScale(d[0]))
      .y0(H - m)
      .y1(d => yScale(d[1]))
      .curve(d3.curveCatmullRom)

    svg.append('path')
      .datum(indexed)
      .attr('d', areaFn)
      .attr('fill', color)
      .attr('opacity', theme.ciOpacity)
  }

  // Line
  const lineFn = d3.line<[number, number]>()
    .x(d => xScale(d[0]))
    .y(d => yScale(d[1]))
    .curve(d3.curveCatmullRom)

  svg.append('path')
    .datum(indexed)
    .attr('d', lineFn)
    .attr('fill', 'none')
    .attr('stroke', color)
    .attr('stroke-width', 1.5)

  // End-point dot
  const lastVal = values[values.length - 1]!
  svg.append('circle')
    .attr('cx', xScale(values.length - 1))
    .attr('cy', yScale(lastVal))
    .attr('r', 3)
    .attr('fill', color)
    .attr('stroke', theme.background)
    .attr('stroke-width', 1.5)
}
