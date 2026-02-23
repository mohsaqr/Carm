/**
 * Smart axis rendering for Carm plots.
 * Auto tick formatting, minimal grid, clean style.
 */

import type * as d3 from 'd3'
import type { CarmTheme } from '../themes/default.js'
import { DEFAULT_THEME } from '../themes/default.js'

type AnyScale = d3.ScaleLinear<number, number> | d3.ScaleBand<string>

/** Render a styled x-axis with label. */
export function renderXAxis(
  g: d3.Selection<SVGGElement, unknown, null, undefined>,
  _scale: AnyScale,
  height: number,
  label: string,
  width: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  g.append('g')
    .attr('transform', `translate(0,${height})`)
    .call((g_: d3.Selection<SVGGElement, unknown, null, undefined>) => {
      g_.selectAll('line').attr('stroke', theme.axisLine)
      g_.selectAll('path').attr('stroke', theme.axisLine)
      g_.selectAll('text')
        .attr('fill', theme.text)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSize)
    })

  // Axis label
  g.append('text')
    .attr('x', width / 2)
    .attr('y', height + 48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(label)
}

/** Render a styled y-axis with label. */
export function renderYAxis(
  g: d3.Selection<SVGGElement, unknown, null, undefined>,
  _scale: d3.ScaleLinear<number, number>,
  height: number,
  label: string,
  theme: CarmTheme = DEFAULT_THEME
): void {
  g.append('g')
    .call(d3Axis => {
      // Will be called with actual d3 in the plot
      void d3Axis
    })

  // Y-axis label
  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -52)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(label)
}

/** Render horizontal grid lines. */
export function renderGridLines(
  g: d3.Selection<SVGGElement, unknown, null, undefined>,
  scale: d3.ScaleLinear<number, number>,
  width: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  const ticks = scale.ticks(6)
  g.selectAll('.grid-line')
    .data(ticks)
    .join('line')
    .attr('class', 'grid-line')
    .attr('x1', 0)
    .attr('x2', width)
    .attr('y1', scale)
    .attr('y2', scale)
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1)
}
