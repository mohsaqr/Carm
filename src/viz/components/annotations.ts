/**
 * Statistical annotation helpers for D3 plots.
 * Renders APA text, regression equations, and stat result boxes.
 */

import type * as d3 from 'd3'
import type { CarmTheme } from '../themes/default.js'
import { DEFAULT_THEME } from '../themes/default.js'

type GSelection = d3.Selection<SVGGElement, unknown, null, undefined>

/** Add a subtitle line below the main title. */
export function addSubtitle(
  svg: d3.Selection<SVGSVGElement, unknown, null, undefined>,
  title: string,
  subtitle: string,
  width: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  // Title
  svg.append('text')
    .attr('x', width / 2)
    .attr('y', 20)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeTitle)
    .attr('font-weight', '600')
    .attr('fill', theme.text)
    .text(title)

  // Subtitle (stat result)
  svg.append('text')
    .attr('x', width / 2)
    .attr('y', 38)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamilyMono)
    .attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textMuted)
    .text(subtitle)
}

/** Add a caption at the bottom of the SVG. */
export function addCaption(
  svg: d3.Selection<SVGSVGElement, unknown, null, undefined>,
  text: string,
  width: number,
  height: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  svg.append('text')
    .attr('x', width / 2)
    .attr('y', height - 6)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted)
    .style('font-style', 'italic')
    .text(text)
}

/** Add a regression equation text on the plot. */
export function addRegressionEquation(
  g: GSelection,
  intercept: number,
  slope: number,
  r2: number,
  x: number,
  y: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  const sign = slope >= 0 ? '+' : '−'
  const eq = `ŷ = ${intercept.toFixed(2)} ${sign} ${Math.abs(slope).toFixed(2)}x, R² = ${r2.toFixed(2)}`
  g.append('text')
    .attr('x', x)
    .attr('y', y)
    .attr('font-family', theme.fontFamilyMono)
    .attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textAnnotation)
    .text(eq)
}

/** Add n= label to a group. */
export function addNLabel(
  g: GSelection,
  n: number,
  x: number,
  y: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  g.append('text')
    .attr('x', x)
    .attr('y', y)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textMuted)
    .text(`n = ${n}`)
}
