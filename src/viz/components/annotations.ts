/**
 * Statistical annotation helpers for D3 plots.
 * Renders APA text, regression equations, and stat result boxes.
 */

import type * as d3 from 'd3'
import type { CarmTheme } from '../themes/default.js'
import { DEFAULT_THEME } from '../themes/default.js'

type SVGSelection = d3.Selection<SVGSVGElement, unknown, null, undefined>
type GSelection = d3.Selection<SVGGElement, unknown, null, undefined>

/**
 * Add a title + statistical subtitle at the top of an SVG.
 *
 * Layout (y positions, px from SVG top):
 *   20 — decorative accent bar (2×20 px, theme color)
 *   22 — title baseline  (16 px bold, near-black)
 *   43 — subtitle baseline (11 px italic sans-serif, slate)
 *
 * Theme marginTop should be ≥ 58 so the plot area starts below the subtitle.
 */
export function addSubtitle(
  svg: SVGSelection,
  title: string,
  subtitle: string,
  _width: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  // Slim accent bar — left edge decoration
  svg.append('rect')
    .attr('x', 20)
    .attr('y', 10)
    .attr('width', 3)
    .attr('height', 22)
    .attr('rx', 1.5)
    .attr('fill', theme.colors[0] ?? '#4e79a7')

  // Title — left-aligned (indented past the accent bar)
  svg.append('text')
    .attr('x', 30)
    .attr('y', 26)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeTitle)
    .attr('font-weight', '700')
    .attr('letter-spacing', '-0.3')
    .attr('fill', theme.text)
    .text(title)

  // Subtitle — left-aligned, italic, regular (not mono), slate color
  if (subtitle) {
    svg.append('text')
      .attr('x', 30)
      .attr('y', 45)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('font-style', 'italic')
      .attr('fill', theme.textAnnotation)
      .text(subtitle)
  }
}

/**
 * Add an italic caption at the bottom-left of the SVG.
 * Used for data source, method notes, sample size.
 */
export function addCaption(
  svg: SVGSelection,
  text: string,
  _width: number,
  height: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  svg.append('text')
    .attr('x', 20)
    .attr('y', height - 8)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted)
    .style('font-style', 'italic')
    .text(text)
}

/**
 * Add a regression equation text annotation on the plot area.
 * Uses monospace font so numbers align cleanly.
 */
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
  const eq = `ŷ = ${intercept.toFixed(2)} ${sign} ${Math.abs(slope).toFixed(2)}x   R² = ${r2.toFixed(3)}`

  // Pill background
  const padX = 8, padY = 4
  const textW = eq.length * 6.2 + padX * 2
  const textH = 14 + padY * 2
  g.append('rect')
    .attr('x', x - padX)
    .attr('y', y - textH + padY)
    .attr('width', textW)
    .attr('height', textH)
    .attr('rx', 4)
    .attr('fill', theme.surface)
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1)

  g.append('text')
    .attr('x', x)
    .attr('y', y)
    .attr('font-family', theme.fontFamilyMono)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textAnnotation)
    .text(eq)
}

/**
 * Add an n= label below a group's x-position.
 */
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
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted)
    .text(`n = ${n}`)
}

/**
 * Add a stat annotation box (pill) directly on the plot area.
 * E.g. AUC = 0.82, r = .91, p < .001
 */
export function addStatBadge(
  g: GSelection,
  lines: string[],
  x: number,
  y: number,
  theme: CarmTheme = DEFAULT_THEME
): void {
  const lineH = 14
  const padX = 10, padY = 6
  const maxLen = Math.max(...lines.map(l => l.length))
  const bw = maxLen * 6.4 + padX * 2
  const bh = lines.length * lineH + padY * 2

  g.append('rect')
    .attr('x', x)
    .attr('y', y)
    .attr('width', bw)
    .attr('height', bh)
    .attr('rx', 5)
    .attr('fill', theme.surface)
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1)
    .attr('opacity', 0.92)

  lines.forEach((line, i) => {
    g.append('text')
      .attr('x', x + padX)
      .attr('y', y + padY + (i + 1) * lineH - 2)
      .attr('font-family', theme.fontFamilyMono)
      .attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.textAnnotation)
      .text(line)
  })
}
