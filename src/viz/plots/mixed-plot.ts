/**
 * LMM visualization: caterpillar plot of BLUPs + variance component chart.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import type { LMMResult } from '../../core/types.js'

export interface MixedPlotConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export function renderMixedPlot(
  container: HTMLElement,
  result: LMMResult,
  blups: ReadonlyArray<{ group: string | number; blup: number }>,
  config: MixedPlotConfig = {}
): void {
  import('d3').then(d3 => renderMixedD3(d3, container, result, blups, config))
}

function renderMixedD3(
  d3: typeof D3,
  container: HTMLElement,
  result: LMMResult,
  blups: ReadonlyArray<{ group: string | number; blup: number }>,
  config: MixedPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 600
  const H = config.height ?? Math.max(blups.length * 20 + 150, 300)
  const margin = { top: theme.marginTop, right: 80, bottom: theme.marginBottom, left: 100 }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Random Effects (BLUPs)', result.formatted, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const sorted = [...blups].sort((a, b) => a.blup - b.blup)
  const groupLabels = sorted.map(b => String(b.group))
  const blupValues = sorted.map(b => b.blup)

  const xExt = d3.extent(blupValues) as [number, number]
  const xPad = Math.max(Math.abs(xExt[0]), Math.abs(xExt[1])) * 0.2
  const xScale = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, width]).nice()
  const yScale = d3.scaleBand<string>().domain(groupLabels).range([height, 0]).padding(0.3)

  // Zero line
  g.append('line').attr('x1', xScale(0)).attr('x2', xScale(0))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.axisLine).attr('stroke-dasharray', '4,2').attr('stroke-width', 1.5)

  // Grid
  g.selectAll('.grid').data(xScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', d => xScale(d)).attr('x2', d => xScale(d))
    .attr('y1', 0).attr('y2', height)
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // BLUP dots
  sorted.forEach((b, _i) => {
    const cy = (yScale(String(b.group)) ?? 0) + yScale.bandwidth() / 2
    const color = b.blup >= 0 ? getColor(0, theme) : getColor(5, theme)
    g.append('circle')
      .attr('cx', xScale(b.blup)).attr('cy', cy)
      .attr('r', 4).attr('fill', color)
    // Whisker from 0
    g.append('line')
      .attr('x1', xScale(0)).attr('x2', xScale(b.blup))
      .attr('y1', cy).attr('y2', cy)
      .attr('stroke', color).attr('stroke-width', 1.5).attr('opacity', 0.5)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text('BLUP (random intercept)')

  // ICC annotation
  svg.append('text').attr('x', W - 10).attr('y', H - 10)
    .attr('text-anchor', 'end').attr('font-family', theme.fontFamilyMono).attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textMuted).text(`ICC = ${result.icc.toFixed(3)}`)

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
