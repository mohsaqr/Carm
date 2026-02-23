/**
 * Area chart with optional stacking.
 * Non-stacked: overlapping semi-transparent areas + lines per series.
 * Stacked: cumulative areas per series.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface AreaChartConfig {
  readonly stacked?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface AreaSeries {
  readonly label: string
  readonly x: readonly number[]
  readonly y: readonly number[]
}

export interface AreaChartData {
  readonly series: readonly AreaSeries[]
  readonly testResult?: StatResult
}

export function renderAreaChart(
  container: HTMLElement,
  data: AreaChartData,
  config: AreaChartConfig = {}
): void {
  import('d3').then(d3 => renderAreaChartD3(d3, container, data, config))
}

interface XYPoint {
  readonly xVal: number
  readonly yVal: number
}

interface StackedPoint {
  readonly xVal: number
  readonly top: number
  readonly bot: number
}

function renderAreaChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: AreaChartData,
  config: AreaChartConfig
): void {
  if (data.series.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Area Chart', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allX = data.series.flatMap(s => [...s.x])
  const xMin = Math.min(...allX)
  const xMax = Math.max(...allX)
  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width])

  const stacked = config.stacked ?? false

  let yMax: number
  if (stacked) {
    const allXSorted = [...new Set(allX)].sort((a, b) => a - b)
    const rowSums = allXSorted.map(xVal =>
      data.series.reduce((sum, s) => {
        const idx = s.x.indexOf(xVal)
        return sum + (idx >= 0 ? (s.y[idx] ?? 0) : 0)
      }, 0)
    )
    yMax = Math.max(...rowSums)
  } else {
    yMax = Math.max(...data.series.flatMap(s => [...s.y]))
  }

  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([height, 0]).nice()

  // Grid
  g.selectAll<SVGLineElement, number>('.grid').data(yScale.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  if (stacked) {
    const allXSorted = [...new Set(allX)].sort((a, b) => a - b)
    const baselines = new Map<number, number>()
    allXSorted.forEach(xVal => baselines.set(xVal, 0))

    data.series.forEach((series, si) => {
      const color = getColor(si, theme)

      const points: StackedPoint[] = series.x.map((xVal, xi) => ({
        xVal,
        top: (baselines.get(xVal) ?? 0) + (series.y[xi] ?? 0),
        bot: baselines.get(xVal) ?? 0
      }))

      const areaFn = d3.area<StackedPoint>()
        .x(d => xScale(d.xVal))
        .y0(d => yScale(d.bot))
        .y1(d => yScale(d.top))
        .curve(d3.curveMonotoneX)

      g.append('path').datum(points)
        .attr('d', areaFn)
        .attr('fill', color)
        .attr('opacity', 0.75)
        .attr('stroke', color)
        .attr('stroke-width', 1.5)

      series.x.forEach((xVal, xi) => {
        baselines.set(xVal, (baselines.get(xVal) ?? 0) + (series.y[xi] ?? 0))
      })
    })
  } else {
    data.series.forEach((series, si) => {
      const color = getColor(si, theme)

      // Build typed points array with guaranteed number (filter out bad indices)
      const points: XYPoint[] = series.x
        .map((xVal, xi) => ({ xVal, yVal: series.y[xi] }))
        .filter((p): p is XYPoint => p.yVal !== undefined)

      const areaFn = d3.area<XYPoint>()
        .x(d => xScale(d.xVal))
        .y0(yScale(0))
        .y1(d => yScale(d.yVal))
        .curve(d3.curveMonotoneX)

      const lineFn = d3.line<XYPoint>()
        .x(d => xScale(d.xVal))
        .y(d => yScale(d.yVal))
        .curve(d3.curveMonotoneX)

      g.append('path').datum(points)
        .attr('d', areaFn)
        .attr('fill', color)
        .attr('opacity', 0.22)

      g.append('path').datum(points)
        .attr('d', lineFn)
        .attr('fill', 'none')
        .attr('stroke', color)
        .attr('stroke-width', 2.5)

      g.selectAll<SVGCircleElement, XYPoint>(`.dot-s${si}`)
        .data(points).join('circle')
        .attr('class', `dot-s${si}`)
        .attr('cx', d => xScale(d.xVal))
        .attr('cy', d => yScale(d.yVal))
        .attr('r', 3.5)
        .attr('fill', color)
        .attr('opacity', 0.7)
        .on('mouseover', function(event: MouseEvent, d: XYPoint) {
          showTooltip(event, [
            formatTooltipRow('Series', series.label),
            formatTooltipRow('x', d.xVal.toFixed(3)),
            formatTooltipRow('y', d.yVal.toFixed(3))
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(8))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'x')

  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Value')

  // Legend
  data.series.forEach((s, i) => {
    const lx = margin.left + i * 130
    const ly = H - 18
    svg.append('rect').attr('x', lx).attr('y', ly - 9).attr('width', 12).attr('height', 12)
      .attr('fill', getColor(i, theme)).attr('rx', 2)
    svg.append('text').attr('x', lx + 16).attr('y', ly)
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(s.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
