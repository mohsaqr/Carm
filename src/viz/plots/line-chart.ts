/**
 * Multi-line chart with optional area fill.
 * Each series has its own line + color. Points at each data location.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface LineChartConfig {
  readonly showArea?: boolean
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface LineChartSeries {
  readonly label: string
  readonly x: readonly number[]
  readonly y: readonly number[]
}

export interface LineChartData {
  readonly series: readonly LineChartSeries[]
  readonly testResult?: StatResult
}

export function renderLineChart(
  container: HTMLElement,
  data: LineChartData,
  config: LineChartConfig = {}
): void {
  import('d3').then(d3 => renderLineChartD3(d3, container, data, config))
}

function renderLineChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: LineChartData,
  config: LineChartConfig
): void {
  if (data.series.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 420
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom
  const showArea = config.showArea ?? false

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Line Chart', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allX = data.series.flatMap(s => [...s.x])
  const allY = data.series.flatMap(s => [...s.y])
  const [xMin, xMax] = d3.extent(allX) as [number, number]
  const [yMin, yMax] = d3.extent(allY) as [number, number]
  const xPad = (xMax - xMin) * 0.05
  const yPad = (yMax - yMin) * 0.1

  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()
  const yScale = d3.scaleLinear().domain([Math.min(0, yMin - yPad), yMax + yPad]).range([height, 0]).nice()

  // Grid lines
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  data.series.forEach((s, si) => {
    if (s.x.length === 0 || s.y.length === 0) return
    const color = getColor(si, theme)
    const pts = s.x.map((x, i) => [x, s.y[i] ?? 0] as [number, number])
      .sort((a, b) => a[0] - b[0])

    // Area fill
    if (showArea) {
      const areaGen = d3.area<[number, number]>()
        .x(d => xScale(d[0])).y0(yScale(0)).y1(d => yScale(d[1])).curve(d3.curveCatmullRom)
      g.append('path').datum(pts).attr('d', areaGen)
        .attr('fill', color).attr('opacity', theme.ciOpacity)
    }

    // Line
    const lineGen = d3.line<[number, number]>()
      .x(d => xScale(d[0])).y(d => yScale(d[1])).curve(d3.curveCatmullRom)
    g.append('path').datum(pts).attr('d', lineGen)
      .attr('fill', 'none').attr('stroke', color).attr('stroke-width', 2)

    // Points
    pts.forEach(([x, y]) => {
      g.append('circle')
        .attr('cx', xScale(x)).attr('cy', yScale(y)).attr('r', 4)
        .attr('fill', color).attr('stroke', theme.background).attr('stroke-width', 1.5)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Series', s.label),
            formatTooltipRow('x', x.toFixed(3)),
            formatTooltipRow('y', y.toFixed(3)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'x')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'y')

  // Legend
  if (data.series.length > 1) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme)
      const lx = width - 130
      const ly = i * 18 + 8
      g.append('line').attr('x1', lx).attr('x2', lx + 20).attr('y1', ly).attr('y2', ly)
        .attr('stroke', color).attr('stroke-width', 2)
      g.append('circle').attr('cx', lx + 10).attr('cy', ly).attr('r', 3)
        .attr('fill', color).attr('stroke', theme.background).attr('stroke-width', 1)
      g.append('text').attr('x', lx + 24).attr('y', ly + 4)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text).text(s.label)
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
