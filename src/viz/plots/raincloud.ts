/**
 * Raincloud plot: half-violin + jitter + box.
 * Premium visualization for showing full distribution shape.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption, addNLabel } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import { quantile } from '../../core/math.js'
import type { GroupData, StatResult } from '../../core/types.js'

export interface RaincloudConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface RaincloudData {
  readonly groups: readonly GroupData[]
  readonly testResult?: StatResult
}

export function renderRaincloud(
  container: HTMLElement,
  data: RaincloudData,
  config: RaincloudConfig = {}
): void {
  import('d3').then(d3 => renderRaincloudD3(d3, container, data, config))
}

function renderRaincloudD3(
  d3: typeof D3,
  container: HTMLElement,
  data: RaincloudData,
  config: RaincloudConfig
): void {
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

  addSubtitle(svg, config.title ?? 'Raincloud Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const labels = data.groups.map(gr => gr.label)
  const allValues = data.groups.flatMap(gr => [...gr.values])
  const [yMin, yMax] = d3.extent(allValues) as [number, number]
  const yPad = (yMax - yMin) * 0.1

  const xScale = d3.scaleBand<string>().domain(labels).range([0, width]).padding(0.3)
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice()

  // Grid
  g.selectAll('.grid').data(yScale.ticks(6)).join('line')
    .attr('class', 'grid').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  data.groups.forEach((gr, gi) => {
    const color = getColor(gi, theme)
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2
    const bw = (yMax - yMin) / 20
    const violinR = xScale.bandwidth() * 0.35

    // KDE
    const kdePoints = kdeEstimate(gr.values, yScale.ticks(60), bw)
    const maxD = Math.max(...kdePoints.map(p => p[1]))
    const dScale = maxD > 0 ? violinR / maxD : 1

    // Half violin (right side only — this is the "cloud")
    const areaFn = d3.area<[number, number]>()
      .x0(cx).x1(d => cx + d[1] * dScale)
      .y(d => yScale(d[0])).curve(d3.curveCatmullRom)

    g.append('path').datum(kdePoints).attr('d', areaFn)
      .attr('fill', color).attr('opacity', theme.violinOpacity)
      .attr('stroke', color).attr('stroke-width', 1)

    // Box (narrow, slightly left of center)
    const q1 = quantile(gr.values, 0.25)
    const med = quantile(gr.values, 0.5)
    const q3 = quantile(gr.values, 0.75)
    const iqr = q3 - q1
    const wLo = Math.min(...gr.values.filter(v => v >= q1 - 1.5 * iqr))
    const wHi = Math.max(...gr.values.filter(v => v <= q3 + 1.5 * iqr))
    const boxW = xScale.bandwidth() * 0.12
    const bx = cx - violinR * 0.5  // shift box left

    g.append('line').attr('x1', bx).attr('x2', bx)
      .attr('y1', yScale(wLo)).attr('y2', yScale(wHi))
      .attr('stroke', color).attr('stroke-width', 1.5)
    g.append('rect')
      .attr('x', bx - boxW / 2).attr('width', boxW)
      .attr('y', yScale(q3)).attr('height', yScale(q1) - yScale(q3))
      .attr('fill', theme.background).attr('stroke', color).attr('stroke-width', 2)
    g.append('line').attr('x1', bx - boxW / 2).attr('x2', bx + boxW / 2)
      .attr('y1', yScale(med)).attr('y2', yScale(med))
      .attr('stroke', color).attr('stroke-width', 2.5)

    // Jitter (the "rain" — left of box)
    const jx = bx - boxW - 4
    const seed = gi * 99991
    gr.values.forEach((v, vi) => {
      const jitter = (pseudoRnd(seed + vi) - 0.5) * xScale.bandwidth() * 0.12
      g.append('circle')
        .attr('cx', jx + jitter).attr('cy', yScale(v))
        .attr('r', 2.5).attr('fill', color).attr('opacity', theme.pointOpacity)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [formatTooltipRow('Group', gr.label), formatTooltipRow('Value', v.toFixed(3))].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    })

    addNLabel(g, gr.values.length, cx, height + 60, theme)
  })

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? '')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle').attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'Value')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

function kdeEstimate(data: readonly number[], xPoints: number[], bw: number): Array<[number, number]> {
  return xPoints.map(x => {
    const d = data.reduce((s, xi) => s + Math.exp(-0.5 * ((x - xi) / bw) ** 2) / (Math.sqrt(2 * Math.PI) * bw), 0) / data.length
    return [x, d] as [number, number]
  })
}

function pseudoRnd(seed: number): number {
  const x = Math.sin(seed) * 10000
  return x - Math.floor(x)
}
