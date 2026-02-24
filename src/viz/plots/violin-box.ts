/**
 * Violin + Box + Jitter combo plot (ggbetweenstats style).
 * Shows distribution shape (violin), summary (box), individual points (jitter),
 * and significance brackets between groups.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption, addNLabel } from '../components/annotations.js'
import { renderBrackets } from '../components/brackets.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult, PairwiseResult, GroupData } from '../../core/types.js'
import { mean as _mean, sd as _sd, quantile } from '../../core/math.js'

export interface ViolinBoxConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showJitter?: boolean
  readonly showBrackets?: boolean
  readonly significantBracketsOnly?: boolean
  readonly jitterWidth?: number
  readonly violinBandwidth?: number
  readonly showN?: boolean
  readonly showMean?: boolean
  readonly showMedian?: boolean
  readonly numericP?: boolean
}

export interface ViolinBoxData {
  readonly groups: readonly GroupData[]
  readonly testResult?: StatResult
  readonly pairwise?: readonly PairwiseResult[]
}

/**
 * Render a violin + box + jitter plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - group data + optional stat results
 * @param config - visual configuration
 */
export function renderViolinBox(
  container: HTMLElement,
  data: ViolinBoxData,
  config: ViolinBoxConfig = {}
): void {
  // Load d3 dynamically to work as peer dependency
  import('d3').then(d3 => renderViolinBoxD3(d3, container, data, config))
}

function renderViolinBoxD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ViolinBoxData,
  config: ViolinBoxConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = {
    top: theme.marginTop + (data.pairwise?.length ?? 0) * 20,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft,
  }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  // Clear container
  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  // Title + subtitle
  if (config.title || data.testResult) {
    addSubtitle(
      svg,
      config.title ?? 'Group Comparison',
      data.testResult?.formatted ?? '',
      W,
      theme
    )
  }

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`)

  // Scales
  const labels = data.groups.map(gr => gr.label)
  const xScale = d3.scaleBand<string>()
    .domain(labels)
    .range([0, width])
    .padding(0.2)

  const allValues = data.groups.flatMap(gr => [...gr.values])
  const [yMin, yMax] = d3.extent(allValues) as [number, number]
  const yPad = (yMax - yMin) * 0.1
  const yScale = d3.scaleLinear()
    .domain([yMin - yPad, yMax + yPad])
    .range([height, 0])
    .nice()

  // Grid lines
  g.selectAll('.grid')
    .data(yScale.ticks(6))
    .join('line')
    .attr('class', 'grid')
    .attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1)

  // Axes
  g.append('g')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(xScale))
    .selectAll('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  g.append('g')
    .call(d3.axisLeft(yScale).ticks(6))
    .selectAll('text')
    .attr('fill', theme.text)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)

  // Axis labels
  g.append('text')
    .attr('x', width / 2)
    .attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.xLabel ?? '')

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.yLabel ?? 'Value')

  // Draw each group: violin → box → jitter
  data.groups.forEach((gr, gi) => {
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2
    const color = getColor(gi, theme)
    const bw = config.violinBandwidth ?? (yMax - yMin) / 20

    // Violin (KDE approximation)
    const kdePoints = kernelDensityEstimate(gr.values, yScale.ticks(60), bw)
    const maxDensity = Math.max(...kdePoints.map(p => p[1]))
    const violinWidth = xScale.bandwidth() * 0.45
    const violinScale = maxDensity > 0 ? violinWidth / maxDensity : 1

    const areaFn = d3.area<[number, number]>()
      .x0(d => cx - d[1] * violinScale)
      .x1(d => cx + d[1] * violinScale)
      .y(d => yScale(d[0]))
      .curve(d3.curveCatmullRom)

    g.append('path')
      .datum(kdePoints)
      .attr('d', areaFn)
      .attr('fill', color)
      .attr('opacity', theme.violinOpacity)
      .attr('stroke', color)
      .attr('stroke-width', 1)

    // Box plot stats
    const q1 = quantile(gr.values, 0.25)
    const med = quantile(gr.values, 0.5)
    const q3 = quantile(gr.values, 0.75)
    const iqr = q3 - q1
    const whiskerLo = Math.min(...gr.values.filter(v => v >= q1 - 1.5 * iqr))
    const whiskerHi = Math.max(...gr.values.filter(v => v <= q3 + 1.5 * iqr))
    const boxW = xScale.bandwidth() * 0.2

    // Whiskers
    g.append('line').attr('x1', cx).attr('x2', cx)
      .attr('y1', yScale(whiskerLo)).attr('y2', yScale(q1))
      .attr('stroke', color).attr('stroke-width', 1.5)
    g.append('line').attr('x1', cx).attr('x2', cx)
      .attr('y1', yScale(q3)).attr('y2', yScale(whiskerHi))
      .attr('stroke', color).attr('stroke-width', 1.5)

    // Box
    g.append('rect')
      .attr('x', cx - boxW / 2).attr('width', boxW)
      .attr('y', yScale(q3)).attr('height', yScale(q1) - yScale(q3))
      .attr('fill', theme.background)
      .attr('stroke', color)
      .attr('stroke-width', 2)

    // Median line + value
    if (config.showMedian !== false) {
      g.append('line').attr('x1', cx - boxW / 2).attr('x2', cx + boxW / 2)
        .attr('y1', yScale(med)).attr('y2', yScale(med))
        .attr('stroke', color).attr('stroke-width', 2.5)
      g.append('text')
        .attr('x', cx + boxW / 2 + 4).attr('y', yScale(med) + 3.5)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textAnnotation).text(med.toFixed(2))
    }

    // Mean diamond marker + value
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length
      const my = yScale(groupMean)
      const ds = 5
      g.append('polygon')
        .attr('points', `${cx},${my - ds} ${cx + ds},${my} ${cx},${my + ds} ${cx - ds},${my}`)
        .attr('fill', 'white').attr('stroke', color).attr('stroke-width', 1.5)
      g.append('text')
        .attr('x', cx + ds + 4).attr('y', my + 3.5)
        .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
        .attr('fill', theme.textAnnotation).text(groupMean.toFixed(2))
    }

    // Jitter points
    if (config.showJitter !== false) {
      const jw = (config.jitterWidth ?? 0.15) * xScale.bandwidth()
      const seed = gi * 12345
      gr.values.forEach((v, vi) => {
        const jx = cx + (pseudoRandom(seed + vi) - 0.5) * jw
        g.append('circle')
          .attr('cx', jx).attr('cy', yScale(v))
          .attr('r', 3)
          .attr('fill', color)
          .attr('opacity', theme.pointOpacity)
          .attr('stroke', 'none')
          .on('mouseover', (event: MouseEvent) => {
            showTooltip(event, [
              formatTooltipRow('Group', gr.label),
              formatTooltipRow('Value', v.toFixed(3)),
            ].join(''), theme)
          })
          .on('mouseout', hideTooltip)
      })
    }

    // n label
    if (config.showN !== false) {
      addNLabel(g, gr.values.length, cx, height + 60, theme)
    }
  })

  // Significance brackets
  if (config.showBrackets !== false && data.pairwise && data.pairwise.length > 0) {
    const posMap = new Map(
      data.groups.map(gr => [gr.label, (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2])
    )
    renderBrackets(g, data.pairwise, {
      groupPositions: posMap,
      yBase: 0,
      bracketHeight: 22,
      significantOnly: config.significantBracketsOnly ?? false,
      ...(config.numericP !== undefined && { numericP: config.numericP }),
    }, theme)
  }

  // Caption
  if (config.caption) {
    addCaption(svg, config.caption, W, H, theme)
  }
}

// ─── KDE ─────────────────────────────────────────────────────────────────

function kernelDensityEstimate(data: readonly number[], xPoints: number[], bandwidth: number): Array<[number, number]> {
  return xPoints.map(x => {
    const density = data.reduce((sum, xi) => sum + gaussianKernel((x - xi) / bandwidth), 0) / (data.length * bandwidth)
    return [x, density] as [number, number]
  })
}

function gaussianKernel(u: number): number {
  return Math.exp(-0.5 * u * u) / Math.sqrt(2 * Math.PI)
}

function pseudoRandom(seed: number): number {
  const x = Math.sin(seed) * 10000
  return x - Math.floor(x)
}
