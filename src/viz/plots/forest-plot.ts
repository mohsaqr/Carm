/**
 * Forest plot: per-study CI lines with square points + pooled diamond.
 * Vertical reference line at 0. Numeric estimates shown on the right.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ForestPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface ForestStudy {
  readonly label: string
  readonly estimate: number
  readonly ciLow: number
  readonly ciHigh: number
  readonly weight?: number
}

export interface ForestPooled {
  readonly estimate: number
  readonly ciLow: number
  readonly ciHigh: number
}

export interface ForestPlotData {
  readonly studies: readonly ForestStudy[]
  readonly pooled?: ForestPooled
  readonly testResult?: StatResult
}

export function renderForestPlot(
  container: HTMLElement,
  data: ForestPlotData,
  config: ForestPlotConfig = {}
): void {
  import('d3').then(d3 => renderForestPlotD3(d3, container, data, config))
}

function renderForestPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ForestPlotData,
  config: ForestPlotConfig
): void {
  if (data.studies.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const totalRows = data.studies.length + (data.pooled != null ? 2 : 0)
  const rowH = 28
  const extraH = theme.marginTop + theme.marginBottom + 20
  const H = config.height ?? Math.max(totalRows * rowH + extraH, 320)
  const W = config.width ?? Math.max(container.clientWidth || 700, 500)

  const labelColW = 160
  const estColW = 140
  const marginLeft = labelColW + 8
  const marginRight = estColW + 16
  const margin = { top: theme.marginTop, right: marginRight, bottom: theme.marginBottom, left: marginLeft }
  const width = W - margin.left - margin.right
  const height = totalRows * rowH

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H + margin.top + margin.bottom).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Forest Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const allLow = data.studies.map(s => s.ciLow).concat(data.pooled != null ? [data.pooled.ciLow] : [])
  const allHigh = data.studies.map(s => s.ciHigh).concat(data.pooled != null ? [data.pooled.ciHigh] : [])
  const xMin = Math.min(...allLow)
  const xMax = Math.max(...allHigh)
  const xPad = (xMax - xMin) * 0.12
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice()

  const maxWeight = data.studies.reduce((m, s) => Math.max(m, s.weight ?? 1), 0)

  // Reference line at 0
  const [domainMin, domainMax] = xScale.domain() as [number, number]
  if (domainMin <= 0 && domainMax >= 0) {
    g.append('line')
      .attr('x1', xScale(0)).attr('x2', xScale(0))
      .attr('y1', 0).attr('y2', height)
      .attr('stroke', theme.axisLine)
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '4,3')
  }

  // Column headers
  svg.append('text')
    .attr('x', margin.left - 8).attr('y', margin.top - 6)
    .attr('text-anchor', 'end')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
    .attr('font-weight', '600').attr('fill', theme.textMuted)
    .text('Study')
  svg.append('text')
    .attr('x', W - marginRight + 8).attr('y', margin.top - 6)
    .attr('text-anchor', 'start')
    .attr('font-family', theme.fontFamilyMono ?? theme.fontFamily).attr('font-size', theme.fontSizeSmall)
    .attr('font-weight', '600').attr('fill', theme.textMuted)
    .text('Estimate [95% CI]')

  const studyColor = getColor(0, theme)

  data.studies.forEach((study, i) => {
    const cy = i * rowH + rowH / 2

    // CI line
    g.append('line')
      .attr('x1', xScale(study.ciLow)).attr('x2', xScale(study.ciHigh))
      .attr('y1', cy).attr('y2', cy)
      .attr('stroke', studyColor).attr('stroke-width', 1.5)

    // CI caps
    const capH = 5
    ;[study.ciLow, study.ciHigh].forEach(xv => {
      g.append('line')
        .attr('x1', xScale(xv)).attr('x2', xScale(xv))
        .attr('y1', cy - capH).attr('y2', cy + capH)
        .attr('stroke', studyColor).attr('stroke-width', 1.5)
    })

    // Square point
    const sqHalf = study.weight != null && maxWeight > 0
      ? 4 + (study.weight / maxWeight) * 6
      : 6
    g.append('rect')
      .attr('x', xScale(study.estimate) - sqHalf / 2)
      .attr('y', cy - sqHalf / 2)
      .attr('width', sqHalf).attr('height', sqHalf)
      .attr('fill', studyColor)
      .on('mouseover', function(event: MouseEvent) {
        showTooltip(event, [
          formatTooltipRow('Study', study.label),
          formatTooltipRow('Estimate', study.estimate.toFixed(3)),
          formatTooltipRow('95% CI', `[${study.ciLow.toFixed(3)}, ${study.ciHigh.toFixed(3)}]`),
          ...(study.weight != null ? [formatTooltipRow('Weight', study.weight.toFixed(2))] : [])
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Label
    svg.append('text')
      .attr('x', margin.left - 8).attr('y', margin.top + cy)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
      .attr('fill', theme.text).text(study.label)

    // Numeric estimate
    const estStr = `${study.estimate.toFixed(2)} [${study.ciLow.toFixed(2)}, ${study.ciHigh.toFixed(2)}]`
    svg.append('text')
      .attr('x', W - marginRight + 8).attr('y', margin.top + cy)
      .attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamilyMono ?? theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(estStr)
  })

  // Pooled diamond
  if (data.pooled != null) {
    const pooled = data.pooled
    const sepY = data.studies.length * rowH
    g.append('line')
      .attr('x1', 0).attr('x2', width)
      .attr('y1', sepY).attr('y2', sepY)
      .attr('stroke', theme.axisLine).attr('stroke-width', 1)

    const pooledY = sepY + rowH * 1.0
    const diamondColor = getColor(4, theme)
    const dLeft = xScale(pooled.ciLow)
    const dRight = xScale(pooled.ciHigh)
    const dMid = xScale(pooled.estimate)
    const dH = 9

    const diamondPath = [
      `M ${dMid},${pooledY - dH}`,
      `L ${dRight},${pooledY}`,
      `L ${dMid},${pooledY + dH}`,
      `L ${dLeft},${pooledY}`,
      'Z'
    ].join(' ')

    g.append('path').attr('d', diamondPath)
      .attr('fill', diamondColor).attr('opacity', 0.88)
      .attr('stroke', diamondColor).attr('stroke-width', 1)
      .on('mouseover', function(event: MouseEvent) {
        showTooltip(event, [
          formatTooltipRow('Pooled estimate', pooled.estimate.toFixed(3)),
          formatTooltipRow('95% CI', `[${pooled.ciLow.toFixed(3)}, ${pooled.ciHigh.toFixed(3)}]`)
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    svg.append('text')
      .attr('x', margin.left - 8).attr('y', margin.top + pooledY)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
      .attr('font-weight', '600').attr('fill', theme.text).text('Pooled')

    const pEstStr = `${pooled.estimate.toFixed(2)} [${pooled.ciLow.toFixed(2)}, ${pooled.ciHigh.toFixed(2)}]`
    svg.append('text')
      .attr('x', W - marginRight + 8).attr('y', margin.top + pooledY)
      .attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamilyMono ?? theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('font-weight', '600').attr('fill', theme.text).text(pEstStr)
  }

  // X axis
  g.append('g').attr('transform', `translate(0,${height + 8})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'Effect Size')

  if (config.caption) addCaption(svg, config.caption, W, H + margin.top + margin.bottom, theme)
}
