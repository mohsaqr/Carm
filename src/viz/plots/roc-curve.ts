/**
 * ROC curve: TPR vs FPR with AUC annotation.
 * Includes diagonal reference line, filled area under curve.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ROCCurveConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface ROCCurveData {
  readonly fpr: readonly number[]
  readonly tpr: readonly number[]
  readonly auc: number
  readonly testResult?: StatResult
}

export function renderROCCurve(
  container: HTMLElement,
  data: ROCCurveData,
  config: ROCCurveConfig = {}
): void {
  import('d3').then(d3 => renderROCCurveD3(d3, container, data, config))
}

interface ROCPoint {
  readonly fprVal: number
  readonly tprVal: number
}

function renderROCCurveD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ROCCurveData,
  config: ROCCurveConfig
): void {
  if (data.fpr.length === 0 || data.tpr.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 520, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  const aucStr = `AUC = ${data.auc.toFixed(4)}`
  const subtitle = data.testResult?.formatted
    ? `${data.testResult.formatted}   ${aucStr}`
    : aucStr
  addSubtitle(svg, config.title ?? 'ROC Curve', subtitle, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width])
  const yScale = d3.scaleLinear().domain([0, 1]).range([height, 0])

  // Grid
  g.selectAll<SVGLineElement, number>('.grid-h').data(yScale.ticks(5)).join('line')
    .attr('class', 'grid-h').attr('x1', 0).attr('x2', width)
    .attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  // Diagonal reference line
  g.append('line')
    .attr('x1', xScale(0)).attr('y1', yScale(0))
    .attr('x2', xScale(1)).attr('y2', yScale(1))
    .attr('stroke', theme.axisLine)
    .attr('stroke-width', 1.5)
    .attr('stroke-dasharray', '6,4')

  // Sort by FPR and build typed pairs
  const pairs: ROCPoint[] = data.fpr
    .map((fprVal, i) => ({ fprVal, tprVal: data.tpr[i] ?? 0 }))
    .sort((a, b) => a.fprVal - b.fprVal)

  const curveColor = getColor(0, theme)

  // Area under ROC curve
  const areaFn = d3.area<ROCPoint>()
    .x(d => xScale(d.fprVal))
    .y0(yScale(0))
    .y1(d => yScale(d.tprVal))
    .curve(d3.curveLinear)

  g.append('path').datum(pairs)
    .attr('d', areaFn)
    .attr('fill', curveColor)
    .attr('opacity', theme.ciOpacity)

  // ROC line
  const lineFn = d3.line<ROCPoint>()
    .x(d => xScale(d.fprVal))
    .y(d => yScale(d.tprVal))
    .curve(d3.curveLinear)

  g.append('path').datum(pairs)
    .attr('d', lineFn)
    .attr('fill', 'none')
    .attr('stroke', curveColor)
    .attr('stroke-width', 2.5)

  // Interactive dots
  g.selectAll<SVGCircleElement, ROCPoint>('.roc-dot')
    .data(pairs).join('circle')
    .attr('class', 'roc-dot')
    .attr('cx', d => xScale(d.fprVal))
    .attr('cy', d => yScale(d.tprVal))
    .attr('r', 3)
    .attr('fill', curveColor)
    .attr('opacity', 0.6)
    .on('mouseover', function(event: MouseEvent, d: ROCPoint) {
      showTooltip(event, [
        formatTooltipRow('FPR', d.fprVal.toFixed(4)),
        formatTooltipRow('TPR', d.tprVal.toFixed(4)),
        formatTooltipRow('Specificity', (1 - d.fprVal).toFixed(4)),
        formatTooltipRow('Sensitivity', d.tprVal.toFixed(4))
      ].join(''), theme)
    })
    .on('mouseout', hideTooltip)

  // AUC annotation box
  g.append('rect')
    .attr('x', 8).attr('y', 8)
    .attr('width', 130).attr('height', 28)
    .attr('fill', theme.surface)
    .attr('stroke', theme.gridLine)
    .attr('rx', 4).attr('opacity', 0.9)

  g.append('text')
    .attr('x', 16).attr('y', 26)
    .attr('font-family', theme.fontFamilyMono ?? theme.fontFamily)
    .attr('font-size', theme.fontSize + 1)
    .attr('font-weight', '600')
    .attr('fill', curveColor)
    .text(`AUC = ${data.auc.toFixed(4)}`)

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.xLabel ?? 'False Positive Rate')

  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
    .attr('fill', theme.text).text(config.yLabel ?? 'True Positive Rate')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
