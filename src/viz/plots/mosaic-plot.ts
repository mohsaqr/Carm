/**
 * Mosaic plot — column widths proportional to column marginals, row heights
 * proportional to within-column conditional proportions.
 * Cells are coloured by Pearson residual: blue = positive, red = negative,
 * near-zero = neutral grey.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface MosaicPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface MosaicPlotData {
  /** table[i][j] = count for row i, column j */
  readonly table: readonly (readonly number[])[]
  readonly rowLabels: readonly string[]
  readonly colLabels: readonly string[]
  readonly testResult?: StatResult
}

/**
 * Render a mosaic plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - contingency table + row/col labels + optional stat result
 * @param config - visual configuration
 */
export function renderMosaicPlot(
  container: HTMLElement,
  data: MosaicPlotData,
  config: MosaicPlotConfig = {}
): void {
  import('d3').then(d3 => renderMosaicPlotD3(d3, container, data, config))
}

function renderMosaicPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: MosaicPlotData,
  config: MosaicPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight + 80,   // legend space
    bottom: theme.marginBottom,
    left: theme.marginLeft,
  }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Mosaic Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const nRows = data.table.length
  const nCols = data.table[0]?.length ?? 0
  if (nRows === 0 || nCols === 0) return

  // Marginals
  const colSums = Array.from({ length: nCols }, (_, j) =>
    data.table.reduce((s, row) => s + (row[j] ?? 0), 0)
  )
  const grandTotal = colSums.reduce((a, b) => a + b, 0)
  const rowSums = data.table.map(row => row.reduce((a, b) => a + b, 0))

  // Expected values for Pearson residuals
  const expected = data.table.map((row, i) =>
    row.map((_, j) => (rowSums[i]! * colSums[j]!) / grandTotal)
  )

  // Pearson residuals
  const residuals = data.table.map((row, i) =>
    row.map((v, j) => {
      const e = expected[i]![j]!
      return e > 0 ? (v - e) / Math.sqrt(e) : 0
    })
  )
  const maxResid = Math.max(1, ...residuals.flat().map(Math.abs))

  // Colour scale: diverging blue–grey–red
  const colorScale = (residual: number): string => {
    const t = residual / maxResid   // [-1, 1]
    if (t > 0) {
      // positive → blue
      const intensity = Math.min(t, 1)
      const r = Math.round(255 * (1 - intensity * 0.65))
      const gv = Math.round(255 * (1 - intensity * 0.5))
      const b = 255
      return `rgb(${r},${gv},${b})`
    } else if (t < 0) {
      // negative → red
      const intensity = Math.min(-t, 1)
      const r = 255
      const gv = Math.round(255 * (1 - intensity * 0.65))
      const bv = Math.round(255 * (1 - intensity * 0.65))
      return `rgb(${r},${gv},${bv})`
    }
    return theme.gridLine
  }

  // Column x-positions (cumulative proportional widths)
  const colWidths = colSums.map(s => (s / grandTotal) * width)
  const colX: number[] = []
  colWidths.reduce((acc, w, i) => { colX[i] = acc; return acc + w }, 0)

  const gap = 2

  // Draw cells
  data.table.forEach((row, ri) => {
    row.forEach((count, ci) => {
      const colSum = colSums[ci] ?? 1
      const cellH = colSum > 0 ? (count / colSum) * height : 0

      // Row y-position within column: cumulative from top
      const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[ci] ?? 0), 0)
      const cellY = colSum > 0 ? (rowsAbove / colSum) * height : 0

      const cx = colX[ci] ?? 0
      const cw = colWidths[ci] ?? 0
      const resid = residuals[ri]![ci]!
      const fillColor = colorScale(resid)

      g.append('rect')
        .attr('x', cx + gap / 2)
        .attr('y', cellY + gap / 2)
        .attr('width', Math.max(cw - gap, 0))
        .attr('height', Math.max(cellH - gap, 0))
        .attr('fill', fillColor)
        .attr('stroke', theme.background)
        .attr('stroke-width', 1)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Row', data.rowLabels[ri] ?? `Row ${ri}`),
            formatTooltipRow('Column', data.colLabels[ci] ?? `Col ${ci}`),
            formatTooltipRow('Count', count),
            formatTooltipRow('Pearson residual', resid.toFixed(3)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)

      // Cell label if large enough
      if (cw > 28 && cellH > 14) {
        g.append('text')
          .attr('x', cx + cw / 2)
          .attr('y', cellY + cellH / 2 + 4)
          .attr('text-anchor', 'middle')
          .attr('font-family', theme.fontFamily)
          .attr('font-size', theme.fontSizeSmall - 1)
          .attr('fill', Math.abs(resid) > maxResid * 0.5 ? '#fff' : theme.text)
          .text(String(count))
      }
    })
  })

  // Column labels (x-axis)
  colX.forEach((x, ci) => {
    const cw = colWidths[ci] ?? 0
    g.append('text')
      .attr('x', x + cw / 2)
      .attr('y', height + 18)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(data.colLabels[ci] ?? `Col ${ci}`)
  })

  // Row labels (y-axis) — positioned at midpoint of each row stripe in the first column
  data.table.forEach((row, ri) => {
    const colSum = colSums[0] ?? 1
    const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[0] ?? 0), 0)
    const cellY = colSum > 0 ? (rowsAbove / colSum) * height : 0
    const cellH = colSum > 0 ? ((row[0] ?? 0) / colSum) * height : 0
    g.append('text')
      .attr('x', -6)
      .attr('y', cellY + cellH / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(data.rowLabels[ri] ?? `Row ${ri}`)
  })

  // Axis labels
  g.append('text')
    .attr('x', width / 2).attr('y', height + 40)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.xLabel ?? '')

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2).attr('y', -50)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSize)
    .attr('fill', theme.text)
    .text(config.yLabel ?? '')

  // Pearson residual legend bar (right side)
  const lgX = width + 16
  const lgH = 120
  const lgW = 14
  const lgY = height / 2 - lgH / 2
  const steps = 20
  Array.from({ length: steps }, (_, i) => {
    const t = 1 - (i / (steps - 1)) * 2   // 1 down to -1
    g.append('rect')
      .attr('x', lgX)
      .attr('y', lgY + (i / steps) * lgH)
      .attr('width', lgW)
      .attr('height', lgH / steps + 0.5)
      .attr('fill', colorScale(t * maxResid))
  })
  g.append('text').attr('x', lgX + lgW + 4).attr('y', lgY + 4)
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted).text(`+${maxResid.toFixed(1)}`)
  g.append('text').attr('x', lgX + lgW + 4).attr('y', lgY + lgH / 2 + 4)
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted).text('0')
  g.append('text').attr('x', lgX + lgW + 4).attr('y', lgY + lgH + 4)
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted).text(`-${maxResid.toFixed(1)}`)
  g.append('text').attr('x', lgX + lgW / 2).attr('y', lgY - 8)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall - 1)
    .attr('fill', theme.textMuted).text('Residual')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
