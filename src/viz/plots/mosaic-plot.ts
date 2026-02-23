/**
 * Mosaic plot — column widths proportional to column marginals, row heights
 * proportional to within-column conditional proportions.
 *
 * Cells coloured by discretised Pearson residual (RdBu ColorBrewer bins):
 *   > +4   dark blue     #2166ac
 *  +2..+4  medium blue   #74add1
 *   0..+2  light blue    #d1e5f0
 *  -2.. 0  light peach   #fddbc7
 *  -4..-2  salmon        #d6604d
 *   < -4   dark red      #b2182b
 *
 * Borders:
 *   |residual| ≥ 4  → thick solid (2.5 px)
 *   |residual| ≥ 2  → dashed (stroke-dasharray 4 2)
 *   otherwise       → thin solid (0.5 px, background colour)
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

// ─── Discrete RdBu bins (matches R vcd / ggmosaic convention) ─────────────

interface ResidBin {
  readonly lo: number       // inclusive lower bound
  readonly hi: number       // exclusive upper bound (Infinity for last)
  readonly fill: string
  readonly label: string
  readonly borderStyle: 'none' | 'dashed' | 'thick'
}

const RESID_BINS: readonly ResidBin[] = [
  { lo:  4,    hi:  Infinity, fill: '#2166ac', label: '> 4',   borderStyle: 'thick'  },
  { lo:  2,    hi:  4,        fill: '#74add1', label: '2 : 4', borderStyle: 'dashed' },
  { lo:  0,    hi:  2,        fill: '#d1e5f0', label: '0 : 2', borderStyle: 'none'   },
  { lo: -2,    hi:  0,        fill: '#fddbc7', label: '−2 : 0',borderStyle: 'none'   },
  { lo: -4,    hi: -2,        fill: '#d6604d', label: '−4 : −2',borderStyle: 'dashed'},
  { lo: -Infinity, hi: -4,   fill: '#b2182b', label: '< −4',  borderStyle: 'thick'  },
]

function binForResidual(r: number): ResidBin {
  return RESID_BINS.find(b => r >= b.lo && r < b.hi) ?? RESID_BINS[2]!
}

/** White text on dark bins, dark text on light bins. */
function textColor(fill: string): string {
  const dark = ['#2166ac', '#b2182b', '#d6604d']
  return dark.includes(fill) ? '#ffffff' : '#212529'
}

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
  const legendW = 110
  const W = config.width ?? Math.max(container.clientWidth || 620, 420)
  const H = config.height ?? 480
  const margin = {
    top:    theme.marginTop,
    right:  theme.marginRight + legendW,
    bottom: theme.marginBottom,
    left:   theme.marginLeft,
  }
  const width  = W - margin.left - margin.right
  const height = H - margin.top  - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Mosaic Plot', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const nRows = data.table.length
  const nCols = data.table[0]?.length ?? 0
  if (nRows === 0 || nCols === 0) return

  // ── Marginals ─────────────────────────────────────────────────────────────
  const colSums = Array.from({ length: nCols }, (_, j) =>
    data.table.reduce((s, row) => s + (row[j] ?? 0), 0)
  )
  const grandTotal = colSums.reduce((a, b) => a + b, 0)
  const rowSums    = data.table.map(row => row.reduce((a, b) => a + b, 0))

  // ── Pearson residuals ─────────────────────────────────────────────────────
  const residuals = data.table.map((row, i) =>
    row.map((v, j) => {
      const e = (rowSums[i]! * colSums[j]!) / grandTotal
      return e > 0 ? (v - e) / Math.sqrt(e) : 0
    })
  )

  // ── Column x-positions (cumulative proportional widths) ───────────────────
  const gap = 3
  const colWidths = colSums.map(s => (s / grandTotal) * width)
  const colX: number[] = []
  colWidths.reduce((acc, w, i) => { colX[i] = acc; return acc + w }, 0)

  // ── Draw cells ────────────────────────────────────────────────────────────
  data.table.forEach((row, ri) => {
    row.forEach((count, ci) => {
      const colSum = colSums[ci] ?? 1
      const cellH  = colSum > 0 ? (count / colSum) * height : 0
      const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[ci] ?? 0), 0)
      const cellY  = colSum > 0 ? (rowsAbove / colSum) * height : 0
      const cx     = colX[ci] ?? 0
      const cw     = colWidths[ci] ?? 0
      const resid  = residuals[ri]![ci]!
      const bin    = binForResidual(resid)

      const rx = cx + gap / 2
      const ry = cellY + gap / 2
      const rw = Math.max(cw - gap, 0)
      const rh = Math.max(cellH - gap, 0)

      // Background rect
      g.append('rect')
        .attr('x', rx).attr('y', ry)
        .attr('width', rw).attr('height', rh)
        .attr('fill', bin.fill)
        .attr('stroke', 'none')

      // Border overlay (significance indicator)
      if (bin.borderStyle !== 'none') {
        const isThick = bin.borderStyle === 'thick'
        g.append('rect')
          .attr('x', rx).attr('y', ry)
          .attr('width', rw).attr('height', rh)
          .attr('fill', 'none')
          .attr('stroke', isThick ? (resid > 0 ? '#08519c' : '#a50f15') : (resid > 0 ? '#4292c6' : '#cb181d'))
          .attr('stroke-width', isThick ? 2.5 : 1.5)
          .attr('stroke-dasharray', isThick ? 'none' : '5,2')

      }

      // Invisible hit-area for tooltip (covers full cell including gaps)
      g.append('rect')
        .attr('x', rx).attr('y', ry)
        .attr('width', rw).attr('height', rh)
        .attr('fill', 'transparent')
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Row',    data.rowLabels[ri] ?? `Row ${ri}`),
            formatTooltipRow('Column', data.colLabels[ci] ?? `Col ${ci}`),
            formatTooltipRow('Count',  count),
            formatTooltipRow('Residual', resid.toFixed(2)),
            formatTooltipRow('Bin',    bin.label),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)

      // Count label if cell is large enough
      if (rw > 30 && rh > 16) {
        g.append('text')
          .attr('x', rx + rw / 2)
          .attr('y', ry + rh / 2 + 4)
          .attr('text-anchor', 'middle')
          .attr('font-family', theme.fontFamily)
          .attr('font-size', Math.min(theme.fontSizeSmall, rh * 0.38, rw * 0.22))
          .attr('fill', textColor(bin.fill))
          .attr('pointer-events', 'none')
          .text(String(count))
      }
    })
  })

  // ── Column labels (x-axis) ────────────────────────────────────────────────
  colX.forEach((x, ci) => {
    const cw = colWidths[ci] ?? 0
    g.append('text')
      .attr('x', x + cw / 2).attr('y', height + 18)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(data.colLabels[ci] ?? `Col ${ci}`)
  })

  // ── Row labels (y-axis) — midpoint of each row in the first column ─────────
  data.table.forEach((row, ri) => {
    const colSum = colSums[0] ?? 1
    const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[0] ?? 0), 0)
    const cellY = colSum > 0 ? (rowsAbove / colSum) * height : 0
    const cellH = colSum > 0 ? ((row[0] ?? 0) / colSum) * height : 0
    g.append('text')
      .attr('x', -8).attr('y', cellY + cellH / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(data.rowLabels[ri] ?? `Row ${ri}`)
  })

  // ── Axis labels ───────────────────────────────────────────────────────────
  if (config.xLabel) {
    g.append('text')
      .attr('x', width / 2).attr('y', height + 44)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.xLabel)
  }
  if (config.yLabel) {
    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -height / 2).attr('y', -52)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.yLabel)
  }

  // ── Discrete legend ───────────────────────────────────────────────────────
  const lgX  = width + 20
  const lgY0 = height / 2 - (RESID_BINS.length * 22) / 2
  const swatchS = 14

  // Legend title
  g.append('text')
    .attr('x', lgX).attr('y', lgY0 - 14)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('font-weight', '600')
    .attr('fill', theme.text)
    .text('Standardised')
  g.append('text')
    .attr('x', lgX).attr('y', lgY0 - 2)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('font-weight', '600')
    .attr('fill', theme.text)
    .text('residual')

  RESID_BINS.forEach((bin, bi) => {
    const ly = lgY0 + bi * 22

    // Colour swatch
    g.append('rect')
      .attr('x', lgX).attr('y', ly)
      .attr('width', swatchS).attr('height', swatchS)
      .attr('fill', bin.fill)
      .attr('rx', 2)

    // Significance border on swatch
    if (bin.borderStyle !== 'none') {
      g.append('rect')
        .attr('x', lgX).attr('y', ly)
        .attr('width', swatchS).attr('height', swatchS)
        .attr('fill', 'none').attr('rx', 2)
        .attr('stroke', bin.fill === '#2166ac' || bin.fill === '#74add1' ? '#4292c6' : '#cb181d')
        .attr('stroke-width', bin.borderStyle === 'thick' ? 2 : 1.5)
        .attr('stroke-dasharray', bin.borderStyle === 'thick' ? 'none' : '4,2')
    }

    // Label
    g.append('text')
      .attr('x', lgX + swatchS + 6).attr('y', ly + swatchS - 3)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text)
      .text(bin.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
