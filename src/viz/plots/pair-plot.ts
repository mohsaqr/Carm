/**
 * Scatter matrix (pairs plot) — n×n grid of mini-plots.
 * Off-diagonal cells: scatter plot of variable i vs variable j.
 * Diagonal cells: histogram of variable i.
 * Each cell is a self-contained SVG <g> fitting in the available space.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface PairPlotConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface PairPlotData {
  /** data[varIndex][obsIndex] — each inner array is one variable's observations */
  readonly data: readonly (readonly number[])[]
  readonly labels: readonly string[]
  readonly testResult?: StatResult
}

/**
 * Render a scatter matrix (pairs plot).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - variable arrays + labels + optional stat result
 * @param config - visual configuration
 */
export function renderPairPlot(
  container: HTMLElement,
  data: PairPlotData,
  config: PairPlotConfig = {}
): void {
  import('d3').then(d3 => renderPairPlotD3(d3, container, data, config))
}

function renderPairPlotD3(
  d3: typeof D3,
  container: HTMLElement,
  data: PairPlotData,
  config: PairPlotConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? Math.max(W, 480)

  // Outer margins just for title/caption
  const outerMarginTop = 52
  const outerMarginBottom = config.caption ? 28 : 16
  const outerMarginSide = 16

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Scatter Matrix', data.testResult?.formatted ?? '', W, theme)

  const n = data.labels.length
  if (n < 2) return

  const gridW = W - outerMarginSide * 2
  const gridH = H - outerMarginTop - outerMarginBottom
  const cellW = gridW / n
  const cellH = gridH / n
  const cellPad = 6   // inner padding per cell

  // Pre-compute per-variable scales (shared across the whole column/row)
  const scales = data.data.map(vals => {
    const [lo, hi] = d3.extent(vals) as [number, number]
    const pad = ((hi - lo) * 0.08) || 0.5
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([cellPad, cellW - cellPad]).nice()
  })

  const scalesY = data.data.map(vals => {
    const [lo, hi] = d3.extent(vals) as [number, number]
    const pad = ((hi - lo) * 0.08) || 0.5
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([cellH - cellPad, cellPad]).nice()
  })

  const gridG = svg.append('g')
    .attr('transform', `translate(${outerMarginSide},${outerMarginTop})`)

  // Draw n×n cells
  for (let row = 0; row < n; row++) {
    for (let col = 0; col < n; col++) {
      const cellX = col * cellW
      const cellY = row * cellH
      const cellG = gridG.append('g')
        .attr('transform', `translate(${cellX},${cellY})`)

      // Cell background
      cellG.append('rect')
        .attr('width', cellW).attr('height', cellH)
        .attr('fill', row === col ? theme.surface : theme.background)
        .attr('stroke', theme.gridLine).attr('stroke-width', 0.5)

      if (row === col) {
        // Diagonal: histogram of variable i
        drawHistogramCell(d3, cellG, data.data[row]!, scales[col]!, cellW, cellH, cellPad, getColor(row, theme), theme)

        // Variable label centred
        cellG.append('text')
          .attr('x', cellW / 2).attr('y', cellPad + 10)
          .attr('text-anchor', 'middle')
          .attr('font-family', theme.fontFamily)
          .attr('font-size', Math.min(theme.fontSizeSmall, cellW / 8))
          .attr('font-weight', '600')
          .attr('fill', theme.text)
          .text(data.labels[row] ?? '')
      } else {
        // Off-diagonal: scatter plot of col (x) vs row (y)
        const xVals = data.data[col]!
        const yVals = data.data[row]!
        const xSc = scales[col]!
        const ySc = scalesY[row]!

        const n_obs = Math.min(xVals.length, yVals.length)
        for (let i = 0; i < n_obs; i++) {
          const xv = xVals[i]!
          const yv = yVals[i]!
          cellG.append('circle')
            .attr('cx', xSc(xv))
            .attr('cy', ySc(yv))
            .attr('r', Math.min(2.5, cellW / 40))
            .attr('fill', getColor(col, theme))
            .attr('opacity', theme.pointOpacity)
            .on('mouseover', (event: MouseEvent) => {
              showTooltip(event, [
                formatTooltipRow(data.labels[col] ?? `Var ${col}`, xv.toFixed(3)),
                formatTooltipRow(data.labels[row] ?? `Var ${row}`, yv.toFixed(3)),
              ].join(''), theme)
            })
            .on('mouseout', hideTooltip)
        }

        // Pearson correlation label
        const r = pearsonR(xVals.slice(0, n_obs), yVals.slice(0, n_obs))
        cellG.append('text')
          .attr('x', cellW - cellPad - 2).attr('y', cellH - cellPad - 2)
          .attr('text-anchor', 'end')
          .attr('font-family', theme.fontFamilyMono)
          .attr('font-size', Math.min(theme.fontSizeSmall - 1, cellW / 9))
          .attr('fill', theme.textMuted)
          .text(`r=${r.toFixed(2)}`)
      }
    }
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

function drawHistogramCell(
  d3: typeof D3,
  g: D3.Selection<SVGGElement, unknown, null, undefined>,
  values: readonly number[],
  xScale: D3.ScaleLinear<number, number>,
  _cellW: number,
  cellH: number,
  pad: number,
  color: string,
  _theme: CarmTheme
): void {
  const bins = d3.bin()
    .domain(xScale.domain() as [number, number])
    .thresholds(8)([...values])

  const maxCount = Math.max(...bins.map(b => b.length))
  const yScale = d3.scaleLinear().domain([0, maxCount]).range([cellH - pad - 14, pad + 14])

  bins.forEach(bin => {
    const x0 = xScale(bin.x0 ?? 0)
    const x1 = xScale(bin.x1 ?? 0)
    const bw = Math.max(x1 - x0 - 1, 1)
    g.append('rect')
      .attr('x', x0)
      .attr('y', yScale(bin.length))
      .attr('width', bw)
      .attr('height', Math.max(cellH - pad - 14 - yScale(bin.length), 0))
      .attr('fill', color)
      .attr('opacity', 0.6)
  })
}

function pearsonR(xs: readonly number[], ys: readonly number[]): number {
  const n = xs.length
  if (n < 2) return 0
  const mx = xs.reduce((a, b) => a + b, 0) / n
  const my = ys.reduce((a, b) => a + b, 0) / n
  let num = 0, sdx = 0, sdy = 0
  for (let i = 0; i < n; i++) {
    const dx = (xs[i]! - mx)
    const dy = (ys[i]! - my)
    num += dx * dy
    sdx += dx * dx
    sdy += dy * dy
  }
  const denom = Math.sqrt(sdx * sdy)
  return denom === 0 ? 0 : num / denom
}
