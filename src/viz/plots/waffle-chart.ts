/**
 * Waffle chart — 10×10 grid of small squares where each square represents 1%
 * of the total. Squares are coloured by slice proportionally (floor-rounded,
 * remainder assigned to largest slice). Legend below shows colour, label, %.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface WaffleChartConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface WaffleChartData {
  readonly slices: readonly { readonly label: string; readonly value: number }[]
  readonly testResult?: StatResult
}

/**
 * Render a waffle chart (10×10 proportional grid).
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - named slices with values + optional stat result
 * @param config - visual configuration
 */
export function renderWaffleChart(
  container: HTMLElement,
  data: WaffleChartData,
  config: WaffleChartConfig = {}
): void {
  import('d3').then(d3 => renderWaffleChartD3(d3, container, data, config))
}

function renderWaffleChartD3(
  d3: typeof D3,
  container: HTMLElement,
  data: WaffleChartData,
  config: WaffleChartConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 400, 300)
  const H = config.height ?? 420

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Waffle Chart', data.testResult?.formatted ?? '', W, theme)

  if (data.slices.length === 0) return

  const total = data.slices.reduce((s, sl) => s + sl.value, 0)
  if (total === 0) return

  // Allocate 100 cells: floor then assign remainder to largest remainders
  const pcts = data.slices.map(sl => (sl.value / total) * 100)
  const floors = pcts.map(p => Math.floor(p))
  const remainders = pcts.map((p, i) => p - floors[i]!)
  let remaining = 100 - floors.reduce((a, b) => a + b, 0)

  // Sort by remainder descending to assign leftover cells
  const sortedIdx = Array.from({ length: data.slices.length }, (_, i) => i)
    .sort((a, b) => remainders[b]! - remainders[a]!)

  const cells = [...floors]
  for (let i = 0; i < remaining; i++) {
    cells[sortedIdx[i % sortedIdx.length]!]! += 1
  }

  // Build flat cell array: 100 entries, each holding the slice index
  const cellColors: number[] = []
  cells.forEach((count, si) => {
    for (let c = 0; c < count; c++) cellColors.push(si)
  })

  // Grid geometry
  const legendH = Math.ceil(data.slices.length / 3) * 22 + 12
  const gridAreaH = H - theme.marginTop - legendH - 16
  const gridAreaW = W - theme.marginLeft - theme.marginRight
  const cellSize = Math.min(Math.floor(gridAreaH / 10), Math.floor(gridAreaW / 10))
  const gap = Math.max(2, Math.floor(cellSize / 8))
  const gridW = cellSize * 10 + gap * 9
  const gridH = cellSize * 10 + gap * 9
  const originX = (W - gridW) / 2
  const originY = theme.marginTop + 8

  const g = svg.append('g')

  // Draw 10×10 grid (row 0 = top)
  for (let row = 0; row < 10; row++) {
    for (let col = 0; col < 10; col++) {
      const idx = row * 10 + col
      const si = cellColors[idx] ?? 0
      const color = getColor(si, theme)
      const x = originX + col * (cellSize + gap)
      const y = originY + row * (cellSize + gap)
      const sliceData = data.slices[si]!

      g.append('rect')
        .attr('x', x).attr('y', y)
        .attr('width', cellSize).attr('height', cellSize)
        .attr('fill', color)
        .attr('opacity', 0.85)
        .attr('rx', Math.max(1, cellSize / 8))
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Slice', sliceData.label),
            formatTooltipRow('Value', sliceData.value.toFixed(2)),
            formatTooltipRow('Share', `${((sliceData.value / total) * 100).toFixed(1)}%`),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)
    }
  }

  // Legend — up to 3 items per row
  const lgTop = originY + gridH + 20
  const lgColW = Math.floor(W / 3)

  data.slices.forEach((sl, si) => {
    const row = Math.floor(si / 3)
    const col = si % 3
    const lgX = col * lgColW + 12
    const lgY = lgTop + row * 22
    const color = getColor(si, theme)
    const pct = ((sl.value / total) * 100).toFixed(1)

    g.append('rect')
      .attr('x', lgX).attr('y', lgY)
      .attr('width', 12).attr('height', 12)
      .attr('fill', color).attr('opacity', 0.85)
      .attr('rx', 2)

    const labelText = sl.label.length > 14 ? sl.label.slice(0, 13) + '…' : sl.label
    g.append('text')
      .attr('x', lgX + 16).attr('y', lgY + 10)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text)
      .text(`${labelText} (${pct}%)`)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
