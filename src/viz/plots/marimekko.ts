/**
 * Marimekko (Mekko) chart — a 100% stacked bar chart where column WIDTHS are
 * proportional to each column's total value, so both axes encode data.
 *
 * - Column width  = that category's total as a fraction of the grand total.
 * - Within each column, stacked segments sum to 100%, coloured by series.
 * - Cells receive percentage (or raw value) labels when large enough.
 * - A "spine" bar below the x-axis shows relative column widths at a glance.
 *
 * Usage:
 *   renderMarimekko(container, data, config)
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

// ─── Public interfaces ────────────────────────────────────────────────────────

export interface MarimekkoData {
  /** Labels for each column (X variable). */
  readonly categories: readonly string[]
  /** One series per stacked segment. values[j] corresponds to categories[j]. */
  readonly series: readonly {
    readonly label: string
    readonly values: readonly number[]
  }[]
  readonly testResult?: StatResult
}

export interface MarimekkoConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  /** Show percentage labels inside cells when the cell is large enough. Default: true. */
  readonly showPercentLabels?: boolean
  /** Show raw values instead of percentages in cell labels. Default: false. */
  readonly showValueLabels?: boolean
}

// ─── Entry point ──────────────────────────────────────────────────────────────

export function renderMarimekko(
  container: HTMLElement,
  data: MarimekkoData,
  config: MarimekkoConfig = {}
): void {
  import('d3').then(d3 => renderMarimekkoD3(d3, container, data, config))
}

// ─── D3 render ────────────────────────────────────────────────────────────────

function renderMarimekkoD3(
  d3: typeof D3,
  container: HTMLElement,
  data: MarimekkoData,
  config: MarimekkoConfig
): void {
  if (data.categories.length === 0 || data.series.length === 0) return

  const theme           = config.theme ?? DEFAULT_THEME
  const showPct         = config.showPercentLabels !== false
  const showVal         = config.showValueLabels === true

  // ── Layout constants ───────────────────────────────────────────────────────
  const legendW  = 130          // reserved width on the right for the legend
  const spineH   = 6            // height of the column-width spine bar
  const spineGap = 44           // vertical offset below x-axis for spine
  const xLabelY  = spineGap + spineH + 18  // y offset for optional x-axis label

  const W = config.width  ?? Math.max(container.clientWidth || 640, 420)
  const H = config.height ?? 500

  const margin = {
    top:    theme.marginTop,
    right:  theme.marginRight + legendW,
    bottom: theme.marginBottom + spineH + 20, // extra room for spine + x label
    left:   theme.marginLeft,
  }
  const width  = W - margin.left - margin.right
  const height = H - margin.top  - margin.bottom

  // ── Column totals and grand total ──────────────────────────────────────────
  const nCols = data.categories.length
  const colTotals: number[] = Array.from({ length: nCols }, (_, j) =>
    data.series.reduce((sum, s) => sum + (s.values[j] ?? 0), 0)
  )
  const grandTotal = colTotals.reduce((a, b) => a + b, 0)
  if (grandTotal === 0) return

  // ── Column widths and x positions (cumulative) ─────────────────────────────
  const gap = 3
  const colWidths: number[] = colTotals.map(t => (t / grandTotal) * width)
  const colX: number[] = []
  colWidths.reduce((acc, w, i) => { colX[i] = acc; return acc + w }, 0)

  // ── Bootstrap SVG ─────────────────────────────────────────────────────────
  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width',   W)
    .attr('height',  H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(
    svg,
    config.title ?? 'Marimekko Chart',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // ── Y-axis gridlines (0 / 25 / 50 / 75 / 100 %) ───────────────────────────
  const yPcts = [0, 25, 50, 75, 100]
  yPcts.forEach(pct => {
    const yPos = height * (1 - pct / 100)

    // Dashed gridline
    g.append('line')
      .attr('x1', 0).attr('x2', width)
      .attr('y1', yPos).attr('y2', yPos)
      .attr('stroke', theme.gridLine)
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '3,3')

    // Tick label
    g.append('text')
      .attr('x', -8)
      .attr('y', yPos + 4)
      .attr('text-anchor', 'end')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted)
      .text(`${pct}%`)
  })

  // Y-axis label (left, rotated)
  if (config.yLabel) {
    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -height / 2)
      .attr('y', -52)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.yLabel)
  } else {
    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -height / 2)
      .attr('y', -52)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.textMuted)
      .text('Percentage (%)')
  }

  // ── Draw cells ─────────────────────────────────────────────────────────────
  data.categories.forEach((cat, ci) => {
    const colTotal = colTotals[ci] ?? 0
    if (colTotal === 0) return

    const cx = colX[ci] ?? 0
    const cw = colWidths[ci] ?? 0

    // Accumulate cumulative proportion from bottom (series drawn bottom→top)
    let cumPct = 0

    data.series.forEach((s, si) => {
      const rawVal  = s.values[ci] ?? 0
      const cellPct = rawVal / colTotal          // proportion within column (0–1)
      const cellH   = cellPct * height

      const rx = cx + gap / 2
      const ry = height - (cumPct + cellPct) * height
      const rw = Math.max(cw - gap, 0)
      const rh = Math.max(cellH - gap / 2, 0)

      const color = getColor(si, theme)

      // Cell rectangle
      g.append('rect')
        .attr('x', rx).attr('y', ry)
        .attr('width', rw).attr('height', rh)
        .attr('fill', color)
        .attr('opacity', theme.violinOpacity)

      // Cell label (% or raw value) — only if large enough
      const labelStr = showVal
        ? String(rawVal)
        : `${(cellPct * 100).toFixed(1)}%`

      if (rh > 18 && rw > 30 && (showPct || showVal)) {
        const isDark = isDarkColor(color)
        g.append('text')
          .attr('x', rx + rw / 2)
          .attr('y', ry + rh / 2 + 4)
          .attr('text-anchor', 'middle')
          .attr('font-family', theme.fontFamilyMono)
          .attr('font-size', Math.min(theme.fontSizeSmall, rh * 0.38, rw * 0.18))
          .attr('fill', isDark ? '#ffffff' : '#212529')
          .attr('pointer-events', 'none')
          .text(labelStr)
      }

      // Invisible hit-area for tooltip
      const colPct = (colTotal / grandTotal * 100).toFixed(1)
      g.append('rect')
        .attr('x', rx).attr('y', ry)
        .attr('width', rw).attr('height', rh)
        .attr('fill', 'transparent')
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Category',     cat),
            formatTooltipRow('Series',       s.label),
            formatTooltipRow('Value',        rawVal),
            formatTooltipRow('% of column',  `${(cellPct * 100).toFixed(1)}%`),
            formatTooltipRow('Column total', `${colTotal} (${colPct}% of total)`),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)

      cumPct += cellPct
    })
  })

  // ── X-axis: column labels + (N=xxx) totals ─────────────────────────────────
  data.categories.forEach((cat, ci) => {
    const cx = colX[ci] ?? 0
    const cw = colWidths[ci] ?? 0
    const midX = cx + cw / 2

    // Category label
    g.append('text')
      .attr('x', midX)
      .attr('y', height + 18)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(cat)

    // N label in muted colour
    g.append('text')
      .attr('x', midX)
      .attr('y', height + 30)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.textMuted)
      .text(`(N=${colTotals[ci] ?? 0})`)
  })

  // ── Spine bar (column-width indicator) ─────────────────────────────────────
  data.categories.forEach((_, ci) => {
    const cx = colX[ci] ?? 0
    const cw = colWidths[ci] ?? 0
    g.append('rect')
      .attr('x', cx + gap / 2)
      .attr('y', height + spineGap)
      .attr('width', Math.max(cw - gap, 0))
      .attr('height', spineH)
      .attr('fill', getColor(0, theme))
      .attr('opacity', 0.4)
      .attr('rx', 1)
  })

  // ── Optional X-axis label (below spine) ────────────────────────────────────
  if (config.xLabel) {
    g.append('text')
      .attr('x', width / 2)
      .attr('y', height + xLabelY)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.xLabel)
  }

  // ── Right-side legend ──────────────────────────────────────────────────────
  const swatchS = 12
  const lgX     = width + 20
  const lgY0    = Math.max(0, height / 2 - (data.series.length * 20) / 2)

  // Legend title
  g.append('text')
    .attr('x', lgX)
    .attr('y', lgY0 - 12)
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall - 1)
    .attr('font-weight', '600')
    .attr('fill', theme.text)
    .text('Series')

  data.series.forEach((s, si) => {
    const color = getColor(si, theme)
    const ly    = lgY0 + si * 20

    g.append('rect')
      .attr('x', lgX).attr('y', ly)
      .attr('width', swatchS).attr('height', swatchS)
      .attr('fill', color)
      .attr('opacity', theme.violinOpacity)
      .attr('rx', 2)

    g.append('text')
      .attr('x', lgX + swatchS + 6)
      .attr('y', ly + swatchS - 2)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall - 1)
      .attr('fill', theme.text)
      .text(s.label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/**
 * Returns true if a hex color is perceptually dark (luminance < 0.35),
 * so overlaid text should be white rather than near-black.
 * Supports 6-digit hex strings (e.g. '#4e79a7').
 */
function isDarkColor(hex: string): boolean {
  const r = parseInt(hex.slice(1, 3), 16) / 255
  const g = parseInt(hex.slice(3, 5), 16) / 255
  const b = parseInt(hex.slice(5, 7), 16) / 255
  // Relative luminance (WCAG 2.1 formula, approximate linearisation)
  const lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
  return lum < 0.35
}
