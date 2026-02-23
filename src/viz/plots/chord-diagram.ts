/**
 * Chord diagram: shows flows between groups as ribbons inside a circular layout.
 * Each group occupies an arc on the outer ring; ribbons inside represent flows.
 * Based on the D3 chord layout (d3.chord / d3.ribbon / d3.arc).
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ChordDiagramData {
  readonly matrix: readonly (readonly number[])[]
  readonly labels: readonly string[]
  readonly testResult?: StatResult
}

export interface ChordDiagramConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly innerRadius?: number
  readonly padAngle?: number
}

export function renderChordDiagram(
  container: HTMLElement,
  data: ChordDiagramData,
  config: ChordDiagramConfig = {}
): void {
  import('d3').then(d3 => renderChordDiagramD3(d3, container, data, config))
}

function renderChordDiagramD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ChordDiagramData,
  config: ChordDiagramConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? Math.max(container.clientHeight || 500, 400)

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .style('background', theme.background)

  addSubtitle(
    svg,
    config.title ?? 'Chord Diagram',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  // Reserve space for subtitle at top (theme.marginTop) and caption at bottom
  const captionH = config.caption ? 24 : 0
  const availH = H - theme.marginTop - captionH - 8

  const cx = W / 2
  const cy = theme.marginTop + availH / 2

  const minDim = Math.min(W, availH)
  const outerRadius = (config.innerRadius != null)
    ? config.innerRadius + 20
    : minDim / 2 * 0.82
  const innerRadius = (config.innerRadius != null)
    ? config.innerRadius
    : outerRadius - 20

  const padAngle = config.padAngle ?? 0.02

  // Build a mutable copy of the matrix for d3.chord (it expects number[][])
  const mutableMatrix = data.matrix.map(row => [...row])

  const chordLayout = d3.chord()
    .padAngle(padAngle)
    .sortSubgroups(d3.descending)

  const chords = chordLayout(mutableMatrix)

  const arcGenerator = d3.arc<D3.ChordGroup>()
    .innerRadius(innerRadius)
    .outerRadius(outerRadius)

  // Use generic overload so TypeScript understands ChordSubgroup as the subgroup datum.
  // The startAngle/endAngle accessors come directly from ChordSubgroup; radius is the
  // constant innerRadius set via the accessor form.
  const ribbonGenerator = d3.ribbon<D3.Chord, D3.ChordSubgroup>()
    .startAngle((sg: D3.ChordSubgroup) => sg.startAngle)
    .endAngle((sg: D3.ChordSubgroup) => sg.endAngle)
    .radius((_sg: D3.ChordSubgroup) => innerRadius)

  const g = svg.append('g')
    .attr('transform', `translate(${cx},${cy})`)

  // ── Ribbons (render before arcs so arcs sit on top) ──────────────────────
  // Sort chords ascending by value so smaller ribbons render on top
  const sortedChords = [...chords].sort(
    (a, b) => b.source.value + b.target.value - (a.source.value + a.target.value)
  )

  g.selectAll<SVGPathElement, D3.Chord>('.ribbon')
    .data(sortedChords)
    .join('path')
    .attr('class', 'ribbon')
    .attr('d', d => ribbonGenerator(d) ?? '')
    .attr('fill', d => getColor(d.source.index, theme))
    .attr('opacity', 0.65)
    .attr('stroke', theme.background)
    .attr('stroke-width', 0.5)
    .on('mouseover', (event: MouseEvent, d: D3.Chord) => {
      const srcLabel = data.labels[d.source.index] ?? `Group ${d.source.index}`
      const tgtLabel = data.labels[d.target.index] ?? `Group ${d.target.index}`
      const content = [
        formatTooltipRow('From', srcLabel),
        formatTooltipRow('To', tgtLabel),
        formatTooltipRow('Flow', d.source.value),
      ].join('')
      showTooltip(event, content, theme)
    })
    .on('mouseout', hideTooltip)

  // ── Group arcs ───────────────────────────────────────────────────────────
  const groupG = g.selectAll<SVGGElement, D3.ChordGroup>('.group')
    .data(chords.groups)
    .join('g')
    .attr('class', 'group')

  groupG.append('path')
    .attr('d', d => arcGenerator(d) ?? '')
    .attr('fill', d => getColor(d.index, theme))
    .attr('stroke', theme.background)
    .attr('stroke-width', 1)
    .on('mouseover', (event: MouseEvent, d: D3.ChordGroup) => {
      const label = data.labels[d.index] ?? `Group ${d.index}`
      // Compute total outflow and inflow from the original matrix
      const outflow = data.matrix[d.index]?.reduce((s, v) => s + v, 0) ?? 0
      const inflow = data.matrix.reduce((s, row) => s + (row[d.index] ?? 0), 0)
      const content = [
        formatTooltipRow('Group', label),
        formatTooltipRow('Outflow', outflow),
        formatTooltipRow('Inflow', inflow),
      ].join('')
      showTooltip(event, content, theme)
    })
    .on('mouseout', hideTooltip)

  // ── Tick marks on group arcs ─────────────────────────────────────────────
  // Compute total sum to derive tick step (every ~5% of total)
  const total = data.matrix.reduce(
    (s, row) => s + row.reduce((rs, v) => rs + v, 0),
    0
  )
  const tickStep = total > 0 ? computeTickStep(total) : 1

  chords.groups.forEach(group => {
    const groupTotal = data.matrix[group.index]?.reduce((s, v) => s + v, 0) ?? 0
    const tickCount = Math.floor(groupTotal / tickStep)

    for (let ti = 0; ti <= tickCount; ti++) {
      const tickValue = ti * tickStep
      // Angle for this tick within the group's arc
      const angleRange = group.endAngle - group.startAngle - padAngle
      const tickAngle = group.startAngle + padAngle / 2 +
        (groupTotal > 0 ? (tickValue / groupTotal) * angleRange : 0)

      const sinA = Math.sin(tickAngle - Math.PI / 2)
      const cosA = Math.cos(tickAngle - Math.PI / 2)

      g.append('line')
        .attr('x1', innerRadius * cosA)
        .attr('y1', innerRadius * sinA)
        .attr('x2', (outerRadius + 8) * cosA)
        .attr('y2', (outerRadius + 8) * sinA)
        .attr('stroke', theme.axisLine)
        .attr('stroke-width', 1)

      // Tick label for non-zero ticks
      if (ti > 0) {
        g.append('text')
          .attr('x', (outerRadius + 14) * cosA)
          .attr('y', (outerRadius + 14) * sinA)
          .attr('text-anchor', tickAngle > Math.PI ? 'end' : 'start')
          .attr('dominant-baseline', 'middle')
          .attr('font-family', theme.fontFamily)
          .attr('font-size', theme.fontSizeSmall - 2)
          .attr('fill', theme.textMuted)
          .text(formatTickValue(tickValue))
      }
    }
  })

  // ── Group labels ─────────────────────────────────────────────────────────
  const labelRadius = outerRadius + 32

  chords.groups.forEach(group => {
    const angle = (group.startAngle + group.endAngle) / 2
    // D3 chord angles are measured from 12 o'clock, clockwise
    const labelAngle = angle - Math.PI / 2
    const lx = labelRadius * Math.cos(labelAngle)
    const ly = labelRadius * Math.sin(labelAngle)
    const rightSide = Math.cos(labelAngle) > 0

    g.append('text')
      .attr('x', lx)
      .attr('y', ly)
      .attr('text-anchor', rightSide ? 'start' : 'end')
      .attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('font-weight', '600')
      .attr('fill', theme.text)
      .text(data.labels[group.index] ?? `Group ${group.index}`)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

/** Compute a round tick step so there are roughly 5–10 ticks across the total. */
function computeTickStep(total: number): number {
  const rough = total / 8
  const magnitude = Math.pow(10, Math.floor(Math.log10(rough)))
  const norm = rough / magnitude
  const step = norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10
  return step * magnitude
}

/** Format a tick value — suppress decimals for whole numbers. */
function formatTickValue(v: number): string {
  if (v >= 1000) return `${(v / 1000).toFixed(v % 1000 === 0 ? 0 : 1)}k`
  return Number.isInteger(v) ? String(v) : v.toFixed(1)
}
