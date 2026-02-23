/**
 * Arc diagram: nodes on a horizontal baseline connected by arcs above them.
 * Short connections produce shallow arcs; distant connections produce tall arcs.
 * Useful for visualising network topology when node order carries meaning.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface ArcDiagramNode {
  readonly id: string
  readonly label?: string
  readonly group?: string
}

export interface ArcDiagramEdge {
  readonly source: string
  readonly target: string
  readonly value?: number
}

export interface ArcDiagramData {
  readonly nodes: readonly ArcDiagramNode[]
  readonly edges: readonly ArcDiagramEdge[]
  readonly testResult?: StatResult
}

export interface ArcDiagramConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly nodeRadius?: number
  readonly sortByGroup?: boolean
}

export function renderArcDiagram(
  container: HTMLElement,
  data: ArcDiagramData,
  config: ArcDiagramConfig = {}
): void {
  import('d3').then(d3 => renderArcDiagramD3(d3, container, data, config))
}

function renderArcDiagramD3(
  d3: typeof D3,
  container: HTMLElement,
  data: ArcDiagramData,
  config: ArcDiagramConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 700, 400)
  const H = config.height ?? 420
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft,
  }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  const nodeRadius = config.nodeRadius ?? 6
  const sortByGroup = config.sortByGroup ?? true

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .style('background', theme.background)

  addSubtitle(
    svg,
    config.title ?? 'Arc Diagram',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // ── Determine node order ──────────────────────────────────────────────────
  // Collect unique groups in the order they first appear, then sort nodes
  const uniqueGroups: string[] = []
  data.nodes.forEach(n => {
    const grp = n.group ?? '__none__'
    if (!uniqueGroups.includes(grp)) uniqueGroups.push(grp)
  })

  const groupIndex = (grp: string | undefined): number =>
    uniqueGroups.indexOf(grp ?? '__none__')

  const orderedNodes: ArcDiagramNode[] = sortByGroup
    ? [...data.nodes].sort((a, b) => {
        const gi = groupIndex(a.group) - groupIndex(b.group)
        if (gi !== 0) return gi
        // Within the same group, preserve original order deterministically
        return data.nodes.indexOf(a) - data.nodes.indexOf(b)
      })
    : [...data.nodes]

  // ── Positions ─────────────────────────────────────────────────────────────
  // Nodes sit on a horizontal baseline
  const nodeY = height * 0.65

  const n = orderedNodes.length
  // Evenly space nodes; for a single node, place it at center
  const xStep = n > 1 ? width / (n - 1) : width / 2
  const xOffset = n > 1 ? 0 : width / 2

  const nodeX = new Map<string, number>()
  orderedNodes.forEach((node, i) => {
    nodeX.set(node.id, xOffset + i * xStep)
  })

  // ── Degree map ────────────────────────────────────────────────────────────
  const degree = new Map<string, number>()
  data.nodes.forEach(node => degree.set(node.id, 0))
  data.edges.forEach(edge => {
    degree.set(edge.source, (degree.get(edge.source) ?? 0) + 1)
    degree.set(edge.target, (degree.get(edge.target) ?? 0) + 1)
  })

  // ── Baseline grid line ────────────────────────────────────────────────────
  g.append('line')
    .attr('x1', 0)
    .attr('x2', width)
    .attr('y1', nodeY)
    .attr('y2', nodeY)
    .attr('stroke', theme.gridLine)
    .attr('stroke-width', 1.5)

  // ── Arcs ─────────────────────────────────────────────────────────────────
  data.edges.forEach(edge => {
    const x1 = nodeX.get(edge.source)
    const x2 = nodeX.get(edge.target)
    if (x1 == null || x2 == null || edge.source === edge.target) return

    const srcNode = data.nodes.find(nd => nd.id === edge.source)
    const srcGroupIdx = groupIndex(srcNode?.group)
    const arcColor = getColor(srcGroupIdx, theme)

    // Arc height proportional to horizontal distance between nodes
    const arcHeight = (Math.abs(x2 - x1) / width) * (height * 0.55)

    const mx = (x1 + x2) / 2
    const my = nodeY - arcHeight

    const rawStroke = 1 + (edge.value ?? 1) * 0.5
    const strokeWidth = Math.min(rawStroke, 4)

    const pathD = `M ${x1},${nodeY} Q ${mx},${my} ${x2},${nodeY}`

    g.append('path')
      .attr('d', pathD)
      .attr('fill', 'none')
      .attr('stroke', arcColor)
      .attr('stroke-width', strokeWidth)
      .attr('opacity', 0.55)
      .on('mouseover', (event: MouseEvent) => {
        const srcLabel =
          data.nodes.find(nd => nd.id === edge.source)?.label ?? edge.source
        const tgtLabel =
          data.nodes.find(nd => nd.id === edge.target)?.label ?? edge.target
        const content = [
          formatTooltipRow('From', srcLabel),
          formatTooltipRow('To', tgtLabel),
          ...(edge.value != null ? [formatTooltipRow('Value', edge.value)] : []),
        ].join('')
        showTooltip(event, content, theme)
      })
      .on('mouseout', hideTooltip)
  })

  // ── Nodes ────────────────────────────────────────────────────────────────
  const rotateLabelThreshold = 8

  orderedNodes.forEach(node => {
    const x = nodeX.get(node.id) ?? 0
    const grpIdx = groupIndex(node.group)
    const color = getColor(grpIdx, theme)
    const displayLabel = node.label ?? node.id
    const deg = degree.get(node.id) ?? 0

    // Node circle
    g.append('circle')
      .attr('cx', x)
      .attr('cy', nodeY)
      .attr('r', nodeRadius)
      .attr('fill', color)
      .attr('stroke', theme.background)
      .attr('stroke-width', 1.5)
      .on('mouseover', (event: MouseEvent) => {
        const groupName = node.group ?? 'none'
        const content = [
          formatTooltipRow('Node', displayLabel),
          formatTooltipRow('Group', groupName),
          formatTooltipRow('Degree', deg),
        ].join('')
        showTooltip(event, content, theme)
      })
      .on('mouseout', hideTooltip)

    // Label below node
    if (n <= rotateLabelThreshold) {
      // Straight label, centered below node
      g.append('text')
        .attr('x', x)
        .attr('y', nodeY + nodeRadius + 14)
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'hanging')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSize)
        .attr('fill', theme.text)
        .text(displayLabel)
    } else {
      // Rotated 45° to avoid overlap when many nodes
      g.append('text')
        .attr('transform', `translate(${x},${nodeY + nodeRadius + 8}) rotate(45)`)
        .attr('text-anchor', 'start')
        .attr('dominant-baseline', 'middle')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text)
        .text(displayLabel)
    }
  })

  // ── x/y axis labels (optional) ───────────────────────────────────────────
  if (config.xLabel) {
    g.append('text')
      .attr('x', width / 2)
      .attr('y', height + 44)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.xLabel)
  }

  if (config.yLabel) {
    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -height / 2)
      .attr('y', -48)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.yLabel)
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
