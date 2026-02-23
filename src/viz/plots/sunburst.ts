/**
 * Sunburst chart — radial hierarchy visualization where each ring represents
 * one level of depth. Arc angles are proportional to node values. Center hole
 * shows total. Depth-1 nodes are listed in a right-side legend.
 * Inspired by ggstatsplot philosophy: every plot tells the full statistical story.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface SunburstNode {
  readonly name: string
  readonly value?: number
  readonly children?: readonly SunburstNode[]
}

export interface SunburstData {
  readonly root: SunburstNode
  readonly testResult?: StatResult
}

export interface SunburstConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly innerRadius?: number   // fraction of outerRadius — 0 = full pie, default 0.25
  readonly maxDepth?: number      // max rings to render, default 4
}

export function renderSunburst(
  container: HTMLElement,
  data: SunburstData,
  config: SunburstConfig = {}
): void {
  import('d3').then(d3 => renderSunburstD3(d3, container, data, config))
}

// Internal type that carries partition layout coordinates on each node datum
interface PartitionDatum {
  readonly name: string
  readonly value?: number
  readonly children?: readonly SunburstNode[]
}

function renderSunburstD3(
  d3: typeof D3,
  container: HTMLElement,
  data: SunburstData,
  config: SunburstConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 640, 400)
  const H = config.height ?? 520
  const maxDepth = config.maxDepth ?? 4
  const innerRadiusFrac = config.innerRadius ?? 0.25

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('background', theme.background)

  addSubtitle(
    svg,
    config.title ?? 'Sunburst Chart',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  // Reserve right margin for legend, caption at bottom
  const captionH = config.caption ? 24 : 0
  const legendW = 150
  const availW = W - legendW - 8
  const availH = H - theme.marginTop - theme.marginBottom - captionH

  // Center the sunburst in the available area
  const cx = theme.marginLeft + (availW - theme.marginLeft) / 2
  const cy = theme.marginTop + availH / 2

  const outerRadius = Math.min(availW - theme.marginLeft, availH) / 2 - 4
  const innerRadius = outerRadius * innerRadiusFrac

  if (outerRadius <= 0) return

  // Build hierarchy
  const hierarchyRoot = d3.hierarchy<PartitionDatum>(
    data.root as PartitionDatum,
    (d: PartitionDatum) => d.children as PartitionDatum[] | undefined
  ).sum((d: PartitionDatum) => d.value ?? 0)

  const total = hierarchyRoot.value ?? 0

  // Partition layout — size([2π, outerRadius]) maps x to angles, y to radii
  const partitionLayout = d3.partition<PartitionDatum>().size([2 * Math.PI, outerRadius])
  const partitioned = partitionLayout(hierarchyRoot)

  // Collect all nodes as HierarchyRectangularNode (partition sets x0/x1/y0/y1)
  const nodes = partitioned.descendants() as D3.HierarchyRectangularNode<PartitionDatum>[]

  // Build parent-index map for color variation: maps a node to its index among siblings
  const siblingIndexMap = new Map<D3.HierarchyRectangularNode<PartitionDatum>, number>()
  partitioned.eachBefore((node) => {
    const rect = node as D3.HierarchyRectangularNode<PartitionDatum>
    const parent = rect.parent as D3.HierarchyRectangularNode<PartitionDatum> | null
    if (parent !== null) {
      const siblings = (parent.children ?? []) as D3.HierarchyRectangularNode<PartitionDatum>[]
      const idx = siblings.indexOf(rect)
      siblingIndexMap.set(rect, idx >= 0 ? idx : 0)
    }
  })

  // Arc generator — innerRadius/outerRadius come from the partition y0/y1 values
  const arcGen = d3.arc<D3.HierarchyRectangularNode<PartitionDatum>>()
    .startAngle(d => d.x0)
    .endAngle(d => d.x1)
    .innerRadius(d => Math.max(d.y0, innerRadius))
    .outerRadius(d => d.y1)
    .padAngle(0.005)
    .padRadius(outerRadius / 2)

  const g = svg.append('g').attr('transform', `translate(${cx},${cy})`)

  // Depth-1 names for the legend
  const depth1Nodes = nodes.filter(n => n.depth === 1)

  // Draw arcs
  nodes.forEach(node => {
    if (node.depth === 0) return                   // root — skip (center hole)
    if (node.depth > maxDepth) return              // beyond max depth
    if (node.y0 < innerRadius) return              // inside the hole

    const sibIdx = siblingIndexMap.get(node) ?? 0
    const colorIdx = (node.depth - 1 + sibIdx) % 8
    const fillColor = getColor(colorIdx, theme)
    const opacity = Math.max(0.4, 1 - node.depth * 0.15)

    const arcAngle = node.x1 - node.x0
    const arcWidth = node.y1 - node.y0
    const nodeValue = node.value ?? 0
    const parentValue = (node.parent?.value ?? 0)
    const pctOfParent = parentValue > 0 ? ((nodeValue / parentValue) * 100).toFixed(1) : '—'

    g.append('path')
      .datum(node)
      .attr('d', arcGen)
      .attr('fill', fillColor)
      .attr('opacity', opacity)
      .attr('stroke', theme.background)
      .attr('stroke-width', 1.5)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('Name', node.data.name),
          formatTooltipRow('Value', nodeValue.toLocaleString()),
          formatTooltipRow('% of parent', pctOfParent + '%'),
          formatTooltipRow('Depth', String(node.depth)),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Text label — only if arc is wide enough
    if (arcAngle > 0.15 && arcWidth > 20) {
      const midAngle = (node.x0 + node.x1) / 2
      const midRadius = (Math.max(node.y0, innerRadius) + node.y1) / 2
      const midAngleDeg = (midAngle * 180 / Math.PI) - 90
      // Flip text on the left half so it reads correctly
      const flip = midAngle > Math.PI
      const rotateDeg = flip ? midAngleDeg + 180 : midAngleDeg

      const maxChars = Math.max(2, Math.floor((arcAngle * midRadius) / 7))
      const label = node.data.name.length > maxChars
        ? node.data.name.slice(0, maxChars - 1) + '…'
        : node.data.name

      g.append('text')
        .attr('transform', `rotate(${rotateDeg}) translate(${midRadius},0) rotate(${flip ? 180 : 0})`)
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'middle')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', Math.min(theme.fontSizeSmall, arcWidth * 0.5))
        .attr('fill', '#fff')
        .attr('pointer-events', 'none')
        .text(label)
    }
  })

  // Center label — total value + "total"
  if (innerRadius > 12) {
    g.append('text')
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'auto')
      .attr('y', -4)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeTitle)
      .attr('font-weight', '700')
      .attr('fill', theme.text)
      .text(total.toLocaleString())

    g.append('text')
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'hanging')
      .attr('y', 6)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted)
      .text('total')
  }

  // Legend — depth-1 nodes on the right side
  const legendX = cx + outerRadius + 16
  const legendStartY = theme.marginTop + 8
  const lineH = 20

  depth1Nodes.forEach((node, i) => {
    const sibIdx = siblingIndexMap.get(node) ?? 0
    const colorIdx = sibIdx % 8
    const fillColor = getColor(colorIdx, theme)
    const ly = legendStartY + i * lineH

    if (ly + lineH > H - theme.marginBottom - captionH) return

    svg.append('rect')
      .attr('x', legendX)
      .attr('y', ly)
      .attr('width', 12)
      .attr('height', 12)
      .attr('fill', fillColor)
      .attr('rx', 2)

    const maxLabelChars = Math.max(3, Math.floor((W - legendX - 28) / 7))
    const rawLabel = node.data.name
    const label = rawLabel.length > maxLabelChars
      ? rawLabel.slice(0, maxLabelChars - 1) + '…'
      : rawLabel

    svg.append('text')
      .attr('x', legendX + 16)
      .attr('y', ly + 10)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .text(label)
  })

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
