/**
 * Treemap — rectangle area proportional to value, laid out using D3's
 * squarified treemap algorithm. Cells are labelled (truncated when small)
 * and coloured by group when provided, otherwise by index.
 * White borders separate cells. Tooltip shows label, value, and %.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface TreemapConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface TreemapChild {
  readonly label: string
  readonly value: number
  readonly group?: string
}

export interface TreemapData {
  readonly children: readonly TreemapChild[]
  readonly testResult?: StatResult
}

/**
 * Render a treemap using D3's squarified layout.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - leaf nodes with labels, values, optional groups + stat result
 * @param config - visual configuration
 */
export function renderTreemap(
  container: HTMLElement,
  data: TreemapData,
  config: TreemapConfig = {}
): void {
  import('d3').then(d3 => renderTreemapD3(d3, container, data, config))
}

function renderTreemapD3(
  d3: typeof D3,
  container: HTMLElement,
  data: TreemapData,
  config: TreemapConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
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

  addSubtitle(svg, config.title ?? 'Treemap', data.testResult?.formatted ?? '', W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  if (data.children.length === 0) return

  // Determine group colour mapping
  const uniqueGroups = Array.from(new Set(data.children.map(c => c.group ?? ''))).filter(Boolean)
  const groupColorMap = new Map(uniqueGroups.map((grp, i) => [grp, getColor(i, theme)]))

  const total = data.children.reduce((s, c) => s + c.value, 0)

  // Build hierarchy and layout
  type HierarchyDatum = { readonly children: readonly TreemapChild[] }
  const root = (d3.hierarchy as unknown as (datum: HierarchyDatum) => D3.HierarchyNode<HierarchyDatum>)(
    { children: data.children }
  )
    .sum((d: unknown) => {
      const node = d as Partial<TreemapChild>
      return node.value ?? 0
    })
    .sort((a, b) => (b.value ?? 0) - (a.value ?? 0))

  const treemapLayout = d3.treemap<HierarchyDatum>()
    .size([width, height])
    .paddingInner(2)
    .paddingOuter(0)
    .round(true)

  treemapLayout(root)

  // Draw leaves
  const leaves = root.leaves() as unknown as Array<D3.HierarchyRectangularNode<HierarchyDatum>>

  leaves.forEach((leaf, li) => {
    const childData = leaf.data as unknown as TreemapChild
    const x0 = leaf.x0
    const y0 = leaf.y0
    const x1 = leaf.x1
    const y1 = leaf.y1
    const cellW = x1 - x0
    const cellH = y1 - y0

    const color = childData.group && groupColorMap.has(childData.group)
      ? groupColorMap.get(childData.group)!
      : getColor(li, theme)

    const pct = total > 0 ? (childData.value / total) * 100 : 0

    const cell = g.append('g').attr('transform', `translate(${x0},${y0})`)

    cell.append('rect')
      .attr('width', cellW)
      .attr('height', cellH)
      .attr('fill', color)
      .attr('opacity', 0.75)
      .attr('stroke', theme.background)
      .attr('stroke-width', 2)
      .attr('rx', 2)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('Label', childData.label),
          formatTooltipRow('Value', childData.value.toFixed(2)),
          formatTooltipRow('Share', `${pct.toFixed(1)}%`),
          ...(childData.group ? [formatTooltipRow('Group', childData.group)] : []),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)

    // Label: only if cell is large enough
    if (cellW > 30 && cellH > 16) {
      const maxChars = Math.max(3, Math.floor(cellW / 7))
      const label = childData.label.length > maxChars
        ? childData.label.slice(0, maxChars - 1) + '…'
        : childData.label

      cell.append('text')
        .attr('x', 5).attr('y', 14)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', Math.min(theme.fontSizeSmall, cellH / 3))
        .attr('font-weight', '600')
        .attr('fill', '#fff')
        .attr('pointer-events', 'none')
        .text(label)

      if (cellH > 28) {
        cell.append('text')
          .attr('x', 5).attr('y', 26)
          .attr('font-family', theme.fontFamily)
          .attr('font-size', Math.min(theme.fontSizeSmall - 2, cellH / 4))
          .attr('fill', 'rgba(255,255,255,0.75)')
          .attr('pointer-events', 'none')
          .text(`${pct.toFixed(1)}%`)
      }
    }
  })

  // Group legend if groups are present
  if (uniqueGroups.length > 0) {
    uniqueGroups.forEach((grp, i) => {
      const lgX = 0
      const lgY = height + 24 + i * 18
      if (lgY + 12 > H - margin.top - 4) return  // no overflow
      g.append('rect')
        .attr('x', lgX).attr('y', lgY)
        .attr('width', 12).attr('height', 12)
        .attr('fill', groupColorMap.get(grp) ?? theme.gridLine)
        .attr('rx', 2)
      g.append('text')
        .attr('x', lgX + 16).attr('y', lgY + 10)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text)
        .text(grp)
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
