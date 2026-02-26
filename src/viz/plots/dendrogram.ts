/**
 * Dendrogram visualization for hierarchical agglomerative clustering.
 * Renders a U-shaped elbow dendrogram with cluster coloring, cut line,
 * leaf labels, and hover tooltips.
 *
 * Uses manual coordinate computation (not d3.hierarchy) because dendrograms
 * need continuous y-axis (merge height) and specific DFS leaf ordering,
 * which d3.cluster/d3.tree cannot represent correctly.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { HACMerge } from '../../stats/clustering.js'

// ─── Data & Config ───────────────────────────────────────────────────────

export interface DendrogramData {
  readonly merges: readonly HACMerge[]
  readonly heights: readonly number[]
  readonly dendrogramOrder: readonly number[]
  readonly labels: readonly number[]
  readonly k: number
  readonly linkage: string
  readonly copheneticCorrelation: number
  readonly observationLabels?: readonly string[]
}

export interface DendrogramConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly orientation?: 'vertical' | 'horizontal'
  readonly showCutLine?: boolean
  readonly showLabels?: boolean
  readonly colorSubtrees?: boolean
}

// ─── Public entry point ──────────────────────────────────────────────────

/**
 * Render a dendrogram into the given container.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data - merge data from runHierarchical + cutTree
 * @param config - visual configuration
 */
export function renderDendrogram(
  container: HTMLElement,
  data: DendrogramData,
  config: DendrogramConfig = {}
): void {
  import('d3').then(d3 => renderDendrogramD3(d3, container, data, config))
}

// ─── Core renderer ───────────────────────────────────────────────────────

function renderDendrogramD3(
  d3: typeof D3,
  container: HTMLElement,
  data: DendrogramData,
  config: DendrogramConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 600
  const H = config.height ?? 480
  const showCut = config.showCutLine !== false
  const showLabels = config.showLabels !== false
  const colorSubs = config.colorSubtrees !== false

  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: showLabels ? Math.max(theme.marginBottom, 100) : theme.marginBottom,
    left: theme.marginLeft,
  }
  const plotW = W - margin.left - margin.right
  const plotH = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  const { merges, heights, dendrogramOrder, labels, k } = data
  const n = merges.length + 1

  // ── Node coordinates ───────────────────────────────────────────────────

  // nodeX / nodeY for all 2n-1 nodes (0..n-1 are leaves, n..2n-2 are internal)
  const nodeX = new Float64Array(2 * n - 1)
  const nodeY = new Float64Array(2 * n - 1)

  // Leaf x positions: index in dendrogramOrder → uniform spacing
  const leafPos = new Float64Array(n)
  for (let i = 0; i < dendrogramOrder.length; i++) {
    leafPos[dendrogramOrder[i]!] = i
  }
  for (let i = 0; i < n; i++) {
    nodeX[i] = leafPos[i]!
    nodeY[i] = 0
  }

  // Internal node positions from merges
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s]!
    const internalIdx = n + s
    nodeX[internalIdx] = (nodeX[m.a]! + nodeX[m.b]!) / 2
    nodeY[internalIdx] = heights[s]!
  }

  // ── Cluster assignment per node (for subtree coloring) ─────────────────

  const nodeCluster = new Int32Array(2 * n - 1)
  // Leaves: use labels directly
  for (let i = 0; i < n; i++) {
    nodeCluster[i] = labels[i]!
  }
  // Internal: same cluster if both children share one, else -1 (mixed)
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s]!
    const cA = nodeCluster[m.a]!
    const cB = nodeCluster[m.b]!
    nodeCluster[n + s] = cA === cB ? cA : -1
  }

  // ── Subtree sizes (for tooltip) ────────────────────────────────────────

  const subtreeSize = new Float64Array(2 * n - 1)
  for (let i = 0; i < n; i++) subtreeSize[i] = 1
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s]!
    subtreeSize[n + s] = subtreeSize[m.a]! + subtreeSize[m.b]!
  }

  // ── Scales ─────────────────────────────────────────────────────────────

  const maxHeight = Math.max(...heights) * 1.05
  const xScale = d3.scaleLinear().domain([0, n - 1]).range([0, plotW])
  const yScale = d3.scaleLinear().domain([0, maxHeight]).range([plotH, 0])

  // ── SVG setup ──────────────────────────────────────────────────────────

  const svg = d3.select(container)
    .append('svg')
    .attr('width', W)
    .attr('height', H)
    .style('background', theme.background)

  const subtitle = `${data.linkage} linkage, K = ${k}, cophenetic r = ${data.copheneticCorrelation.toFixed(3)}`
  addSubtitle(svg, config.title ?? 'Dendrogram', subtitle, W, theme)

  if (config.caption) {
    addCaption(svg, config.caption, W, H, theme)
  }

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`)

  // ── Y-axis (height) ────────────────────────────────────────────────────

  const yAxis = d3.axisLeft(yScale).ticks(6)
  g.append('g')
    .attr('class', 'y-axis')
    .call(yAxis)
    .call(sel => {
      sel.select('.domain').attr('stroke', theme.axisLine)
      sel.selectAll('.tick line').attr('stroke', theme.gridLine)
      sel.selectAll('.tick text')
        .attr('fill', theme.textMuted)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
    })

  // Horizontal gridlines
  g.append('g')
    .attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(6).tickSize(-plotW).tickFormat(() => ''))
    .call(sel => {
      sel.select('.domain').remove()
      sel.selectAll('.tick line')
        .attr('stroke', theme.gridLine)
        .attr('stroke-opacity', 0.5)
    })

  // Y-axis label
  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -plotH / 2)
    .attr('y', -margin.left + 16)
    .attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily)
    .attr('font-size', theme.fontSizeSmall)
    .attr('fill', theme.textMuted)
    .text('Merge Height')

  // ── U-shaped elbow links ───────────────────────────────────────────────

  for (let s = 0; s < merges.length; s++) {
    const m = merges[s]!
    const internalIdx = n + s

    const ax = xScale(nodeX[m.a]!)
    const ay = yScale(nodeY[m.a]!)
    const bx = xScale(nodeX[m.b]!)
    const by = yScale(nodeY[m.b]!)
    const my = yScale(nodeY[internalIdx]!)

    // U-shape: from A up to merge height, across, down to B
    const pathD = `M${ax},${ay} V${my} H${bx} V${by}`

    // Link color: cluster color if pure subtree, neutral gray if mixed
    const cluster = nodeCluster[internalIdx]!
    const color = colorSubs && cluster >= 0
      ? getColor(cluster, theme)
      : theme.axisLine

    const sizeA = subtreeSize[m.a]!
    const sizeB = subtreeSize[m.b]!
    const sizeC = subtreeSize[internalIdx]!
    const mergeH = heights[s]!

    g.append('path')
      .attr('d', pathD)
      .attr('fill', 'none')
      .attr('stroke', color)
      .attr('stroke-width', 1.5)
      .attr('stroke-linejoin', 'round')
      .style('cursor', 'pointer')
      .on('mouseenter', (event: MouseEvent) => {
        const content = [
          formatTooltipRow('Merge height', mergeH.toFixed(3)),
          formatTooltipRow('Left subtree', String(sizeA)),
          formatTooltipRow('Right subtree', String(sizeB)),
          formatTooltipRow('Combined', String(sizeC)),
        ].join('')
        showTooltip(event, content, theme)
      })
      .on('mousemove', (event: MouseEvent) => {
        const content = [
          formatTooltipRow('Merge height', mergeH.toFixed(3)),
          formatTooltipRow('Left subtree', String(sizeA)),
          formatTooltipRow('Right subtree', String(sizeB)),
          formatTooltipRow('Combined', String(sizeC)),
        ].join('')
        showTooltip(event, content, theme)
      })
      .on('mouseleave', () => hideTooltip())
  }

  // ── Cut line ───────────────────────────────────────────────────────────

  if (showCut && k > 1 && k <= n) {
    // Cut between the (n-k)-th and (n-k+1)-th merge heights
    // merges are sorted by height ascending: index n-k-1 is last merge kept, n-k is first merge cut
    const belowIdx = n - k - 1
    const aboveIdx = n - k
    const hBelow = belowIdx >= 0 ? heights[belowIdx]! : 0
    const hAbove = aboveIdx < heights.length ? heights[aboveIdx]! : maxHeight
    const cutY = yScale((hBelow + hAbove) / 2)

    g.append('line')
      .attr('x1', -8)
      .attr('x2', plotW + 8)
      .attr('y1', cutY)
      .attr('y2', cutY)
      .attr('stroke', theme.textMuted)
      .attr('stroke-width', 1.5)
      .attr('stroke-dasharray', '6,4')
      .attr('opacity', 0.7)

    g.append('text')
      .attr('x', plotW + 4)
      .attr('y', cutY - 6)
      .attr('font-family', theme.fontFamilyMono)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.textMuted)
      .attr('text-anchor', 'end')
      .text(`K = ${k}`)
  }

  // ── Leaf labels ────────────────────────────────────────────────────────

  if (showLabels) {
    const obsLabels = data.observationLabels
    // For large n, thin out labels to prevent overlap
    const step = n > 80 ? Math.ceil(n / 40) : 1

    for (let i = 0; i < dendrogramOrder.length; i++) {
      if (i % step !== 0) continue
      const leafIdx = dendrogramOrder[i]!
      const lx = xScale(i)
      const ly = plotH + 8

      const label = obsLabels?.[leafIdx] ?? String(leafIdx + 1)
      const cluster = labels[leafIdx]!

      g.append('text')
        .attr('x', lx)
        .attr('y', ly)
        .attr('transform', `rotate(-60, ${lx}, ${ly})`)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', Math.min(theme.fontSizeSmall, Math.max(8, plotW / n * 0.8)))
        .attr('fill', colorSubs ? getColor(cluster, theme) : theme.textMuted)
        .attr('text-anchor', 'end')
        .attr('dominant-baseline', 'auto')
        .text(label)
    }
  }
}
