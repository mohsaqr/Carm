/**
 * Hierarchical edge bundling (Holten 2006).
 *
 * Leaf nodes are arranged on a circle. Non-leaf nodes form the hierarchy used
 * to bundle edges — they are not drawn. Edges between leaf nodes are drawn as
 * cubic B-spline paths routed through the common ancestor path in the
 * hierarchy, which causes edges that share ancestry to visually cluster
 * ("bundle") together.
 *
 * Layout:
 *   d3.stratify builds the hierarchy from the flat node list.
 *   d3.cluster assigns angular positions (x) and radial depth (y).
 *   Only leaf nodes (no children) are shown as circles with radial labels.
 *   Each edge walks the path [source → LCA → target] through the tree and
 *   passes those anchor points to d3.lineRadial with d3.curveBundle.beta(β).
 *
 * Colour:
 *   Nodes are coloured by their `group` field if provided, otherwise by their
 *   immediate parent id. The same colour mapping is reused for edges (source
 *   colour at low opacity so overlapping edges remain readable).
 *
 * Interactivity:
 *   Hovering a node highlights all its edges (full opacity, thicker stroke)
 *   and dims all other edges. Clicking a node locks the highlight; clicking
 *   again unlocks it.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

// ─── Public interfaces ────────────────────────────────────────────────────────

export interface BundleNode {
  readonly id: string
  readonly parent: string   // id of parent node; root has parent = ""
  readonly label?: string
  readonly group?: string   // optional group for colouring
}

export interface BundleEdge {
  readonly source: string   // leaf node id
  readonly target: string   // leaf node id
}

export interface EdgeBundlingData {
  readonly nodes: readonly BundleNode[]
  readonly edges: readonly BundleEdge[]
  readonly testResult?: StatResult
}

export interface EdgeBundlingConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly bundlingStrength?: number  // beta for d3.curveBundle, 0–1, default 0.85
  readonly nodeRadius?: number        // leaf node dot radius, default 4
}

// ─── Public entry-point ───────────────────────────────────────────────────────

/**
 * Render a hierarchical edge bundling diagram.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes (with parent hierarchy), edges between leaves, optional stat result
 * @param config    - visual configuration
 */
export function renderEdgeBundling(
  container: HTMLElement,
  data: EdgeBundlingData,
  config: EdgeBundlingConfig = {}
): void {
  import('d3').then(d3 => renderEdgeBundlingD3(d3, container, data, config))
}

// ─── Internal types ───────────────────────────────────────────────────────────

// d3.stratify does not expose a generic-typed HierarchyNode with .path() in
// all versions, so we use a minimal typed wrapper that we can annotate.
interface RadialHierarchyNode {
  data: BundleNode
  x: number          // angle in degrees (from d3.cluster)
  y: number          // radius (from d3.cluster)
  children?: RadialHierarchyNode[]
  parent: RadialHierarchyNode | null
  // ancestors() is present on all d3-hierarchy nodes
  ancestors(): RadialHierarchyNode[]
}

// ─── D3 render function ───────────────────────────────────────────────────────

function renderEdgeBundlingD3(
  d3: typeof D3,
  container: HTMLElement,
  data: EdgeBundlingData,
  config: EdgeBundlingConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 640, 400)
  const H = config.height ?? Math.max(W, 480)
  const beta = config.bundlingStrength ?? 0.85
  const nodeR = config.nodeRadius ?? 4

  // Radial layout fills the inner square minus generous padding for labels
  const labelPad = 90
  const radius = Math.min(W, H) / 2 - labelPad

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
    config.title ?? 'Edge Bundling',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  // Centre the radial layout in the SVG
  const cx = W / 2
  const cy = H / 2 + (theme.marginTop / 2)

  const g = svg.append('g').attr('transform', `translate(${cx},${cy})`)

  if (data.nodes.length === 0) return

  // ── 1. Colour mapping ────────────────────────────────────────────────────────
  // Use group if provided, otherwise parent id. Collect unique values and map
  // them to palette colours.
  const colorKey = (n: BundleNode): string => n.group ?? n.parent ?? ''
  const uniqueColorKeys = Array.from(new Set(data.nodes.map(colorKey)))
  const colorKeyMap = new Map(uniqueColorKeys.map((k, i) => [k, getColor(i, theme)]))
  const nodeColor = (n: BundleNode): string => colorKeyMap.get(colorKey(n)) ?? getColor(0, theme)

  // ── 2. Build d3 hierarchy ────────────────────────────────────────────────────
  // d3.stratify requires a root node whose parentId returns null/undefined.
  // Our root node has parent = "".
  type RawNode = { id: string; parentId: string | undefined; data: BundleNode }
  const rawNodes: RawNode[] = data.nodes.map(n => ({
    id: n.id,
    parentId: n.parent === '' ? undefined : n.parent,
    data: n,
  }))

  let root: RadialHierarchyNode
  try {
    const stratify = (d3.stratify as unknown as (opts: {
      id: (d: RawNode) => string
      parentId: (d: RawNode) => string | undefined
    }) => (data: RawNode[]) => RadialHierarchyNode)({
      id: (d: RawNode) => d.id,
      parentId: (d: RawNode) => d.parentId,
    })
    root = stratify(rawNodes)
  } catch {
    // If stratify fails (bad hierarchy), render empty
    g.append('text')
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.textMuted)
      .text('Invalid hierarchy — check node parent ids.')
    return
  }

  // ── 3. Cluster layout: leaf nodes on a circle ───────────────────────────────
  // d3.cluster assigns .x (angle in [0, 360)) and .y (depth in [0, radius]).
  // We fix inner nodes at y=0 and leaves at y=radius.
  interface D3ClusterLayout {
    size(sz: [number, number]): D3ClusterLayout
    (root: RadialHierarchyNode): RadialHierarchyNode
  }
  const cluster = (d3.cluster as unknown as () => D3ClusterLayout)()
  cluster.size([360, radius])
  cluster(root)

  // ── 4. Index nodes by id ────────────────────────────────────────────────────
  const nodeIndex = new Map<string, RadialHierarchyNode>()
  const allNodes: RadialHierarchyNode[] = []

  // Walk the hierarchy using a queue (no for loops — use array methods)
  ;(function collectNodes(node: RadialHierarchyNode): void {
    allNodes.push(node)
    nodeIndex.set(node.data.id, node)
    node.children?.forEach(collectNodes)
  })(root)

  const leafNodes = allNodes.filter(n => !n.children || n.children.length === 0)

  // ── 5. LCA path helper ──────────────────────────────────────────────────────
  // Returns the path [a, ...up_to_LCA, ...down_to_b] as an array of nodes.
  function lcaPath(
    a: RadialHierarchyNode,
    b: RadialHierarchyNode
  ): RadialHierarchyNode[] {
    const ancestorsA = a.ancestors()          // [a, parent(a), ..., root]
    const ancestorsB = b.ancestors()          // [b, parent(b), ..., root]
    const ancestorSetA = new Set(ancestorsA.map(n => n.data.id))

    // Find LCA: first ancestor of b that is also an ancestor of a
    const lca = ancestorsB.find(n => ancestorSetA.has(n.data.id))
    if (!lca) return [a, b]

    // Path from a up to LCA (inclusive on both ends)
    const lcaIdx = ancestorsA.findIndex(n => n.data.id === lca.data.id)
    const pathUp = ancestorsA.slice(0, lcaIdx + 1)

    // Path from LCA down to b (exclusive LCA, inclusive b)
    const lcaIdxB = ancestorsB.findIndex(n => n.data.id === lca.data.id)
    const pathDown = ancestorsB.slice(0, lcaIdxB).reverse()

    return [...pathUp, ...pathDown]
  }

  // ── 6. Line generator: radial + cubic bundle ────────────────────────────────
  interface D3LineRadial {
    curve(c: unknown): D3LineRadial
    angle(fn: (d: RadialHierarchyNode) => number): D3LineRadial
    radius(fn: (d: RadialHierarchyNode) => number): D3LineRadial
    (pts: RadialHierarchyNode[]): string | null
  }
  const bundleCurve = (d3.curveBundle as unknown as { beta: (b: number) => unknown }).beta(beta)
  const lineRadial = (d3.lineRadial as unknown as () => D3LineRadial)()
    .curve(bundleCurve)
    .angle((d: RadialHierarchyNode) => (d.x * Math.PI) / 180)
    .radius((d: RadialHierarchyNode) => d.y)

  // ── 7. Draw edges ────────────────────────────────────────────────────────────
  const edgeGroup = g.append('g').attr('class', 'edges')
  const nodeGroup = g.append('g').attr('class', 'nodes')

  // Build edge paths and store them for hover interaction
  type EdgePath = {
    srcId: string
    tgtId: string
    el: D3.Selection<SVGPathElement, unknown, null, undefined>
    srcColor: string
  }
  const edgePaths: EdgePath[] = []

  data.edges.forEach(e => {
    const srcNode = nodeIndex.get(e.source)
    const tgtNode = nodeIndex.get(e.target)
    if (!srcNode || !tgtNode) return

    const path = lcaPath(srcNode, tgtNode)
    const pathStr = (lineRadial as unknown as (pts: RadialHierarchyNode[]) => string | null)(path)
    if (!pathStr) return

    const color = nodeColor(srcNode.data)

    const el = edgeGroup.append('path')
      .attr('d', pathStr)
      .attr('fill', 'none')
      .attr('stroke', color)
      .attr('stroke-width', 1)
      .attr('opacity', 0.28)

    edgePaths.push({ srcId: e.source, tgtId: e.target, el, srcColor: color })
  })

  // ── 8. Draw leaf nodes and labels ────────────────────────────────────────────
  let lockedNodeId: string | null = null

  leafNodes.forEach(n => {
    const angleRad = (n.x * Math.PI) / 180
    const nx = Math.sin(angleRad) * n.y
    const ny = -Math.cos(angleRad) * n.y
    const color = nodeColor(n.data)
    const label = n.data.label ?? n.data.id

    // Node circle
    nodeGroup.append('circle')
      .attr('cx', nx)
      .attr('cy', ny)
      .attr('r', nodeR)
      .attr('fill', color)
      .attr('stroke', theme.background)
      .attr('stroke-width', 1.5)
      .attr('cursor', 'pointer')
      .on('mouseover', (event: MouseEvent) => {
        if (lockedNodeId !== null) return
        highlightNode(n.data.id)
        showTooltip(event, [
          formatTooltipRow('Node', label),
          formatTooltipRow('Group', n.data.group ?? n.data.parent ?? '—'),
          formatTooltipRow(
            'Connections',
            data.edges.filter(e => e.source === n.data.id || e.target === n.data.id).length.toString()
          ),
        ].join(''), theme)
      })
      .on('mouseout', () => {
        if (lockedNodeId !== null) return
        resetHighlight()
        hideTooltip()
      })
      .on('click', (event: MouseEvent) => {
        event.stopPropagation()
        if (lockedNodeId === n.data.id) {
          lockedNodeId = null
          resetHighlight()
          hideTooltip()
        } else {
          lockedNodeId = n.data.id
          highlightNode(n.data.id)
          showTooltip(event, [
            formatTooltipRow('Node', label),
            formatTooltipRow('Group', n.data.group ?? n.data.parent ?? '—'),
            formatTooltipRow(
              'Connections',
              data.edges.filter(e => e.source === n.data.id || e.target === n.data.id).length.toString()
            ),
          ].join(''), theme)
        }
      })

    // Radial label: flip text on the left half so it always reads outward
    const labelAngle = n.x              // degrees
    const flipped = labelAngle > 180
    const rotDeg = flipped ? labelAngle - 90 : labelAngle - 90
    const labelR = n.y + nodeR + 5

    const lx = Math.sin(angleRad) * labelR
    const ly = -Math.cos(angleRad) * labelR

    // Rotation: align text along the radial line, flip if on left half
    const textRotate = flipped
      ? `rotate(${rotDeg + 180},${lx.toFixed(2)},${ly.toFixed(2)})`
      : `rotate(${rotDeg},${lx.toFixed(2)},${ly.toFixed(2)})`

    nodeGroup.append('text')
      .attr('x', lx)
      .attr('y', ly)
      .attr('dy', '0.35em')
      .attr('text-anchor', flipped ? 'end' : 'start')
      .attr('transform', textRotate)
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text)
      .attr('pointer-events', 'none')
      .text(label)
  })

  // Click background to unlock
  svg.on('click', () => {
    if (lockedNodeId !== null) {
      lockedNodeId = null
      resetHighlight()
      hideTooltip()
    }
  })

  // ── 9. Highlight helpers ─────────────────────────────────────────────────────
  function highlightNode(nodeId: string): void {
    edgePaths.forEach(ep => {
      const connected = ep.srcId === nodeId || ep.tgtId === nodeId
      ep.el
        .attr('opacity', connected ? 0.9 : 0.05)
        .attr('stroke-width', connected ? 2 : 0.8)
    })
  }

  function resetHighlight(): void {
    edgePaths.forEach(ep => {
      ep.el.attr('opacity', 0.28).attr('stroke-width', 1)
    })
  }

  // ── 10. Colour legend ────────────────────────────────────────────────────────
  if (uniqueColorKeys.filter(k => k !== '').length > 1) {
    const legendKeys = uniqueColorKeys.filter(k => k !== '')
    const lgX = -W / 2 + (theme.marginLeft)
    const lgY = H / 2 - theme.marginBottom - (legendKeys.length * 18)

    legendKeys.forEach((k, i) => {
      g.append('rect')
        .attr('x', lgX)
        .attr('y', lgY + i * 18)
        .attr('width', 10)
        .attr('height', 10)
        .attr('fill', colorKeyMap.get(k) ?? getColor(i, theme))
        .attr('rx', 2)
      g.append('text')
        .attr('x', lgX + 14)
        .attr('y', lgY + i * 18 + 9)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('fill', theme.text)
        .text(k)
    })
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
