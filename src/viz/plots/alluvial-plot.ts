/**
 * Alluvial (Sankey-like) plot — nodes stacked by stage, flows as cubic bezier
 * ribbons. Every stage is a vertical column; node height is proportional to
 * total flow value. Ribbons connect adjacent stages and are coloured by source
 * node colour at reduced opacity so layering is visible.
 *
 * Layout:
 *   stages are evenly spaced across the plot width.
 *   within each stage nodes are sorted largest→smallest, stacked top→bottom
 *   with a configurable gap between them.
 *   each node tracks separate y-offsets for outgoing and incoming ribbons so
 *   ribbons stack neatly without crossing inside a node.
 *
 * Colour: all unique node labels are collected up-front and assigned a fixed
 * colour from the theme palette so the same category always has the same
 * colour regardless of which stage it appears in.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

// ─── Public interfaces ────────────────────────────────────────────────────────

export interface AlluvialNode {
  readonly id: string
  readonly stage: number     // 0-based column index
  readonly label: string
}

export interface AlluvialFlow {
  readonly source: string    // node id
  readonly target: string    // node id
  readonly value: number
}

export interface AlluvialData {
  readonly nodes: readonly AlluvialNode[]
  readonly flows: readonly AlluvialFlow[]
  readonly stageLabels?: readonly string[]
  readonly testResult?: StatResult
}

export interface AlluvialConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly nodePadding?: number   // vertical gap between nodes in same stage, default 8
  readonly nodeWidth?: number     // width of node rectangles, default 18
}

// ─── Public entry-point ───────────────────────────────────────────────────────

/**
 * Render an alluvial / Sankey-style plot.
 * @param container - HTMLElement to render into (cleared on call)
 * @param data      - nodes, flows, optional stage labels and stat result
 * @param config    - visual configuration
 */
export function renderAlluvialPlot(
  container: HTMLElement,
  data: AlluvialData,
  config: AlluvialConfig = {}
): void {
  import('d3').then(d3 => renderAlluvialD3(d3, container, data, config))
}

// ─── Internal types ───────────────────────────────────────────────────────────

interface LayoutNode {
  readonly id: string
  readonly stage: number
  readonly label: string
  readonly color: string
  readonly totalValue: number   // sum of all flows in + out (used for sizing)
  readonly inValue: number      // sum of incoming flows (used for height in non-first stages)
  readonly outValue: number     // sum of outgoing flows
  x: number                     // left edge of node rect
  y: number                     // top edge of node rect
  readonly height: number       // pixel height
  outOffset: number             // running y-offset for stacking outgoing ribbons
  inOffset: number              // running y-offset for stacking incoming ribbons
}

interface LayoutFlow {
  readonly sourceNode: LayoutNode
  readonly targetNode: LayoutNode
  readonly value: number
  readonly ribbonHeight: number
  readonly sourceY: number      // top of this ribbon's slot on source node
  readonly targetY: number      // top of this ribbon's slot on target node
}

// ─── D3 render function ───────────────────────────────────────────────────────

function renderAlluvialD3(
  d3: typeof D3,
  container: HTMLElement,
  data: AlluvialData,
  config: AlluvialConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 700, 400)
  const H = config.height ?? 520
  const nodePadding = config.nodePadding ?? 8
  const nodeWidth = config.nodeWidth ?? 18

  // When stage headers are shown, add 24 px so they don't overlap the subtitle
  // (subtitle baseline is at SVG y = 45; headers render at g-coord y = -12, i.e.
  //  SVG y = marginTop - 12; we need marginTop - 12 > 55 → marginTop > 67).
  const hasStageLabels = (data.stageLabels?.length ?? 0) > 0
  const margin = {
    top: theme.marginTop + (hasStageLabels ? 24 : 0),
    right: theme.marginRight + 80,    // extra right room for last stage's labels
    bottom: theme.marginBottom,
    left: theme.marginLeft + 60,      // extra left room for first stage's labels
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

  addSubtitle(
    svg,
    config.title ?? 'Alluvial Plot',
    data.testResult?.formatted ?? '',
    W,
    theme
  )

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  if (data.nodes.length === 0 || data.flows.length === 0) return

  // ── 1. Colour mapping: unique labels get a fixed palette colour ─────────────
  const uniqueLabels = Array.from(new Set(data.nodes.map(n => n.label)))
  const labelColor = new Map(uniqueLabels.map((lbl, i) => [lbl, getColor(i, theme)]))

  // ── 2. Determine stage count and x positions ────────────────────────────────
  const stageIndices = Array.from(new Set(data.nodes.map(n => n.stage))).sort((a, b) => a - b)
  const nStages = stageIndices.length
  const stageX = (stageIdx: number): number =>
    nStages < 2 ? width / 2 : (stageIdx / (nStages - 1)) * (width - nodeWidth)

  // ── 3. Compute per-node flow sums ───────────────────────────────────────────
  const nodeById = new Map(data.nodes.map(n => [n.id, n]))
  const inSum = new Map<string, number>()
  const outSum = new Map<string, number>()

  data.flows.forEach(f => {
    outSum.set(f.source, (outSum.get(f.source) ?? 0) + f.value)
    inSum.set(f.target, (inSum.get(f.target) ?? 0) + f.value)
  })

  // Total value per node: for the first stage, use outgoing; last stage uses incoming;
  // middle stages use max(in, out) so nodes are sized to the larger of the two.
  const firstStage = stageIndices[0] ?? 0
  const lastStage = stageIndices[nStages - 1] ?? 0

  const nodeTotal = (id: string): number => {
    const node = nodeById.get(id)
    if (!node) return 0
    if (node.stage === firstStage) return outSum.get(id) ?? 0
    if (node.stage === lastStage) return inSum.get(id) ?? 0
    return Math.max(outSum.get(id) ?? 0, inSum.get(id) ?? 0)
  }

  // ── 4. Layout nodes within each stage ──────────────────────────────────────
  // Total value across all nodes in a stage determines the height scale.
  const stageNodes = new Map<number, LayoutNode[]>()
  stageIndices.forEach(si => {
    const nodesInStage = data.nodes
      .filter(n => n.stage === si)
      .sort((a, b) => nodeTotal(b.id) - nodeTotal(a.id))

    const totalVal = nodesInStage.reduce((s, n) => s + nodeTotal(n.id), 0)
    const usableHeight = height - nodePadding * Math.max(0, nodesInStage.length - 1)
    const heightScale = totalVal > 0 ? usableHeight / totalVal : 1

    const lnodes: LayoutNode[] = []
    let yCursor = 0

    nodesInStage.forEach(n => {
      const tv = nodeTotal(n.id)
      const nh = Math.max(4, tv * heightScale)
      lnodes.push({
        id: n.id,
        stage: n.stage,
        label: n.label,
        color: labelColor.get(n.label) ?? getColor(0, theme),
        totalValue: tv,
        inValue: inSum.get(n.id) ?? 0,
        outValue: outSum.get(n.id) ?? 0,
        x: stageX(stageIndices.indexOf(si)),
        y: yCursor,
        height: nh,
        outOffset: 0,
        inOffset: 0,
      })
      yCursor += nh + nodePadding
    })

    stageNodes.set(si, lnodes)
  })

  const layoutNodeMap = new Map<string, LayoutNode>()
  stageNodes.forEach(lnodes => lnodes.forEach(ln => layoutNodeMap.set(ln.id, ln)))

  // ── 5. Compute flow layout (ribbon slots) ───────────────────────────────────
  // Scale ribbon height relative to the source node's height / out-value.
  const layoutFlows: LayoutFlow[] = []

  // Sort flows by value desc so large ribbons sit at the top of each node.
  const sortedFlows = [...data.flows].sort((a, b) => b.value - a.value)

  sortedFlows.forEach(f => {
    const src = layoutNodeMap.get(f.source)
    const tgt = layoutNodeMap.get(f.target)
    if (!src || !tgt) return

    const srcHeightScale = src.outValue > 0 ? src.height / src.outValue : 0
    const tgtHeightScale = tgt.inValue > 0 ? tgt.height / tgt.inValue : 0
    const ribbonH = Math.max(1, f.value * Math.min(srcHeightScale, tgtHeightScale))

    const sourceY = src.y + src.outOffset
    const targetY = tgt.y + tgt.inOffset

    layoutFlows.push({ sourceNode: src, targetNode: tgt, value: f.value, ribbonHeight: ribbonH, sourceY, targetY })

    // Advance offsets (mutate LayoutNode offsets — they are tracking state)
    ;(src as { outOffset: number }).outOffset += ribbonH
    ;(tgt as { inOffset: number }).inOffset += ribbonH
  })

  // ── 6. Draw ribbons (below nodes) ───────────────────────────────────────────
  const flowGroup = g.append('g').attr('class', 'flows')
  const nodeGroup = g.append('g').attr('class', 'nodes')

  layoutFlows.forEach(lf => {
    const x1 = lf.sourceNode.x + nodeWidth
    const x2 = lf.targetNode.x
    const cpX1 = x1 + (x2 - x1) * 0.4
    const cpX2 = x1 + (x2 - x1) * 0.6

    const yTop1 = lf.sourceY
    const yBot1 = lf.sourceY + lf.ribbonHeight
    const yTop2 = lf.targetY
    const yBot2 = lf.targetY + lf.ribbonHeight

    // Closed ribbon shape: top bezier forward + bottom bezier backward
    const pathD = [
      `M ${x1.toFixed(2)},${yTop1.toFixed(2)}`,
      `C ${cpX1.toFixed(2)},${yTop1.toFixed(2)} ${cpX2.toFixed(2)},${yTop2.toFixed(2)} ${x2.toFixed(2)},${yTop2.toFixed(2)}`,
      `L ${x2.toFixed(2)},${yBot2.toFixed(2)}`,
      `C ${cpX2.toFixed(2)},${yBot2.toFixed(2)} ${cpX1.toFixed(2)},${yBot1.toFixed(2)} ${x1.toFixed(2)},${yBot1.toFixed(2)}`,
      'Z',
    ].join(' ')

    const srcPct = lf.sourceNode.outValue > 0
      ? ((lf.value / lf.sourceNode.outValue) * 100).toFixed(1)
      : '–'

    flowGroup.append('path')
      .attr('d', pathD)
      .attr('fill', lf.sourceNode.color)
      .attr('opacity', 0.42)
      .attr('stroke', lf.sourceNode.color)
      .attr('stroke-width', 0.5)
      .attr('stroke-opacity', 0.25)
      .on('mouseover', (event: MouseEvent) => {
        showTooltip(event, [
          formatTooltipRow('From', lf.sourceNode.label),
          formatTooltipRow('To', lf.targetNode.label),
          formatTooltipRow('Value', lf.value.toFixed(2)),
          formatTooltipRow('% of source', `${srcPct}%`),
        ].join(''), theme)
      })
      .on('mouseout', hideTooltip)
      // Highlight on hover
      .on('mouseover.highlight', function (event: MouseEvent) {
        d3.select(this as SVGPathElement).attr('opacity', 0.72)
        showTooltip(event, [
          formatTooltipRow('From', lf.sourceNode.label),
          formatTooltipRow('To', lf.targetNode.label),
          formatTooltipRow('Value', lf.value.toFixed(2)),
          formatTooltipRow('% of source', `${srcPct}%`),
        ].join(''), theme)
      })
      .on('mouseout.highlight', function () {
        d3.select(this as SVGPathElement).attr('opacity', 0.42)
        hideTooltip()
      })
  })

  // ── 7. Draw nodes and labels ────────────────────────────────────────────────
  stageIndices.forEach(si => {
    const lnodes = stageNodes.get(si) ?? []
    const isFirst = si === firstStage
    const isLast = si === lastStage

    lnodes.forEach(ln => {
      // Node rectangle
      nodeGroup.append('rect')
        .attr('x', ln.x)
        .attr('y', ln.y)
        .attr('width', nodeWidth)
        .attr('height', ln.height)
        .attr('fill', ln.color)
        .attr('rx', 3)
        .attr('stroke', theme.background)
        .attr('stroke-width', 1.5)
        .on('mouseover', (event: MouseEvent) => {
          showTooltip(event, [
            formatTooltipRow('Label', ln.label),
            formatTooltipRow('Stage', ln.stage.toString()),
            formatTooltipRow('Total value', ln.totalValue.toFixed(2)),
          ].join(''), theme)
        })
        .on('mouseout', hideTooltip)

      // Label placement
      const labelX = isFirst
        ? ln.x - 6                      // left of leftmost stage
        : isLast
          ? ln.x + nodeWidth + 6        // right of rightmost stage
          : ln.x + nodeWidth / 2        // centred above middle stages

      const labelAnchor = isFirst ? 'end' : isLast ? 'start' : 'middle'
      const labelY = ln.y + ln.height / 2 + 4

      nodeGroup.append('text')
        .attr('x', labelX)
        .attr('y', labelY)
        .attr('text-anchor', labelAnchor)
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('font-weight', '500')
        .attr('fill', theme.text)
        .attr('pointer-events', 'none')
        .text(ln.label)
    })
  })

  // ── 8. Stage column headers ─────────────────────────────────────────────────
  if (data.stageLabels) {
    stageIndices.forEach((si, sIdx) => {
      const label = data.stageLabels?.[sIdx] ?? `Stage ${si}`
      const x = stageX(sIdx) + nodeWidth / 2

      nodeGroup.append('text')
        .attr('x', x)
        .attr('y', -12)
        .attr('text-anchor', 'middle')
        .attr('font-family', theme.fontFamily)
        .attr('font-size', theme.fontSizeSmall)
        .attr('font-weight', '700')
        .attr('fill', theme.textMuted)
        .attr('letter-spacing', '0.5')
        .text(label.toUpperCase())
    })
  }

  // ── 9. Optional axis labels ─────────────────────────────────────────────────
  if (config.xLabel) {
    g.append('text')
      .attr('x', width / 2)
      .attr('y', height + 48)
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
      .attr('y', -margin.left + 16)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSize)
      .attr('fill', theme.text)
      .text(config.yLabel)
  }

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
