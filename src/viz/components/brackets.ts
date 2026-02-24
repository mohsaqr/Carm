/**
 * Significance bracket layout for group comparison plots.
 * Renders stacked brackets above groups with p-value labels.
 * Non-overlapping bracket stacking algorithm.
 */

import type * as d3 from 'd3'
import type { PairwiseResult } from '../../core/types.js'
import type { CarmTheme } from '../themes/default.js'
import { DEFAULT_THEME } from '../themes/default.js'

export interface BracketConfig {
  readonly groupPositions: ReadonlyMap<string, number>  // group label → x-center
  readonly yBase: number              // y coordinate of the top of data (plot top)
  readonly bracketHeight: number      // vertical space per bracket level
  readonly significantOnly: boolean   // show only p < 0.05 brackets
  readonly numericP?: boolean         // true = "p = .025", false = "***" stars (default true)
}

/** Format p-value for bracket label. */
function formatBracketP(p: number, numeric: boolean): string {
  if (numeric) {
    if (p < 0.001) return 'p < .001'
    return `p = ${p.toFixed(3).replace(/^0\./, '.')}`
  }
  if (p < 0.001) return '***'
  if (p < 0.01)  return '**'
  if (p < 0.05)  return '*'
  return 'ns'
}

/**
 * Render significance brackets onto an SVG group.
 * Uses greedy level-assignment to avoid bracket overlap.
 */
export function renderBrackets(
  g: d3.Selection<SVGGElement, unknown, null, undefined>,
  comparisons: readonly PairwiseResult[],
  config: BracketConfig,
  theme: CarmTheme = DEFAULT_THEME
): void {
  const toShow = config.significantOnly
    ? comparisons.filter(c => c.significant)
    : comparisons

  if (toShow.length === 0) return

  // Greedy level assignment: assign each bracket to the lowest level that doesn't overlap
  interface Bracket {
    x1: number
    x2: number
    p: number
    label: string
    level: number
  }

  const brackets: Bracket[] = toShow.map(c => ({
    x1: Math.min(
      config.groupPositions.get(c.group1) ?? 0,
      config.groupPositions.get(c.group2) ?? 0
    ),
    x2: Math.max(
      config.groupPositions.get(c.group1) ?? 0,
      config.groupPositions.get(c.group2) ?? 0
    ),
    p: c.pValueAdj,
    label: formatBracketP(c.pValueAdj, config.numericP !== false),
    level: 0,
  })).sort((a, b) => (a.x2 - a.x1) - (b.x2 - b.x1))  // shortest first

  // Assign levels greedily
  const levelMaxX: number[] = []
  for (const bracket of brackets) {
    let level = 0
    while ((levelMaxX[level] ?? -Infinity) >= bracket.x1 - 5) {
      level++
    }
    bracket.level = level
    levelMaxX[level] = bracket.x2
  }

  // Draw brackets
  for (const b of brackets) {
    const y = config.yBase - b.level * config.bracketHeight - 12
    const tipLen = 6
    const color = b.label === 'ns' ? theme.textMuted : theme.text

    // Horizontal line
    g.append('line')
      .attr('x1', b.x1).attr('x2', b.x2)
      .attr('y1', y).attr('y2', y)
      .attr('stroke', color)
      .attr('stroke-width', 1.2)

    // Left tick
    g.append('line')
      .attr('x1', b.x1).attr('x2', b.x1)
      .attr('y1', y).attr('y2', y + tipLen)
      .attr('stroke', color)
      .attr('stroke-width', 1.2)

    // Right tick
    g.append('line')
      .attr('x1', b.x2).attr('x2', b.x2)
      .attr('y1', y).attr('y2', y + tipLen)
      .attr('stroke', color)
      .attr('stroke-width', 1.2)

    // Label
    g.append('text')
      .attr('x', (b.x1 + b.x2) / 2)
      .attr('y', y - 3)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily)
      .attr('font-size', theme.fontSizeSmall)
      .attr('fill', color)
      .text(b.label)
  }
}

/** Total bracket height needed (to set top margin). */
export function totalBracketHeight(
  comparisons: readonly PairwiseResult[],
  significantOnly: boolean,
  bracketHeight: number
): number {
  const n = significantOnly ? comparisons.filter(c => c.significant).length : comparisons.length
  return n * bracketHeight + 20
}
