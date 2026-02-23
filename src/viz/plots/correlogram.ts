/**
 * Correlation heatmap (correlogram) with significance stars and clustering.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { CorrelationMatrix } from '../../stats/correlation.js'

export interface CorrelogramConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
  readonly showValues?: boolean
  readonly showSignificance?: boolean
}

export function renderCorrelogram(
  container: HTMLElement,
  data: CorrelationMatrix,
  config: CorrelogramConfig = {}
): void {
  import('d3').then(d3 => renderCorrelogramD3(d3, container, data, config))
}

function renderCorrelogramD3(
  d3: typeof D3,
  container: HTMLElement,
  data: CorrelationMatrix,
  config: CorrelogramConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const k = data.labels.length
  const cellSize = Math.min(
    Math.floor((config.width ?? 500) / (k + 2)),
    Math.floor((config.height ?? 500) / (k + 2)),
    70
  )
  const W = cellSize * k + 120
  const H = cellSize * k + 100

  container.innerHTML = ''
  applyTheme(container, theme)

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).style('background', theme.background)

  addSubtitle(svg, config.title ?? 'Correlation Matrix', '', W, theme)

  const margin = { top: 60, right: 20, bottom: 20, left: 80 }
  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  // Color scale: blue (−1) → white → red (+1)
  const colorScale = d3.scaleSequential()
    .domain([-1, 1])
    .interpolator(d3.interpolateRdBu)

  // Draw cells
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) {
      const r = data.r[i]?.[j] ?? 0
      const p = data.pValues[i]?.[j] ?? 1

      g.append('rect')
        .attr('x', j * cellSize).attr('y', i * cellSize)
        .attr('width', cellSize - 2).attr('height', cellSize - 2)
        .attr('rx', 3)
        .attr('fill', i === j ? theme.surface : colorScale(r))
        .on('mouseover', (event: MouseEvent) => {
          if (i !== j) {
            showTooltip(event, [
              formatTooltipRow(`${data.labels[i]} × ${data.labels[j]}`, ''),
              formatTooltipRow('r', r.toFixed(3)),
              formatTooltipRow('p', p < 0.001 ? '< .001' : p.toFixed(3)),
            ].join(''), theme)
          }
        })
        .on('mouseout', hideTooltip)

      // Value text
      if (config.showValues !== false && i !== j) {
        g.append('text')
          .attr('x', j * cellSize + cellSize / 2)
          .attr('y', i * cellSize + cellSize / 2 + 1)
          .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
          .attr('font-family', theme.fontFamilyMono).attr('font-size', Math.min(theme.fontSizeSmall, cellSize / 3))
          .attr('fill', Math.abs(r) > 0.6 ? '#fff' : theme.text)
          .text(r.toFixed(2))
      }

      // Significance stars
      if (config.showSignificance !== false && i !== j && !isNaN(p)) {
        const stars = p < 0.001 ? '***' : p < 0.01 ? '**' : p < 0.05 ? '*' : ''
        if (stars) {
          g.append('text')
            .attr('x', j * cellSize + cellSize / 2)
            .attr('y', i * cellSize + cellSize - 4)
            .attr('text-anchor', 'middle')
            .attr('font-size', 8)
            .attr('fill', Math.abs(r) > 0.6 ? '#fff' : theme.textMuted)
            .text(stars)
        }
      }
    }
  }

  // Row labels (left)
  data.labels.forEach((lbl, i) => {
    g.append('text')
      .attr('x', -6).attr('y', i * cellSize + cellSize / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(lbl)
  })

  // Column labels (top)
  data.labels.forEach((lbl, j) => {
    g.append('text')
      .attr('x', j * cellSize + cellSize / 2).attr('y', -8)
      .attr('text-anchor', 'middle')
      .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSizeSmall)
      .attr('fill', theme.text).text(lbl)
  })

  // Color legend (gradient bar)
  const legendW = 80, legendH = 12
  const legendX = W - 120, legendY = H - 35
  const defs = svg.append('defs')
  const gradId = 'corr-grad-' + Math.random().toString(36).slice(2)
  const grad = defs.append('linearGradient').attr('id', gradId)
  const stops = [-1, -0.5, 0, 0.5, 1]
  stops.forEach(v => {
    grad.append('stop')
      .attr('offset', `${((v + 1) / 2 * 100).toFixed(0)}%`)
      .attr('stop-color', colorScale(v))
  })
  svg.append('rect').attr('x', legendX).attr('y', legendY)
    .attr('width', legendW).attr('height', legendH)
    .attr('fill', `url(#${gradId})`)
  svg.append('text').attr('x', legendX).attr('y', legendY + legendH + 10)
    .attr('font-size', 9).attr('fill', theme.textMuted).text('−1')
  svg.append('text').attr('x', legendX + legendW / 2).attr('y', legendY + legendH + 10)
    .attr('text-anchor', 'middle').attr('font-size', 9).attr('fill', theme.textMuted).text('0')
  svg.append('text').attr('x', legendX + legendW).attr('y', legendY + legendH + 10)
    .attr('text-anchor', 'end').attr('font-size', 9).attr('fill', theme.textMuted).text('+1')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}
