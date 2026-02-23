/**
 * Funnel chart: trapezoidal segments stacked vertically.
 * Width proportional to value. Labels on left, value + % drop on right.
 */

import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { showTooltip, hideTooltip, formatTooltipRow } from '../components/tooltip.js'
import type { StatResult } from '../../core/types.js'

export interface FunnelConfig {
  readonly title?: string
  readonly xLabel?: string
  readonly yLabel?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export interface FunnelStage {
  readonly label: string
  readonly value: number
}

export interface FunnelData {
  readonly stages: readonly FunnelStage[]
  readonly testResult?: StatResult
}

export function renderFunnel(
  container: HTMLElement,
  data: FunnelData,
  config: FunnelConfig = {}
): void {
  import('d3').then(d3 => renderFunnelD3(d3, container, data, config))
}

function renderFunnelD3(
  _d3: unknown,
  container: HTMLElement,
  data: FunnelData,
  config: FunnelConfig
): void {
  if (data.stages.length === 0) return

  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? Math.max(container.clientWidth || 600, 400)
  const H = config.height ?? 480
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right
  const height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  // Build SVG manually without D3 (funnel uses SVG path primitives only)
  const svgNS = 'http://www.w3.org/2000/svg'
  const svgEl = document.createElementNS(svgNS, 'svg')
  svgEl.setAttribute('width', String(W))
  svgEl.setAttribute('height', String(H))
  svgEl.style.background = theme.background
  container.appendChild(svgEl)

  // Title
  const titleEl = document.createElementNS(svgNS, 'text')
  titleEl.setAttribute('x', String(W / 2))
  titleEl.setAttribute('y', '20')
  titleEl.setAttribute('text-anchor', 'middle')
  titleEl.setAttribute('font-family', theme.fontFamily)
  titleEl.setAttribute('font-size', String(theme.fontSizeTitle))
  titleEl.setAttribute('font-weight', '600')
  titleEl.setAttribute('fill', theme.text)
  titleEl.textContent = config.title ?? 'Funnel Chart'
  svgEl.appendChild(titleEl)

  const subtitleEl = document.createElementNS(svgNS, 'text')
  subtitleEl.setAttribute('x', String(W / 2))
  subtitleEl.setAttribute('y', '38')
  subtitleEl.setAttribute('text-anchor', 'middle')
  subtitleEl.setAttribute('font-family', theme.fontFamilyMono ?? theme.fontFamily)
  subtitleEl.setAttribute('font-size', String(theme.fontSizeSmall))
  subtitleEl.setAttribute('fill', theme.textMuted)
  subtitleEl.textContent = data.testResult?.formatted ?? ''
  svgEl.appendChild(subtitleEl)

  const maxValue = Math.max(...data.stages.map(s => s.value))
  const n = data.stages.length
  const segH = height / n
  const gap = 2
  const labelColW = 110
  const valueColW = 110
  const funnelW = width - labelColW - valueColW
  const funnelOffsetX = margin.left + labelColW
  const offsetY = margin.top

  data.stages.forEach((stage, i) => {
    const color = getColor(i, theme)
    const topW = (stage.value / maxValue) * funnelW
    const nextStage = i + 1 < data.stages.length ? data.stages[i + 1] : undefined
    const bottomW = nextStage != null ? (nextStage.value / maxValue) * funnelW : topW * 0.6

    const topY = offsetY + i * segH
    const bottomY = offsetY + (i + 1) * segH - gap

    const topLeft = funnelOffsetX + (funnelW - topW) / 2
    const topRight = funnelOffsetX + (funnelW + topW) / 2
    const botLeft = funnelOffsetX + (funnelW - bottomW) / 2
    const botRight = funnelOffsetX + (funnelW + bottomW) / 2

    const pathD = `M ${topLeft},${topY} L ${topRight},${topY} L ${botRight},${bottomY} L ${botLeft},${bottomY} Z`
    const pathEl = document.createElementNS(svgNS, 'path')
    pathEl.setAttribute('d', pathD)
    pathEl.setAttribute('fill', color)
    pathEl.setAttribute('opacity', '0.85')
    pathEl.addEventListener('mouseover', (event: MouseEvent) => {
      const pct = ((stage.value / maxValue) * 100).toFixed(1)
      const prevStage = i > 0 ? data.stages[i - 1] : undefined
      const drop = prevStage != null
        ? (((prevStage.value - stage.value) / prevStage.value) * 100).toFixed(1)
        : '—'
      showTooltip(event, [
        formatTooltipRow('Stage', stage.label),
        formatTooltipRow('Value', stage.value.toFixed(0)),
        formatTooltipRow('% of Max', pct + '%'),
        formatTooltipRow('Drop from prev', drop === '—' ? drop : drop + '%')
      ].join(''), theme)
    })
    pathEl.addEventListener('mouseout', hideTooltip)
    svgEl.appendChild(pathEl)

    const midY = topY + segH / 2

    // Label on left
    const labelEl = document.createElementNS(svgNS, 'text')
    labelEl.setAttribute('x', String(margin.left + labelColW - 8))
    labelEl.setAttribute('y', String(midY))
    labelEl.setAttribute('text-anchor', 'end')
    labelEl.setAttribute('dominant-baseline', 'middle')
    labelEl.setAttribute('font-family', theme.fontFamily)
    labelEl.setAttribute('font-size', String(theme.fontSize))
    labelEl.setAttribute('fill', theme.text)
    labelEl.textContent = stage.label
    svgEl.appendChild(labelEl)

    // Value + % on right
    const pctOfMax = ((stage.value / maxValue) * 100).toFixed(0)
    const prevStage = i > 0 ? data.stages[i - 1] : undefined
    const dropStr = prevStage != null
      ? `\u2212${(((prevStage.value - stage.value) / prevStage.value) * 100).toFixed(0)}%`
      : ''

    const valEl = document.createElementNS(svgNS, 'text')
    valEl.setAttribute('x', String(funnelOffsetX + funnelW + 8))
    valEl.setAttribute('y', String(midY - (dropStr ? 7 : 0)))
    valEl.setAttribute('dominant-baseline', 'middle')
    valEl.setAttribute('font-family', theme.fontFamilyMono ?? theme.fontFamily)
    valEl.setAttribute('font-size', String(theme.fontSize))
    valEl.setAttribute('fill', theme.text)
    valEl.textContent = `${stage.value.toLocaleString()} (${pctOfMax}%)`
    svgEl.appendChild(valEl)

    if (dropStr) {
      const dropEl = document.createElementNS(svgNS, 'text')
      dropEl.setAttribute('x', String(funnelOffsetX + funnelW + 8))
      dropEl.setAttribute('y', String(midY + 9))
      dropEl.setAttribute('dominant-baseline', 'middle')
      dropEl.setAttribute('font-family', theme.fontFamily)
      dropEl.setAttribute('font-size', String(theme.fontSizeSmall))
      dropEl.setAttribute('fill', '#D55E00')
      dropEl.textContent = dropStr
      svgEl.appendChild(dropEl)
    }
  })

  if (config.caption) {
    const capEl = document.createElementNS(svgNS, 'text')
    capEl.setAttribute('x', String(W / 2))
    capEl.setAttribute('y', String(H - 6))
    capEl.setAttribute('text-anchor', 'middle')
    capEl.setAttribute('font-family', theme.fontFamily)
    capEl.setAttribute('font-size', String(theme.fontSizeSmall - 1))
    capEl.setAttribute('fill', theme.textMuted)
    capEl.style.fontStyle = 'italic'
    capEl.textContent = config.caption
    svgEl.appendChild(capEl)
  }
}
