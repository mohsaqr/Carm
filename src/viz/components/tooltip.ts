/**
 * Hover tooltip component.
 * Shows a floating tooltip near the cursor with formatted content.
 */

import type { CarmTheme } from '../themes/default.js'
import { DEFAULT_THEME } from '../themes/default.js'

let tooltipEl: HTMLDivElement | null = null

function ensureTooltip(): HTMLDivElement {
  if (!tooltipEl) {
    tooltipEl = document.createElement('div')
    tooltipEl.id = 'carm-tooltip'
    tooltipEl.style.cssText = `
      position: fixed;
      pointer-events: none;
      z-index: 9999;
      padding: 8px 12px;
      border-radius: 6px;
      font-size: 12px;
      line-height: 1.5;
      max-width: 240px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
      opacity: 0;
      transition: opacity 0.15s ease;
    `
    document.body.appendChild(tooltipEl)
  }
  return tooltipEl
}

export function showTooltip(
  event: MouseEvent,
  content: string,
  theme: CarmTheme = DEFAULT_THEME
): void {
  const el = ensureTooltip()
  el.innerHTML = content
  el.style.background = theme.surface
  el.style.color = theme.text
  el.style.border = `1px solid ${theme.gridLine}`
  el.style.fontFamily = theme.fontFamily

  const x = event.clientX + 14
  const y = event.clientY - 28
  const viewW = window.innerWidth
  const elW = 240

  el.style.left = `${Math.min(x, viewW - elW - 8)}px`
  el.style.top = `${Math.max(y, 8)}px`
  el.style.opacity = '1'
}

export function hideTooltip(): void {
  if (tooltipEl) tooltipEl.style.opacity = '0'
}

export function formatTooltipRow(label: string, value: string | number): string {
  return `<div style="display:flex;justify-content:space-between;gap:12px">
    <span style="opacity:0.7">${label}</span>
    <strong>${typeof value === 'number' ? value.toFixed(3) : value}</strong>
  </div>`
}
