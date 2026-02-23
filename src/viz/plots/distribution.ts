/**
 * Interactive distribution explorer.
 * Shows PDF/CDF for normal, t, chi-square, F, binomial, Poisson distributions.
 * Uses jStat for accurate distribution functions.
 */

import type * as D3 from 'd3'
import { DEFAULT_THEME, applyTheme, getColor } from '../themes/default.js'
import type { CarmTheme } from '../themes/default.js'
import { addSubtitle, addCaption } from '../components/annotations.js'
import { normalCDF } from '../../core/math.js'

export type DistributionName = 'normal' | 't' | 'chi-square' | 'F' | 'uniform' | 'exponential'

export interface DistributionParams {
  readonly distribution: DistributionName
  readonly params: Readonly<Record<string, number>>  // e.g. { mean: 0, sd: 1 }
  readonly showPDF?: boolean
  readonly showCDF?: boolean
  readonly highlightX?: number  // shade area to the left of this x
}

export interface DistributionConfig {
  readonly title?: string
  readonly caption?: string
  readonly width?: number
  readonly height?: number
  readonly theme?: CarmTheme
}

export function renderDistribution(
  container: HTMLElement,
  params: DistributionParams,
  config: DistributionConfig = {}
): void {
  import('d3').then(d3 => renderDistD3(d3, container, params, config))
}

function renderDistD3(
  d3: typeof D3,
  container: HTMLElement,
  params: DistributionParams,
  config: DistributionConfig
): void {
  const theme = config.theme ?? DEFAULT_THEME
  const W = config.width ?? 550, H = config.height ?? 350
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft }
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom

  container.innerHTML = ''
  applyTheme(container, theme)

  // Get PDF and x domain
  const { xMin, xMax, pdf, cdf } = getDistributionFunctions(params)
  const nPoints = 200
  const xs = Array.from({ length: nPoints }, (_, i) => xMin + (xMax - xMin) * i / (nPoints - 1))
  const pdfs = xs.map(x => pdf(x))
  const cdfs = xs.map(x => cdf(x))

  const subtitle = formatDistributionParams(params)
  const svg = d3.select(container).append('svg').attr('width', W).attr('height', H).style('background', theme.background)
  addSubtitle(svg, config.title ?? `${params.distribution} distribution`, subtitle, W, theme)

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`)

  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width])
  const maxPDF = Math.max(...pdfs.filter(isFinite))
  const yScale = d3.scaleLinear().domain([0, maxPDF * 1.1]).range([height, 0])

  // Grid
  g.selectAll('.grid').data(yScale.ticks(5)).join('line')
    .attr('x1', 0).attr('x2', width).attr('y1', d => yScale(d)).attr('y2', d => yScale(d))
    .attr('stroke', theme.gridLine).attr('stroke-width', 1)

  const lineData = xs.map((x, i) => ({ x, y: pdfs[i] ?? 0 })).filter(d => isFinite(d.y))

  // Highlight area
  if (params.highlightX !== undefined) {
    const hx = params.highlightX
    const shadedData = lineData.filter(d => d.x <= hx)
    const area = d3.area<typeof shadedData[0]>()
      .x(d => xScale(d.x)).y0(height).y1(d => yScale(d.y))
    g.append('path').datum(shadedData).attr('d', area)
      .attr('fill', getColor(0, theme)).attr('opacity', 0.25)

    // P-value label
    const p = cdf(hx)
    g.append('text').attr('x', xScale(hx)).attr('y', 20)
      .attr('text-anchor', 'middle').attr('font-family', theme.fontFamilyMono)
      .attr('font-size', theme.fontSizeSmall).attr('fill', theme.textAnnotation)
      .text(`P(X ≤ ${hx.toFixed(2)}) = ${p.toFixed(4)}`)
  }

  // PDF curve
  if (params.showPDF !== false) {
    const line = d3.line<typeof lineData[0]>().x(d => xScale(d.x)).y(d => yScale(d.y)).curve(d3.curveBasis)
    g.append('path').datum(lineData).attr('d', line)
      .attr('fill', 'none').attr('stroke', getColor(0, theme)).attr('stroke-width', 2.5)
  }

  // CDF curve (secondary y-axis, mapped to same pixel space [0, maxPDF])
  if (params.showCDF === true) {
    const cdfData = xs.map((x, i) => ({ x, y: (cdfs[i] ?? 0) * maxPDF })).filter(d => isFinite(d.y))
    const line = d3.line<typeof cdfData[0]>().x(d => xScale(d.x)).y(d => yScale(d.y)).curve(d3.curveBasis)
    g.append('path').datum(cdfData).attr('d', line)
      .attr('fill', 'none').attr('stroke', getColor(1, theme))
      .attr('stroke-width', 2).attr('stroke-dasharray', '6,3')
  }

  // Axes
  g.append('g').attr('transform', `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('fill', theme.text).attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize)

  g.append('text').attr('x', width / 2).attr('y', height + 44).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('x')
  g.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', -48).attr('text-anchor', 'middle')
    .attr('font-family', theme.fontFamily).attr('font-size', theme.fontSize).attr('fill', theme.text).text('Density')

  if (config.caption) addCaption(svg, config.caption, W, H, theme)
}

function getDistributionFunctions(params: DistributionParams): {
  xMin: number; xMax: number
  pdf: (x: number) => number
  cdf: (x: number) => number
} {
  const p = params.params
  switch (params.distribution) {
    case 'normal': {
      const mu = p['mean'] ?? 0, sigma = p['sd'] ?? 1
      return {
        xMin: mu - 4 * sigma, xMax: mu + 4 * sigma,
        pdf: x => Math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * Math.sqrt(2 * Math.PI)),
        cdf: x => normalCDF((x - mu) / sigma),
      }
    }
    case 't': {
      const df = p['df'] ?? 5
      return {
        xMin: -5, xMax: 5,
        pdf: x => {
          const c = Math.exp(lgamma((df + 1) / 2) - lgamma(df / 2)) / Math.sqrt(df * Math.PI)
          return c * Math.pow(1 + x * x / df, -(df + 1) / 2)
        },
        cdf: x => tDistCDFLocal(x, df),
      }
    }
    case 'chi-square': {
      const df = p['df'] ?? 5
      return {
        xMin: 0, xMax: df + 5 * Math.sqrt(2 * df),
        pdf: x => x <= 0 ? 0 : Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - (df / 2) * Math.log(2) - lgamma(df / 2)),
        cdf: x => chiSqCDFLocal(x, df),
      }
    }
    case 'uniform': {
      const a = p['min'] ?? 0, b = p['max'] ?? 1
      return {
        xMin: a - 0.2 * (b - a), xMax: b + 0.2 * (b - a),
        pdf: x => (x >= a && x <= b) ? 1 / (b - a) : 0,
        cdf: x => x < a ? 0 : x > b ? 1 : (x - a) / (b - a),
      }
    }
    case 'exponential': {
      const rate = p['rate'] ?? 1
      return {
        xMin: 0, xMax: 5 / rate,
        pdf: x => x < 0 ? 0 : rate * Math.exp(-rate * x),
        cdf: x => x < 0 ? 0 : 1 - Math.exp(-rate * x),
      }
    }
    case 'F': {
      const df1 = p['df1'] ?? 2, df2 = p['df2'] ?? 10
      return {
        xMin: 0, xMax: 6,
        pdf: x => {
          if (x <= 0) return 0
          const num = Math.pow(df1 * x, df1) * Math.pow(df2, df2)
          const den = Math.pow(df1 * x + df2, df1 + df2)
          return Math.sqrt(num / den) / (x * Math.exp(logBetaLocal(df1 / 2, df2 / 2)))
        },
        cdf: x => {
          if (x <= 0) return 0
          return incompleteBetaLocal(df1 * x / (df1 * x + df2), df1 / 2, df2 / 2)
        },
      }
    }
    default:
      throw new Error(`Unknown distribution: ${params.distribution}`)
  }
}

function formatDistributionParams(params: DistributionParams): string {
  const entries = Object.entries(params.params).map(([k, v]) => `${k} = ${v}`).join(', ')
  return `${params.distribution}(${entries})`
}

// Minimal standalone implementations to avoid circular imports
function lgamma(z: number): number {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
  const x = z - 1
  let sum = c[0]!
  for (let i = 1; i < 9; i++) sum += (c[i] ?? 0) / (x + i)
  const t = x + 7.5
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum)
}

function tDistCDFLocal(t: number, df: number): number {
  const x = df / (df + t * t)
  const p = incompleteBetaLocal(x, df / 2, 0.5) / 2
  return t >= 0 ? 1 - p : p
}

function chiSqCDFLocal(x: number, df: number): number {
  if (x <= 0) return 0
  return incompleteGammaLocal(df / 2, x / 2)
}

function incompleteGammaLocal(a: number, x: number): number {
  if (x < a + 1) {
    let term = 1 / a, sum = term
    for (let n = 1; n < 200; n++) {
      term *= x / (a + n); sum += term
      if (Math.abs(term) < Math.abs(sum) * 3e-7) break
    }
    return sum * Math.exp(-x + a * Math.log(x) - lgamma(a))
  }
  let f = x + 1 - a; const FPMIN = 1e-30
  if (Math.abs(f) < FPMIN) f = FPMIN
  let C = f, D = 0
  for (let i = 1; i <= 200; i++) {
    const an = -i * (i - a), bn = x + 2 * i + 1 - a
    D = bn + an * D; if (Math.abs(D) < FPMIN) D = FPMIN
    C = bn + an / C; if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D; const delta = C * D; f *= delta
    if (Math.abs(delta - 1) < 3e-7) break
  }
  return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) / f
}

function logBetaLocal(a: number, b: number): number {
  return lgamma(a) + lgamma(b) - lgamma(a + b)
}

function incompleteBetaLocal(x: number, a: number, b: number): number {
  if (x <= 0) return 0; if (x >= 1) return 1
  if (x > (a + 1) / (a + b + 2)) return 1 - incompleteBetaLocal(1 - x, b, a)
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - logBetaLocal(a, b)) / a
  const FPMIN = 1e-30
  let f = 1, C = 1, D = 1 - (a + b) * x / (a + 1)
  if (Math.abs(D) < FPMIN) D = FPMIN; D = 1 / D; f = D
  for (let m = 1; m <= 200; m++) {
    let num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
    D = 1 + num * D; if (Math.abs(D) < FPMIN) D = FPMIN
    C = 1 + num / C; if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D; f *= C * D
    num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
    D = 1 + num * D; if (Math.abs(D) < FPMIN) D = FPMIN
    C = 1 + num / C; if (Math.abs(C) < FPMIN) C = FPMIN
    D = 1 / D; const delta = C * D; f *= delta
    if (Math.abs(delta - 1) < 3e-7) break
  }
  return front * f
}
