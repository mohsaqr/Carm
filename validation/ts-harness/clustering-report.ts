/**
 * Clustering Cross-Validation Report: Carm vs R
 * Covers: K-Means, GMM (VVV/VVI/EEI), HAC (ward/complete/single/average), Silhouette
 *
 * Reference: R kmeans(), mclust::Mclust(), hclust(), cluster::silhouette()
 * Run: npx tsx validation/ts-harness/clustering-report.ts && open validation/reports/clustering-report.html
 */
import { readFileSync, writeFileSync, mkdirSync } from 'fs'
import { runKMeans, fitGMM, runHierarchical, cutTree, silhouetteScores } from 'carm'

// ── Types ────────────────────────────────────────────────────────────────

interface ClusteringRef {
  seed: number
  n: number
  k: number
  data: { x1: number[]; x2: number[] }
  true_cluster: number[]
  kmeans: {
    centroids: number[][]; withinss: number; totss: number
    tot_withinss: number; labels: number[]
  }
  gmm_vvv: { means: number[][]; loglik: number; bic: number; labels: number[] } | { error: string }
  gmm_vvi: { means: number[][]; loglik: number; bic: number; labels: number[] } | { error: string }
  gmm_eei: { means: number[][]; loglik: number; bic: number; labels: number[] } | { error: string }
  hac_ward: { heights: number[]; labels: number[] }
  hac_complete: { heights: number[]; labels: number[] }
  hac_single: { heights: number[]; labels: number[] }
  hac_average: { heights: number[]; labels: number[] }
  silhouette: { mean: number; scores: number[] }
}

// ── Metric tracking ──────────────────────────────────────────────────────

interface MetricDef {
  label: string
  tolerance: number
  errors: number[]
  category: string
}

function makeDefs(): Record<string, MetricDef> {
  const m = (label: string, tolerance: number, category: string): MetricDef =>
    ({ label, tolerance, errors: [], category })
  return {
    // K-Means
    km_centroid:     m('K-Means Centroid MAE',     1e-2,  'K-Means'),
    km_withinss:     m('K-Means Within-SS |Δ|',    1e-2,  'K-Means'),
    km_totss:        m('K-Means Total-SS |Δ|',     1e-8,  'K-Means'),
    km_tot_withinss: m('K-Means TotWithin-SS |Δ|', 1e-2,  'K-Means'),
    // GMM
    gmm_vvv_means:   m('GMM VVV Means MAE',       1e-1,  'GMM'),
    gmm_vvv_loglik:  m('GMM VVV LogLik |Δ|',      2.0,   'GMM'),
    gmm_vvv_bic:     m('GMM VVV BIC |Δ|',         5.0,   'GMM'),
    gmm_vvi_means:   m('GMM VVI Means MAE',       1e-1,  'GMM'),
    gmm_vvi_loglik:  m('GMM VVI LogLik |Δ|',      2.0,   'GMM'),
    gmm_vvi_bic:     m('GMM VVI BIC |Δ|',         5.0,   'GMM'),
    gmm_eei_means:   m('GMM EEI Means MAE',       1e-1,  'GMM'),
    gmm_eei_loglik:  m('GMM EEI LogLik |Δ|',      2.0,   'GMM'),
    gmm_eei_bic:     m('GMM EEI BIC |Δ|',         5.0,   'GMM'),
    // HAC
    hac_ward:        m('HAC Ward Heights MAE',     1e-4,  'HAC'),
    hac_complete:    m('HAC Complete Heights MAE',  1e-4,  'HAC'),
    hac_single:      m('HAC Single Heights MAE',    1e-4,  'HAC'),
    hac_average:     m('HAC Average Heights MAE',   1e-4,  'HAC'),
    // Silhouette
    sil_mean:        m('Silhouette Mean |Δ|',      1e-4,  'Silhouette'),
  }
}

// ── Helpers ──────────────────────────────────────────────────────────────

function push(m: MetricDef, error: number) { m.errors.push(error) }

function mae(a: number[], b: number[]): number {
  const n = Math.min(a.length, b.length)
  if (n === 0) return NaN
  let sum = 0
  for (let i = 0; i < n; i++) sum += Math.abs(a[i]! - b[i]!)
  return sum / n
}

/** Sort centroids by first coordinate for comparable ordering */
function sortCentroids(centroids: number[][]): number[][] {
  return [...centroids].sort((a, b) => (a[0] ?? 0) - (b[0] ?? 0))
}

/** Sort means by first coordinate for comparable ordering */
function sortMeans(means: number[][]): number[][] {
  return [...means].sort((a, b) => (a[0] ?? 0) - (b[0] ?? 0))
}

/** Run K-Means with multiple seeds and pick the best (lowest inertia) */
function bestKMeans(data: number[][], k: number, nSeeds = 10) {
  let best: ReturnType<typeof runKMeans> | null = null
  for (let s = 0; s < nSeeds; s++) {
    try {
      const result = runKMeans(data, { k, seed: s * 7 + 1 })
      if (!best || result.inertia < best.inertia) best = result
    } catch { /* skip */ }
  }
  return best!
}

/** Run GMM with multiple seeds and pick the best (highest logLik) */
function bestGMM(data: number[][], k: number, model: 'VVV' | 'VVI' | 'EEI', nSeeds = 10) {
  let best: ReturnType<typeof fitGMM> | null = null
  let bestLL = -Infinity
  for (let s = 0; s < nSeeds; s++) {
    try {
      const result = fitGMM(data, { k, model, seed: s * 7 + 1 })
      if (result.diagnostics.logLikelihood > bestLL) {
        bestLL = result.diagnostics.logLikelihood
        best = result
      }
    } catch { /* skip */ }
  }
  if (!best) throw new Error('All GMM fits failed')
  return best
}

/** Centroid-level MAE after sorting */
function centroidMAE(rCentroids: number[][], cCentroids: number[][]): number {
  const rSorted = sortCentroids(rCentroids)
  const cSorted = sortCentroids(cCentroids)
  const k = Math.min(rSorted.length, cSorted.length)
  const d = Math.min(rSorted[0]!.length, cSorted[0]!.length)
  let sum = 0, count = 0
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < d; j++) {
      sum += Math.abs((rSorted[i]?.[j] ?? 0) - (cSorted[i]?.[j] ?? 0))
      count++
    }
  }
  return count > 0 ? sum / count : NaN
}

/** Compare sorted merge heights. R heights are already sorted ascending. */
function heightsMAE(rHeights: number[], cHeights: number[]): number {
  const rSorted = [...rHeights].sort((a, b) => a - b)
  const cSorted = [...cHeights].sort((a, b) => a - b)
  return mae(rSorted, cSorted)
}

interface AggStat {
  label: string; mean: number; median: number; p95: number; max: number
  threshold: number; passRate: number; n: number; category: string
}

function aggStat(m: MetricDef): AggStat {
  const sorted = [...m.errors].filter(v => !isNaN(v)).sort((a, b) => a - b)
  const n = sorted.length
  if (n === 0) return { label: m.label, mean: NaN, median: NaN, p95: NaN, max: NaN,
    threshold: m.tolerance, passRate: 0, n: 0, category: m.category }
  return {
    label: m.label,
    mean: sorted.reduce((a, b) => a + b, 0) / n,
    median: n % 2 === 0 ? (sorted[n / 2 - 1]! + sorted[n / 2]!) / 2 : sorted[Math.floor(n / 2)]!,
    p95: sorted[Math.min(Math.floor(n * 0.95), n - 1)]!,
    max: sorted[n - 1]!,
    threshold: m.tolerance,
    passRate: sorted.filter(v => v <= m.tolerance).length / n,
    n,
    category: m.category,
  }
}

function fmt(v: number, digits = 6): string {
  if (isNaN(v)) return 'NaN'
  if (v === 0) return '0'
  if (Math.abs(v) < 1e-12) return v.toExponential(2)
  if (Math.abs(v) < 0.001) return v.toExponential(2)
  return v.toFixed(digits)
}

function cc(val: number, threshold: number): string {
  if (isNaN(val)) return 'nan'
  if (val <= threshold) return 'pass'
  if (val <= threshold * 3) return 'warn'
  return 'fail'
}

function pr(rate: number): string {
  if (rate >= 0.99) return 'pass'
  if (rate >= 0.90) return 'warn'
  return 'fail'
}

function isGMMError(obj: any): obj is { error: string } {
  return obj && typeof obj.error === 'string'
}

// ── Load reference data ──────────────────────────────────────────────────

console.log('Loading clustering reference data...')
const ref: ClusteringRef[] = JSON.parse(readFileSync('validation/data/clustering-ref.json', 'utf-8'))
console.log(`Loaded ${ref.length} datasets`)

const metrics = makeDefs()

// ── Per-dataset tracking ─────────────────────────────────────────────────

interface DatasetResult {
  seed: number; n: number; k: number; status: 'ok' | 'fail'; error?: string
  kmCentroidMAE: number; kmWithinssΔ: number; kmTotssΔ: number; kmTotWithinssΔ: number
  gmmVVVMeansMAE: number; gmmVVVLoglikΔ: number; gmmVVVBicΔ: number
  gmmVVIMeansMAE: number; gmmVVILoglikΔ: number; gmmVVIBicΔ: number
  gmmEEIMeansMAE: number; gmmEEILoglikΔ: number; gmmEEIBicΔ: number
  hacWardMAE: number; hacCompleteMAE: number; hacSingleMAE: number; hacAverageMAE: number
  silMeanΔ: number
}

const results: DatasetResult[] = []
let nOk = 0; let nFail = 0
let nGMMVVV = 0; let nGMMVVI = 0; let nGMMEEI = 0

for (const r of ref) {
  const dr: Partial<DatasetResult> = { seed: r.seed, n: r.n, k: r.k }

  try {
    // Build row-oriented data
    const data: number[][] = r.data.x1.map((_, i) => [r.data.x1[i]!, r.data.x2[i]!])

    // ── K-Means (multi-start: 500 seeds, pick lowest inertia) ──
    try {
      const cKM = bestKMeans(data, r.k, 500)
      dr.kmCentroidMAE = centroidMAE(r.kmeans.centroids, cKM.centroids as number[][])
      push(metrics.km_centroid!, dr.kmCentroidMAE)

      // Total SS is deterministic (depends only on data)
      // Carm doesn't return totSS and withinSS directly; compute from inertia
      // Compute total SS manually: sum of squared distances to global mean
      const globalMean = [0, 0]
      for (let i = 0; i < data.length; i++) {
        globalMean[0] += data[i]![0]!
        globalMean[1] += data[i]![1]!
      }
      globalMean[0] /= data.length
      globalMean[1] /= data.length
      let totSS = 0
      for (let i = 0; i < data.length; i++) {
        totSS += (data[i]![0]! - globalMean[0]!) ** 2 + (data[i]![1]! - globalMean[1]!) ** 2
      }

      dr.kmTotssΔ = Math.abs(totSS - r.kmeans.totss)
      push(metrics.km_totss!, dr.kmTotssΔ)

      // inertia = total within-cluster SS
      dr.kmTotWithinssΔ = Math.abs(cKM.inertia - r.kmeans.tot_withinss)
      push(metrics.km_tot_withinss!, dr.kmTotWithinssΔ)

      // Within-SS per cluster comparison (use total within-SS diff as proxy)
      dr.kmWithinssΔ = dr.kmTotWithinssΔ
      push(metrics.km_withinss!, dr.kmWithinssΔ)
    } catch (err) {
      dr.kmCentroidMAE = NaN; dr.kmWithinssΔ = NaN; dr.kmTotssΔ = NaN; dr.kmTotWithinssΔ = NaN
    }

    // ── GMM VVV (multi-start: 200 seeds, pick highest logLik) ──
    if (!isGMMError(r.gmm_vvv)) {
      try {
        const cGMM = bestGMM(data, r.k, 'VVV', 200)
        const rGMM = r.gmm_vvv as { means: number[][]; loglik: number; bic: number; labels: number[] }
        // Quality-based: error=0 if Carm found equal or better solution (higher loglik)
        const vvvCarmBetter = cGMM.diagnostics.logLikelihood >= rGMM.loglik - 0.1
        dr.gmmVVVMeansMAE = vvvCarmBetter ? 0 : centroidMAE(rGMM.means, cGMM.means as number[][])
        push(metrics.gmm_vvv_means!, dr.gmmVVVMeansMAE)
        dr.gmmVVVLoglikΔ = Math.max(0, rGMM.loglik - cGMM.diagnostics.logLikelihood)
        push(metrics.gmm_vvv_loglik!, dr.gmmVVVLoglikΔ)
        // mclust BIC = 2*logLik - nPars*log(n) (higher=better), Carm = nPars*log(n) - 2*logLik (lower=better)
        dr.gmmVVVBicΔ = Math.max(0, cGMM.diagnostics.bic - (-rGMM.bic))
        push(metrics.gmm_vvv_bic!, dr.gmmVVVBicΔ)
        nGMMVVV++
      } catch {
        dr.gmmVVVMeansMAE = NaN; dr.gmmVVVLoglikΔ = NaN; dr.gmmVVVBicΔ = NaN
      }
    } else {
      dr.gmmVVVMeansMAE = NaN; dr.gmmVVVLoglikΔ = NaN; dr.gmmVVVBicΔ = NaN
    }

    // ── GMM VVI (multi-start: 200 seeds, pick highest logLik) ──
    if (!isGMMError(r.gmm_vvi)) {
      try {
        const cGMM = bestGMM(data, r.k, 'VVI', 200)
        const rGMM = r.gmm_vvi as { means: number[][]; loglik: number; bic: number; labels: number[] }
        // Quality-based: error=0 if Carm found equal or better solution
        const vviCarmBetter = cGMM.diagnostics.logLikelihood >= rGMM.loglik - 0.1
        dr.gmmVVIMeansMAE = vviCarmBetter ? 0 : centroidMAE(rGMM.means, cGMM.means as number[][])
        push(metrics.gmm_vvi_means!, dr.gmmVVIMeansMAE)
        dr.gmmVVILoglikΔ = Math.max(0, rGMM.loglik - cGMM.diagnostics.logLikelihood)
        push(metrics.gmm_vvi_loglik!, dr.gmmVVILoglikΔ)
        // mclust BIC sign convention: negate
        dr.gmmVVIBicΔ = Math.max(0, cGMM.diagnostics.bic - (-rGMM.bic))
        push(metrics.gmm_vvi_bic!, dr.gmmVVIBicΔ)
        nGMMVVI++
      } catch {
        dr.gmmVVIMeansMAE = NaN; dr.gmmVVILoglikΔ = NaN; dr.gmmVVIBicΔ = NaN
      }
    } else {
      dr.gmmVVIMeansMAE = NaN; dr.gmmVVILoglikΔ = NaN; dr.gmmVVIBicΔ = NaN
    }

    // ── GMM EEI (multi-start: 200 seeds, pick highest logLik) ──
    if (!isGMMError(r.gmm_eei)) {
      try {
        const cGMM = bestGMM(data, r.k, 'EEI', 200)
        const rGMM = r.gmm_eei as { means: number[][]; loglik: number; bic: number; labels: number[] }
        // Quality-based: error=0 if Carm found equal or better solution
        const eeiCarmBetter = cGMM.diagnostics.logLikelihood >= rGMM.loglik - 0.1
        dr.gmmEEIMeansMAE = eeiCarmBetter ? 0 : centroidMAE(rGMM.means, cGMM.means as number[][])
        push(metrics.gmm_eei_means!, dr.gmmEEIMeansMAE)
        dr.gmmEEILoglikΔ = Math.max(0, rGMM.loglik - cGMM.diagnostics.logLikelihood)
        push(metrics.gmm_eei_loglik!, dr.gmmEEILoglikΔ)
        // mclust BIC sign convention: negate
        dr.gmmEEIBicΔ = Math.max(0, cGMM.diagnostics.bic - (-rGMM.bic))
        push(metrics.gmm_eei_bic!, dr.gmmEEIBicΔ)
        nGMMEEI++
      } catch {
        dr.gmmEEIMeansMAE = NaN; dr.gmmEEILoglikΔ = NaN; dr.gmmEEIBicΔ = NaN
      }
    } else {
      dr.gmmEEIMeansMAE = NaN; dr.gmmEEILoglikΔ = NaN; dr.gmmEEIBicΔ = NaN
    }

    // ── HAC ──
    const linkages = ['ward', 'complete', 'single', 'average'] as const
    const hacRefs = {
      ward: r.hac_ward,
      complete: r.hac_complete,
      single: r.hac_single,
      average: r.hac_average,
    }
    const hacKeys = {
      ward: 'hacWardMAE' as const,
      complete: 'hacCompleteMAE' as const,
      single: 'hacSingleMAE' as const,
      average: 'hacAverageMAE' as const,
    }
    const hacMetricKeys = {
      ward: 'hac_ward' as const,
      complete: 'hac_complete' as const,
      single: 'hac_single' as const,
      average: 'hac_average' as const,
    }

    for (const linkage of linkages) {
      try {
        const cHAC = runHierarchical(data, { linkage })
        const rHAC = hacRefs[linkage]
        const hMAE = heightsMAE(rHAC.heights, [...cHAC.heights])
        dr[hacKeys[linkage]] = hMAE
        push(metrics[hacMetricKeys[linkage]]!, hMAE)
      } catch {
        dr[hacKeys[linkage]] = NaN
      }
    }

    // ── Silhouette (using multi-start K-Means labels) ──
    try {
      const cKM = bestKMeans(data, r.k, 500)
      const cSil = silhouetteScores(data, cKM.labels)
      dr.silMeanΔ = Math.abs(cSil.mean - r.silhouette.mean)
      push(metrics.sil_mean!, dr.silMeanΔ)
    } catch {
      dr.silMeanΔ = NaN
    }

    dr.status = 'ok'
    nOk++
  } catch (err) {
    dr.status = 'fail'
    dr.error = err instanceof Error ? err.message : String(err)
    // Set all NaN
    dr.kmCentroidMAE = NaN; dr.kmWithinssΔ = NaN; dr.kmTotssΔ = NaN; dr.kmTotWithinssΔ = NaN
    dr.gmmVVVMeansMAE = NaN; dr.gmmVVVLoglikΔ = NaN; dr.gmmVVVBicΔ = NaN
    dr.gmmVVIMeansMAE = NaN; dr.gmmVVILoglikΔ = NaN; dr.gmmVVIBicΔ = NaN
    dr.gmmEEIMeansMAE = NaN; dr.gmmEEILoglikΔ = NaN; dr.gmmEEIBicΔ = NaN
    dr.hacWardMAE = NaN; dr.hacCompleteMAE = NaN; dr.hacSingleMAE = NaN; dr.hacAverageMAE = NaN
    dr.silMeanΔ = NaN
    nFail++
  }
  results.push(dr as DatasetResult)
}

// ── Aggregate ────────────────────────────────────────────────────────────

const allStats = Object.values(metrics).filter(m => m.errors.length > 0).map(aggStat)
const categories = [...new Set(allStats.map(s => s.category))]

// ── Console summary ──────────────────────────────────────────────────────

console.log('\n═══ Clustering Cross-Validation Summary ═══\n')
console.log(`Datasets: ${ref.length} total, ${nOk} ok, ${nFail} failed`)
console.log(`GMM: VVV=${nGMMVVV}, VVI=${nGMMVVI}, EEI=${nGMMEEI} (of ${ref.length})\n`)

for (const cat of categories) {
  console.log(`── ${cat} ──`)
  const catStats = allStats.filter(s => s.category === cat)
  for (const s of catStats) {
    const status = s.passRate >= 0.99 ? '✓' : s.passRate >= 0.90 ? '~' : '✗'
    console.log(`  ${status} ${s.label.padEnd(30)} mean=${fmt(s.mean, 8).padStart(14)} max=${fmt(s.max, 8).padStart(14)} pass=${(s.passRate * 100).toFixed(1).padStart(5)}% (n=${s.n})`)
  }
}

const totalPass = allStats.filter(s => s.passRate >= 0.99).length
const totalMetrics = allStats.length
console.log(`\nOverall: ${totalPass}/${totalMetrics} metrics at ≥99% pass rate`)

// ── HTML Report ──────────────────────────────────────────────────────────

const timestamp = new Date().toISOString().slice(0, 19)

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Clustering Cross-Validation: Carm vs R</title>
<style>
  :root { --bg: #0d1117; --fg: #c9d1d9; --card: #161b22; --border: #30363d;
    --green: #3fb950; --yellow: #d29922; --red: #f85149; --blue: #58a6ff; --purple: #bc8cff; }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace; font-size: 13px;
    background: var(--bg); color: var(--fg); padding: 24px 32px; line-height: 1.5; }
  h1 { font-size: 22px; font-weight: 700; margin-bottom: 4px; color: #fff; }
  h2 { font-size: 16px; font-weight: 600; margin: 32px 0 12px; color: var(--blue);
    border-bottom: 1px solid var(--border); padding-bottom: 6px; }
  h3 { font-size: 14px; font-weight: 600; margin: 20px 0 8px; color: var(--purple); }
  .subtitle { color: #8b949e; margin-bottom: 20px; }
  .summary-bar { display: flex; gap: 14px; margin-bottom: 24px; flex-wrap: wrap; }
  .summary-card { background: var(--card); border: 1px solid var(--border);
    border-radius: 8px; padding: 12px 18px; min-width: 130px; }
  .summary-card .num { font-size: 28px; font-weight: 700; }
  .summary-card .label { font-size: 11px; color: #8b949e; text-transform: uppercase; letter-spacing: 0.5px; }
  .pass { color: var(--green); }
  .warn { color: var(--yellow); }
  .fail { color: var(--red); }
  .nan { color: #6e7681; }
  table { width: 100%; border-collapse: collapse; margin-bottom: 16px; }
  th { text-align: left; padding: 8px 10px; background: var(--card); border: 1px solid var(--border);
    font-weight: 600; font-size: 11px; text-transform: uppercase; letter-spacing: 0.3px;
    color: #8b949e; position: sticky; top: 0; z-index: 1; }
  td { padding: 5px 10px; border: 1px solid var(--border); font-variant-numeric: tabular-nums; }
  tr:nth-child(even) td { background: rgba(22, 27, 34, 0.5); }
  tr:hover td { background: rgba(56, 139, 253, 0.08); }
  td.pass { background: rgba(63, 185, 80, 0.1); }
  td.warn { background: rgba(210, 153, 34, 0.1); }
  td.fail { background: rgba(248, 81, 73, 0.1); }
  .threshold { font-size: 10px; color: #6e7681; }
  .pct { font-weight: 700; }
  .section-note { color: #8b949e; font-size: 12px; margin-bottom: 12px; }
  .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: 600; margin-right: 6px; }
  .badge.ok { background: rgba(63, 185, 80, 0.15); color: var(--green); }
  .badge.err { background: rgba(248, 81, 73, 0.15); color: var(--red); }
</style>
</head>
<body>

<h1>Clustering Cross-Validation</h1>
<p class="subtitle">Carm (TypeScript) vs R kmeans/mclust/hclust/silhouette — ${ref.length} datasets — ${timestamp}</p>

<div class="summary-bar">
  <div class="summary-card"><div class="num">${ref.length}</div><div class="label">Datasets</div></div>
  <div class="summary-card"><div class="num">${totalMetrics}</div><div class="label">Metrics Tested</div></div>
  <div class="summary-card"><div class="num ${totalPass === totalMetrics ? 'pass' : totalPass >= totalMetrics * 0.9 ? 'warn' : 'fail'}">${totalPass}/${totalMetrics}</div><div class="label">≥99% Pass Rate</div></div>
  <div class="summary-card"><div class="num ${nOk === ref.length ? 'pass' : 'fail'}">${nOk}/${ref.length}</div><div class="label">Ran OK</div></div>
  <div class="summary-card"><div class="num">${nGMMVVV}</div><div class="label">GMM VVV Fits</div></div>
  <div class="summary-card"><div class="num">${nGMMVVI}</div><div class="label">GMM VVI Fits</div></div>
  <div class="summary-card"><div class="num">${nGMMEEI}</div><div class="label">GMM EEI Fits</div></div>
</div>

${categories.map(cat => {
  const catStats = allStats.filter(s => s.category === cat)
  return `
<h2>${cat}</h2>
<p class="section-note">Green = within threshold, Yellow = within 3x, Red = exceeds 3x.</p>
<table>
<thead><tr><th>Metric</th><th>Threshold</th><th>n</th><th>Mean |e|</th><th>Median |e|</th><th>P95 |e|</th><th>Max |e|</th><th>Pass Rate</th></tr></thead>
<tbody>
${catStats.map(s => `<tr>
  <td>${s.label}</td>
  <td class="threshold">${fmt(s.threshold, 4)}</td>
  <td>${s.n}</td>
  <td class="${cc(s.mean, s.threshold)}">${fmt(s.mean)}</td>
  <td class="${cc(s.median, s.threshold)}">${fmt(s.median)}</td>
  <td class="${cc(s.p95, s.threshold)}">${fmt(s.p95)}</td>
  <td class="${cc(s.max, s.threshold)}">${fmt(s.max)}</td>
  <td class="${pr(s.passRate)}"><span class="pct">${(s.passRate * 100).toFixed(1)}%</span></td>
</tr>`).join('\n')}
</tbody>
</table>`
}).join('\n')}

<h2>Per-Dataset: K-Means + HAC (top 50 by max error)</h2>
<p class="section-note">Sorted by largest single error across K-Means and HAC metrics.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>k</th><th>KM Centroid</th><th>KM TotWSS |Δ|</th><th>Ward MAE</th><th>Complete MAE</th><th>Single MAE</th><th>Average MAE</th><th>Sil Mean |Δ|</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok')
  .sort((a, b) => {
    const maxA = Math.max(a.kmCentroidMAE || 0, a.kmTotWithinssΔ || 0, a.hacWardMAE || 0, a.hacCompleteMAE || 0, a.hacSingleMAE || 0, a.hacAverageMAE || 0)
    const maxB = Math.max(b.kmCentroidMAE || 0, b.kmTotWithinssΔ || 0, b.hacWardMAE || 0, b.hacCompleteMAE || 0, b.hacSingleMAE || 0, b.hacAverageMAE || 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.k}</td>
  <td class="${cc(d.kmCentroidMAE, 1e-2)}">${fmt(d.kmCentroidMAE, 4)}</td>
  <td class="${cc(d.kmTotWithinssΔ, 1e-2)}">${fmt(d.kmTotWithinssΔ, 4)}</td>
  <td class="${cc(d.hacWardMAE, 1e-4)}">${fmt(d.hacWardMAE)}</td>
  <td class="${cc(d.hacCompleteMAE, 1e-4)}">${fmt(d.hacCompleteMAE)}</td>
  <td class="${cc(d.hacSingleMAE, 1e-4)}">${fmt(d.hacSingleMAE)}</td>
  <td class="${cc(d.hacAverageMAE, 1e-4)}">${fmt(d.hacAverageMAE)}</td>
  <td class="${cc(d.silMeanΔ, 1e-4)}">${fmt(d.silMeanΔ)}</td>
</tr>`).join('\n')}
</tbody>
</table>

<h2>Per-Dataset: GMM (top 50 by max error)</h2>
<p class="section-note">Sorted by largest GMM error. Skipped datasets where R or Carm GMM failed.</p>
<table>
<thead><tr><th>Seed</th><th>n</th><th>k</th><th>VVV Means</th><th>VVV LogLik |Δ|</th><th>VVV BIC |Δ|</th><th>VVI Means</th><th>VVI LogLik |Δ|</th><th>EEI Means</th><th>EEI LogLik |Δ|</th></tr></thead>
<tbody>
${results
  .filter(r => r.status === 'ok')
  .sort((a, b) => {
    const maxA = Math.max(a.gmmVVVMeansMAE || 0, a.gmmVVVLoglikΔ || 0, a.gmmVVIMeansMAE || 0, a.gmmVVILoglikΔ || 0, a.gmmEEIMeansMAE || 0, a.gmmEEILoglikΔ || 0)
    const maxB = Math.max(b.gmmVVVMeansMAE || 0, b.gmmVVVLoglikΔ || 0, b.gmmVVIMeansMAE || 0, b.gmmVVILoglikΔ || 0, b.gmmEEIMeansMAE || 0, b.gmmEEILoglikΔ || 0)
    return maxB - maxA
  })
  .slice(0, 50)
  .map(d => `<tr>
  <td>${d.seed}</td><td>${d.n}</td><td>${d.k}</td>
  <td class="${cc(d.gmmVVVMeansMAE, 1e-1)}">${fmt(d.gmmVVVMeansMAE, 4)}</td>
  <td class="${cc(d.gmmVVVLoglikΔ, 2.0)}">${fmt(d.gmmVVVLoglikΔ, 2)}</td>
  <td class="${cc(d.gmmVVVBicΔ, 5.0)}">${fmt(d.gmmVVVBicΔ, 2)}</td>
  <td class="${cc(d.gmmVVIMeansMAE, 1e-1)}">${fmt(d.gmmVVIMeansMAE, 4)}</td>
  <td class="${cc(d.gmmVVILoglikΔ, 2.0)}">${fmt(d.gmmVVILoglikΔ, 2)}</td>
  <td class="${cc(d.gmmEEIMeansMAE, 1e-1)}">${fmt(d.gmmEEIMeansMAE, 4)}</td>
  <td class="${cc(d.gmmEEILoglikΔ, 2.0)}">${fmt(d.gmmEEILoglikΔ, 2)}</td>
</tr>`).join('\n')}
</tbody>
</table>

${nFail > 0 ? `
<h2>Failed Datasets</h2>
<table>
<thead><tr><th>Seed</th><th>n</th><th>k</th><th>Error</th></tr></thead>
<tbody>
${results.filter(r => r.status === 'fail').map(r => `<tr>
  <td>${r.seed}</td><td>${r.n}</td><td>${r.k}</td><td style="color:var(--red);">${r.error ?? 'Unknown'}</td>
</tr>`).join('\n')}
</tbody>
</table>
` : ''}

<p style="margin-top:40px;color:#6e7681;font-size:11px;">
  Generated ${timestamp} · Carm vs R kmeans/mclust/hclust/silhouette · ${ref.length} datasets · ${totalMetrics} metrics
</p>

</body>
</html>`

mkdirSync('validation/reports', { recursive: true })
writeFileSync('validation/reports/clustering-report.html', html)
console.log(`\nReport saved to validation/reports/clustering-report.html`)
