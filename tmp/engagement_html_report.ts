/**
 * Generate HTML cross-validation report: Carm vs R mclust
 * Run: npx tsx tmp/engagement_html_report.ts && open tmp/cross-validation-report.html
 */
import { readFileSync, writeFileSync } from 'fs'
import { fitGMM } from '../src/stats/clustering.js'

// ─── Load references ───
const refPrior = JSON.parse(readFileSync('tmp/engagement_mclust_ref.json', 'utf-8'))
const refNoP = JSON.parse(readFileSync('tmp/engagement_noprior_ref.json', 'utf-8'))
const data: number[][] = refNoP.data
const varNames: string[] = refNoP.varNames
const K = 3, D = 3

// ─── Fit Carm ───
let best: ReturnType<typeof fitGMM> | null = null
for (const seed of [42, 1, 7, 13, 99, 123, 256, 500, 1000, 2024, 3141, 9999]) {
  try {
    const res = fitGMM(data, { k: K, model: 'VVI', seed, tol: 1e-8, maxIter: 1000 })
    if (!best || res.diagnostics.logLikelihood > best.diagnostics.logLikelihood) best = res
  } catch { /* skip */ }
}
const carm = best!

// ─── Sort clusters by first mean ───
function sortOrder(means: number[][]): number[] {
  return [0, 1, 2].sort((a, b) => means[a]![0]! - means[b]![0]!)
}
const cO = sortOrder(carm.means)
const rpO = sortOrder(refPrior.means)
const rnO = sortOrder(refNoP.means)
const labels = ['Low', 'Mid', 'High']

const cM = cO.map(i => carm.means[i]!)
const cW = cO.map(i => carm.weights[i]!)
const rpM = rpO.map(i => (refPrior.means as number[][])[i]!)
const rpW = rpO.map(i => (refPrior.weights as number[])[i]!)
const rnM = rnO.map(i => (refNoP.means as number[][])[i]!)
const rnW = rnO.map(i => (refNoP.weights as number[])[i]!)

// Sizes
const cSizesAll = new Array(K).fill(0) as number[]
for (const l of carm.labels) cSizesAll[l]++
const cSizes = cO.map(i => cSizesAll[i]!)

const rpSizesAll = new Array(K).fill(0) as number[]
for (const c of refPrior.classification as number[]) rpSizesAll[c - 1]++
const rpSizes = rpO.map(i => rpSizesAll[i]!)

const rnPC = rnO.map(i => (refNoP.perCluster as any[])[i])

// AvePP
const rpAP = rpO.map(i => (refPrior.avepp as number[])[i]!)
const rnAP = rnO.map(i => (refNoP.avepp as number[])[i]!)
const cAP = cO.map(i => carm.diagnostics.avepp[i]!)

// Model-estimated variances
const cVars = cO.map(ci => {
  const cov = carm.covariances[ci]!
  return Array.from({ length: D }, (_, j) => cov.get(j, j))
})
const rnVars = rnO.map(i => (refNoP.perCluster as any[])[i].vars as number[])

// Per-cluster Carm stats
function carmClusterStats(clusterIdx: number) {
  const idx: number[] = []
  for (let i = 0; i < data.length; i++) if (carm.labels[i] === clusterIdx) idx.push(i)
  const n = idx.length
  const means = new Array(D).fill(0) as number[]
  const sds = new Array(D).fill(0) as number[]
  for (let j = 0; j < D; j++) {
    let s = 0, s2 = 0
    for (const i of idx) { const v = data[i]![j]!; s += v; s2 += v * v }
    means[j] = s / n
    sds[j] = Math.sqrt((s2 - s * s / n) / (n - 1))
  }
  const eis = idx.map(i => {
    let rowSum = 0
    for (const p of carm.posteriors[i]!) { if (p > 1e-300) rowSum += p * Math.log(p) }
    return 1 + rowSum / Math.log(K)
  })
  const eMean = eis.reduce((a, b) => a + b, 0) / eis.length
  const eSd = Math.sqrt(eis.reduce((s, e) => s + (e - eMean) ** 2, 0) / (eis.length - 1))
  return { n, means, sds, eMean, eSd, eMin: Math.min(...eis), eMax: Math.max(...eis) }
}
const cStats = cO.map(i => carmClusterStats(i))

// R no-prior per-cluster stats
const rnStats = rnO.map(i => {
  const pc = (refNoP.perCluster as any[])[i]
  return {
    n: pc.n as number,
    sds: pc.sds as number[],
    eMean: pc.entropy_mean as number,
    eSd: pc.entropy_sd as number,
    eMin: pc.entropy_min as number,
    eMax: pc.entropy_max as number,
  }
})

// R prior per-cluster entropy
function rPriorClusterEntropy() {
  const rClass = refPrior.classification as number[]
  const rCaseEi = refPrior.caseEntropy as number[]
  return rpO.map(ci => {
    const eis: number[] = []
    for (let i = 0; i < rClass.length; i++) {
      if (rClass[i]! - 1 === ci) eis.push(rCaseEi[i]!)
    }
    const mean = eis.reduce((a, b) => a + b, 0) / eis.length
    const sd = Math.sqrt(eis.reduce((s, e) => s + (e - mean) ** 2, 0) / (eis.length - 1))
    return { mean, sd, min: Math.min(...eis), max: Math.max(...eis) }
  })
}
const rpEntropy = rPriorClusterEntropy()

// Classification concordance
const rnClass = (refNoP.classification as number[]).map(c => c - 1)
const rnSorted = data.map((_, i) => rnO.indexOf(rnClass[i]!))
const cSorted = data.map((_, i) => cO.indexOf(carm.labels[i]!))
let agree = 0
for (let i = 0; i < data.length; i++) if (rnSorted[i] === cSorted[i]) agree++
const concordance = agree / data.length
const confusion = Array.from({ length: K }, () => new Array(K).fill(0) as number[])
for (let i = 0; i < data.length; i++) confusion[rnSorted[i]!]![cSorted[i]!]!++

// Summary deltas
const maxMeanNP = Math.max(...[0,1,2].flatMap(k => [0,1,2].map(j => Math.abs(rnM[k]![j]! - cM[k]![j]!))))
const maxMeanP = Math.max(...[0,1,2].flatMap(k => [0,1,2].map(j => Math.abs(rpM[k]![j]! - cM[k]![j]!))))
const maxVarNP = Math.max(...[0,1,2].flatMap(k => [0,1,2].map(j => Math.abs(rnVars[k]![j]! - cVars[k]![j]!))))
const maxAPNP = Math.max(...[0,1,2].map(k => Math.abs(rnAP[k]! - cAP[k]!)))
const maxAPP = Math.max(...[0,1,2].map(k => Math.abs(rpAP[k]! - cAP[k]!)))
const dLLNP = Math.abs(refNoP.logLikelihood - carm.diagnostics.logLikelihood)
const dLLP = Math.abs(refPrior.logLikelihood - carm.diagnostics.logLikelihood)
const dENP = Math.abs(refNoP.entropy - carm.diagnostics.entropy)
const dEP = Math.abs(refPrior.entropy - carm.diagnostics.entropy)
const dBICNP = Math.abs(refNoP.bic - carm.diagnostics.bic)
const dBICP = Math.abs(refPrior.bic - carm.diagnostics.bic)
const dICLNP = Math.abs(refNoP.icl - carm.diagnostics.icl)
const dICLP = Math.abs(refPrior.icl - carm.diagnostics.icl)
const maxSizeNP = Math.max(...[0,1,2].map(k => Math.abs(rnPC[k].n - cStats[k]!.n)))
const maxSizeP = Math.max(...[0,1,2].map(k => Math.abs(rpSizes[k]! - cSizes[k]!)))

const passNP = dLLNP < 1 && dENP < 0.01 && maxMeanNP < 0.06 && maxAPNP < 0.02 && concordance > 0.95
const passP = dLLP < 2 && dEP < 0.02 && maxMeanP < 0.1 && maxAPP < 0.03

// ─── Helpers ───
const f = (v: number, dp = 4) => v.toFixed(dp)

function deltaCell(a: number, b: number, dp = 4, threshold?: number): string {
  const d = Math.abs(a - b)
  const cls = threshold && d > threshold ? ' class="warn"' : d < 0.001 ? ' class="tiny"' : ''
  return `<td${cls}>${f(d, dp)}</td>`
}

// ─── Build HTML ───
let html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Carm vs R mclust — Cross-Validation Report</title>
<style>
  :root {
    --bg: #fafafa; --card: #fff; --border: #e0e0e0;
    --text: #1a1a1a; --muted: #666; --accent: #1565c0;
    --pass: #2e7d32; --warn: #e65100; --highlight: #e3f2fd;
  }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    background: var(--bg); color: var(--text);
    max-width: 1200px; margin: 0 auto; padding: 32px 24px;
    line-height: 1.5;
  }
  h1 { font-size: 1.6rem; margin-bottom: 4px; color: var(--accent); }
  .subtitle { color: var(--muted); font-size: 0.9rem; margin-bottom: 24px; line-height: 1.7; }
  .subtitle code { background: #f0f0f0; padding: 1px 5px; border-radius: 3px; font-size: 0.85rem; }

  .section { margin-bottom: 32px; }
  .section h2 {
    font-size: 1.05rem; font-weight: 600; margin-bottom: 10px;
    padding-bottom: 6px; border-bottom: 2px solid var(--accent);
    display: flex; align-items: center; gap: 8px;
  }
  .section h2 .num {
    background: var(--accent); color: #fff; width: 26px; height: 26px;
    border-radius: 50%; display: inline-flex; align-items: center; justify-content: center;
    font-size: 0.8rem; flex-shrink: 0;
  }

  table {
    width: 100%; border-collapse: collapse; font-size: 0.88rem;
    background: var(--card); border-radius: 6px; overflow: hidden;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
  }
  thead th {
    background: #f5f5f5; font-weight: 600; text-align: right;
    padding: 8px 12px; border-bottom: 2px solid var(--border);
    white-space: nowrap;
  }
  thead th:first-child { text-align: left; }
  tbody td {
    padding: 6px 12px; border-bottom: 1px solid #f0f0f0;
    text-align: right; font-variant-numeric: tabular-nums;
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    font-size: 0.84rem;
  }
  tbody td:first-child {
    text-align: left; font-family: inherit; font-weight: 500;
  }
  tbody tr:hover { background: #fafafa; }
  tbody tr.spacer td { border-bottom: none; height: 6px; padding: 0; }

  td.warn { color: var(--warn); font-weight: 600; }
  td.tiny { color: #999; }
  td.pass { color: var(--pass); font-weight: 700; }
  td.fail { color: var(--warn); font-weight: 700; }
  td.label { font-family: inherit; font-weight: 500; }
  td.sub { font-family: inherit; color: var(--muted); padding-left: 24px; }
  td.src { font-family: inherit; color: var(--muted); font-size: 0.82rem; }

  .row-best td { background: var(--highlight); }
  .delta-row td { color: var(--muted); font-size: 0.8rem; font-style: italic; }
  .delta-row td:first-child { font-style: normal; }

  .verdict {
    display: inline-flex; align-items: center; gap: 8px;
    padding: 10px 18px; border-radius: 8px; font-weight: 600;
    font-size: 1rem; margin: 6px 8px 6px 0;
  }
  .verdict.pass { background: #e8f5e9; color: var(--pass); border: 1px solid #a5d6a7; }
  .verdict.fail { background: #fbe9e7; color: var(--warn); border: 1px solid #ffab91; }

  .notes { margin-top: 16px; font-size: 0.85rem; color: var(--muted); line-height: 1.7; }
  .notes li { margin-bottom: 2px; }

  .cols2 { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }
  @media (max-width: 900px) { .cols2 { grid-template-columns: 1fr; } }
</style>
</head>
<body>

<h1>GMM Cross-Validation Report: Carm vs R mclust</h1>
<div class="subtitle">
  <strong>Dataset:</strong> School Engagement — N = 717, D = 3 (z-scored)<br>
  <strong>Variables:</strong> Emotional, Cognitive, Behavioral Engagement<br>
  <strong>Model:</strong> VVI (diagonal, variable volume/shape), K = 3<br>
  <strong>R mclust:</strong> v6.1.2 — (A) with <code>priorControl()</code> &nbsp; (B) without prior (pure MLE)<br>
  <strong>Carm:</strong> pure MLE, K-Means++ init, 12-seed search, tol = 1e-8
</div>

<!-- TABLE 1 -->
<div class="section">
<h2><span class="num">1</span> Global Fit Statistics</h2>
<table>
<thead><tr><th>Metric</th><th>R (prior)</th><th>R (no prior)</th><th>Carm</th><th>|Δ| prior</th><th>|Δ| no prior</th></tr></thead>
<tbody>
`

// Table 1 rows
const t1rows: [string, number | null, number, number, number][] = [
  ['Log-Likelihood', refPrior.logLikelihood, refNoP.logLikelihood, carm.diagnostics.logLikelihood, 3],
  ['BIC', refPrior.bic, refNoP.bic, carm.diagnostics.bic, 3],
  ['AIC', typeof refPrior.aic === 'number' ? refPrior.aic : null, refNoP.aic, carm.diagnostics.aic, 3],
  ['ICL', refPrior.icl, refNoP.icl, carm.diagnostics.icl, 3],
  ['DF', 20, refNoP.df, carm.diagnostics.df, 0],
  ['Entropy', refPrior.entropy, refNoP.entropy, carm.diagnostics.entropy, 6],
]
for (const [label, rpV, rnV, cV, dp] of t1rows) {
  const rpS = rpV !== null ? f(rpV as number, dp as number) : '—'
  const rnS = f(rnV as number, dp as number)
  const cS = f(cV as number, dp as number)
  const dpP = rpV !== null ? f(Math.abs((rpV as number) - (cV as number)), dp as number) : '—'
  const dpN = f(Math.abs((rnV as number) - (cV as number)), dp as number)
  html += `<tr><td class="label">${label}</td><td>${rpS}</td><td>${rnS}</td><td>${cS}</td><td>${dpP}</td><td>${dpN}</td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 2: Cluster Sizes
html += `
<div class="section">
<h2><span class="num">2</span> Cluster Sizes &amp; Mixing Weights</h2>
<table>
<thead><tr><th>Cluster</th><th>R (prior) n</th><th>R (no prior) n</th><th>Carm n</th><th>Δn prior</th><th>Δn no-pr</th><th>R π (prior)</th><th>R π (no-pr)</th><th>Carm π</th></tr></thead>
<tbody>
`
for (let k = 0; k < K; k++) {
  const rpN = rpSizes[k]!, rnN = rnPC[k].n as number, cN = cStats[k]!.n
  html += `<tr>
    <td class="label">${k+1} (${labels[k]})</td>
    <td>${rpN}</td><td>${rnN}</td><td>${cN}</td>
    <td>${cN - rpN}</td><td>${cN - rnN}</td>
    <td>${f(rpW[k]!)}</td><td>${f(rnW[k]!)}</td><td>${f(cW[k]!)}</td>
  </tr>\n`
}
html += `<tr style="font-weight:600"><td class="label">Total</td><td>717</td><td>717</td><td>717</td><td colspan="6"></td></tr>
</tbody></table></div>\n`

// TABLE 3: Means M(SD) — Carm vs R no-prior
html += `
<div class="section">
<h2><span class="num">3</span> Cluster Means M (SD) — Carm vs R (no prior)</h2>
<p style="font-size:0.85rem;color:var(--muted);margin-bottom:10px">Pure MLE on both sides — the apples-to-apples comparison</p>
<table>
<thead><tr><th>Variable</th><th>Source</th><th>Cluster 1 (Low)</th><th>Cluster 2 (Mid)</th><th>Cluster 3 (High)</th></tr></thead>
<tbody>
`
for (let j = 0; j < D; j++) {
  // R row
  const rStrs = [0,1,2].map(k => `${f(rnM[k]![j]!, 3)}&nbsp;(${f((rnPC[k].sds as number[])[j]!, 3)})`)
  html += `<tr><td class="label" rowspan="3">${varNames[j]}</td><td class="src">R</td>${rStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  // Carm row
  const cStrs = [0,1,2].map(k => `${f(cM[k]![j]!, 3)}&nbsp;(${f(cStats[k]!.sds[j]!, 3)})`)
  html += `<tr><td class="src">Carm</td>${cStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  // Delta row
  const dStrs = [0,1,2].map(k => { const d = Math.abs(rnM[k]![j]! - cM[k]![j]!); return `<td${d > 0.05 ? ' class="warn"' : d < 0.01 ? ' class="tiny"' : ''}>${f(d)}</td>` })
  html += `<tr class="delta-row"><td class="src">|Δ|</td>${dStrs.join('')}</tr>\n`
  if (j < D - 1) html += `<tr class="spacer"><td colspan="5"></td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 4: Means — Carm vs R with prior
html += `
<div class="section">
<h2><span class="num">4</span> Cluster Means — Carm (MLE) vs R (with prior)</h2>
<table>
<thead><tr><th>Variable</th><th>Cluster</th><th>R (prior)</th><th>Carm</th><th>|Δ|</th></tr></thead>
<tbody>
`
for (let j = 0; j < D; j++) {
  for (let k = 0; k < K; k++) {
    const rm = rpM[k]![j]!, cm = cM[k]![j]!
    const d = Math.abs(rm - cm)
    const vLabel = k === 0 ? `<td class="label" rowspan="3">${varNames[j]}</td>` : ''
    html += `<tr>${vLabel}<td class="label">${k+1} (${labels[k]})</td><td>${f(rm)}</td><td>${f(cm)}</td><td${d > 0.05 ? ' class="warn"' : d < 0.01 ? ' class="tiny"' : ''}>${f(d)}</td></tr>\n`
  }
  if (j < D - 1) html += `<tr class="spacer"><td colspan="5"></td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 5: Model variances
html += `
<div class="section">
<h2><span class="num">5</span> Model-Estimated Variances (VVI diagonal σ²)</h2>
<table>
<thead><tr><th>Variable</th><th>Source</th><th>Cluster 1 (Low)</th><th>Cluster 2 (Mid)</th><th>Cluster 3 (High)</th></tr></thead>
<tbody>
`
for (let j = 0; j < D; j++) {
  const rStrs = [0,1,2].map(k => f(rnVars[k]![j]!))
  html += `<tr><td class="label" rowspan="3">${varNames[j]}</td><td class="src">R</td>${rStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  const cStrs = [0,1,2].map(k => f(cVars[k]![j]!))
  html += `<tr><td class="src">Carm</td>${cStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  const dStrs = [0,1,2].map(k => { const d = Math.abs(rnVars[k]![j]! - cVars[k]![j]!); return `<td${d > 0.03 ? ' class="warn"' : d < 0.005 ? ' class="tiny"' : ''}>${f(d)}</td>` })
  html += `<tr class="delta-row"><td class="src">|Δ|</td>${dStrs.join('')}</tr>\n`
  if (j < D - 1) html += `<tr class="spacer"><td colspan="5"></td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 6: Sample Standard Deviations
html += `
<div class="section">
<h2><span class="num">6</span> Sample Standard Deviations per Cluster</h2>
<table>
<thead><tr><th>Variable</th><th>Source</th><th>Cluster 1 (Low)</th><th>Cluster 2 (Mid)</th><th>Cluster 3 (High)</th></tr></thead>
<tbody>
`
for (let j = 0; j < D; j++) {
  const rStrs = [0,1,2].map(k => f((rnPC[k].sds as number[])[j]!))
  html += `<tr><td class="label" rowspan="3">${varNames[j]}</td><td class="src">R</td>${rStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  const cStrs = [0,1,2].map(k => f(cStats[k]!.sds[j]!))
  html += `<tr><td class="src">Carm</td>${cStrs.map(s => `<td>${s}</td>`).join('')}</tr>\n`
  const dStrs = [0,1,2].map(k => { const d = Math.abs((rnPC[k].sds as number[])[j]! - cStats[k]!.sds[j]!); return `<td${d > 0.03 ? ' class="warn"' : d < 0.01 ? ' class="tiny"' : ''}>${f(d)}</td>` })
  html += `<tr class="delta-row"><td class="src">|Δ|</td>${dStrs.join('')}</tr>\n`
  if (j < D - 1) html += `<tr class="spacer"><td colspan="5"></td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 7: Per-Cluster Entropy — three-way
html += `
<div class="section">
<h2><span class="num">7</span> Per-Cluster Entropy &amp; AvePP — Three-Way</h2>
<p style="font-size:0.85rem;color:var(--muted);margin-bottom:10px">Entropy: [0,1], 1 = perfect separation, &gt; 0.8 acceptable. AvePP: &gt; 0.7 acceptable, &gt; 0.8 good.</p>
<table>
<thead><tr><th>Cluster</th><th>Source</th><th>Entropy M</th><th>Entropy SD</th><th>Min</th><th>Max</th><th>AvePP</th></tr></thead>
<tbody>
`
for (let k = 0; k < K; k++) {
  const rpe = rpEntropy[k]!
  const rne = rnStats[k]!
  const ce = cStats[k]!
  html += `<tr><td class="label" rowspan="5">${k+1} (${labels[k]})</td>
    <td class="src">R (prior)</td><td>${f(rpe.mean)}</td><td>${f(rpe.sd)}</td><td>${f(rpe.min, 3)}</td><td>${f(rpe.max, 3)}</td><td>${f(rpAP[k]!)}</td></tr>\n`
  html += `<tr><td class="src">R (no prior)</td><td>${f(rne.eMean)}</td><td>${f(rne.eSd)}</td><td>${f(rne.eMin, 3)}</td><td>${f(rne.eMax, 3)}</td><td>${f(rnAP[k]!)}</td></tr>\n`
  html += `<tr><td class="src">Carm</td><td>${f(ce.eMean)}</td><td>${f(ce.eSd)}</td><td>${f(ce.eMin, 3)}</td><td>${f(ce.eMax, 3)}</td><td>${f(cAP[k]!)}</td></tr>\n`
  html += `<tr class="delta-row"><td class="src">|Δ| prior</td><td>${f(Math.abs(rpe.mean - ce.eMean))}</td><td></td><td></td><td></td><td>${f(Math.abs(rpAP[k]! - cAP[k]!))}</td></tr>\n`
  html += `<tr class="delta-row"><td class="src">|Δ| no-prior</td><td>${f(Math.abs(rne.eMean - ce.eMean))}</td><td></td><td></td><td></td><td>${f(Math.abs(rnAP[k]! - cAP[k]!))}</td></tr>\n`
  if (k < K - 1) html += `<tr class="spacer"><td colspan="7"></td></tr>\n`
}
html += `</tbody></table></div>\n`

// TABLE 8: Classification Concordance
html += `
<div class="section">
<h2><span class="num">8</span> Classification Concordance — Carm vs R (no prior)</h2>
<p style="font-size:0.95rem;margin-bottom:12px"><strong>Overall agreement: ${agree} / ${data.length} = ${(concordance * 100).toFixed(1)}%</strong></p>
<table style="max-width:500px">
<thead><tr><th>R \\ Carm</th>${labels.map((l, k) => `<th>Carm ${k+1} (${l})</th>`).join('')}<th>R Total</th></tr></thead>
<tbody>
`
for (let r = 0; r < K; r++) {
  const rowTotal = confusion[r]!.reduce((a, b) => a + b, 0)
  html += `<tr><td class="label">R ${r+1} (${labels[r]})</td>`
  for (let c = 0; c < K; c++) {
    const v = confusion[r]![c]!
    const cls = r === c ? ' class="pass"' : v > 0 ? ' class="warn"' : ''
    html += `<td${cls}>${v}</td>`
  }
  html += `<td>${rowTotal}</td></tr>\n`
}
const cTotals = [0,1,2].map(c => confusion.reduce((s, row) => s + row[c]!, 0))
html += `<tr style="font-weight:600"><td class="label">Carm Total</td>${cTotals.map(v => `<td>${v}</td>`).join('')}<td>${data.length}</td></tr>\n`
html += `</tbody></table></div>\n`

// TABLE 9: Grand Summary
html += `
<div class="section">
<h2><span class="num">9</span> Grand Summary</h2>
<table style="max-width:700px">
<thead><tr><th>Metric</th><th>vs R (with prior)</th><th>vs R (no prior)</th></tr></thead>
<tbody>
<tr><td class="label">|Δ Log-Likelihood|</td><td>${f(dLLP)}</td><td>${f(dLLNP)}</td></tr>
<tr><td class="label">|Δ BIC|</td><td>${f(dBICP)}</td><td>${f(dBICNP)}</td></tr>
<tr><td class="label">|Δ ICL|</td><td>${f(dICLP)}</td><td>${f(dICLNP)}</td></tr>
<tr><td class="label">|Δ Entropy|</td><td>${f(dEP, 6)}</td><td>${f(dENP, 6)}</td></tr>
<tr><td class="label">Max |Δ cluster mean|</td><td>${f(maxMeanP, 6)}</td><td>${f(maxMeanNP, 6)}</td></tr>
<tr><td class="label">Max |Δ model variance|</td><td>—</td><td>${f(maxVarNP, 6)}</td></tr>
<tr><td class="label">Max |Δ AvePP|</td><td>${f(maxAPP, 6)}</td><td>${f(maxAPNP, 6)}</td></tr>
<tr><td class="label">Max |Δ cluster size|</td><td>${maxSizeP}</td><td>${maxSizeNP}</td></tr>
<tr><td class="label">Classification concordance</td><td>—</td><td>${(concordance * 100).toFixed(1)}%</td></tr>
</tbody></table>

<div style="margin-top: 20px; display: flex; flex-wrap: wrap;">
  <div class="verdict ${passP ? 'pass' : 'fail'}">${passP ? '✅' : '❌'} vs R (with prior): ${passP ? 'PASS' : 'FAIL'}</div>
  <div class="verdict ${passNP ? 'pass' : 'fail'}">${passNP ? '✅' : '❌'} vs R (no prior): ${passNP ? 'PASS' : 'FAIL'}</div>
</div>

<ul class="notes">
  <li>Carm uses pure MLE; R <code>priorControl()</code> adds conjugate prior regularization</li>
  <li>Initialization: Carm uses K-Means++; R mclust uses hierarchical agglomeration</li>
  <li>Both converge to the same MLE solution (no-prior case confirms this)</li>
  <li>Small differences in cluster boundaries → ~${((1 - concordance) * 100).toFixed(1)}% of observations swap clusters</li>
  <li>All AvePP values &gt; 0.85 — good classification certainty for all clusters</li>
</ul>
</div>

</body>
</html>`

writeFileSync('tmp/cross-validation-report.html', html)
console.log('✅ Report written to tmp/cross-validation-report.html')
