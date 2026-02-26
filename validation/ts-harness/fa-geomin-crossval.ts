#!/usr/bin/env npx tsx
/**
 * Cross-validate Carm geomin rotation against R GPArotation::geominQ
 */
import { readFileSync } from 'fs'
import { runEFA } from 'carm'

// Load synthetic datasets
const datasets = JSON.parse(
  readFileSync('tmp/fa-crossval-data.json', 'utf-8')
).datasets as Array<{
  id: number
  params: { n: number; k: number }
  data: number[][]
  variableNames: string[]
}>

// Load R geomin reference
const rRef = JSON.parse(
  readFileSync('tmp/fa-geomin-ref.json', 'utf-8')
) as Array<{
  id: number; n: number; p: number; k: number
  loadings: number[][]; phi: number[][]; converged: boolean
}>

// Build lookup
const refMap = new Map(rRef.map(r => [r.id, r]))

// Permutation + sign matching
function permutations(n: number): number[][] {
  if (n === 1) return [[0]]
  const result: number[][] = []
  const sub = permutations(n - 1)
  for (let i = 0; i < n; i++) {
    for (const perm of sub) {
      result.push([i, ...perm.map(x => x >= i ? x + 1 : x)])
    }
  }
  return result
}

function bestMAE(cLoad: number[][], rLoad: number[][], p: number, k: number): {
  mae: number; perm: number[]; signs: number[]
} {
  let bestMae = Infinity
  let bestPerm: number[] = []
  let bestSigns: number[] = []

  for (const perm of permutations(k)) {
    for (let signMask = 0; signMask < (1 << k); signMask++) {
      const signs = Array.from({ length: k }, (_, j) =>
        (signMask & (1 << j)) ? -1 : 1
      )
      let sumAE = 0
      for (let i = 0; i < p; i++) {
        for (let j = 0; j < k; j++) {
          sumAE += Math.abs(
            cLoad[i]![perm[j]!]! * signs[j]! - rLoad[i]![j]!
          )
        }
      }
      const mae = sumAE / (p * k)
      if (mae < bestMae) {
        bestMae = mae
        bestPerm = perm
        bestSigns = signs
      }
    }
  }
  return { mae: bestMae, perm: bestPerm, signs: bestSigns }
}

// Run cross-validation
let passed = 0, failed = 0
const threshold = 0.05

for (const ds of datasets) {
  const ref = refMap.get(ds.id)
  if (!ref) continue

  const result = runEFA(ds.data, {
    nFactors: ds.params.k,
    extraction: 'ml',
    rotation: 'geomin',
    geominDelta: 0.01,
    variableNames: ds.variableNames,
  })

  const cLoad = result.loadings.map(r => [...r])
  const { mae, perm, signs } = bestMAE(cLoad, ref.loadings, ref.p, ref.k)

  const pass = mae <= threshold
  if (pass) passed++; else failed++

  const marker = pass ? '✓' : '✗'
  console.log(
    `[${String(ds.id).padStart(3)}] n=${ref.n} p=${ref.p} k=${ref.k} ${marker} loadMAE=${mae.toFixed(4)}`
  )
}

console.log(`\nDone: ${passed} passed, ${failed} failed`)

// Also test real dataset
console.log('\n=== Real dataset geomin ===')
const csv = readFileSync('/Users/mohammedsaqr/Downloads/rraw_dataaw_data.csv', 'utf-8')
const lines = csv.trim().split('\n')
const header = lines[0]!.split(',')
const realData = lines.slice(1).map(l => l.split(',').map(Number))

const realRef = JSON.parse(
  readFileSync('tmp/fa-geomin-real-ref.json', 'utf-8')
) as Record<string, { k: number; loadings: number[][]; phi: number[][] }>

for (const kStr of ['3', '5']) {
  const rr = realRef[kStr]!
  const result = runEFA(realData, {
    nFactors: rr.k,
    extraction: 'ml',
    rotation: 'geomin',
    geominDelta: 0.01,
    variableNames: header,
  })

  const cLoad = result.loadings.map(r => [...r])
  const { mae } = bestMAE(cLoad, rr.loadings, header.length, rr.k)
  const pass = mae <= threshold
  console.log(`  k=${rr.k}: ${pass ? '✓' : '✗'} loadMAE=${mae.toFixed(4)}`)

  // Show first 3 rows comparison
  for (let i = 0; i < 3; i++) {
    console.log(`    ${header[i]}: Carm=[${cLoad[i]!.map(v => v.toFixed(4)).join(', ')}]`)
    console.log(`    ${' '.repeat(header[i]!.length)}: R   =[${rr.loadings[i]!.map(v => v.toFixed(4)).join(', ')}]`)
  }
}
