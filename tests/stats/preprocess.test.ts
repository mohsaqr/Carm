/**
 * Tests for preprocessing module.
 *
 * Cross-validated with R:
 * > scale(data)                       # standardize
 * > scale(data, scale = FALSE)        # center only
 * > log(data)                         # log
 * > sqrt(data)                        # sqrt
 */

import { describe, it, expect } from 'vitest'
import { preprocessData, inverseTransform } from '../../src/stats/preprocess.js'

// ─── Test data ───────────────────────────────────────────────────────────

const DATA = [
  [1, 10, 100],
  [2, 20, 200],
  [3, 30, 300],
  [4, 40, 400],
  [5, 50, 500],
]

// Cross-validated with R:
// > d <- matrix(c(1,2,3,4,5, 10,20,30,40,50, 100,200,300,400,500), ncol=3)
// > scale(d, center=TRUE, scale=FALSE)
// > scale(d)
// > colMeans(d)     # 3, 30, 300
// > apply(d, 2, sd) # 1.581139, 15.81139, 158.1139

describe('preprocessData', () => {
  it('method=none passes data through unchanged', () => {
    const result = preprocessData(DATA, { method: 'none' })
    expect(result.data).toEqual(DATA)
    expect(result.method).toBe('none')
    expect(result.centered).toBe(false)
    expect(result.scaled).toBe(false)
    expect(result.colMeans[0]).toBeCloseTo(3, 10)
    expect(result.colMeans[1]).toBeCloseTo(30, 10)
    expect(result.colMeans[2]).toBeCloseTo(300, 10)
  })

  it('method=center matches R scale(x, scale=FALSE)', () => {
    // R: scale(d, scale=FALSE)
    // > [1] -2, -1, 0, 1, 2 for col1
    const result = preprocessData(DATA, { method: 'center' })
    expect(result.centered).toBe(true)
    expect(result.scaled).toBe(false)
    expect(result.data[0]![0]).toBeCloseTo(-2, 10) // 1 - 3
    expect(result.data[2]![0]).toBeCloseTo(0, 10)  // 3 - 3
    expect(result.data[4]![0]).toBeCloseTo(2, 10)  // 5 - 3
    expect(result.data[0]![1]).toBeCloseTo(-20, 10) // 10 - 30
    expect(result.data[0]![2]).toBeCloseTo(-200, 10) // 100 - 300
  })

  it('method=standardize matches R scale(x)', () => {
    // R: scale(d)
    // col1: sd=1.581139, so (1-3)/1.581139 = -1.264911, etc.
    const result = preprocessData(DATA, { method: 'standardize' })
    expect(result.centered).toBe(true)
    expect(result.scaled).toBe(true)

    // R: (1-3)/sd(c(1,2,3,4,5)) = -2/1.581139 = -1.264911
    expect(result.data[0]![0]).toBeCloseTo(-1.264911, 5)
    expect(result.data[2]![0]).toBeCloseTo(0, 10)
    expect(result.data[4]![0]).toBeCloseTo(1.264911, 5)

    // All columns should have same standardized values (since they're proportional)
    expect(result.data[0]![1]).toBeCloseTo(-1.264911, 5)
    expect(result.data[0]![2]).toBeCloseTo(-1.264911, 5)
  })

  it('method=log transforms correctly', () => {
    const result = preprocessData(DATA, { method: 'log' })
    expect(result.data[0]![0]).toBeCloseTo(Math.log(1), 10)
    expect(result.data[1]![0]).toBeCloseTo(Math.log(2), 10)
    expect(result.data[0]![1]).toBeCloseTo(Math.log(10), 10)
    expect(result.data[0]![2]).toBeCloseTo(Math.log(100), 10)
    expect(result.method).toBe('log')
  })

  it('method=sqrt transforms correctly', () => {
    const result = preprocessData(DATA, { method: 'sqrt' })
    expect(result.data[0]![0]).toBeCloseTo(1, 10)
    expect(result.data[1]![0]).toBeCloseTo(Math.sqrt(2), 10)
    expect(result.data[0]![1]).toBeCloseTo(Math.sqrt(10), 10)
    expect(result.data[0]![2]).toBeCloseTo(10, 10)
    expect(result.method).toBe('sqrt')
  })

  it('log throws on non-positive values', () => {
    expect(() => preprocessData([[0, 1], [1, 2]], { method: 'log' })).toThrow()
    expect(() => preprocessData([[-1, 1], [1, 2]], { method: 'log' })).toThrow()
  })

  it('sqrt throws on negative values', () => {
    expect(() => preprocessData([[-1, 1], [1, 2]], { method: 'sqrt' })).toThrow()
  })

  it('sqrt accepts zero values', () => {
    const result = preprocessData([[0, 1], [4, 9]], { method: 'sqrt' })
    expect(result.data[0]![0]).toBeCloseTo(0, 10)
    expect(result.data[1]![0]).toBeCloseTo(2, 10)
    expect(result.data[1]![1]).toBeCloseTo(3, 10)
  })

  it('handles zero-variance columns (SD=1)', () => {
    const constData = [[5, 1], [5, 2], [5, 3]]
    const result = preprocessData(constData, { method: 'standardize' })
    // Col 0 is constant (5,5,5) → centered to (0,0,0), SD=1 so no scaling
    expect(result.data[0]![0]).toBeCloseTo(0, 10)
    expect(result.data[1]![0]).toBeCloseTo(0, 10)
    expect(result.colSDs[0]).toBe(1) // zero-variance → replaced with 1
  })

  it('handles single row (SD=0 for all)', () => {
    const singleRow = [[1, 2, 3]]
    const result = preprocessData(singleRow, { method: 'standardize' })
    // SD = 0 for all columns (n=1), replaced with 1
    expect(result.data[0]![0]).toBeCloseTo(0, 10)  // centered: 1 - 1 = 0
    expect(result.data[0]![1]).toBeCloseTo(0, 10)
  })

  it('throws on empty data', () => {
    expect(() => preprocessData([])).toThrow()
  })

  it('default method is none', () => {
    const result = preprocessData(DATA)
    expect(result.method).toBe('none')
    expect(result.data).toEqual(DATA)
  })

  it('colMeans and colSDs are computed correctly for R scale()', () => {
    // R: colMeans(d) = c(3, 30, 300)
    // R: apply(d, 2, sd) = c(1.581139, 15.81139, 158.1139)
    const result = preprocessData(DATA, { method: 'standardize' })
    expect(result.colMeans[0]).toBeCloseTo(3, 10)
    expect(result.colMeans[1]).toBeCloseTo(30, 10)
    expect(result.colMeans[2]).toBeCloseTo(300, 10)
    expect(result.colSDs[0]).toBeCloseTo(1.581139, 4)
    expect(result.colSDs[1]).toBeCloseTo(15.81139, 3)
    expect(result.colSDs[2]).toBeCloseTo(158.1139, 2)
  })
})

describe('inverseTransform', () => {
  it('round-trip center: inverse(preprocess(x)) ≈ x', () => {
    const pp = preprocessData(DATA, { method: 'center' })
    const restored = inverseTransform(pp.data, pp)
    for (let i = 0; i < DATA.length; i++) {
      for (let j = 0; j < DATA[0]!.length; j++) {
        expect(restored[i]![j]).toBeCloseTo(DATA[i]![j]!, 10)
      }
    }
  })

  it('round-trip standardize: inverse(preprocess(x)) ≈ x', () => {
    const pp = preprocessData(DATA, { method: 'standardize' })
    const restored = inverseTransform(pp.data, pp)
    for (let i = 0; i < DATA.length; i++) {
      for (let j = 0; j < DATA[0]!.length; j++) {
        expect(restored[i]![j]).toBeCloseTo(DATA[i]![j]!, 8)
      }
    }
  })

  it('round-trip log: inverse(preprocess(x)) ≈ x', () => {
    const pp = preprocessData(DATA, { method: 'log' })
    const restored = inverseTransform(pp.data, pp)
    for (let i = 0; i < DATA.length; i++) {
      for (let j = 0; j < DATA[0]!.length; j++) {
        expect(restored[i]![j]).toBeCloseTo(DATA[i]![j]!, 8)
      }
    }
  })

  it('round-trip sqrt: inverse(preprocess(x)) ≈ x', () => {
    const pp = preprocessData(DATA, { method: 'sqrt' })
    const restored = inverseTransform(pp.data, pp)
    for (let i = 0; i < DATA.length; i++) {
      for (let j = 0; j < DATA[0]!.length; j++) {
        expect(restored[i]![j]).toBeCloseTo(DATA[i]![j]!, 8)
      }
    }
  })

  it('round-trip none: pass-through', () => {
    const pp = preprocessData(DATA, { method: 'none' })
    const restored = inverseTransform(pp.data, pp)
    expect(restored).toEqual(DATA)
  })
})
