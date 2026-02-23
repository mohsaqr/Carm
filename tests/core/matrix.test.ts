/**
 * Tests for src/core/matrix.ts
 * Cross-validated against R:
 * > A <- matrix(c(1,3,2,4), nrow=2); solve(A)
 * > A %*% B
 */
import { describe, it, expect } from 'vitest'
import { Matrix, solveLinear } from '../../src/core/matrix.js'
import gt from '../fixtures/ground_truth.json'

describe('Matrix', () => {
  describe('construction', () => {
    it('creates from 2D array', () => {
      const m = Matrix.fromArray([[1, 2], [3, 4]])
      expect(m.rows).toBe(2)
      expect(m.cols).toBe(2)
      expect(m.get(0, 0)).toBe(1)
      expect(m.get(1, 1)).toBe(4)
    })

    it('throws on mismatched row lengths', () => {
      expect(() => Matrix.fromArray([[1, 2], [3]])).toThrow()
    })

    it('creates identity matrix', () => {
      const I = Matrix.identity(3)
      expect(I.get(0, 0)).toBe(1)
      expect(I.get(0, 1)).toBe(0)
      expect(I.get(2, 2)).toBe(1)
    })
  })

  describe('transpose', () => {
    it('transposes 2×3 to 3×2', () => {
      const m = Matrix.fromArray([[1, 2, 3], [4, 5, 6]])
      const t = m.transpose()
      expect(t.rows).toBe(3)
      expect(t.cols).toBe(2)
      expect(t.get(0, 0)).toBe(1)
      expect(t.get(2, 1)).toBe(6)
    })
  })

  describe('multiply', () => {
    it('multiplies 2×2 matrices (R ground truth)', () => {
      const A = Matrix.fromArray(gt.matrix.multiply_2x2.A as number[][])
      const B = Matrix.fromArray(gt.matrix.multiply_2x2.B as number[][])
      const C = A.multiply(B)
      const expected = gt.matrix.multiply_2x2.result as number[][]
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 2; j++) {
          expect(C.get(i, j)).toBeCloseTo(expected[i]![j]!, 8)
        }
      }
    })

    it('throws on dimension mismatch', () => {
      const A = Matrix.fromArray([[1, 2], [3, 4]])
      const B = Matrix.fromArray([[1, 2, 3]])
      expect(() => A.multiply(B)).toThrow()
    })
  })

  describe('inverse', () => {
    it('inverts 2×2 matrix (R ground truth)', () => {
      const A = Matrix.fromArray(gt.matrix.inverse_2x2.A as number[][])
      const Ainv = A.inverse()
      const expected = gt.matrix.inverse_2x2.result as number[][]
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 2; j++) {
          expect(Ainv.get(i, j)).toBeCloseTo(expected[i]![j]!, 6)
        }
      }
    })

    it('inverts 3×3 matrix (R ground truth)', () => {
      const A = Matrix.fromArray(gt.matrix.inverse_3x3.A as number[][])
      const Ainv = A.inverse()
      const expected = gt.matrix.inverse_3x3.result as number[][]
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(Ainv.get(i, j)).toBeCloseTo(expected[i]![j]!, 4)
        }
      }
    })

    it('A * inv(A) ≈ I', () => {
      const A = Matrix.fromArray([[4, 7], [2, 6]])
      const I = A.multiply(A.inverse())
      expect(I.get(0, 0)).toBeCloseTo(1, 10)
      expect(I.get(0, 1)).toBeCloseTo(0, 10)
      expect(I.get(1, 0)).toBeCloseTo(0, 10)
      expect(I.get(1, 1)).toBeCloseTo(1, 10)
    })

    it('throws on singular matrix', () => {
      const singular = Matrix.fromArray([[1, 2], [2, 4]])
      expect(() => singular.inverse()).toThrow()
    })
  })

  describe('cholesky', () => {
    it('decomposes SPD matrix (R ground truth)', () => {
      const A = Matrix.fromArray(gt.matrix.cholesky_3x3.A as number[][])
      const L = A.cholesky()
      const expected = gt.matrix.cholesky_3x3.L as number[][]
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(L.get(i, j)).toBeCloseTo(expected[i]![j]!, 8)
        }
      }
    })

    it('L * L^T ≈ A', () => {
      const A = Matrix.fromArray([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])
      const L = A.cholesky()
      const LLt = L.multiply(L.transpose())
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(LLt.get(i, j)).toBeCloseTo(A.get(i, j), 8)
        }
      }
    })

    it('throws on non-SPD matrix', () => {
      const nonSPD = Matrix.fromArray([[1, 2], [2, 1]])  // not positive definite
      expect(() => nonSPD.cholesky()).toThrow()
    })
  })

  describe('logDet', () => {
    it('computes log determinant of SPD matrix', () => {
      // det(A) = 4*6 - 7*2 = 10, logDet = log(10) ≈ 2.302585
      const A = Matrix.fromArray([[4, 7], [2, 6]])
      // This is not SPD (not symmetric), so use a symmetric one
      const B = Matrix.fromArray([[4, 2], [2, 3]])
      // det = 4*3 - 2*2 = 8, logDet = log(8) ≈ 2.079442
      expect(B.logDet()).toBeCloseTo(Math.log(8), 5)
    })
  })

  describe('SVD', () => {
    it('A ≈ U * diag(S) * V^T', () => {
      const A = Matrix.fromArray([[1, 2], [3, 4], [5, 6]])
      const { U, S, V } = A.svd()
      // Reconstruct: U * diag(S) * V^T
      const Diag = Matrix.fromArray(
        Array.from({ length: U.cols }, (_, i) =>
          Array.from({ length: V.cols }, (_, j) => (i === j ? (S[i] ?? 0) : 0))
        )
      )
      const reconstructed = U.multiply(Diag).multiply(V.transpose())
      for (let i = 0; i < A.rows; i++) {
        for (let j = 0; j < A.cols; j++) {
          expect(reconstructed.get(i, j)).toBeCloseTo(A.get(i, j), 4)
        }
      }
    })

    it('singular values are non-negative and descending', () => {
      const A = Matrix.fromArray([[1, 2, 3], [4, 5, 6]])
      const { S } = A.svd()
      expect(S.length).toBe(2)
      for (const s of S) expect(s).toBeGreaterThanOrEqual(0)
      expect(S[0]!).toBeGreaterThanOrEqual(S[1]!)
    })
  })

  describe('eigen', () => {
    it('computes eigenvalues of symmetric matrix', () => {
      // A = [[2, 1], [1, 2]], eigenvalues = 3, 1
      const A = Matrix.fromArray([[2, 1], [1, 2]])
      const { values } = A.eigen()
      expect(values[0]).toBeCloseTo(3, 5)
      expect(values[1]).toBeCloseTo(1, 5)
    })
  })

  describe('solveLinear', () => {
    it('solves Ax = b', () => {
      // A x = b → x
      const A = Matrix.fromArray([[2, 1], [5, 3]])
      const b = [4, 7]
      const x = solveLinear(A, b)
      // x = [5, -6]
      expect(x[0]).toBeCloseTo(5, 5)
      expect(x[1]).toBeCloseTo(-6, 5)
    })
  })
})
