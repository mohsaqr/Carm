/**
 * Matrix class for Carm.
 * Implements multiply, transpose, inverse (via Cholesky or LU), log-determinant, and SVD.
 * Pure computation — no DOM, no D3, no side effects.
 */

export class Matrix {
  readonly rows: number
  readonly cols: number
  private readonly _data: readonly number[]  // row-major flat array

  constructor(rows: number, cols: number, data?: readonly number[]) {
    this.rows = rows
    this.cols = cols
    this._data = data ?? new Array<number>(rows * cols).fill(0)
    if (this._data.length !== rows * cols) {
      throw new Error(`Matrix data length ${this._data.length} does not match ${rows}×${cols}`)
    }
  }

  /** Build from 2-D array (row-major). */
  static fromArray(arr: readonly (readonly number[])[]): Matrix {
    const rows = arr.length
    if (rows === 0) throw new Error('Matrix cannot have 0 rows')
    const cols = arr[0]!.length
    const data: number[] = []
    for (const row of arr) {
      if (row.length !== cols) throw new Error('All rows must have equal length')
      for (const v of row) data.push(v)
    }
    return new Matrix(rows, cols, data)
  }

  /** Identity matrix of size n. */
  static identity(n: number): Matrix {
    const data = new Array<number>(n * n).fill(0)
    for (let i = 0; i < n; i++) data[i * n + i] = 1
    return new Matrix(n, n, data)
  }

  /** Zero matrix. */
  static zeros(rows: number, cols: number): Matrix {
    return new Matrix(rows, cols, new Array<number>(rows * cols).fill(0))
  }

  /** Get element at (i, j) — 0-indexed. */
  get(i: number, j: number): number {
    const v = this._data[i * this.cols + j]
    if (v === undefined) throw new Error(`Index (${i},${j}) out of bounds for ${this.rows}×${this.cols}`)
    return v
  }

  /** Return 2-D array representation. */
  toArray(): number[][] {
    return Array.from({ length: this.rows }, (_, i) =>
      Array.from({ length: this.cols }, (_, j) => this.get(i, j))
    )
  }

  /** Return flat row-major copy. */
  toFlat(): number[] {
    return [...this._data]
  }

  /** Matrix transpose. */
  transpose(): Matrix {
    const data: number[] = new Array(this.rows * this.cols)
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        data[j * this.rows + i] = this.get(i, j)
      }
    }
    return new Matrix(this.cols, this.rows, data)
  }

  /** Matrix multiplication: this × other. */
  multiply(other: Matrix): Matrix {
    if (this.cols !== other.rows) {
      throw new Error(`Dimension mismatch: ${this.rows}×${this.cols} × ${other.rows}×${other.cols}`)
    }
    const data: number[] = new Array(this.rows * other.cols).fill(0)
    for (let i = 0; i < this.rows; i++) {
      for (let k = 0; k < this.cols; k++) {
        const aik = this.get(i, k)
        for (let j = 0; j < other.cols; j++) {
          data[i * other.cols + j]! += aik * other.get(k, j)
        }
      }
    }
    return new Matrix(this.rows, other.cols, data)
  }

  /** Scalar multiplication. */
  scale(s: number): Matrix {
    return new Matrix(this.rows, this.cols, this._data.map(v => v * s))
  }

  /** Element-wise add. */
  add(other: Matrix): Matrix {
    if (this.rows !== other.rows || this.cols !== other.cols) {
      throw new Error('Matrix dimensions must match for addition')
    }
    return new Matrix(this.rows, this.cols, this._data.map((v, i) => v + (other._data[i] ?? 0)))
  }

  /** Element-wise subtract. */
  subtract(other: Matrix): Matrix {
    if (this.rows !== other.rows || this.cols !== other.cols) {
      throw new Error('Matrix dimensions must match for subtraction')
    }
    return new Matrix(this.rows, this.cols, this._data.map((v, i) => v - (other._data[i] ?? 0)))
  }

  /**
   * Cholesky decomposition for symmetric positive-definite matrices.
   * Returns lower-triangular L such that this = L * L^T.
   * Throws if matrix is not SPD.
   */
  cholesky(): Matrix {
    if (this.rows !== this.cols) throw new Error('Cholesky requires square matrix')
    const n = this.rows
    const L: number[] = new Array(n * n).fill(0)

    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = this.get(i, j)
        for (let k = 0; k < j; k++) {
          sum -= (L[i * n + k] ?? 0) * (L[j * n + k] ?? 0)
        }
        if (i === j) {
          if (sum <= 0) throw new Error(`Matrix is not positive definite (diagonal became ${sum} at position ${i})`)
          L[i * n + j] = Math.sqrt(sum)
        } else {
          const diag = L[j * n + j]
          if (!diag || diag === 0) throw new Error('Zero diagonal in Cholesky')
          L[i * n + j] = sum / diag
        }
      }
    }
    return new Matrix(n, n, L)
  }

  /**
   * Log-determinant via Cholesky: log|A| = 2 * Σ log(L_ii).
   * Only valid for symmetric positive-definite matrices.
   */
  logDet(): number {
    const L = this.cholesky()
    let logdet = 0
    for (let i = 0; i < this.rows; i++) {
      const diag = L.get(i, i)
      if (diag <= 0) throw new Error('Non-positive diagonal in Cholesky')
      logdet += Math.log(diag)
    }
    return 2 * logdet
  }

  /**
   * Inverse via LU decomposition with partial pivoting.
   * Works for any non-singular square matrix.
   * Formula: Doolittle LU, then forward/back substitution for each column of I.
   */
  inverse(): Matrix {
    if (this.rows !== this.cols) throw new Error('Inverse requires square matrix')
    const n = this.rows
    // Build augmented matrix [A | I]
    const aug: number[] = new Array(n * 2 * n)
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        aug[i * 2 * n + j] = this.get(i, j)
        aug[i * 2 * n + n + j] = i === j ? 1 : 0
      }
    }

    // Gaussian elimination with partial pivoting
    for (let col = 0; col < n; col++) {
      // Find pivot
      let maxRow = col
      let maxVal = Math.abs(aug[col * 2 * n + col] ?? 0)
      for (let row = col + 1; row < n; row++) {
        const v = Math.abs(aug[row * 2 * n + col] ?? 0)
        if (v > maxVal) { maxVal = v; maxRow = row }
      }
      if (maxVal < 1e-12) throw new Error('Matrix is singular or near-singular')

      // Swap rows
      if (maxRow !== col) {
        for (let j = 0; j < 2 * n; j++) {
          const tmp = aug[col * 2 * n + j] ?? 0
          aug[col * 2 * n + j] = aug[maxRow * 2 * n + j] ?? 0
          aug[maxRow * 2 * n + j] = tmp
        }
      }

      const pivot = aug[col * 2 * n + col] ?? 0
      // Eliminate below
      for (let row = 0; row < n; row++) {
        if (row === col) continue
        const factor = (aug[row * 2 * n + col] ?? 0) / pivot
        for (let j = 0; j < 2 * n; j++) {
          aug[row * 2 * n + j]! -= factor * (aug[col * 2 * n + j] ?? 0)
        }
      }
      // Scale pivot row
      for (let j = 0; j < 2 * n; j++) {
        aug[col * 2 * n + j]! /= pivot
      }
    }

    // Extract right half (inverse)
    const inv: number[] = new Array(n * n)
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        inv[i * n + j] = aug[i * 2 * n + n + j] ?? 0
      }
    }
    return new Matrix(n, n, inv)
  }

  /**
   * Singular Value Decomposition: A = U · S · V^T
   * Returns { U, S (diagonal values), V }.
   * Algorithm: Golub-Reinsch (one-sided Jacobi for small matrices).
   * Reference: Golub & Van Loan, "Matrix Computations", 4th ed., Algorithm 8.6.2
   */
  svd(): { U: Matrix; S: number[]; V: Matrix } {
    // Use Jacobi one-sided SVD
    const m = this.rows
    const n = this.cols
    // Work on a copy of A as column matrix V^T = I
    const a = this.toArray()
    const v: number[][] = Array.from({ length: n }, (_, i) =>
      Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
    )

    const MAX_ITER = 200 * n * n
    for (let iter = 0; iter < MAX_ITER; iter++) {
      let converged = true
      for (let p = 0; p < n - 1; p++) {
        for (let q = p + 1; q < n; q++) {
          // Compute alpha, beta, gamma for Jacobi rotation
          let alpha = 0, beta = 0, gamma = 0
          for (let i = 0; i < m; i++) {
            alpha += (a[i]![p] ?? 0) ** 2
            beta  += (a[i]![q] ?? 0) ** 2
            gamma += (a[i]![p] ?? 0) * (a[i]![q] ?? 0)
          }
          if (Math.abs(gamma) < 1e-15 * Math.sqrt(alpha * beta)) continue
          converged = false
          const zeta = (beta - alpha) / (2 * gamma)
          const t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta))
          const c = 1 / Math.sqrt(1 + t * t)
          const s = t * c

          // Update columns p and q of a
          for (let i = 0; i < m; i++) {
            const ap = a[i]![p] ?? 0
            const aq = a[i]![q] ?? 0
            a[i]![p] = c * ap + s * aq
            a[i]![q] = -s * ap + c * aq
          }
          // Update V
          for (let i = 0; i < n; i++) {
            const vp = v[i]![p] ?? 0
            const vq = v[i]![q] ?? 0
            v[i]![p] = c * vp + s * vq
            v[i]![q] = -s * vp + c * vq
          }
        }
      }
      if (converged) break
    }

    // Compute singular values and normalize columns of a → U
    const singularValues = Array.from({ length: n }, (_, j) => {
      let sum = 0
      for (let i = 0; i < m; i++) sum += (a[i]![j] ?? 0) ** 2
      return Math.sqrt(sum)
    })

    // Build U (m×n), normalize each column
    const uData: number[][] = Array.from({ length: m }, () => new Array<number>(n).fill(0))
    for (let j = 0; j < n; j++) {
      const sv = singularValues[j] ?? 0
      for (let i = 0; i < m; i++) {
        uData[i]![j] = sv > 1e-15 ? (a[i]![j] ?? 0) / sv : 0
      }
    }

    // Sort by descending singular value, truncate to k = min(m, n)
    const order = singularValues.map((_, i) => i).sort((a, b) => (singularValues[b] ?? 0) - (singularValues[a] ?? 0))
    const k = Math.min(m, n)
    const S = order.slice(0, k).map(i => singularValues[i] ?? 0)
    const Uarr: number[][] = Array.from({ length: m }, (_, i) => order.slice(0, k).map(j => uData[i]![j] ?? 0))
    const Varr: number[][] = Array.from({ length: n }, (_, i) => order.slice(0, k).map(j => v[i]![j] ?? 0))

    return {
      U: Matrix.fromArray(Uarr),
      S,
      V: Matrix.fromArray(Varr),
    }
  }

  /**
   * Pseudo-inverse via SVD: A+ = V · S^{-1} · U^T
   */
  pseudoInverse(tol = 1e-10): Matrix {
    const { U, S, V } = this.svd()
    const maxS = Math.max(...S)
    const threshold = tol * maxS
    const SInv = S.map(s => (s > threshold ? 1 / s : 0))

    // V * diag(SInv)
    const VSInv = Matrix.fromArray(
      Array.from({ length: V.rows }, (_, i) =>
        Array.from({ length: V.cols }, (_, j) => V.get(i, j) * (SInv[j] ?? 0))
      )
    )
    return VSInv.multiply(U.transpose())
  }

  /** Trace (sum of diagonal elements). */
  trace(): number {
    const n = Math.min(this.rows, this.cols)
    let t = 0
    for (let i = 0; i < n; i++) t += this.get(i, i)
    return t
  }

  /** Extract diagonal as array. */
  diagonal(): number[] {
    const n = Math.min(this.rows, this.cols)
    return Array.from({ length: n }, (_, i) => this.get(i, i))
  }

  /** Column vector as Matrix from array. */
  static colVec(data: readonly number[]): Matrix {
    return new Matrix(data.length, 1, [...data])
  }

  /** Row vector as Matrix from array. */
  static rowVec(data: readonly number[]): Matrix {
    return new Matrix(1, data.length, [...data])
  }

  /**
   * Eigendecomposition for symmetric matrices via Jacobi iterations.
   * Returns { values, vectors } where vectors are columns of the eigenvector matrix.
   * Reference: Golub & Van Loan, Algorithm 8.4.1
   */
  eigen(): { values: number[]; vectors: Matrix } {
    if (this.rows !== this.cols) throw new Error('eigen() requires square matrix')
    const n = this.rows
    const a = this.toArray()
    const q: number[][] = Array.from({ length: n }, (_, i) =>
      Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
    )

    for (let iter = 0; iter < 100 * n * n; iter++) {
      // Find largest off-diagonal element
      let p = 0, r = 1
      let maxVal = Math.abs(a[0]![1] ?? 0)
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          if (Math.abs(a[i]![j] ?? 0) > maxVal) {
            maxVal = Math.abs(a[i]![j] ?? 0)
            p = i; r = j
          }
        }
      }
      if (maxVal < 1e-12) break

      const app = a[p]![p] ?? 0, arr = a[r]![r] ?? 0, apr = a[p]![r] ?? 0
      const theta = 0.5 * Math.atan2(2 * apr, app - arr)
      const c = Math.cos(theta), s = Math.sin(theta)

      // Update a
      const newApp = c * c * app + 2 * s * c * apr + s * s * arr
      const newArr = s * s * app - 2 * s * c * apr + c * c * arr
      a[p]![p] = newApp; a[r]![r] = newArr; a[p]![r] = 0; a[r]![p] = 0

      for (let i = 0; i < n; i++) {
        if (i === p || i === r) continue
        const aip = a[i]![p] ?? 0, air = a[i]![r] ?? 0
        a[i]![p] = c * aip + s * air
        a[p]![i] = a[i]![p]!
        a[i]![r] = -s * aip + c * air
        a[r]![i] = a[i]![r]!
      }

      // Update eigenvectors
      for (let i = 0; i < n; i++) {
        const qip = q[i]![p] ?? 0, qir = q[i]![r] ?? 0
        q[i]![p] = c * qip + s * qir
        q[i]![r] = -s * qip + c * qir
      }
    }

    const values = Array.from({ length: n }, (_, i) => a[i]![i] ?? 0)
    // Sort descending
    const order = values.map((_, i) => i).sort((a, b) => values[b]! - values[a]!)
    return {
      values: order.map(i => values[i]!),
      vectors: Matrix.fromArray(Array.from({ length: n }, (_, i) => order.map(j => q[i]![j] ?? 0))),
    }
  }
}

/** Solve linear system A·x = b using the (already computed) inverse.  */
export function solveLinear(A: Matrix, b: readonly number[]): number[] {
  const inv = A.inverse()
  const x = inv.multiply(Matrix.colVec(b))
  return Array.from({ length: b.length }, (_, i) => x.get(i, 0))
}
