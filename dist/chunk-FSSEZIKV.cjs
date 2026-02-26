'use strict';

// src/core/matrix.ts
var Matrix = class _Matrix {
  rows;
  cols;
  _data;
  // row-major flat array
  constructor(rows, cols, data) {
    this.rows = rows;
    this.cols = cols;
    this._data = data ?? new Array(rows * cols).fill(0);
    if (this._data.length !== rows * cols) {
      throw new Error(`Matrix data length ${this._data.length} does not match ${rows}\xD7${cols}`);
    }
  }
  /** Build from 2-D array (row-major). */
  static fromArray(arr) {
    const rows = arr.length;
    if (rows === 0) throw new Error("Matrix cannot have 0 rows");
    const cols = arr[0].length;
    const data = [];
    for (const row of arr) {
      if (row.length !== cols) throw new Error("All rows must have equal length");
      for (const v of row) data.push(v);
    }
    return new _Matrix(rows, cols, data);
  }
  /** Identity matrix of size n. */
  static identity(n) {
    const data = new Array(n * n).fill(0);
    for (let i = 0; i < n; i++) data[i * n + i] = 1;
    return new _Matrix(n, n, data);
  }
  /** Zero matrix. */
  static zeros(rows, cols) {
    return new _Matrix(rows, cols, new Array(rows * cols).fill(0));
  }
  /** Get element at (i, j) — 0-indexed. */
  get(i, j) {
    const v = this._data[i * this.cols + j];
    if (v === void 0) throw new Error(`Index (${i},${j}) out of bounds for ${this.rows}\xD7${this.cols}`);
    return v;
  }
  /** Return 2-D array representation. */
  toArray() {
    return Array.from(
      { length: this.rows },
      (_, i) => Array.from({ length: this.cols }, (_2, j) => this.get(i, j))
    );
  }
  /** Return flat row-major copy. */
  toFlat() {
    return [...this._data];
  }
  /** Matrix transpose. */
  transpose() {
    const data = new Array(this.rows * this.cols);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        data[j * this.rows + i] = this.get(i, j);
      }
    }
    return new _Matrix(this.cols, this.rows, data);
  }
  /** Matrix multiplication: this × other. */
  multiply(other) {
    if (this.cols !== other.rows) {
      throw new Error(`Dimension mismatch: ${this.rows}\xD7${this.cols} \xD7 ${other.rows}\xD7${other.cols}`);
    }
    const data = new Array(this.rows * other.cols).fill(0);
    for (let i = 0; i < this.rows; i++) {
      for (let k = 0; k < this.cols; k++) {
        const aik = this.get(i, k);
        for (let j = 0; j < other.cols; j++) {
          data[i * other.cols + j] += aik * other.get(k, j);
        }
      }
    }
    return new _Matrix(this.rows, other.cols, data);
  }
  /** Scalar multiplication. */
  scale(s) {
    return new _Matrix(this.rows, this.cols, this._data.map((v) => v * s));
  }
  /** Element-wise add. */
  add(other) {
    if (this.rows !== other.rows || this.cols !== other.cols) {
      throw new Error("Matrix dimensions must match for addition");
    }
    return new _Matrix(this.rows, this.cols, this._data.map((v, i) => v + (other._data[i] ?? 0)));
  }
  /** Element-wise subtract. */
  subtract(other) {
    if (this.rows !== other.rows || this.cols !== other.cols) {
      throw new Error("Matrix dimensions must match for subtraction");
    }
    return new _Matrix(this.rows, this.cols, this._data.map((v, i) => v - (other._data[i] ?? 0)));
  }
  /**
   * Cholesky decomposition for symmetric positive-definite matrices.
   * Returns lower-triangular L such that this = L * L^T.
   * Throws if matrix is not SPD.
   */
  cholesky() {
    if (this.rows !== this.cols) throw new Error("Cholesky requires square matrix");
    const n = this.rows;
    const L = new Array(n * n).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = this.get(i, j);
        for (let k = 0; k < j; k++) {
          sum -= (L[i * n + k] ?? 0) * (L[j * n + k] ?? 0);
        }
        if (i === j) {
          if (sum <= 0) throw new Error(`Matrix is not positive definite (diagonal became ${sum} at position ${i})`);
          L[i * n + j] = Math.sqrt(sum);
        } else {
          const diag = L[j * n + j];
          if (!diag || diag === 0) throw new Error("Zero diagonal in Cholesky");
          L[i * n + j] = sum / diag;
        }
      }
    }
    return new _Matrix(n, n, L);
  }
  /**
   * Log-determinant via Cholesky: log|A| = 2 * Σ log(L_ii).
   * Only valid for symmetric positive-definite matrices.
   */
  logDet() {
    const L = this.cholesky();
    let logdet = 0;
    for (let i = 0; i < this.rows; i++) {
      const diag = L.get(i, i);
      if (diag <= 0) throw new Error("Non-positive diagonal in Cholesky");
      logdet += Math.log(diag);
    }
    return 2 * logdet;
  }
  /**
   * Inverse via LU decomposition with partial pivoting.
   * Works for any non-singular square matrix.
   * Formula: Doolittle LU, then forward/back substitution for each column of I.
   */
  inverse() {
    if (this.rows !== this.cols) throw new Error("Inverse requires square matrix");
    const n = this.rows;
    const aug = new Array(n * 2 * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        aug[i * 2 * n + j] = this.get(i, j);
        aug[i * 2 * n + n + j] = i === j ? 1 : 0;
      }
    }
    for (let col = 0; col < n; col++) {
      let maxRow = col;
      let maxVal = Math.abs(aug[col * 2 * n + col] ?? 0);
      for (let row = col + 1; row < n; row++) {
        const v = Math.abs(aug[row * 2 * n + col] ?? 0);
        if (v > maxVal) {
          maxVal = v;
          maxRow = row;
        }
      }
      if (maxVal < 1e-12) throw new Error("Matrix is singular or near-singular");
      if (maxRow !== col) {
        for (let j = 0; j < 2 * n; j++) {
          const tmp = aug[col * 2 * n + j] ?? 0;
          aug[col * 2 * n + j] = aug[maxRow * 2 * n + j] ?? 0;
          aug[maxRow * 2 * n + j] = tmp;
        }
      }
      const pivot = aug[col * 2 * n + col] ?? 0;
      for (let row = 0; row < n; row++) {
        if (row === col) continue;
        const factor = (aug[row * 2 * n + col] ?? 0) / pivot;
        for (let j = 0; j < 2 * n; j++) {
          aug[row * 2 * n + j] -= factor * (aug[col * 2 * n + j] ?? 0);
        }
      }
      for (let j = 0; j < 2 * n; j++) {
        aug[col * 2 * n + j] /= pivot;
      }
    }
    const inv = new Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        inv[i * n + j] = aug[i * 2 * n + n + j] ?? 0;
      }
    }
    return new _Matrix(n, n, inv);
  }
  /**
   * Singular Value Decomposition: A = U · S · V^T
   * Returns { U, S (diagonal values), V }.
   * Algorithm: Golub-Reinsch (one-sided Jacobi for small matrices).
   * Reference: Golub & Van Loan, "Matrix Computations", 4th ed., Algorithm 8.6.2
   */
  svd() {
    const m = this.rows;
    const n = this.cols;
    const a = this.toArray();
    const v = Array.from(
      { length: n },
      (_, i) => Array.from({ length: n }, (_2, j) => i === j ? 1 : 0)
    );
    const MAX_ITER = 200 * n * n;
    for (let iter = 0; iter < MAX_ITER; iter++) {
      let converged = true;
      for (let p = 0; p < n - 1; p++) {
        for (let q = p + 1; q < n; q++) {
          let alpha = 0, beta = 0, gamma = 0;
          for (let i = 0; i < m; i++) {
            alpha += (a[i][p] ?? 0) ** 2;
            beta += (a[i][q] ?? 0) ** 2;
            gamma += (a[i][p] ?? 0) * (a[i][q] ?? 0);
          }
          if (Math.abs(gamma) < 1e-15 * Math.sqrt(alpha * beta)) continue;
          converged = false;
          const zeta = (beta - alpha) / (2 * gamma);
          const t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
          const c = 1 / Math.sqrt(1 + t * t);
          const s = t * c;
          for (let i = 0; i < m; i++) {
            const ap = a[i][p] ?? 0;
            const aq = a[i][q] ?? 0;
            a[i][p] = c * ap + s * aq;
            a[i][q] = -s * ap + c * aq;
          }
          for (let i = 0; i < n; i++) {
            const vp = v[i][p] ?? 0;
            const vq = v[i][q] ?? 0;
            v[i][p] = c * vp + s * vq;
            v[i][q] = -s * vp + c * vq;
          }
        }
      }
      if (converged) break;
    }
    const singularValues = Array.from({ length: n }, (_, j) => {
      let sum = 0;
      for (let i = 0; i < m; i++) sum += (a[i][j] ?? 0) ** 2;
      return Math.sqrt(sum);
    });
    const uData = Array.from({ length: m }, () => new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const sv = singularValues[j] ?? 0;
      for (let i = 0; i < m; i++) {
        uData[i][j] = sv > 1e-15 ? (a[i][j] ?? 0) / sv : 0;
      }
    }
    const order = singularValues.map((_, i) => i).sort((a2, b) => (singularValues[b] ?? 0) - (singularValues[a2] ?? 0));
    const k = Math.min(m, n);
    const S = order.slice(0, k).map((i) => singularValues[i] ?? 0);
    const Uarr = Array.from({ length: m }, (_, i) => order.slice(0, k).map((j) => uData[i][j] ?? 0));
    const Varr = Array.from({ length: n }, (_, i) => order.slice(0, k).map((j) => v[i][j] ?? 0));
    return {
      U: _Matrix.fromArray(Uarr),
      S,
      V: _Matrix.fromArray(Varr)
    };
  }
  /**
   * Pseudo-inverse via SVD: A+ = V · S^{-1} · U^T
   */
  pseudoInverse(tol = 1e-10) {
    const { U, S, V } = this.svd();
    const maxS = Math.max(...S);
    const threshold = tol * maxS;
    const SInv = S.map((s) => s > threshold ? 1 / s : 0);
    const VSInv = _Matrix.fromArray(
      Array.from(
        { length: V.rows },
        (_, i) => Array.from({ length: V.cols }, (_2, j) => V.get(i, j) * (SInv[j] ?? 0))
      )
    );
    return VSInv.multiply(U.transpose());
  }
  /** Trace (sum of diagonal elements). */
  trace() {
    const n = Math.min(this.rows, this.cols);
    let t = 0;
    for (let i = 0; i < n; i++) t += this.get(i, i);
    return t;
  }
  /** Extract diagonal as array. */
  diagonal() {
    const n = Math.min(this.rows, this.cols);
    return Array.from({ length: n }, (_, i) => this.get(i, i));
  }
  /** Column vector as Matrix from array. */
  static colVec(data) {
    return new _Matrix(data.length, 1, [...data]);
  }
  /** Row vector as Matrix from array. */
  static rowVec(data) {
    return new _Matrix(1, data.length, [...data]);
  }
  /**
   * Eigendecomposition for symmetric matrices via Jacobi iterations.
   * Returns { values, vectors } where vectors are columns of the eigenvector matrix.
   * Reference: Golub & Van Loan, Algorithm 8.4.1
   */
  eigen() {
    if (this.rows !== this.cols) throw new Error("eigen() requires square matrix");
    const n = this.rows;
    const a = this.toArray();
    const q = Array.from(
      { length: n },
      (_, i) => Array.from({ length: n }, (_2, j) => i === j ? 1 : 0)
    );
    for (let iter = 0; iter < 100 * n * n; iter++) {
      let p = 0, r = 1;
      let maxVal = Math.abs(a[0][1] ?? 0);
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          if (Math.abs(a[i][j] ?? 0) > maxVal) {
            maxVal = Math.abs(a[i][j] ?? 0);
            p = i;
            r = j;
          }
        }
      }
      if (maxVal < 1e-12) break;
      const app = a[p][p] ?? 0, arr = a[r][r] ?? 0, apr = a[p][r] ?? 0;
      const theta = 0.5 * Math.atan2(2 * apr, app - arr);
      const c = Math.cos(theta), s = Math.sin(theta);
      const newApp = c * c * app + 2 * s * c * apr + s * s * arr;
      const newArr = s * s * app - 2 * s * c * apr + c * c * arr;
      a[p][p] = newApp;
      a[r][r] = newArr;
      a[p][r] = 0;
      a[r][p] = 0;
      for (let i = 0; i < n; i++) {
        if (i === p || i === r) continue;
        const aip = a[i][p] ?? 0, air = a[i][r] ?? 0;
        a[i][p] = c * aip + s * air;
        a[p][i] = a[i][p];
        a[i][r] = -s * aip + c * air;
        a[r][i] = a[i][r];
      }
      for (let i = 0; i < n; i++) {
        const qip = q[i][p] ?? 0, qir = q[i][r] ?? 0;
        q[i][p] = c * qip + s * qir;
        q[i][r] = -s * qip + c * qir;
      }
    }
    const values = Array.from({ length: n }, (_, i) => a[i][i] ?? 0);
    const order = values.map((_, i) => i).sort((a2, b) => values[b] - values[a2]);
    return {
      values: order.map((i) => values[i]),
      vectors: _Matrix.fromArray(Array.from({ length: n }, (_, i) => order.map((j) => q[i][j] ?? 0)))
    };
  }
};
function solveLinear(A, b) {
  const inv = A.inverse();
  const x = inv.multiply(Matrix.colVec(b));
  return Array.from({ length: b.length }, (_, i) => x.get(i, 0));
}

exports.Matrix = Matrix;
exports.solveLinear = solveLinear;
//# sourceMappingURL=chunk-FSSEZIKV.cjs.map
//# sourceMappingURL=chunk-FSSEZIKV.cjs.map