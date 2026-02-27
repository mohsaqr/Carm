"use strict";
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};
var __copyProps = (to, from, except, desc) => {
  if (from && typeof from === "object" || typeof from === "function") {
    for (let key of __getOwnPropNames(from))
      if (!__hasOwnProp.call(to, key) && key !== except)
        __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
  }
  return to;
};
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// src/core/index.ts
var core_exports = {};
__export(core_exports, {
  Matrix: () => Matrix,
  adjustPValues: () => adjustPValues,
  betaFn: () => betaFn,
  chiSqCDF: () => chiSqCDF,
  chiSqPValue: () => chiSqPValue,
  chiSqQuantile: () => chiSqQuantile,
  clamp: () => clamp,
  cov: () => cov,
  fDistCDF: () => fDistCDF,
  fDistPValue: () => fDistPValue,
  formatANOVA: () => formatANOVA,
  formatCFAFit: () => formatCFAFit,
  formatCI: () => formatCI,
  formatChiSq: () => formatChiSq,
  formatCorrelation: () => formatCorrelation,
  formatDF: () => formatDF,
  formatKruskalWallis: () => formatKruskalWallis,
  formatLMM: () => formatLMM,
  formatMannWhitney: () => formatMannWhitney,
  formatP: () => formatP,
  formatRMANOVA: () => formatRMANOVA,
  formatRegression: () => formatRegression,
  formatStat: () => formatStat,
  formatTTest: () => formatTTest,
  gamma: () => gamma,
  incompleteBeta: () => incompleteBeta,
  incompleteGamma: () => incompleteGamma,
  interpretCohensD: () => interpretCohensD,
  interpretCramerV: () => interpretCramerV,
  interpretEffect: () => interpretEffect,
  interpretEtaSq: () => interpretEtaSq,
  interpretR: () => interpretR,
  logBeta: () => logBeta,
  logGamma: () => logGamma,
  mean: () => mean,
  median: () => median,
  nelderMead: () => nelderMead,
  normalCDF: () => normalCDF,
  normalQuantile: () => normalQuantile,
  normalSurvival: () => normalSurvival,
  pKendallExact: () => pKendallExact,
  pSpearmanExact: () => pSpearmanExact,
  pValueStudentizedRangeApprox: () => pValueStudentizedRangeApprox,
  ptukeyApprox: () => ptukeyApprox,
  quantile: () => quantile,
  rank: () => rank,
  roundTo: () => roundTo,
  sd: () => sd,
  se: () => se,
  solveLinear: () => solveLinear,
  sortAsc: () => sortAsc,
  ss: () => ss,
  tDistCDF: () => tDistCDF,
  tDistPValue: () => tDistPValue,
  tDistQuantile: () => tDistQuantile,
  variance: () => variance
});
module.exports = __toCommonJS(core_exports);

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
          let alpha = 0, beta = 0, gamma2 = 0;
          for (let i = 0; i < m; i++) {
            alpha += (a[i][p] ?? 0) ** 2;
            beta += (a[i][q] ?? 0) ** 2;
            gamma2 += (a[i][p] ?? 0) * (a[i][q] ?? 0);
          }
          if (Math.abs(gamma2) < 1e-15 * Math.sqrt(alpha * beta)) continue;
          converged = false;
          const zeta = (beta - alpha) / (2 * gamma2);
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

// src/core/math.ts
function logGamma(z) {
  if (z <= 0) throw new Error(`logGamma: z must be positive, got ${z}`);
  if (z < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
  }
  const g = 7;
  const c = [
    0.9999999999998099,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9984369578019572e-21,
    15056327351493116e-23
  ];
  const x = z - 1;
  let sum = c[0];
  for (let i = 1; i < g + 2; i++) sum += (c[i] ?? 0) / (x + i);
  const t = x + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum);
}
function gamma(z) {
  return Math.exp(logGamma(z));
}
function logBeta(a, b) {
  return logGamma(a) + logGamma(b) - logGamma(a + b);
}
function betaFn(a, b) {
  return Math.exp(logBeta(a, b));
}
function incompleteBeta(x, a, b) {
  if (x < 0 || x > 1) throw new Error(`incompleteBeta: x must be in [0,1], got ${x}`);
  if (x === 0) return 0;
  if (x === 1) return 1;
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - incompleteBeta(1 - x, b, a);
  }
  const lbeta_ab = logBeta(a, b);
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lbeta_ab) / a;
  const MAX_ITER = 200;
  const EPS = 3e-7;
  const FPMIN = 1e-30;
  let f = 1;
  let C = 1;
  let D = 1 - (a + b) * x / (a + 1);
  if (Math.abs(D) < FPMIN) D = FPMIN;
  D = 1 / D;
  f = D;
  for (let m = 1; m <= MAX_ITER; m++) {
    let num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
    D = 1 + num * D;
    if (Math.abs(D) < FPMIN) D = FPMIN;
    C = 1 + num / C;
    if (Math.abs(C) < FPMIN) C = FPMIN;
    D = 1 / D;
    f *= C * D;
    num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
    D = 1 + num * D;
    if (Math.abs(D) < FPMIN) D = FPMIN;
    C = 1 + num / C;
    if (Math.abs(C) < FPMIN) C = FPMIN;
    D = 1 / D;
    const delta = C * D;
    f *= delta;
    if (Math.abs(delta - 1) < EPS) break;
  }
  return front * f;
}
function incompleteGamma(a, x) {
  if (x < 0) throw new Error(`incompleteGamma: x must be \u2265 0, got ${x}`);
  if (x === 0) return 0;
  if (x < a + 1) {
    let term = 1 / a;
    let sum = term;
    for (let n = 1; n < 200; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 3e-7) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
  } else {
    return 1 - incompleteGammaComplement(a, x);
  }
}
function incompleteGammaComplement(a, x) {
  const FPMIN = 1e-30;
  let f = x + 1 - a;
  if (Math.abs(f) < FPMIN) f = FPMIN;
  let C = f, D = 0;
  for (let i = 1; i <= 200; i++) {
    const an = -i * (i - a);
    const bn = x + 2 * i + 1 - a;
    D = bn + an * D;
    if (Math.abs(D) < FPMIN) D = FPMIN;
    C = bn + an / C;
    if (Math.abs(C) < FPMIN) C = FPMIN;
    D = 1 / D;
    const delta = C * D;
    f *= delta;
    if (Math.abs(delta - 1) < 3e-7) break;
  }
  return Math.exp(-x + a * Math.log(x) - logGamma(a)) / f;
}
function tDistPValue(t, df) {
  const x = df / (df + t * t);
  const p = incompleteBeta(x, df / 2, 0.5);
  return p;
}
function tDistCDF(t, df) {
  const x = df / (df + t * t);
  const p = 0.5 * incompleteBeta(x, df / 2, 0.5);
  return t >= 0 ? 1 - p : p;
}
function tDistQuantile(p, df) {
  if (p <= 0 || p >= 1) throw new Error(`tDistQuantile: p must be in (0,1), got ${p}`);
  let lo = -50, hi = 50;
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (tDistCDF(mid, df) < p) lo = mid;
    else hi = mid;
    if (hi - lo < 1e-10) break;
  }
  return (lo + hi) / 2;
}
function fDistCDF(f, df1, df2) {
  if (f <= 0) return 0;
  const x = df1 * f / (df1 * f + df2);
  return incompleteBeta(x, df1 / 2, df2 / 2);
}
function fDistPValue(f, df1, df2) {
  return 1 - fDistCDF(f, df1, df2);
}
function chiSqCDF(x, df) {
  if (x <= 0) return 0;
  return incompleteGamma(df / 2, x / 2);
}
function chiSqPValue(x, df) {
  return 1 - chiSqCDF(x, df);
}
function chiSqQuantile(p, df) {
  if (p <= 0) return 0;
  if (p >= 1) return Infinity;
  let lo = 0, hi = df + 10 * Math.sqrt(2 * df);
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (chiSqCDF(mid, df) < p) lo = mid;
    else hi = mid;
    if (hi - lo < 1e-10) break;
  }
  return (lo + hi) / 2;
}
function normalCDF(z) {
  return 0.5 * (1 + erf(z / Math.SQRT2));
}
function normalSurvival(z) {
  if (z >= 0) {
    return 0.5 * erfc(z / Math.SQRT2);
  }
  return 1 - 0.5 * erfc(-z / Math.SQRT2);
}
function pKendallExact(q, n) {
  q = Math.floor(q + 1e-7);
  const maxT = n * (n - 1) / 2;
  if (q < 0) return 0;
  if (q > maxT) return 1;
  const cache = /* @__PURE__ */ new Map();
  function ckendall(k, nn) {
    const u = nn * (nn - 1) / 2;
    if (k < 0 || k > u) return 0;
    if (nn === 1) return k === 0 ? 1 : 0;
    if (!cache.has(nn)) {
      const arr = new Float64Array(u + 1);
      for (let i = 0; i <= u; i++) arr[i] = -1;
      cache.set(nn, arr);
    }
    const w = cache.get(nn);
    if (w[k] < 0) {
      let s = 0;
      for (let i = 0; i < nn; i++) s += ckendall(k - i, nn - 1);
      w[k] = s;
    }
    return w[k];
  }
  let p = 0;
  for (let j = 0; j <= q; j++) p += ckendall(j, n);
  let nfact = 1;
  for (let i = 2; i <= n; i++) nfact *= i;
  return p / nfact;
}
function pSpearmanExact(D, n, lowerTail) {
  const n3 = n * (n * n - 1) / 3;
  if (n <= 1) return lowerTail ? 0 : 1;
  if (D <= 0) return lowerTail ? 0 : 1;
  if (D > n3) return lowerTail ? 1 : 0;
  if (n <= 9) {
    const perm = new Array(n);
    for (let i = 0; i < n; i++) perm[i] = i + 1;
    let nfac = 1;
    for (let i = 1; i <= n; i++) nfac *= i;
    let ifr = 0;
    for (let m = 0; m < nfac; m++) {
      let ise = 0;
      for (let i = 0; i < n; i++) {
        const d = i + 1 - perm[i];
        ise += d * d;
      }
      if (D <= ise) ifr++;
      let n1 = n;
      do {
        const mt = perm[0];
        for (let i = 1; i < n1; i++) perm[i - 1] = perm[i];
        n1--;
        perm[n1] = mt;
        if (mt !== n1 + 1 || n1 <= 1) break;
      } while (true);
    }
    return (lowerTail ? nfac - ifr : ifr) / nfac;
  }
  const c1 = 0.2274, c2 = 0.2531, c3 = 0.1745, c4 = 0.0758;
  const c5 = 0.1033, c6 = 0.3932, c7 = 0.0879, c8 = 0.0151;
  const c9 = 72e-4, c10 = 0.0831, c11 = 0.0131, c12 = 46e-5;
  const yn = n;
  const b = 1 / yn;
  const x = (6 * (D - 1) * b / (yn * yn - 1) - 1) * Math.sqrt(yn - 1);
  let y = x * x;
  const u = x * b * (c1 + b * (c2 + c3 * b) + y * (-c4 + b * (c5 + c6 * b) - y * b * (c7 + c8 * b - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))));
  y = x * x;
  const correction = u / Math.exp(y / 2);
  let pv;
  if (lowerTail) {
    pv = -correction + normalCDF(x);
  } else {
    pv = correction + normalSurvival(x);
  }
  if (pv < 0) pv = 0;
  if (pv > 1) pv = 1;
  return pv;
}
function normalQuantile(p) {
  if (p <= 0 || p >= 1) throw new Error(`normalQuantile: p must be in (0,1), got ${p}`);
  const a1 = -39.69683028665376, a2 = 220.9460984245205;
  const a3 = -275.9285104469687, a4 = 138.357751867269;
  const a5 = -30.66479806614716, a6 = 2.506628277459239;
  const b1 = -54.47609879822406, b2 = 161.5858368580409;
  const b3 = -155.6989798598866, b4 = 66.80131188771972;
  const b5 = -13.28068155288572;
  const c1 = -0.007784894002430293, c2 = -0.3223964580411365;
  const c3 = -2.400758277161838, c4 = -2.549732539343734;
  const c5 = 4.374664141464968, c6 = 2.938163982698783;
  const d1 = 0.007784695709041462, d2 = 0.3224671290700398;
  const d3 = 2.445134137142996, d4 = 3.754408661907416;
  const pLow = 0.02425, pHigh = 1 - pLow;
  if (p < pLow) {
    const q = Math.sqrt(-2 * Math.log(p));
    return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  } else if (p <= pHigh) {
    const q = p - 0.5, r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  }
}
function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const result = 1 - poly * Math.exp(-x * x);
  return x >= 0 ? result : -result;
}
function erfc(x) {
  const ax = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * ax);
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const val = poly * Math.exp(-ax * ax);
  return x >= 0 ? val : 2 - val;
}
function nelderMead(fn, x0, opts = {}) {
  const n = x0.length;
  const maxIter = opts.maxIter ?? 1e3 * n;
  const tol = opts.tol ?? 1e-8;
  const alpha = opts.alpha ?? 1;
  const beta = opts.beta ?? 0.5;
  const gam = opts.gamma ?? 2;
  const delta = opts.delta ?? 0.5;
  const simplex = Array.from({ length: n + 1 }, (_, i) => {
    const v = [...x0];
    if (i > 0) v[i - 1] = (v[i - 1] ?? 0) + (Math.abs(v[i - 1] ?? 0) > 1e-8 ? 0.05 * (v[i - 1] ?? 0) : 25e-5);
    return v;
  });
  let fvals = simplex.map((v) => fn(v));
  let iter = 0;
  let converged = false;
  for (; iter < maxIter; iter++) {
    const order = fvals.map((_, i) => i).sort((a, b) => fvals[a] - fvals[b]);
    const sorted = order.map((i) => simplex[i]);
    fvals = order.map((i) => fvals[i]);
    for (let i = 0; i <= n; i++) simplex[i] = sorted[i];
    const fRange = (fvals[n] ?? 0) - (fvals[0] ?? 0);
    if (fRange < tol) {
      converged = true;
      break;
    }
    const centroid = Array.from(
      { length: n },
      (_, j) => simplex.slice(0, n).reduce((s, v) => s + (v[j] ?? 0), 0) / n
    );
    const worst = simplex[n];
    const reflected = centroid.map((c, j) => c + alpha * (c - (worst[j] ?? 0)));
    const fr = fn(reflected);
    if (fr < (fvals[0] ?? 0)) {
      const expanded = centroid.map((c, j) => c + gam * (reflected[j] - c));
      const fe = fn(expanded);
      simplex[n] = fe < fr ? expanded : reflected;
      fvals[n] = fe < fr ? fe : fr;
    } else if (fr < (fvals[n - 1] ?? 0)) {
      simplex[n] = reflected;
      fvals[n] = fr;
    } else {
      const contracted = centroid.map((c, j) => c + beta * ((worst[j] ?? 0) - c));
      const fc = fn(contracted);
      if (fc < (fvals[n] ?? 0)) {
        simplex[n] = contracted;
        fvals[n] = fc;
      } else {
        const best = simplex[0];
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, j) => (best[j] ?? 0) + delta * (v - (best[j] ?? 0)));
          fvals[i] = fn(simplex[i]);
        }
      }
    }
  }
  const order0 = fvals.map((_, i) => i).sort((a, b) => fvals[a] - fvals[b]);
  return {
    x: simplex[order0[0]],
    fval: fvals[order0[0]],
    iterations: iter,
    converged
  };
}
function adjustPValues(pValues, method) {
  const n = pValues.length;
  if (n === 0) return [];
  if (method === "none") return [...pValues];
  if (method === "bonferroni") {
    return pValues.map((p) => Math.min(1, p * n));
  }
  if (method === "holm") {
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => a.p - b.p);
    const adj = new Array(n);
    let running = 0;
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k];
      running = Math.max(running, p * (n - k));
      adj[i] = Math.min(1, running);
    }
    return adj;
  }
  if (method === "BH") {
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p);
    const adj = new Array(n);
    let minSoFar = Infinity;
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k];
      const rank2 = n - k;
      const adjusted = p * n / rank2;
      minSoFar = Math.min(minSoFar, adjusted);
      adj[i] = Math.min(1, minSoFar);
    }
    return adj;
  }
  if (method === "BY") {
    const cm = pValues.reduce((s, _, k) => s + 1 / (k + 1), 0);
    const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p);
    const adj = new Array(n);
    let minSoFar = Infinity;
    for (let k = 0; k < n; k++) {
      const { p, i } = order[k];
      const rank2 = n - k;
      const adjusted = p * n * cm / rank2;
      minSoFar = Math.min(minSoFar, adjusted);
      adj[i] = Math.min(1, minSoFar);
    }
    return adj;
  }
  throw new Error(`Unknown p-adjust method: ${method}`);
}
function mean(x) {
  if (x.length === 0) throw new Error("mean: empty array");
  return x.reduce((s, v) => s + Number(v), 0) / x.length;
}
function variance(x) {
  if (x.length < 2) throw new Error("variance: need at least 2 observations");
  const m = mean(x);
  if (!isFinite(m)) return NaN;
  return x.reduce((s, v) => s + (v - m) ** 2, 0) / (x.length - 1);
}
function sd(x) {
  return Math.sqrt(variance(x));
}
function se(x) {
  return sd(x) / Math.sqrt(x.length);
}
function median(x) {
  if (x.length === 0) throw new Error("median: empty array");
  const sorted = [...x].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
}
function quantile(x, p) {
  if (p < 0 || p > 1) throw new Error(`quantile: p must be in [0,1], got ${p}`);
  const sorted = [...x].sort((a, b) => a - b);
  const n = sorted.length;
  const h = (n - 1) * p;
  const lo = Math.floor(h);
  const hi = Math.ceil(h);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
}
function sortAsc(x) {
  return [...x].sort((a, b) => a - b);
}
function ss(x) {
  const m = mean(x);
  return x.reduce((s, v) => s + (v - m) ** 2, 0);
}
function rank(x) {
  const n = x.length;
  const order = x.map((v, i2) => ({ v, i: i2 })).sort((a, b) => a.v - b.v);
  const ranks = new Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    while (j < n && order[j].v === order[i].v) j++;
    const avgRank = (i + j - 1) / 2 + 1;
    for (let k = i; k < j; k++) ranks[order[k].i] = avgRank;
    i = j;
  }
  return ranks;
}
function cov(x, y) {
  if (x.length !== y.length) throw new Error("cov: arrays must have equal length");
  if (x.length < 2) throw new Error("cov: need at least 2 observations");
  const mx = mean(x), my = mean(y);
  return x.reduce((s, v, i) => s + (v - mx) * ((y[i] ?? 0) - my), 0) / (x.length - 1);
}
function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}
function roundTo(v, n) {
  const factor = 10 ** n;
  return Math.round(v * factor) / factor;
}
function normalPDF(z) {
  return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
}
function ptukeyApprox(q, k, df) {
  if (q <= 0) return 0;
  if (!isFinite(q) || isNaN(q)) return NaN;
  if (k < 2) return NaN;
  if (df > 5e3) return ptukeyInfDf(q, k);
  let wStep;
  if (df <= 100) wStep = 1;
  else if (df <= 800) wStep = 0.5;
  else if (df <= 5e3) wStep = 0.25;
  else wStep = 0.125;
  const halfDf = df / 2;
  const rate = df / 4;
  const logNorm = halfDf * Math.log(rate) - logGamma(halfDf);
  const glPts = 16;
  const gl = gaussLegendreRef(glPts);
  let result = 0;
  let vLo = 0;
  for (let interval = 0; interval < 200; interval++) {
    const vHi = vLo + wStep;
    const mid = (vLo + vHi) / 2;
    const half = (vHi - vLo) / 2;
    let subResult = 0;
    for (let i = 0; i < glPts; i++) {
      const v = mid + half * gl.nodes[i];
      const w = half * gl.weights[i];
      if (v <= 0) continue;
      const logDens = logNorm + (halfDf - 1) * Math.log(v) - rate * v;
      const dens = Math.exp(logDens);
      if (dens < 1e-300) continue;
      const scaledQ = q * Math.sqrt(v / 2);
      const pInf = ptukeyInfDf(scaledQ, k);
      subResult += w * dens * pInf;
    }
    result += subResult;
    vLo = vHi;
    if (vLo > 2 && Math.abs(subResult) < 1e-10) break;
  }
  return Math.max(0, Math.min(1, result));
}
function ptukeyInfDf(q, k) {
  if (q <= 0) return 0;
  const lo = -8;
  const hi = 8;
  const npts = 64;
  const { nodes, weights } = gaussLegendreAB(npts, lo, hi);
  let result = 0;
  for (let i = 0; i < npts; i++) {
    const z = nodes[i];
    const w = weights[i];
    const phi = normalPDF(z);
    const diff = normalCDF(z + q) - normalCDF(z);
    if (diff <= 0) continue;
    result += w * phi * Math.pow(diff, k - 1);
  }
  return Math.max(0, Math.min(1, k * result));
}
function pValueStudentizedRangeApprox(q, k, df) {
  return Math.max(0, Math.min(1, 1 - ptukeyApprox(q, k, df)));
}
function gaussLegendreAB(n, a, b) {
  const ref = gaussLegendreRef(n);
  const mid = (a + b) / 2;
  const half = (b - a) / 2;
  return {
    nodes: ref.nodes.map((t) => mid + half * t),
    weights: ref.weights.map((w) => half * w)
  };
}
function gaussLegendreRef(n) {
  const nodes = [];
  const weights = [];
  const m = Math.ceil(n / 2);
  for (let i = 0; i < m; i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));
    let pp = 0;
    for (let iter = 0; iter < 100; iter++) {
      let p0 = 1;
      let p1 = z;
      for (let j = 2; j <= n; j++) {
        const p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j;
        p0 = p1;
        p1 = p2;
      }
      pp = n * (z * p1 - p0) / (z * z - 1);
      const zOld = z;
      z = z - p1 / pp;
      if (Math.abs(z - zOld) < 1e-15) break;
    }
    nodes.push(z);
    nodes.push(-z);
    const w = 2 / ((1 - z * z) * pp * pp);
    weights.push(w);
    weights.push(w);
  }
  if (n % 2 === 1) {
    nodes.pop();
    weights.pop();
  }
  return { nodes, weights };
}

// src/core/apa.ts
function formatP(p) {
  if (p < 1e-3) return "p < .001";
  const rounded = roundTo(p, 3).toFixed(3).replace("0.", ".");
  return `p = ${rounded}`;
}
function formatStat(v, decimals = 2) {
  return roundTo(v, decimals).toFixed(decimals);
}
function formatCI(ci, decimals = 2) {
  return `[${roundTo(ci[0], decimals).toFixed(decimals)}, ${roundTo(ci[1], decimals).toFixed(decimals)}]`;
}
function formatDF(df) {
  if (typeof df === "number") {
    return Number.isInteger(df) ? String(df) : df.toFixed(2);
  }
  const [df1, df2] = df;
  const f1 = Number.isInteger(df1) ? String(df1) : df1.toFixed(2);
  const f2 = Number.isInteger(df2) ? String(df2) : df2.toFixed(2);
  return `${f1}, ${f2}`;
}
function formatTTest(t, df, pValue, effectSize, effectName, ci, ciLevel = 0.95) {
  const ciPct = Math.round(ciLevel * 100);
  return `t(${formatDF(df)}) = ${formatStat(t)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}, ${ciPct}% CI ${formatCI(ci)}`;
}
function formatANOVA(F, df1, df2, pValue, effectSize, effectName = "\u03B7\xB2") {
  return `F(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`;
}
function formatRMANOVA(F, df1, df2, pValue, effectSize, effectName = "\u03B7\xB2_p", correction) {
  const prefix = correction === "Greenhouse-Geisser" ? "F_GG" : correction === "Huynh-Feldt" ? "F_HF" : "F";
  return `${prefix}(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`;
}
function formatChiSq(chiSq, df, pValue, effectSize, effectName = "V") {
  return `\u03C7\xB2(${formatDF(df)}) = ${formatStat(chiSq)}, ${formatP(pValue)}, ${effectName} = ${formatStat(effectSize)}`;
}
function formatCorrelation(r, df, pValue, ci, method = "r", ciLevel = 0.95) {
  const ciPct = Math.round(ciLevel * 100);
  return `${method}(${formatDF(df)}) = ${formatStat(r)}, ${formatP(pValue)}, ${ciPct}% CI ${formatCI(ci)}`;
}
function formatRegression(r2, adjR2, F, df1, df2, pValue) {
  return `R\xB2 = ${formatStat(r2)}, adj. R\xB2 = ${formatStat(adjR2)}, F(${formatDF([df1, df2])}) = ${formatStat(F)}, ${formatP(pValue)}`;
}
function formatMannWhitney(W, pValue, r) {
  return `W = ${formatStat(W, 0)}, ${formatP(pValue)}, r = ${formatStat(r)}`;
}
function formatKruskalWallis(H, df, pValue, etaSq) {
  return `H(${formatDF(df)}) = ${formatStat(H)}, ${formatP(pValue)}, \u03B7\xB2_H = ${formatStat(etaSq)}`;
}
function formatLMM(icc, aic, bic, logLik) {
  return `ICC = ${formatStat(icc)}, AIC = ${formatStat(aic, 1)}, BIC = ${formatStat(bic, 1)}, logLik = ${formatStat(logLik, 1)}`;
}
function interpretCohensD(d) {
  const abs = Math.abs(d);
  if (abs < 0.2) return "negligible";
  if (abs < 0.5) return "small";
  if (abs < 0.8) return "medium";
  return "large";
}
function interpretEtaSq(eta) {
  if (eta < 0.01) return "negligible";
  if (eta < 0.06) return "small";
  if (eta < 0.14) return "medium";
  return "large";
}
function interpretR(r) {
  const abs = Math.abs(r);
  if (abs < 0.1) return "negligible";
  if (abs < 0.3) return "small";
  if (abs < 0.5) return "medium";
  if (abs < 0.7) return "large";
  return "very large";
}
function interpretCramerV(v, df) {
  const small = df === 1 ? 0.1 : df === 2 ? 0.07 : 0.06;
  const medium = df === 1 ? 0.3 : df === 2 ? 0.21 : 0.17;
  if (v < small) return "negligible";
  if (v < medium) return "small";
  if (v < medium * 1.5) return "medium";
  return "large";
}
function formatCFAFit(fit) {
  return `\u03C7\xB2(${formatDF(fit.df)}) = ${formatStat(fit.chiSq)}, ${formatP(fit.pValue)}, RMSEA = ${formatStat(fit.rmsea, 3)} ${formatCI(fit.rmseaCI, 3)}, CFI = ${formatStat(fit.cfi, 3)}, TLI = ${formatStat(fit.tli, 3)}, SRMR = ${formatStat(fit.srmr, 3)}`;
}
function interpretEffect(value, thresholds) {
  const abs = Math.abs(value);
  if (abs < thresholds[0]) return "negligible";
  if (abs < thresholds[1]) return "small";
  if (abs < thresholds[2]) return "medium";
  return "large";
}
//# sourceMappingURL=index.cjs.map
