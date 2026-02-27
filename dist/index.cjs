"use strict";
var __create = Object.create;
var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __getProtoOf = Object.getPrototypeOf;
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
var __toESM = (mod, isNodeMode, target) => (target = mod != null ? __create(__getProtoOf(mod)) : {}, __copyProps(
  // If the importer is in node compatibility mode or this is not an ESM
  // file that has been converted to a CommonJS file using a Babel-
  // compatible transform (i.e. "__esModule" has not been set), then set
  // "default" to the CommonJS "module.exports" for node compatibility.
  isNodeMode || !mod || !mod.__esModule ? __defProp(target, "default", { value: mod, enumerable: true }) : target,
  mod
));
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// src/index.ts
var index_exports = {};
__export(index_exports, {
  CARM_PALETTE: () => CARM_PALETTE,
  DARK_THEME: () => DARK_THEME,
  DEFAULT_THEME: () => DEFAULT_THEME,
  Matrix: () => Matrix,
  OKABE_ITO: () => OKABE_ITO,
  addCaption: () => addCaption,
  addNLabel: () => addNLabel,
  addRegressionEquation: () => addRegressionEquation,
  addStatBadge: () => addStatBadge,
  addSubtitle: () => addSubtitle,
  adjustPValues: () => adjustPValues,
  analyze: () => analyze,
  applyTheme: () => applyTheme,
  betaFn: () => betaFn,
  chiSqCDF: () => chiSqCDF,
  chiSqPValue: () => chiSqPValue,
  chiSqQuantile: () => chiSqQuantile,
  chiSquareTest: () => chiSquareTest,
  ciMean: () => ciMean,
  clamp: () => clamp,
  cohensD: () => cohensD,
  cohensDCI: () => cohensDCI,
  cohensDPaired: () => cohensDPaired,
  computeBLUPs: () => computeBLUPs,
  contingencyTable: () => contingencyTable,
  correlationMatrix: () => correlationMatrix,
  cov: () => cov,
  cutTree: () => cutTree,
  cutTreeHeight: () => cutTreeHeight,
  describe: () => describe,
  detectFieldType: () => detectFieldType,
  dunnTest: () => dunnTest,
  etaSquared: () => etaSquared,
  etaSquaredKW: () => etaSquaredKW,
  euclideanDistMatrix: () => euclideanDistMatrix,
  exportPNG: () => exportPNG,
  exportSVG: () => exportSVG,
  fDistCDF: () => fDistCDF,
  fDistPValue: () => fDistPValue,
  findBestGMM: () => findBestGMM,
  fisherExactTest: () => fisherExactTest,
  fitGMM: () => fitGMM,
  fitGMMRange: () => fitGMMRange,
  fitKMeansRange: () => fitKMeansRange,
  fitLCA: () => fitLCA,
  fitLTA: () => fitLTA,
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
  formatRegression: () => formatRegression,
  formatStat: () => formatStat,
  formatTTest: () => formatTTest,
  formatTooltipRow: () => formatTooltipRow,
  frequencyTable: () => frequencyTable,
  friedmanTest: () => friedmanTest,
  gamesHowell: () => gamesHowell,
  gamma: () => gamma,
  getColor: () => getColor,
  goodnessOfFit: () => goodnessOfFit,
  hedgesG: () => hedgesG,
  hideTooltip: () => hideTooltip,
  incompleteBeta: () => incompleteBeta,
  incompleteGamma: () => incompleteGamma,
  interpretCohensD: () => interpretCohensD,
  interpretCramerV: () => interpretCramerV,
  interpretEffect: () => interpretEffect,
  interpretEtaSq: () => interpretEtaSq,
  interpretR: () => interpretR,
  inverseTransform: () => inverseTransform,
  kDistancePlot: () => kDistancePlot,
  kendallTau: () => kendallTau,
  kruskalWallis: () => kruskalWallis,
  kurtosis: () => kurtosis,
  linearRegression: () => linearRegression,
  logBeta: () => logBeta,
  logGamma: () => logGamma,
  logisticRegression: () => logisticRegression,
  mannWhitneyU: () => mannWhitneyU,
  mean: () => mean,
  median: () => median,
  multipleRegression: () => multipleRegression,
  nelderMead: () => nelderMead,
  normalCDF: () => normalCDF,
  normalQuantile: () => normalQuantile,
  normalSurvival: () => normalSurvival,
  omegaSquared: () => omegaSquared,
  oneWayANOVA: () => oneWayANOVA,
  pKendallExact: () => pKendallExact,
  pSpearmanExact: () => pSpearmanExact,
  pValueStudentizedRangeApprox: () => pValueStudentizedRangeApprox,
  partialCorrelation: () => partialCorrelation,
  pearsonCorrelation: () => pearsonCorrelation,
  phiCoefficient: () => phiCoefficient,
  polynomialRegression: () => polynomialRegression,
  predictGMM: () => predictGMM,
  predictKMeans: () => predictKMeans,
  preprocessData: () => preprocessData,
  ptukeyApprox: () => ptukeyApprox,
  quantile: () => quantile,
  rank: () => rank,
  rankBiserial: () => rankBiserial,
  rankBiserialWilcoxon: () => rankBiserialWilcoxon,
  regressionDiagnostics: () => regressionDiagnostics,
  renderAlluvialPlot: () => renderAlluvialPlot,
  renderArcDiagram: () => renderArcDiagram,
  renderAreaChart: () => renderAreaChart,
  renderBarStats: () => renderBarStats,
  renderBoxplot: () => renderBoxplot,
  renderBrackets: () => renderBrackets,
  renderBubbleChart: () => renderBubbleChart,
  renderChordDiagram: () => renderChordDiagram,
  renderCoefPlot: () => renderCoefPlot,
  renderCorrelogram: () => renderCorrelogram,
  renderDendrogram: () => renderDendrogram,
  renderDensity: () => renderDensity,
  renderDistribution: () => renderDistribution,
  renderDotPlot: () => renderDotPlot,
  renderEdgeBundling: () => renderEdgeBundling,
  renderFAPlot: () => renderFAPlot,
  renderForestPlot: () => renderForestPlot,
  renderFunnel: () => renderFunnel,
  renderGridLines: () => renderGridLines,
  renderGroupedBar: () => renderGroupedBar,
  renderHistogram: () => renderHistogram,
  renderLineChart: () => renderLineChart,
  renderLollipop: () => renderLollipop,
  renderMarimekko: () => renderMarimekko,
  renderMixedPlot: () => renderMixedPlot,
  renderMosaicPlot: () => renderMosaicPlot,
  renderPCAPlot: () => renderPCAPlot,
  renderPairPlot: () => renderPairPlot,
  renderParallelCoords: () => renderParallelCoords,
  renderPareto: () => renderPareto,
  renderPieChart: () => renderPieChart,
  renderQQPlot: () => renderQQPlot,
  renderROCCurve: () => renderROCCurve,
  renderRadarChart: () => renderRadarChart,
  renderRaincloud: () => renderRaincloud,
  renderResidualPanel: () => renderResidualPanel,
  renderScatterStats: () => renderScatterStats,
  renderSparkline: () => renderSparkline,
  renderStripPlot: () => renderStripPlot,
  renderSunburst: () => renderSunburst,
  renderSwarmPlot: () => renderSwarmPlot,
  renderTreemap: () => renderTreemap,
  renderViolinBox: () => renderViolinBox,
  renderWaffleChart: () => renderWaffleChart,
  renderXAxis: () => renderXAxis,
  renderYAxis: () => renderYAxis,
  roundTo: () => roundTo,
  runCFA: () => runCFA,
  runDBSCAN: () => runDBSCAN,
  runEFA: () => runEFA,
  runFADiagnostics: () => runFADiagnostics,
  runHierarchical: () => runHierarchical,
  runKMeans: () => runKMeans,
  runLMM: () => runLMM,
  runPCA: () => runPCA,
  screeData: () => screeData,
  sd: () => sd,
  se: () => se,
  shapiroWilk: () => shapiroWilk,
  showTooltip: () => showTooltip,
  silhouetteScores: () => silhouetteScores,
  skewness: () => skewness,
  solveLinear: () => solveLinear,
  sortAsc: () => sortAsc,
  spearmanCorrelation: () => spearmanCorrelation,
  ss: () => ss,
  tDistCDF: () => tDistCDF,
  tDistPValue: () => tDistPValue,
  tDistQuantile: () => tDistQuantile,
  tTestIndependent: () => tTestIndependent,
  tTestPaired: () => tTestPaired,
  themeColorScale: () => themeColorScale,
  totalBracketHeight: () => totalBracketHeight,
  trimmedMean: () => trimmedMean,
  tukeyHSD: () => tukeyHSD,
  variance: () => variance,
  varimaxRotation: () => varimaxRotation,
  wilcoxonSignedRank: () => wilcoxonSignedRank
});
module.exports = __toCommonJS(index_exports);

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

// src/stats/descriptive.ts
function mode(x) {
  const counts = /* @__PURE__ */ new Map();
  for (const v of x) counts.set(v, (counts.get(v) ?? 0) + 1);
  const maxCount = Math.max(...counts.values());
  const modes = [];
  for (const [v, c] of counts) {
    if (c === maxCount) modes.push(v);
  }
  return modes.sort((a, b) => a - b);
}
function trimmedMean(x, alpha = 0.05) {
  const sorted = sortAsc(x);
  const n = sorted.length;
  const trim = Math.floor(n * alpha);
  const trimmed = sorted.slice(trim, n - trim);
  if (trimmed.length === 0) throw new Error("trimmedMean: too much trimming, no data remains");
  return mean(trimmed);
}
function skewness(x) {
  const n = x.length;
  if (n < 3) throw new Error("skewness: need at least 3 observations");
  const m = mean(x);
  const s = sd(x);
  if (s === 0) return 0;
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 3, 0);
  return n / ((n - 1) * (n - 2)) * sum;
}
function kurtosis(x) {
  const n = x.length;
  if (n < 4) throw new Error("kurtosis: need at least 4 observations");
  const m = mean(x);
  const s = sd(x);
  if (s === 0) return 0;
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 4, 0);
  const k1 = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * sum;
  const k2 = 3 * (n - 1) ** 2 / ((n - 2) * (n - 3));
  return k1 - k2;
}
function ciMean(x, ciLevel = 0.95) {
  const n = x.length;
  if (n < 2) throw new Error("ciMean: need at least 2 observations");
  const m = mean(x);
  const s = se(x);
  const t = tDistQuantile(1 - (1 - ciLevel) / 2, n - 1);
  return [m - t * s, m + t * s];
}
var SW_C1 = [0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056];
var SW_C2 = [0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633];
function swPoly(c, x) {
  let r = 0;
  for (let i = c.length - 1; i >= 0; i--) r = r * x + (c[i] ?? 0);
  return r;
}
function shapiroWilk(x) {
  const n = x.length;
  if (n < 3) throw new Error("shapiroWilk: need at least 3 observations");
  if (n > 5e3) throw new Error("shapiroWilk: n > 5000 not supported");
  const sorted = sortAsc(x);
  const nn2 = Math.floor(n / 2);
  const a = new Array(nn2).fill(0);
  if (n === 3) {
    a[0] = Math.SQRT1_2;
  } else if (n === 4) {
    a[0] = 0.6872;
    a[1] = 0.1677;
  } else if (n === 5) {
    a[0] = 0.6646;
    a[1] = 0.2413;
  } else {
    for (let i = 0; i < nn2; i++) {
      a[i] = normalQuantile((i + 1 - 0.375) / (n + 0.25));
    }
    let summ2 = 0;
    for (let i = 0; i < nn2; i++) summ2 += a[i] * a[i];
    summ2 *= 2;
    const ssumm2 = Math.sqrt(summ2);
    const rsn = 1 / Math.sqrt(n);
    const a0orig = a[0];
    const a1orig = a[1];
    const a1corr = swPoly(SW_C1, rsn) - a0orig / ssumm2;
    if (n > 5) {
      const a2corr = -a1orig / ssumm2 + swPoly(SW_C2, rsn);
      const num = summ2 - 2 * a0orig * a0orig - 2 * a1orig * a1orig;
      const den = 1 - 2 * a1corr * a1corr - 2 * a2corr * a2corr;
      const fac = num > 0 && den > 0 ? Math.sqrt(num / den) : 1;
      a[0] = a1corr;
      a[1] = a2corr;
      for (let i = 2; i < nn2; i++) a[i] = (a[i] ?? 0) / -fac;
    } else {
      const num = summ2 - 2 * a0orig * a0orig;
      const den = 1 - 2 * a1corr * a1corr;
      const fac = num > 0 && den > 0 ? Math.sqrt(num / den) : 1;
      a[0] = a1corr;
      for (let i = 1; i < nn2; i++) a[i] = (a[i] ?? 0) / -fac;
    }
  }
  let w1 = 0;
  for (let i = 0; i < nn2; i++) {
    w1 += (a[i] ?? 0) * ((sorted[n - 1 - i] ?? 0) - (sorted[i] ?? 0));
  }
  const sst = variance(sorted) * (n - 1);
  const W = sst > 0 ? Math.min(1, w1 * w1 / sst) : 1;
  const pValue = shapiroWilkPValue(W, n);
  return { statistic: roundTo(W, 4), pValue: roundTo(pValue, 4) };
}
function polyAsc(c, x) {
  let r = c[c.length - 1];
  for (let i = c.length - 2; i >= 0; i--) r = r * x + c[i];
  return r;
}
var SW_G = [-2.273, 0.459];
var SW_C3 = [0.544, -0.39978, 0.025054, -6714e-7];
var SW_C4 = [1.3822, -0.77857, 0.062767, -20322e-7];
var SW_C5 = [-1.5861, -0.31082, -0.083751, 38915e-7];
var SW_C6 = [-0.4803, -0.082676, 30302e-7];
function shapiroWilkPValue(W, n) {
  if (n === 3) {
    const pi6 = 6 / Math.PI;
    const stqr = normalQuantile(0.75);
    const pw = pi6 * (Math.asin(Math.sqrt(W)) - stqr);
    return Math.max(0, Math.min(1, pw));
  }
  const y0 = Math.log(1 - W);
  let z;
  if (n <= 11) {
    const gamma2 = polyAsc(SW_G, n);
    if (y0 >= gamma2) return 0;
    const y = -Math.log(gamma2 - y0);
    const mu = polyAsc(SW_C3, n);
    const sigma = Math.exp(polyAsc(SW_C4, n));
    z = (y - mu) / sigma;
  } else {
    const xx = Math.log(n);
    const mu = polyAsc(SW_C5, xx);
    const sigma = Math.exp(polyAsc(SW_C6, xx));
    z = (y0 - mu) / sigma;
  }
  return Math.max(0, Math.min(1, 1 - normalCDF(z)));
}
function describe(x, ciLevel = 0.95) {
  if (x.length === 0) throw new Error("describe: empty array");
  const n = x.length;
  const m = mean(x);
  const med = median(x);
  const modes = mode(x);
  const tm = trimmedMean(x, 0.05);
  const s = sd(x);
  const sem = se(x);
  const v = variance(x);
  const sorted = sortAsc(x);
  const mn = sorted[0];
  const mx = sorted[n - 1];
  const q1 = quantile(x, 0.25);
  const q3 = quantile(x, 0.75);
  const iqr = q3 - q1;
  const skew = n >= 3 ? skewness(x) : 0;
  const kurt = n >= 4 ? kurtosis(x) : 0;
  const ci = ciMean(x, ciLevel);
  const sw = n >= 3 ? shapiroWilk(x) : { statistic: NaN, pValue: NaN };
  const ciPct = Math.round(ciLevel * 100);
  const formatted = [
    `n = ${n}`,
    `M = ${roundTo(m, 2)}, SD = ${roundTo(s, 2)}, SE = ${roundTo(sem, 2)}`,
    `Mdn = ${roundTo(med, 2)}, IQR = ${roundTo(iqr, 2)}`,
    `Skew = ${roundTo(skew, 2)}, Kurt = ${roundTo(kurt, 2)}`,
    `${ciPct}% CI [${roundTo(ci[0], 2)}, ${roundTo(ci[1], 2)}]`,
    n >= 3 ? `Shapiro-Wilk W = ${roundTo(sw.statistic, 3)}, ${formatP(sw.pValue)}` : ""
  ].filter(Boolean).join("; ");
  return {
    n,
    mean: roundTo(m, 6),
    median: roundTo(med, 6),
    mode: modes,
    trimmedMean: roundTo(tm, 6),
    sd: roundTo(s, 6),
    se: roundTo(sem, 6),
    variance: roundTo(v, 6),
    min: mn,
    max: mx,
    range: mx - mn,
    iqr: roundTo(iqr, 6),
    q1: roundTo(q1, 6),
    q3: roundTo(q3, 6),
    skewness: roundTo(skew, 6),
    kurtosis: roundTo(kurt, 6),
    ci,
    ciLevel,
    shapiroWilk: sw,
    formatted
  };
}

// src/stats/effect-size.ts
function cohensD(x1, x2) {
  if (x1.length < 2 || x2.length < 2) throw new Error("cohensD: need at least 2 observations per group");
  const n1 = x1.length, n2 = x2.length;
  const m1 = mean(x1), m2 = mean(x2);
  const v1 = variance(x1), v2 = variance(x2);
  const sdPooled = Math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2));
  if (sdPooled === 0) return { value: 0, name: "Cohen's d", interpretation: "negligible" };
  const d = (m1 - m2) / sdPooled;
  return {
    value: d,
    name: "Cohen's d",
    interpretation: interpretCohensD(d)
  };
}
function cohensDPaired(diffs) {
  if (diffs.length < 2) throw new Error("cohensDPaired: need at least 2 differences");
  const m = mean(diffs);
  const s = sd(diffs);
  const d = s === 0 ? 0 : m / s;
  return {
    value: d,
    name: "Cohen's d",
    interpretation: interpretCohensD(d)
  };
}
function hedgesG(x1, x2) {
  const d = cohensD(x1, x2);
  const df = x1.length + x2.length - 2;
  const J = 1 - 3 / (4 * df - 1);
  const g = d.value * J;
  return {
    value: g,
    name: "Hedges' g",
    interpretation: interpretCohensD(g)
  };
}
function etaSquared(ssBetween, ssTotal) {
  if (ssTotal <= 0) return { value: 0, name: "\u03B7\xB2", interpretation: "negligible" };
  const eta2 = Math.max(0, Math.min(1, ssBetween / ssTotal));
  return {
    value: eta2,
    name: "\u03B7\xB2",
    interpretation: interpretEtaSq(eta2)
  };
}
function omegaSquared(ssBetween, ssTotal, dfBetween, msWithin) {
  const denom = ssTotal + msWithin;
  if (denom <= 0) return { value: 0, name: "\u03C9\xB2", interpretation: "negligible" };
  const omega2 = Math.max(0, (ssBetween - dfBetween * msWithin) / denom);
  return {
    value: omega2,
    name: "\u03C9\xB2",
    interpretation: interpretEtaSq(omega2)
  };
}
function rankBiserial(U, n1, n2) {
  const r = 1 - 2 * U / (n1 * n2);
  return {
    value: r,
    name: "r (rank-biserial)",
    interpretation: interpretR(r)
  };
}
function rankBiserialWilcoxon(T, n) {
  const maxT = n * (n + 1) / 2;
  const r = maxT > 0 ? T / maxT * 2 - 1 : 0;
  return {
    value: r,
    name: "r (rank-biserial)",
    interpretation: interpretR(r)
  };
}
function etaSquaredKW(H, k, n) {
  const eta2 = Math.max(0, (H - k + 1) / (n - k));
  return {
    value: eta2,
    name: "\u03B7\xB2_H",
    interpretation: interpretEtaSq(eta2)
  };
}
function cohensDCI(d, n1, n2, ciLevel = 0.95) {
  const sampleSE = Math.sqrt((n1 + n2) / (n1 * n2) + d * d / (2 * (n1 + n2)));
  const z = normalQuantileInline(1 - (1 - ciLevel) / 2);
  return [d - z * sampleSE, d + z * sampleSE];
}
function normalQuantileInline(p) {
  const a = [2.515517, 0.802853, 0.010328];
  const b = [1.432788, 0.189269, 1308e-6];
  const t = Math.sqrt(-2 * Math.log(p <= 0.5 ? p : 1 - p));
  const num = a[0] + a[1] * t + a[2] * t * t;
  const den = 1 + b[0] * t + b[1] * t * t + b[2] * t * t * t;
  const x = t - num / den;
  return p <= 0.5 ? -x : x;
}

// src/stats/frequency.ts
function frequencyTable(data) {
  const counts = /* @__PURE__ */ new Map();
  for (const v of data) counts.set(v, (counts.get(v) ?? 0) + 1);
  const total = data.length;
  const sorted = [...counts.entries()].sort((a, b) => {
    const av = typeof a[0] === "number" ? a[0] : String(a[0]);
    const bv = typeof b[0] === "number" ? b[0] : String(b[0]);
    return av < bv ? -1 : av > bv ? 1 : 0;
  });
  let cumulative = 0;
  return sorted.map(([value, count]) => {
    const relative = count / total;
    cumulative += relative;
    return { value, count, relative, cumulative: Math.min(1, cumulative) };
  });
}
function contingencyTable(group1, group2) {
  if (group1.length !== group2.length) throw new Error("contingencyTable: arrays must have equal length");
  const rowSet = [...new Set(group1)].sort();
  const colSet = [...new Set(group2)].sort();
  const table = Array.from(
    { length: rowSet.length },
    () => new Array(colSet.length).fill(0)
  );
  for (let i = 0; i < group1.length; i++) {
    const r = rowSet.indexOf(group1[i]);
    const c = colSet.indexOf(group2[i]);
    table[r][c]++;
  }
  return { table, rowLabels: rowSet, colLabels: colSet };
}
function expectedCounts(observed) {
  const R = observed.length;
  const C = observed[0].length;
  const rowTotals = observed.map((row) => row.reduce((s, v) => s + v, 0));
  const colTotals = Array.from(
    { length: C },
    (_, j) => observed.reduce((s, row) => s + (row[j] ?? 0), 0)
  );
  const total = rowTotals.reduce((s, v) => s + v, 0);
  return Array.from(
    { length: R },
    (_, i) => Array.from({ length: C }, (_2, j) => rowTotals[i] * (colTotals[j] ?? 0) / total)
  );
}
function chiSquareTest(observed, yatesCorrection = false) {
  const R = observed.length;
  if (R < 2) throw new Error("chiSquareTest: need at least 2 rows");
  const C = observed[0].length;
  if (C < 2) throw new Error("chiSquareTest: need at least 2 columns");
  const expected = expectedCounts(observed);
  const n = observed.reduce((s, row) => s + row.reduce((a, v) => a + v, 0), 0);
  let chiSq = 0;
  for (let i = 0; i < R; i++) {
    for (let j = 0; j < C; j++) {
      const o = observed[i][j] ?? 0;
      const e = expected[i][j] ?? 0;
      if (e < 5) {
      }
      const num = yatesCorrection ? Math.max(0, Math.abs(o - e) - 0.5) : o - e;
      chiSq += num * num / e;
    }
  }
  const df = (R - 1) * (C - 1);
  const pValue = chiSqPValue(chiSq, df);
  const minDim = Math.min(R, C) - 1;
  const cramersV = Math.sqrt(chiSq / (n * Math.max(1, minDim)));
  const effectSize = {
    value: cramersV,
    name: "Cram\xE9r's V",
    interpretation: interpretCramerV(cramersV, df)
  };
  const ci = [NaN, NaN];
  const formatted = formatChiSq(chiSq, df, pValue, cramersV, "V");
  return {
    testName: "Pearson's \u03C7\xB2",
    statistic: chiSq,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel: 0.95,
    n,
    formatted,
    table: frequencyTable(observed.flatMap(
      (row, i) => row.flatMap((count, j) => new Array(count).fill(`${i},${j}`))
    )),
    expectedCounts: expected
  };
}
function fisherExactTest(a, b, c, d) {
  if (a < 0 || b < 0 || c < 0 || d < 0) throw new Error("fisherExactTest: all cells must be non-negative");
  const n = a + b + c + d;
  const r1 = a + b;
  const c1 = a + c;
  const c2 = b + d;
  const pObserved = hypergeomPMF(a, r1, c1, n);
  let pValue = 0;
  const kMin = Math.max(0, r1 - c2);
  const kMax = Math.min(r1, c1);
  for (let k = kMin; k <= kMax; k++) {
    const p = hypergeomPMF(k, r1, c1, n);
    if (p <= pObserved + 1e-10) pValue += p;
  }
  pValue = Math.min(1, pValue);
  const oddsRatio = b === 0 || c === 0 ? (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)) : a * d / (b * c);
  const se2 = Math.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5));
  const z = normalQuantile(0.975);
  const logOR = Math.log(oddsRatio);
  const ci = [Math.exp(logOR - z * se2), Math.exp(logOR + z * se2)];
  const effectSize = {
    value: oddsRatio,
    name: "Odds ratio",
    interpretation: oddsRatio >= 3 ? "large" : oddsRatio >= 1.5 ? "medium" : "small"
  };
  return {
    testName: "Fisher's Exact Test",
    statistic: oddsRatio,
    df: 1,
    pValue,
    effectSize,
    ci,
    ciLevel: 0.95,
    n,
    formatted: `OR = ${oddsRatio.toFixed(2)}, ${formatP(pValue)}, 95% CI [${ci[0].toFixed(2)}, ${ci[1].toFixed(2)}]`
  };
}
function hypergeomPMF(k, n, K, N) {
  return Math.exp(
    logCombination(K, k) + logCombination(N - K, n - k) - logCombination(N, n)
  );
}
function logCombination(n, k) {
  if (k < 0 || k > n) return -Infinity;
  if (k === 0 || k === n) return 0;
  return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
}
function logFactorial(n) {
  if (n <= 1) return 0;
  let result = 0;
  for (let i = 2; i <= n; i++) result += Math.log(i);
  return result;
}
function phiCoefficient(a, b, c, d) {
  const denom = Math.sqrt((a + b) * (c + d) * (a + c) * (b + d));
  return denom === 0 ? 0 : (a * d - b * c) / denom;
}
function goodnessOfFit(observed, expected) {
  const n = observed.reduce((s, v) => s + v, 0);
  const k = observed.length;
  if (k < 2) throw new Error("goodnessOfFit: need at least 2 categories");
  const expCounts = expected ? expected.map((p) => p * n) : new Array(k).fill(n / k);
  let chiSq = 0;
  for (let i = 0; i < k; i++) {
    const o = observed[i] ?? 0;
    const e = expCounts[i] ?? 0;
    if (e <= 0) throw new Error(`goodnessOfFit: expected count must be > 0 at index ${i}`);
    chiSq += (o - e) ** 2 / e;
  }
  const df = k - 1;
  const pValue = chiSqPValue(chiSq, df);
  const w = Math.sqrt(chiSq / n);
  return {
    testName: "Chi-square goodness-of-fit",
    statistic: chiSq,
    df,
    pValue,
    effectSize: {
      value: w,
      name: "Cohen's w",
      interpretation: w < 0.1 ? "negligible" : w < 0.3 ? "small" : w < 0.5 ? "medium" : "large"
    },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: formatChiSq(chiSq, df, pValue, w, "w")
  };
}

// src/stats/comparison.ts
function tTestIndependent(x1, x2, equalVariances = false, ciLevel = 0.95, alternative = "two.sided") {
  if (x1.length < 2 || x2.length < 2) throw new Error("tTestIndependent: need at least 2 per group");
  const n1 = x1.length, n2 = x2.length;
  const m1 = mean(x1), m2 = mean(x2);
  const v1 = variance(x1), v2 = variance(x2);
  let df;
  let se2;
  if (equalVariances) {
    const spSq = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2);
    se2 = Math.sqrt(spSq * (1 / n1 + 1 / n2));
    df = n1 + n2 - 2;
  } else {
    se2 = Math.sqrt(v1 / n1 + v2 / n2);
    const num = (v1 / n1 + v2 / n2) ** 2;
    const den = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1);
    df = den > 0 ? num / den : n1 + n2 - 2;
  }
  const t = se2 === 0 ? 0 : (m1 - m2) / se2;
  const pFull = tDistPValue(t, df);
  const pValue = alternative === "two.sided" ? pFull : alternative === "less" ? t < 0 ? pFull / 2 : 1 - pFull / 2 : t > 0 ? pFull / 2 : 1 - pFull / 2;
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const diff = m1 - m2;
  const ci = [diff - tCrit * se2, diff + tCrit * se2];
  const effectSize = cohensD(x1, x2);
  const formatted = formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel);
  return {
    testName: equalVariances ? "Student's t-test" : "Welch's t-test",
    statistic: t,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n: n1 + n2,
    formatted
  };
}
function tTestPaired(x1, x2, ciLevel = 0.95) {
  if (x1.length !== x2.length) throw new Error("tTestPaired: arrays must have equal length");
  if (x1.length < 2) throw new Error("tTestPaired: need at least 2 pairs");
  const diffs = x1.map((v, i) => v - (x2[i] ?? 0));
  const n = diffs.length;
  const mDiff = mean(diffs);
  const seDiff = se(diffs);
  const df = n - 1;
  const t = seDiff === 0 ? 0 : mDiff / seDiff;
  const pValue = tDistPValue(t, df);
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const ci = [mDiff - tCrit * seDiff, mDiff + tCrit * seDiff];
  const effectSize = cohensDPaired(diffs);
  const formatted = formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel);
  return {
    testName: "Paired t-test",
    statistic: t,
    df,
    pValue,
    effectSize,
    ci,
    ciLevel,
    n,
    formatted
  };
}
function oneWayANOVA(groups) {
  if (groups.length < 2) throw new Error("oneWayANOVA: need at least 2 groups");
  const k = groups.length;
  const allValues = groups.flatMap((g) => [...g.values]);
  const n = allValues.length;
  const grandMean = mean(allValues);
  let ssBetween = 0, ssWithin = 0;
  const groupStats = groups.map((g) => {
    const gm = mean(g.values);
    const gn = g.values.length;
    ssBetween += gn * (gm - grandMean) ** 2;
    ssWithin += g.values.reduce((s, v) => s + (v - gm) ** 2, 0);
    return { label: g.label, n: gn, mean: gm, sd: sd(g.values) };
  });
  const ssTotal = ssBetween + ssWithin;
  const dfBetween = k - 1;
  const dfWithin = n - k;
  const msBetween = ssBetween / dfBetween;
  const msWithin = ssWithin / dfWithin;
  const F = msWithin === 0 ? Infinity : msBetween / msWithin;
  const pValue = fDistPValue(F, dfBetween, dfWithin);
  const omega = omegaSquared(ssBetween, ssTotal, dfBetween, msWithin);
  const formatted = formatANOVA(F, dfBetween, dfWithin, pValue, omega.value, "\u03C9\xB2");
  return {
    testName: "One-way ANOVA",
    statistic: F,
    df: [dfBetween, dfWithin],
    pValue,
    effectSize: omega,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted,
    groups: groupStats,
    ssBetween,
    ssWithin,
    ssTotal,
    msBetween,
    msWithin,
    dfBetween,
    dfWithin
  };
}
function mannWhitneyU(x1, x2, alternative = "two.sided") {
  if (x1.length < 1 || x2.length < 1) throw new Error("mannWhitneyU: empty group");
  const n1 = x1.length, n2 = x2.length;
  const combined = [
    ...x1.map((v, i) => ({ v, group: 1, i })),
    ...x2.map((v, i) => ({ v, group: 2, i }))
  ].sort((a, b) => a.v - b.v);
  const allCombined = rank(combined.map((d) => d.v));
  let R1 = 0;
  for (let i = 0; i < combined.length; i++) {
    if (combined[i].group === 1) R1 += allCombined[i];
  }
  const U1 = R1 - n1 * (n1 + 1) / 2;
  const N = n1 + n2;
  const muU = n1 * n2 / 2;
  const allVals = combined.map((d) => d.v);
  const tieCountsU = /* @__PURE__ */ new Map();
  for (const v of allVals) tieCountsU.set(v, (tieCountsU.get(v) ?? 0) + 1);
  let tieCorr = 0;
  for (const t of tieCountsU.values()) tieCorr += t * t * t - t;
  const varU = n1 * n2 / 12 * (N + 1 - tieCorr / (N * (N - 1)));
  const z = varU === 0 ? 0 : (U1 - muU) / Math.sqrt(varU);
  const zAbs = Math.abs(z);
  const pNormal = 2 * (1 - normalCDFInline(zAbs));
  const pValue = alternative === "two.sided" ? pNormal : alternative === "less" ? z < 0 ? pNormal / 2 : 1 - pNormal / 2 : z > 0 ? pNormal / 2 : 1 - pNormal / 2;
  const effect = rankBiserial(U1, n1, n2);
  const formatted = formatMannWhitney(U1, pValue, effect.value);
  return {
    testName: "Mann-Whitney U",
    statistic: U1,
    df: 0,
    pValue: Math.min(1, Math.max(0, pValue)),
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n: n1 + n2,
    formatted
  };
}
function normalCDFInline(z) {
  const x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const erf2 = 1 - poly * Math.exp(-x * x);
  return 0.5 * (1 + (z >= 0 ? erf2 : -erf2));
}
function wilcoxonExactP(W, n) {
  const maxW = n * (n + 1) / 2;
  let dist = new Array(maxW + 1).fill(0);
  dist[0] = 1;
  for (let k = 1; k <= n; k++) {
    const newDist = [...dist];
    for (let w = k; w <= maxW; w++) {
      newDist[w] += dist[w - k];
    }
    dist = newDist;
  }
  const total = Math.pow(2, n);
  let cumP = 0;
  for (let w = 0; w <= W && w <= maxW; w++) {
    cumP += (dist[w] ?? 0) / total;
  }
  return Math.min(1, 2 * cumP);
}
function wilcoxonSignedRank(x1, x2) {
  if (x1.length !== x2.length) throw new Error("wilcoxonSignedRank: arrays must match length");
  const diffs = x1.map((v, i) => v - (x2[i] ?? 0)).filter((d) => d !== 0);
  const n = diffs.length;
  if (n === 0) throw new Error("wilcoxonSignedRank: no non-zero differences");
  const absDiffs = diffs.map(Math.abs);
  const ranks_ = rank(absDiffs);
  let Wplus = 0, Wminus = 0;
  for (let i = 0; i < n; i++) {
    if ((diffs[i] ?? 0) > 0) Wplus += ranks_[i];
    else Wminus += ranks_[i];
  }
  const W = Math.min(Wplus, Wminus);
  let pValue;
  if (n <= 20) {
    pValue = wilcoxonExactP(W, n);
  } else {
    const muW = n * (n + 1) / 4;
    const varW = n * (n + 1) * (2 * n + 1) / 24;
    const z = varW === 0 ? 0 : (W + 0.5 - muW) / Math.sqrt(varW);
    pValue = 2 * (1 - normalCDFInline(Math.abs(z)));
  }
  const effect = { value: Wplus / (n * (n + 1) / 2) * 2 - 1, name: "r (rank-biserial)", interpretation: "small" };
  return {
    testName: "Wilcoxon Signed-Rank",
    statistic: W,
    df: 0,
    pValue: Math.min(1, Math.max(0, pValue)),
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `V = ${W}, ${formatP(pValue)}, r = ${effect.value.toFixed(2)}`
  };
}
function kruskalWallis(groups) {
  if (groups.length < 2) throw new Error("kruskalWallis: need at least 2 groups");
  const k = groups.length;
  const allValues = groups.flatMap((g) => [...g.values]);
  const n = allValues.length;
  const allRanks_ = rank(allValues);
  let offset = 0;
  let H = 0;
  for (const g of groups) {
    const gn = g.values.length;
    let Rj = 0;
    for (let i = 0; i < gn; i++) Rj += allRanks_[offset + i];
    H += Rj * Rj / gn;
    offset += gn;
  }
  H = 12 / (n * (n + 1)) * H - 3 * (n + 1);
  const tieCounts = /* @__PURE__ */ new Map();
  for (const v of allValues) tieCounts.set(v, (tieCounts.get(v) ?? 0) + 1);
  let C = 0;
  for (const t of tieCounts.values()) C += t * t * t - t;
  const correction = 1 - C / (n * n * n - n);
  if (correction > 0) H /= correction;
  const df = k - 1;
  const pValue = chiSqPValue(H, df);
  const effect = etaSquaredKW(H, k, n);
  const formatted = formatKruskalWallis(H, df, pValue, effect.value);
  return {
    testName: "Kruskal-Wallis",
    statistic: H,
    df,
    pValue,
    effectSize: effect,
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted
  };
}
function friedmanTest(data) {
  const n = data.length;
  if (n < 2) throw new Error("friedmanTest: need at least 2 subjects");
  const k = data[0].length;
  if (k < 2) throw new Error("friedmanTest: need at least 2 conditions");
  let sumRankSq = 0;
  for (const row of data) {
    const rowRanks = rank(row);
    for (const r of rowRanks) sumRankSq += r * r;
  }
  const colRankSums = Array.from(
    { length: k },
    (_, j) => data.reduce((s, row) => {
      const rowRanks = rank(row);
      return s + (rowRanks[j] ?? 0);
    }, 0)
  );
  const Rj2sum = colRankSums.reduce((s, r) => s + r * r, 0);
  const chi2 = 12 / (n * k * (k + 1)) * Rj2sum - 3 * n * (k + 1);
  const df = k - 1;
  const pValue = chiSqPValue(chi2, df);
  const w = chi2 / (n * (k - 1));
  return {
    testName: "Friedman Test",
    statistic: chi2,
    df,
    pValue,
    effectSize: {
      value: w,
      name: "Kendall's W",
      interpretation: w < 0.1 ? "negligible" : w < 0.3 ? "small" : w < 0.5 ? "medium" : "large"
    },
    ci: [NaN, NaN],
    ciLevel: 0.95,
    n,
    formatted: `\u03C7\xB2_F(${df}) = ${chi2.toFixed(2)}, ${formatP(pValue)}, W = ${w.toFixed(2)}`
  };
}

// src/stats/post-hoc.ts
function tukeyHSD(groups, msWithin, dfWithin, ciLevel = 0.95) {
  const results = [];
  const alpha = 1 - ciLevel;
  const k = groups.length;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const g1 = groups[i];
      const g2 = groups[j];
      const n1 = g1.values.length;
      const n2 = g2.values.length;
      const m1 = mean(g1.values);
      const m2 = mean(g2.values);
      const diff = m1 - m2;
      const se2 = Math.sqrt(msWithin / 2 * (1 / n1 + 1 / n2));
      const q = se2 === 0 ? 0 : Math.abs(diff) / se2;
      const pValue = pValueStudentizedRangeApprox(q, k, dfWithin);
      const tCrit = tDistQuantile(1 - alpha / (k * (k - 1)), dfWithin);
      const ciHalf = tCrit * se2;
      const ci = [diff - ciHalf, diff + ciHalf];
      results.push({
        group1: g1.label,
        group2: g2.label,
        meanDiff: diff,
        se: se2,
        statistic: q,
        pValue,
        pValueAdj: pValue,
        // already family-wise corrected
        ci,
        significant: pValue < alpha
      });
    }
  }
  return results;
}
function gamesHowell(groups, ciLevel = 0.95) {
  const results = [];
  const alpha = 1 - ciLevel;
  const k = groups.length;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const g1 = groups[i];
      const g2 = groups[j];
      const n1 = g1.values.length;
      const n2 = g2.values.length;
      const m1 = mean(g1.values);
      const m2 = mean(g2.values);
      const v1 = variance(g1.values);
      const v2 = variance(g2.values);
      const diff = m1 - m2;
      const se2 = Math.sqrt((v1 / n1 + v2 / n2) / 2);
      const q = se2 === 0 ? 0 : Math.abs(diff) / se2;
      const dfNum = (v1 / n1 + v2 / n2) ** 2;
      const dfDen = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1);
      const df = dfDen > 0 ? dfNum / dfDen : n1 + n2 - 2;
      const pValue = pValueStudentizedRangeApprox(q, k, df);
      const tCrit = tDistQuantile(1 - alpha / (k * (k - 1)), df);
      const ciHalf = tCrit * se2;
      results.push({
        group1: g1.label,
        group2: g2.label,
        meanDiff: diff,
        se: se2,
        statistic: q,
        pValue,
        pValueAdj: pValue,
        ci: [diff - ciHalf, diff + ciHalf],
        significant: pValue < alpha
      });
    }
  }
  return results;
}
function dunnTest(groups, method = "bonferroni") {
  const k = groups.length;
  const allValues = groups.flatMap((g) => [...g.values]);
  const n = allValues.length;
  const allRanks = rank(allValues);
  let offset = 0;
  const groupRankMeans = [];
  const groupNs = [];
  for (const g of groups) {
    const gn = g.values.length;
    let Rj = 0;
    for (let i = 0; i < gn; i++) Rj += allRanks[offset + i];
    groupRankMeans.push(Rj / gn);
    groupNs.push(gn);
    offset += gn;
  }
  const tieCounts = /* @__PURE__ */ new Map();
  for (const v of allValues) tieCounts.set(v, (tieCounts.get(v) ?? 0) + 1);
  let C = 0;
  for (const t of tieCounts.values()) C += t * t * t - t;
  const tieAdj = C / (12 * (n - 1));
  const baseVar = n * (n + 1) / 12 - tieAdj;
  const rawPValues = [];
  const rawZs = [];
  const rawSEs = [];
  const pairs = [];
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const diff = (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0);
      const se2 = Math.sqrt(
        baseVar * (1 / (groupNs[i] ?? 1) + 1 / (groupNs[j] ?? 1))
      );
      const z = se2 === 0 ? 0 : diff / se2;
      const p = normalSurvival(Math.abs(z));
      rawPValues.push(p);
      rawZs.push(z);
      rawSEs.push(se2);
      pairs.push({ i, j });
    }
  }
  const adjPValues = method === "holm" ? dunnStyleHolm(rawPValues) : method === "BH" ? dunnStyleBH(rawPValues) : adjustPValues(rawPValues, method);
  return pairs.map(({ i, j }, idx) => ({
    group1: groups[i].label,
    group2: groups[j].label,
    meanDiff: (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0),
    se: rawSEs[idx],
    statistic: rawZs[idx],
    pValue: rawPValues[idx],
    pValueAdj: adjPValues[idx],
    ci: [NaN, NaN],
    significant: (adjPValues[idx] ?? 1) < 0.05
  }));
}
function dunnStyleHolm(pValues) {
  const m = pValues.length;
  if (m === 0) return [];
  const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => a.p - b.p);
  const adj = new Array(m);
  for (let k = 0; k < m; k++) {
    const { p, i } = order[k];
    adj[i] = Math.min(1, p * (m - k));
  }
  return adj;
}
function dunnStyleBH(pValues) {
  const m = pValues.length;
  if (m === 0) return [];
  const order = pValues.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p);
  const adj = new Array(m);
  for (let k = 0; k < m; k++) {
    const { p, i } = order[k];
    const rank2 = m - k;
    adj[i] = Math.min(1, p * m / rank2);
  }
  return adj;
}

// src/stats/correlation.ts
function pearsonCorrelation(x, y, ciLevel = 0.95) {
  if (x.length !== y.length) throw new Error("pearsonCorrelation: arrays must have equal length");
  const n = x.length;
  if (n < 3) throw new Error("pearsonCorrelation: need at least 3 observations");
  const sdX = sd(x), sdY = sd(y);
  if (sdX === 0 || sdY === 0) throw new Error("pearsonCorrelation: zero variance in input");
  const r = cov(x, y) / (sdX * sdY);
  const rClamped = Math.max(-1, Math.min(1, r));
  const df = n - 2;
  const t = Math.abs(rClamped) === 1 ? Infinity : rClamped * Math.sqrt(df / (1 - rClamped * rClamped));
  const pValue = tDistPValue(t, df);
  const ci = fisherZCI(rClamped, n, ciLevel);
  return {
    testName: "Pearson r",
    statistic: roundTo(rClamped, 4),
    df,
    pValue: roundTo(pValue, 4),
    effectSize: {
      value: rClamped,
      name: "Pearson r",
      interpretation: interpretR(rClamped)
    },
    ci,
    ciLevel,
    n,
    formatted: formatCorrelation(rClamped, df, pValue, ci, "r", ciLevel)
  };
}
function fisherZCI(r, n, ciLevel) {
  const z = Math.log((1 + r) / (1 - r)) / 2;
  const se2 = 1 / Math.sqrt(n - 3);
  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2);
  const lo = z - zCrit * se2;
  const hi = z + zCrit * se2;
  return [
    Math.tanh(lo),
    Math.tanh(hi)
  ];
}
function spearmanCorrelation(x, y, ciLevel = 0.95) {
  if (x.length !== y.length) throw new Error("spearmanCorrelation: arrays must have equal length");
  const n = x.length;
  if (n < 3) throw new Error("spearmanCorrelation: need at least 3 observations");
  const rx = rank(x), ry = rank(y);
  const rhoResult = pearsonCorrelation(rx, ry, ciLevel);
  const rho = rhoResult.statistic;
  const hasTies = new Set(rx).size < n || new Set(ry).size < n;
  let pValue = rhoResult.pValue;
  if (!hasTies && n <= 1290) {
    let D = 0;
    for (let i = 0; i < n; i++) {
      const d = rx[i] - ry[i];
      D += d * d;
    }
    const q = (n * n * n - n) / 6;
    let p;
    if (D > q) {
      p = pSpearmanExact(Math.round(D), n, false);
    } else {
      p = pSpearmanExact(Math.round(D) + 2, n, true);
    }
    pValue = Math.min(2 * p, 1);
  }
  return {
    ...rhoResult,
    pValue,
    testName: "Spearman's \u03C1",
    effectSize: {
      ...rhoResult.effectSize,
      name: "Spearman's \u03C1"
    },
    formatted: formatCorrelation(rho, typeof rhoResult.df === "number" ? rhoResult.df : 0, pValue, rhoResult.ci, "\u03C1", ciLevel)
  };
}
function kendallTau(x, y, ciLevel = 0.95) {
  if (x.length !== y.length) throw new Error("kendallTau: arrays must have equal length");
  const n = x.length;
  if (n < 3) throw new Error("kendallTau: need at least 3 observations");
  let concordant = 0, discordant = 0;
  let tiesX = 0, tiesY = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dx = (x[i] ?? 0) - (x[j] ?? 0);
      const dy = (y[i] ?? 0) - (y[j] ?? 0);
      const sign = Math.sign(dx * dy);
      if (sign > 0) concordant++;
      else if (sign < 0) discordant++;
      if (dx === 0) tiesX++;
      if (dy === 0) tiesY++;
    }
  }
  const n2 = n * (n - 1) / 2;
  const tau = (concordant - discordant) / Math.sqrt((n2 - tiesX) * (n2 - tiesY));
  const xCounts = /* @__PURE__ */ new Map();
  const yCounts = /* @__PURE__ */ new Map();
  for (let i = 0; i < n; i++) {
    const xv = x[i];
    const yv = y[i];
    xCounts.set(xv, (xCounts.get(xv) ?? 0) + 1);
    yCounts.set(yv, (yCounts.get(yv) ?? 0) + 1);
  }
  const v0 = n * (n - 1) * (2 * n + 5);
  let vt = 0;
  let vu = 0;
  let v1t = 0;
  let v1u = 0;
  let v2t = 0;
  let v2u = 0;
  for (const t of xCounts.values()) {
    vt += t * (t - 1) * (2 * t + 5);
    v1t += t * (t - 1);
    v2t += t * (t - 1) * (t - 2);
  }
  for (const u of yCounts.values()) {
    vu += u * (u - 1) * (2 * u + 5);
    v1u += u * (u - 1);
    v2u += u * (u - 1) * (u - 2);
  }
  const sigma2 = (v0 - vt - vu) / 18 + v1t * v1u / (2 * n * (n - 1)) + v2t * v2u / (9 * n * (n - 1) * (n - 2));
  const z = (concordant - discordant) / Math.sqrt(sigma2);
  const hasTies = tiesX > 0 || tiesY > 0;
  let pValue;
  if (!hasTies && n < 50) {
    const S = concordant - discordant;
    const T = (S + n2) / 2;
    const pLower = pKendallExact(Math.round(T), n);
    const pUpper = 1 - pKendallExact(Math.round(T) - 1, n);
    pValue = Math.min(1, 2 * Math.min(pLower, pUpper));
  } else {
    pValue = 2 * normalSurvival(Math.abs(z));
  }
  const ci = fisherZCI(tau, n, ciLevel);
  const df = n - 2;
  return {
    testName: "Kendall's \u03C4",
    statistic: roundTo(tau, 4),
    df,
    pValue: roundTo(pValue, 4),
    effectSize: {
      value: tau,
      name: "Kendall's \u03C4",
      interpretation: interpretR(tau)
    },
    ci,
    ciLevel,
    n,
    formatted: formatCorrelation(tau, df, pValue, ci, "\u03C4", ciLevel)
  };
}
function partialCorrelation(x, y, controls) {
  if (x.length !== y.length) throw new Error("partialCorrelation: arrays must have equal length");
  const xRes = residualize(x, controls);
  const yRes = residualize(y, controls);
  const result = pearsonCorrelation(xRes, yRes);
  const n = x.length;
  const q = controls.length;
  const correctedDf = n - 2 - q;
  if (correctedDf < 1) return result;
  const r = result.effectSize.value;
  const t = Math.abs(r) === 1 ? Infinity : r * Math.sqrt(correctedDf / (1 - r * r));
  const pValue = tDistPValue(t, correctedDf);
  return {
    ...result,
    df: correctedDf,
    pValue: roundTo(pValue, 4),
    formatted: formatCorrelation(result.statistic, correctedDf, pValue, result.ci, "r", result.ciLevel)
  };
}
function residualize(y, predictors) {
  if (predictors.length === 0) return [...y];
  const n = y.length;
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((p) => p[i] ?? 0)])
  );
  const Xt = X.transpose();
  const XtX = Xt.multiply(X);
  const XtY = Xt.multiply(Matrix.colVec(y));
  const beta = XtX.inverse().multiply(XtY);
  const fitted = X.multiply(beta);
  return Array.from({ length: n }, (_, i) => (y[i] ?? 0) - fitted.get(i, 0));
}
function correlationMatrix(data, labels, method = "pearson") {
  const k = data.length;
  if (k < 2) throw new Error("correlationMatrix: need at least 2 variables");
  const n = data[0].length;
  const corrFn = method === "pearson" ? pearsonCorrelation : method === "spearman" ? spearmanCorrelation : kendallTau;
  const r = Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      if (i === j) return 1;
      if (j < i) return 0;
      try {
        return corrFn(data[i], data[j]).effectSize.value;
      } catch {
        return NaN;
      }
    })
  );
  const pValues = Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      if (i === j) return NaN;
      if (j < i) return 0;
      try {
        return corrFn(data[i], data[j]).pValue;
      } catch {
        return NaN;
      }
    })
  );
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < i; j++) {
      r[i][j] = r[j][i];
      pValues[i][j] = pValues[j][i];
    }
  }
  return {
    r,
    pValues,
    n,
    labels: labels ?? Array.from({ length: k }, (_, i) => `Var${i + 1}`)
  };
}

// src/stats/regression.ts
function fitOLS(X, y, coefNames, ciLevel = 0.95) {
  const n = y.length;
  const p = X.cols;
  const Xt = X.transpose();
  const XtX = Xt.multiply(X);
  const XtY = Xt.multiply(Matrix.colVec(y));
  const XtXInv = XtX.inverse();
  const betaM = XtXInv.multiply(XtY);
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0));
  const fitted = Array.from({ length: n }, (_, i) => {
    let val = 0;
    for (let j = 0; j < p; j++) val += X.get(i, j) * (beta[j] ?? 0);
    return val;
  });
  const residuals = y.map((v, i) => v - (fitted[i] ?? 0));
  const yMean = mean(y);
  const ss_res = residuals.reduce((s, r) => s + r * r, 0);
  const ss_tot = y.reduce((s, v) => s + (v - yMean) ** 2, 0);
  const r2 = ss_tot > 0 ? Math.max(0, 1 - ss_res / ss_tot) : 0;
  const adjR2 = 1 - (1 - r2) * (n - 1) / (n - p);
  const dfRes = n - p;
  if (dfRes <= 0) throw new Error("fitOLS: not enough degrees of freedom");
  const sigma2 = ss_res / dfRes;
  const covBeta = XtXInv.scale(sigma2);
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, dfRes);
  const coefficients = beta.map((b, i) => {
    const se2 = Math.sqrt(Math.max(0, covBeta.get(i, i)));
    const t = se2 === 0 ? 0 : b / se2;
    const pVal = tDistPValue(t, dfRes);
    const ci = [b - tCrit * se2, b + tCrit * se2];
    return {
      name: coefNames[i] ?? `\u03B2${i}`,
      estimate: roundTo(b, 6),
      se: roundTo(se2, 6),
      tValue: roundTo(t, 4),
      pValue: roundTo(pVal, 4),
      ci
    };
  });
  const dfModel = p - 1;
  const ss_reg = ss_tot - ss_res;
  const F = sigma2 === 0 || dfModel === 0 ? 0 : ss_reg / dfModel / sigma2;
  const fPValue = fDistPValue(F, dfModel, dfRes);
  const rssSafe = Math.max(ss_res, 1e-15);
  const logLik = -n / 2 * (Math.log(2 * Math.PI) + Math.log(rssSafe / n) + 1);
  const aic = -2 * logLik + 2 * (p + 1);
  const bic = -2 * logLik + Math.log(n) * (p + 1);
  const formatted = formatRegression(r2, adjR2, F, dfModel, dfRes, fPValue);
  return {
    coefficients,
    r2: roundTo(r2, 6),
    adjR2: roundTo(adjR2, 6),
    fStatistic: roundTo(F, 4),
    fDf: [dfModel, dfRes],
    fPValue: roundTo(fPValue, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    residuals,
    fitted,
    n,
    formatted
  };
}
function linearRegression(x, y, ciLevel = 0.95) {
  if (x.length !== y.length) throw new Error("linearRegression: arrays must have equal length");
  if (x.length < 3) throw new Error("linearRegression: need at least 3 observations");
  const n = x.length;
  const X = Matrix.fromArray(Array.from({ length: n }, (_, i) => [1, x[i] ?? 0]));
  return fitOLS(X, y, ["(Intercept)", "x"], ciLevel);
}
function multipleRegression(y, predictors, ciLevel = 0.95) {
  if (predictors.length === 0) throw new Error("multipleRegression: need at least 1 predictor");
  const n = y.length;
  for (const p of predictors) {
    if (p.values.length !== n) throw new Error(`multipleRegression: predictor '${p.name}' length mismatch`);
  }
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((p) => p.values[i] ?? 0)])
  );
  const names = ["(Intercept)", ...predictors.map((p) => p.name)];
  return fitOLS(X, y, names, ciLevel);
}
function polynomialRegression(x, y, degree, ciLevel = 0.95) {
  if (degree < 1) throw new Error("polynomialRegression: degree must be \u2265 1");
  if (x.length !== y.length) throw new Error("polynomialRegression: arrays must match length");
  const n = x.length;
  const X = Matrix.fromArray(
    Array.from(
      { length: n },
      (_, i) => [1, ...Array.from({ length: degree }, (_2, d) => (x[i] ?? 0) ** (d + 1))]
    )
  );
  const names = ["(Intercept)", ...Array.from({ length: degree }, (_, d) => `x^${d + 1}`)];
  return fitOLS(X, y, names, ciLevel);
}
function logisticRegression(y, predictors, ciLevel = 0.95, maxIter = 100, tol = 1e-8) {
  for (const v of y) {
    if (v !== 0 && v !== 1) throw new Error("logisticRegression: y must be 0 or 1");
  }
  const n = y.length;
  const p = predictors.length + 1;
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((pr) => pr.values[i] ?? 0)])
  );
  const names = ["(Intercept)", ...predictors.map((pr) => pr.name)];
  const computeEta = (b) => Array.from({ length: n }, (_, i) => {
    let v = 0;
    for (let j = 0; j < p; j++) v += X.get(i, j) * (b[j] ?? 0);
    return v;
  });
  const EPS_MU = 1e-15;
  const logistic = (e) => {
    const mu2 = 1 / (1 + Math.exp(-Math.min(700, Math.max(-700, e))));
    return Math.min(1 - EPS_MU, Math.max(EPS_MU, mu2));
  };
  const computeDeviance = (yArr, muArr) => {
    let dev = 0;
    for (let i = 0; i < n; i++) {
      const yi = yArr[i] ?? 0;
      const mi = muArr[i];
      dev += -2 * (yi * Math.log(Math.max(1e-15, mi)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - mi)));
    }
    return dev;
  };
  let beta = new Array(p).fill(0);
  const pMeanInit = Math.min(1 - 1e-6, Math.max(1e-6, mean([...y])));
  beta[0] = Math.log(pMeanInit / (1 - pMeanInit));
  let prevDeviance = Infinity;
  for (let iter = 0; iter < maxIter; iter++) {
    const eta2 = computeEta(beta);
    const mu2 = eta2.map(logistic);
    const w2 = mu2.map((m) => Math.max(1e-10, m * (1 - m)));
    const dev = computeDeviance(y, mu2);
    if (iter > 0 && Math.abs(dev - prevDeviance) / (0.1 + Math.abs(dev)) < tol) break;
    prevDeviance = dev;
    const z = Array.from({ length: n }, (_, i) => eta2[i] + ((y[i] ?? 0) - mu2[i]) / w2[i]);
    const sqrtW2 = w2.map(Math.sqrt);
    const Xw2 = Matrix.fromArray(
      Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_2, j) => X.get(i, j) * sqrtW2[i]))
    );
    const zw = z.map((zi, i) => zi * sqrtW2[i]);
    try {
      const Xwt = Xw2.transpose();
      const XwtXw = Xwt.multiply(Xw2);
      const XwtZw = Xwt.multiply(Matrix.colVec(zw));
      const betaNewM = XwtXw.inverse().multiply(XwtZw);
      const betaNew = Array.from({ length: p }, (_, j) => betaNewM.get(j, 0));
      let accepted = false;
      let stepScale = 1;
      for (let half = 0; half < 10; half++) {
        const betaTry = beta.map((b, j) => b + stepScale * (betaNew[j] - b));
        const etaTry = computeEta(betaTry);
        const muTry = etaTry.map(logistic);
        const trialDev = computeDeviance(y, muTry);
        if (isFinite(trialDev) && (trialDev <= dev + 1e-8 || half === 9)) {
          beta = betaTry;
          accepted = true;
          break;
        }
        stepScale *= 0.5;
      }
      if (!accepted) break;
    } catch {
      break;
    }
  }
  const eta = computeEta(beta);
  const mu = eta.map(logistic);
  const w = mu.map((m) => Math.max(1e-10, m * (1 - m)));
  const sqrtW = w.map(Math.sqrt);
  const Xw = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_2, j) => X.get(i, j) * sqrtW[i]))
  );
  let covBeta;
  try {
    covBeta = Xw.transpose().multiply(Xw).inverse();
  } catch {
    covBeta = Matrix.identity(p);
  }
  const zCrit = normalQuantile(1 - (1 - ciLevel) / 2);
  const coefficients = beta.map((b, i) => {
    const se2 = Math.sqrt(Math.max(0, covBeta.get(i, i)));
    const zStat = se2 === 0 ? 0 : b / se2;
    const pVal = 2 * (1 - normCDFLocal(Math.abs(zStat)));
    const ci = [b - zCrit * se2, b + zCrit * se2];
    return {
      name: names[i] ?? `\u03B2${i}`,
      estimate: roundTo(b, 10),
      se: roundTo(se2, 10),
      tValue: roundTo(zStat, 6),
      pValue: roundTo(pVal, 6),
      ci
    };
  });
  const logLik = mu.reduce((s, m, i) => {
    const yi = y[i] ?? 0;
    return s + yi * Math.log(Math.max(1e-15, m)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - m));
  }, 0);
  const pMeanRaw = mean([...y]);
  const pMean = Math.min(1 - 1e-12, Math.max(1e-12, pMeanRaw));
  const nullLogLik = n * (pMean * Math.log(Math.max(1e-15, pMean)) + (1 - pMean) * Math.log(Math.max(1e-15, 1 - pMean)));
  const r2 = Math.abs(nullLogLik) < 1e-12 ? NaN : 1 - logLik / nullLogLik;
  const aic = -2 * logLik + 2 * p;
  const bic = -2 * logLik + Math.log(n) * p;
  const residuals = y.map((v, i) => (v ?? 0) - (mu[i] ?? 0));
  return {
    coefficients,
    r2: roundTo(r2, 6),
    adjR2: roundTo(r2, 6),
    // McFadden's for logistic
    fStatistic: NaN,
    fDf: [p - 1, n - p],
    fPValue: NaN,
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    residuals,
    fitted: mu,
    n,
    formatted: `McFadden R\xB2 = ${roundTo(r2, 3)}, AIC = ${roundTo(aic, 1)}`
  };
}
function normCDFLocal(z) {
  const x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const erf2 = 1 - poly * Math.exp(-x * x);
  return 0.5 * (1 + (z >= 0 ? erf2 : -erf2));
}
function regressionDiagnostics(result, predictors) {
  const n = result.n;
  const p = result.coefficients.length;
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((pr) => pr.values[i] ?? 0)])
  );
  const Xt = X.transpose();
  let XtXInv;
  try {
    XtXInv = Xt.multiply(X).inverse();
  } catch {
    XtXInv = Matrix.identity(p);
  }
  const hat = X.multiply(XtXInv).multiply(Xt);
  const leverage = Array.from({ length: n }, (_, i) => hat.get(i, i));
  const mse = result.residuals.reduce((s, r) => s + r * r, 0) / (n - p);
  const standardizedResiduals = result.residuals.map((r, i) => {
    const denom = Math.sqrt(mse * (1 - (leverage[i] ?? 0)));
    return denom === 0 ? 0 : r / denom;
  });
  const cooksDistance = result.residuals.map((r, i) => {
    const h = leverage[i] ?? 0;
    return r * r * h / (p * mse * (1 - h) ** 2);
  });
  const vif = predictors.map((_, j) => {
    const otherPreds = predictors.filter((__, k) => k !== j);
    if (otherPreds.length === 0) return 1;
    const xj = predictors[j].values;
    const others = otherPreds.map((p2) => ({ name: p2.name, values: p2.values }));
    try {
      const res = multipleRegression(xj, others);
      return 1 / Math.max(1e-10, 1 - res.r2);
    } catch {
      return NaN;
    }
  });
  return { leverage, cooksDistance, standardizedResiduals, vif };
}

// src/stats/preprocess.ts
function preprocessData(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("preprocessData: data cannot be empty");
  const d = data[0].length;
  const method = options?.method ?? "none";
  if (method === "log") {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < d; j++) {
        if (data[i][j] <= 0) {
          throw new Error(`preprocessData: log requires all values > 0, found ${data[i][j]} at row ${i}, col ${j}`);
        }
      }
    }
    const transformed = data.map((row) => row.map((v) => Math.log(v)));
    const colMeans2 = computeColMeans(transformed, n, d);
    const colSDs2 = computeColSDs(transformed, colMeans2, n, d);
    return { data: transformed, colMeans: colMeans2, colSDs: colSDs2, method, centered: false, scaled: false };
  }
  if (method === "sqrt") {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < d; j++) {
        if (data[i][j] < 0) {
          throw new Error(`preprocessData: sqrt requires all values >= 0, found ${data[i][j]} at row ${i}, col ${j}`);
        }
      }
    }
    const transformed = data.map((row) => row.map((v) => Math.sqrt(v)));
    const colMeans2 = computeColMeans(transformed, n, d);
    const colSDs2 = computeColSDs(transformed, colMeans2, n, d);
    return { data: transformed, colMeans: colMeans2, colSDs: colSDs2, method, centered: false, scaled: false };
  }
  const colMeans = computeColMeans(data, n, d);
  const colSDs = computeColSDs(data, colMeans, n, d);
  if (method === "none") {
    return { data, colMeans, colSDs, method, centered: false, scaled: false };
  }
  if (method === "center") {
    const centered = data.map(
      (row) => row.map((v, j) => v - colMeans[j])
    );
    return { data: centered, colMeans, colSDs, method, centered: true, scaled: false };
  }
  const safeSDs = colSDs.map((s) => s === 0 ? 1 : s);
  const standardized = data.map(
    (row) => row.map((v, j) => (v - colMeans[j]) / safeSDs[j])
  );
  return { data: standardized, colMeans, colSDs: safeSDs, method, centered: true, scaled: true };
}
function inverseTransform(data, params) {
  const { method, colMeans, colSDs } = params;
  if (method === "none") return data;
  if (method === "log") {
    return data.map((row) => row.map((v) => Math.exp(v)));
  }
  if (method === "sqrt") {
    return data.map((row) => row.map((v) => v * v));
  }
  if (method === "center") {
    return data.map(
      (row) => row.map((v, j) => v + colMeans[j])
    );
  }
  return data.map(
    (row) => row.map((v, j) => v * colSDs[j] + colMeans[j])
  );
}
function computeColMeans(data, n, d) {
  return Array.from({ length: d }, (_, j) => {
    let sum = 0;
    for (let i = 0; i < n; i++) sum += data[i][j];
    return sum / n;
  });
}
function computeColSDs(data, colMeans, n, d) {
  if (n < 2) return new Array(d).fill(0);
  return Array.from({ length: d }, (_, j) => {
    let ss2 = 0;
    for (let i = 0; i < n; i++) {
      const diff = data[i][j] - colMeans[j];
      ss2 += diff * diff;
    }
    return Math.sqrt(ss2 / (n - 1));
  });
}

// src/stats/pca.ts
function runPCA(data, nComponents, scale = true) {
  if (data.length < 2) throw new Error("runPCA: need at least 2 observations");
  const n = data.length;
  const k = data[0].length;
  if (k < 2) throw new Error("runPCA: need at least 2 variables");
  const pp = preprocessData(data, { method: scale ? "standardize" : "center" });
  const X = Matrix.fromArray(pp.data);
  const Xt = X.transpose();
  const R = Xt.multiply(X).scale(1 / (n - 1));
  const { values: eigVals, vectors: eigVecs } = R.eigen();
  const nc = nComponents ?? Math.min(n - 1, k);
  const eigenvalues = eigVals.slice(0, nc);
  let totalVar = 0;
  for (let i = 0; i < k; i++) totalVar += R.get(i, i);
  const loadings = Array.from(
    { length: k },
    (_, varIdx) => Array.from({ length: nc }, (_2, compIdx) => eigVecs.get(varIdx, compIdx))
  );
  const Vk = Matrix.fromArray(
    Array.from({ length: k }, (_, i) => Array.from({ length: nc }, (_2, j) => eigVecs.get(i, j)))
  );
  const scoresM = X.multiply(Vk);
  const scores = Array.from(
    { length: n },
    (_, i) => Array.from({ length: nc }, (_2, j) => scoresM.get(i, j))
  );
  const varianceExplained = eigenvalues.map((e) => totalVar > 0 ? e / totalVar : 0);
  const cumulativeVariance = varianceExplained.reduce((acc, v, i) => {
    acc.push((acc[i - 1] ?? 0) + v);
    return acc;
  }, []);
  return {
    loadings,
    scores,
    eigenvalues,
    varianceExplained,
    cumulativeVariance,
    nComponents: nc
  };
}
function varimaxRotation(loadings, maxIter = 1e3, tol = 1e-5, normalize = true) {
  const p = loadings.length;
  const nc = loadings[0].length;
  if (nc < 2) {
    return {
      rotatedLoadings: loadings.map((row) => [...row]),
      rotationMatrix: [[1]]
    };
  }
  let X = loadings.map((row) => [...row]);
  let sc = null;
  if (normalize) {
    sc = X.map((row) => Math.sqrt(row.reduce((s, v) => s + v * v, 0)));
    X = X.map((row, i) => {
      const h = sc[i];
      return h > 0 ? row.map((v) => v / h) : [...row];
    });
  }
  let T = Array.from(
    { length: nc },
    (_, i) => Array.from({ length: nc }, (_2, j) => i === j ? 1 : 0)
  );
  let d = 0;
  for (let iter = 0; iter < maxIter; iter++) {
    const Z2 = X.map((row) => {
      const out = new Array(nc).fill(0);
      for (let j = 0; j < nc; j++) {
        let s = 0;
        for (let l = 0; l < nc; l++) s += row[l] * T[l][j];
        out[j] = s;
      }
      return out;
    });
    const colSumSq = new Array(nc).fill(0);
    for (let i = 0; i < p; i++) {
      for (let j = 0; j < nc; j++) colSumSq[j] += Z2[i][j] * Z2[i][j];
    }
    for (let j = 0; j < nc; j++) colSumSq[j] /= p;
    const target = Z2.map(
      (row) => row.map((v, j) => v * v * v - v * colSumSq[j])
    );
    const B = Array.from({ length: nc }, (_, i) => {
      const out = new Array(nc).fill(0);
      for (let j = 0; j < nc; j++) {
        let s = 0;
        for (let r = 0; r < p; r++) s += X[r][i] * target[r][j];
        out[j] = s;
      }
      return out;
    });
    const Bm = Matrix.fromArray(B);
    const BtB = Bm.transpose().multiply(Bm);
    const { values: eigVals, vectors: eigVecs } = BtB.eigen();
    let dNew = 0;
    const invSqrt = eigVals.map((lambda) => {
      const sv = Math.sqrt(Math.max(0, lambda));
      dNew += sv;
      return sv > 1e-15 ? 1 / sv : 0;
    });
    const invSqrtBtB = Array.from({ length: nc }, (_, i) => {
      const row = new Array(nc).fill(0);
      for (let j = 0; j < nc; j++) {
        let s = 0;
        for (let l = 0; l < nc; l++) s += eigVecs.get(i, l) * invSqrt[l] * eigVecs.get(j, l);
        row[j] = s;
      }
      return row;
    });
    const Tnew = Array.from({ length: nc }, (_, i) => {
      const row = new Array(nc).fill(0);
      for (let j = 0; j < nc; j++) {
        let s = 0;
        for (let l = 0; l < nc; l++) s += B[i][l] * invSqrtBtB[l][j];
        row[j] = s;
      }
      return row;
    });
    T = Tnew;
    const dPast = d;
    d = dNew;
    if (d < dPast * (1 + tol)) break;
  }
  let Z = X.map((row) => {
    const out = new Array(nc).fill(0);
    for (let j = 0; j < nc; j++) {
      let s = 0;
      for (let l = 0; l < nc; l++) s += row[l] * T[l][j];
      out[j] = s;
    }
    return out;
  });
  if (normalize && sc) {
    Z = Z.map((row, i) => {
      const h = sc[i];
      return row.map((v) => v * h);
    });
  }
  return { rotatedLoadings: Z, rotationMatrix: T };
}
function screeData(pca) {
  return {
    components: Array.from({ length: pca.nComponents }, (_, i) => i + 1),
    eigenvalues: pca.eigenvalues,
    varianceExplained: pca.varianceExplained,
    cumulativeVariance: pca.cumulativeVariance
  };
}

// src/stats/mixed.ts
function remlProfileLogLik(logPsi, y, X, Z) {
  const psi = Math.exp(logPsi);
  const n = y.length;
  const q = Z.cols;
  const p = X.cols;
  const ZtZ = Z.transpose().multiply(Z);
  const Dmat = ZtZ.add(Matrix.identity(q).scale(1 / psi));
  let DInv;
  let logDetD;
  try {
    DInv = Dmat.inverse();
    logDetD = Dmat.logDet();
  } catch {
    return { negLogLik: Infinity, sigmae2: 0, sigmab2: 0 };
  }
  const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose());
  const VpsiInv = Matrix.fromArray(
    Array.from(
      { length: n },
      (_, i) => Array.from(
        { length: n },
        (_2, j) => (i === j ? 1 : 0) - ZDinvZt.get(i, j)
      )
    )
  );
  const logDetVpsi = q * Math.log(psi) + logDetD;
  const XtVinv = X.transpose().multiply(VpsiInv);
  const XtVinvX = XtVinv.multiply(X);
  let XtVinvXInv;
  let logDetXVX;
  try {
    XtVinvXInv = XtVinvX.inverse();
    logDetXVX = XtVinvX.logDet();
  } catch {
    return { negLogLik: Infinity, sigmae2: 0, sigmab2: 0 };
  }
  const XtVinvY = XtVinv.multiply(Matrix.colVec(y));
  const beta = XtVinvXInv.multiply(XtVinvY);
  const Xbeta = X.multiply(beta);
  const e = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0));
  const eM = Matrix.colVec(e);
  const quadForm = eM.transpose().multiply(VpsiInv).multiply(eM).get(0, 0);
  const sigmae2 = Math.max(1e-8, quadForm / (n - p));
  const sigmab2 = psi * sigmae2;
  const reml = -0.5 * ((n - p) * Math.log(sigmae2) + logDetVpsi + logDetXVX);
  return { negLogLik: -reml, sigmae2, sigmab2 };
}
function runLMM(input) {
  const { outcome: y, fixedPredictors, groupId, ciLevel = 0.95 } = input;
  const n = y.length;
  if (n < 5) throw new Error("runLMM: need at least 5 observations");
  if (groupId.length !== n) throw new Error("runLMM: groupId must have same length as outcome");
  const groupLevels = [...new Set(groupId)];
  const nGroups = groupLevels.length;
  if (nGroups < 2) throw new Error("runLMM: need at least 2 groups");
  const predNames = Object.keys(fixedPredictors);
  const p = predNames.length + 1;
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [
      1,
      ...predNames.map((name) => (fixedPredictors[name] ?? [])[i] ?? 0)
    ])
  );
  const Z = Matrix.fromArray(
    Array.from(
      { length: n },
      (_, i) => groupLevels.map((g) => groupId[i] === g ? 1 : 0)
    )
  );
  const objFn = (theta) => remlProfileLogLik(theta[0] ?? 0, y, X, Z).negLogLik;
  const starts = [-4, -2, 0, 2, 4];
  let optResult = nelderMead(objFn, [starts[0]], { maxIter: 1e3, tol: 1e-8 });
  for (let si = 1; si < starts.length; si++) {
    const cand = nelderMead(objFn, [starts[si]], { maxIter: 1e3, tol: 1e-8 });
    if (cand.fval < optResult.fval) optResult = cand;
  }
  const finalModel = remlProfileLogLik(optResult.x[0] ?? 0, y, X, Z);
  const sigmab2 = finalModel.sigmab2;
  const sigmae2 = finalModel.sigmae2;
  const scale = sigmab2 / sigmae2;
  const ZtZ = Z.transpose().multiply(Z);
  let VinvScaled;
  if (scale < 1e-10) {
    VinvScaled = Matrix.identity(n);
  } else {
    const Dmat = ZtZ.add(Matrix.identity(nGroups).scale(1 / scale));
    let DInv;
    try {
      DInv = Dmat.inverse();
      const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose());
      VinvScaled = Matrix.fromArray(
        Array.from(
          { length: n },
          (_, i) => Array.from(
            { length: n },
            (_2, j) => (i === j ? 1 : 0) - ZDinvZt.get(i, j)
          )
        )
      );
    } catch {
      VinvScaled = Matrix.identity(n);
    }
  }
  const Vinv = VinvScaled.scale(1 / sigmae2);
  const Xt = X.transpose();
  const XtVinv = Xt.multiply(Vinv);
  const XtVinvX = XtVinv.multiply(X);
  let XtVinvXInv;
  try {
    XtVinvXInv = XtVinvX.inverse();
  } catch {
    XtVinvXInv = Matrix.identity(p);
  }
  const XtVinvY = XtVinv.multiply(Matrix.colVec([...y]));
  const betaM = XtVinvXInv.multiply(XtVinvY);
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0));
  const df = Math.max(1, n - p - nGroups + 1);
  const tCrit = tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const covBeta = XtVinvXInv.scale(sigmae2);
  const fixedEffectNames = ["(Intercept)", ...predNames];
  const fixedEffects = beta.map((b, i) => {
    const seVal = Math.sqrt(Math.max(0, covBeta.get(i, i)));
    const t = seVal === 0 ? 0 : b / seVal;
    const pVal = tDistPValue(t, df);
    return {
      name: fixedEffectNames[i] ?? `\u03B2${i}`,
      estimate: roundTo(b, 6),
      se: roundTo(seVal, 6),
      tValue: roundTo(t, 4),
      pValue: roundTo(pVal, 4),
      ci: [roundTo(b - tCrit * seVal, 6), roundTo(b + tCrit * seVal, 6)]
    };
  });
  const icc = sigmab2 / (sigmab2 + sigmae2);
  const remlConst = 0.5 * (n - p) * (1 + Math.log(2 * Math.PI));
  const logLik = -finalModel.negLogLik - remlConst;
  const aic = -2 * logLik + 2 * (p + 2);
  const bic = -2 * logLik + Math.log(n) * (p + 2);
  const formatted = formatLMM(icc, aic, bic, logLik);
  return {
    fixedEffects,
    varianceComponents: {
      intercept: roundTo(sigmab2, 6),
      residual: roundTo(sigmae2, 6)
    },
    icc: roundTo(icc, 6),
    logLik: roundTo(logLik, 4),
    aic: roundTo(aic, 2),
    bic: roundTo(bic, 2),
    nObs: n,
    nGroups,
    formatted
  };
}
function computeBLUPs(input, result) {
  const { outcome: y, fixedPredictors, groupId } = input;
  const n = y.length;
  const groupLevels = [...new Set(groupId)];
  const predNames = Object.keys(fixedPredictors);
  const sigmab2 = result.varianceComponents.intercept;
  const sigmae2 = result.varianceComponents.residual;
  const X = Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predNames.map((name) => (fixedPredictors[name] ?? [])[i] ?? 0)])
  );
  const beta = result.fixedEffects.map((fe) => fe.estimate);
  const Xbeta = X.multiply(Matrix.colVec(beta));
  const residuals = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0));
  const psi = sigmab2 / sigmae2;
  return groupLevels.map((g) => {
    const indices = Array.from({ length: n }, (_, i) => i).filter((i) => groupId[i] === g);
    const sumResid = indices.reduce((s, i) => s + (residuals[i] ?? 0), 0);
    const nj = indices.length;
    const blup = psi / (1 + psi * nj) * sumResid;
    return { group: g, blup: roundTo(blup, 6) };
  });
}

// src/stats/analyze.ts
var DEFAULTS = {
  ciLevel: 0.95,
  paired: false,
  pAdjMethod: "holm",
  forceTest: void 0,
  equalVariances: false,
  normalityAlpha: 0.05
};
function detectFieldType(values) {
  if (values.length === 0) return "categorical";
  const allFiniteNumber = values.every((v) => typeof v === "number" && isFinite(v));
  if (allFiniteNumber) {
    const nums = values;
    const unique2 = new Set(nums);
    if (unique2.size === 2) {
      const sorted = [...unique2].sort((a, b) => a - b);
      if (sorted[0] === 0 && sorted[1] === 1) return "binary";
    }
    return "numeric";
  }
  const unique = new Set(values);
  return unique.size === 2 ? "binary" : "categorical";
}
function splitGroups(outcome, labels) {
  if (outcome.length !== labels.length) {
    throw new Error(
      `analyze: outcome length (${outcome.length}) !== predictor length (${labels.length})`
    );
  }
  const map = /* @__PURE__ */ new Map();
  for (let i = 0; i < labels.length; i++) {
    const key = String(labels[i]);
    if (!map.has(key)) map.set(key, []);
    map.get(key).push(outcome[i]);
  }
  return [...map.entries()].sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0).map(([label, values]) => ({ label, values }));
}
function checkNormality(groups, alpha = 0.05) {
  const results = groups.map((g) => {
    const n = g.values.length;
    if (n < 3 || n > 50) {
      return { group: g.label, W: 1, p: 1 };
    }
    const sw = shapiroWilk(g.values);
    return { group: g.label, W: sw.statistic, p: sw.pValue };
  });
  const allNormal = results.every((r) => r.p >= alpha);
  return { allNormal, results };
}
function selectTest(outcomeFt, predictorFt, groups, opts) {
  if (opts.forceTest) return opts.forceTest;
  if (outcomeFt === "numeric") {
    if (predictorFt === void 0) return "describe-only";
    if (predictorFt === "binary") {
      const { allNormal } = checkNormality(groups, opts.normalityAlpha);
      if (opts.paired) {
        return allNormal ? "t-test-paired" : "wilcoxon";
      }
      return allNormal ? "t-test-independent" : "mann-whitney";
    }
    if (predictorFt === "categorical") {
      const { allNormal } = checkNormality(groups, opts.normalityAlpha);
      return allNormal ? "one-way-anova" : "kruskal-wallis";
    }
    throw new Error(
      `analyze: unsupported predictor type '${predictorFt}' for numeric outcome. Use 'binary' or 'categorical'.`
    );
  }
  return "chi-square";
}
function analyze(outcome, predictor, opts) {
  const options = { ...DEFAULTS, ...opts };
  if (outcome.type === "binary" || outcome.type === "categorical") {
    if (predictor === void 0) {
      throw new Error(
        `analyze: categorical/binary outcome requires a predictor for frequency tests`
      );
    }
    const groupOutcome = outcome;
    const groupPredictor = predictor;
    const { table } = contingencyTable(
      groupOutcome.values,
      groupPredictor.values
    );
    const testName = options.forceTest ?? "chi-square";
    let result2;
    if (testName === "fisher") {
      if (table.length !== 2 || table[0].length !== 2) {
        throw new Error(`analyze: Fisher's exact test requires a 2\xD72 table`);
      }
      result2 = fisherExactTest(table[0][0], table[0][1], table[1][0], table[1][1]);
    } else {
      if (!options.forceTest && table.length === 2 && table[0].length === 2) {
        const n = table[0][0] + table[0][1] + table[1][0] + table[1][1];
        const r0 = table[0][0] + table[0][1];
        const r1 = n - r0;
        const c0 = table[0][0] + table[1][0];
        const c1 = n - c0;
        const anyLowExpected = r0 * c0 / n < 5 || r0 * c1 / n < 5 || r1 * c0 / n < 5 || r1 * c1 / n < 5;
        if (anyLowExpected) {
          result2 = fisherExactTest(table[0][0], table[0][1], table[1][0], table[1][1]);
          return {
            test: result2.testName,
            outcome: outcome.name,
            predictor: predictor.name,
            result: result2
          };
        }
      }
      result2 = chiSquareTest(table);
    }
    return {
      test: result2.testName,
      outcome: outcome.name,
      predictor: predictor.name,
      result: result2
    };
  }
  const numericOutcome = outcome;
  const numericValues = numericOutcome.values;
  if (predictor === void 0) {
    const desc = describe(numericValues, options.ciLevel);
    const dummyResult = {
      testName: "Descriptive statistics",
      statistic: NaN,
      df: 0,
      pValue: NaN,
      effectSize: { value: NaN, name: "none", interpretation: "negligible" },
      ci: [NaN, NaN],
      ciLevel: options.ciLevel,
      n: numericValues.length,
      formatted: desc.formatted
    };
    return {
      test: "descriptive",
      outcome: outcome.name,
      result: dummyResult,
      descriptives: [desc]
    };
  }
  const groups = splitGroups(numericValues, predictor.values);
  const normalityCheck = checkNormality(groups, options.normalityAlpha);
  const selectedTest = selectTest(outcome.type, predictor.type, groups, options);
  let result;
  let posthoc;
  if (selectedTest === "t-test-independent" || selectedTest === "mann-whitney") {
    if (groups.length !== 2) {
      throw new Error(
        `analyze: '${selectedTest}' requires exactly 2 groups, got ${groups.length}`
      );
    }
    const [g1, g2] = groups;
    if (selectedTest === "t-test-independent") {
      result = tTestIndependent(g1.values, g2.values, options.equalVariances, options.ciLevel);
    } else {
      result = mannWhitneyU(g1.values, g2.values);
    }
  } else if (selectedTest === "t-test-paired" || selectedTest === "wilcoxon") {
    if (groups.length !== 2) {
      throw new Error(
        `analyze: '${selectedTest}' requires exactly 2 groups, got ${groups.length}`
      );
    }
    const [g1, g2] = groups;
    if (g1.values.length !== g2.values.length) {
      throw new Error(
        `analyze: paired=true but group sizes are unequal (${g1.label}: ${g1.values.length}, ${g2.label}: ${g2.values.length})`
      );
    }
    if (selectedTest === "t-test-paired") {
      result = tTestPaired(g1.values, g2.values, options.ciLevel);
    } else {
      result = wilcoxonSignedRank(g1.values, g2.values);
    }
  } else if (selectedTest === "one-way-anova") {
    const anovaResult = oneWayANOVA(groups);
    result = anovaResult;
    posthoc = tukeyHSD(groups, anovaResult.msWithin, anovaResult.dfWithin, options.ciLevel);
  } else if (selectedTest === "kruskal-wallis") {
    result = kruskalWallis(groups);
    posthoc = dunnTest(groups, options.pAdjMethod);
  } else if (selectedTest === "chi-square" || selectedTest === "fisher") {
    throw new Error(
      `analyze: '${selectedTest}' is not applicable to numeric outcomes`
    );
  } else {
    throw new Error(`analyze: unknown test '${selectedTest}'`);
  }
  const descriptives = groups.map(
    (g) => describe(g.values, options.ciLevel)
  );
  return {
    test: result.testName,
    outcome: outcome.name,
    predictor: predictor.name,
    result,
    descriptives,
    // Only include optional fields when they have a value (exactOptionalPropertyTypes)
    ...posthoc !== void 0 && { posthoc },
    normality: normalityCheck.results
  };
}

// src/stats/clustering.ts
var LOG_2PI = Math.log(2 * Math.PI);
var MIN_PROB = 1e-300;
var PRNG = class {
  state;
  constructor(seed) {
    this.state = seed >>> 0;
  }
  next() {
    this.state = this.state + 2654435769 | 0;
    let t = this.state ^ this.state >>> 16;
    t = Math.imul(t, 569420461);
    t = t ^ t >>> 15;
    t = Math.imul(t, 1935289751);
    t = t ^ t >>> 15;
    return (t >>> 0) / 4294967296;
  }
};
function logSumExp(arr) {
  let max = -Infinity;
  for (let i = 0; i < arr.length; i++) {
    if (arr[i] > max) max = arr[i];
  }
  if (max === -Infinity) return -Infinity;
  let sum = 0;
  for (let i = 0; i < arr.length; i++) sum += Math.exp(arr[i] - max);
  return max + Math.log(sum);
}
function computeRawEntropy(resp, k) {
  let ent = 0;
  for (let i = 0; i < resp.length; i++) {
    for (let j = 0; j < k; j++) {
      const z = resp[i][j];
      if (z > MIN_PROB) ent -= z * Math.log(z);
    }
  }
  return ent;
}
function computeNormalizedEntropy(resp, k) {
  if (k <= 1) return 1;
  const rawE = computeRawEntropy(resp, k);
  const n = resp.length;
  const denom = n * Math.log(k);
  return denom > 0 ? 1 - rawE / denom : 1;
}
function computeAvePP(resp, k) {
  const sums = new Float64Array(k);
  const counts = new Float64Array(k);
  for (let i = 0; i < resp.length; i++) {
    let maxP = -1, best = 0;
    for (let j = 0; j < k; j++) {
      if (resp[i][j] > maxP) {
        maxP = resp[i][j];
        best = j;
      }
    }
    sums[best] += maxP;
    counts[best]++;
  }
  return Array.from(sums).map((s, i) => counts[i] > 0 ? s / counts[i] : 0);
}
function kMeansPlusPlus(data, k, rng) {
  const n = data.length;
  const d = data[0].length;
  const means = [[...data[Math.floor(rng.next() * n)]]];
  const dists = new Float64Array(n).fill(Infinity);
  for (let j = 1; j < k; j++) {
    const lastMean = means[j - 1];
    let sumSqDist = 0;
    for (let i = 0; i < n; i++) {
      const pt = data[i];
      let dSq = 0;
      for (let dim = 0; dim < d; dim++) {
        const diff = pt[dim] - lastMean[dim];
        dSq += diff * diff;
      }
      if (dSq < dists[i]) dists[i] = dSq;
      sumSqDist += dists[i];
    }
    let target = rng.next() * sumSqDist;
    let cumulative = 0;
    for (let i = 0; i < n; i++) {
      cumulative += dists[i];
      if (cumulative >= target) {
        means.push([...data[i]]);
        break;
      }
    }
    if (means.length <= j) means.push([...data[n - 1]]);
  }
  return means;
}
function mvnLogPdf(x, mu, uFlat, eigenvals, d) {
  let mahal = 0;
  let logDet = 0;
  for (let i = 0; i < d; i++) {
    let yi = 0;
    for (let j = 0; j < d; j++) {
      yi += uFlat[j * d + i] * (x[j] - mu[j]);
    }
    const eig = eigenvals[i];
    mahal += yi * yi / eig;
    logDet += Math.log(eig);
  }
  return -0.5 * (d * LOG_2PI + logDet + mahal);
}
function fitGMM(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("fitGMM: data cannot be empty");
  const d = data[0].length;
  const k = options.k;
  if (k < 1) throw new Error("fitGMM: k must be >= 1");
  if (k > n) throw new Error("fitGMM: k cannot exceed n");
  const modelType = options.model ?? "VVV";
  const rng = new PRNG(options.seed ?? 42);
  const tol = options.tol ?? 1e-6;
  const maxIter = options.maxIter ?? 200;
  const regCovar = options.regCovar ?? 1e-6;
  const means = kMeansPlusPlus(data, k, rng);
  const weights = new Float64Array(k).fill(1 / k);
  const resp = Array.from({ length: n }, () => new Float64Array(k));
  const I = Matrix.identity(d);
  const iFlat = I.toFlat();
  const uFlats = Array.from({ length: k }, () => [...iFlat]);
  const eigenvals = Array.from({ length: k }, () => new Array(d).fill(1));
  let globalVar = regCovar;
  const colMeans = Array.from({ length: d }, (_, j) => {
    let s = 0;
    for (let i = 0; i < n; i++) s += data[i][j];
    return s / n;
  });
  for (let i = 0; i < n; i++) {
    for (let dim = 0; dim < d; dim++) {
      globalVar += (data[i][dim] - colMeans[dim]) ** 2 / (n * d);
    }
  }
  for (let j = 0; j < k; j++) eigenvals[j].fill(globalVar);
  let prevLogL = -Infinity;
  let converged = false;
  let iter = 0;
  for (; iter < maxIter; iter++) {
    let currentLogL = 0;
    const logLiks = new Float64Array(k);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) {
        logLiks[j] = Math.log(Math.max(weights[j], MIN_PROB)) + mvnLogPdf(data[i], means[j], uFlats[j], eigenvals[j], d);
      }
      const marg = logSumExp(logLiks);
      currentLogL += marg;
      for (let j = 0; j < k; j++) {
        resp[i][j] = Math.exp(logLiks[j] - marg);
      }
    }
    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true;
      break;
    }
    prevLogL = currentLogL;
    const Nk = new Float64Array(k);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) Nk[j] += resp[i][j];
    }
    const needPool = modelType === "EEE" || modelType === "EEI" || modelType === "EII";
    const pool = needPool ? Array.from({ length: d }, () => new Array(d).fill(0)) : null;
    for (let j = 0; j < k; j++) {
      const nk = Math.max(Nk[j], MIN_PROB);
      weights[j] = Nk[j] / n;
      const mu = means[j];
      mu.fill(0);
      for (let i = 0; i < n; i++) {
        const r = resp[i][j];
        for (let dim = 0; dim < d; dim++) mu[dim] += r * data[i][dim];
      }
      for (let dim = 0; dim < d; dim++) mu[dim] /= nk;
      const cov2 = Array.from({ length: d }, () => new Array(d).fill(0));
      for (let i = 0; i < n; i++) {
        const r = resp[i][j];
        for (let r_idx = 0; r_idx < d; r_idx++) {
          const dr = data[i][r_idx] - mu[r_idx];
          for (let c_idx = r_idx; c_idx < d; c_idx++) {
            const val = r * dr * (data[i][c_idx] - mu[c_idx]);
            cov2[r_idx][c_idx] += val;
          }
        }
      }
      for (let r_idx = 0; r_idx < d; r_idx++) {
        for (let c_idx = r_idx; c_idx < d; c_idx++) {
          cov2[r_idx][c_idx] /= nk;
          cov2[c_idx][r_idx] = cov2[r_idx][c_idx];
        }
        cov2[r_idx][r_idx] += regCovar;
      }
      if (pool) {
        const w = weights[j];
        for (let r_idx = 0; r_idx < d; r_idx++) {
          for (let c_idx = 0; c_idx < d; c_idx++) {
            pool[r_idx][c_idx] += w * cov2[r_idx][c_idx];
          }
        }
      }
      if (modelType === "VVV") {
        const { values, vectors } = Matrix.fromArray(cov2).eigen();
        eigenvals[j] = values.map((v) => Math.max(v, regCovar));
        uFlats[j] = vectors.toFlat();
      } else if (modelType === "VVI") {
        for (let dim = 0; dim < d; dim++) eigenvals[j][dim] = cov2[dim][dim];
        uFlats[j] = [...iFlat];
      } else if (modelType === "VII") {
        let trace = 0;
        for (let dim = 0; dim < d; dim++) trace += cov2[dim][dim];
        eigenvals[j].fill(trace / d);
        uFlats[j] = [...iFlat];
      }
    }
    if (pool) {
      if (modelType === "EEE") {
        const { values, vectors } = Matrix.fromArray(pool).eigen();
        const sharedEig = values.map((v) => Math.max(v, regCovar));
        const sharedU = vectors.toFlat();
        for (let j = 0; j < k; j++) {
          eigenvals[j] = [...sharedEig];
          uFlats[j] = [...sharedU];
        }
      } else if (modelType === "EEI") {
        const sharedDiag = new Array(d);
        for (let dim = 0; dim < d; dim++) sharedDiag[dim] = Math.max(pool[dim][dim], regCovar);
        for (let j = 0; j < k; j++) {
          eigenvals[j] = [...sharedDiag];
          uFlats[j] = [...iFlat];
        }
      } else if (modelType === "EII") {
        let trace = 0;
        for (let dim = 0; dim < d; dim++) trace += pool[dim][dim];
        const scalar = Math.max(trace / d, regCovar);
        for (let j = 0; j < k; j++) {
          eigenvals[j].fill(scalar);
          uFlats[j] = [...iFlat];
        }
      }
    }
  }
  const labels = new Array(n);
  for (let i = 0; i < n; i++) {
    let maxP = -1, best = 0;
    for (let j = 0; j < k; j++) {
      if (resp[i][j] > maxP) {
        maxP = resp[i][j];
        best = j;
      }
    }
    labels[i] = best;
  }
  const covariances = eigenvals.map((eig, j) => {
    const U = new Matrix(d, d, uFlats[j]);
    const D = Matrix.fromArray(Array.from(
      { length: d },
      (_, r) => Array.from({ length: d }, (_2, c) => r === c ? eig[r] : 0)
    ));
    return U.multiply(D).multiply(U.transpose());
  });
  const logL = prevLogL === -Infinity ? 0 : prevLogL;
  const dfMap = {
    "VVV": k * d * (d + 1) / 2,
    "EEE": d * (d + 1) / 2,
    "VVI": k * d,
    "EEI": d,
    "VII": k,
    "EII": 1
  };
  const df = k - 1 + k * d + dfMap[modelType];
  const rawEntropy = computeRawEntropy(resp, k);
  const entropy = computeNormalizedEntropy(resp, k);
  const bic = df * Math.log(n) - 2 * logL;
  const aic = 2 * df - 2 * logL;
  const icl = bic + 2 * rawEntropy;
  return {
    weights: Array.from(weights),
    means: means.map((m) => [...m]),
    covariances,
    posteriors: resp.map((r) => Array.from(r)),
    labels,
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(resp, k),
      formatted: `GMM (K = ${k}, ${modelType}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}, AvePP = [${computeAvePP(resp, k).map((v) => roundTo(v, 2)).join(", ")}]`
    }
  };
}
function predictGMM(data, result) {
  const k = result.weights.length;
  const d = result.means[0].length;
  const uFlats = result.covariances.map((cov2) => {
    const { vectors } = cov2.eigen();
    return vectors.toFlat();
  });
  const eigVals = result.covariances.map((cov2) => {
    const { values } = cov2.eigen();
    return values.map((v) => Math.max(v, 1e-12));
  });
  const n = data.length;
  const labels = new Array(n);
  const posteriors = new Array(n);
  const logLiks = new Float64Array(k);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < k; j++) {
      logLiks[j] = Math.log(Math.max(result.weights[j], MIN_PROB)) + mvnLogPdf(data[i], result.means[j], uFlats[j], eigVals[j], d);
    }
    const marg = logSumExp(logLiks);
    const post = new Array(k);
    let maxP = -1, best = 0;
    for (let j = 0; j < k; j++) {
      post[j] = Math.exp(logLiks[j] - marg);
      if (post[j] > maxP) {
        maxP = post[j];
        best = j;
      }
    }
    labels[i] = best;
    posteriors[i] = post;
  }
  return { labels, posteriors };
}
function findBestGMM(data, kRange = [1, 2, 3, 4, 5], models = ["VVV", "EEE", "VVI", "EEI", "VII", "EII"]) {
  let best = null;
  for (const k of kRange) {
    for (const model of models) {
      try {
        const res = fitGMM(data, { k, model });
        if (!best || res.diagnostics.bic < best.diagnostics.bic) best = res;
      } catch {
      }
    }
  }
  if (!best) throw new Error("findBestGMM: all model fits failed");
  return best;
}
function fitLCA(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("fitLCA: data cannot be empty");
  const m = data[0].length;
  const k = options.k;
  if (k < 1) throw new Error("fitLCA: k must be >= 1");
  const rng = new PRNG(options.seed ?? 42);
  const tol = options.tol ?? 1e-6;
  const maxIter = options.maxIter ?? 200;
  const rho = Array.from(
    { length: k },
    () => Array.from({ length: m }, () => 0.1 + rng.next() * 0.8)
  );
  const weights = new Float64Array(k).fill(1 / k);
  const resp = Array.from({ length: n }, () => new Float64Array(k));
  let prevLogL = -Infinity;
  let converged = false;
  let iter = 0;
  for (; iter < maxIter; iter++) {
    let currentLogL = 0;
    const logLiks = new Float64Array(k);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) {
        let ll = Math.log(Math.max(weights[j], MIN_PROB));
        for (let d = 0; d < m; d++) {
          const r = Math.max(Math.min(rho[j][d], 1 - 1e-12), 1e-12);
          ll += data[i][d] === 1 ? Math.log(r) : Math.log(1 - r);
        }
        logLiks[j] = ll;
      }
      const marg = logSumExp(logLiks);
      currentLogL += marg;
      for (let j = 0; j < k; j++) resp[i][j] = Math.exp(logLiks[j] - marg);
    }
    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true;
      break;
    }
    prevLogL = currentLogL;
    const Nk = new Float64Array(k);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < k; j++) Nk[j] += resp[i][j];
    }
    for (let j = 0; j < k; j++) {
      weights[j] = Nk[j] / n;
      const sumW = Nk[j];
      for (let d = 0; d < m; d++) {
        let sumX = 0;
        for (let i = 0; i < n; i++) {
          if (data[i][d] === 1) sumX += resp[i][j];
        }
        rho[j][d] = Math.max(Math.min(sumX / sumW, 1 - 1e-10), 1e-10);
      }
    }
  }
  const labels = new Array(n);
  for (let i = 0; i < n; i++) {
    let maxP = -1, best = 0;
    for (let j = 0; j < k; j++) {
      if (resp[i][j] > maxP) {
        maxP = resp[i][j];
        best = j;
      }
    }
    labels[i] = best;
  }
  const logL = prevLogL === -Infinity ? 0 : prevLogL;
  const df = k - 1 + k * m;
  const rawEntropy = computeRawEntropy(resp, k);
  const entropy = computeNormalizedEntropy(resp, k);
  const bic = df * Math.log(n) - 2 * logL;
  const aic = 2 * df - 2 * logL;
  const icl = bic + 2 * rawEntropy;
  return {
    rho: rho.map((r) => [...r]),
    priorWeights: Array.from(weights),
    posteriors: resp.map((r) => Array.from(r)),
    labels,
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(resp, k),
      formatted: `LCA (K = ${k}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}`
    }
  };
}
function logEmission(x, rho_s, m) {
  let ll = 0;
  for (let d = 0; d < m; d++) {
    const r = Math.max(Math.min(rho_s[d], 1 - 1e-12), 1e-12);
    ll += x[d] === 1 ? Math.log(r) : Math.log(1 - r);
  }
  return ll;
}
function fitLTA(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("fitLTA: data cannot be empty");
  const T = data[0].length;
  if (T < 2) throw new Error("fitLTA: requires at least 2 timepoints");
  const m = data[0][0].length;
  const k = options.k;
  if (k < 2) throw new Error("fitLTA: k must be >= 2");
  const rng = new PRNG(options.seed ?? 42);
  const tol = options.tol ?? 1e-6;
  const maxIter = options.maxIter ?? 200;
  const pi = new Float64Array(k).fill(1 / k);
  const tau = Array.from(
    { length: k },
    (_, j) => Array.from({ length: k }, (_2, l) => j === l ? 0.8 : 0.2 / (k - 1))
  );
  const rho = Array.from(
    { length: k },
    () => Array.from({ length: m }, () => 0.1 + rng.next() * 0.8)
  );
  let prevLogL = -Infinity;
  let converged = false;
  let iter = 0;
  let finalGamma = [];
  for (; iter < maxIter; iter++) {
    let currentLogL = 0;
    const piAcc = new Float64Array(k);
    const tauNum = Array.from({ length: k }, () => new Array(k).fill(0));
    const tauDen = new Float64Array(k);
    const rhoNum = Array.from({ length: k }, () => new Array(m).fill(0));
    const rhoDen = new Float64Array(k);
    const iterGamma = new Array(n);
    for (let i = 0; i < n; i++) {
      const logB = Array.from(
        { length: T },
        (_, t) => Array.from({ length: k }, (_2, s) => logEmission(data[i][t], rho[s], m))
      );
      const logAlpha = Array.from({ length: T }, () => new Array(k));
      for (let s = 0; s < k; s++) {
        logAlpha[0][s] = Math.log(Math.max(pi[s], MIN_PROB)) + logB[0][s];
      }
      for (let t = 1; t < T; t++) {
        for (let s = 0; s < k; s++) {
          const trans = new Float64Array(k);
          for (let p = 0; p < k; p++) {
            trans[p] = logAlpha[t - 1][p] + Math.log(Math.max(tau[p][s], MIN_PROB));
          }
          logAlpha[t][s] = logB[t][s] + logSumExp(trans);
        }
      }
      const subjLL = logSumExp(logAlpha[T - 1]);
      currentLogL += subjLL;
      const logBeta2 = Array.from({ length: T }, () => new Array(k));
      for (let s = 0; s < k; s++) logBeta2[T - 1][s] = 0;
      for (let t = T - 2; t >= 0; t--) {
        for (let s = 0; s < k; s++) {
          const combined = new Float64Array(k);
          for (let next = 0; next < k; next++) {
            combined[next] = Math.log(Math.max(tau[s][next], MIN_PROB)) + logB[t + 1][next] + logBeta2[t + 1][next];
          }
          logBeta2[t][s] = logSumExp(combined);
        }
      }
      const gamma2 = Array.from({ length: T }, () => new Array(k));
      for (let t = 0; t < T; t++) {
        const logRow = new Float64Array(k);
        for (let s = 0; s < k; s++) logRow[s] = logAlpha[t][s] + logBeta2[t][s];
        const den = logSumExp(logRow);
        for (let s = 0; s < k; s++) gamma2[t][s] = Math.exp(logRow[s] - den);
      }
      iterGamma[i] = gamma2;
      for (let t = 0; t < T - 1; t++) {
        const logXi = new Float64Array(k * k);
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            logXi[j * k + l] = logAlpha[t][j] + Math.log(Math.max(tau[j][l], MIN_PROB)) + logB[t + 1][l] + logBeta2[t + 1][l];
          }
        }
        const xiDen = logSumExp(logXi);
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            tauNum[j][l] += Math.exp(logXi[j * k + l] - xiDen);
          }
        }
      }
      for (let s = 0; s < k; s++) piAcc[s] += gamma2[0][s];
      for (let t = 0; t < T - 1; t++) {
        for (let s = 0; s < k; s++) tauDen[s] += gamma2[t][s];
      }
      for (let t = 0; t < T; t++) {
        for (let s = 0; s < k; s++) {
          rhoDen[s] += gamma2[t][s];
          for (let d = 0; d < m; d++) {
            if (data[i][t][d] === 1) rhoNum[s][d] += gamma2[t][s];
          }
        }
      }
    }
    finalGamma = iterGamma;
    if (Math.abs(currentLogL - prevLogL) < tol) {
      converged = true;
      break;
    }
    prevLogL = currentLogL;
    for (let s = 0; s < k; s++) {
      pi[s] = (piAcc[s] + 1) / (n + k);
    }
    for (let j = 0; j < k; j++) {
      const den = tauDen[j];
      for (let l = 0; l < k; l++) {
        tau[j][l] = (tauNum[j][l] + 0.1) / (den + 0.1 * k);
      }
    }
    for (let s = 0; s < k; s++) {
      const den = rhoDen[s];
      for (let d = 0; d < m; d++) {
        rho[s][d] = Math.max(Math.min(rhoNum[s][d] / den, 1 - 1e-10), 1e-10);
      }
    }
  }
  const trajectories = data.map((subj) => {
    const vt = Array.from({ length: T }, () => new Array(k));
    const ptr = Array.from({ length: T }, () => new Array(k));
    for (let s = 0; s < k; s++) {
      vt[0][s] = Math.log(Math.max(pi[s], MIN_PROB)) + logEmission(subj[0], rho[s], m);
    }
    for (let t = 1; t < T; t++) {
      for (let s = 0; s < k; s++) {
        const emit = logEmission(subj[t], rho[s], m);
        let maxVal = -Infinity, bestP = 0;
        for (let p = 0; p < k; p++) {
          const sc = vt[t - 1][p] + Math.log(Math.max(tau[p][s], MIN_PROB));
          if (sc > maxVal) {
            maxVal = sc;
            bestP = p;
          }
        }
        vt[t][s] = emit + maxVal;
        ptr[t][s] = bestP;
      }
    }
    const path = new Array(T);
    let maxFinal = -Infinity, bestFinal = 0;
    for (let s = 0; s < k; s++) {
      if (vt[T - 1][s] > maxFinal) {
        maxFinal = vt[T - 1][s];
        bestFinal = s;
      }
    }
    path[T - 1] = bestFinal;
    for (let t = T - 2; t >= 0; t--) {
      path[t] = ptr[t + 1][path[t + 1]];
    }
    return path;
  });
  const logL = prevLogL === -Infinity ? 0 : prevLogL;
  const df = k - 1 + k * (k - 1) + k * m;
  const flatGamma = [];
  for (let i = 0; i < n; i++) {
    for (let t = 0; t < T; t++) {
      flatGamma.push(finalGamma[i][t]);
    }
  }
  const rawEntropy = computeRawEntropy(flatGamma, k);
  const entropy = computeNormalizedEntropy(flatGamma, k);
  const bic = df * Math.log(n) - 2 * logL;
  const aic = 2 * df - 2 * logL;
  const icl = bic + 2 * rawEntropy;
  return {
    pi: Array.from(pi),
    tau: tau.map((row) => [...row]),
    rho: rho.map((r) => [...r]),
    trajectories,
    posteriors: finalGamma.map((subj) => subj.map((t) => [...t])),
    diagnostics: {
      converged,
      iterations: iter,
      logLikelihood: logL,
      df,
      aic,
      bic,
      icl,
      entropy,
      avepp: computeAvePP(flatGamma, k),
      formatted: `LTA (K = ${k}, T = ${T}): BIC = ${roundTo(bic, 1)}, AIC = ${roundTo(aic, 1)}, LL = ${roundTo(logL, 1)}`
    }
  };
}
function runKMeans(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("runKMeans: data cannot be empty");
  const d = data[0].length;
  const k = options.k;
  if (k < 1) throw new Error("runKMeans: k must be >= 1");
  if (k > n) throw new Error("runKMeans: k cannot exceed n");
  const rng = new PRNG(options.seed ?? 42);
  const maxIter = options.maxIter ?? 300;
  const tol = options.tol ?? 1e-6;
  const centroids = kMeansPlusPlus(data, k, rng);
  const labels = new Array(n).fill(0);
  let converged = false;
  let iter = 0;
  let inertia = 0;
  for (; iter < maxIter; iter++) {
    const nextCentroids = Array.from({ length: k }, () => new Array(d).fill(0));
    const counts = new Array(k).fill(0);
    inertia = 0;
    for (let i = 0; i < n; i++) {
      let minD = Infinity, bestK = 0;
      for (let j = 0; j < k; j++) {
        let dist = 0;
        for (let dim = 0; dim < d; dim++) {
          const diff = data[i][dim] - centroids[j][dim];
          dist += diff * diff;
        }
        if (dist < minD) {
          minD = dist;
          bestK = j;
        }
      }
      labels[i] = bestK;
      counts[bestK]++;
      inertia += minD;
      for (let dim = 0; dim < d; dim++) nextCentroids[bestK][dim] += data[i][dim];
    }
    let shiftSq = 0;
    for (let j = 0; j < k; j++) {
      if (counts[j] === 0) {
        let maxDist = -1, farIdx = 0;
        for (let i = 0; i < n; i++) {
          let di = 0;
          for (let dim = 0; dim < d; dim++) {
            const diff = data[i][dim] - centroids[labels[i]][dim];
            di += diff * diff;
          }
          if (di > maxDist) {
            maxDist = di;
            farIdx = i;
          }
        }
        for (let dim = 0; dim < d; dim++) {
          const updated = data[farIdx][dim];
          const diff = updated - centroids[j][dim];
          shiftSq += diff * diff;
          centroids[j][dim] = updated;
        }
      } else {
        for (let dim = 0; dim < d; dim++) {
          const updated = nextCentroids[j][dim] / counts[j];
          const diff = updated - centroids[j][dim];
          shiftSq += diff * diff;
          centroids[j][dim] = updated;
        }
      }
    }
    if (shiftSq < tol) {
      converged = true;
      break;
    }
  }
  return {
    centroids: centroids.map((c) => [...c]),
    labels,
    inertia,
    converged,
    iterations: iter
  };
}
function fitGMMRange(data, kRange, model = "VVV") {
  const entries = [];
  for (const k of kRange) {
    try {
      const result = fitGMM(data, { k, model, seed: 42 });
      entries.push({ k, model, result });
    } catch {
    }
  }
  if (entries.length === 0) throw new Error("fitGMMRange: all fits failed");
  return entries.sort((a, b) => a.k - b.k);
}
function fitKMeansRange(data, kRange) {
  const entries = [];
  for (const k of kRange) {
    try {
      const result = runKMeans(data, { k, seed: 42 });
      entries.push({ k, result });
    } catch {
    }
  }
  if (entries.length === 0) throw new Error("fitKMeansRange: all fits failed");
  return entries.sort((a, b) => a.k - b.k);
}
function predictKMeans(data, centroids) {
  const k = centroids.length;
  const d = centroids[0].length;
  return data.map((pt) => {
    let minD = Infinity, best = 0;
    for (let j = 0; j < k; j++) {
      let dist = 0;
      for (let dim = 0; dim < d; dim++) {
        const diff = pt[dim] - centroids[j][dim];
        dist += diff * diff;
      }
      if (dist < minD) {
        minD = dist;
        best = j;
      }
    }
    return best;
  });
}
function euclideanDistMatrix(data) {
  const n = data.length;
  const d = data[0]?.length ?? 0;
  const dist = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let sq = 0;
      for (let dim = 0; dim < d; dim++) {
        const diff = data[i][dim] - data[j][dim];
        sq += diff * diff;
      }
      const val = Math.sqrt(sq);
      dist[i * n + j] = val;
      dist[j * n + i] = val;
    }
  }
  return dist;
}
function silhouetteScores(data, labels) {
  const n = data.length;
  const dist = euclideanDistMatrix(data);
  const labelSet = /* @__PURE__ */ new Set();
  for (let i = 0; i < n; i++) {
    if (labels[i] >= 0) labelSet.add(labels[i]);
  }
  const uniqueLabels = [...labelSet].sort((a, b) => a - b);
  const nClusters = uniqueLabels.length;
  if (nClusters < 2) {
    return { scores: new Array(n).fill(0), mean: 0 };
  }
  const scores = new Array(n).fill(NaN);
  let sumSil = 0;
  let countNonNoise = 0;
  for (let i = 0; i < n; i++) {
    const ci = labels[i];
    if (ci < 0) continue;
    let aSum = 0, aCount = 0;
    for (let j = 0; j < n; j++) {
      if (j === i || labels[j] !== ci) continue;
      aSum += dist[i * n + j];
      aCount++;
    }
    const ai = aCount > 0 ? aSum / aCount : 0;
    let bi = Infinity;
    for (const cl of uniqueLabels) {
      if (cl === ci) continue;
      let bSum = 0, bCount = 0;
      for (let j = 0; j < n; j++) {
        if (labels[j] !== cl) continue;
        bSum += dist[i * n + j];
        bCount++;
      }
      if (bCount > 0) {
        const meanDist = bSum / bCount;
        if (meanDist < bi) bi = meanDist;
      }
    }
    const maxAB = Math.max(ai, bi);
    scores[i] = maxAB > 0 ? (bi - ai) / maxAB : 0;
    sumSil += scores[i];
    countNonNoise++;
  }
  return {
    scores,
    mean: countNonNoise > 0 ? sumSil / countNonNoise : 0
  };
}
function runDBSCAN(data, options) {
  const n = data.length;
  if (n === 0) throw new Error("runDBSCAN: data cannot be empty");
  const { eps, minPts } = options;
  if (eps <= 0) throw new Error("runDBSCAN: eps must be > 0");
  if (minPts < 1) throw new Error("runDBSCAN: minPts must be >= 1");
  const dist = euclideanDistMatrix(data);
  const neighbors = new Array(n);
  for (let i = 0; i < n; i++) {
    const nb = [];
    for (let j = 0; j < n; j++) {
      if (dist[i * n + j] <= eps) nb.push(j);
    }
    neighbors[i] = nb;
  }
  const isCore = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    if (neighbors[i].length >= minPts) isCore[i] = 1;
  }
  const labels = new Int32Array(n).fill(-1);
  let clusterId = 0;
  for (let i = 0; i < n; i++) {
    if (!isCore[i] || labels[i] !== -1) continue;
    const queue = [i];
    labels[i] = clusterId;
    let head = 0;
    while (head < queue.length) {
      const current = queue[head];
      head++;
      for (const nb of neighbors[current]) {
        if (labels[nb] === -1) {
          labels[nb] = clusterId;
          if (isCore[nb]) queue.push(nb);
        }
      }
    }
    clusterId++;
  }
  const pointTypes = new Array(n);
  for (let i = 0; i < n; i++) {
    if (isCore[i]) {
      pointTypes[i] = "core";
    } else if (labels[i] >= 0) {
      pointTypes[i] = "border";
    } else {
      pointTypes[i] = "noise";
    }
  }
  const nClusters = clusterId;
  const clusterSizes = new Array(nClusters).fill(0);
  let nNoise = 0;
  for (let i = 0; i < n; i++) {
    if (labels[i] >= 0) clusterSizes[labels[i]]++;
    else nNoise++;
  }
  const sil = nClusters >= 2 ? silhouetteScores(data, Array.from(labels)) : { scores: new Array(n).fill(NaN), mean: NaN };
  const formatted = `DBSCAN (eps = ${roundTo(eps, 3)}, minPts = ${minPts}): ${nClusters} clusters, ${nNoise} noise points, silhouette = ${isNaN(sil.mean) ? "N/A" : roundTo(sil.mean, 3)}`;
  return {
    labels: Array.from(labels),
    pointTypes,
    nClusters,
    nNoise,
    clusterSizes,
    silhouette: sil,
    formatted
  };
}
function kDistancePlot(data, k) {
  const n = data.length;
  if (k < 1 || k >= n) throw new Error(`kDistancePlot: k must be in [1, n-1], got ${k}`);
  const dist = euclideanDistMatrix(data);
  const kDists = new Array(n);
  for (let i = 0; i < n; i++) {
    const dists = new Array(n - 1);
    let idx = 0;
    for (let j = 0; j < n; j++) {
      if (j !== i) dists[idx++] = dist[i * n + j];
    }
    dists.sort((a, b) => a - b);
    kDists[i] = dists[k - 1];
  }
  return kDists.sort((a, b) => a - b);
}
function runHierarchical(data, options) {
  const n = data.length;
  if (n < 2) throw new Error("runHierarchical: need at least 2 observations");
  const linkage = options?.linkage ?? "ward";
  const origDist = euclideanDistMatrix(data);
  const useSquared = linkage === "ward";
  const D = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const val = useSquared ? origDist[i * n + j] * origDist[i * n + j] : origDist[i * n + j];
      D[i * n + j] = val;
      D[j * n + i] = val;
    }
  }
  const sizes = new Float64Array(2 * n - 1);
  for (let i = 0; i < n; i++) sizes[i] = 1;
  const active = /* @__PURE__ */ new Set();
  for (let i = 0; i < n; i++) active.add(i);
  const merges = [];
  const mergeHeights = [];
  const distMap = /* @__PURE__ */ new Map();
  function getDist(a, b) {
    if (a < n && b < n) return D[a * n + b];
    const key = a < b ? `${a},${b}` : `${b},${a}`;
    return distMap.get(key) ?? Infinity;
  }
  function setDist(a, b, val) {
    if (a < n && b < n) {
      D[a * n + b] = val;
      D[b * n + a] = val;
    } else {
      const key = a < b ? `${a},${b}` : `${b},${a}`;
      distMap.set(key, val);
    }
  }
  let nextCluster = n;
  for (let step = 0; step < n - 1; step++) {
    let minDist = Infinity;
    let mergeA = -1, mergeB = -1;
    const activeArr = [...active];
    for (let ii = 0; ii < activeArr.length; ii++) {
      const ci = activeArr[ii];
      for (let jj = ii + 1; jj < activeArr.length; jj++) {
        const cj = activeArr[jj];
        const d = getDist(ci, cj);
        if (d < minDist) {
          minDist = d;
          mergeA = ci;
          mergeB = cj;
        }
      }
    }
    const height = useSquared ? Math.sqrt(minDist) : minDist;
    merges.push({ a: mergeA, b: mergeB, height });
    mergeHeights.push(height);
    const newCluster = nextCluster++;
    sizes[newCluster] = sizes[mergeA] + sizes[mergeB];
    active.delete(mergeA);
    active.delete(mergeB);
    const ni = sizes[mergeA];
    const nj = sizes[mergeB];
    for (const ck of active) {
      const nk = sizes[ck];
      const dik = getDist(mergeA, ck);
      const djk = getDist(mergeB, ck);
      const dij = getDist(mergeA, mergeB);
      let newDist;
      if (linkage === "single") {
        newDist = 0.5 * dik + 0.5 * djk - 0.5 * Math.abs(dik - djk);
      } else if (linkage === "complete") {
        newDist = 0.5 * dik + 0.5 * djk + 0.5 * Math.abs(dik - djk);
      } else if (linkage === "average") {
        const ai = ni / (ni + nj);
        const aj = nj / (ni + nj);
        newDist = ai * dik + aj * djk;
      } else {
        const nt = ni + nj + nk;
        newDist = (ni + nk) / nt * dik + (nj + nk) / nt * djk - nk / nt * dij;
      }
      setDist(newCluster, ck, newDist);
    }
    active.add(newCluster);
  }
  const order = buildLeafOrder(merges, n);
  const cophCorr = computeCopheneticCorrelation(merges, n, origDist);
  const formatted = `HAC (${linkage}): ${n} observations, cophenetic r = ${roundTo(cophCorr, 3)}`;
  return {
    merges,
    heights: mergeHeights,
    order,
    copheneticCorrelation: cophCorr,
    formatted
  };
}
function cutTree(result, k) {
  const nMerges = result.merges.length;
  const n = nMerges + 1;
  if (k < 1 || k > n) throw new Error(`cutTree: k must be in [1, ${n}], got ${k}`);
  const parent = new Int32Array(2 * n - 1);
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i;
  function find(x) {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }
  const nMergesToApply = n - k;
  for (let s = 0; s < nMergesToApply; s++) {
    const merge = result.merges[s];
    const newCluster = n + s;
    parent[find(merge.a)] = newCluster;
    parent[find(merge.b)] = newCluster;
  }
  const rootToLabel = /* @__PURE__ */ new Map();
  let nextLabel = 0;
  const labels = new Array(n);
  for (let i = 0; i < n; i++) {
    const root = find(i);
    let label = rootToLabel.get(root);
    if (label === void 0) {
      label = nextLabel++;
      rootToLabel.set(root, label);
    }
    labels[i] = label;
  }
  return labels;
}
function cutTreeHeight(result, h) {
  const nMerges = result.merges.length;
  const n = nMerges + 1;
  const parent = new Int32Array(2 * n - 1);
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i;
  function find(x) {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }
  for (let s = 0; s < nMerges; s++) {
    const merge = result.merges[s];
    if (merge.height > h) break;
    const newCluster = n + s;
    parent[find(merge.a)] = newCluster;
    parent[find(merge.b)] = newCluster;
  }
  const rootToLabel = /* @__PURE__ */ new Map();
  let nextLabel = 0;
  const labels = new Array(n);
  for (let i = 0; i < n; i++) {
    const root = find(i);
    let label = rootToLabel.get(root);
    if (label === void 0) {
      label = nextLabel++;
      rootToLabel.set(root, label);
    }
    labels[i] = label;
  }
  return labels;
}
function buildLeafOrder(merges, n) {
  const children = /* @__PURE__ */ new Map();
  for (let s = 0; s < merges.length; s++) {
    const merge = merges[s];
    const nodeId = n + s;
    children.set(nodeId, [merge.a, merge.b]);
  }
  const root = n + merges.length - 1;
  const order = [];
  function dfs(node) {
    const ch = children.get(node);
    if (ch) {
      dfs(ch[0]);
      dfs(ch[1]);
    } else {
      order.push(node);
    }
  }
  dfs(root);
  return order;
}
function computeCopheneticCorrelation(merges, n, origDist) {
  const parent = new Int32Array(2 * n - 1);
  for (let i = 0; i < 2 * n - 1; i++) parent[i] = i;
  function find(x) {
    while (parent[x] !== x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }
  const members = /* @__PURE__ */ new Map();
  for (let i = 0; i < n; i++) members.set(i, [i]);
  const nPairs = n * (n - 1) / 2;
  const cophDist = new Float64Array(nPairs);
  const origFlat = new Float64Array(nPairs);
  let idx = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      origFlat[idx] = origDist[i * n + j];
      idx++;
    }
  }
  for (let s = 0; s < merges.length; s++) {
    const merge = merges[s];
    const newCluster = n + s;
    const rootA = find(merge.a);
    const rootB = find(merge.b);
    const membersA = members.get(rootA) ?? [];
    const membersB = members.get(rootB) ?? [];
    for (const a of membersA) {
      for (const b of membersB) {
        const i = Math.min(a, b);
        const j = Math.max(a, b);
        const flatIdx = i * n - i * (i + 1) / 2 + (j - i - 1);
        cophDist[flatIdx] = merge.height;
      }
    }
    parent[rootA] = newCluster;
    parent[rootB] = newCluster;
    members.set(newCluster, [...membersA, ...membersB]);
    members.delete(rootA);
    members.delete(rootB);
  }
  return pearsonR(origFlat, cophDist, nPairs);
}
function pearsonR(x, y, n) {
  let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sx += x[i];
    sy += y[i];
    sxx += x[i] * x[i];
    syy += y[i] * y[i];
    sxy += x[i] * y[i];
  }
  const num = n * sxy - sx * sy;
  const den = Math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
  return den > 0 ? num / den : 0;
}

// src/stats/factor-analysis.ts
var PRNG2 = class {
  state;
  constructor(seed) {
    this.state = seed >>> 0;
  }
  next() {
    this.state = this.state + 2654435769 | 0;
    let t = this.state ^ this.state >>> 16;
    t = Math.imul(t, 569420461);
    t = t ^ t >>> 15;
    t = Math.imul(t, 1935289751);
    t = t ^ t >>> 15;
    return (t >>> 0) / 4294967296;
  }
};
function prngNormal(rng) {
  let u = 0, v = 0;
  while (u === 0) u = rng.next();
  while (v === 0) v = rng.next();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function randomOrthogonalMatrix(k, rng) {
  const M = Array.from(
    { length: k },
    () => Array.from({ length: k }, () => prngNormal(rng))
  );
  const Q = M.map((r) => [...r]);
  const Rdiag = new Float64Array(k);
  for (let j = 0; j < k; j++) {
    for (let prev = 0; prev < j; prev++) {
      let dot = 0;
      for (let i = 0; i < k; i++) dot += Q[i][j] * Q[i][prev];
      for (let i = 0; i < k; i++) Q[i][j] = Q[i][j] - dot * Q[i][prev];
    }
    let norm = 0;
    for (let i = 0; i < k; i++) norm += Q[i][j] ** 2;
    norm = Math.sqrt(norm);
    Rdiag[j] = norm;
    if (norm > 1e-15) {
      for (let i = 0; i < k; i++) Q[i][j] = Q[i][j] / norm;
    }
  }
  for (let j = 0; j < k; j++) {
    if (Rdiag[j] < 0) {
      for (let i = 0; i < k; i++) Q[i][j] = -Q[i][j];
    }
  }
  return Q;
}
function computeCorrelationMatrix(data, n, d) {
  const means = new Float64Array(d);
  const sds = new Float64Array(d);
  for (let i = 0; i < n; i++) {
    const row = data[i];
    for (let j = 0; j < d; j++) means[j] = means[j] + row[j] / n;
  }
  for (let i = 0; i < n; i++) {
    const row = data[i];
    for (let j = 0; j < d; j++) sds[j] = sds[j] + (row[j] - means[j]) ** 2;
  }
  for (let j = 0; j < d; j++) sds[j] = Math.sqrt(sds[j] / (n - 1));
  const R = Array.from({ length: d }, () => new Array(d).fill(0));
  for (let r = 0; r < d; r++) {
    R[r][r] = 1;
    for (let c = r + 1; c < d; c++) {
      let sum = 0;
      const sdR = sds[r] || 1;
      const sdC = sds[c] || 1;
      for (let i = 0; i < n; i++) {
        sum += (data[i][r] - means[r]) / sdR * ((data[i][c] - means[c]) / sdC);
      }
      const val = sum / (n - 1);
      R[r][c] = val;
      R[c][r] = val;
    }
  }
  return Matrix.fromArray(R);
}
function computeFit(S, Sigma, n, d, nFreeParams, nFactors) {
  let logDetSigma, logDetS;
  try {
    logDetSigma = Sigma.logDet();
  } catch {
    const ev = Sigma.eigen().values;
    logDetSigma = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0);
  }
  try {
    logDetS = S.logDet();
  } catch {
    const ev = S.eigen().values;
    logDetS = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0);
  }
  let traceVal;
  try {
    traceVal = Sigma.inverse().multiply(S).trace();
  } catch {
    traceVal = Sigma.pseudoInverse().multiply(S).trace();
  }
  const F_ml = Math.max(0, logDetSigma + traceVal - logDetS - d);
  const bartlettN = nFactors !== void 0 ? n - 1 - (2 * d + 5) / 6 - 2 * nFactors / 3 : n - 1;
  const chiSq = Math.max(0, bartlettN * F_ml);
  const totalElements = d * (d + 1) / 2;
  const df = Math.max(0, totalElements - nFreeParams);
  const pValue = df > 0 ? chiSqPValue(chiSq, df) : 1;
  let logDetDiagS = 0;
  let traceNull = 0;
  for (let i = 0; i < d; i++) {
    const diagVal = Math.max(S.get(i, i), 1e-15);
    logDetDiagS += Math.log(diagVal);
    traceNull += S.get(i, i) / diagVal;
  }
  const diagArr = Array.from(
    { length: d },
    (_, i) => Array.from({ length: d }, (_2, j) => i === j ? S.get(i, i) : 0)
  );
  const diagMat = Matrix.fromArray(diagArr);
  let traceNullFull;
  try {
    traceNullFull = diagMat.inverse().multiply(S).trace();
  } catch {
    traceNullFull = d;
  }
  const F_null = Math.max(0, logDetDiagS + traceNullFull - logDetS - d);
  const bartlettNNull = nFactors !== void 0 ? n - 1 - (2 * d + 5) / 6 : n - 1;
  const chiSqNull = bartlettNNull * F_null;
  const dfNull = d * (d - 1) / 2;
  const ncp = Math.max(chiSq - df, 0);
  const rmsea = df > 0 ? Math.sqrt(ncp / (df * (n - 1))) : 0;
  let rmseaLo = 0, rmseaHi = rmsea * 2;
  if (df > 0) {
    const ncpLo = Math.max(chiSq - df - 1.645 * Math.sqrt(2 * df), 0);
    const ncpHi = Math.max(chiSq - df + 1.645 * Math.sqrt(2 * df), 0);
    rmseaLo = Math.sqrt(Math.max(ncpLo / (df * (n - 1)), 0));
    rmseaHi = Math.sqrt(ncpHi / (df * (n - 1)));
  }
  const ncpNull = Math.max(chiSqNull - dfNull, 0);
  const cfi = ncpNull > 0 ? Math.max(0, Math.min(1, 1 - ncp / ncpNull)) : 1;
  const tli = df > 0 && dfNull > 0 ? (chiSqNull / dfNull - chiSq / df) / (chiSqNull / dfNull - 1) : 1;
  let srmrSum = 0;
  let srmrCount = 0;
  for (let i = 0; i < d; i++) {
    for (let j = 0; j <= i; j++) {
      const sij = S.get(i, j);
      const sigij = Sigma.get(i, j);
      const denom = Math.sqrt(S.get(i, i) * S.get(j, j));
      const r_obs = denom > 0 ? sij / denom : 0;
      const r_imp = denom > 0 ? sigij / denom : 0;
      srmrSum += (r_obs - r_imp) ** 2;
      srmrCount++;
    }
  }
  const srmr = Math.sqrt(srmrSum / srmrCount);
  const aic = chiSq + 2 * nFreeParams;
  const bic = chiSq + nFreeParams * Math.log(n);
  return {
    chiSq,
    df,
    pValue,
    rmsea,
    rmseaCI: [rmseaLo, rmseaHi],
    cfi,
    tli,
    srmr,
    aic,
    bic
  };
}
function extractPAF(R, k, maxIter, tol) {
  const d = R.rows;
  const h2 = new Float64Array(d);
  try {
    const invR = R.inverse();
    for (let i = 0; i < d; i++) h2[i] = Math.max(0.01, 1 - 1 / Math.max(invR.get(i, i), 1e-12));
  } catch {
    for (let i = 0; i < d; i++) h2[i] = 0.5;
  }
  const loadings = Array.from({ length: d }, () => new Array(k).fill(0));
  for (let iter = 0; iter < maxIter; iter++) {
    const adjR = R.toArray();
    for (let i = 0; i < d; i++) adjR[i][i] = h2[i];
    const adjM = Matrix.fromArray(adjR);
    const { values, vectors } = adjM.eigen();
    const oldH2 = new Float64Array(h2);
    for (let f = 0; f < k; f++) {
      const eigenVal = Math.max(values[f], 0);
      const scale = Math.sqrt(eigenVal);
      for (let r = 0; r < d; r++) {
        loadings[r][f] = vectors.get(r, f) * scale;
      }
    }
    let maxDelta = 0;
    for (let r = 0; r < d; r++) {
      let sum = 0;
      for (let f = 0; f < k; f++) sum += loadings[r][f] ** 2;
      h2[r] = Math.max(1e-3, Math.min(0.9999, sum));
      maxDelta = Math.max(maxDelta, Math.abs(h2[r] - oldH2[r]));
    }
    if (maxDelta < tol) break;
  }
  for (let f = 0; f < k; f++) {
    let maxAbs = 0, maxIdx = 0;
    for (let i = 0; i < d; i++) {
      const a = Math.abs(loadings[i][f]);
      if (a > maxAbs) {
        maxAbs = a;
        maxIdx = i;
      }
    }
    if (loadings[maxIdx][f] < 0) {
      for (let i = 0; i < d; i++) loadings[i][f] = -loadings[i][f];
    }
  }
  return { loadings, communalities: h2 };
}
function extractML(R, k, maxIter, tol) {
  const d = R.rows;
  const Rarr = R.toArray();
  const Theta = new Float64Array(d);
  try {
    const invR = R.inverse();
    const factor = 1 - 0.5 * k / d;
    for (let i = 0; i < d; i++) Theta[i] = Math.max(5e-3, Math.min(0.995, factor / Math.max(invR.get(i, i), 1e-12)));
  } catch {
    for (let i = 0; i < d; i++) Theta[i] = 0.5;
  }
  function extractLoadingsFromTheta(theta) {
    const scaledR = Array.from(
      { length: d },
      (_, i) => Array.from({ length: d }, (_2, j) => {
        const si = 1 / Math.sqrt(Math.max(theta[i], 1e-12));
        const sj = 1 / Math.sqrt(Math.max(theta[j], 1e-12));
        return Rarr[i][j] * si * sj;
      })
    );
    const { values, vectors } = Matrix.fromArray(scaledR).eigen();
    const L2 = Array.from({ length: d }, () => new Array(k).fill(0));
    for (let f = 0; f < k; f++) {
      const ev = Math.max(values[f] - 1, 0);
      const scale = Math.sqrt(ev);
      for (let i = 0; i < d; i++) {
        L2[i][f] = Math.sqrt(Math.max(theta[i], 1e-12)) * vectors.get(i, f) * scale;
      }
    }
    return L2;
  }
  const lr0 = 0.02;
  const lrMin = 1e-3;
  let L = extractLoadingsFromTheta(Theta);
  const vTheta = new Float64Array(d);
  const momentum = 0.8;
  for (let iter = 0; iter < maxIter; iter++) {
    const lr = lrMin + (lr0 - lrMin) * 0.5 * (1 + Math.cos(Math.PI * iter / maxIter));
    L = extractLoadingsFromTheta(Theta);
    const sigmaArr = Array.from(
      { length: d },
      (_, i) => Array.from({ length: d }, (_2, j) => {
        let sum = 0;
        for (let f = 0; f < k; f++) sum += L[i][f] * L[j][f];
        return sum + (i === j ? Theta[i] : 0);
      })
    );
    const Sigma = Matrix.fromArray(sigmaArr);
    let invSigma;
    try {
      invSigma = Sigma.inverse();
    } catch {
      invSigma = Sigma.pseudoInverse();
    }
    const Delta = invSigma.multiply(Sigma.subtract(R)).multiply(invSigma);
    let maxGrad = 0;
    for (let i = 0; i < d; i++) {
      const grad = Delta.get(i, i);
      maxGrad = Math.max(maxGrad, Math.abs(grad));
      const v = momentum * vTheta[i] - lr * grad;
      vTheta[i] = v;
      Theta[i] = Math.max(5e-3, Math.min(0.995, Theta[i] + v));
    }
    if (maxGrad < tol) break;
  }
  function concentratedML(x) {
    const scaledR = Array.from(
      { length: d },
      (_, i) => Array.from({ length: d }, (_2, j) => {
        const psi_i = Math.max(x[i], 5e-3);
        const psi_j = Math.max(x[j], 5e-3);
        return Rarr[i][j] / (Math.sqrt(psi_i) * Math.sqrt(psi_j));
      })
    );
    const { values } = Matrix.fromArray(scaledR).eigen();
    let sum = 0;
    for (let j = k; j < d; j++) {
      const ev = Math.max(values[j], 1e-15);
      sum += Math.log(ev) - ev;
    }
    let penalty = 0;
    for (let i = 0; i < d; i++) {
      if (x[i] < 5e-3 || x[i] > 0.995) penalty += 1e3;
    }
    return -sum + k - d + penalty;
  }
  const psi0 = Array.from(Theta);
  const nmResult = nelderMead(concentratedML, psi0, {
    maxIter: 5e3 * d,
    tol: 1e-10
  });
  const finalTheta = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    finalTheta[i] = Math.max(5e-3, Math.min(0.995, nmResult.x[i]));
  }
  L = extractLoadingsFromTheta(finalTheta);
  for (let f = 0; f < k; f++) {
    let maxAbs = 0, maxIdx = 0;
    for (let i = 0; i < d; i++) {
      const a = Math.abs(L[i][f]);
      if (a > maxAbs) {
        maxAbs = a;
        maxIdx = i;
      }
    }
    if (L[maxIdx][f] < 0) {
      for (let i = 0; i < d; i++) L[i][f] = -L[i][f];
    }
  }
  const h2 = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    let sum = 0;
    for (let f = 0; f < k; f++) sum += L[i][f] ** 2;
    h2[i] = Math.max(1e-3, Math.min(0.9999, sum));
  }
  return { loadings: L, communalities: h2 };
}
function rotateVarimax(L, maxIter, tol) {
  const p = L.length;
  const k = L[0].length;
  if (k < 2) return { rotated: L.map((r) => [...r]), T: [[1]] };
  const sc = new Float64Array(p);
  for (let i = 0; i < p; i++) {
    let ss2 = 0;
    for (let j = 0; j < k; j++) ss2 += L[i][j] ** 2;
    sc[i] = Math.sqrt(ss2 || 1e-15);
  }
  const x = L.map((row, i) => row.map((v) => v / sc[i]));
  let T = Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
  );
  let dPast = 0;
  for (let iter = 0; iter < Math.min(maxIter, 1e3); iter++) {
    const z = Array.from(
      { length: p },
      (_, i) => Array.from({ length: k }, (_2, j) => {
        let s = 0;
        for (let m = 0; m < k; m++) s += x[i][m] * T[m][j];
        return s;
      })
    );
    const csz2 = new Float64Array(k);
    for (let j = 0; j < k; j++) {
      let s = 0;
      for (let i = 0; i < p; i++) s += z[i][j] ** 2;
      csz2[j] = s / p;
    }
    const B = Array.from(
      { length: k },
      (_, r) => Array.from({ length: k }, (_2, c) => {
        let s = 0;
        for (let i = 0; i < p; i++) {
          const zij = z[i][c];
          s += x[i][r] * (zij * zij * zij - zij * csz2[c]);
        }
        return s;
      })
    );
    const BtB = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => {
        let s = 0;
        for (let m = 0; m < k; m++) s += B[m][i] * B[m][j];
        return s;
      })
    );
    const { values: sigma2, vectors: Vmat } = Matrix.fromArray(BtB).eigen();
    const svals = new Float64Array(k);
    for (let j = 0; j < k; j++) svals[j] = Math.sqrt(Math.max(sigma2[j], 0));
    const Varr = Vmat.toArray();
    const BV = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => {
        let s = 0;
        for (let m = 0; m < k; m++) s += B[i][m] * Varr[m][j];
        return s;
      })
    );
    for (let j = 0; j < k; j++) {
      const invS = svals[j] > 1e-15 ? 1 / svals[j] : 0;
      for (let i = 0; i < k; i++) BV[i][j] = BV[i][j] * invS;
    }
    T = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => {
        let s = 0;
        for (let m = 0; m < k; m++) s += BV[i][m] * Varr[j][m];
        return s;
      })
    );
    let dNew = 0;
    for (let j = 0; j < k; j++) dNew += svals[j];
    if (dNew < dPast * (1 + tol)) break;
    dPast = dNew;
  }
  const rot = Array.from(
    { length: p },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      let s = 0;
      for (let m = 0; m < k; m++) s += x[i][m] * T[m][j];
      return s * sc[i];
    })
  );
  return { rotated: rot, T };
}
function criterionOblimin(L, gamma2) {
  const p = L.length, k = L[0].length;
  const Gq = Array.from({ length: p }, () => new Array(k).fill(0));
  let f = 0;
  for (let i = 0; i < p; i++) {
    let rowSumSq = 0;
    for (let m = 0; m < k; m++) rowSumSq += L[i][m] ** 2;
    for (let j = 0; j < k; j++) {
      const lij = L[i][j];
      const otherSumSq = rowSumSq - lij * lij;
      Gq[i][j] = lij * (otherSumSq - gamma2 / p * rowSumSq);
      f += lij * lij * otherSumSq;
    }
  }
  f /= 4;
  return { f, Gq };
}
function criterionGeomin(L, delta) {
  const p = L.length, k = L[0].length;
  const Gq = Array.from({ length: p }, () => new Array(k).fill(0));
  let f = 0;
  for (let i = 0; i < p; i++) {
    let logSum = 0;
    for (let j = 0; j < k; j++) logSum += Math.log(L[i][j] ** 2 + delta);
    const pro = Math.exp(logSum / k);
    f += pro;
    for (let j = 0; j < k; j++) {
      Gq[i][j] = 2 / k * L[i][j] / (L[i][j] ** 2 + delta) * pro;
    }
  }
  return { f, Gq };
}
function gpfOblq(A, criterion, maxIter, tol, Tinit) {
  const k = A[0].length;
  const Amat = Matrix.fromArray(A);
  let Tarr = Tinit ? Tinit.map((r) => [...r]) : Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
  );
  let Tmat = Matrix.fromArray(Tarr);
  let Larr;
  if (Tinit) {
    let invT;
    try {
      invT = Tmat.inverse();
    } catch {
      invT = Tmat.pseudoInverse();
    }
    Larr = Amat.multiply(invT.transpose()).toArray();
  } else {
    Larr = A.map((r) => [...r]);
  }
  let vgQ = criterion(Larr);
  let f = vgQ.f;
  let G = computeGPFoblqGradient(Larr, vgQ.Gq, Tmat);
  let al = 1;
  for (let iter = 0; iter < maxIter; iter++) {
    const colSumTG = new Float64Array(k);
    for (let j = 0; j < k; j++) {
      let s3 = 0;
      for (let i = 0; i < k; i++) s3 += Tarr[i][j] * G[i][j];
      colSumTG[j] = s3;
    }
    const Gp = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => G[i][j] - Tarr[i][j] * colSumTG[j])
    );
    let s2 = 0;
    for (let i = 0; i < k; i++)
      for (let j = 0; j < k; j++) s2 += Gp[i][j] ** 2;
    const s = Math.sqrt(s2);
    if (s < tol) break;
    al = 2 * al;
    let Tmatt = Tarr;
    let TmattMat = Tmat;
    let Lnew = Larr;
    let vgQnew = vgQ;
    for (let lsIter = 0; lsIter <= 10; lsIter++) {
      const X = Array.from(
        { length: k },
        (_, i) => Array.from({ length: k }, (_2, j) => Tarr[i][j] - al * Gp[i][j])
      );
      for (let j = 0; j < k; j++) {
        let ss2 = 0;
        for (let i = 0; i < k; i++) ss2 += X[i][j] ** 2;
        const invNorm = 1 / Math.sqrt(ss2 || 1e-15);
        for (let i = 0; i < k; i++) X[i][j] = X[i][j] * invNorm;
      }
      const Ttmat = Matrix.fromArray(X);
      let invTt;
      try {
        invTt = Ttmat.inverse();
      } catch {
        invTt = Ttmat.pseudoInverse();
      }
      Lnew = Amat.multiply(invTt.transpose()).toArray();
      vgQnew = criterion(Lnew);
      Tmatt = X;
      TmattMat = Ttmat;
      const improvement = f - vgQnew.f;
      if (improvement > 0.5 * s2 * al) break;
      al = al / 2;
    }
    Tarr = Tmatt;
    Tmat = TmattMat;
    Larr = Lnew;
    f = vgQnew.f;
    G = computeGPFoblqGradient(Lnew, vgQnew.Gq, TmattMat);
  }
  const Phi = Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      let s = 0;
      for (let l = 0; l < k; l++) s += Tarr[l][i] * Tarr[l][j];
      return s;
    })
  );
  return { rotated: Larr, T: Tarr, Phi, f };
}
function computeGPFoblqGradient(L, Gq, Tmat) {
  const Lmat = Matrix.fromArray(L);
  const GqMat = Matrix.fromArray(Gq);
  let invT;
  try {
    invT = Tmat.inverse();
  } catch {
    invT = Tmat.pseudoInverse();
  }
  const inner = Lmat.transpose().multiply(GqMat).multiply(invT);
  const G = Array.from(
    { length: inner.rows },
    (_, i) => Array.from({ length: inner.cols }, (_2, j) => -inner.get(j, i))
  );
  return G;
}
function rotatePromax(L, power, maxIter, tol) {
  const d = L.length;
  const k = L[0].length;
  const h2 = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    let sum = 0;
    for (let j = 0; j < k; j++) sum += L[i][j] * L[i][j];
    h2[i] = sum;
  }
  const sqrtH2 = new Float64Array(d);
  for (let i = 0; i < d; i++) sqrtH2[i] = Math.sqrt(Math.max(h2[i], 1e-15));
  const L_norm = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => L[i][j] / sqrtH2[i])
  );
  const { rotated: vari, T: T_varimax } = rotateVarimax(L_norm, maxIter, tol);
  const Q = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      const val = vari[i][j];
      return Math.sign(val) * Math.abs(val) ** power;
    })
  );
  const Vmat = Matrix.fromArray(vari);
  const Qmat = Matrix.fromArray(Q);
  let VtV_inv;
  try {
    VtV_inv = Vmat.transpose().multiply(Vmat).inverse();
  } catch {
    VtV_inv = Vmat.transpose().multiply(Vmat).pseudoInverse();
  }
  const U = VtV_inv.multiply(Vmat.transpose()).multiply(Qmat);
  const Uarr = U.toArray();
  const UtU = U.transpose().multiply(U);
  let UtU_inv;
  try {
    UtU_inv = UtU.inverse();
  } catch {
    UtU_inv = UtU.pseudoInverse();
  }
  for (let j = 0; j < k; j++) {
    const d_j = Math.max(UtU_inv.get(j, j), 1e-15);
    const scale = Math.sqrt(d_j);
    for (let i = 0; i < k; i++) Uarr[i][j] = Uarr[i][j] * scale;
  }
  const Umat = Matrix.fromArray(Uarr);
  const rotatedNorm = Vmat.multiply(Umat).toArray();
  const rotated = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => rotatedNorm[i][j] * sqrtH2[i])
  );
  const TvarMat = Matrix.fromArray(T_varimax);
  const Tcompound = TvarMat.multiply(Umat);
  const Tarr = Tcompound.toArray();
  let invT;
  try {
    invT = Tcompound.inverse();
  } catch {
    invT = Tcompound.pseudoInverse();
  }
  const PhiMat = invT.multiply(invT.transpose());
  const Phi = PhiMat.toArray();
  return { rotated, T: Tarr, Phi };
}
function applyRotation(loadings, method, maxIter, tol, geominDelta = 0.01, randomStarts = 1, seed = 42) {
  const k = loadings[0].length;
  if (method === "none") {
    const Phi = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
    );
    return { rotated: loadings.map((r) => [...r]), Phi };
  }
  if (method === "varimax") {
    const { rotated } = rotateVarimax(loadings, maxIter, tol);
    const Phi = Array.from(
      { length: k },
      (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
    );
    return { rotated, Phi };
  }
  if (method === "oblimin" || method === "quartimin") {
    const gamma2 = method === "quartimin" ? 0 : 0;
    const criterionFn = (L) => criterionOblimin(L, gamma2);
    return gpfOblqWithRandomStarts(loadings, criterionFn, maxIter, tol, k, randomStarts, seed);
  }
  if (method === "geomin") {
    const criterionFn = (L) => criterionGeomin(L, geominDelta);
    return gpfOblqWithRandomStarts(loadings, criterionFn, maxIter, tol, k, randomStarts, seed);
  }
  if (method === "promax") {
    const { rotated, Phi } = rotatePromax(loadings, 4, maxIter, tol);
    return { rotated, Phi };
  }
  throw new Error(`runEFA: unknown rotation method '${method}'`);
}
function gpfOblqWithRandomStarts(loadings, criterionFn, maxIter, tol, k, randomStarts, seed) {
  let best = gpfOblq(loadings, criterionFn, maxIter, tol);
  if (randomStarts > 1) {
    const rng = new PRNG2(seed);
    for (let s = 1; s < randomStarts; s++) {
      const Trand = randomOrthogonalMatrix(k, rng);
      const result = gpfOblq(loadings, criterionFn, maxIter, tol, Trand);
      if (result.f < best.f) {
        best = result;
      }
    }
  }
  return { rotated: best.rotated, Phi: best.Phi };
}
function velicerMAP(R) {
  const d = R.rows;
  const { values, vectors } = R.eigen();
  let bestK = 0;
  let minMap = Infinity;
  for (let k = 0; k < d - 1; k++) {
    const Larr = Array.from(
      { length: d },
      (_, i) => Array.from(
        { length: k + 1 },
        (_2, j) => vectors.get(i, j) * Math.sqrt(Math.max(values[j], 0))
      )
    );
    const Lmat = Matrix.fromArray(Larr);
    const R_star = R.subtract(Lmat.multiply(Lmat.transpose()));
    const C_star = Array.from(
      { length: d },
      (_, i) => Array.from({ length: d }, (_2, j) => {
        const denom = Math.sqrt(Math.abs(R_star.get(i, i)) * Math.abs(R_star.get(j, j)));
        return denom > 1e-15 ? R_star.get(i, j) / denom : i === j ? 1 : 0;
      })
    );
    let sumSq = 0;
    for (let r = 0; r < d; r++) {
      for (let c = 0; c < r; c++) sumSq += C_star[r][c] ** 2;
    }
    const mapVal = sumSq / (d * (d - 1) / 2);
    if (mapVal < minMap) {
      minMap = mapVal;
      bestK = k + 1;
    }
  }
  return bestK;
}
function parallelAnalysis(observedEigenvalues, n, d, iterations, rng) {
  const allEigens = Array.from({ length: d }, () => new Array(iterations).fill(0));
  for (let iter = 0; iter < iterations; iter++) {
    const randomData = Array.from(
      { length: n },
      () => Array.from({ length: d }, () => prngNormal(rng))
    );
    const randR = computeCorrelationMatrix(randomData, n, d);
    const randEig = randR.eigen().values;
    for (let i = 0; i < d; i++) {
      allEigens[i][iter] = randEig[i];
    }
  }
  const simulated = allEigens.map((eigArray) => {
    const sorted = [...eigArray].sort((a, b) => a - b);
    const idx = Math.floor(0.95 * iterations);
    return sorted[Math.min(idx, iterations - 1)];
  });
  let suggested = 0;
  for (let i = 0; i < d; i++) {
    if (observedEigenvalues[i] > simulated[i]) suggested++;
    else break;
  }
  return { simulated, suggested: Math.max(1, suggested) };
}
function computeKMOBartlett(R, n, d) {
  let invR;
  try {
    invR = R.inverse();
  } catch {
    invR = R.pseudoInverse();
  }
  const S2diag = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    S2diag[i] = 1 / Math.max(invR.get(i, i), 1e-12);
  }
  const antiImage = Array.from(
    { length: d },
    (_, i) => Array.from({ length: d }, (_2, j) => {
      if (i === j) return 1;
      return -invR.get(i, j) / Math.sqrt(Math.max(invR.get(i, i) * invR.get(j, j), 1e-12));
    })
  );
  let rSum = 0, qSum = 0;
  const msaItems = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    let rIdx = 0, qIdx = 0;
    for (let j = 0; j < d; j++) {
      if (i === j) continue;
      rIdx += R.get(i, j) ** 2;
      qIdx += antiImage[i][j] ** 2;
    }
    msaItems[i] = rIdx / Math.max(rIdx + qIdx, 1e-12);
    rSum += rIdx;
    qSum += qIdx;
  }
  const kmo = rSum / Math.max(rSum + qSum, 1e-12);
  let logDetR;
  try {
    logDetR = R.logDet();
  } catch {
    const ev = R.eigen().values;
    logDetR = ev.reduce((s, v) => s + Math.log(Math.max(v, 1e-15)), 0);
  }
  const chiSq = -(n - 1 - (2 * d + 5) / 6) * logDetR;
  const df = d * (d - 1) / 2;
  const pValue = chiSqPValue(Math.max(chiSq, 0), df);
  return {
    kmo,
    kmoPerItem: Array.from(msaItems),
    bartlett: {
      chiSq: Math.max(chiSq, 0),
      df,
      pValue
    }
  };
}
function computeImpliedCov(L, Phi, Theta) {
  const d = L.length;
  const k = L[0].length;
  const sigma = Array.from(
    { length: d },
    (_, i) => Array.from({ length: d }, (_2, j) => {
      let sum = 0;
      for (let f1 = 0; f1 < k; f1++) {
        for (let f2 = 0; f2 < k; f2++) {
          sum += L[i][f1] * Phi[f1][f2] * L[j][f2];
        }
      }
      return sum + (i === j ? Theta[i] : 0);
    })
  );
  return Matrix.fromArray(sigma);
}
function cfaOptimize(S, model, d, maxIter, tol) {
  const factors = Object.keys(model);
  const k = factors.length;
  const L = Array.from({ length: d }, () => new Array(k).fill(0));
  const Theta = new Float64Array(d).fill(0.5);
  const Phi = Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
  );
  factors.forEach((f, c) => {
    const items = model[f];
    for (const r of items) {
      if (r < d) L[r][c] = 0.7;
    }
  });
  let converged = false;
  let iter = 0;
  const c1 = 1e-4;
  function objective(Lc, Phic, Thetac) {
    const Sigma = computeImpliedCov(Lc, Phic, Thetac);
    let logDetSigma;
    try {
      logDetSigma = Sigma.logDet();
    } catch {
      return 1e10;
    }
    let logDetS;
    try {
      logDetS = S.logDet();
    } catch {
      return 1e10;
    }
    let trVal;
    try {
      trVal = Sigma.inverse().multiply(S).trace();
    } catch {
      return 1e10;
    }
    return Math.max(0, logDetSigma + trVal - logDetS - d);
  }
  for (iter = 0; iter < maxIter; iter++) {
    const Sigma = computeImpliedCov(L, Phi, Theta);
    let invSigma;
    try {
      invSigma = Sigma.inverse();
    } catch {
      try {
        invSigma = Sigma.pseudoInverse();
      } catch {
        break;
      }
    }
    const Delta = invSigma.multiply(Sigma.subtract(S)).multiply(invSigma);
    const LMat = Matrix.fromArray(L);
    const PhiMat = Matrix.fromArray(Phi);
    const gL = Delta.multiply(LMat).multiply(PhiMat).scale(2);
    const gPhi = LMat.transpose().multiply(Delta).multiply(LMat);
    const currentLoss = objective(L, Phi, Theta);
    let alpha = 1;
    let armijoSatisfied = false;
    while (!armijoSatisfied && alpha > 1e-6) {
      let maxGrad = 0;
      let dirDotGrad = 0;
      const nextL = Array.from({ length: d }, () => new Array(k).fill(0));
      factors.forEach((f, c) => {
        for (const r of model[f]) {
          if (r >= d) continue;
          const grad = gL.get(r, c);
          maxGrad = Math.max(maxGrad, Math.abs(grad));
          nextL[r][c] = L[r][c] - alpha * grad;
          dirDotGrad += grad * -grad;
        }
      });
      const nextTheta = new Float64Array(d);
      for (let i = 0; i < d; i++) {
        const grad = Delta.get(i, i);
        maxGrad = Math.max(maxGrad, Math.abs(grad));
        nextTheta[i] = Math.max(1e-3, Theta[i] - alpha * grad);
        dirDotGrad += grad * -grad;
      }
      const nextPhi = Array.from(
        { length: k },
        (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
      );
      for (let r = 0; r < k; r++) {
        for (let c = 0; c < r; c++) {
          const grad = gPhi.get(r, c) + gPhi.get(c, r);
          maxGrad = Math.max(maxGrad, Math.abs(grad));
          const updated = Math.max(-0.99, Math.min(0.99, Phi[r][c] - alpha * grad));
          nextPhi[r][c] = updated;
          nextPhi[c][r] = updated;
          dirDotGrad += grad * -grad;
        }
      }
      if (maxGrad < tol) {
        converged = true;
        break;
      }
      const nextLoss = objective(nextL, nextPhi, nextTheta);
      if (nextLoss <= currentLoss + c1 * alpha * dirDotGrad) {
        armijoSatisfied = true;
        for (let i = 0; i < d; i++) {
          for (let j = 0; j < k; j++) L[i][j] = nextL[i][j];
          Theta[i] = nextTheta[i];
        }
        for (let i = 0; i < k; i++) {
          for (let j = 0; j < k; j++) Phi[i][j] = nextPhi[i][j];
        }
      } else {
        alpha *= 0.5;
      }
    }
    if (converged) break;
    if (!armijoSatisfied) break;
  }
  return { L, Theta, Phi, converged, iterations: iter };
}
function cfaStandardErrors(L, Phi, Theta, S, model, n, d) {
  const factors = Object.keys(model);
  const k = factors.length;
  const params = [];
  factors.forEach((f, c) => {
    for (const r of model[f]) {
      if (r < d) params.push({ type: "loading", i: r, j: c });
    }
  });
  for (let i = 0; i < d; i++) params.push({ type: "theta", i, j: 0 });
  for (let r = 0; r < k; r++) {
    for (let c = 0; c < r; c++) params.push({ type: "phi", i: r, j: c });
  }
  const nParams = params.length;
  function getParam(idx) {
    const p = params[idx];
    if (p.type === "loading") return L[p.i][p.j];
    if (p.type === "theta") return Theta[p.i];
    return Phi[p.i][p.j];
  }
  function setParam(idx, val) {
    const p = params[idx];
    if (p.type === "loading") L[p.i][p.j] = val;
    else if (p.type === "theta") Theta[p.i] = val;
    else {
      Phi[p.i][p.j] = val;
      Phi[p.j][p.i] = val;
    }
  }
  function obj() {
    const Sigma = computeImpliedCov(L, Phi, Theta);
    try {
      const logDetSig = Sigma.logDet();
      const logDetS = S.logDet();
      const tr = Sigma.inverse().multiply(S).trace();
      return Math.max(0, logDetSig + tr - logDetS - d);
    } catch {
      return 1e10;
    }
  }
  const h = 1e-4;
  const hessian = Array.from(
    { length: nParams },
    () => new Array(nParams).fill(0)
  );
  for (let i = 0; i < nParams; i++) {
    for (let j = i; j < nParams; j++) {
      const vi = getParam(i);
      const vj = getParam(j);
      if (i === j) {
        setParam(i, vi + h);
        const fPlus = obj();
        setParam(i, vi - h);
        const fMinus = obj();
        setParam(i, vi);
        const f0 = obj();
        hessian[i][i] = (fPlus - 2 * f0 + fMinus) / (h * h);
      } else {
        setParam(i, vi + h);
        setParam(j, vj + h);
        const fpp = obj();
        setParam(i, vi + h);
        setParam(j, vj - h);
        const fpm = obj();
        setParam(i, vi - h);
        setParam(j, vj + h);
        const fmp = obj();
        setParam(i, vi - h);
        setParam(j, vj - h);
        const fmm = obj();
        setParam(i, vi);
        setParam(j, vj);
        hessian[i][j] = (fpp - fpm - fmp + fmm) / (4 * h * h);
        hessian[j][i] = hessian[i][j];
      }
    }
  }
  const infoArr = hessian.map(
    (row) => row.map((v) => (n - 1) / 2 * v)
  );
  let covMat;
  try {
    covMat = Matrix.fromArray(infoArr).inverse();
  } catch {
    try {
      covMat = Matrix.fromArray(infoArr).pseudoInverse();
    } catch {
      const loadingSE2 = Array.from({ length: d }, () => new Array(k).fill(0.05));
      const thetaSE2 = new Float64Array(d).fill(0.05);
      const phiSE2 = Array.from({ length: k }, () => new Array(k).fill(0.05));
      return { loadingSE: loadingSE2, thetaSE: thetaSE2, phiSE: phiSE2 };
    }
  }
  const loadingSE = Array.from({ length: d }, () => new Array(k).fill(0));
  const thetaSE = new Float64Array(d);
  const phiSE = Array.from({ length: k }, () => new Array(k).fill(0));
  for (let idx = 0; idx < nParams; idx++) {
    const variance2 = Math.max(covMat.get(idx, idx), 0);
    const se2 = Math.sqrt(variance2);
    const p = params[idx];
    if (p.type === "loading") loadingSE[p.i][p.j] = se2;
    else if (p.type === "theta") thetaSE[p.i] = se2;
    else phiSE[p.i][p.j] = se2;
  }
  return { loadingSE, thetaSE, phiSE };
}
function runFADiagnostics(data, options) {
  const n = data.length;
  if (n < 3) throw new Error("runFADiagnostics: need at least 3 observations");
  const d = data[0].length;
  if (d < 2) throw new Error("runFADiagnostics: need at least 2 variables");
  const seed = options?.seed ?? 42;
  const iterations = options?.parallelIterations ?? 100;
  const rng = new PRNG2(seed);
  const R = computeCorrelationMatrix(data, n, d);
  const { kmo, kmoPerItem, bartlett } = computeKMOBartlett(R, n, d);
  const eigenvalues = R.eigen().values;
  const { simulated, suggested: parallelSuggested } = parallelAnalysis(
    eigenvalues,
    n,
    d,
    iterations,
    rng
  );
  const mapSuggested = velicerMAP(R);
  return {
    kmo,
    kmoPerItem,
    bartlett,
    mapSuggested,
    parallelEigenvalues: [...eigenvalues],
    parallelSimulated: [...simulated],
    parallelSuggested
  };
}
function runEFA(data, options) {
  const n = data.length;
  if (n < 3) throw new Error("runEFA: need at least 3 observations");
  const d = data[0].length;
  if (d < 2) throw new Error("runEFA: need at least 2 variables");
  const extraction = options?.extraction ?? "ml";
  const rotation = options?.rotation ?? "promax";
  const geominDelta = options?.geominDelta ?? 0.01;
  const maxIter = options?.maxIter ?? 1e3;
  const tol = options?.tol ?? 1e-6;
  const seed = options?.seed ?? 42;
  const randomStarts = options?.randomStarts ?? 50;
  const R = computeCorrelationMatrix(data, n, d);
  const eigenvalues = R.eigen().values;
  let nFactors = options?.nFactors;
  if (nFactors === void 0) {
    const rng = new PRNG2(seed);
    const { suggested } = parallelAnalysis(eigenvalues, n, d, 100, rng);
    nFactors = suggested;
  }
  if (nFactors < 1) throw new Error("runEFA: nFactors must be at least 1");
  if (nFactors >= d) throw new Error("runEFA: nFactors must be less than number of variables");
  const extracted = extraction === "ml" ? extractML(R, nFactors, maxIter, tol) : extractPAF(R, nFactors, maxIter, tol);
  const { rotated, Phi } = applyRotation(
    extracted.loadings,
    rotation,
    maxIter,
    tol,
    geominDelta,
    randomStarts,
    seed
  );
  const communalities = new Float64Array(d);
  const uniqueness = new Float64Array(d);
  for (let i = 0; i < d; i++) {
    let comm = 0;
    for (let f1 = 0; f1 < nFactors; f1++) {
      for (let f2 = 0; f2 < nFactors; f2++) {
        comm += rotated[i][f1] * Phi[f1][f2] * rotated[i][f2];
      }
    }
    communalities[i] = Math.max(1e-3, Math.min(0.9999, comm));
    uniqueness[i] = 1 - communalities[i];
  }
  const thetaArr = new Float64Array(d);
  for (let i = 0; i < d; i++) thetaArr[i] = Math.max(1e-3, uniqueness[i]);
  const impliedSigma = computeImpliedCov(rotated, Phi, thetaArr);
  const nLoadingParams = d * nFactors;
  const nUniqueParams = d;
  const nRotationConstraints = nFactors * (nFactors - 1) / 2;
  const nFreeParams = nLoadingParams + nUniqueParams - nRotationConstraints;
  const fit = computeFit(R, impliedSigma, n, d, nFreeParams, nFactors);
  const stdLoadings = Array.from(
    { length: d },
    (_, i) => Array.from({ length: nFactors }, (_2, j) => {
      const sigmaII = impliedSigma.get(i, i);
      return rotated[i][j] * Math.sqrt(Phi[j][j]) / Math.sqrt(Math.max(sigmaII, 1e-12));
    })
  );
  const variableNames = options?.variableNames ?? Array.from({ length: d }, (_, i) => `V${i + 1}`);
  const factorNames = Array.from({ length: nFactors }, (_, i) => `F${i + 1}`);
  const formatted = `EFA (${extraction}/${rotation}): ${formatCFAFit(fit)}`;
  return {
    loadings: rotated,
    standardizedLoadings: stdLoadings,
    uniqueness: Array.from(uniqueness),
    communalities: Array.from(communalities),
    factorCorrelations: Phi,
    fit,
    eigenvalues: [...eigenvalues],
    nFactors,
    rotation,
    extraction,
    variableNames,
    factorNames,
    formatted
  };
}
function runCFA(data, model, options) {
  const n = data.length;
  if (n < 3) throw new Error("runCFA: need at least 3 observations");
  const d = data[0].length;
  if (d < 2) throw new Error("runCFA: need at least 2 variables");
  const factors = Object.keys(model);
  const k = factors.length;
  if (k < 1) throw new Error("runCFA: model must specify at least 1 factor");
  const maxIter = options?.maxIter ?? 1e3;
  const tolVal = options?.tol ?? 1e-6;
  const S = computeCorrelationMatrix(data, n, d);
  const { L, Theta, Phi } = cfaOptimize(S, model, d, maxIter, tolVal);
  const impliedSigma = computeImpliedCov(L, Phi, Theta);
  let nLoadingParams = 0;
  factors.forEach((f) => {
    nLoadingParams += model[f].length;
  });
  const nUniqueParams = d;
  const nCovParams = k * (k - 1) / 2;
  const nFreeParams = nLoadingParams + nUniqueParams + nCovParams;
  const fit = computeFit(S, impliedSigma, n, d, nFreeParams);
  const Lcopy = L.map((r) => [...r]);
  const Phicopy = Phi.map((r) => [...r]);
  const Thetacopy = new Float64Array(Theta);
  const { loadingSE, thetaSE, phiSE } = cfaStandardErrors(
    Lcopy,
    Phicopy,
    Thetacopy,
    S,
    model,
    n,
    d
  );
  const paramLoadings = factors.map(
    (f, c) => model[f].map((r) => {
      if (r >= d) {
        return { estimate: 0, se: 0, z: 0, pValue: 1, stdAll: 0 };
      }
      const est = L[r][c];
      const se2 = Math.max(loadingSE[r][c], 1e-6);
      const z = est / se2;
      const pValue = 2 * (1 - normalCDF(Math.abs(z)));
      const sigmaII = Math.max(impliedSigma.get(r, r), 1e-12);
      const stdAll = est * Math.sqrt(Phi[c][c]) / Math.sqrt(sigmaII);
      return { estimate: roundTo(est, 4), se: roundTo(se2, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(stdAll, 4) };
    })
  );
  const paramUniqueness = Array.from(Theta).map((est, i) => {
    const se2 = Math.max(thetaSE[i], 1e-6);
    const z = est / se2;
    const pValue = 2 * (1 - normalCDF(Math.abs(z)));
    const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12);
    const stdAll = est / sigmaII;
    return { estimate: roundTo(est, 4), se: roundTo(se2, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(stdAll, 4) };
  });
  const paramCov = Array.from(
    { length: k },
    (_, r) => Array.from({ length: k }, (_2, c) => {
      if (r === c) {
        return { estimate: 1, se: 0, z: Infinity, pValue: 0, stdAll: 1 };
      }
      const est = Phi[r][c];
      const se2 = Math.max(phiSE[r][c] || phiSE[c][r], 1e-6);
      const z = est / se2;
      const pValue = 2 * (1 - normalCDF(Math.abs(z)));
      return { estimate: roundTo(est, 4), se: roundTo(se2, 4), z: roundTo(z, 3), pValue: roundTo(pValue, 4), stdAll: roundTo(est, 4) };
    })
  );
  const communalities = Array.from(Theta).map((t) => roundTo(1 - t, 4));
  const uniquenessArr = Array.from(Theta).map((t) => roundTo(t, 4));
  const stdLoadings = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12);
      return roundTo(L[i][j] * Math.sqrt(Phi[j][j]) / Math.sqrt(sigmaII), 4);
    })
  );
  const eigenvalues = S.eigen().values.map((v) => roundTo(v, 4));
  const variableNames = options?.variableNames ?? Array.from({ length: d }, (_, i) => `V${i + 1}`);
  const factorNames = options?.factorNames ?? factors;
  const formatted = `CFA: ${formatCFAFit(fit)}`;
  return {
    loadings: L,
    standardizedLoadings: stdLoadings,
    uniqueness: uniquenessArr,
    communalities,
    factorCorrelations: Phi,
    fit,
    eigenvalues,
    nFactors: k,
    rotation: "none",
    extraction: "ml",
    variableNames,
    factorNames,
    formatted,
    parameterEstimates: {
      loadings: paramLoadings,
      uniquenesses: paramUniqueness,
      factorCovariances: paramCov
    },
    model
  };
}

// src/viz/themes/default.ts
var CARM_PALETTE = [
  "#4e79a7",
  // cornflower blue
  "#f28e2b",
  // amber orange
  "#e15759",
  // soft crimson
  "#76b7b2",
  // dusty teal
  "#59a14f",
  // forest green
  "#af7aa1",
  // dusty mauve
  "#ff9da7",
  // rose pink
  "#9c755f"
  // warm umber
];
var OKABE_ITO = CARM_PALETTE;
var DEFAULT_THEME = {
  background: "#ffffff",
  surface: "#f8f9fa",
  text: "#1a1a2e",
  textMuted: "#6c757d",
  textAnnotation: "#495057",
  gridLine: "#eaeef3",
  axisLine: "#c4cdd6",
  colors: CARM_PALETTE,
  fontFamily: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
  fontFamilyMono: "'JetBrains Mono', 'Fira Code', 'Cascadia Code', Menlo, Consolas, monospace",
  fontSize: 12,
  fontSizeSmall: 11,
  fontSizeTitle: 16,
  marginTop: 58,
  marginRight: 32,
  marginBottom: 84,
  marginLeft: 64,
  pointOpacity: 0.55,
  violinOpacity: 0.72,
  ciOpacity: 0.15
};
var DARK_THEME = {
  ...DEFAULT_THEME,
  background: "#14142b",
  surface: "#1e1e3f",
  text: "#e9ecef",
  textMuted: "#9aa5b1",
  textAnnotation: "#ced4da",
  gridLine: "#252548",
  axisLine: "#3d3d6b"
};
function applyTheme(container, theme = DEFAULT_THEME) {
  const vars = {
    "--js-bg": theme.background,
    "--js-surface": theme.surface,
    "--js-text": theme.text,
    "--js-text-muted": theme.textMuted,
    "--js-text-annotation": theme.textAnnotation,
    "--js-grid": theme.gridLine,
    "--js-axis": theme.axisLine,
    "--js-font": theme.fontFamily,
    "--js-font-mono": theme.fontFamilyMono,
    "--js-font-size": `${theme.fontSize}px`,
    "--js-font-size-sm": `${theme.fontSizeSmall}px`,
    "--js-font-size-title": `${theme.fontSizeTitle}px`
  };
  theme.colors.forEach((c, i) => {
    vars[`--js-color-${i}`] = c;
  });
  for (const [k, v] of Object.entries(vars)) {
    container.style.setProperty(k, v);
  }
  container.style.background = theme.background;
}
function getColor(index, theme = DEFAULT_THEME) {
  return theme.colors[index % theme.colors.length];
}
function themeColorScale(theme = DEFAULT_THEME) {
  return (i) => getColor(i, theme);
}

// src/viz/components/tooltip.ts
var tooltipEl = null;
function ensureTooltip() {
  if (!tooltipEl) {
    tooltipEl = document.createElement("div");
    tooltipEl.id = "carm-tooltip";
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
    `;
    document.body.appendChild(tooltipEl);
  }
  return tooltipEl;
}
function showTooltip(event, content, theme = DEFAULT_THEME) {
  const el = ensureTooltip();
  el.innerHTML = content;
  el.style.background = theme.surface;
  el.style.color = theme.text;
  el.style.border = `1px solid ${theme.gridLine}`;
  el.style.fontFamily = theme.fontFamily;
  const x = event.clientX + 14;
  const y = event.clientY - 28;
  const viewW = window.innerWidth;
  const elW = 240;
  el.style.left = `${Math.min(x, viewW - elW - 8)}px`;
  el.style.top = `${Math.max(y, 8)}px`;
  el.style.opacity = "1";
}
function hideTooltip() {
  if (tooltipEl) tooltipEl.style.opacity = "0";
}
function formatTooltipRow(label, value) {
  return `<div style="display:flex;justify-content:space-between;gap:12px">
    <span style="opacity:0.7">${label}</span>
    <strong>${typeof value === "number" ? value.toFixed(3) : value}</strong>
  </div>`;
}

// src/viz/components/annotations.ts
function addSubtitle(svg, title, subtitle, _width, theme = DEFAULT_THEME) {
  svg.append("rect").attr("x", 20).attr("y", 10).attr("width", 3).attr("height", 22).attr("rx", 1.5).attr("fill", theme.colors[0] ?? "#4e79a7");
  svg.append("text").attr("class", "plot-title").attr("x", 30).attr("y", 26).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeTitle).attr("font-weight", "700").attr("letter-spacing", "-0.3").attr("fill", theme.text).text(title);
  if (subtitle) {
    svg.append("text").attr("class", "plot-subtitle").attr("x", 30).attr("y", 45).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-style", "italic").attr("fill", theme.textAnnotation).text(subtitle);
  }
}
function addCaption(svg, text, _width, height, theme = DEFAULT_THEME) {
  svg.append("text").attr("x", 20).attr("y", height - 8).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).style("font-style", "italic").text(text);
}
function addRegressionEquation(g, intercept, slope, r2, x, y, theme = DEFAULT_THEME) {
  const sign = slope >= 0 ? "+" : "\u2212";
  const eq = `\u0177 = ${intercept.toFixed(2)} ${sign} ${Math.abs(slope).toFixed(2)}x   R\xB2 = ${r2.toFixed(3)}`;
  const padX = 8, padY = 4;
  const textW = eq.length * 6.2 + padX * 2;
  const textH = 14 + padY * 2;
  g.append("rect").attr("x", x - padX).attr("y", y - textH + padY).attr("width", textW).attr("height", textH).attr("rx", 4).attr("fill", theme.surface).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("text").attr("x", x).attr("y", y).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(eq);
}
function addNLabel(g, n, x, y, theme = DEFAULT_THEME) {
  g.append("text").attr("x", x).attr("y", y).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).text(`n = ${n}`);
}
function addStatBadge(g, lines, x, y, theme = DEFAULT_THEME) {
  const lineH = 14;
  const padX = 10, padY = 6;
  const maxLen = Math.max(...lines.map((l) => l.length));
  const bw = maxLen * 6.4 + padX * 2;
  const bh = lines.length * lineH + padY * 2;
  g.append("rect").attr("x", x).attr("y", y).attr("width", bw).attr("height", bh).attr("rx", 5).attr("fill", theme.surface).attr("stroke", theme.gridLine).attr("stroke-width", 1).attr("opacity", 0.92);
  lines.forEach((line, i) => {
    g.append("text").attr("x", x + padX).attr("y", y + padY + (i + 1) * lineH - 2).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(line);
  });
}

// src/viz/components/brackets.ts
function formatBracketP(p, numeric) {
  if (numeric) {
    if (p < 1e-3) return "p < .001";
    return `p = ${p.toFixed(3).replace(/^0\./, ".")}`;
  }
  if (p < 1e-3) return "***";
  if (p < 0.01) return "**";
  if (p < 0.05) return "*";
  return "ns";
}
function renderBrackets(g, comparisons, config, theme = DEFAULT_THEME) {
  const toShow = config.significantOnly ? comparisons.filter((c) => c.significant) : comparisons;
  if (toShow.length === 0) return;
  const brackets = toShow.map((c) => ({
    x1: Math.min(
      config.groupPositions.get(c.group1) ?? 0,
      config.groupPositions.get(c.group2) ?? 0
    ),
    x2: Math.max(
      config.groupPositions.get(c.group1) ?? 0,
      config.groupPositions.get(c.group2) ?? 0
    ),
    p: c.pValueAdj,
    label: formatBracketP(c.pValueAdj, config.numericP !== false),
    level: 0
  })).sort((a, b) => a.x2 - a.x1 - (b.x2 - b.x1));
  const levelMaxX = [];
  for (const bracket of brackets) {
    let level = 0;
    while ((levelMaxX[level] ?? -Infinity) >= bracket.x1 - 5) {
      level++;
    }
    bracket.level = level;
    levelMaxX[level] = bracket.x2;
  }
  for (const b of brackets) {
    const y = config.yBase - b.level * config.bracketHeight - 12;
    const tipLen = 6;
    const color = b.label === "ns" ? theme.textMuted : theme.text;
    g.append("line").attr("x1", b.x1).attr("x2", b.x2).attr("y1", y).attr("y2", y).attr("stroke", color).attr("stroke-width", 1.2);
    g.append("line").attr("x1", b.x1).attr("x2", b.x1).attr("y1", y).attr("y2", y + tipLen).attr("stroke", color).attr("stroke-width", 1.2);
    g.append("line").attr("x1", b.x2).attr("x2", b.x2).attr("y1", y).attr("y2", y + tipLen).attr("stroke", color).attr("stroke-width", 1.2);
    g.append("text").attr("x", (b.x1 + b.x2) / 2).attr("y", y - 3).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", color).text(b.label);
  }
}
function totalBracketHeight(comparisons, significantOnly, bracketHeight) {
  const n = significantOnly ? comparisons.filter((c) => c.significant).length : comparisons.length;
  return n * bracketHeight + 20;
}

// src/viz/components/axis.ts
function renderXAxis(g, _scale, height, label, width, theme = DEFAULT_THEME) {
  g.append("g").attr("transform", `translate(0,${height})`).call((g_) => {
    g_.selectAll("line").attr("stroke", theme.axisLine);
    g_.selectAll("path").attr("stroke", theme.axisLine);
    g_.selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  });
  g.append("text").attr("x", width / 2).attr("y", height + 48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(label);
}
function renderYAxis(g, _scale, height, label, theme = DEFAULT_THEME) {
  g.append("g").call((d3Axis) => {
    void d3Axis;
  });
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(label);
}
function renderGridLines(g, scale, width, theme = DEFAULT_THEME) {
  const ticks = scale.ticks(6);
  g.selectAll(".grid-line").data(ticks).join("line").attr("class", "grid-line").attr("x1", 0).attr("x2", width).attr("y1", scale).attr("y2", scale).attr("stroke", theme.gridLine).attr("stroke-width", 1);
}

// src/viz/plots/violin-box.ts
function renderViolinBox(container, data, config = {}) {
  import("d3").then((d3) => renderViolinBoxD3(d3, container, data, config));
}
function renderViolinBoxD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = {
    top: theme.marginTop + (data.pairwise?.length ?? 0) * 20,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  if (config.title || data.testResult) {
    addSubtitle(
      svg,
      config.title ?? "Group Comparison",
      data.testResult?.formatted ?? "",
      W,
      theme
    );
  }
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const labels = data.groups.map((gr) => gr.label);
  const xScale = d3.scaleBand().domain(labels).range([0, width]).padding(0.2);
  const allValues = data.groups.flatMap((gr) => [...gr.values]);
  const [yMin, yMax] = d3.extent(allValues);
  const yPad = (yMax - yMin) * 0.1;
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  data.groups.forEach((gr, gi) => {
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2;
    const color = getColor(gi, theme);
    const bw = config.violinBandwidth ?? (yMax - yMin) / 20;
    const kdePoints = kernelDensityEstimate(gr.values, yScale.ticks(60), bw);
    const maxDensity = Math.max(...kdePoints.map((p) => p[1]));
    const violinWidth = xScale.bandwidth() * 0.45;
    const violinScale = maxDensity > 0 ? violinWidth / maxDensity : 1;
    const areaFn = d3.area().x0((d) => cx - d[1] * violinScale).x1((d) => cx + d[1] * violinScale).y((d) => yScale(d[0])).curve(d3.curveCatmullRom);
    g.append("path").datum(kdePoints).attr("d", areaFn).attr("fill", color).attr("opacity", theme.violinOpacity).attr("stroke", color).attr("stroke-width", 1);
    const q1 = quantile(gr.values, 0.25);
    const med = quantile(gr.values, 0.5);
    const q3 = quantile(gr.values, 0.75);
    const iqr = q3 - q1;
    const whiskerLo = Math.min(...gr.values.filter((v) => v >= q1 - 1.5 * iqr));
    const whiskerHi = Math.max(...gr.values.filter((v) => v <= q3 + 1.5 * iqr));
    const boxW = xScale.bandwidth() * 0.2;
    g.append("line").attr("x1", cx).attr("x2", cx).attr("y1", yScale(whiskerLo)).attr("y2", yScale(q1)).attr("stroke", color).attr("stroke-width", 1.5);
    g.append("line").attr("x1", cx).attr("x2", cx).attr("y1", yScale(q3)).attr("y2", yScale(whiskerHi)).attr("stroke", color).attr("stroke-width", 1.5);
    g.append("rect").attr("x", cx - boxW / 2).attr("width", boxW).attr("y", yScale(q3)).attr("height", yScale(q1) - yScale(q3)).attr("fill", theme.background).attr("stroke", color).attr("stroke-width", 2);
    if (config.showMedian !== false) {
      g.append("line").attr("x1", cx - boxW / 2).attr("x2", cx + boxW / 2).attr("y1", yScale(med)).attr("y2", yScale(med)).attr("stroke", color).attr("stroke-width", 2.5);
      g.append("text").attr("x", cx + boxW / 2 + 4).attr("y", yScale(med) + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(med.toFixed(2));
    }
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length;
      const my = yScale(groupMean);
      const ds = 5;
      g.append("polygon").attr("points", `${cx},${my - ds} ${cx + ds},${my} ${cx},${my + ds} ${cx - ds},${my}`).attr("fill", "white").attr("stroke", color).attr("stroke-width", 1.5);
      const medY = yScale(med);
      const tooClose = config.showMedian !== false && Math.abs(my - medY) < 14;
      if (!tooClose) {
        g.append("text").attr("x", cx + ds + 4).attr("y", my + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(groupMean.toFixed(2));
      }
    }
    if (config.showJitter !== false) {
      const jw = (config.jitterWidth ?? 0.15) * xScale.bandwidth();
      const seed = gi * 12345;
      gr.values.forEach((v, vi) => {
        const jx = cx + (pseudoRandom(seed + vi) - 0.5) * jw;
        g.append("circle").attr("cx", jx).attr("cy", yScale(v)).attr("r", 3).attr("fill", color).attr("opacity", theme.pointOpacity).attr("stroke", "none").on("mouseover", (event) => {
          showTooltip(event, [
            formatTooltipRow("Group", gr.label),
            formatTooltipRow("Value", v.toFixed(3))
          ].join(""), theme);
        }).on("mouseout", hideTooltip);
      });
    }
    if (config.showN !== false) {
      addNLabel(g, gr.values.length, cx, height + 60, theme);
    }
  });
  if (config.showBrackets !== false && data.pairwise && data.pairwise.length > 0) {
    const posMap = new Map(
      data.groups.map((gr) => [gr.label, (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2])
    );
    renderBrackets(g, data.pairwise, {
      groupPositions: posMap,
      yBase: 0,
      bracketHeight: 22,
      significantOnly: config.significantBracketsOnly ?? false,
      ...config.numericP !== void 0 && { numericP: config.numericP }
    }, theme);
  }
  if (config.caption) {
    addCaption(svg, config.caption, W, H, theme);
  }
}
function kernelDensityEstimate(data, xPoints, bandwidth) {
  return xPoints.map((x) => {
    const density = data.reduce((sum, xi) => sum + gaussianKernel((x - xi) / bandwidth), 0) / (data.length * bandwidth);
    return [x, density];
  });
}
function gaussianKernel(u) {
  return Math.exp(-0.5 * u * u) / Math.sqrt(2 * Math.PI);
}
function pseudoRandom(seed) {
  const x = Math.sin(seed) * 1e4;
  return x - Math.floor(x);
}

// src/viz/plots/scatter-stats.ts
function renderScatterStats(container, data, config = {}) {
  import("d3").then((d3) => renderScatterD3(d3, container, data, config));
}
function renderScatterD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight + 40, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Scatter Plot", data.correlationResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const [xMin, xMax] = d3.extent(data.x);
  const [yMin, yMax] = d3.extent(data.y);
  const xPad = (xMax - xMin) * 0.05;
  const yPad = (yMax - yMin) * 0.05;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid-h").data(yScale.ticks(6)).join("line").attr("class", "grid-h").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "X");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Y");
  if (data.regressionResult) {
    const coefs = data.regressionResult.coefficients;
    const b0 = coefs[0]?.estimate ?? 0;
    const b1 = coefs[1]?.estimate ?? 1;
    const r2 = data.regressionResult.r2;
    const xRange = xScale.domain();
    const lineData = [
      [xRange[0], b0 + b1 * xRange[0]],
      [xRange[1], b0 + b1 * xRange[1]]
    ];
    if (config.showCI !== false) {
      const n = data.x.length;
      const se_reg = Math.sqrt(data.regressionResult.residuals.reduce((s, r) => s + r * r, 0) / (n - 2));
      const xMean = data.x.reduce((s, v) => s + v, 0) / n;
      const sxx = data.x.reduce((s, v) => s + (v - xMean) ** 2, 0);
      const ciPoints = xScale.ticks(50).map((x) => {
        const yHat = b0 + b1 * x;
        const h = 1 / n + (x - xMean) ** 2 / sxx;
        const halfCI = 1.96 * se_reg * Math.sqrt(h);
        return { x, lo: yHat - halfCI, hi: yHat + halfCI };
      });
      const area = d3.area().x((d) => xScale(d.x)).y0((d) => yScale(d.lo)).y1((d) => yScale(d.hi)).curve(d3.curveBasis);
      g.append("path").datum(ciPoints).attr("d", area).attr("fill", getColor(0, theme)).attr("opacity", theme.ciOpacity);
    }
    g.append("line").attr("x1", xScale(lineData[0][0])).attr("x2", xScale(lineData[1][0])).attr("y1", yScale(lineData[0][1])).attr("y2", yScale(lineData[1][1])).attr("stroke", getColor(0, theme)).attr("stroke-width", 2);
    if (config.showEquation !== false) {
      addRegressionEquation(g, b0, b1, r2, 10, 20, theme);
    }
  }
  const color = getColor(0, theme);
  data.x.forEach((xi, i) => {
    const yi = data.y[i] ?? 0;
    g.append("circle").attr("cx", xScale(xi)).attr("cy", yScale(yi)).attr("r", config.pointSize ?? 4).attr("fill", color).attr("opacity", theme.pointOpacity).attr("stroke", theme.background).attr("stroke-width", 0.5).on("mouseover", (event) => {
      const rows = [formatTooltipRow(config.xLabel ?? "X", xi.toFixed(3)), formatTooltipRow(config.yLabel ?? "Y", yi.toFixed(3))];
      if (data.labels?.[i]) rows.unshift(`<strong>${data.labels[i]}</strong>`);
      showTooltip(event, rows.join(""), theme);
    }).on("mouseout", hideTooltip);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/histogram.ts
function renderHistogram(container, data, config = {}) {
  import("d3").then((d3) => renderHistogramD3(d3, container, data, config));
}
function renderHistogramD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 400;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Distribution",
    data.descriptives ? `M = ${data.descriptives.mean.toFixed(2)}, SD = ${data.descriptives.sd.toFixed(2)}, n = ${data.descriptives.n}` : "",
    W,
    theme
  );
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const [xMin, xMax] = d3.extent(data.values);
  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]).nice();
  const nBins = config.bins ?? Math.ceil(Math.sqrt(data.values.length));
  const histogram = d3.bin().domain(xScale.domain()).thresholds(nBins);
  const bins = histogram(data.values);
  const maxCount = Math.max(...bins.map((b) => b.length));
  const yScale = d3.scaleLinear().domain([0, maxCount * 1.1]).range([height, 0]);
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const barColor = config.color ?? getColor(0, theme);
  g.selectAll(".bar").data(bins).join("rect").attr("class", "bar").attr("x", (b) => xScale(b.x0)).attr("y", (b) => yScale(b.length)).attr("width", (b) => Math.max(0, xScale(b.x1) - xScale(b.x0) - 1)).attr("height", (b) => height - yScale(b.length)).attr("fill", barColor).attr("opacity", 0.7).attr("stroke", theme.background).attr("stroke-width", 0.5);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Value");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Count");
  if (config.showNormalCurve !== false && data.descriptives) {
    const m = data.descriptives.mean;
    const s = data.descriptives.sd;
    const n = data.descriptives.n;
    const binWidth = bins[0]?.x1 - bins[0]?.x0;
    const scale_ = n * binWidth;
    const xTicks = d3.range(xMin, xMax, (xMax - xMin) / 100);
    const normalPoints = xTicks.map((x) => ({
      x,
      y: scale_ * Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI))
    }));
    const line = d3.line().x((d) => xScale(d.x)).y((d) => yScale(d.y)).curve(d3.curveBasis);
    g.append("path").datum(normalPoints).attr("d", line).attr("fill", "none").attr("stroke", getColor(4, theme)).attr("stroke-width", 2).attr("stroke-dasharray", "6,3");
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/bar-stats.ts
function renderBarStats(container, data, config = {}) {
  import("d3").then((d3) => renderBarD3(d3, container, data, config));
}
function renderBarD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 400;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Frequency", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const categories = data.rows.map((r) => String(r.value));
  const maxCount = Math.max(...data.rows.map((r) => r.count));
  const xScale = d3.scaleBand().domain(categories).range([0, width]).padding(0.25);
  const yScale = d3.scaleLinear().domain([0, maxCount * 1.15]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.selectAll(".bar").data(data.rows).join("rect").attr("class", "bar").attr("x", (row) => xScale(String(row.value)) ?? 0).attr("y", (row) => yScale(row.count)).attr("width", xScale.bandwidth()).attr("height", (row) => height - yScale(row.count)).attr("fill", (_, i) => getColor(i, theme)).attr("rx", 2);
  if (config.showCounts !== false) {
    g.selectAll(".label-count").data(data.rows).join("text").attr("class", "label-count").attr("x", (row) => (xScale(String(row.value)) ?? 0) + xScale.bandwidth() / 2).attr("y", (row) => yScale(row.count) - 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text((row) => `${row.count}`);
  }
  if (config.showPercentages !== false) {
    g.selectAll(".label-pct").data(data.rows).join("text").attr("class", "label-pct").attr("x", (row) => (xScale(String(row.value)) ?? 0) + xScale.bandwidth() / 2).attr("y", (row) => yScale(row.count) - 16).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text((row) => `(${(row.relative * 100).toFixed(1)}%)`);
  }
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Category");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Count");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/correlogram.ts
function renderCorrelogram(container, data, config = {}) {
  import("d3").then((d3) => renderCorrelogramD3(d3, container, data, config));
}
function renderCorrelogramD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const k = data.labels.length;
  const cellSize = Math.min(
    Math.floor((config.width ?? 500) / (k + 2)),
    Math.floor((config.height ?? 500) / (k + 2)),
    70
  );
  const W = cellSize * k + 120;
  const H = cellSize * k + 100;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Correlation Matrix", "", W, theme);
  const margin = { top: 60, right: 20, bottom: 20, left: 80 };
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu);
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) {
      const r = data.r[i]?.[j] ?? 0;
      const p = data.pValues[i]?.[j] ?? 1;
      g.append("rect").attr("x", j * cellSize).attr("y", i * cellSize).attr("width", cellSize - 2).attr("height", cellSize - 2).attr("rx", 3).attr("fill", i === j ? theme.surface : colorScale(r)).on("mouseover", (event) => {
        if (i !== j) {
          showTooltip(event, [
            formatTooltipRow(`${data.labels[i]} \xD7 ${data.labels[j]}`, ""),
            formatTooltipRow("r", r.toFixed(3)),
            formatTooltipRow("p", p < 1e-3 ? "< .001" : p.toFixed(3))
          ].join(""), theme);
        }
      }).on("mouseout", hideTooltip);
      if (config.showValues !== false && i !== j) {
        g.append("text").attr("x", j * cellSize + cellSize / 2).attr("y", i * cellSize + cellSize / 2 + 1).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", Math.min(theme.fontSizeSmall, cellSize / 3)).attr("fill", Math.abs(r) > 0.6 ? "#fff" : theme.text).text(r.toFixed(2));
      }
      if (config.showSignificance !== false && i !== j && !isNaN(p)) {
        const stars = p < 1e-3 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : "";
        if (stars) {
          g.append("text").attr("x", j * cellSize + cellSize / 2).attr("y", i * cellSize + cellSize - 4).attr("text-anchor", "middle").attr("font-size", 8).attr("fill", Math.abs(r) > 0.6 ? "#fff" : theme.textMuted).text(stars);
        }
      }
    }
  }
  data.labels.forEach((lbl, i) => {
    g.append("text").attr("x", -6).attr("y", i * cellSize + cellSize / 2).attr("text-anchor", "end").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(lbl);
  });
  data.labels.forEach((lbl, j) => {
    g.append("text").attr("x", j * cellSize + cellSize / 2).attr("y", -8).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(lbl);
  });
  if (config.showLegend !== false) {
    const legendW = 80, legendH = 12;
    const legendX = W - 120, legendY = H - 35;
    const defs = svg.append("defs");
    const gradId = "corr-grad-" + Math.random().toString(36).slice(2);
    const grad = defs.append("linearGradient").attr("id", gradId);
    const stops = [-1, -0.5, 0, 0.5, 1];
    stops.forEach((v) => {
      grad.append("stop").attr("offset", `${((v + 1) / 2 * 100).toFixed(0)}%`).attr("stop-color", colorScale(v));
    });
    svg.append("rect").attr("x", legendX).attr("y", legendY).attr("width", legendW).attr("height", legendH).attr("fill", `url(#${gradId})`);
    svg.append("text").attr("x", legendX).attr("y", legendY + legendH + 10).attr("font-size", 9).attr("fill", theme.textMuted).text("\u22121");
    svg.append("text").attr("x", legendX + legendW / 2).attr("y", legendY + legendH + 10).attr("text-anchor", "middle").attr("font-size", 9).attr("fill", theme.textMuted).text("0");
    svg.append("text").attr("x", legendX + legendW).attr("y", legendY + legendH + 10).attr("text-anchor", "end").attr("font-size", 9).attr("fill", theme.textMuted).text("+1");
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/coef-plot.ts
function renderCoefPlot(container, coefficients, config = {}) {
  import("d3").then((d3) => renderCoefD3(d3, container, coefficients, config));
}
function renderCoefD3(d3, container, coefficientsAll, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const coefs = config.excludeIntercept !== false ? coefficientsAll.filter((c) => c.name !== "(Intercept)") : coefficientsAll;
  const k = coefs.length;
  const W = config.width ?? 500;
  const H = config.height ?? Math.max(k * 40 + 100, 200);
  const margin = { top: theme.marginTop, right: theme.marginRight + 60, bottom: theme.marginBottom, left: 120 };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Coefficient Plot", "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allX = coefs.flatMap((c) => [c.ci[0], c.ci[1], c.estimate]);
  const [xMin, xMax] = d3.extent(allX);
  const xPad = (xMax - xMin) * 0.1;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleBand().domain([...coefs].reverse().map((c) => c.name)).range([0, height]).padding(0.3);
  if (config.showZeroLine !== false) {
    const x0 = xScale(0);
    g.append("line").attr("x1", x0).attr("x2", x0).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-dasharray", "4,2").attr("stroke-width", 1);
  }
  g.selectAll(".grid").data(xScale.ticks(6)).join("line").attr("class", "grid").attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d)).attr("y1", 0).attr("y2", height).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  coefs.forEach((coef, _i) => {
    const cy = (yScale(coef.name) ?? 0) + yScale.bandwidth() / 2;
    const significant = coef.pValue < 0.05;
    const color = significant ? getColor(0, theme) : theme.textMuted;
    g.append("line").attr("x1", xScale(coef.ci[0])).attr("x2", xScale(coef.ci[1])).attr("y1", cy).attr("y2", cy).attr("stroke", color).attr("stroke-width", 2.5).attr("opacity", 0.6);
    g.append("circle").attr("cx", xScale(coef.estimate)).attr("cy", cy).attr("r", 5).attr("fill", color);
    const pLabel = coef.pValue < 1e-3 ? "p < .001" : `p = ${coef.pValue.toFixed(3)}`;
    g.append("text").attr("x", width + 5).attr("y", cy + 4).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).text(pLabel);
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Estimate (95% CI)");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/qq-plot.ts
function renderQQPlot(container, values, config = {}) {
  import("d3").then((d3) => renderQQD3(d3, container, values, config));
}
function renderQQD3(d3, container, values, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 400;
  const H = config.height ?? 400;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const n = values.length;
  const sorted = sortAsc(values);
  const mean_ = sorted.reduce((s, v) => s + v, 0) / n;
  const sd_ = Math.sqrt(sorted.reduce((s, v) => s + (v - mean_) ** 2, 0) / (n - 1));
  const points = sorted.map((y, i) => {
    const p = (i + 1 - 0.375) / (n + 0.25);
    const x = normalQuantile(Math.max(1e-4, Math.min(0.9999, p)));
    return { x, y };
  });
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "QQ Plot", "Normal probability plot", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const xExtent = d3.extent(points.map((p) => p.x));
  const yExtent = d3.extent(points.map((p) => p.y));
  const xPad = (xExtent[1] - xExtent[0]) * 0.05;
  const yPad = (yExtent[1] - yExtent[0]) * 0.05;
  const xScale = d3.scaleLinear().domain([xExtent[0] - xPad, xExtent[1] + xPad]).range([0, width]).nice();
  const yScale = d3.scaleLinear().domain([yExtent[0] - yPad, yExtent[1] + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const xDom = xScale.domain();
  const xd0 = xDom[0] ?? 0, xd1 = xDom[1] ?? 1;
  g.append("line").attr("x1", xScale(xd0)).attr("x2", xScale(xd1)).attr("y1", yScale(mean_ + sd_ * xd0)).attr("y2", yScale(mean_ + sd_ * xd1)).attr("stroke", getColor(4, theme)).attr("stroke-width", 2).attr("stroke-dasharray", "6,3");
  if (config.showCI !== false) {
    const CI_MULT = 1.36 / Math.sqrt(n);
    const bandData = points.map((p) => ({
      x: p.x,
      lo: mean_ + sd_ * p.x - CI_MULT * sd_ * Math.sqrt(n),
      hi: mean_ + sd_ * p.x + CI_MULT * sd_ * Math.sqrt(n)
    }));
    const area = d3.area().x((d) => xScale(d.x)).y0((d) => yScale(d.lo)).y1((d) => yScale(d.hi));
    g.append("path").datum(bandData).attr("d", area).attr("fill", getColor(4, theme)).attr("opacity", theme.ciOpacity);
  }
  g.selectAll(".qq-point").data(points).join("circle").attr("class", "qq-point").attr("cx", (p) => xScale(p.x)).attr("cy", (p) => yScale(p.y)).attr("r", 3).attr("fill", getColor(0, theme)).attr("opacity", theme.pointOpacity);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Theoretical Quantiles");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Sample Quantiles");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/residual-panel.ts
function renderResidualPanel(container, result, leverage, config = {}) {
  import("d3").then((d3) => renderResidualD3(d3, container, result, leverage, config));
}
function renderResidualD3(d3, container, result, leverage, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 700;
  const H = config.height ?? 600;
  const panelW = W / 2, panelH = H / 2;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  svg.append("text").attr("x", W / 2).attr("y", 18).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeTitle).attr("font-weight", "600").attr("fill", theme.text).text(config.title ?? "Regression Diagnostics");
  const color = getColor(0, theme);
  const panels = [
    { title: "Residuals vs Fitted", row: 0, col: 0 },
    { title: "Normal Q-Q", row: 0, col: 1 },
    { title: "Scale-Location", row: 1, col: 0 },
    { title: "Residuals vs Leverage", row: 1, col: 1 }
  ];
  const margin = { top: 40, right: 20, bottom: 50, left: 50 };
  const pw = panelW - margin.left - margin.right;
  const ph = panelH - margin.top - margin.bottom;
  panels.forEach((panel, idx) => {
    const ox = panel.col * panelW + margin.left;
    const oy = panel.row * panelH + margin.top + 25;
    const g = svg.append("g").attr("transform", `translate(${ox},${oy})`);
    g.append("text").attr("x", pw / 2).attr("y", -18).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("font-weight", "600").attr("fill", theme.text).text(panel.title);
    let xData;
    let yData;
    if (idx === 0) {
      xData = result.fitted;
      yData = result.residuals;
    } else if (idx === 1) {
      const sorted = sortAsc(result.residuals);
      const n = sorted.length;
      yData = sorted;
      xData = sorted.map((_, i) => {
        const p = (i + 1 - 0.375) / (n + 0.25);
        return normalQuantile(Math.max(1e-4, Math.min(0.9999, p)));
      });
    } else if (idx === 2) {
      xData = result.fitted;
      yData = result.residuals.map((r) => Math.sqrt(Math.abs(r)));
    } else {
      xData = leverage;
      yData = result.residuals;
    }
    const xExt = d3.extent([...xData]);
    const yExt = d3.extent([...yData]);
    const xPad = (xExt[1] - xExt[0]) * 0.05;
    const yPad = (yExt[1] - yExt[0]) * 0.05;
    const xs = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, pw]).nice();
    const ys = d3.scaleLinear().domain([yExt[0] - yPad, yExt[1] + yPad]).range([ph, 0]).nice();
    g.selectAll(".grid").data(ys.ticks(4)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", pw).attr("y1", (d) => ys(d)).attr("y2", (d) => ys(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
    if (idx === 0 || idx === 2) {
      g.append("line").attr("x1", 0).attr("x2", pw).attr("y1", ys(0)).attr("y2", ys(0)).attr("stroke", theme.axisLine).attr("stroke-dasharray", "4,2");
    }
    if (idx === 1) {
      const mean_ = mean(result.residuals);
      const sd_ = Math.sqrt(result.residuals.reduce((s, r) => s + (r - mean_) ** 2, 0) / (result.residuals.length - 1));
      const xd = xs.domain();
      const xd0 = xd[0] ?? 0, xd1 = xd[1] ?? 1;
      g.append("line").attr("x1", xs(xd0)).attr("x2", xs(xd1)).attr("y1", ys(mean_ + sd_ * xd0)).attr("y2", ys(mean_ + sd_ * xd1)).attr("stroke", getColor(4, theme)).attr("stroke-width", 1.5).attr("stroke-dasharray", "5,2");
    }
    xData.forEach((xi, i) => {
      g.append("circle").attr("cx", xs(xi)).attr("cy", ys(yData[i] ?? 0)).attr("r", 2.5).attr("fill", color).attr("opacity", theme.pointOpacity);
    });
    g.append("g").attr("transform", `translate(0,${ph})`).call(d3.axisBottom(xs).ticks(4)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall);
    g.append("g").call(d3.axisLeft(ys).ticks(4)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall);
  });
}

// src/viz/plots/raincloud.ts
function renderRaincloud(container, data, config = {}) {
  import("d3").then((d3) => renderRaincloudD3(d3, container, data, config));
}
function renderRaincloudD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Raincloud Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const labels = data.groups.map((gr) => gr.label);
  const allValues = data.groups.flatMap((gr) => [...gr.values]);
  const [yMin, yMax] = d3.extent(allValues);
  const yPad = (yMax - yMin) * 0.1;
  const xScale = d3.scaleBand().domain(labels).range([0, width]).padding(0.3);
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  data.groups.forEach((gr, gi) => {
    const color = getColor(gi, theme);
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2;
    const bw = (yMax - yMin) / 20;
    const violinR = xScale.bandwidth() * 0.35;
    const kdePoints = kdeEstimate(gr.values, yScale.ticks(60), bw);
    const maxD = Math.max(...kdePoints.map((p) => p[1]));
    const dScale = maxD > 0 ? violinR / maxD : 1;
    const areaFn = d3.area().x0(cx).x1((d) => cx + d[1] * dScale).y((d) => yScale(d[0])).curve(d3.curveCatmullRom);
    g.append("path").datum(kdePoints).attr("d", areaFn).attr("fill", color).attr("opacity", theme.violinOpacity).attr("stroke", color).attr("stroke-width", 1);
    const q1 = quantile(gr.values, 0.25);
    const med = quantile(gr.values, 0.5);
    const q3 = quantile(gr.values, 0.75);
    const iqr = q3 - q1;
    const wLo = Math.min(...gr.values.filter((v) => v >= q1 - 1.5 * iqr));
    const wHi = Math.max(...gr.values.filter((v) => v <= q3 + 1.5 * iqr));
    const boxW = xScale.bandwidth() * 0.12;
    const bx = cx - violinR * 0.5;
    g.append("line").attr("x1", bx).attr("x2", bx).attr("y1", yScale(wLo)).attr("y2", yScale(wHi)).attr("stroke", color).attr("stroke-width", 1.5);
    g.append("rect").attr("x", bx - boxW / 2).attr("width", boxW).attr("y", yScale(q3)).attr("height", yScale(q1) - yScale(q3)).attr("fill", theme.background).attr("stroke", color).attr("stroke-width", 2);
    g.append("line").attr("x1", bx - boxW / 2).attr("x2", bx + boxW / 2).attr("y1", yScale(med)).attr("y2", yScale(med)).attr("stroke", color).attr("stroke-width", 2.5);
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length;
      const my = yScale(groupMean);
      const ds = 5;
      g.append("polygon").attr("points", `${bx},${my - ds} ${bx + ds},${my} ${bx},${my + ds} ${bx - ds},${my}`).attr("fill", "white").attr("stroke", color).attr("stroke-width", 1.5);
      g.append("text").attr("x", bx + ds + 4).attr("y", my + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(groupMean.toFixed(2));
    }
    if (config.showJitter !== false) {
      const jx = bx - boxW - 4;
      const seed = gi * 99991;
      gr.values.forEach((v, vi) => {
        const jitter = (pseudoRnd(seed + vi) - 0.5) * xScale.bandwidth() * 0.12;
        g.append("circle").attr("cx", jx + jitter).attr("cy", yScale(v)).attr("r", 2.5).attr("fill", color).attr("opacity", theme.pointOpacity).on("mouseover", (event) => {
          showTooltip(event, [formatTooltipRow("Group", gr.label), formatTooltipRow("Value", v.toFixed(3))].join(""), theme);
        }).on("mouseout", hideTooltip);
      });
    }
    if (config.showN !== false) {
      addNLabel(g, gr.values.length, cx, height + 60, theme);
    }
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function kdeEstimate(data, xPoints, bw) {
  return xPoints.map((x) => {
    const d = data.reduce((s, xi) => s + Math.exp(-0.5 * ((x - xi) / bw) ** 2) / (Math.sqrt(2 * Math.PI) * bw), 0) / data.length;
    return [x, d];
  });
}
function pseudoRnd(seed) {
  const x = Math.sin(seed) * 1e4;
  return x - Math.floor(x);
}

// src/viz/plots/mixed-plot.ts
function renderMixedPlot(container, result, blups, config = {}) {
  import("d3").then((d3) => renderMixedD3(d3, container, result, blups, config));
}
function renderMixedD3(d3, container, result, blups, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 600;
  const H = config.height ?? Math.max(blups.length * 20 + 150, 300);
  const margin = { top: theme.marginTop, right: 80, bottom: theme.marginBottom, left: 100 };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Random Effects (BLUPs)", result.formatted, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const sorted = [...blups].sort((a, b) => a.blup - b.blup);
  const groupLabels = sorted.map((b) => String(b.group));
  const blupValues = sorted.map((b) => b.blup);
  const xExt = d3.extent(blupValues);
  const xPad = Math.max(Math.abs(xExt[0]), Math.abs(xExt[1])) * 0.2;
  const xScale = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, width]).nice();
  const yScale = d3.scaleBand().domain(groupLabels).range([height, 0]).padding(0.3);
  g.append("line").attr("x1", xScale(0)).attr("x2", xScale(0)).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-dasharray", "4,2").attr("stroke-width", 1.5);
  g.selectAll(".grid").data(xScale.ticks(5)).join("line").attr("class", "grid").attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d)).attr("y1", 0).attr("y2", height).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  sorted.forEach((b, _i) => {
    const cy = (yScale(String(b.group)) ?? 0) + yScale.bandwidth() / 2;
    const color = b.blup >= 0 ? getColor(0, theme) : getColor(5, theme);
    g.append("circle").attr("cx", xScale(b.blup)).attr("cy", cy).attr("r", 4).attr("fill", color);
    g.append("line").attr("x1", xScale(0)).attr("x2", xScale(b.blup)).attr("y1", cy).attr("y2", cy).attr("stroke", color).attr("stroke-width", 1.5).attr("opacity", 0.5);
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("BLUP (random intercept)");
  svg.append("text").attr("x", W - 10).attr("y", H - 10).attr("text-anchor", "end").attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text(`ICC = ${result.icc.toFixed(3)}`);
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/pca-plot.ts
function renderPCAPlot(container, pca, config = {}) {
  import("d3").then((d3) => {
    const type = config.type ?? "biplot";
    if (type === "biplot") renderBiplot(d3, container, pca, config);
    else if (type === "scree") renderScree(d3, container, pca, config);
    else renderLoadingsHeatmap(d3, container, pca, config);
  });
}
function renderBiplot(d3, container, pca, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 500, H = config.height ?? 500;
  const margin = { top: theme.marginTop, right: 60, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const pct1 = ((pca.varianceExplained[0] ?? 0) * 100).toFixed(1);
  const pct2 = ((pca.varianceExplained[1] ?? 0) * 100).toFixed(1);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "PCA Biplot", `PC1 (${pct1}%) \xD7 PC2 (${pct2}%)`, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const scores = pca.scores;
  const xs = scores.map((s) => s[0] ?? 0), ys = scores.map((s) => s[1] ?? 0);
  const xExt = d3.extent(xs), yExt = d3.extent(ys);
  const xPad = (xExt[1] - xExt[0]) * 0.1, yPad = (yExt[1] - yExt[0]) * 0.1;
  const xScale = d3.scaleLinear().domain([xExt[0] - xPad, xExt[1] + xPad]).range([0, width]).nice();
  const yScale = d3.scaleLinear().domain([yExt[0] - yPad, yExt[1] + yPad]).range([height, 0]).nice();
  g.selectAll(".grid-h").data(yScale.ticks(5)).join("line").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", xScale(0)).attr("x2", xScale(0)).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-dasharray", "4,2");
  g.append("line").attr("x1", 0).attr("x2", width).attr("y1", yScale(0)).attr("y2", yScale(0)).attr("stroke", theme.axisLine).attr("stroke-dasharray", "4,2");
  scores.forEach((s, i) => {
    g.append("circle").attr("cx", xScale(s[0] ?? 0)).attr("cy", yScale(s[1] ?? 0)).attr("r", 3.5).attr("fill", getColor(0, theme)).attr("opacity", theme.pointOpacity);
    if (config.observationLabels?.[i]) {
      g.append("text").attr("x", xScale(s[0] ?? 0) + 5).attr("y", yScale(s[1] ?? 0) + 4).attr("font-size", theme.fontSizeSmall - 2).attr("fill", theme.textMuted).text(config.observationLabels[i]);
    }
  });
  const scale_ = Math.min(width, height) * 0.4;
  pca.loadings.forEach((loading, vi) => {
    const lx = (loading[0] ?? 0) * scale_, ly = (loading[1] ?? 0) * scale_;
    const x0 = xScale(0), y0 = yScale(0);
    const color = getColor(vi % 8 + 1, theme);
    g.append("line").attr("x1", x0).attr("y1", y0).attr("x2", x0 + lx).attr("y2", y0 - ly).attr("stroke", color).attr("stroke-width", 1.5).attr("marker-end", "url(#arrow)");
    g.append("text").attr("x", x0 + lx + 5).attr("y", y0 - ly + 4).attr("font-size", theme.fontSizeSmall).attr("fill", color).text(config.variableLabels?.[vi] ?? `Var${vi + 1}`);
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(`PC1 (${pct1}%)`);
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(`PC2 (${pct2}%)`);
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderScree(d3, container, pca, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 500, H = config.height ?? 350;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Scree Plot", "Eigenvalue by component", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const comps = pca.eigenvalues.map((_, i) => i + 1);
  const xScale = d3.scaleBand().domain(comps).range([0, width]).padding(0.3);
  const yScale = d3.scaleLinear().domain([0, Math.max(...pca.eigenvalues) * 1.1]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", 0).attr("x2", width).attr("y1", yScale(1)).attr("y2", yScale(1)).attr("stroke", theme.axisLine).attr("stroke-dasharray", "6,3").attr("stroke-width", 1.5);
  pca.eigenvalues.forEach((ev, i) => {
    const x = xScale(i + 1) ?? 0;
    g.append("rect").attr("x", x).attr("y", yScale(ev)).attr("width", xScale.bandwidth()).attr("height", height - yScale(ev)).attr("fill", getColor(i, theme)).attr("opacity", 0.8).attr("rx", 2);
    g.append("text").attr("x", x + xScale.bandwidth() / 2).attr("y", yScale(ev) - 4).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(ev.toFixed(2));
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Component");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Eigenvalue");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderLoadingsHeatmap(d3, container, pca, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const nVars = pca.loadings.length;
  const nc = pca.nComponents;
  const cellSize = 50;
  const W = config.width ?? nc * cellSize + 120, H = config.height ?? nVars * cellSize + 100;
  const margin = { top: 60, right: 20, bottom: 40, left: 100 };
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "PCA Loadings", "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu);
  pca.loadings.forEach((row, vi) => {
    row.forEach((val, ci) => {
      g.append("rect").attr("x", ci * cellSize).attr("y", vi * cellSize).attr("width", cellSize - 2).attr("height", cellSize - 2).attr("rx", 3).attr("fill", colorScale(val));
      g.append("text").attr("x", ci * cellSize + cellSize / 2).attr("y", vi * cellSize + cellSize / 2 + 4).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall - 1).attr("fill", Math.abs(val) > 0.6 ? "#fff" : theme.text).text(val.toFixed(2));
    });
    const label = config.variableLabels?.[vi] ?? `Var${vi + 1}`;
    g.append("text").attr("x", -6).attr("y", vi * cellSize + cellSize / 2 + 4).attr("text-anchor", "end").attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label);
  });
  Array.from({ length: nc }, (_, ci) => {
    g.append("text").attr("x", ci * cellSize + cellSize / 2).attr("y", -8).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(`PC${ci + 1}`);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/distribution.ts
function renderDistribution(container, params, config = {}) {
  import("d3").then((d3) => renderDistD3(d3, container, params, config));
}
function renderDistD3(d3, container, params, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 550, H = config.height ?? 350;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right, height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const { xMin, xMax, pdf, cdf } = getDistributionFunctions(params);
  const nPoints = 200;
  const xs = Array.from({ length: nPoints }, (_, i) => xMin + (xMax - xMin) * i / (nPoints - 1));
  const pdfs = xs.map((x) => pdf(x));
  const cdfs = xs.map((x) => cdf(x));
  const subtitle = formatDistributionParams(params);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? `${params.distribution} distribution`, subtitle, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]);
  const maxPDF = Math.max(...pdfs.filter(isFinite));
  const yScale = d3.scaleLinear().domain([0, maxPDF * 1.1]).range([height, 0]);
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const lineData = xs.map((x, i) => ({ x, y: pdfs[i] ?? 0 })).filter((d) => isFinite(d.y));
  if (params.highlightX !== void 0) {
    const hx = params.highlightX;
    const shadedData = lineData.filter((d) => d.x <= hx);
    const area = d3.area().x((d) => xScale(d.x)).y0(height).y1((d) => yScale(d.y));
    g.append("path").datum(shadedData).attr("d", area).attr("fill", getColor(0, theme)).attr("opacity", 0.25);
    const p = cdf(hx);
    g.append("text").attr("x", xScale(hx)).attr("y", 20).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textAnnotation).text(`P(X \u2264 ${hx.toFixed(2)}) = ${p.toFixed(4)}`);
  }
  if (params.showPDF !== false) {
    const line = d3.line().x((d) => xScale(d.x)).y((d) => yScale(d.y)).curve(d3.curveBasis);
    g.append("path").datum(lineData).attr("d", line).attr("fill", "none").attr("stroke", getColor(0, theme)).attr("stroke-width", 2.5);
  }
  if (params.showCDF === true) {
    const cdfData = xs.map((x, i) => ({ x, y: (cdfs[i] ?? 0) * maxPDF })).filter((d) => isFinite(d.y));
    const line = d3.line().x((d) => xScale(d.x)).y((d) => yScale(d.y)).curve(d3.curveBasis);
    g.append("path").datum(cdfData).attr("d", line).attr("fill", "none").attr("stroke", getColor(1, theme)).attr("stroke-width", 2).attr("stroke-dasharray", "6,3");
  }
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("x");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Density");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function getDistributionFunctions(params) {
  const p = params.params;
  switch (params.distribution) {
    case "normal": {
      const mu = p["mean"] ?? 0, sigma = p["sd"] ?? 1;
      return {
        xMin: mu - 4 * sigma,
        xMax: mu + 4 * sigma,
        pdf: (x) => Math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * Math.sqrt(2 * Math.PI)),
        cdf: (x) => normalCDF((x - mu) / sigma)
      };
    }
    case "t": {
      const df = p["df"] ?? 5;
      return {
        xMin: -5,
        xMax: 5,
        pdf: (x) => {
          const c = Math.exp(lgamma((df + 1) / 2) - lgamma(df / 2)) / Math.sqrt(df * Math.PI);
          return c * Math.pow(1 + x * x / df, -(df + 1) / 2);
        },
        cdf: (x) => tDistCDFLocal(x, df)
      };
    }
    case "chi-square": {
      const df = p["df"] ?? 5;
      return {
        xMin: 0,
        xMax: df + 5 * Math.sqrt(2 * df),
        pdf: (x) => x <= 0 ? 0 : Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - df / 2 * Math.log(2) - lgamma(df / 2)),
        cdf: (x) => chiSqCDFLocal(x, df)
      };
    }
    case "uniform": {
      const a = p["min"] ?? 0, b = p["max"] ?? 1;
      return {
        xMin: a - 0.2 * (b - a),
        xMax: b + 0.2 * (b - a),
        pdf: (x) => x >= a && x <= b ? 1 / (b - a) : 0,
        cdf: (x) => x < a ? 0 : x > b ? 1 : (x - a) / (b - a)
      };
    }
    case "exponential": {
      const rate = p["rate"] ?? 1;
      return {
        xMin: 0,
        xMax: 5 / rate,
        pdf: (x) => x < 0 ? 0 : rate * Math.exp(-rate * x),
        cdf: (x) => x < 0 ? 0 : 1 - Math.exp(-rate * x)
      };
    }
    case "F": {
      const df1 = p["df1"] ?? 2, df2 = p["df2"] ?? 10;
      return {
        xMin: 0,
        xMax: 6,
        pdf: (x) => {
          if (x <= 0) return 0;
          const num = Math.pow(df1 * x, df1) * Math.pow(df2, df2);
          const den = Math.pow(df1 * x + df2, df1 + df2);
          return Math.sqrt(num / den) / (x * Math.exp(logBetaLocal(df1 / 2, df2 / 2)));
        },
        cdf: (x) => {
          if (x <= 0) return 0;
          return incompleteBetaLocal(df1 * x / (df1 * x + df2), df1 / 2, df2 / 2);
        }
      };
    }
    default:
      throw new Error(`Unknown distribution: ${params.distribution}`);
  }
}
function formatDistributionParams(params) {
  const entries = Object.entries(params.params).map(([k, v]) => `${k} = ${v}`).join(", ");
  return `${params.distribution}(${entries})`;
}
function lgamma(z) {
  const c = [0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9984369578019572e-21, 15056327351493116e-23];
  const x = z - 1;
  let sum = c[0];
  for (let i = 1; i < 9; i++) sum += (c[i] ?? 0) / (x + i);
  const t = x + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum);
}
function tDistCDFLocal(t, df) {
  const x = df / (df + t * t);
  const p = incompleteBetaLocal(x, df / 2, 0.5) / 2;
  return t >= 0 ? 1 - p : p;
}
function chiSqCDFLocal(x, df) {
  if (x <= 0) return 0;
  return incompleteGammaLocal(df / 2, x / 2);
}
function incompleteGammaLocal(a, x) {
  if (x < a + 1) {
    let term = 1 / a, sum = term;
    for (let n = 1; n < 200; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 3e-7) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
  }
  let f = x + 1 - a;
  const FPMIN = 1e-30;
  if (Math.abs(f) < FPMIN) f = FPMIN;
  let C = f, D = 0;
  for (let i = 1; i <= 200; i++) {
    const an = -i * (i - a), bn = x + 2 * i + 1 - a;
    D = bn + an * D;
    if (Math.abs(D) < FPMIN) D = FPMIN;
    C = bn + an / C;
    if (Math.abs(C) < FPMIN) C = FPMIN;
    D = 1 / D;
    const delta = C * D;
    f *= delta;
    if (Math.abs(delta - 1) < 3e-7) break;
  }
  return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) / f;
}
function logBetaLocal(a, b) {
  return lgamma(a) + lgamma(b) - lgamma(a + b);
}
function incompleteBetaLocal(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  if (x > (a + 1) / (a + b + 2)) return 1 - incompleteBetaLocal(1 - x, b, a);
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - logBetaLocal(a, b)) / a;
  const FPMIN = 1e-30;
  let f = 1, C = 1, D = 1 - (a + b) * x / (a + 1);
  if (Math.abs(D) < FPMIN) D = FPMIN;
  D = 1 / D;
  f = D;
  for (let m = 1; m <= 200; m++) {
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
    if (Math.abs(delta - 1) < 3e-7) break;
  }
  return front * f;
}

// src/viz/plots/density.ts
function renderDensity(container, data, config = {}) {
  import("d3").then((d3) => renderDensityD3(d3, container, data, config));
}
function silvermanBW(values) {
  const n = values.length;
  if (n < 2) return 1;
  const m = values.reduce((s, v) => s + v, 0) / n;
  const std = Math.sqrt(values.reduce((s, v) => s + (v - m) ** 2, 0) / (n - 1));
  return 1.06 * std * Math.pow(n, -0.2);
}
function kdeGaussian(values, xPoints, bw) {
  return xPoints.map((x) => {
    const density = values.reduce(
      (s, xi) => s + Math.exp(-0.5 * ((x - xi) / bw) ** 2) / (Math.sqrt(2 * Math.PI) * bw),
      0
    ) / values.length;
    return [x, density];
  });
}
function renderDensityD3(d3, container, data, config) {
  if (data.series.length === 0 || data.series.every((s) => s.values.length === 0)) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 420;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const showRug = config.showRug ?? true;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Density Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allValues = data.series.flatMap((s) => [...s.values]);
  const [xMin, xMax] = d3.extent(allValues);
  const xPad = (xMax - xMin) * 0.1;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const xTicks = xScale.ticks(100);
  const seriesDensities = data.series.map((s) => {
    const bw = config.bandwidth ?? silvermanBW(s.values);
    return kdeGaussian(s.values, xTicks, bw);
  });
  const maxDensity = Math.max(...seriesDensities.flatMap((pts) => pts.map((p) => p[1])));
  const rugHeight = showRug ? height - 12 : height;
  const yScale = d3.scaleLinear().domain([0, maxDensity * 1.1]).range([rugHeight, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  seriesDensities.forEach((pts, i) => {
    const color = getColor(i, theme);
    const lineGen = d3.line().x((d) => xScale(d[0])).y((d) => yScale(d[1])).curve(d3.curveCatmullRom);
    const areaGen = d3.area().x((d) => xScale(d[0])).y0(yScale(0)).y1((d) => yScale(d[1])).curve(d3.curveCatmullRom);
    g.append("path").datum(pts).attr("d", areaGen).attr("fill", color).attr("opacity", theme.ciOpacity);
    g.append("path").datum(pts).attr("d", lineGen).attr("fill", "none").attr("stroke", color).attr("stroke-width", 2);
  });
  if (showRug) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme);
      const rugY = height - 10;
      s.values.forEach((v) => {
        g.append("line").attr("x1", xScale(v)).attr("x2", xScale(v)).attr("y1", rugY).attr("y2", rugY + 6).attr("stroke", color).attr("stroke-width", 1).attr("opacity", 0.5).on("mouseover", (event) => {
          showTooltip(
            event,
            [formatTooltipRow("Series", s.label), formatTooltipRow("Value", v.toFixed(3))].join(""),
            theme
          );
        }).on("mouseout", hideTooltip);
      });
    });
  }
  g.append("g").attr("transform", `translate(0,${rugHeight})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Density");
  if (config.showLegend !== false && data.series.length > 1) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme);
      const lx = width - 120;
      const ly = i * 18 + 8;
      g.append("line").attr("x1", lx).attr("x2", lx + 20).attr("y1", ly).attr("y2", ly).attr("stroke", color).attr("stroke-width", 2);
      g.append("text").attr("x", lx + 24).attr("y", ly + 4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(s.label);
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/boxplot.ts
function renderBoxplot(container, data, config = {}) {
  import("d3").then((d3) => renderBoxplotD3(d3, container, data, config));
}
function computeBoxStats(values) {
  const q1 = quantile(values, 0.25);
  const med = quantile(values, 0.5);
  const q3 = quantile(values, 0.75);
  const iqr = q3 - q1;
  const fence_lo = q1 - 1.5 * iqr;
  const fence_hi = q3 + 1.5 * iqr;
  const whisker_lo = Math.min(...values.filter((v) => v >= fence_lo));
  const whisker_hi = Math.max(...values.filter((v) => v <= fence_hi));
  const outliers = values.filter((v) => v < fence_lo || v > fence_hi);
  return { q1, med, q3, whisker_lo, whisker_hi, outliers };
}
function renderBoxplotD3(d3, container, data, config) {
  if (data.groups.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 440;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Box Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const labels = data.groups.map((gr) => gr.label);
  const allValues = data.groups.flatMap((gr) => [...gr.values]);
  const [yMin, yMax] = d3.extent(allValues);
  const yPad = (yMax - yMin) * 0.1;
  const xScale = d3.scaleBand().domain(labels).range([0, width]).padding(0.3);
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  data.groups.forEach((gr, gi) => {
    if (gr.values.length === 0) return;
    const color = getColor(gi, theme);
    const bx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2;
    const boxW = xScale.bandwidth() * 0.5;
    const { q1, med, q3, whisker_lo, whisker_hi, outliers } = computeBoxStats(gr.values);
    g.append("line").attr("x1", bx).attr("x2", bx).attr("y1", yScale(whisker_lo)).attr("y2", yScale(q1)).attr("stroke", color).attr("stroke-width", 1.5).attr("stroke-dasharray", "3,2");
    g.append("line").attr("x1", bx).attr("x2", bx).attr("y1", yScale(q3)).attr("y2", yScale(whisker_hi)).attr("stroke", color).attr("stroke-width", 1.5).attr("stroke-dasharray", "3,2");
    const capW = boxW * 0.4;
    g.append("line").attr("x1", bx - capW / 2).attr("x2", bx + capW / 2).attr("y1", yScale(whisker_lo)).attr("y2", yScale(whisker_lo)).attr("stroke", color).attr("stroke-width", 1.5);
    g.append("line").attr("x1", bx - capW / 2).attr("x2", bx + capW / 2).attr("y1", yScale(whisker_hi)).attr("y2", yScale(whisker_hi)).attr("stroke", color).attr("stroke-width", 1.5);
    g.append("rect").attr("x", bx - boxW / 2).attr("width", boxW).attr("y", yScale(q3)).attr("height", Math.max(1, yScale(q1) - yScale(q3))).attr("fill", color).attr("opacity", theme.violinOpacity).attr("stroke", color).attr("stroke-width", 2).attr("rx", 2).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("Group", gr.label),
        formatTooltipRow("Median", med.toFixed(3)),
        formatTooltipRow("Q1", q1.toFixed(3)),
        formatTooltipRow("Q3", q3.toFixed(3))
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    if (config.showMedian !== false) {
      g.append("line").attr("x1", bx - boxW / 2).attr("x2", bx + boxW / 2).attr("y1", yScale(med)).attr("y2", yScale(med)).attr("stroke", theme.background).attr("stroke-width", 2.5);
      g.append("text").attr("x", bx + boxW / 2 + 4).attr("y", yScale(med) + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(med.toFixed(2));
    }
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length;
      const my = yScale(groupMean);
      const ds = 5;
      g.append("polygon").attr("points", `${bx},${my - ds} ${bx + ds},${my} ${bx},${my + ds} ${bx - ds},${my}`).attr("fill", "white").attr("stroke", color).attr("stroke-width", 1.5);
      g.append("text").attr("x", bx + ds + 4).attr("y", my + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(groupMean.toFixed(2));
    }
    if (config.showOutliers !== false) {
      outliers.forEach((v) => {
        g.append("circle").attr("cx", bx).attr("cy", yScale(v)).attr("r", 3.5).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.5).on("mouseover", (event) => {
          showTooltip(event, [
            formatTooltipRow("Group", gr.label),
            formatTooltipRow("Outlier", v.toFixed(3))
          ].join(""), theme);
        }).on("mouseout", hideTooltip);
      });
    }
    if (config.showN !== false) {
      g.append("text").attr("x", bx).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text(`n = ${gr.values.length}`);
    }
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 60).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/lollipop.ts
function renderLollipop(container, data, config = {}) {
  import("d3").then((d3) => renderLollipopD3(d3, container, data, config));
}
function renderLollipopD3(d3, container, data, config) {
  if (data.labels.length === 0 || data.values.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const itemCount = data.labels.length;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? Math.max(200, itemCount * 36 + theme.marginTop + theme.marginBottom);
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft + 40 };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Lollipop Chart", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  let indices = data.labels.map((_, i) => i);
  if (config.sorted) {
    indices = [...indices].sort((a, b) => (data.values[b] ?? 0) - (data.values[a] ?? 0));
  }
  const orderedLabels = indices.map((i) => data.labels[i] ?? "");
  const orderedValues = indices.map((i) => data.values[i] ?? 0);
  const [vMin, vMax] = d3.extent(orderedValues);
  const xMin = Math.min(0, vMin);
  const xMax = Math.max(0, vMax);
  const xPad = (xMax - xMin) * 0.1;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleBand().domain(orderedLabels).range([0, height]).padding(0.35);
  g.selectAll(".grid").data(xScale.ticks(5)).join("line").attr("class", "grid").attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d)).attr("y1", 0).attr("y2", height).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", xScale(0)).attr("x2", xScale(0)).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-width", 1.5);
  orderedLabels.forEach((label, i) => {
    const val = orderedValues[i] ?? 0;
    const cy = (yScale(label) ?? 0) + yScale.bandwidth() / 2;
    const color = getColor(0, theme);
    g.append("line").attr("x1", xScale(0)).attr("x2", xScale(val)).attr("y1", cy).attr("y2", cy).attr("stroke", color).attr("stroke-width", 2).attr("opacity", 0.7);
    g.append("circle").attr("cx", xScale(val)).attr("cy", cy).attr("r", 6).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("Label", label),
        formatTooltipRow("Value", val.toFixed(3))
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Value");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -margin.left + 12).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/dot-plot.ts
function renderDotPlot(container, data, config = {}) {
  import("d3").then((d3) => renderDotPlotD3(d3, container, data, config));
}
function renderDotPlotD3(d3, container, data, config) {
  if (data.labels.length === 0 || data.values.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const itemCount = data.labels.length;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? Math.max(200, itemCount * 36 + theme.marginTop + theme.marginBottom);
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft + 40 };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const isPaired = data.group2 != null && data.group2.length === data.labels.length;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Dot Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allValues = [...data.values];
  if (isPaired) allValues.push(...data.group2);
  const [vMin, vMax] = d3.extent(allValues);
  const xPad = (vMax - vMin) * 0.1;
  const xScale = d3.scaleLinear().domain([vMin - xPad, vMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleBand().domain([...data.labels]).range([0, height]).padding(0.35);
  g.selectAll(".grid").data(xScale.ticks(5)).join("line").attr("class", "grid").attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d)).attr("y1", 0).attr("y2", height).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const color1 = getColor(0, theme);
  const color2 = getColor(1, theme);
  data.labels.forEach((label, i) => {
    const v1 = data.values[i] ?? 0;
    const v2 = isPaired ? data.group2[i] ?? 0 : null;
    const cy = (yScale(label) ?? 0) + yScale.bandwidth() / 2;
    if (isPaired && v2 !== null) {
      g.append("line").attr("x1", xScale(v1)).attr("x2", xScale(v2)).attr("y1", cy).attr("y2", cy).attr("stroke", theme.gridLine).attr("stroke-width", 2);
    }
    g.append("circle").attr("cx", xScale(v1)).attr("cy", cy).attr("r", 6).attr("fill", color1).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("Label", label),
        formatTooltipRow(data.group1Label ?? "Value", v1.toFixed(3))
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    if (isPaired && v2 !== null) {
      g.append("circle").attr("cx", xScale(v2)).attr("cy", cy).attr("r", 6).attr("fill", color2).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Label", label),
          formatTooltipRow(data.group2Label ?? "Group 2", v2.toFixed(3))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    }
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Value");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -margin.left + 12).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "");
  if (isPaired) {
    const legendX = width - 140;
    const label1 = data.group1Label ?? "Group 1";
    const label2 = data.group2Label ?? "Group 2";
    g.append("circle").attr("cx", legendX + 6).attr("cy", 8).attr("r", 5).attr("fill", color1);
    g.append("text").attr("x", legendX + 16).attr("y", 12).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label1);
    g.append("circle").attr("cx", legendX + 6).attr("cy", 26).attr("r", 5).attr("fill", color2);
    g.append("text").attr("x", legendX + 16).attr("y", 30).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label2);
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/grouped-bar.ts
function renderGroupedBar(container, data, config = {}) {
  import("d3").then((d3) => renderGroupedBarD3(d3, container, data, config));
}
function renderGroupedBarD3(d3, container, data, config) {
  if (data.categories.length === 0 || data.series.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const chartType = config.type ?? "grouped";
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 440;
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? (chartType === "stacked" ? "Stacked Bar Chart" : "Grouped Bar Chart"),
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const categories = [...data.categories];
  const seriesLabels = data.series.map((s) => s.label);
  const rows = categories.map((cat, ci) => {
    const obj = { cat };
    data.series.forEach((s) => {
      obj[s.label] = s.values[ci] ?? 0;
    });
    return obj;
  });
  let yMax;
  if (chartType === "stacked") {
    yMax = Math.max(...rows.map(
      (row) => data.series.reduce((sum, s) => sum + row[s.label], 0)
    ));
  } else {
    yMax = Math.max(...data.series.flatMap((s) => [...s.values].map((v) => Math.abs(v))));
  }
  const xOuter = d3.scaleBand().domain(categories).range([0, width]).padding(0.2);
  const xInner = d3.scaleBand().domain(seriesLabels).range([0, xOuter.bandwidth()]).padding(0.05);
  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  if (chartType === "stacked") {
    const stackGen = d3.stack().keys(seriesLabels).value((row, key) => row[key] ?? 0);
    const stackedData = stackGen(rows);
    stackedData.forEach((layer, si) => {
      const color = getColor(si, theme);
      layer.forEach((seg, ci) => {
        const cat = categories[ci] ?? "";
        const bx = xOuter(cat) ?? 0;
        const bw = xOuter.bandwidth();
        g.append("rect").attr("x", bx).attr("width", bw).attr("y", yScale(seg[1])).attr("height", Math.max(0, yScale(seg[0]) - yScale(seg[1]))).attr("fill", color).attr("opacity", theme.violinOpacity).attr("stroke", theme.background).attr("stroke-width", 0.5).on("mouseover", (event) => {
          showTooltip(event, [
            formatTooltipRow("Category", cat),
            formatTooltipRow("Series", seriesLabels[si] ?? ""),
            formatTooltipRow("Value", (seg[1] - seg[0]).toFixed(3))
          ].join(""), theme);
        }).on("mouseout", hideTooltip);
      });
    });
  } else {
    categories.forEach((cat) => {
      const bxOuter = xOuter(cat) ?? 0;
      data.series.forEach((s, si) => {
        const color = getColor(si, theme);
        const val = s.values[categories.indexOf(cat)] ?? 0;
        const bx = bxOuter + (xInner(s.label) ?? 0);
        const bw = xInner.bandwidth();
        g.append("rect").attr("x", bx).attr("width", bw).attr("y", yScale(Math.max(0, val))).attr("height", Math.abs(yScale(0) - yScale(val))).attr("fill", color).attr("opacity", theme.violinOpacity).attr("stroke", color).attr("stroke-width", 0.5).attr("rx", 2).on("mouseover", (event) => {
          showTooltip(event, [
            formatTooltipRow("Category", cat),
            formatTooltipRow("Series", s.label),
            formatTooltipRow("Value", val.toFixed(3))
          ].join(""), theme);
        }).on("mouseout", hideTooltip);
      });
    });
  }
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xOuter)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  data.series.forEach((s, i) => {
    const color = getColor(i, theme);
    const lx = width - 130;
    const ly = i * 18 + 8;
    g.append("rect").attr("x", lx).attr("y", ly - 9).attr("width", 14).attr("height", 10).attr("fill", color).attr("opacity", theme.violinOpacity).attr("rx", 2);
    g.append("text").attr("x", lx + 18).attr("y", ly).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(s.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/line-chart.ts
function renderLineChart(container, data, config = {}) {
  import("d3").then((d3) => renderLineChartD3(d3, container, data, config));
}
function renderLineChartD3(d3, container, data, config) {
  if (data.series.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 420;
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const showArea = config.showArea ?? false;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Line Chart", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allX = data.series.flatMap((s) => [...s.x]);
  const allY = data.series.flatMap((s) => [...s.y]);
  const [xMin, xMax] = d3.extent(allX);
  const [yMin, yMax] = d3.extent(allY);
  const xPad = (xMax - xMin) * 0.05;
  const yPad = (yMax - yMin) * 0.1;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleLinear().domain([Math.min(0, yMin - yPad), yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  data.series.forEach((s, si) => {
    if (s.x.length === 0 || s.y.length === 0) return;
    const color = getColor(si, theme);
    const pts = s.x.map((x, i) => [x, s.y[i] ?? 0]).sort((a, b) => a[0] - b[0]);
    if (showArea) {
      const areaGen = d3.area().x((d) => xScale(d[0])).y0(yScale(0)).y1((d) => yScale(d[1])).curve(d3.curveCatmullRom);
      g.append("path").datum(pts).attr("d", areaGen).attr("fill", color).attr("opacity", theme.ciOpacity);
    }
    const lineGen = d3.line().x((d) => xScale(d[0])).y((d) => yScale(d[1])).curve(d3.curveCatmullRom);
    g.append("path").datum(pts).attr("d", lineGen).attr("fill", "none").attr("stroke", color).attr("stroke-width", 2);
    pts.forEach(([x, y]) => {
      g.append("circle").attr("cx", xScale(x)).attr("cy", yScale(y)).attr("r", 4).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Series", s.label),
          formatTooltipRow("x", x.toFixed(3)),
          formatTooltipRow("y", y.toFixed(3))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    });
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "x");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "y");
  if (data.series.length > 1) {
    data.series.forEach((s, i) => {
      const color = getColor(i, theme);
      const lx = width - 130;
      const ly = i * 18 + 8;
      g.append("line").attr("x1", lx).attr("x2", lx + 20).attr("y1", ly).attr("y2", ly).attr("stroke", color).attr("stroke-width", 2);
      g.append("circle").attr("cx", lx + 10).attr("cy", ly).attr("r", 3).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1);
      g.append("text").attr("x", lx + 24).attr("y", ly + 4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(s.label);
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/bubble-chart.ts
function renderBubbleChart(container, data, config = {}) {
  import("d3").then((d3) => renderBubbleChartD3(d3, container, data, config));
}
function renderBubbleChartD3(d3, container, data, config) {
  if (data.points.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 460;
  const margin = { top: theme.marginTop, right: theme.marginRight + 20, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Bubble Chart", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const xVals = data.points.map((p) => p.x);
  const yVals = data.points.map((p) => p.y);
  const rVals = data.points.map((p) => p.r);
  const [xMin, xMax] = d3.extent(xVals);
  const [yMin, yMax] = d3.extent(yVals);
  const [rMin, rMax] = d3.extent(rVals);
  const xPad = (xMax - xMin) * 0.1;
  const yPad = (yMax - yMin) * 0.1;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  const rScale = rMax === rMin ? (_) => 12 : d3.scaleSqrt().domain([rMin, rMax]).range([4, 30]);
  const groups = [...new Set(data.points.map((p) => p.group ?? "__default__"))];
  const groupIndex = (grp) => groups.indexOf(grp);
  const hasGroups = !(groups.length === 1 && groups[0] === "__default__");
  g.selectAll(".gridH").data(yScale.ticks(5)).join("line").attr("class", "gridH").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const sorted = [...data.points].sort((a, b) => b.r - a.r);
  sorted.forEach((pt) => {
    const grp = pt.group ?? "__default__";
    const color = getColor(groupIndex(grp), theme);
    const cx = xScale(pt.x);
    const cy = yScale(pt.y);
    const radius = rScale(pt.r);
    g.append("circle").attr("cx", cx).attr("cy", cy).attr("r", radius).attr("fill", color).attr("opacity", theme.pointOpacity).attr("stroke", color).attr("stroke-width", 1.5).attr("stroke-opacity", 0.8).on("mouseover", (event) => {
      const rows = [
        formatTooltipRow("x", pt.x.toFixed(3)),
        formatTooltipRow("y", pt.y.toFixed(3)),
        formatTooltipRow("size", pt.r.toFixed(3))
      ];
      if (pt.label) rows.unshift(formatTooltipRow("Label", pt.label));
      if (hasGroups && pt.group) rows.push(formatTooltipRow("Group", pt.group));
      showTooltip(event, rows.join(""), theme);
    }).on("mouseout", hideTooltip);
    if (pt.label && radius > 16) {
      g.append("text").attr("x", cx).attr("y", cy + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.background).attr("pointer-events", "none").text(pt.label.slice(0, 6));
    }
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "x");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "y");
  if (hasGroups) {
    groups.forEach((grp, i) => {
      const color = getColor(i, theme);
      const lx = width - 130;
      const ly = i * 18 + 8;
      g.append("circle").attr("cx", lx + 6).attr("cy", ly).attr("r", 5).attr("fill", color).attr("opacity", theme.pointOpacity);
      g.append("text").attr("x", lx + 16).attr("y", ly + 4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(grp);
    });
  }
  if (rMax > rMin) {
    const sizeLegendVals = [rMin, (rMin + rMax) / 2, rMax];
    const slx = 16;
    let sly = height - 70;
    g.append("text").attr("x", slx).attr("y", sly - 6).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text("Size:");
    sizeLegendVals.forEach((rv) => {
      const rad = rScale(rv);
      g.append("circle").attr("cx", slx + 16).attr("cy", sly + rad).attr("r", rad).attr("fill", "none").attr("stroke", theme.textMuted).attr("stroke-width", 1);
      g.append("text").attr("x", slx + 36).attr("y", sly + rad + 4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).text(rv.toFixed(1));
      sly += rad * 2 + 10;
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/pareto.ts
function renderPareto(container, data, config = {}) {
  import("d3").then((d3) => renderParetoD3(d3, container, data, config));
}
function renderParetoD3(d3, container, data, config) {
  if (data.labels.length === 0 || data.values.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: 64, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Pareto Chart", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const indices = Array.from({ length: data.labels.length }, (_, i) => i).sort((a, b) => (data.values[b] ?? 0) - (data.values[a] ?? 0));
  const total = indices.reduce((s, i) => s + (data.values[i] ?? 0), 0);
  let running = 0;
  const points = indices.map((origIdx, sortPos) => {
    const val = data.values[origIdx] ?? 0;
    const lbl = data.labels[origIdx] ?? "";
    running += val;
    return {
      label: lbl,
      value: val,
      cumPct: running / total * 100,
      index: sortPos
    };
  });
  const sortedLabels = points.map((p) => p.label);
  const xScale = d3.scaleBand().domain(sortedLabels).range([0, width]).padding(0.25);
  const yScaleBar = d3.scaleLinear().domain([0, (points[0]?.value ?? 0) * 1.1]).range([height, 0]).nice();
  const yScalePct = d3.scaleLinear().domain([0, 100]).range([height, 0]);
  g.selectAll(".grid").data(yScaleBar.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScaleBar(d)).attr("y2", (d) => yScaleBar(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  const barColor = getColor(0, theme);
  const lineColor = getColor(1, theme);
  g.selectAll(".bar").data(points).join("rect").attr("class", "bar").attr("x", (d) => xScale(d.label) ?? 0).attr("y", (d) => yScaleBar(d.value)).attr("width", xScale.bandwidth()).attr("height", (d) => height - yScaleBar(d.value)).attr("fill", barColor).attr("opacity", 0.85).on("mouseover", function(event, d) {
    showTooltip(event, [
      formatTooltipRow("Category", d.label),
      formatTooltipRow("Value", d.value.toFixed(2)),
      formatTooltipRow("Cumulative %", d.cumPct.toFixed(1) + "%")
    ].join(""), theme);
  }).on("mouseout", hideTooltip);
  const lineGen = d3.line().x((d) => (xScale(d.label) ?? 0) + xScale.bandwidth() / 2).y((d) => yScalePct(d.cumPct)).curve(d3.curveMonotoneX);
  g.append("path").datum(points).attr("d", lineGen).attr("fill", "none").attr("stroke", lineColor).attr("stroke-width", 2.5);
  g.selectAll(".cum-dot").data(points).join("circle").attr("class", "cum-dot").attr("cx", (d) => (xScale(d.label) ?? 0) + xScale.bandwidth() / 2).attr("cy", (d) => yScalePct(d.cumPct)).attr("r", 4).attr("fill", lineColor).attr("stroke", theme.background).attr("stroke-width", 1.5);
  g.append("line").attr("x1", 0).attr("x2", width).attr("y1", yScalePct(80)).attr("y2", yScalePct(80)).attr("stroke", "#D55E00").attr("stroke-width", 1.5).attr("stroke-dasharray", "6,4");
  g.append("text").attr("x", width + 4).attr("y", yScalePct(80) + 4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", "#D55E00").text("80%");
  const leftAxis = g.append("g").call(d3.axisLeft(yScaleBar).ticks(6));
  leftAxis.selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  const rightAxis = g.append("g").attr("transform", `translate(${width},0)`).call(d3.axisRight(yScalePct).ticks(5).tickFormat((d) => `${d}%`));
  rightAxis.selectAll("text").attr("fill", lineColor).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  const bottomAxis = g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale));
  bottomAxis.selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("text-anchor", "end").attr("transform", "rotate(-35)");
  g.append("text").attr("x", width / 2).attr("y", height + 60).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Frequency");
  g.append("text").attr("transform", `translate(${width + 56},${height / 2}) rotate(90)`).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", lineColor).text("Cumulative %");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/funnel.ts
function renderFunnel(container, data, config = {}) {
  import("d3").then((d3) => renderFunnelD3(d3, container, data, config));
}
function renderFunnelD3(_d3, container, data, config) {
  if (data.stages.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svgNS = "http://www.w3.org/2000/svg";
  const svgEl = document.createElementNS(svgNS, "svg");
  svgEl.setAttribute("width", String(W));
  svgEl.setAttribute("height", String(H));
  svgEl.style.background = theme.background;
  container.appendChild(svgEl);
  const titleEl = document.createElementNS(svgNS, "text");
  titleEl.setAttribute("x", String(W / 2));
  titleEl.setAttribute("y", "20");
  titleEl.setAttribute("text-anchor", "middle");
  titleEl.setAttribute("font-family", theme.fontFamily);
  titleEl.setAttribute("font-size", String(theme.fontSizeTitle));
  titleEl.setAttribute("font-weight", "600");
  titleEl.setAttribute("fill", theme.text);
  titleEl.textContent = config.title ?? "Funnel Chart";
  svgEl.appendChild(titleEl);
  const subtitleEl = document.createElementNS(svgNS, "text");
  subtitleEl.setAttribute("x", String(W / 2));
  subtitleEl.setAttribute("y", "38");
  subtitleEl.setAttribute("text-anchor", "middle");
  subtitleEl.setAttribute("font-family", theme.fontFamilyMono ?? theme.fontFamily);
  subtitleEl.setAttribute("font-size", String(theme.fontSizeSmall));
  subtitleEl.setAttribute("fill", theme.textMuted);
  subtitleEl.textContent = data.testResult?.formatted ?? "";
  svgEl.appendChild(subtitleEl);
  const maxValue = Math.max(...data.stages.map((s) => s.value));
  const n = data.stages.length;
  const segH = height / n;
  const gap = 2;
  const labelColW = 110;
  const valueColW = 110;
  const funnelW = width - labelColW - valueColW;
  const funnelOffsetX = margin.left + labelColW;
  const offsetY = margin.top;
  data.stages.forEach((stage, i) => {
    const color = getColor(i, theme);
    const topW = stage.value / maxValue * funnelW;
    const nextStage = i + 1 < data.stages.length ? data.stages[i + 1] : void 0;
    const bottomW = nextStage != null ? nextStage.value / maxValue * funnelW : topW * 0.6;
    const topY = offsetY + i * segH;
    const bottomY = offsetY + (i + 1) * segH - gap;
    const topLeft = funnelOffsetX + (funnelW - topW) / 2;
    const topRight = funnelOffsetX + (funnelW + topW) / 2;
    const botLeft = funnelOffsetX + (funnelW - bottomW) / 2;
    const botRight = funnelOffsetX + (funnelW + bottomW) / 2;
    const pathD = `M ${topLeft},${topY} L ${topRight},${topY} L ${botRight},${bottomY} L ${botLeft},${bottomY} Z`;
    const pathEl = document.createElementNS(svgNS, "path");
    pathEl.setAttribute("d", pathD);
    pathEl.setAttribute("fill", color);
    pathEl.setAttribute("opacity", "0.85");
    pathEl.addEventListener("mouseover", (event) => {
      const pct = (stage.value / maxValue * 100).toFixed(1);
      const prevStage2 = i > 0 ? data.stages[i - 1] : void 0;
      const drop = prevStage2 != null ? ((prevStage2.value - stage.value) / prevStage2.value * 100).toFixed(1) : "\u2014";
      showTooltip(event, [
        formatTooltipRow("Stage", stage.label),
        formatTooltipRow("Value", stage.value.toFixed(0)),
        formatTooltipRow("% of Max", pct + "%"),
        formatTooltipRow("Drop from prev", drop === "\u2014" ? drop : drop + "%")
      ].join(""), theme);
    });
    pathEl.addEventListener("mouseout", hideTooltip);
    svgEl.appendChild(pathEl);
    const midY = topY + segH / 2;
    const labelEl = document.createElementNS(svgNS, "text");
    labelEl.setAttribute("x", String(margin.left + labelColW - 8));
    labelEl.setAttribute("y", String(midY));
    labelEl.setAttribute("text-anchor", "end");
    labelEl.setAttribute("dominant-baseline", "middle");
    labelEl.setAttribute("font-family", theme.fontFamily);
    labelEl.setAttribute("font-size", String(theme.fontSize));
    labelEl.setAttribute("fill", theme.text);
    labelEl.textContent = stage.label;
    svgEl.appendChild(labelEl);
    const pctOfMax = (stage.value / maxValue * 100).toFixed(0);
    const prevStage = i > 0 ? data.stages[i - 1] : void 0;
    const dropStr = prevStage != null ? `\u2212${((prevStage.value - stage.value) / prevStage.value * 100).toFixed(0)}%` : "";
    const valEl = document.createElementNS(svgNS, "text");
    valEl.setAttribute("x", String(funnelOffsetX + funnelW + 8));
    valEl.setAttribute("y", String(midY - (dropStr ? 7 : 0)));
    valEl.setAttribute("dominant-baseline", "middle");
    valEl.setAttribute("font-family", theme.fontFamilyMono ?? theme.fontFamily);
    valEl.setAttribute("font-size", String(theme.fontSize));
    valEl.setAttribute("fill", theme.text);
    valEl.textContent = `${stage.value.toLocaleString()} (${pctOfMax}%)`;
    svgEl.appendChild(valEl);
    if (dropStr) {
      const dropEl = document.createElementNS(svgNS, "text");
      dropEl.setAttribute("x", String(funnelOffsetX + funnelW + 8));
      dropEl.setAttribute("y", String(midY + 9));
      dropEl.setAttribute("dominant-baseline", "middle");
      dropEl.setAttribute("font-family", theme.fontFamily);
      dropEl.setAttribute("font-size", String(theme.fontSizeSmall));
      dropEl.setAttribute("fill", "#D55E00");
      dropEl.textContent = dropStr;
      svgEl.appendChild(dropEl);
    }
  });
  if (config.caption) {
    const capEl = document.createElementNS(svgNS, "text");
    capEl.setAttribute("x", String(W / 2));
    capEl.setAttribute("y", String(H - 6));
    capEl.setAttribute("text-anchor", "middle");
    capEl.setAttribute("font-family", theme.fontFamily);
    capEl.setAttribute("font-size", String(theme.fontSizeSmall - 1));
    capEl.setAttribute("fill", theme.textMuted);
    capEl.style.fontStyle = "italic";
    capEl.textContent = config.caption;
    svgEl.appendChild(capEl);
  }
}

// src/viz/plots/pie-chart.ts
function renderPieChart(container, data, config = {}) {
  import("d3").then((d3) => renderPieChartD3(d3, container, data, config));
}
function renderPieChartD3(d3, container, data, config) {
  if (data.slices.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Pie Chart", data.testResult?.formatted ?? "", W, theme);
  const legendH = Math.ceil(data.slices.length / 3) * 22;
  const plotH = H - margin.top - margin.bottom - legendH;
  const cx = W / 2;
  const cy = margin.top + plotH / 2;
  const outerR = Math.min(W / 2, plotH / 2) * 0.72;
  const innerR = config.donut ?? false ? outerR * 0.55 : 0;
  const total = data.slices.reduce((s, sl) => s + sl.value, 0);
  const pieGen = d3.pie().value((d) => d.value).sort(null);
  const arcGen = d3.arc().innerRadius(innerR).outerRadius(outerR);
  const labelArc = d3.arc().innerRadius(outerR * 1.15).outerRadius(outerR * 1.15);
  const polyArc = d3.arc().innerRadius(outerR * 0.85).outerRadius(outerR * 0.85);
  const arcs = pieGen(data.slices);
  const g = svg.append("g").attr("transform", `translate(${cx},${cy})`);
  arcs.forEach((arc, i) => {
    const color = getColor(i, theme);
    const pct = (arc.data.value / total * 100).toFixed(1);
    g.append("path").attr("d", arcGen(arc) ?? "").attr("fill", color).attr("opacity", 0.88).attr("stroke", theme.background).attr("stroke-width", 2).on("mouseover", function(event) {
      showTooltip(event, [
        formatTooltipRow("Label", arc.data.label),
        formatTooltipRow("Value", arc.data.value.toFixed(2)),
        formatTooltipRow("Percentage", pct + "%")
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    const midAngle = (arc.startAngle + arc.endAngle) / 2;
    const isRight = midAngle < Math.PI;
    const labelPt = labelArc.centroid(arc);
    const polyPt = polyArc.centroid(arc);
    const endX = isRight ? outerR * 1.35 : -outerR * 1.35;
    if (arc.data.value / total > 0.03) {
      g.append("polyline").attr("points", [polyPt, labelPt, [endX, labelPt[1]]].map((p) => p.join(",")).join(" ")).attr("fill", "none").attr("stroke", theme.textMuted).attr("stroke-width", 1);
      const displayStr = config.showPercentages ?? true ? `${arc.data.label} (${pct}%)` : arc.data.label;
      g.append("text").attr("x", isRight ? endX + 4 : endX - 4).attr("y", labelPt[1]).attr("dominant-baseline", "middle").attr("text-anchor", isRight ? "start" : "end").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(displayStr);
    }
  });
  if (config.donut ?? false) {
    g.append("text").attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize + 2).attr("font-weight", "600").attr("fill", theme.text).text(`n = ${total.toFixed(0)}`);
  }
  const legendY = margin.top + plotH + 12;
  const colW = W / 3;
  data.slices.forEach((sl, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const lx = col * colW + 16;
    const ly = legendY + row * 22;
    svg.append("rect").attr("x", lx).attr("y", ly).attr("width", 12).attr("height", 12).attr("fill", getColor(i, theme)).attr("rx", 2);
    svg.append("text").attr("x", lx + 16).attr("y", ly + 9).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(sl.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/area-chart.ts
function renderAreaChart(container, data, config = {}) {
  import("d3").then((d3) => renderAreaChartD3(d3, container, data, config));
}
function renderAreaChartD3(d3, container, data, config) {
  if (data.series.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Area Chart", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allX = data.series.flatMap((s) => [...s.x]);
  const xMin = Math.min(...allX);
  const xMax = Math.max(...allX);
  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]);
  const stacked = config.stacked ?? false;
  let yMax;
  if (stacked) {
    const allXSorted = [...new Set(allX)].sort((a, b) => a - b);
    const rowSums = allXSorted.map(
      (xVal) => data.series.reduce((sum, s) => {
        const idx = s.x.indexOf(xVal);
        return sum + (idx >= 0 ? s.y[idx] ?? 0 : 0);
      }, 0)
    );
    yMax = Math.max(...rowSums);
  } else {
    yMax = Math.max(...data.series.flatMap((s) => [...s.y]));
  }
  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  if (stacked) {
    const allXSorted = [...new Set(allX)].sort((a, b) => a - b);
    const baselines = /* @__PURE__ */ new Map();
    allXSorted.forEach((xVal) => baselines.set(xVal, 0));
    data.series.forEach((series, si) => {
      const color = getColor(si, theme);
      const points = series.x.map((xVal, xi) => ({
        xVal,
        top: (baselines.get(xVal) ?? 0) + (series.y[xi] ?? 0),
        bot: baselines.get(xVal) ?? 0
      }));
      const areaFn = d3.area().x((d) => xScale(d.xVal)).y0((d) => yScale(d.bot)).y1((d) => yScale(d.top)).curve(d3.curveMonotoneX);
      g.append("path").datum(points).attr("d", areaFn).attr("fill", color).attr("opacity", 0.75).attr("stroke", color).attr("stroke-width", 1.5);
      series.x.forEach((xVal, xi) => {
        baselines.set(xVal, (baselines.get(xVal) ?? 0) + (series.y[xi] ?? 0));
      });
    });
  } else {
    data.series.forEach((series, si) => {
      const color = getColor(si, theme);
      const points = series.x.map((xVal, xi) => ({ xVal, yVal: series.y[xi] })).filter((p) => p.yVal !== void 0);
      const areaFn = d3.area().x((d) => xScale(d.xVal)).y0(yScale(0)).y1((d) => yScale(d.yVal)).curve(d3.curveMonotoneX);
      const lineFn = d3.line().x((d) => xScale(d.xVal)).y((d) => yScale(d.yVal)).curve(d3.curveMonotoneX);
      g.append("path").datum(points).attr("d", areaFn).attr("fill", color).attr("opacity", 0.22);
      g.append("path").datum(points).attr("d", lineFn).attr("fill", "none").attr("stroke", color).attr("stroke-width", 2.5);
      g.selectAll(`.dot-s${si}`).data(points).join("circle").attr("class", `dot-s${si}`).attr("cx", (d) => xScale(d.xVal)).attr("cy", (d) => yScale(d.yVal)).attr("r", 3.5).attr("fill", color).attr("opacity", 0.7).on("mouseover", function(event, d) {
        showTooltip(event, [
          formatTooltipRow("Series", series.label),
          formatTooltipRow("x", d.xVal.toFixed(3)),
          formatTooltipRow("y", d.yVal.toFixed(3))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    });
  }
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(8)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "x");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  data.series.forEach((s, i) => {
    const lx = margin.left + i * 130;
    const ly = H - 18;
    svg.append("rect").attr("x", lx).attr("y", ly - 9).attr("width", 12).attr("height", 12).attr("fill", getColor(i, theme)).attr("rx", 2);
    svg.append("text").attr("x", lx + 16).attr("y", ly).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(s.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/forest-plot.ts
function renderForestPlot(container, data, config = {}) {
  import("d3").then((d3) => renderForestPlotD3(d3, container, data, config));
}
function renderForestPlotD3(d3, container, data, config) {
  if (data.studies.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const totalRows = data.studies.length + (data.pooled != null ? 2 : 0);
  const rowH = 28;
  const extraH = theme.marginTop + theme.marginBottom + 20;
  const H = config.height ?? Math.max(totalRows * rowH + extraH, 320);
  const W = config.width ?? Math.max(container.clientWidth || 700, 500);
  const labelColW = 160;
  const estColW = 140;
  const marginLeft = labelColW + 8;
  const marginRight = estColW + 16;
  const margin = { top: theme.marginTop, right: marginRight, bottom: theme.marginBottom, left: marginLeft };
  const width = W - margin.left - margin.right;
  const height = totalRows * rowH;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H + margin.top + margin.bottom).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Forest Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const allLow = data.studies.map((s) => s.ciLow).concat(data.pooled != null ? [data.pooled.ciLow] : []);
  const allHigh = data.studies.map((s) => s.ciHigh).concat(data.pooled != null ? [data.pooled.ciHigh] : []);
  const xMin = Math.min(...allLow);
  const xMax = Math.max(...allHigh);
  const xPad = (xMax - xMin) * 0.12;
  const xScale = d3.scaleLinear().domain([xMin - xPad, xMax + xPad]).range([0, width]).nice();
  const maxWeight = data.studies.reduce((m, s) => Math.max(m, s.weight ?? 1), 0);
  const [domainMin, domainMax] = xScale.domain();
  if (domainMin <= 0 && domainMax >= 0) {
    g.append("line").attr("x1", xScale(0)).attr("x2", xScale(0)).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-width", 1).attr("stroke-dasharray", "4,3");
  }
  svg.append("text").attr("x", margin.left - 8).attr("y", margin.top - 6).attr("text-anchor", "end").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.textMuted).text("Study");
  svg.append("text").attr("x", W - marginRight + 8).attr("y", margin.top - 6).attr("text-anchor", "start").attr("font-family", theme.fontFamilyMono ?? theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.textMuted).text("Estimate [95% CI]");
  const studyColor = getColor(0, theme);
  data.studies.forEach((study, i) => {
    const cy = i * rowH + rowH / 2;
    g.append("line").attr("x1", xScale(study.ciLow)).attr("x2", xScale(study.ciHigh)).attr("y1", cy).attr("y2", cy).attr("stroke", studyColor).attr("stroke-width", 1.5);
    const capH = 5;
    [study.ciLow, study.ciHigh].forEach((xv) => {
      g.append("line").attr("x1", xScale(xv)).attr("x2", xScale(xv)).attr("y1", cy - capH).attr("y2", cy + capH).attr("stroke", studyColor).attr("stroke-width", 1.5);
    });
    const sqHalf = study.weight != null && maxWeight > 0 ? 4 + study.weight / maxWeight * 6 : 6;
    g.append("rect").attr("x", xScale(study.estimate) - sqHalf / 2).attr("y", cy - sqHalf / 2).attr("width", sqHalf).attr("height", sqHalf).attr("fill", studyColor).on("mouseover", function(event) {
      showTooltip(event, [
        formatTooltipRow("Study", study.label),
        formatTooltipRow("Estimate", study.estimate.toFixed(3)),
        formatTooltipRow("95% CI", `[${study.ciLow.toFixed(3)}, ${study.ciHigh.toFixed(3)}]`),
        ...study.weight != null ? [formatTooltipRow("Weight", study.weight.toFixed(2))] : []
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    svg.append("text").attr("x", margin.left - 8).attr("y", margin.top + cy).attr("text-anchor", "end").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(study.label);
    const estStr = `${study.estimate.toFixed(2)} [${study.ciLow.toFixed(2)}, ${study.ciHigh.toFixed(2)}]`;
    svg.append("text").attr("x", W - marginRight + 8).attr("y", margin.top + cy).attr("dominant-baseline", "middle").attr("font-family", theme.fontFamilyMono ?? theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(estStr);
  });
  if (data.pooled != null) {
    const pooled = data.pooled;
    const sepY = data.studies.length * rowH;
    g.append("line").attr("x1", 0).attr("x2", width).attr("y1", sepY).attr("y2", sepY).attr("stroke", theme.axisLine).attr("stroke-width", 1);
    const pooledY = sepY + rowH * 1;
    const diamondColor = getColor(4, theme);
    const dLeft = xScale(pooled.ciLow);
    const dRight = xScale(pooled.ciHigh);
    const dMid = xScale(pooled.estimate);
    const dH = 9;
    const diamondPath = [
      `M ${dMid},${pooledY - dH}`,
      `L ${dRight},${pooledY}`,
      `L ${dMid},${pooledY + dH}`,
      `L ${dLeft},${pooledY}`,
      "Z"
    ].join(" ");
    g.append("path").attr("d", diamondPath).attr("fill", diamondColor).attr("opacity", 0.88).attr("stroke", diamondColor).attr("stroke-width", 1).on("mouseover", function(event) {
      showTooltip(event, [
        formatTooltipRow("Pooled estimate", pooled.estimate.toFixed(3)),
        formatTooltipRow("95% CI", `[${pooled.ciLow.toFixed(3)}, ${pooled.ciHigh.toFixed(3)}]`)
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    svg.append("text").attr("x", margin.left - 8).attr("y", margin.top + pooledY).attr("text-anchor", "end").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("font-weight", "600").attr("fill", theme.text).text("Pooled");
    const pEstStr = `${pooled.estimate.toFixed(2)} [${pooled.ciLow.toFixed(2)}, ${pooled.ciHigh.toFixed(2)}]`;
    svg.append("text").attr("x", W - marginRight + 8).attr("y", margin.top + pooledY).attr("dominant-baseline", "middle").attr("font-family", theme.fontFamilyMono ?? theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.text).text(pEstStr);
  }
  g.append("g").attr("transform", `translate(0,${height + 8})`).call(d3.axisBottom(xScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "Effect Size");
  if (config.caption) addCaption(svg, config.caption, W, H + margin.top + margin.bottom, theme);
}

// src/viz/plots/roc-curve.ts
function renderROCCurve(container, data, config = {}) {
  import("d3").then((d3) => renderROCCurveD3(d3, container, data, config));
}
function renderROCCurveD3(d3, container, data, config) {
  if (data.fpr.length === 0 || data.tpr.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 520, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  const aucStr = `AUC = ${data.auc.toFixed(4)}`;
  const subtitle = data.testResult?.formatted ? `${data.testResult.formatted}   ${aucStr}` : aucStr;
  addSubtitle(svg, config.title ?? "ROC Curve", subtitle, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width]);
  const yScale = d3.scaleLinear().domain([0, 1]).range([height, 0]);
  g.selectAll(".grid-h").data(yScale.ticks(5)).join("line").attr("class", "grid-h").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", xScale(0)).attr("y1", yScale(0)).attr("x2", xScale(1)).attr("y2", yScale(1)).attr("stroke", theme.axisLine).attr("stroke-width", 1.5).attr("stroke-dasharray", "6,4");
  const pairs = data.fpr.map((fprVal, i) => ({ fprVal, tprVal: data.tpr[i] ?? 0 })).sort((a, b) => a.fprVal - b.fprVal);
  const curveColor = getColor(0, theme);
  const areaFn = d3.area().x((d) => xScale(d.fprVal)).y0(yScale(0)).y1((d) => yScale(d.tprVal)).curve(d3.curveLinear);
  g.append("path").datum(pairs).attr("d", areaFn).attr("fill", curveColor).attr("opacity", theme.ciOpacity);
  const lineFn = d3.line().x((d) => xScale(d.fprVal)).y((d) => yScale(d.tprVal)).curve(d3.curveLinear);
  g.append("path").datum(pairs).attr("d", lineFn).attr("fill", "none").attr("stroke", curveColor).attr("stroke-width", 2.5);
  g.selectAll(".roc-dot").data(pairs).join("circle").attr("class", "roc-dot").attr("cx", (d) => xScale(d.fprVal)).attr("cy", (d) => yScale(d.tprVal)).attr("r", 3).attr("fill", curveColor).attr("opacity", 0.6).on("mouseover", function(event, d) {
    showTooltip(event, [
      formatTooltipRow("FPR", d.fprVal.toFixed(4)),
      formatTooltipRow("TPR", d.tprVal.toFixed(4)),
      formatTooltipRow("Specificity", (1 - d.fprVal).toFixed(4)),
      formatTooltipRow("Sensitivity", d.tprVal.toFixed(4))
    ].join(""), theme);
  }).on("mouseout", hideTooltip);
  g.append("rect").attr("x", 8).attr("y", 8).attr("width", 130).attr("height", 28).attr("fill", theme.surface).attr("stroke", theme.gridLine).attr("rx", 4).attr("opacity", 0.9);
  g.append("text").attr("x", 16).attr("y", 26).attr("font-family", theme.fontFamilyMono ?? theme.fontFamily).attr("font-size", theme.fontSize + 1).attr("font-weight", "600").attr("fill", curveColor).text(`AUC = ${data.auc.toFixed(4)}`);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "False Positive Rate");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "True Positive Rate");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/strip-plot.ts
function renderStripPlot(container, data, config = {}) {
  import("d3").then((d3) => renderStripPlotD3(d3, container, data, config));
}
function renderStripPlotD3(d3, container, data, config) {
  if (data.groups.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Strip Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const labels = data.groups.map((gr) => gr.label);
  const allValues = data.groups.flatMap((gr) => [...gr.values]);
  const yMinRaw = Math.min(...allValues);
  const yMaxRaw = Math.max(...allValues);
  const yPad = (yMaxRaw - yMinRaw) * 0.1;
  const xScale = d3.scaleBand().domain(labels).range([0, width]).padding(0.3);
  const yScale = d3.scaleLinear().domain([yMinRaw - yPad, yMaxRaw + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  data.groups.forEach((gr, gi) => {
    if (gr.values.length === 0) return;
    const color = getColor(gi, theme);
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2;
    const jitterWidth = xScale.bandwidth() * 0.38;
    const seed = gi * 77773;
    const mean2 = gr.values.reduce((s, v) => s + v, 0) / gr.values.length;
    gr.values.forEach((v, vi) => {
      const jitter = (pseudoRnd2(seed + vi * 13 + gi) - 0.5) * jitterWidth;
      g.append("circle").attr("cx", cx + jitter).attr("cy", yScale(v)).attr("r", 3.5).attr("fill", color).attr("opacity", theme.pointOpacity).on("mouseover", function(event) {
        showTooltip(event, [
          formatTooltipRow("Group", gr.label),
          formatTooltipRow("Value", v.toFixed(4)),
          formatTooltipRow("Group mean", mean2.toFixed(4))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    });
    const lineHalfW = xScale.bandwidth() * 0.32;
    g.append("line").attr("x1", cx - lineHalfW).attr("x2", cx + lineHalfW).attr("y1", yScale(mean2)).attr("y2", yScale(mean2)).attr("stroke", color).attr("stroke-width", 2.5).attr("opacity", 0.9);
    if (config.showN !== false) {
      g.append("text").attr("x", cx).attr("y", height + 52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text(`n = ${gr.values.length}`);
    }
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function pseudoRnd2(seed) {
  const x = Math.sin(seed) * 1e4;
  return x - Math.floor(x);
}

// src/viz/plots/swarm-plot.ts
function renderSwarmPlot(container, data, config = {}) {
  import("d3").then((d3) => renderSwarmPlotD3(d3, container, data, config));
}
function renderSwarmPlotD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const r = config.pointRadius ?? 4;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Beeswarm Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const labels = data.groups.map((gr) => gr.label);
  const xScale = d3.scaleBand().domain(labels).range([0, width]).padding(0.3);
  const allValues = data.groups.flatMap((gr) => [...gr.values]);
  const [yMin, yMax] = d3.extent(allValues);
  const yPad = (yMax - yMin) * 0.08 || 1;
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(6)).join("line").attr("class", "grid").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).call((ax) => ax.select(".domain").attr("stroke", theme.axisLine)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(6)).call((ax) => ax.select(".domain").attr("stroke", theme.axisLine)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel ?? "");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel ?? "Value");
  data.groups.forEach((gr, gi) => {
    const color = getColor(gi, theme);
    const cx = (xScale(gr.label) ?? 0) + xScale.bandwidth() / 2;
    const maxHalfWidth = xScale.bandwidth() / 2 - r;
    const positions = beeswarm([...gr.values], yScale, r, maxHalfWidth);
    positions.forEach(({ value, xOffset }) => {
      g.append("circle").attr("cx", cx + xOffset).attr("cy", yScale(value)).attr("r", r).attr("fill", color).attr("opacity", theme.pointOpacity).attr("stroke", theme.background).attr("stroke-width", 0.5).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Group", gr.label),
          formatTooltipRow("Value", value.toFixed(4))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    });
    if (config.showMean) {
      const groupMean = gr.values.reduce((s, v) => s + v, 0) / gr.values.length;
      const my = yScale(groupMean);
      const ds = 5;
      g.append("polygon").attr("points", `${cx},${my - ds} ${cx + ds},${my} ${cx},${my + ds} ${cx - ds},${my}`).attr("fill", "white").attr("stroke", color).attr("stroke-width", 1.5);
      g.append("text").attr("x", cx + ds + 4).attr("y", my + 3.5).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textAnnotation).text(groupMean.toFixed(2));
    }
    if (config.showN !== false) {
      svg.append("text").attr("x", margin.left + cx).attr("y", H - theme.marginBottom + 28).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text(`n = ${gr.values.length}`);
    }
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function beeswarm(values, yScale, r, maxHalf) {
  const sorted = values.slice().sort((a, b) => a - b);
  const placed = [];
  const diameter = r * 2 + 0.5;
  const result = sorted.map((value) => {
    const y = yScale(value);
    let xOffset = 0;
    let placed_flag = false;
    for (let attempt = 0; attempt <= Math.ceil(maxHalf / diameter) + 1; attempt++) {
      const candidates = attempt === 0 ? [0] : [attempt * diameter, -attempt * diameter];
      for (const xOff of candidates) {
        if (Math.abs(xOff) > maxHalf) continue;
        const overlaps = placed.some((p) => {
          const dx = xOff - p.xOffset;
          const dy = y - p.y;
          return Math.sqrt(dx * dx + dy * dy) < diameter;
        });
        if (!overlaps) {
          xOffset = xOff;
          placed_flag = true;
          break;
        }
      }
      if (placed_flag) break;
    }
    placed.push({ y, xOffset });
    return { value, xOffset };
  });
  return result;
}

// src/viz/plots/mosaic-plot.ts
var RESID_BINS = [
  { lo: 4, hi: Infinity, fill: "#2166ac", label: "> 4", borderStyle: "thick" },
  { lo: 2, hi: 4, fill: "#74add1", label: "2 : 4", borderStyle: "dashed" },
  { lo: 0, hi: 2, fill: "#d1e5f0", label: "0 : 2", borderStyle: "none" },
  { lo: -2, hi: 0, fill: "#fddbc7", label: "\u22122 : 0", borderStyle: "none" },
  { lo: -4, hi: -2, fill: "#d6604d", label: "\u22124 : \u22122", borderStyle: "dashed" },
  { lo: -Infinity, hi: -4, fill: "#b2182b", label: "< \u22124", borderStyle: "thick" }
];
function binForResidual(r) {
  return RESID_BINS.find((b) => r >= b.lo && r < b.hi) ?? RESID_BINS[2];
}
function textColor(fill) {
  const dark = ["#2166ac", "#b2182b", "#d6604d"];
  return dark.includes(fill) ? "#ffffff" : "#212529";
}
function renderMosaicPlot(container, data, config = {}) {
  import("d3").then((d3) => renderMosaicPlotD3(d3, container, data, config));
}
function renderMosaicPlotD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const legendW = 110;
  const W = config.width ?? Math.max(container.clientWidth || 620, 420);
  const H = config.height ?? 480;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight + legendW,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Mosaic Plot", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const nRows = data.table.length;
  const nCols = data.table[0]?.length ?? 0;
  if (nRows === 0 || nCols === 0) return;
  const colSums = Array.from(
    { length: nCols },
    (_, j) => data.table.reduce((s, row) => s + (row[j] ?? 0), 0)
  );
  const grandTotal = colSums.reduce((a, b) => a + b, 0);
  const rowSums = data.table.map((row) => row.reduce((a, b) => a + b, 0));
  const residuals = data.table.map(
    (row, i) => row.map((v, j) => {
      const e = rowSums[i] * colSums[j] / grandTotal;
      return e > 0 ? (v - e) / Math.sqrt(e) : 0;
    })
  );
  const gap = 3;
  const colWidths = colSums.map((s) => s / grandTotal * width);
  const colX = [];
  colWidths.reduce((acc, w, i) => {
    colX[i] = acc;
    return acc + w;
  }, 0);
  data.table.forEach((row, ri) => {
    row.forEach((count, ci) => {
      const colSum = colSums[ci] ?? 1;
      const cellH = colSum > 0 ? count / colSum * height : 0;
      const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[ci] ?? 0), 0);
      const cellY = colSum > 0 ? rowsAbove / colSum * height : 0;
      const cx = colX[ci] ?? 0;
      const cw = colWidths[ci] ?? 0;
      const resid = residuals[ri][ci];
      const bin = binForResidual(resid);
      const rx = cx + gap / 2;
      const ry = cellY + gap / 2;
      const rw = Math.max(cw - gap, 0);
      const rh = Math.max(cellH - gap, 0);
      g.append("rect").attr("x", rx).attr("y", ry).attr("width", rw).attr("height", rh).attr("fill", bin.fill).attr("stroke", "none");
      if (bin.borderStyle !== "none") {
        const isThick = bin.borderStyle === "thick";
        g.append("rect").attr("x", rx).attr("y", ry).attr("width", rw).attr("height", rh).attr("fill", "none").attr("stroke", isThick ? resid > 0 ? "#08519c" : "#a50f15" : resid > 0 ? "#4292c6" : "#cb181d").attr("stroke-width", isThick ? 2.5 : 1.5).attr("stroke-dasharray", isThick ? "none" : "5,2");
      }
      g.append("rect").attr("x", rx).attr("y", ry).attr("width", rw).attr("height", rh).attr("fill", "transparent").on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Row", data.rowLabels[ri] ?? `Row ${ri}`),
          formatTooltipRow("Column", data.colLabels[ci] ?? `Col ${ci}`),
          formatTooltipRow("Count", count),
          formatTooltipRow("Residual", resid.toFixed(2)),
          formatTooltipRow("Bin", bin.label)
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
      if (rw > 30 && rh > 16) {
        g.append("text").attr("x", rx + rw / 2).attr("y", ry + rh / 2 + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall, rh * 0.38, rw * 0.22)).attr("fill", textColor(bin.fill)).attr("pointer-events", "none").text(String(count));
      }
    });
  });
  colX.forEach((x, ci) => {
    const cw = colWidths[ci] ?? 0;
    g.append("text").attr("x", x + cw / 2).attr("y", height + 18).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(data.colLabels[ci] ?? `Col ${ci}`);
  });
  data.table.forEach((row, ri) => {
    const colSum = colSums[0] ?? 1;
    const rowsAbove = data.table.slice(0, ri).reduce((s, r2) => s + (r2[0] ?? 0), 0);
    const cellY = colSum > 0 ? rowsAbove / colSum * height : 0;
    const cellH = colSum > 0 ? (row[0] ?? 0) / colSum * height : 0;
    g.append("text").attr("x", -8).attr("y", cellY + cellH / 2 + 4).attr("text-anchor", "end").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(data.rowLabels[ri] ?? `Row ${ri}`);
  });
  if (config.xLabel) {
    g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel);
  }
  if (config.yLabel) {
    g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel);
  }
  const lgX = width + 20;
  const lgY0 = height / 2 - RESID_BINS.length * 22 / 2;
  const swatchS = 14;
  g.append("text").attr("x", lgX).attr("y", lgY0 - 14).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("font-weight", "600").attr("fill", theme.text).text("Standardised");
  g.append("text").attr("x", lgX).attr("y", lgY0 - 2).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("font-weight", "600").attr("fill", theme.text).text("residual");
  RESID_BINS.forEach((bin, bi) => {
    const ly = lgY0 + bi * 22;
    g.append("rect").attr("x", lgX).attr("y", ly).attr("width", swatchS).attr("height", swatchS).attr("fill", bin.fill).attr("rx", 2);
    if (bin.borderStyle !== "none") {
      g.append("rect").attr("x", lgX).attr("y", ly).attr("width", swatchS).attr("height", swatchS).attr("fill", "none").attr("rx", 2).attr("stroke", bin.fill === "#2166ac" || bin.fill === "#74add1" ? "#4292c6" : "#cb181d").attr("stroke-width", bin.borderStyle === "thick" ? 2 : 1.5).attr("stroke-dasharray", bin.borderStyle === "thick" ? "none" : "4,2");
    }
    g.append("text").attr("x", lgX + swatchS + 6).attr("y", ly + swatchS - 3).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(bin.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/pair-plot.ts
function renderPairPlot(container, data, config = {}) {
  import("d3").then((d3) => renderPairPlotD3(d3, container, data, config));
}
function renderPairPlotD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? Math.max(W, 480);
  const outerMarginTop = 52;
  const outerMarginBottom = config.caption ? 28 : 16;
  const outerMarginSide = 16;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Scatter Matrix", data.testResult?.formatted ?? "", W, theme);
  const n = data.labels.length;
  if (n < 2) return;
  const gridW = W - outerMarginSide * 2;
  const gridH = H - outerMarginTop - outerMarginBottom;
  const cellW = gridW / n;
  const cellH = gridH / n;
  const cellPad = 6;
  const scales = data.data.map((vals) => {
    const [lo, hi] = d3.extent(vals);
    const pad = (hi - lo) * 0.08 || 0.5;
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([cellPad, cellW - cellPad]).nice();
  });
  const scalesY = data.data.map((vals) => {
    const [lo, hi] = d3.extent(vals);
    const pad = (hi - lo) * 0.08 || 0.5;
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([cellH - cellPad, cellPad]).nice();
  });
  const gridG = svg.append("g").attr("transform", `translate(${outerMarginSide},${outerMarginTop})`);
  for (let row = 0; row < n; row++) {
    for (let col = 0; col < n; col++) {
      const cellX = col * cellW;
      const cellY = row * cellH;
      const cellG = gridG.append("g").attr("transform", `translate(${cellX},${cellY})`);
      cellG.append("rect").attr("width", cellW).attr("height", cellH).attr("fill", row === col ? theme.surface : theme.background).attr("stroke", theme.gridLine).attr("stroke-width", 0.5);
      if (row === col) {
        drawHistogramCell(d3, cellG, data.data[row], scales[col], cellW, cellH, cellPad, getColor(row, theme), theme);
        cellG.append("text").attr("x", cellW / 2).attr("y", cellPad + 10).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall, cellW / 8)).attr("font-weight", "600").attr("fill", theme.text).text(data.labels[row] ?? "");
      } else {
        const xVals = data.data[col];
        const yVals = data.data[row];
        const xSc = scales[col];
        const ySc = scalesY[row];
        const n_obs = Math.min(xVals.length, yVals.length);
        for (let i = 0; i < n_obs; i++) {
          const xv = xVals[i];
          const yv = yVals[i];
          cellG.append("circle").attr("cx", xSc(xv)).attr("cy", ySc(yv)).attr("r", Math.min(2.5, cellW / 40)).attr("fill", getColor(col, theme)).attr("opacity", theme.pointOpacity).on("mouseover", (event) => {
            showTooltip(event, [
              formatTooltipRow(data.labels[col] ?? `Var ${col}`, xv.toFixed(3)),
              formatTooltipRow(data.labels[row] ?? `Var ${row}`, yv.toFixed(3))
            ].join(""), theme);
          }).on("mouseout", hideTooltip);
        }
        const r = pearsonR2(xVals.slice(0, n_obs), yVals.slice(0, n_obs));
        cellG.append("text").attr("x", cellW - cellPad - 2).attr("y", cellH - cellPad - 2).attr("text-anchor", "end").attr("font-family", theme.fontFamilyMono).attr("font-size", Math.min(theme.fontSizeSmall - 1, cellW / 9)).attr("fill", theme.textMuted).text(`r=${r.toFixed(2)}`);
      }
    }
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function drawHistogramCell(d3, g, values, xScale, _cellW, cellH, pad, color, _theme) {
  const bins = d3.bin().domain(xScale.domain()).thresholds(8)([...values]);
  const maxCount = Math.max(...bins.map((b) => b.length));
  const yScale = d3.scaleLinear().domain([0, maxCount]).range([cellH - pad - 14, pad + 14]);
  bins.forEach((bin) => {
    const x0 = xScale(bin.x0 ?? 0);
    const x1 = xScale(bin.x1 ?? 0);
    const bw = Math.max(x1 - x0 - 1, 1);
    g.append("rect").attr("x", x0).attr("y", yScale(bin.length)).attr("width", bw).attr("height", Math.max(cellH - pad - 14 - yScale(bin.length), 0)).attr("fill", color).attr("opacity", 0.6);
  });
}
function pearsonR2(xs, ys) {
  const n = xs.length;
  if (n < 2) return 0;
  const mx = xs.reduce((a, b) => a + b, 0) / n;
  const my = ys.reduce((a, b) => a + b, 0) / n;
  let num = 0, sdx = 0, sdy = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx;
    const dy = ys[i] - my;
    num += dx * dy;
    sdx += dx * dx;
    sdy += dy * dy;
  }
  const denom = Math.sqrt(sdx * sdy);
  return denom === 0 ? 0 : num / denom;
}

// src/viz/plots/radar-chart.ts
function renderRadarChart(container, data, config = {}) {
  import("d3").then((d3) => renderRadarChartD3(d3, container, data, config));
}
function renderRadarChartD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const levels = config.levels ?? 5;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Radar Chart", data.testResult?.formatted ?? "", W, theme);
  const nAxes = data.axes.length;
  if (nAxes < 3 || data.series.length === 0) return;
  const legendWidth = 120;
  const chartCentreX = (W - legendWidth) / 2;
  const chartCentreY = (H - theme.marginTop - 32) / 2 + theme.marginTop + 8;
  const radius = Math.min((W - legendWidth) / 2, (H - theme.marginTop - 48) / 2) * 0.72;
  const g = svg.append("g");
  const axisMin = Array.from(
    { length: nAxes },
    (_, ai) => Math.min(...data.series.map((s) => s.values[ai] ?? 0))
  );
  const axisMax = Array.from(
    { length: nAxes },
    (_, ai) => Math.max(...data.series.map((s) => s.values[ai] ?? 0))
  );
  const normalize = (val, ai) => {
    const lo = axisMin[ai];
    const hi = axisMax[ai];
    return hi === lo ? 0.5 : (val - lo) / (hi - lo);
  };
  const angleOf = (i) => i / nAxes * 2 * Math.PI - Math.PI / 2;
  const pt = (t, ai) => {
    const a = angleOf(ai);
    return [
      chartCentreX + radius * t * Math.cos(a),
      chartCentreY + radius * t * Math.sin(a)
    ];
  };
  for (let lvl = 1; lvl <= levels; lvl++) {
    const t = lvl / levels;
    const ringPts = Array.from({ length: nAxes }, (_, ai) => pt(t, ai));
    g.append("polygon").attr("points", ringPts.map(([x, y]) => `${x},${y}`).join(" ")).attr("fill", "none").attr("stroke", theme.gridLine).attr("stroke-width", 1);
  }
  Array.from({ length: nAxes }, (_, ai) => {
    const [x2, y2] = pt(1, ai);
    g.append("line").attr("x1", chartCentreX).attr("y1", chartCentreY).attr("x2", x2).attr("y2", y2).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  });
  data.axes.forEach((label, ai) => {
    const angle = angleOf(ai);
    const labelDist = radius * 1.15;
    const lx = chartCentreX + labelDist * Math.cos(angle);
    const ly = chartCentreY + labelDist * Math.sin(angle);
    let anchor = "middle";
    if (Math.cos(angle) > 0.1) anchor = "start";
    else if (Math.cos(angle) < -0.1) anchor = "end";
    g.append("text").attr("x", lx).attr("y", ly + 4).attr("text-anchor", anchor).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label);
  });
  data.series.forEach((series, si) => {
    const color = getColor(si, theme);
    const polyPts = data.axes.map((_, ai) => {
      const t = normalize(series.values[ai] ?? 0, ai);
      return pt(t, ai);
    });
    g.append("polygon").attr("points", polyPts.map(([x, y]) => `${x},${y}`).join(" ")).attr("fill", color).attr("fill-opacity", theme.ciOpacity + 0.05).attr("stroke", color).attr("stroke-width", 2);
    polyPts.forEach(([px, py], ai) => {
      g.append("circle").attr("cx", px).attr("cy", py).attr("r", 4).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Series", series.label),
          formatTooltipRow("Axis", data.axes[ai] ?? `Axis ${ai}`),
          formatTooltipRow("Value", (series.values[ai] ?? 0).toFixed(3))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    });
  });
  const lgX = W - legendWidth + 8;
  data.series.forEach((series, si) => {
    const color = getColor(si, theme);
    const lgY = theme.marginTop + si * 22;
    g.append("rect").attr("x", lgX).attr("y", lgY).attr("width", 12).attr("height", 12).attr("fill", color).attr("opacity", 0.8).attr("rx", 2);
    g.append("text").attr("x", lgX + 16).attr("y", lgY + 10).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(series.label.length > 14 ? series.label.slice(0, 13) + "\u2026" : series.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/parallel-coords.ts
function renderParallelCoords(container, data, config = {}) {
  import("d3").then((d3) => renderParallelCoordsD3(d3, container, data, config));
}
function renderParallelCoordsD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 700, 400);
  const H = config.height ?? 480;
  const margin = {
    top: theme.marginTop + 16,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Parallel Coordinates", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const nAxes = data.axes.length;
  if (nAxes < 2 || data.rows.length === 0) return;
  const axisX = (i) => i / (nAxes - 1) * width;
  const yScales = data.axes.map((_, ai) => {
    const vals = data.rows.map((row) => row[ai] ?? 0);
    const [lo, hi] = d3.extent(vals);
    const pad = (hi - lo) * 0.05 || 0.5;
    return d3.scaleLinear().domain([lo - pad, hi + pad]).range([height, 0]).nice();
  });
  data.rows.forEach((row, ri) => {
    const groupIdx = data.groups?.[ri] ?? 0;
    const color = getColor(groupIdx, theme);
    const pathData = data.axes.map((_, ai) => {
      const val = row[ai] ?? 0;
      return [axisX(ai), yScales[ai](val)];
    });
    g.append("path").attr("d", lineThrough(pathData)).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.2).attr("opacity", Math.max(0.1, Math.min(0.55, 30 / data.rows.length))).on("mouseover", (event) => {
      const tooltipRows = data.axes.map(
        (label, ai) => formatTooltipRow(label, (row[ai] ?? 0).toFixed(3))
      );
      showTooltip(event, tooltipRows.join(""), theme);
    }).on("mouseout", hideTooltip);
  });
  data.axes.forEach((label, ai) => {
    const x = axisX(ai);
    const ySc = yScales[ai];
    g.append("line").attr("x1", x).attr("x2", x).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-width", 1.5);
    const ticks = ySc.ticks(5);
    ticks.forEach((tick) => {
      g.append("line").attr("x1", x - 4).attr("x2", x + 4).attr("y1", ySc(tick)).attr("y2", ySc(tick)).attr("stroke", theme.axisLine).attr("stroke-width", 1);
      g.append("text").attr("x", x - 6).attr("y", ySc(tick) + 4).attr("text-anchor", "end").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).text(tick.toFixed(1));
    });
    g.append("text").attr("x", x).attr("y", -10).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.text).text(label);
  });
  if (data.groups) {
    const uniqueGroups = Array.from(new Set(data.groups)).sort((a, b) => a - b);
    uniqueGroups.forEach((grp, i) => {
      const lgX = 0;
      const lgY = height + 32 + i * 18;
      g.append("rect").attr("x", lgX).attr("y", lgY).attr("width", 10).attr("height", 10).attr("fill", getColor(grp, theme)).attr("rx", 2);
      g.append("text").attr("x", lgX + 14).attr("y", lgY + 9).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(`Group ${grp}`);
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function lineThrough(pts) {
  if (pts.length === 0) return "";
  return pts.map(([x, y], i) => `${i === 0 ? "M" : "L"}${x.toFixed(2)},${y.toFixed(2)}`).join(" ");
}

// src/viz/plots/treemap.ts
function renderTreemap(container, data, config = {}) {
  import("d3").then((d3) => renderTreemapD3(d3, container, data, config));
}
function renderTreemapD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? 480;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Treemap", data.testResult?.formatted ?? "", W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  if (data.children.length === 0) return;
  const uniqueGroups = Array.from(new Set(data.children.map((c) => c.group ?? ""))).filter(Boolean);
  const groupColorMap = new Map(uniqueGroups.map((grp, i) => [grp, getColor(i, theme)]));
  const total = data.children.reduce((s, c) => s + c.value, 0);
  const root = d3.hierarchy(
    { children: data.children }
  ).sum((d) => {
    const node = d;
    return node.value ?? 0;
  }).sort((a, b) => (b.value ?? 0) - (a.value ?? 0));
  const treemapLayout = d3.treemap().size([width, height]).paddingInner(2).paddingOuter(0).round(true);
  treemapLayout(root);
  const leaves = root.leaves();
  leaves.forEach((leaf, li) => {
    const childData = leaf.data;
    const x0 = leaf.x0;
    const y0 = leaf.y0;
    const x1 = leaf.x1;
    const y1 = leaf.y1;
    const cellW = x1 - x0;
    const cellH = y1 - y0;
    const color = childData.group && groupColorMap.has(childData.group) ? groupColorMap.get(childData.group) : getColor(li, theme);
    const pct = total > 0 ? childData.value / total * 100 : 0;
    const cell = g.append("g").attr("transform", `translate(${x0},${y0})`);
    cell.append("rect").attr("width", cellW).attr("height", cellH).attr("fill", color).attr("opacity", 0.75).attr("stroke", theme.background).attr("stroke-width", 2).attr("rx", 2).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("Label", childData.label),
        formatTooltipRow("Value", childData.value.toFixed(2)),
        formatTooltipRow("Share", `${pct.toFixed(1)}%`),
        ...childData.group ? [formatTooltipRow("Group", childData.group)] : []
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    if (cellW > 30 && cellH > 16) {
      const maxChars = Math.max(3, Math.floor(cellW / 7));
      const label = childData.label.length > maxChars ? childData.label.slice(0, maxChars - 1) + "\u2026" : childData.label;
      cell.append("text").attr("x", 5).attr("y", 14).attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall, cellH / 3)).attr("font-weight", "600").attr("fill", "#fff").attr("pointer-events", "none").text(label);
      if (cellH > 28) {
        cell.append("text").attr("x", 5).attr("y", 26).attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall - 2, cellH / 4)).attr("fill", "rgba(255,255,255,0.75)").attr("pointer-events", "none").text(`${pct.toFixed(1)}%`);
      }
    }
  });
  if (uniqueGroups.length > 0) {
    uniqueGroups.forEach((grp, i) => {
      const lgX = 0;
      const lgY = height + 24 + i * 18;
      if (lgY + 12 > H - margin.top - 4) return;
      g.append("rect").attr("x", lgX).attr("y", lgY).attr("width", 12).attr("height", 12).attr("fill", groupColorMap.get(grp) ?? theme.gridLine).attr("rx", 2);
      g.append("text").attr("x", lgX + 16).attr("y", lgY + 10).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(grp);
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/waffle-chart.ts
function renderWaffleChart(container, data, config = {}) {
  import("d3").then((d3) => renderWaffleChartD3(d3, container, data, config));
}
function renderWaffleChartD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 400, 300);
  const H = config.height ?? 420;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Waffle Chart", data.testResult?.formatted ?? "", W, theme);
  if (data.slices.length === 0) return;
  const total = data.slices.reduce((s, sl) => s + sl.value, 0);
  if (total === 0) return;
  const pcts = data.slices.map((sl) => sl.value / total * 100);
  const floors = pcts.map((p) => Math.floor(p));
  const remainders = pcts.map((p, i) => p - floors[i]);
  let remaining = 100 - floors.reduce((a, b) => a + b, 0);
  const sortedIdx = Array.from({ length: data.slices.length }, (_, i) => i).sort((a, b) => remainders[b] - remainders[a]);
  const cells = [...floors];
  for (let i = 0; i < remaining; i++) {
    cells[sortedIdx[i % sortedIdx.length]] += 1;
  }
  const cellColors = [];
  cells.forEach((count, si) => {
    for (let c = 0; c < count; c++) cellColors.push(si);
  });
  const legendH = Math.ceil(data.slices.length / 3) * 22 + 12;
  const gridAreaH = H - theme.marginTop - legendH - 16;
  const gridAreaW = W - theme.marginLeft - theme.marginRight;
  const cellSize = Math.min(Math.floor(gridAreaH / 10), Math.floor(gridAreaW / 10));
  const gap = Math.max(2, Math.floor(cellSize / 8));
  const gridW = cellSize * 10 + gap * 9;
  const gridH = cellSize * 10 + gap * 9;
  const originX = (W - gridW) / 2;
  const originY = theme.marginTop + 8;
  const g = svg.append("g");
  for (let row = 0; row < 10; row++) {
    for (let col = 0; col < 10; col++) {
      const idx = row * 10 + col;
      const si = cellColors[idx] ?? 0;
      const color = getColor(si, theme);
      const x = originX + col * (cellSize + gap);
      const y = originY + row * (cellSize + gap);
      const sliceData = data.slices[si];
      g.append("rect").attr("x", x).attr("y", y).attr("width", cellSize).attr("height", cellSize).attr("fill", color).attr("opacity", 0.85).attr("rx", Math.max(1, cellSize / 8)).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Slice", sliceData.label),
          formatTooltipRow("Value", sliceData.value.toFixed(2)),
          formatTooltipRow("Share", `${(sliceData.value / total * 100).toFixed(1)}%`)
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
    }
  }
  const lgTop = originY + gridH + 20;
  const lgColW = Math.floor(W / 3);
  data.slices.forEach((sl, si) => {
    const row = Math.floor(si / 3);
    const col = si % 3;
    const lgX = col * lgColW + 12;
    const lgY = lgTop + row * 22;
    const color = getColor(si, theme);
    const pct = (sl.value / total * 100).toFixed(1);
    g.append("rect").attr("x", lgX).attr("y", lgY).attr("width", 12).attr("height", 12).attr("fill", color).attr("opacity", 0.85).attr("rx", 2);
    const labelText = sl.label.length > 14 ? sl.label.slice(0, 13) + "\u2026" : sl.label;
    g.append("text").attr("x", lgX + 16).attr("y", lgY + 10).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(`${labelText} (${pct}%)`);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/sparkline.ts
function renderSparkline(container, data, config = {}) {
  import("d3").then((d3) => renderSparklineD3(d3, container, data, config));
}
function renderSparklineD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 200, 80);
  const H = config.height ?? 80;
  const m = 8;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  const values = [...data.values];
  if (values.length < 2) {
    if (values.length === 1) {
      const color2 = getColor(0, theme);
      svg.append("circle").attr("cx", W / 2).attr("cy", H / 2).attr("r", 3).attr("fill", color2);
    }
    return;
  }
  const color = getColor(0, theme);
  const [yMin, yMax] = d3.extent(values);
  const yPad = (yMax - yMin) * 0.1 || 1;
  const xScale = d3.scaleLinear().domain([0, values.length - 1]).range([m, W - m]);
  const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([H - m, m]);
  const indexed = values.map((v, i) => [i, v]);
  if (config.showArea !== false) {
    const areaFn = d3.area().x((d) => xScale(d[0])).y0(H - m).y1((d) => yScale(d[1])).curve(d3.curveCatmullRom);
    svg.append("path").datum(indexed).attr("d", areaFn).attr("fill", color).attr("opacity", theme.ciOpacity);
  }
  const lineFn = d3.line().x((d) => xScale(d[0])).y((d) => yScale(d[1])).curve(d3.curveCatmullRom);
  svg.append("path").datum(indexed).attr("d", lineFn).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.5);
  const lastVal = values[values.length - 1];
  svg.append("circle").attr("cx", xScale(values.length - 1)).attr("cy", yScale(lastVal)).attr("r", 3).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5);
}

// src/viz/plots/sunburst.ts
function renderSunburst(container, data, config = {}) {
  import("d3").then((d3) => renderSunburstD3(d3, container, data, config));
}
function renderSunburstD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 640, 400);
  const H = config.height ?? 520;
  const maxDepth = config.maxDepth ?? 4;
  const innerRadiusFrac = config.innerRadius ?? 0.25;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Sunburst Chart",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const captionH = config.caption ? 24 : 0;
  const legendW = 150;
  const availW = W - legendW - 8;
  const availH = H - theme.marginTop - theme.marginBottom - captionH;
  const cx = theme.marginLeft + (availW - theme.marginLeft) / 2;
  const cy = theme.marginTop + availH / 2;
  const outerRadius = Math.min(availW - theme.marginLeft, availH) / 2 - 4;
  const innerRadius = outerRadius * innerRadiusFrac;
  if (outerRadius <= 0) return;
  const hierarchyRoot = d3.hierarchy(
    data.root,
    (d) => d.children
  ).sum((d) => d.value ?? 0);
  const total = hierarchyRoot.value ?? 0;
  const partitionLayout = d3.partition().size([2 * Math.PI, outerRadius]);
  const partitioned = partitionLayout(hierarchyRoot);
  const nodes = partitioned.descendants();
  const siblingIndexMap = /* @__PURE__ */ new Map();
  partitioned.eachBefore((node) => {
    const rect = node;
    const parent = rect.parent;
    if (parent !== null) {
      const siblings = parent.children ?? [];
      const idx = siblings.indexOf(rect);
      siblingIndexMap.set(rect, idx >= 0 ? idx : 0);
    }
  });
  const arcGen = d3.arc().startAngle((d) => d.x0).endAngle((d) => d.x1).innerRadius((d) => Math.max(d.y0, innerRadius)).outerRadius((d) => d.y1).padAngle(5e-3).padRadius(outerRadius / 2);
  const g = svg.append("g").attr("transform", `translate(${cx},${cy})`);
  const depth1Nodes = nodes.filter((n) => n.depth === 1);
  nodes.forEach((node) => {
    if (node.depth === 0) return;
    if (node.depth > maxDepth) return;
    if (node.y0 < innerRadius) return;
    const sibIdx = siblingIndexMap.get(node) ?? 0;
    const colorIdx = (node.depth - 1 + sibIdx) % 8;
    const fillColor = getColor(colorIdx, theme);
    const opacity = Math.max(0.4, 1 - node.depth * 0.15);
    const arcAngle = node.x1 - node.x0;
    const arcWidth = node.y1 - node.y0;
    const nodeValue = node.value ?? 0;
    const parentValue = node.parent?.value ?? 0;
    const pctOfParent = parentValue > 0 ? (nodeValue / parentValue * 100).toFixed(1) : "\u2014";
    g.append("path").datum(node).attr("d", arcGen).attr("fill", fillColor).attr("opacity", opacity).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("Name", node.data.name),
        formatTooltipRow("Value", nodeValue.toLocaleString()),
        formatTooltipRow("% of parent", pctOfParent + "%"),
        formatTooltipRow("Depth", String(node.depth))
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    if (arcAngle > 0.15 && arcWidth > 20) {
      const midAngle = (node.x0 + node.x1) / 2;
      const midRadius = (Math.max(node.y0, innerRadius) + node.y1) / 2;
      const midAngleDeg = midAngle * 180 / Math.PI - 90;
      const flip = midAngle > Math.PI;
      const rotateDeg = flip ? midAngleDeg + 180 : midAngleDeg;
      const maxChars = Math.max(2, Math.floor(arcAngle * midRadius / 7));
      const label = node.data.name.length > maxChars ? node.data.name.slice(0, maxChars - 1) + "\u2026" : node.data.name;
      g.append("text").attr("transform", `rotate(${rotateDeg}) translate(${midRadius},0) rotate(${flip ? 180 : 0})`).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall, arcWidth * 0.5)).attr("fill", "#fff").attr("pointer-events", "none").text(label);
    }
  });
  if (innerRadius > 12) {
    g.append("text").attr("text-anchor", "middle").attr("dominant-baseline", "auto").attr("y", -4).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeTitle).attr("font-weight", "700").attr("fill", theme.text).text(total.toLocaleString());
    g.append("text").attr("text-anchor", "middle").attr("dominant-baseline", "hanging").attr("y", 6).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text("total");
  }
  const legendX = cx + outerRadius + 16;
  const legendStartY = theme.marginTop + 8;
  const lineH = 20;
  depth1Nodes.forEach((node, i) => {
    const sibIdx = siblingIndexMap.get(node) ?? 0;
    const colorIdx = sibIdx % 8;
    const fillColor = getColor(colorIdx, theme);
    const ly = legendStartY + i * lineH;
    if (ly + lineH > H - theme.marginBottom - captionH) return;
    svg.append("rect").attr("x", legendX).attr("y", ly).attr("width", 12).attr("height", 12).attr("fill", fillColor).attr("rx", 2);
    const maxLabelChars = Math.max(3, Math.floor((W - legendX - 28) / 7));
    const rawLabel = node.data.name;
    const label = rawLabel.length > maxLabelChars ? rawLabel.slice(0, maxLabelChars - 1) + "\u2026" : rawLabel;
    svg.append("text").attr("x", legendX + 16).attr("y", ly + 10).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/marimekko.ts
function renderMarimekko(container, data, config = {}) {
  import("d3").then((d3) => renderMarimekkoD3(d3, container, data, config));
}
function renderMarimekkoD3(d3, container, data, config) {
  if (data.categories.length === 0 || data.series.length === 0) return;
  const theme = config.theme ?? DEFAULT_THEME;
  const showPct = config.showPercentLabels !== false;
  const showVal = config.showValueLabels === true;
  const legendW = 130;
  const spineH = 6;
  const spineGap = 44;
  const xLabelY = spineGap + spineH + 18;
  const W = config.width ?? Math.max(container.clientWidth || 640, 420);
  const H = config.height ?? 500;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight + legendW,
    bottom: theme.marginBottom + spineH + 20,
    // extra room for spine + x label
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const nCols = data.categories.length;
  const colTotals = Array.from(
    { length: nCols },
    (_, j) => data.series.reduce((sum, s) => sum + (s.values[j] ?? 0), 0)
  );
  const grandTotal = colTotals.reduce((a, b) => a + b, 0);
  if (grandTotal === 0) return;
  const gap = 3;
  const colWidths = colTotals.map((t) => t / grandTotal * width);
  const colX = [];
  colWidths.reduce((acc, w, i) => {
    colX[i] = acc;
    return acc + w;
  }, 0);
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Marimekko Chart",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const yPcts = [0, 25, 50, 75, 100];
  yPcts.forEach((pct) => {
    const yPos = height * (1 - pct / 100);
    g.append("line").attr("x1", 0).attr("x2", width).attr("y1", yPos).attr("y2", yPos).attr("stroke", theme.gridLine).attr("stroke-width", 1).attr("stroke-dasharray", "3,3");
    g.append("text").attr("x", -8).attr("y", yPos + 4).attr("text-anchor", "end").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text(`${pct}%`);
  });
  if (config.yLabel) {
    g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel);
  } else {
    g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.textMuted).text("Percentage (%)");
  }
  data.categories.forEach((cat, ci) => {
    const colTotal = colTotals[ci] ?? 0;
    if (colTotal === 0) return;
    const cx = colX[ci] ?? 0;
    const cw = colWidths[ci] ?? 0;
    let cumPct = 0;
    data.series.forEach((s, si) => {
      const rawVal = s.values[ci] ?? 0;
      const cellPct = rawVal / colTotal;
      const cellH = cellPct * height;
      const rx = cx + gap / 2;
      const ry = height - (cumPct + cellPct) * height;
      const rw = Math.max(cw - gap, 0);
      const rh = Math.max(cellH - gap / 2, 0);
      const color = getColor(si, theme);
      g.append("rect").attr("x", rx).attr("y", ry).attr("width", rw).attr("height", rh).attr("fill", color).attr("opacity", theme.violinOpacity);
      const labelStr = showVal ? String(rawVal) : `${(cellPct * 100).toFixed(1)}%`;
      if (rh > 18 && rw > 30 && (showPct || showVal)) {
        const isDark = isDarkColor(color);
        g.append("text").attr("x", rx + rw / 2).attr("y", ry + rh / 2 + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", Math.min(theme.fontSizeSmall, rh * 0.38, rw * 0.18)).attr("fill", isDark ? "#ffffff" : "#212529").attr("pointer-events", "none").text(labelStr);
      }
      const colPct = (colTotal / grandTotal * 100).toFixed(1);
      g.append("rect").attr("x", rx).attr("y", ry).attr("width", rw).attr("height", rh).attr("fill", "transparent").on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Category", cat),
          formatTooltipRow("Series", s.label),
          formatTooltipRow("Value", rawVal),
          formatTooltipRow("% of column", `${(cellPct * 100).toFixed(1)}%`),
          formatTooltipRow("Column total", `${colTotal} (${colPct}% of total)`)
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
      cumPct += cellPct;
    });
  });
  data.categories.forEach((cat, ci) => {
    const cx = colX[ci] ?? 0;
    const cw = colWidths[ci] ?? 0;
    const midX = cx + cw / 2;
    g.append("text").attr("x", midX).attr("y", height + 18).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(cat);
    g.append("text").attr("x", midX).attr("y", height + 30).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.textMuted).text(`(N=${colTotals[ci] ?? 0})`);
  });
  data.categories.forEach((_, ci) => {
    const cx = colX[ci] ?? 0;
    const cw = colWidths[ci] ?? 0;
    g.append("rect").attr("x", cx + gap / 2).attr("y", height + spineGap).attr("width", Math.max(cw - gap, 0)).attr("height", spineH).attr("fill", getColor(0, theme)).attr("opacity", 0.4).attr("rx", 1);
  });
  if (config.xLabel) {
    g.append("text").attr("x", width / 2).attr("y", height + xLabelY).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel);
  }
  const swatchS = 12;
  const lgX = width + 20;
  const lgY0 = Math.max(0, height / 2 - data.series.length * 20 / 2);
  g.append("text").attr("x", lgX).attr("y", lgY0 - 12).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("font-weight", "600").attr("fill", theme.text).text("Series");
  data.series.forEach((s, si) => {
    const color = getColor(si, theme);
    const ly = lgY0 + si * 20;
    g.append("rect").attr("x", lgX).attr("y", ly).attr("width", swatchS).attr("height", swatchS).attr("fill", color).attr("opacity", theme.violinOpacity).attr("rx", 2);
    g.append("text").attr("x", lgX + swatchS + 6).attr("y", ly + swatchS - 2).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(s.label);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function isDarkColor(hex) {
  const r = parseInt(hex.slice(1, 3), 16) / 255;
  const g = parseInt(hex.slice(3, 5), 16) / 255;
  const b = parseInt(hex.slice(5, 7), 16) / 255;
  const lum = 0.2126 * r + 0.7152 * g + 0.0722 * b;
  return lum < 0.35;
}

// src/viz/plots/chord-diagram.ts
function renderChordDiagram(container, data, config = {}) {
  import("d3").then((d3) => renderChordDiagramD3(d3, container, data, config));
}
function renderChordDiagramD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 600, 400);
  const H = config.height ?? Math.max(container.clientHeight || 500, 400);
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Chord Diagram",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const captionH = config.caption ? 24 : 0;
  const availH = H - theme.marginTop - captionH - 8;
  const cx = W / 2;
  const cy = theme.marginTop + availH / 2;
  const minDim = Math.min(W, availH);
  const outerRadius = config.innerRadius != null ? config.innerRadius + 20 : minDim / 2 * 0.82;
  const innerRadius = config.innerRadius != null ? config.innerRadius : outerRadius - 20;
  const padAngle = config.padAngle ?? 0.02;
  const mutableMatrix = data.matrix.map((row) => [...row]);
  const chordLayout = d3.chord().padAngle(padAngle).sortSubgroups(d3.descending);
  const chords = chordLayout(mutableMatrix);
  const arcGenerator = d3.arc().innerRadius(innerRadius).outerRadius(outerRadius);
  const ribbonGenerator = d3.ribbon().startAngle((sg) => sg.startAngle).endAngle((sg) => sg.endAngle).radius((_sg) => innerRadius);
  const g = svg.append("g").attr("transform", `translate(${cx},${cy})`);
  const sortedChords = [...chords].sort(
    (a, b) => b.source.value + b.target.value - (a.source.value + a.target.value)
  );
  g.selectAll(".ribbon").data(sortedChords).join("path").attr("class", "ribbon").attr("d", (d) => ribbonGenerator(d) ?? "").attr("fill", (d) => getColor(d.source.index, theme)).attr("opacity", 0.65).attr("stroke", theme.background).attr("stroke-width", 0.5).on("mouseover", (event, d) => {
    const srcLabel = data.labels[d.source.index] ?? `Group ${d.source.index}`;
    const tgtLabel = data.labels[d.target.index] ?? `Group ${d.target.index}`;
    const content = [
      formatTooltipRow("From", srcLabel),
      formatTooltipRow("To", tgtLabel),
      formatTooltipRow("Flow", d.source.value)
    ].join("");
    showTooltip(event, content, theme);
  }).on("mouseout", hideTooltip);
  const groupG = g.selectAll(".group").data(chords.groups).join("g").attr("class", "group");
  groupG.append("path").attr("d", (d) => arcGenerator(d) ?? "").attr("fill", (d) => getColor(d.index, theme)).attr("stroke", theme.background).attr("stroke-width", 1).on("mouseover", (event, d) => {
    const label = data.labels[d.index] ?? `Group ${d.index}`;
    const outflow = data.matrix[d.index]?.reduce((s, v) => s + v, 0) ?? 0;
    const inflow = data.matrix.reduce((s, row) => s + (row[d.index] ?? 0), 0);
    const content = [
      formatTooltipRow("Group", label),
      formatTooltipRow("Outflow", outflow),
      formatTooltipRow("Inflow", inflow)
    ].join("");
    showTooltip(event, content, theme);
  }).on("mouseout", hideTooltip);
  const total = data.matrix.reduce(
    (s, row) => s + row.reduce((rs, v) => rs + v, 0),
    0
  );
  const tickStep = total > 0 ? computeTickStep(total) : 1;
  chords.groups.forEach((group) => {
    const groupTotal = data.matrix[group.index]?.reduce((s, v) => s + v, 0) ?? 0;
    const tickCount = Math.floor(groupTotal / tickStep);
    for (let ti = 0; ti <= tickCount; ti++) {
      const tickValue = ti * tickStep;
      const angleRange = group.endAngle - group.startAngle - padAngle;
      const tickAngle = group.startAngle + padAngle / 2 + (groupTotal > 0 ? tickValue / groupTotal * angleRange : 0);
      const sinA = Math.sin(tickAngle - Math.PI / 2);
      const cosA = Math.cos(tickAngle - Math.PI / 2);
      g.append("line").attr("x1", innerRadius * cosA).attr("y1", innerRadius * sinA).attr("x2", (outerRadius + 8) * cosA).attr("y2", (outerRadius + 8) * sinA).attr("stroke", theme.axisLine).attr("stroke-width", 1);
      if (ti > 0) {
        g.append("text").attr("x", (outerRadius + 14) * cosA).attr("y", (outerRadius + 14) * sinA).attr("text-anchor", tickAngle > Math.PI ? "end" : "start").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall - 2).attr("fill", theme.textMuted).text(formatTickValue(tickValue));
      }
    }
  });
  const labelRadius = outerRadius + 32;
  chords.groups.forEach((group) => {
    const angle = (group.startAngle + group.endAngle) / 2;
    const labelAngle = angle - Math.PI / 2;
    const lx = labelRadius * Math.cos(labelAngle);
    const ly = labelRadius * Math.sin(labelAngle);
    const rightSide = Math.cos(labelAngle) > 0;
    g.append("text").attr("x", lx).attr("y", ly).attr("text-anchor", rightSide ? "start" : "end").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("font-weight", "600").attr("fill", theme.text).text(data.labels[group.index] ?? `Group ${group.index}`);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function computeTickStep(total) {
  const rough = total / 8;
  const magnitude = Math.pow(10, Math.floor(Math.log10(rough)));
  const norm = rough / magnitude;
  const step = norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10;
  return step * magnitude;
}
function formatTickValue(v) {
  if (v >= 1e3) return `${(v / 1e3).toFixed(v % 1e3 === 0 ? 0 : 1)}k`;
  return Number.isInteger(v) ? String(v) : v.toFixed(1);
}

// src/viz/plots/arc-diagram.ts
function renderArcDiagram(container, data, config = {}) {
  import("d3").then((d3) => renderArcDiagramD3(d3, container, data, config));
}
function renderArcDiagramD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 700, 400);
  const H = config.height ?? 420;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: theme.marginBottom,
    left: theme.marginLeft
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  const nodeRadius = config.nodeRadius ?? 6;
  const sortByGroup = config.sortByGroup ?? true;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Arc Diagram",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const uniqueGroups = [];
  data.nodes.forEach((n2) => {
    const grp = n2.group ?? "__none__";
    if (!uniqueGroups.includes(grp)) uniqueGroups.push(grp);
  });
  const groupIndex = (grp) => uniqueGroups.indexOf(grp ?? "__none__");
  const orderedNodes = sortByGroup ? [...data.nodes].sort((a, b) => {
    const gi = groupIndex(a.group) - groupIndex(b.group);
    if (gi !== 0) return gi;
    return data.nodes.indexOf(a) - data.nodes.indexOf(b);
  }) : [...data.nodes];
  const nodeY = height * 0.65;
  const n = orderedNodes.length;
  const xStep = n > 1 ? width / (n - 1) : width / 2;
  const xOffset = n > 1 ? 0 : width / 2;
  const nodeX = /* @__PURE__ */ new Map();
  orderedNodes.forEach((node, i) => {
    nodeX.set(node.id, xOffset + i * xStep);
  });
  const degree = /* @__PURE__ */ new Map();
  data.nodes.forEach((node) => degree.set(node.id, 0));
  data.edges.forEach((edge) => {
    degree.set(edge.source, (degree.get(edge.source) ?? 0) + 1);
    degree.set(edge.target, (degree.get(edge.target) ?? 0) + 1);
  });
  g.append("line").attr("x1", 0).attr("x2", width).attr("y1", nodeY).attr("y2", nodeY).attr("stroke", theme.gridLine).attr("stroke-width", 1.5);
  data.edges.forEach((edge) => {
    const x1 = nodeX.get(edge.source);
    const x2 = nodeX.get(edge.target);
    if (x1 == null || x2 == null || edge.source === edge.target) return;
    const srcNode = data.nodes.find((nd) => nd.id === edge.source);
    const srcGroupIdx = groupIndex(srcNode?.group);
    const arcColor = getColor(srcGroupIdx, theme);
    const arcHeight = Math.abs(x2 - x1) / width * (height * 0.55);
    const mx = (x1 + x2) / 2;
    const my = nodeY - arcHeight;
    const rawStroke = 1 + (edge.value ?? 1) * 0.5;
    const strokeWidth = Math.min(rawStroke, 4);
    const pathD = `M ${x1},${nodeY} Q ${mx},${my} ${x2},${nodeY}`;
    g.append("path").attr("d", pathD).attr("fill", "none").attr("stroke", arcColor).attr("stroke-width", strokeWidth).attr("opacity", 0.55).on("mouseover", (event) => {
      const srcLabel = data.nodes.find((nd) => nd.id === edge.source)?.label ?? edge.source;
      const tgtLabel = data.nodes.find((nd) => nd.id === edge.target)?.label ?? edge.target;
      const content = [
        formatTooltipRow("From", srcLabel),
        formatTooltipRow("To", tgtLabel),
        ...edge.value != null ? [formatTooltipRow("Value", edge.value)] : []
      ].join("");
      showTooltip(event, content, theme);
    }).on("mouseout", hideTooltip);
  });
  const rotateLabelThreshold = 8;
  orderedNodes.forEach((node) => {
    const x = nodeX.get(node.id) ?? 0;
    const grpIdx = groupIndex(node.group);
    const color = getColor(grpIdx, theme);
    const displayLabel = node.label ?? node.id;
    const deg = degree.get(node.id) ?? 0;
    g.append("circle").attr("cx", x).attr("cy", nodeY).attr("r", nodeRadius).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
      const groupName = node.group ?? "none";
      const content = [
        formatTooltipRow("Node", displayLabel),
        formatTooltipRow("Group", groupName),
        formatTooltipRow("Degree", deg)
      ].join("");
      showTooltip(event, content, theme);
    }).on("mouseout", hideTooltip);
    if (n <= rotateLabelThreshold) {
      g.append("text").attr("x", x).attr("y", nodeY + nodeRadius + 14).attr("text-anchor", "middle").attr("dominant-baseline", "hanging").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(displayLabel);
    } else {
      g.append("text").attr("transform", `translate(${x},${nodeY + nodeRadius + 8}) rotate(45)`).attr("text-anchor", "start").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(displayLabel);
    }
  });
  if (config.xLabel) {
    g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel);
  }
  if (config.yLabel) {
    g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel);
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/alluvial-plot.ts
function renderAlluvialPlot(container, data, config = {}) {
  import("d3").then((d3) => renderAlluvialD3(d3, container, data, config));
}
function renderAlluvialD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 700, 400);
  const H = config.height ?? 520;
  const nodePadding = config.nodePadding ?? 8;
  const nodeWidth = config.nodeWidth ?? 18;
  const hasStageLabels = (data.stageLabels?.length ?? 0) > 0;
  const margin = {
    top: theme.marginTop + (hasStageLabels ? 24 : 0),
    right: theme.marginRight + 80,
    // extra right room for last stage's labels
    bottom: theme.marginBottom,
    left: theme.marginLeft + 60
    // extra left room for first stage's labels
  };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Alluvial Plot",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  if (data.nodes.length === 0 || data.flows.length === 0) return;
  const uniqueLabels = Array.from(new Set(data.nodes.map((n) => n.label)));
  const labelColor = new Map(uniqueLabels.map((lbl, i) => [lbl, getColor(i, theme)]));
  const stageIndices = Array.from(new Set(data.nodes.map((n) => n.stage))).sort((a, b) => a - b);
  const nStages = stageIndices.length;
  const stageX = (stageIdx) => nStages < 2 ? width / 2 : stageIdx / (nStages - 1) * (width - nodeWidth);
  const nodeById = new Map(data.nodes.map((n) => [n.id, n]));
  const inSum = /* @__PURE__ */ new Map();
  const outSum = /* @__PURE__ */ new Map();
  data.flows.forEach((f) => {
    outSum.set(f.source, (outSum.get(f.source) ?? 0) + f.value);
    inSum.set(f.target, (inSum.get(f.target) ?? 0) + f.value);
  });
  const firstStage = stageIndices[0] ?? 0;
  const lastStage = stageIndices[nStages - 1] ?? 0;
  const nodeTotal = (id) => {
    const node = nodeById.get(id);
    if (!node) return 0;
    if (node.stage === firstStage) return outSum.get(id) ?? 0;
    if (node.stage === lastStage) return inSum.get(id) ?? 0;
    return Math.max(outSum.get(id) ?? 0, inSum.get(id) ?? 0);
  };
  const stageNodes = /* @__PURE__ */ new Map();
  stageIndices.forEach((si) => {
    const nodesInStage = data.nodes.filter((n) => n.stage === si).sort((a, b) => nodeTotal(b.id) - nodeTotal(a.id));
    const totalVal = nodesInStage.reduce((s, n) => s + nodeTotal(n.id), 0);
    const usableHeight = height - nodePadding * Math.max(0, nodesInStage.length - 1);
    const heightScale = totalVal > 0 ? usableHeight / totalVal : 1;
    const lnodes = [];
    let yCursor = 0;
    nodesInStage.forEach((n) => {
      const tv = nodeTotal(n.id);
      const nh = Math.max(4, tv * heightScale);
      lnodes.push({
        id: n.id,
        stage: n.stage,
        label: n.label,
        color: labelColor.get(n.label) ?? getColor(0, theme),
        totalValue: tv,
        inValue: inSum.get(n.id) ?? 0,
        outValue: outSum.get(n.id) ?? 0,
        x: stageX(stageIndices.indexOf(si)),
        y: yCursor,
        height: nh,
        outOffset: 0,
        inOffset: 0
      });
      yCursor += nh + nodePadding;
    });
    stageNodes.set(si, lnodes);
  });
  const layoutNodeMap = /* @__PURE__ */ new Map();
  stageNodes.forEach((lnodes) => lnodes.forEach((ln) => layoutNodeMap.set(ln.id, ln)));
  const layoutFlows = [];
  const sortedFlows = [...data.flows].sort((a, b) => b.value - a.value);
  sortedFlows.forEach((f) => {
    const src = layoutNodeMap.get(f.source);
    const tgt = layoutNodeMap.get(f.target);
    if (!src || !tgt) return;
    const srcHeightScale = src.outValue > 0 ? src.height / src.outValue : 0;
    const tgtHeightScale = tgt.inValue > 0 ? tgt.height / tgt.inValue : 0;
    const ribbonH = Math.max(1, f.value * Math.min(srcHeightScale, tgtHeightScale));
    const sourceY = src.y + src.outOffset;
    const targetY = tgt.y + tgt.inOffset;
    layoutFlows.push({ sourceNode: src, targetNode: tgt, value: f.value, ribbonHeight: ribbonH, sourceY, targetY });
    src.outOffset += ribbonH;
    tgt.inOffset += ribbonH;
  });
  const flowGroup = g.append("g").attr("class", "flows");
  const nodeGroup = g.append("g").attr("class", "nodes");
  layoutFlows.forEach((lf) => {
    const x1 = lf.sourceNode.x + nodeWidth;
    const x2 = lf.targetNode.x;
    const cpX1 = x1 + (x2 - x1) * 0.4;
    const cpX2 = x1 + (x2 - x1) * 0.6;
    const yTop1 = lf.sourceY;
    const yBot1 = lf.sourceY + lf.ribbonHeight;
    const yTop2 = lf.targetY;
    const yBot2 = lf.targetY + lf.ribbonHeight;
    const pathD = [
      `M ${x1.toFixed(2)},${yTop1.toFixed(2)}`,
      `C ${cpX1.toFixed(2)},${yTop1.toFixed(2)} ${cpX2.toFixed(2)},${yTop2.toFixed(2)} ${x2.toFixed(2)},${yTop2.toFixed(2)}`,
      `L ${x2.toFixed(2)},${yBot2.toFixed(2)}`,
      `C ${cpX2.toFixed(2)},${yBot2.toFixed(2)} ${cpX1.toFixed(2)},${yBot1.toFixed(2)} ${x1.toFixed(2)},${yBot1.toFixed(2)}`,
      "Z"
    ].join(" ");
    const srcPct = lf.sourceNode.outValue > 0 ? (lf.value / lf.sourceNode.outValue * 100).toFixed(1) : "\u2013";
    flowGroup.append("path").attr("d", pathD).attr("fill", lf.sourceNode.color).attr("opacity", 0.42).attr("stroke", lf.sourceNode.color).attr("stroke-width", 0.5).attr("stroke-opacity", 0.25).on("mouseover", (event) => {
      showTooltip(event, [
        formatTooltipRow("From", lf.sourceNode.label),
        formatTooltipRow("To", lf.targetNode.label),
        formatTooltipRow("Value", lf.value.toFixed(2)),
        formatTooltipRow("% of source", `${srcPct}%`)
      ].join(""), theme);
    }).on("mouseout", hideTooltip).on("mouseover.highlight", function(event) {
      d3.select(this).attr("opacity", 0.72);
      showTooltip(event, [
        formatTooltipRow("From", lf.sourceNode.label),
        formatTooltipRow("To", lf.targetNode.label),
        formatTooltipRow("Value", lf.value.toFixed(2)),
        formatTooltipRow("% of source", `${srcPct}%`)
      ].join(""), theme);
    }).on("mouseout.highlight", function() {
      d3.select(this).attr("opacity", 0.42);
      hideTooltip();
    });
  });
  stageIndices.forEach((si) => {
    const lnodes = stageNodes.get(si) ?? [];
    const isFirst = si === firstStage;
    const isLast = si === lastStage;
    lnodes.forEach((ln) => {
      nodeGroup.append("rect").attr("x", ln.x).attr("y", ln.y).attr("width", nodeWidth).attr("height", ln.height).attr("fill", ln.color).attr("rx", 3).attr("stroke", theme.background).attr("stroke-width", 1.5).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow("Label", ln.label),
          formatTooltipRow("Stage", ln.stage.toString()),
          formatTooltipRow("Total value", ln.totalValue.toFixed(2))
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
      const labelX = isFirst ? ln.x - 6 : isLast ? ln.x + nodeWidth + 6 : ln.x + nodeWidth / 2;
      const labelAnchor = isFirst ? "end" : isLast ? "start" : "middle";
      const labelY = ln.y + ln.height / 2 + 4;
      nodeGroup.append("text").attr("x", labelX).attr("y", labelY).attr("text-anchor", labelAnchor).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "500").attr("fill", theme.text).attr("pointer-events", "none").text(ln.label);
    });
  });
  if (data.stageLabels) {
    stageIndices.forEach((si, sIdx) => {
      const label = data.stageLabels?.[sIdx] ?? `Stage ${si}`;
      const x = stageX(sIdx) + nodeWidth / 2;
      nodeGroup.append("text").attr("x", x).attr("y", -12).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("font-weight", "700").attr("fill", theme.textMuted).attr("letter-spacing", "0.5").text(label.toUpperCase());
    });
  }
  if (config.xLabel) {
    g.append("text").attr("x", width / 2).attr("y", height + 48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.xLabel);
  }
  if (config.yLabel) {
    g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -margin.left + 16).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(config.yLabel);
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/edge-bundling.ts
function renderEdgeBundling(container, data, config = {}) {
  import("d3").then((d3) => renderEdgeBundlingD3(d3, container, data, config));
}
function renderEdgeBundlingD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? Math.max(container.clientWidth || 640, 400);
  const H = config.height ?? Math.max(W, 480);
  const beta = config.bundlingStrength ?? 0.85;
  const nodeR = config.nodeRadius ?? 4;
  const labelPad = 90;
  const radius = Math.min(W, H) / 2 - labelPad;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`).style("background", theme.background);
  addSubtitle(
    svg,
    config.title ?? "Edge Bundling",
    data.testResult?.formatted ?? "",
    W,
    theme
  );
  const cx = W / 2;
  const cy = H / 2 + theme.marginTop / 2;
  const g = svg.append("g").attr("transform", `translate(${cx},${cy})`);
  if (data.nodes.length === 0) return;
  const colorKey = (n) => n.group ?? n.parent ?? "";
  const uniqueColorKeys = Array.from(new Set(data.nodes.map(colorKey)));
  const colorKeyMap = new Map(uniqueColorKeys.map((k, i) => [k, getColor(i, theme)]));
  const nodeColor = (n) => colorKeyMap.get(colorKey(n)) ?? getColor(0, theme);
  const rawNodes = data.nodes.map((n) => ({
    id: n.id,
    parentId: n.parent === "" ? void 0 : n.parent,
    data: n
  }));
  let root;
  try {
    const stratify = d3.stratify({
      id: (d) => d.id,
      parentId: (d) => d.parentId
    });
    root = stratify(rawNodes);
  } catch {
    g.append("text").attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.textMuted).text("Invalid hierarchy \u2014 check node parent ids.");
    return;
  }
  const cluster = d3.cluster();
  cluster.size([360, radius]);
  cluster(root);
  const nodeIndex = /* @__PURE__ */ new Map();
  const allNodes = [];
  (function collectNodes(node) {
    allNodes.push(node);
    nodeIndex.set(node.data.id, node);
    node.children?.forEach(collectNodes);
  })(root);
  const leafNodes = allNodes.filter((n) => !n.children || n.children.length === 0);
  function lcaPath(a, b) {
    const ancestorsA = a.ancestors();
    const ancestorsB = b.ancestors();
    const ancestorSetA = new Set(ancestorsA.map((n) => n.data.id));
    const lca = ancestorsB.find((n) => ancestorSetA.has(n.data.id));
    if (!lca) return [a, b];
    const lcaIdx = ancestorsA.findIndex((n) => n.data.id === lca.data.id);
    const pathUp = ancestorsA.slice(0, lcaIdx + 1);
    const lcaIdxB = ancestorsB.findIndex((n) => n.data.id === lca.data.id);
    const pathDown = ancestorsB.slice(0, lcaIdxB).reverse();
    return [...pathUp, ...pathDown];
  }
  const bundleCurve = d3.curveBundle.beta(beta);
  const lineRadial = d3.lineRadial().curve(bundleCurve).angle((d) => d.x * Math.PI / 180).radius((d) => d.y);
  const edgeGroup = g.append("g").attr("class", "edges");
  const nodeGroup = g.append("g").attr("class", "nodes");
  const edgePaths = [];
  data.edges.forEach((e) => {
    const srcNode = nodeIndex.get(e.source);
    const tgtNode = nodeIndex.get(e.target);
    if (!srcNode || !tgtNode) return;
    const path = lcaPath(srcNode, tgtNode);
    const pathStr = lineRadial(path);
    if (!pathStr) return;
    const color = nodeColor(srcNode.data);
    const el = edgeGroup.append("path").attr("d", pathStr).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1).attr("opacity", 0.28);
    edgePaths.push({ srcId: e.source, tgtId: e.target, el, srcColor: color });
  });
  let lockedNodeId = null;
  leafNodes.forEach((n) => {
    const angleRad = n.x * Math.PI / 180;
    const nx = Math.sin(angleRad) * n.y;
    const ny = -Math.cos(angleRad) * n.y;
    const color = nodeColor(n.data);
    const label = n.data.label ?? n.data.id;
    nodeGroup.append("circle").attr("cx", nx).attr("cy", ny).attr("r", nodeR).attr("fill", color).attr("stroke", theme.background).attr("stroke-width", 1.5).attr("cursor", "pointer").on("mouseover", (event) => {
      if (lockedNodeId !== null) return;
      highlightNode(n.data.id);
      showTooltip(event, [
        formatTooltipRow("Node", label),
        formatTooltipRow("Group", n.data.group ?? n.data.parent ?? "\u2014"),
        formatTooltipRow(
          "Connections",
          data.edges.filter((e) => e.source === n.data.id || e.target === n.data.id).length.toString()
        )
      ].join(""), theme);
    }).on("mouseout", () => {
      if (lockedNodeId !== null) return;
      resetHighlight();
      hideTooltip();
    }).on("click", (event) => {
      event.stopPropagation();
      if (lockedNodeId === n.data.id) {
        lockedNodeId = null;
        resetHighlight();
        hideTooltip();
      } else {
        lockedNodeId = n.data.id;
        highlightNode(n.data.id);
        showTooltip(event, [
          formatTooltipRow("Node", label),
          formatTooltipRow("Group", n.data.group ?? n.data.parent ?? "\u2014"),
          formatTooltipRow(
            "Connections",
            data.edges.filter((e) => e.source === n.data.id || e.target === n.data.id).length.toString()
          )
        ].join(""), theme);
      }
    });
    const labelAngle = n.x;
    const flipped = labelAngle > 180;
    const rotDeg = flipped ? labelAngle - 90 : labelAngle - 90;
    const labelR = n.y + nodeR + 5;
    const lx = Math.sin(angleRad) * labelR;
    const ly = -Math.cos(angleRad) * labelR;
    const textRotate = flipped ? `rotate(${rotDeg + 180},${lx.toFixed(2)},${ly.toFixed(2)})` : `rotate(${rotDeg},${lx.toFixed(2)},${ly.toFixed(2)})`;
    nodeGroup.append("text").attr("x", lx).attr("y", ly).attr("dy", "0.35em").attr("text-anchor", flipped ? "end" : "start").attr("transform", textRotate).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).attr("pointer-events", "none").text(label);
  });
  svg.on("click", () => {
    if (lockedNodeId !== null) {
      lockedNodeId = null;
      resetHighlight();
      hideTooltip();
    }
  });
  function highlightNode(nodeId) {
    edgePaths.forEach((ep) => {
      const connected = ep.srcId === nodeId || ep.tgtId === nodeId;
      ep.el.attr("opacity", connected ? 0.9 : 0.05).attr("stroke-width", connected ? 2 : 0.8);
    });
  }
  function resetHighlight() {
    edgePaths.forEach((ep) => {
      ep.el.attr("opacity", 0.28).attr("stroke-width", 1);
    });
  }
  if (uniqueColorKeys.filter((k) => k !== "").length > 1) {
    const legendKeys = uniqueColorKeys.filter((k) => k !== "");
    const lgX = -W / 2 + theme.marginLeft;
    const lgY = H / 2 - theme.marginBottom - legendKeys.length * 18;
    legendKeys.forEach((k, i) => {
      g.append("rect").attr("x", lgX).attr("y", lgY + i * 18).attr("width", 10).attr("height", 10).attr("fill", colorKeyMap.get(k) ?? getColor(i, theme)).attr("rx", 2);
      g.append("text").attr("x", lgX + 14).attr("y", lgY + i * 18 + 9).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(k);
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/plots/dendrogram.ts
function renderDendrogram(container, data, config = {}) {
  import("d3").then((d3) => renderDendrogramD3(d3, container, data, config));
}
function renderDendrogramD3(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 600;
  const H = config.height ?? 480;
  const showCut = config.showCutLine !== false;
  const showLabels = config.showLabels !== false;
  const colorSubs = config.colorSubtrees !== false;
  const margin = {
    top: theme.marginTop,
    right: theme.marginRight,
    bottom: showLabels ? Math.max(theme.marginBottom, 100) : theme.marginBottom,
    left: theme.marginLeft
  };
  const plotW = W - margin.left - margin.right;
  const plotH = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const { merges, heights, dendrogramOrder, labels, k } = data;
  const n = merges.length + 1;
  const nodeX = new Float64Array(2 * n - 1);
  const nodeY = new Float64Array(2 * n - 1);
  const leafPos = new Float64Array(n);
  for (let i = 0; i < dendrogramOrder.length; i++) {
    leafPos[dendrogramOrder[i]] = i;
  }
  for (let i = 0; i < n; i++) {
    nodeX[i] = leafPos[i];
    nodeY[i] = 0;
  }
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s];
    const internalIdx = n + s;
    nodeX[internalIdx] = (nodeX[m.a] + nodeX[m.b]) / 2;
    nodeY[internalIdx] = heights[s];
  }
  const nodeCluster = new Int32Array(2 * n - 1);
  for (let i = 0; i < n; i++) {
    nodeCluster[i] = labels[i];
  }
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s];
    const cA = nodeCluster[m.a];
    const cB = nodeCluster[m.b];
    nodeCluster[n + s] = cA === cB ? cA : -1;
  }
  const subtreeSize = new Float64Array(2 * n - 1);
  for (let i = 0; i < n; i++) subtreeSize[i] = 1;
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s];
    subtreeSize[n + s] = subtreeSize[m.a] + subtreeSize[m.b];
  }
  const maxHeight = Math.max(...heights) * 1.05;
  const xScale = d3.scaleLinear().domain([0, n - 1]).range([0, plotW]);
  const yScale = d3.scaleLinear().domain([0, maxHeight]).range([plotH, 0]);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  const subtitle = `${data.linkage} linkage, K = ${k}, cophenetic r = ${data.copheneticCorrelation.toFixed(3)}`;
  addSubtitle(svg, config.title ?? "Dendrogram", subtitle, W, theme);
  if (config.caption) {
    addCaption(svg, config.caption, W, H, theme);
  }
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const yAxis = d3.axisLeft(yScale).ticks(6);
  g.append("g").attr("class", "y-axis").call(yAxis).call((sel) => {
    sel.select(".domain").attr("stroke", theme.axisLine);
    sel.selectAll(".tick line").attr("stroke", theme.gridLine);
    sel.selectAll(".tick text").attr("fill", theme.textMuted).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall);
  });
  g.append("g").attr("class", "grid").call(d3.axisLeft(yScale).ticks(6).tickSize(-plotW).tickFormat(() => "")).call((sel) => {
    sel.select(".domain").remove();
    sel.selectAll(".tick line").attr("stroke", theme.gridLine).attr("stroke-opacity", 0.5);
  });
  g.append("text").attr("transform", "rotate(-90)").attr("x", -plotH / 2).attr("y", -margin.left + 16).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).text("Merge Height");
  for (let s = 0; s < merges.length; s++) {
    const m = merges[s];
    const internalIdx = n + s;
    const ax = xScale(nodeX[m.a]);
    const ay = yScale(nodeY[m.a]);
    const bx = xScale(nodeX[m.b]);
    const by = yScale(nodeY[m.b]);
    const my = yScale(nodeY[internalIdx]);
    const pathD = `M${ax},${ay} V${my} H${bx} V${by}`;
    const cluster = nodeCluster[internalIdx];
    const color = colorSubs && cluster >= 0 ? getColor(cluster, theme) : theme.axisLine;
    const sizeA = subtreeSize[m.a];
    const sizeB = subtreeSize[m.b];
    const sizeC = subtreeSize[internalIdx];
    const mergeH = heights[s];
    g.append("path").attr("d", pathD).attr("fill", "none").attr("stroke", color).attr("stroke-width", 1.5).attr("stroke-linejoin", "round").style("cursor", "pointer").on("mouseenter", (event) => {
      const content = [
        formatTooltipRow("Merge height", mergeH.toFixed(3)),
        formatTooltipRow("Left subtree", String(sizeA)),
        formatTooltipRow("Right subtree", String(sizeB)),
        formatTooltipRow("Combined", String(sizeC))
      ].join("");
      showTooltip(event, content, theme);
    }).on("mousemove", (event) => {
      const content = [
        formatTooltipRow("Merge height", mergeH.toFixed(3)),
        formatTooltipRow("Left subtree", String(sizeA)),
        formatTooltipRow("Right subtree", String(sizeB)),
        formatTooltipRow("Combined", String(sizeC))
      ].join("");
      showTooltip(event, content, theme);
    }).on("mouseleave", () => hideTooltip());
  }
  if (showCut && k > 1 && k <= n) {
    const belowIdx = n - k - 1;
    const aboveIdx = n - k;
    const hBelow = belowIdx >= 0 ? heights[belowIdx] : 0;
    const hAbove = aboveIdx < heights.length ? heights[aboveIdx] : maxHeight;
    const cutY = yScale((hBelow + hAbove) / 2);
    g.append("line").attr("x1", -8).attr("x2", plotW + 8).attr("y1", cutY).attr("y2", cutY).attr("stroke", theme.textMuted).attr("stroke-width", 1.5).attr("stroke-dasharray", "6,4").attr("opacity", 0.7);
    g.append("text").attr("x", plotW + 4).attr("y", cutY - 6).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textMuted).attr("text-anchor", "end").text(`K = ${k}`);
  }
  if (showLabels) {
    const obsLabels = data.observationLabels;
    const step = n > 80 ? Math.ceil(n / 40) : 1;
    for (let i = 0; i < dendrogramOrder.length; i++) {
      if (i % step !== 0) continue;
      const leafIdx = dendrogramOrder[i];
      const lx = xScale(i);
      const ly = plotH + 8;
      const label = obsLabels?.[leafIdx] ?? String(leafIdx + 1);
      const cluster = labels[leafIdx];
      g.append("text").attr("x", lx).attr("y", ly).attr("transform", `rotate(-60, ${lx}, ${ly})`).attr("font-family", theme.fontFamily).attr("font-size", Math.min(theme.fontSizeSmall, Math.max(8, plotW / n * 0.8))).attr("fill", colorSubs ? getColor(cluster, theme) : theme.textMuted).attr("text-anchor", "end").attr("dominant-baseline", "auto").text(label);
    }
  }
}

// src/viz/plots/fa-plot.ts
function renderFAPlot(container, data, config = {}) {
  import("d3").then((d3) => {
    const type = config.type ?? "loadings";
    if (type === "scree") renderScree2(d3, container, data, config);
    else if (type === "loadings") renderLoadings(d3, container, data, config);
    else if (type === "path") renderPath(d3, container, data, config);
    else if (type === "communality") renderCommunality(d3, container, data, config);
    else if (type === "factor-correlation") renderFactorCorrelation(d3, container, data, config);
    else if (type === "fit-indices") renderFitIndices(d3, container, data, config);
  });
}
function isCFA(data) {
  return "parameterEstimates" in data;
}
function renderScree2(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 520, H = config.height ?? 380;
  const margin = { top: theme.marginTop, right: theme.marginRight, bottom: theme.marginBottom, left: theme.marginLeft };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  const diag = config.diagnostics;
  const subtitleParts = [];
  if (diag) {
    subtitleParts.push(`Parallel Analysis: ${diag.parallelSuggested} factors`);
    subtitleParts.push(`MAP: ${diag.mapSuggested} factors`);
  }
  addSubtitle(svg, config.title ?? "Scree Plot with Parallel Analysis", subtitleParts.join("; "), W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const eigenvalues = data.eigenvalues;
  const comps = eigenvalues.map((_, i) => i + 1);
  const maxEv = Math.max(...eigenvalues, ...diag?.parallelSimulated ?? [], 1);
  const xScale = d3.scaleBand().domain(comps).range([0, width]).padding(0.25);
  const yScale = d3.scaleLinear().domain([0, maxEv * 1.15]).range([height, 0]).nice();
  g.selectAll(".grid").data(yScale.ticks(5)).join("line").attr("x1", 0).attr("x2", width).attr("y1", (d) => yScale(d)).attr("y2", (d) => yScale(d)).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", 0).attr("x2", width).attr("y1", yScale(1)).attr("y2", yScale(1)).attr("stroke", theme.axisLine).attr("stroke-dasharray", "6,3").attr("stroke-width", 1.5);
  g.append("text").attr("x", width + 4).attr("y", yScale(1) + 4).attr("font-size", 9).attr("fill", theme.textMuted).text("Kaiser");
  eigenvalues.forEach((ev, i) => {
    const x = xScale(i + 1) ?? 0;
    const barColor = i < data.nFactors ? getColor(0, theme) : theme.gridLine;
    g.append("rect").attr("x", x).attr("y", yScale(ev)).attr("width", xScale.bandwidth()).attr("height", height - yScale(ev)).attr("fill", barColor).attr("opacity", 0.8).attr("rx", 2).on("mouseover", (event) => {
      const pctVar = eigenvalues.length > 0 ? (ev / eigenvalues.reduce((s, v) => s + v, 0) * 100).toFixed(1) : "0";
      showTooltip(event, [
        formatTooltipRow("Component", `${i + 1}`),
        formatTooltipRow("Eigenvalue", ev),
        formatTooltipRow("% Variance", `${pctVar}%`)
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    g.append("text").attr("x", x + xScale.bandwidth() / 2).attr("y", yScale(ev) - 4).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(ev.toFixed(2));
  });
  if (config.showParallelAnalysis !== false && diag?.parallelSimulated) {
    const sim = diag.parallelSimulated;
    const lineData = comps.map((c, i) => ({
      x: (xScale(c) ?? 0) + xScale.bandwidth() / 2,
      y: yScale(sim[i] ?? 0)
    }));
    const line = d3.line().x((d) => d.x).y((d) => d.y);
    g.append("path").datum(lineData).attr("d", line).attr("fill", "none").attr("stroke", getColor(2, theme)).attr("stroke-width", 2).attr("stroke-dasharray", "5,3");
    lineData.forEach((pt) => {
      g.append("circle").attr("cx", pt.x).attr("cy", pt.y).attr("r", 3).attr("fill", getColor(2, theme));
    });
    const legY = 8;
    g.append("line").attr("x1", width - 120).attr("x2", width - 90).attr("y1", legY).attr("y2", legY).attr("stroke", getColor(2, theme)).attr("stroke-width", 2).attr("stroke-dasharray", "5,3");
    g.append("text").attr("x", width - 86).attr("y", legY + 4).attr("font-size", 9).attr("fill", theme.textMuted).text("Simulated (95th)");
  }
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Factor");
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -48).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Eigenvalue");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderLoadings(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const threshold = config.loadingThreshold ?? 0.3;
  const nItems = data.loadings.length;
  const nFactors = data.nFactors;
  const itemOrder = Array.from({ length: nItems }, (_, i) => i);
  itemOrder.sort((a, b) => {
    const rowA = data.loadings[a];
    const rowB = data.loadings[b];
    const maxFactorA = rowA.indexOf(Math.max(...rowA.map(Math.abs)));
    const maxFactorB = rowB.indexOf(Math.max(...rowB.map(Math.abs)));
    if (maxFactorA !== maxFactorB) return maxFactorA - maxFactorB;
    return Math.max(...rowB.map(Math.abs)) - Math.max(...rowA.map(Math.abs));
  });
  const cellW = 56, cellH = 28;
  const showComm = true;
  const extraCols = showComm ? 1 : 0;
  const W = config.width ?? (nFactors + extraCols) * cellW + 140;
  const H = config.height ?? nItems * cellH + 100;
  const margin = { top: 60, right: 20, bottom: 30, left: 110 };
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Factor Loadings", data.formatted, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu);
  const commScale = d3.scaleSequential().domain([0, 1]).interpolator(d3.interpolateGreens);
  const varLabels = config.variableLabels ?? data.variableNames;
  const facLabels = config.factorLabels ?? data.factorNames;
  facLabels.forEach((label, ci) => {
    g.append("text").attr("x", ci * cellW + cellW / 2).attr("y", -8).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.text).text(label);
  });
  if (showComm) {
    g.append("text").attr("x", nFactors * cellW + cellW / 2).attr("y", -8).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.text).text("h\xB2");
  }
  itemOrder.forEach((origIdx, rowIdx) => {
    const row = data.loadings[origIdx];
    const label = varLabels[origIdx] ?? `V${origIdx + 1}`;
    g.append("text").attr("x", -6).attr("y", rowIdx * cellH + cellH / 2 + 4).attr("text-anchor", "end").attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(label);
    row.forEach((val, ci) => {
      const dimmed = Math.abs(val) < threshold;
      g.append("rect").attr("x", ci * cellW).attr("y", rowIdx * cellH).attr("width", cellW - 2).attr("height", cellH - 2).attr("rx", 3).attr("fill", colorScale(val)).attr("opacity", dimmed ? 0.3 : 1).on("mouseover", (event) => {
        showTooltip(event, [
          formatTooltipRow(label, ""),
          formatTooltipRow(facLabels[ci] ?? `F${ci + 1}`, val)
        ].join(""), theme);
      }).on("mouseout", hideTooltip);
      g.append("text").attr("x", ci * cellW + cellW / 2).attr("y", rowIdx * cellH + cellH / 2 + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", Math.min(10, cellW / 5)).attr("fill", Math.abs(val) > 0.6 ? "#fff" : theme.text).attr("opacity", dimmed ? 0.5 : 1).text(val.toFixed(2));
    });
    if (showComm) {
      const comm = data.communalities[origIdx] ?? 0;
      g.append("rect").attr("x", nFactors * cellW).attr("y", rowIdx * cellH).attr("width", cellW - 2).attr("height", cellH - 2).attr("rx", 3).attr("fill", commScale(comm));
      g.append("text").attr("x", nFactors * cellW + cellW / 2).attr("y", rowIdx * cellH + cellH / 2 + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", Math.min(10, cellW / 5)).attr("fill", comm > 0.6 ? "#fff" : theme.text).text(comm.toFixed(2));
    }
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function resolvePathStyle(s) {
  return {
    factorRx: s?.factorRx ?? 56,
    factorRy: s?.factorRy ?? 28,
    factorStroke: s?.factorStroke ?? 2.5,
    factorFontSize: s?.factorFontSize ?? 14,
    factorFontWeight: s?.factorFontWeight ?? "800",
    factorGlow: s?.factorGlow ?? true,
    factorGlowRadius: s?.factorGlowRadius ?? 6,
    factorHighlight: s?.factorHighlight ?? true,
    itemWidth: s?.itemWidth ?? 96,
    itemHeight: s?.itemHeight ?? 30,
    itemRadius: s?.itemRadius ?? 5,
    itemStroke: s?.itemStroke ?? 1.2,
    itemFontSize: s?.itemFontSize ?? 11,
    itemFontWeight: s?.itemFontWeight ?? "600",
    itemAccentWidth: s?.itemAccentWidth ?? 4,
    itemGradient: s?.itemGradient ?? true,
    errorRx: s?.errorRx ?? 14,
    errorRy: s?.errorRy ?? 12,
    errorStroke: s?.errorStroke ?? 1.2,
    errorFontSize: s?.errorFontSize ?? 9,
    errorSelfLoop: s?.errorSelfLoop ?? true,
    errorSelfLoopSize: s?.errorSelfLoopSize ?? 12,
    arrowMinWidth: s?.arrowMinWidth ?? 0.8,
    arrowMaxWidth: s?.arrowMaxWidth ?? 2.5,
    arrowMinOpacity: s?.arrowMinOpacity ?? 0.4,
    arrowMaxOpacity: s?.arrowMaxOpacity ?? 0.95,
    arrowMarkerSize: s?.arrowMarkerSize ?? 6,
    arrowCurvature: s?.arrowCurvature ?? 0.4,
    loadingFontSize: s?.loadingFontSize ?? 9.5,
    loadingFontWeight: s?.loadingFontWeight ?? "700",
    loadingPillRadius: s?.loadingPillRadius ?? 7,
    loadingLabelOffset: s?.loadingLabelOffset ?? -6,
    loadingPosition: s?.loadingPosition ?? 0.6,
    halo: s?.halo ?? true,
    haloColor: s?.haloColor ?? "#ffffff",
    haloWidth: s?.haloWidth ?? 3.5,
    covColor: s?.covColor ?? "#8892a0",
    covMinWidth: s?.covMinWidth ?? 1,
    covMaxWidth: s?.covMaxWidth ?? 3,
    covMarkerSize: s?.covMarkerSize ?? 6,
    covLabelFontSize: s?.covLabelFontSize ?? 10,
    covNestSpacing: s?.covNestSpacing ?? 20,
    errorArrowWidth: s?.errorArrowWidth ?? 0.8,
    errorArrowMarkerSize: s?.errorArrowMarkerSize ?? 5,
    itemGap: s?.itemGap ?? 6,
    factorGroupGap: s?.factorGroupGap ?? 32,
    arrowSpan: s?.arrowSpan ?? 120,
    errorSpan: s?.errorSpan ?? 36,
    covArcReserve: s?.covArcReserve ?? 70,
    topPadding: s?.topPadding ?? 74,
    bottomPadding: s?.bottomPadding ?? 44,
    rightPadding: s?.rightPadding ?? 50,
    factorColors: s?.factorColors ?? [],
    fitCardWidth: s?.fitCardWidth ?? 210,
    fitCardHeight: s?.fitCardHeight ?? 118,
    crossLoadingThreshold: s?.crossLoadingThreshold ?? 0.2,
    crossLoadingDash: s?.crossLoadingDash ?? "5,4",
    crossLoadingOpacity: s?.crossLoadingOpacity ?? 0.3,
    showGroupBrackets: s?.showGroupBrackets ?? true,
    groupBracketOpacity: s?.groupBracketOpacity ?? 0.12,
    showShadows: s?.showShadows ?? true
  };
}
function fColor(fi, ps, theme) {
  if (ps.factorColors.length > 0) return ps.factorColors[fi % ps.factorColors.length];
  return getColor(fi, theme);
}
function renderPath(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const ps = resolvePathStyle(config.pathStyle);
  const showErrors = config.showErrorTerms !== false;
  const showFit = config.showFitBox !== false;
  const nFactors = data.nFactors;
  const nItems = data.loadings.length;
  const varLabels = config.variableLabels ?? data.variableNames;
  const facLabels = config.factorLabels ?? data.factorNames;
  const hasCFA = isCFA(data);
  const factorItems = Array.from({ length: nFactors }, () => []);
  if (hasCFA) {
    const cfaModel = data.model;
    for (let fi = 0; fi < nFactors; fi++) {
      const key = facLabels[fi];
      const items = cfaModel[key];
      if (items) {
        for (let j = 0; j < items.length; j++) factorItems[fi].push(items[j]);
      }
    }
  } else {
    for (let i = 0; i < nItems; i++) {
      const row = data.loadings[i];
      let maxF = 0, maxVal = 0;
      for (let f = 0; f < nFactors; f++) {
        if (Math.abs(row[f]) > maxVal) {
          maxVal = Math.abs(row[f]);
          maxF = f;
        }
      }
      factorItems[maxF].push(i);
    }
  }
  const groupHeights = factorItems.map(
    (items) => Math.max(items.length * (ps.itemHeight + ps.itemGap) - ps.itemGap, ps.factorRy * 2)
  );
  const totalGroupH = groupHeights.reduce((s, h) => s + h, 0) + (nFactors - 1) * ps.factorGroupGap;
  const bottomPad = ps.bottomPadding;
  const factorCx = ps.covArcReserve + ps.factorRx + 16;
  const itemLeft = factorCx + ps.factorRx + ps.arrowSpan;
  const itemRight = itemLeft + ps.itemWidth;
  const errorCx = showErrors ? itemRight + ps.errorSpan + ps.errorRx : itemRight + 10;
  const selfLoopExtra = showErrors && ps.errorSelfLoop ? ps.errorSelfLoopSize + 30 : 0;
  const W = config.width ?? errorCx + (showErrors ? ps.errorRx : 0) + selfLoopExtra + ps.rightPadding;
  const H = config.height ?? ps.topPadding + totalGroupH + bottomPad;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background).attr("xmlns", "http://www.w3.org/2000/svg");
  addSubtitle(
    svg,
    config.title ?? (hasCFA ? "Confirmatory Factor Analysis \u2014 Path Diagram" : "Exploratory Factor Analysis \u2014 Path Diagram"),
    data.formatted,
    W,
    theme
  );
  const defs = svg.append("defs");
  const g = svg.append("g");
  if (ps.factorGlow) {
    const glowFilter = defs.append("filter").attr("id", "fa-glow").attr("x", "-40%").attr("y", "-40%").attr("width", "180%").attr("height", "180%");
    glowFilter.append("feGaussianBlur").attr("in", "SourceGraphic").attr("stdDeviation", ps.factorGlowRadius).attr("result", "blur");
    glowFilter.append("feColorMatrix").attr("in", "blur").attr("type", "saturate").attr("values", "0.35").attr("result", "desatBlur");
    const glowMerge = glowFilter.append("feMerge");
    glowMerge.append("feMergeNode").attr("in", "desatBlur");
    glowMerge.append("feMergeNode").attr("in", "SourceGraphic");
  }
  if (ps.showShadows) {
    const shadowMed = defs.append("filter").attr("id", "fa-shadow-med").attr("x", "-20%").attr("y", "-15%").attr("width", "140%").attr("height", "140%");
    shadowMed.append("feDropShadow").attr("dx", 0).attr("dy", 2).attr("stdDeviation", 4).attr("flood-color", "rgba(0,0,0,0.10)");
    const shadowSm = defs.append("filter").attr("id", "fa-shadow-sm").attr("x", "-10%").attr("y", "-10%").attr("width", "120%").attr("height", "130%");
    shadowSm.append("feDropShadow").attr("dx", 0).attr("dy", 1).attr("stdDeviation", 2).attr("flood-color", "rgba(0,0,0,0.07)");
  }
  for (let f = 0; f < nFactors; f++) {
    const color = fColor(f, ps, theme);
    const grad = defs.append("radialGradient").attr("id", `fa-fgrad-${f}`).attr("cx", "35%").attr("cy", "28%").attr("r", "70%");
    grad.append("stop").attr("offset", "0%").attr("stop-color", "#fff").attr("stop-opacity", 0.55);
    grad.append("stop").attr("offset", "35%").attr("stop-color", color).attr("stop-opacity", 0.18);
    grad.append("stop").attr("offset", "100%").attr("stop-color", color).attr("stop-opacity", 0.06);
    if (ps.itemGradient) {
      const lg = defs.append("linearGradient").attr("id", `fa-igrad-${f}`).attr("x1", "0%").attr("y1", "0%").attr("x2", "100%").attr("y2", "0%");
      lg.append("stop").attr("offset", "0%").attr("stop-color", color).attr("stop-opacity", 0.08);
      lg.append("stop").attr("offset", "100%").attr("stop-color", color).attr("stop-opacity", 0);
    }
  }
  const mw = ps.arrowMarkerSize, mh = ps.arrowMarkerSize * 0.7;
  for (let f = 0; f < nFactors; f++) {
    const color = fColor(f, ps, theme);
    defs.append("marker").attr("id", `fa-arr-${f}`).attr("viewBox", "0 0 10 7").attr("refX", 9.5).attr("refY", 3.5).attr("markerWidth", mw).attr("markerHeight", mh).attr("orient", "auto").append("path").attr("d", "M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z").attr("fill", color);
  }
  defs.append("marker").attr("id", "fa-arr-ns").attr("viewBox", "0 0 10 7").attr("refX", 9.5).attr("refY", 3.5).attr("markerWidth", mw).attr("markerHeight", mh).attr("orient", "auto").append("path").attr("d", "M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z").attr("fill", theme.textMuted);
  const cmw = ps.covMarkerSize, cmh = ps.covMarkerSize * 0.7;
  defs.append("marker").attr("id", "fa-cov-s").attr("viewBox", "0 0 10 7").attr("refX", 1).attr("refY", 3.5).attr("markerWidth", cmw).attr("markerHeight", cmh).attr("orient", "auto").append("path").attr("d", "M10,0.5 L1,3.5 L10,6.5 L8,3.5 Z").attr("fill", ps.covColor);
  defs.append("marker").attr("id", "fa-cov-e").attr("viewBox", "0 0 10 7").attr("refX", 9).attr("refY", 3.5).attr("markerWidth", cmw).attr("markerHeight", cmh).attr("orient", "auto").append("path").attr("d", "M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z").attr("fill", ps.covColor);
  const emw = ps.errorArrowMarkerSize, emh = ps.errorArrowMarkerSize * 0.7;
  defs.append("marker").attr("id", "fa-arr-err").attr("viewBox", "0 0 10 7").attr("refX", 9).attr("refY", 3.5).attr("markerWidth", emw).attr("markerHeight", emh).attr("orient", "auto").append("path").attr("d", "M0,0.5 L9,3.5 L0,6.5 L2,3.5 Z").attr("fill", theme.axisLine);
  const factorCy = [];
  const itemCy = new Array(nItems).fill(0);
  let yAccum = ps.topPadding;
  for (let fi = 0; fi < nFactors; fi++) {
    const items = factorItems[fi];
    const gh = groupHeights[fi];
    factorCy.push(yAccum + gh / 2);
    for (let j = 0; j < items.length; j++) {
      itemCy[items[j]] = yAccum + j * (ps.itemHeight + ps.itemGap) + ps.itemHeight / 2;
    }
    yAccum += gh + ps.factorGroupGap;
  }
  if (ps.showGroupBrackets) {
    for (let fi = 0; fi < nFactors; fi++) {
      const items = factorItems[fi];
      if (items.length < 2) continue;
      const firstY = itemCy[items[0]];
      const lastY = itemCy[items[items.length - 1]];
      g.append("line").attr("x1", itemLeft - 6).attr("y1", firstY).attr("x2", itemLeft - 6).attr("y2", lastY).attr("stroke", fColor(fi, ps, theme)).attr("stroke-width", 1.5).attr("opacity", ps.groupBracketOpacity).attr("stroke-linecap", "round");
    }
  }
  const Phi = data.factorCorrelations;
  for (let fi = 0; fi < nFactors; fi++) {
    for (let fj = fi + 1; fj < nFactors; fj++) {
      const corr = Phi[fi][fj];
      if (Math.abs(corr) < 0.01) continue;
      const y1 = factorCy[fi], y2 = factorCy[fj];
      const span = Math.abs(y2 - y1);
      const nestLevel = fj - fi;
      const arcX = factorCx - ps.factorRx - 14 - nestLevel * ps.covNestSpacing - span * 0.04;
      const absCorr = Math.abs(corr);
      const strokeW = ps.covMinWidth + absCorr * (ps.covMaxWidth - ps.covMinWidth);
      const path = `M${factorCx - ps.factorRx - 2},${y1} C${arcX},${y1} ${arcX},${y2} ${factorCx - ps.factorRx - 2},${y2}`;
      g.append("path").attr("d", path).attr("fill", "none").attr("stroke", ps.covColor).attr("stroke-width", strokeW).attr("stroke-linecap", "round").attr("marker-start", "url(#fa-cov-s)").attr("marker-end", "url(#fa-cov-e)");
      const pillW = 38, pillH = 17;
      const t = 0.5, t1 = 0.5;
      const labelX = t1 * t1 * t1 * (factorCx - ps.factorRx - 2) + 3 * t1 * t1 * t * arcX + 3 * t1 * t * t * arcX + t * t * t * (factorCx - ps.factorRx - 2);
      const labelY = t1 * t1 * t1 * y1 + 3 * t1 * t1 * t * y1 + 3 * t1 * t * t * y2 + t * t * t * y2;
      g.append("rect").attr("x", labelX - pillW / 2).attr("y", labelY - pillH / 2).attr("width", pillW).attr("height", pillH).attr("rx", 8).attr("fill", theme.background).attr("stroke", theme.gridLine).attr("stroke-width", 0.7);
      g.append("text").attr("x", labelX).attr("y", labelY + 3.5).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", ps.covLabelFontSize).attr("font-weight", "600").attr("fill", ps.covColor).text(corr.toFixed(2));
    }
  }
  for (let fi = 0; fi < nFactors; fi++) {
    const items = factorItems[fi];
    const color = fColor(fi, ps, theme);
    for (let idx = 0; idx < items.length; idx++) {
      const itemIdx = items[idx];
      if (itemIdx >= nItems) continue;
      const loading = data.loadings[itemIdx][fi];
      if (Math.abs(loading) < 0.01) continue;
      let sig = true;
      let stars = "";
      if (hasCFA) {
        const pe = data.parameterEstimates.loadings[fi]?.[idx];
        if (pe) {
          sig = pe.pValue < 0.05;
          stars = pe.pValue < 1e-3 ? "***" : pe.pValue < 0.01 ? "**" : pe.pValue < 0.05 ? "*" : "";
        }
      }
      const x1 = factorCx + ps.factorRx + 2;
      const y1 = factorCy[fi];
      const x2 = itemLeft - 2;
      const y2 = itemCy[itemIdx];
      const dy = y2 - y1;
      const curv = ps.arrowCurvature;
      const cp1x = x1 + ps.arrowSpan * curv, cp1y = y1 + dy * 0.05;
      const cp2x = x2 - ps.arrowSpan * (curv * 0.85), cp2y = y2;
      const absL = Math.abs(loading);
      const strokeW = ps.arrowMinWidth + absL * (ps.arrowMaxWidth - ps.arrowMinWidth);
      const opacity = ps.arrowMinOpacity + absL * (ps.arrowMaxOpacity - ps.arrowMinOpacity);
      const arrowColor = sig ? color : theme.textMuted;
      const markerId = sig ? `fa-arr-${fi}` : "fa-arr-ns";
      g.append("path").attr("d", `M${x1},${y1} C${cp1x},${cp1y} ${cp2x},${cp2y} ${x2},${y2}`).attr("fill", "none").attr("stroke", arrowColor).attr("stroke-width", strokeW).attr("stroke-linecap", "round").attr("opacity", opacity).attr("marker-end", `url(#${markerId})`);
      const mt = ps.loadingPosition, mt1 = 1 - mt;
      const lx = mt1 * mt1 * mt1 * x1 + 3 * mt1 * mt1 * mt * cp1x + 3 * mt1 * mt * mt * cp2x + mt * mt * mt * x2;
      const ly = mt1 * mt1 * mt1 * y1 + 3 * mt1 * mt1 * mt * cp1y + 3 * mt1 * mt * mt * cp2y + mt * mt * mt * y2;
      const labelText = `${loading.toFixed(2)}${stars}`;
      const lbl = g.append("text").attr("x", lx).attr("y", ly + ps.loadingLabelOffset).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", ps.loadingFontSize).attr("font-weight", ps.loadingFontWeight).attr("fill", arrowColor).text(labelText);
      if (ps.halo) lbl.attr("paint-order", "stroke").attr("stroke", ps.haloColor).attr("stroke-width", ps.haloWidth).attr("stroke-linejoin", "round");
    }
  }
  if (!hasCFA) {
    for (let i = 0; i < nItems; i++) {
      const row = data.loadings[i];
      let primaryF = 0, maxVal = 0;
      for (let f = 0; f < nFactors; f++) {
        if (Math.abs(row[f]) > maxVal) {
          maxVal = Math.abs(row[f]);
          primaryF = f;
        }
      }
      for (let f = 0; f < nFactors; f++) {
        if (f === primaryF) continue;
        if (Math.abs(row[f]) < ps.crossLoadingThreshold) continue;
        const x1 = factorCx + ps.factorRx + 2, y1 = factorCy[f];
        const x2 = itemLeft - 2, y2 = itemCy[i];
        const dy = y2 - y1;
        const curv = ps.arrowCurvature;
        g.append("path").attr("d", `M${x1},${y1} C${x1 + ps.arrowSpan * curv},${y1 + dy * 0.05} ${x2 - ps.arrowSpan * curv * 0.85},${y2} ${x2},${y2}`).attr("fill", "none").attr("stroke", theme.textMuted).attr("stroke-width", 0.9).attr("stroke-dasharray", ps.crossLoadingDash).attr("opacity", ps.crossLoadingOpacity);
      }
    }
  }
  for (let fi = 0; fi < nFactors; fi++) {
    const cy = factorCy[fi];
    const color = fColor(fi, ps, theme);
    if (ps.factorGlow) {
      g.append("ellipse").attr("cx", factorCx).attr("cy", cy).attr("rx", ps.factorRx + 4).attr("ry", ps.factorRy + 3).attr("fill", color).attr("opacity", 0.06).attr("filter", "url(#fa-glow)");
    }
    if (ps.showShadows) {
      g.append("ellipse").attr("cx", factorCx).attr("cy", cy).attr("rx", ps.factorRx).attr("ry", ps.factorRy).attr("fill", "none").attr("filter", "url(#fa-shadow-med)").attr("stroke", "transparent");
    }
    g.append("ellipse").attr("cx", factorCx).attr("cy", cy).attr("rx", ps.factorRx).attr("ry", ps.factorRy).attr("fill", `url(#fa-fgrad-${fi})`).attr("stroke", color).attr("stroke-width", ps.factorStroke);
    if (ps.factorHighlight && ps.factorRx > 20 && ps.factorRy > 15) {
      g.append("ellipse").attr("cx", factorCx).attr("cy", cy - 4).attr("rx", ps.factorRx - 10).attr("ry", ps.factorRy - 10).attr("fill", "none").attr("stroke", "#fff").attr("stroke-width", 1).attr("opacity", 0.25);
    }
    g.append("text").attr("x", factorCx).attr("y", cy + 1).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", ps.factorFontSize).attr("font-weight", ps.factorFontWeight).attr("letter-spacing", "-0.4").attr("fill", theme.text).text(facLabels[fi] ?? `F${fi + 1}`);
  }
  for (let i = 0; i < nItems; i++) {
    const cy = itemCy[i];
    let primaryFactor = 0, primVal = 0;
    for (let f = 0; f < nFactors; f++) {
      if (Math.abs(data.loadings[i][f]) > primVal) {
        primVal = Math.abs(data.loadings[i][f]);
        primaryFactor = f;
      }
    }
    const accentColor = fColor(primaryFactor, ps, theme);
    if (ps.showShadows) {
      g.append("rect").attr("x", itemLeft).attr("y", cy - ps.itemHeight / 2).attr("width", ps.itemWidth).attr("height", ps.itemHeight).attr("rx", ps.itemRadius).attr("fill", theme.surface).attr("filter", "url(#fa-shadow-sm)");
    }
    g.append("rect").attr("x", itemLeft).attr("y", cy - ps.itemHeight / 2).attr("width", ps.itemWidth).attr("height", ps.itemHeight).attr("rx", ps.itemRadius).attr("fill", ps.itemGradient ? `url(#fa-igrad-${primaryFactor})` : theme.background).attr("stroke", accentColor).attr("stroke-width", ps.itemStroke).attr("stroke-opacity", 0.5);
    if (ps.itemAccentWidth > 0) {
      g.append("rect").attr("x", itemLeft).attr("y", cy - ps.itemHeight / 2).attr("width", ps.itemAccentWidth).attr("height", ps.itemHeight).attr("rx", Math.min(ps.itemAccentWidth, ps.itemRadius)).attr("fill", accentColor).attr("opacity", 0.85);
    }
    g.append("text").attr("x", itemLeft + ps.itemWidth / 2 + ps.itemAccentWidth / 2).attr("y", cy + 1).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamily).attr("font-size", ps.itemFontSize).attr("font-weight", ps.itemFontWeight).attr("fill", theme.text).text(varLabels[i] ?? `V${i + 1}`);
    g.append("rect").attr("x", itemLeft).attr("y", cy - ps.itemHeight / 2).attr("width", ps.itemWidth).attr("height", ps.itemHeight).attr("fill", "transparent").attr("cursor", "default").on("mouseover", (event) => {
      const rows = facLabels.map(
        (fl, fi) => formatTooltipRow(fl, data.loadings[i][fi] ?? 0)
      ).join("");
      showTooltip(event, [
        formatTooltipRow(varLabels[i] ?? `V${i + 1}`, ""),
        rows,
        formatTooltipRow("h\xB2", data.communalities[i] ?? 0)
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
  }
  if (showErrors) {
    for (let i = 0; i < nItems; i++) {
      const cy = itemCy[i];
      const uniq = data.uniqueness[i] ?? 0;
      const errArrX1 = errorCx - ps.errorRx - 1;
      const errArrX2 = itemRight + 2;
      g.append("line").attr("x1", errArrX1).attr("y1", cy).attr("x2", errArrX2).attr("y2", cy).attr("stroke", theme.axisLine).attr("stroke-width", ps.errorArrowWidth).attr("marker-end", "url(#fa-arr-err)");
      const errMidX = (errArrX1 + errArrX2) / 2;
      const oneLbl = g.append("text").attr("x", errMidX).attr("y", cy - 5).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", ps.errorFontSize).attr("font-weight", "600").attr("fill", theme.textMuted).text("1");
      if (ps.halo) oneLbl.attr("paint-order", "stroke").attr("stroke", ps.haloColor).attr("stroke-width", ps.haloWidth).attr("stroke-linejoin", "round");
      g.append("ellipse").attr("cx", errorCx).attr("cy", cy).attr("rx", ps.errorRx).attr("ry", ps.errorRy).attr("fill", theme.surface).attr("stroke", theme.axisLine).attr("stroke-width", ps.errorStroke);
      if (ps.errorSelfLoop) {
        const loopRight = errorCx + ps.errorRx + ps.errorSelfLoopSize;
        g.append("path").attr("d", `M${errorCx + ps.errorRx},${cy - 3} C${loopRight},${cy - 14} ${loopRight},${cy + 14} ${errorCx + ps.errorRx},${cy + 3}`).attr("fill", "none").attr("stroke", theme.axisLine).attr("stroke-width", ps.errorArrowWidth).attr("marker-end", "url(#fa-arr-err)");
        const errLbl = g.append("text").attr("x", loopRight + 2).attr("y", cy + ps.errorFontSize * 0.35).attr("text-anchor", "start").attr("font-family", theme.fontFamilyMono).attr("font-size", ps.errorFontSize - 0.5).attr("font-weight", "500").attr("fill", theme.textMuted).text(uniq.toFixed(2));
        if (ps.halo) errLbl.attr("paint-order", "stroke").attr("stroke", ps.haloColor).attr("stroke-width", ps.haloWidth).attr("stroke-linejoin", "round");
      }
      g.append("text").attr("x", errorCx).attr("y", cy + 1).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", ps.errorFontSize).attr("font-weight", "600").attr("fill", theme.textMuted).text(`e${i + 1}`);
    }
  }
  if (showFit) {
    const fit = data.fit;
    const fitColor = (good, ok) => good ? "#2e7d32" : ok ? "#e65100" : "#c62828";
    const table = d3.select(container).append("table").style("border-collapse", "collapse").style("margin", "12px auto 4px").style("font-family", theme.fontFamilyMono).style("font-size", "12px").style("color", theme.text).style("min-width", "420px");
    const thead = table.append("thead");
    const headerRow = thead.append("tr").style("border-bottom", `2px solid ${theme.text}`);
    const headers = ["Index", "Value", "90% CI", "Verdict"];
    headers.forEach((h) => {
      headerRow.append("th").style("padding", "6px 14px").style("text-align", h === "Index" ? "left" : "right").style("font-weight", "700").style("font-size", "11px").style("letter-spacing", "0.5px").style("color", theme.textMuted).text(h);
    });
    const tbody = table.append("tbody");
    const fitRows = [
      {
        label: "\u03C7\xB2",
        value: `${formatStat(fit.chiSq)} (df = ${fit.df})`,
        ci: "",
        verdict: `${formatP(fit.pValue)}`,
        color: theme.textAnnotation
      },
      {
        label: "RMSEA",
        value: formatStat(fit.rmsea, 3),
        ci: `[${formatStat(fit.rmseaCI[0], 3)}, ${formatStat(fit.rmseaCI[1], 3)}]`,
        verdict: fit.rmsea <= 0.05 ? "Good" : fit.rmsea <= 0.08 ? "Acceptable" : "Poor",
        color: fitColor(fit.rmsea <= 0.05, fit.rmsea <= 0.08)
      },
      {
        label: "CFI",
        value: formatStat(fit.cfi, 3),
        ci: "",
        verdict: fit.cfi >= 0.95 ? "Good" : fit.cfi >= 0.9 ? "Acceptable" : "Poor",
        color: fitColor(fit.cfi >= 0.95, fit.cfi >= 0.9)
      },
      {
        label: "TLI",
        value: formatStat(fit.tli, 3),
        ci: "",
        verdict: fit.tli >= 0.95 ? "Good" : fit.tli >= 0.9 ? "Acceptable" : "Poor",
        color: fitColor(fit.tli >= 0.95, fit.tli >= 0.9)
      },
      {
        label: "SRMR",
        value: formatStat(fit.srmr, 3),
        ci: "",
        verdict: fit.srmr <= 0.05 ? "Good" : fit.srmr <= 0.08 ? "Acceptable" : "Poor",
        color: fitColor(fit.srmr <= 0.05, fit.srmr <= 0.08)
      },
      {
        label: "AIC",
        value: formatStat(fit.aic, 1),
        ci: "",
        verdict: "",
        color: theme.textAnnotation
      },
      {
        label: "BIC",
        value: formatStat(fit.bic, 1),
        ci: "",
        verdict: "",
        color: theme.textAnnotation
      }
    ];
    fitRows.forEach((row, i) => {
      const tr = tbody.append("tr").style("border-bottom", i < fitRows.length - 1 ? `1px solid ${theme.gridLine}` : "none");
      tr.append("td").style("padding", "5px 14px").style("font-weight", "600").style("white-space", "nowrap").text(row.label);
      tr.append("td").style("padding", "5px 14px").style("text-align", "right").style("font-weight", "500").text(row.value);
      tr.append("td").style("padding", "5px 14px").style("text-align", "right").style("color", theme.textMuted).style("font-size", "11px").text(row.ci);
      const verdictCell = tr.append("td").style("padding", "5px 14px").style("text-align", "right").style("font-weight", "700").style("color", row.color);
      if (row.verdict) {
        verdictCell.append("span").style("display", "inline-block").style("width", "8px").style("height", "8px").style("border-radius", "50%").style("background", row.color).style("margin-right", "6px").style("vertical-align", "middle");
        verdictCell.append("span").style("vertical-align", "middle").text(row.verdict);
      }
    });
  }
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderCommunality(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const nItems = data.communalities.length;
  const W = config.width ?? 500;
  const H = config.height ?? Math.max(nItems * 28 + 120, 250);
  const margin = { top: theme.marginTop, right: 70, bottom: theme.marginBottom, left: 110 };
  const width = W - margin.left - margin.right;
  const height = H - margin.top - margin.bottom;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Communalities", data.formatted, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const varLabels = config.variableLabels ?? data.variableNames;
  const indices = Array.from({ length: nItems }, (_, i) => i);
  indices.sort((a, b) => (data.communalities[b] ?? 0) - (data.communalities[a] ?? 0));
  const names = indices.map((i) => varLabels[i] ?? `V${i + 1}`);
  const values = indices.map((i) => data.communalities[i] ?? 0);
  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width]).nice();
  const yScale = d3.scaleBand().domain(names).range([0, height]).padding(0.25);
  const commColor = (v) => {
    if (v < 0.3) return "#e15759";
    if (v < 0.5) return "#f28e2b";
    return "#59a14f";
  };
  g.selectAll(".grid").data(xScale.ticks(5)).join("line").attr("x1", (d) => xScale(d)).attr("x2", (d) => xScale(d)).attr("y1", 0).attr("y2", height).attr("stroke", theme.gridLine).attr("stroke-width", 1);
  g.append("line").attr("x1", xScale(0.4)).attr("x2", xScale(0.4)).attr("y1", 0).attr("y2", height).attr("stroke", theme.axisLine).attr("stroke-dasharray", "6,3").attr("stroke-width", 1.5);
  values.forEach((val, i) => {
    const y = yScale(names[i]) ?? 0;
    g.append("rect").attr("x", 0).attr("y", y).attr("width", xScale(val)).attr("height", yScale.bandwidth()).attr("fill", commColor(val)).attr("opacity", 0.8).attr("rx", 2).on("mouseover", (event) => {
      const origIdx = indices[i];
      showTooltip(event, [
        formatTooltipRow("Variable", varLabels[origIdx] ?? `V${origIdx + 1}`),
        formatTooltipRow("Communality", val),
        formatTooltipRow("Uniqueness", data.uniqueness[origIdx] ?? 0)
      ].join(""), theme);
    }).on("mouseout", hideTooltip);
    g.append("text").attr("x", xScale(val) + 4).attr("y", y + yScale.bandwidth() / 2 + 4).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall - 1).attr("fill", theme.text).text(val.toFixed(2));
  });
  g.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(5)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("g").call(d3.axisLeft(yScale)).selectAll("text").attr("fill", theme.text).attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize);
  g.append("text").attr("x", width / 2).attr("y", height + 44).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text("Communality (h\xB2)");
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderFactorCorrelation(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const Phi = data.factorCorrelations;
  const k = data.nFactors;
  const facLabels = config.factorLabels ?? data.factorNames;
  let isOrthogonal = true;
  for (let i = 0; i < k && isOrthogonal; i++) {
    for (let j = 0; j < k && isOrthogonal; j++) {
      if (i !== j && Math.abs(Phi[i][j]) > 1e-3) isOrthogonal = false;
    }
  }
  const cellSize = 60;
  const W = config.width ?? k * cellSize + 140;
  const H = config.height ?? k * cellSize + 120;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Factor Correlations (\u03A6)", data.rotation === "varimax" ? "Orthogonal rotation \u2014 \u03A6 = I" : "", W, theme);
  if (isOrthogonal) {
    svg.append("text").attr("x", W / 2).attr("y", H / 2).attr("text-anchor", "middle").attr("font-size", 14).attr("fill", theme.textMuted).text("Orthogonal rotation: all factor correlations = 0");
    return;
  }
  const margin = { top: 60, right: 20, bottom: 20, left: 80 };
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const colorScale = d3.scaleSequential().domain([-1, 1]).interpolator(d3.interpolateRdBu);
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) {
      const val = Phi[i][j];
      g.append("rect").attr("x", j * cellSize).attr("y", i * cellSize).attr("width", cellSize - 2).attr("height", cellSize - 2).attr("rx", 3).attr("fill", i === j ? theme.surface : colorScale(val)).on("mouseover", (event) => {
        if (i !== j) {
          showTooltip(event, [
            formatTooltipRow(`${facLabels[i]} \xD7 ${facLabels[j]}`, ""),
            formatTooltipRow("r", val)
          ].join(""), theme);
        }
      }).on("mouseout", hideTooltip);
      if (i !== j) {
        g.append("text").attr("x", j * cellSize + cellSize / 2).attr("y", i * cellSize + cellSize / 2 + 4).attr("text-anchor", "middle").attr("font-family", theme.fontFamilyMono).attr("font-size", 12).attr("fill", Math.abs(val) > 0.6 ? "#fff" : theme.text).text(val.toFixed(2));
      }
    }
  }
  facLabels.forEach((lbl, i) => {
    g.append("text").attr("x", -6).attr("y", i * cellSize + cellSize / 2 + 4).attr("text-anchor", "end").attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(lbl);
    g.append("text").attr("x", i * cellSize + cellSize / 2).attr("y", -8).attr("text-anchor", "middle").attr("font-size", theme.fontSizeSmall).attr("fill", theme.text).text(lbl);
  });
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}
function renderFitIndices(d3, container, data, config) {
  const theme = config.theme ?? DEFAULT_THEME;
  const W = config.width ?? 500, H = config.height ?? 340;
  const margin = { top: theme.marginTop, right: 40, bottom: 50, left: 80 };
  const width = W - margin.left - margin.right;
  container.innerHTML = "";
  applyTheme(container, theme);
  const svg = d3.select(container).append("svg").attr("width", W).attr("height", H).style("background", theme.background);
  addSubtitle(svg, config.title ?? "Model Fit Indices", data.formatted, W, theme);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const fit = data.fit;
  const gauges = [
    { label: "RMSEA", value: fit.rmsea, range: [0, 0.15], good: 0.05, ok: 0.08, reversed: true, ci: fit.rmseaCI },
    { label: "CFI", value: fit.cfi, range: [0.8, 1], good: 0.95, ok: 0.9, reversed: false },
    { label: "TLI", value: fit.tli, range: [0.8, 1], good: 0.95, ok: 0.9, reversed: false },
    { label: "SRMR", value: fit.srmr, range: [0, 0.12], good: 0.05, ok: 0.08, reversed: true }
  ];
  const barH = 24;
  const gap = 50;
  const startY = 10;
  gauges.forEach((gauge, i) => {
    const y = startY + i * gap;
    const xScale = d3.scaleLinear().domain([gauge.range[0], gauge.range[1]]).range([0, width]).clamp(true);
    if (gauge.reversed) {
      g.append("rect").attr("x", 0).attr("y", y).attr("width", xScale(gauge.good)).attr("height", barH).attr("fill", "#59a14f").attr("opacity", 0.15).attr("rx", 4);
      g.append("rect").attr("x", xScale(gauge.good)).attr("y", y).attr("width", xScale(gauge.ok) - xScale(gauge.good)).attr("height", barH).attr("fill", "#f28e2b").attr("opacity", 0.15);
      g.append("rect").attr("x", xScale(gauge.ok)).attr("y", y).attr("width", width - xScale(gauge.ok)).attr("height", barH).attr("fill", "#e15759").attr("opacity", 0.15).attr("rx", 4);
    } else {
      g.append("rect").attr("x", 0).attr("y", y).attr("width", xScale(gauge.ok)).attr("height", barH).attr("fill", "#e15759").attr("opacity", 0.15).attr("rx", 4);
      g.append("rect").attr("x", xScale(gauge.ok)).attr("y", y).attr("width", xScale(gauge.good) - xScale(gauge.ok)).attr("height", barH).attr("fill", "#f28e2b").attr("opacity", 0.15);
      g.append("rect").attr("x", xScale(gauge.good)).attr("y", y).attr("width", width - xScale(gauge.good)).attr("height", barH).attr("fill", "#59a14f").attr("opacity", 0.15).attr("rx", 4);
    }
    g.append("line").attr("x1", xScale(gauge.good)).attr("x2", xScale(gauge.good)).attr("y1", y).attr("y2", y + barH).attr("stroke", theme.axisLine).attr("stroke-dasharray", "2,2").attr("stroke-width", 1);
    g.append("line").attr("x1", xScale(gauge.ok)).attr("x2", xScale(gauge.ok)).attr("y1", y).attr("y2", y + barH).attr("stroke", theme.axisLine).attr("stroke-dasharray", "2,2").attr("stroke-width", 1);
    if ("ci" in gauge && gauge.ci) {
      const ciLo = gauge.ci[0];
      const ciHi = gauge.ci[1];
      g.append("rect").attr("x", xScale(ciLo)).attr("y", y + barH / 2 - 3).attr("width", xScale(ciHi) - xScale(ciLo)).attr("height", 6).attr("fill", getColor(0, theme)).attr("opacity", 0.3).attr("rx", 2);
    }
    const valColor = gauge.reversed ? gauge.value <= gauge.good ? "#59a14f" : gauge.value <= gauge.ok ? "#f28e2b" : "#e15759" : gauge.value >= gauge.good ? "#59a14f" : gauge.value >= gauge.ok ? "#f28e2b" : "#e15759";
    g.append("circle").attr("cx", xScale(gauge.value)).attr("cy", y + barH / 2).attr("r", 6).attr("fill", valColor).attr("stroke", "#fff").attr("stroke-width", 2);
    g.append("text").attr("x", -6).attr("y", y + barH / 2 + 4).attr("text-anchor", "end").attr("font-size", theme.fontSizeSmall).attr("font-weight", "600").attr("fill", theme.text).text(gauge.label);
    g.append("text").attr("x", width + 6).attr("y", y + barH / 2 + 4).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", valColor).attr("font-weight", "600").text(gauge.value.toFixed(3));
  });
  const textY = startY + gauges.length * gap + 10;
  g.append("text").attr("x", 0).attr("y", textY).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textAnnotation).text(`\u03C7\xB2(${fit.df}) = ${formatStat(fit.chiSq)}, ${formatP(fit.pValue)}`);
  g.append("text").attr("x", 0).attr("y", textY + 16).attr("font-family", theme.fontFamilyMono).attr("font-size", theme.fontSizeSmall).attr("fill", theme.textAnnotation).text(`AIC = ${formatStat(fit.aic, 1)}, BIC = ${formatStat(fit.bic, 1)}`);
  if (config.caption) addCaption(svg, config.caption, W, H, theme);
}

// src/viz/export.ts
function exportSVG(container, filename = "carm-plot.svg") {
  const svgEl = container.querySelector("svg");
  if (!svgEl) throw new Error("exportSVG: no SVG element found in container");
  if (!svgEl.hasAttribute("xmlns")) {
    svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  }
  const serializer = new XMLSerializer();
  const svgStr = serializer.serializeToString(svgEl);
  const blob = new Blob([svgStr], { type: "image/svg+xml;charset=utf-8" });
  triggerDownload(URL.createObjectURL(blob), filename);
}
function exportPNG(container, filename = "carm-plot.png", dpi = 300) {
  const svgEl = container.querySelector("svg");
  if (!svgEl) throw new Error("exportPNG: no SVG element found in container");
  const svgW = svgEl.viewBox.baseVal.width || svgEl.clientWidth;
  const svgH = svgEl.viewBox.baseVal.height || svgEl.clientHeight;
  const scale = dpi / 96;
  const canvas = document.createElement("canvas");
  canvas.width = svgW * scale;
  canvas.height = svgH * scale;
  const ctx = canvas.getContext("2d");
  ctx.scale(scale, scale);
  if (!svgEl.hasAttribute("xmlns")) svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const svgStr = new XMLSerializer().serializeToString(svgEl);
  const blob = new Blob([svgStr], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => {
      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      canvas.toBlob((pngBlob) => {
        if (!pngBlob) {
          reject(new Error("Canvas toBlob failed"));
          return;
        }
        triggerDownload(URL.createObjectURL(pngBlob), filename);
        resolve();
      }, "image/png");
    };
    img.onerror = () => {
      URL.revokeObjectURL(url);
      reject(new Error("SVG to Image conversion failed"));
    };
    img.src = url;
  });
}
function triggerDownload(url, filename) {
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1e3);
}
//# sourceMappingURL=index.cjs.map
