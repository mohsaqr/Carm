'use strict';

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

exports.adjustPValues = adjustPValues;
exports.betaFn = betaFn;
exports.chiSqCDF = chiSqCDF;
exports.chiSqPValue = chiSqPValue;
exports.chiSqQuantile = chiSqQuantile;
exports.clamp = clamp;
exports.cov = cov;
exports.fDistCDF = fDistCDF;
exports.fDistPValue = fDistPValue;
exports.gamma = gamma;
exports.incompleteBeta = incompleteBeta;
exports.incompleteGamma = incompleteGamma;
exports.logBeta = logBeta;
exports.logGamma = logGamma;
exports.mean = mean;
exports.median = median;
exports.nelderMead = nelderMead;
exports.normalCDF = normalCDF;
exports.normalQuantile = normalQuantile;
exports.quantile = quantile;
exports.rank = rank;
exports.roundTo = roundTo;
exports.sd = sd;
exports.se = se;
exports.sortAsc = sortAsc;
exports.ss = ss;
exports.tDistCDF = tDistCDF;
exports.tDistPValue = tDistPValue;
exports.tDistQuantile = tDistQuantile;
exports.variance = variance;
//# sourceMappingURL=chunk-AVJH2SAO.cjs.map
//# sourceMappingURL=chunk-AVJH2SAO.cjs.map