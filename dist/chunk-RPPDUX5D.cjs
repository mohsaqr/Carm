'use strict';

var chunkFSSEZIKV_cjs = require('./chunk-FSSEZIKV.cjs');
var chunkMX4OLB7V_cjs = require('./chunk-MX4OLB7V.cjs');

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
  const sorted = chunkMX4OLB7V_cjs.sortAsc(x);
  const n = sorted.length;
  const trim = Math.floor(n * alpha);
  const trimmed = sorted.slice(trim, n - trim);
  if (trimmed.length === 0) throw new Error("trimmedMean: too much trimming, no data remains");
  return chunkMX4OLB7V_cjs.mean(trimmed);
}
function skewness(x) {
  const n = x.length;
  if (n < 3) throw new Error("skewness: need at least 3 observations");
  const m = chunkMX4OLB7V_cjs.mean(x);
  const s = chunkMX4OLB7V_cjs.sd(x);
  if (s === 0) return 0;
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 3, 0);
  return n / ((n - 1) * (n - 2)) * sum;
}
function kurtosis(x) {
  const n = x.length;
  if (n < 4) throw new Error("kurtosis: need at least 4 observations");
  const m = chunkMX4OLB7V_cjs.mean(x);
  const s = chunkMX4OLB7V_cjs.sd(x);
  if (s === 0) return 0;
  const sum = x.reduce((acc, v) => acc + ((v - m) / s) ** 4, 0);
  const k1 = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * sum;
  const k2 = 3 * (n - 1) ** 2 / ((n - 2) * (n - 3));
  return k1 - k2;
}
function ciMean(x, ciLevel = 0.95) {
  const n = x.length;
  if (n < 2) throw new Error("ciMean: need at least 2 observations");
  const m = chunkMX4OLB7V_cjs.mean(x);
  const s = chunkMX4OLB7V_cjs.se(x);
  const t = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, n - 1);
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
  const sorted = chunkMX4OLB7V_cjs.sortAsc(x);
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
      a[i] = chunkMX4OLB7V_cjs.normalQuantile((i + 1 - 0.375) / (n + 0.25));
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
  const sst = chunkMX4OLB7V_cjs.variance(sorted) * (n - 1);
  const W = sst > 0 ? Math.min(1, w1 * w1 / sst) : 1;
  const pValue = shapiroWilkPValue(W, n);
  return { statistic: chunkMX4OLB7V_cjs.roundTo(W, 4), pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4) };
}
function polynomialEval(coeffs, x) {
  return coeffs.reduce((acc, c) => acc * x + c, 0);
}
function shapiroWilkPValue(W, n) {
  let y = Math.log(1 - W);
  let z;
  if (n <= 11) {
    const gamma = polynomialEval([0.459, -2.273], 1 / n);
    if (y >= gamma) return 5e-7;
    y = -Math.log(gamma - y);
    const mu = polynomialEval([-1.2725, 1.0521, -0.0895], 1 / n);
    const sigma = Math.exp(polynomialEval([-6714e-7, 0.025054, -0.6714, 0.724], 1 / n));
    z = (y - mu) / sigma;
  } else {
    const mu = polynomialEval([38915e-7, -0.083751, -0.31082, -1.5861], Math.log(n));
    const sigma = Math.exp(polynomialEval([30302e-7, -0.082676, -0.4803], Math.log(n)));
    z = (y - mu) / sigma;
  }
  return Math.max(0, Math.min(1, 1 - chunkMX4OLB7V_cjs.normalCDF(z)));
}
function describe(x, ciLevel = 0.95) {
  if (x.length === 0) throw new Error("describe: empty array");
  const n = x.length;
  const m = chunkMX4OLB7V_cjs.mean(x);
  const med = chunkMX4OLB7V_cjs.median(x);
  const modes = mode(x);
  const tm = trimmedMean(x, 0.05);
  const s = chunkMX4OLB7V_cjs.sd(x);
  const sem = chunkMX4OLB7V_cjs.se(x);
  const v = chunkMX4OLB7V_cjs.variance(x);
  const sorted = chunkMX4OLB7V_cjs.sortAsc(x);
  const mn = sorted[0];
  const mx = sorted[n - 1];
  const q1 = chunkMX4OLB7V_cjs.quantile(x, 0.25);
  const q3 = chunkMX4OLB7V_cjs.quantile(x, 0.75);
  const iqr = q3 - q1;
  const skew = n >= 3 ? skewness(x) : 0;
  const kurt = n >= 4 ? kurtosis(x) : 0;
  const ci = ciMean(x, ciLevel);
  const sw = n >= 3 ? shapiroWilk(x) : { statistic: NaN, pValue: NaN };
  const ciPct = Math.round(ciLevel * 100);
  const formatted = [
    `n = ${n}`,
    `M = ${chunkMX4OLB7V_cjs.roundTo(m, 2)}, SD = ${chunkMX4OLB7V_cjs.roundTo(s, 2)}, SE = ${chunkMX4OLB7V_cjs.roundTo(sem, 2)}`,
    `Mdn = ${chunkMX4OLB7V_cjs.roundTo(med, 2)}, IQR = ${chunkMX4OLB7V_cjs.roundTo(iqr, 2)}`,
    `Skew = ${chunkMX4OLB7V_cjs.roundTo(skew, 2)}, Kurt = ${chunkMX4OLB7V_cjs.roundTo(kurt, 2)}`,
    `${ciPct}% CI [${chunkMX4OLB7V_cjs.roundTo(ci[0], 2)}, ${chunkMX4OLB7V_cjs.roundTo(ci[1], 2)}]`,
    n >= 3 ? `Shapiro-Wilk W = ${chunkMX4OLB7V_cjs.roundTo(sw.statistic, 3)}, ${chunkMX4OLB7V_cjs.formatP(sw.pValue)}` : ""
  ].filter(Boolean).join("; ");
  return {
    n,
    mean: chunkMX4OLB7V_cjs.roundTo(m, 6),
    median: chunkMX4OLB7V_cjs.roundTo(med, 6),
    mode: modes,
    trimmedMean: chunkMX4OLB7V_cjs.roundTo(tm, 6),
    sd: chunkMX4OLB7V_cjs.roundTo(s, 6),
    se: chunkMX4OLB7V_cjs.roundTo(sem, 6),
    variance: chunkMX4OLB7V_cjs.roundTo(v, 6),
    min: mn,
    max: mx,
    range: mx - mn,
    iqr: chunkMX4OLB7V_cjs.roundTo(iqr, 6),
    q1: chunkMX4OLB7V_cjs.roundTo(q1, 6),
    q3: chunkMX4OLB7V_cjs.roundTo(q3, 6),
    skewness: chunkMX4OLB7V_cjs.roundTo(skew, 6),
    kurtosis: chunkMX4OLB7V_cjs.roundTo(kurt, 6),
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
  const m1 = chunkMX4OLB7V_cjs.mean(x1), m2 = chunkMX4OLB7V_cjs.mean(x2);
  const v1 = chunkMX4OLB7V_cjs.variance(x1), v2 = chunkMX4OLB7V_cjs.variance(x2);
  const sdPooled = Math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2));
  if (sdPooled === 0) return { value: 0, name: "Cohen's d", interpretation: "negligible" };
  const d = (m1 - m2) / sdPooled;
  return {
    value: d,
    name: "Cohen's d",
    interpretation: chunkMX4OLB7V_cjs.interpretCohensD(d)
  };
}
function cohensDPaired(diffs) {
  if (diffs.length < 2) throw new Error("cohensDPaired: need at least 2 differences");
  const m = chunkMX4OLB7V_cjs.mean(diffs);
  const s = chunkMX4OLB7V_cjs.sd(diffs);
  const d = s === 0 ? 0 : m / s;
  return {
    value: d,
    name: "Cohen's d",
    interpretation: chunkMX4OLB7V_cjs.interpretCohensD(d)
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
    interpretation: chunkMX4OLB7V_cjs.interpretCohensD(g)
  };
}
function etaSquared(ssBetween, ssTotal) {
  if (ssTotal <= 0) return { value: 0, name: "\u03B7\xB2", interpretation: "negligible" };
  const eta2 = Math.max(0, Math.min(1, ssBetween / ssTotal));
  return {
    value: eta2,
    name: "\u03B7\xB2",
    interpretation: chunkMX4OLB7V_cjs.interpretEtaSq(eta2)
  };
}
function omegaSquared(ssBetween, ssTotal, dfBetween, msWithin) {
  const denom = ssTotal + msWithin;
  if (denom <= 0) return { value: 0, name: "\u03C9\xB2", interpretation: "negligible" };
  const omega2 = Math.max(0, (ssBetween - dfBetween * msWithin) / denom);
  return {
    value: omega2,
    name: "\u03C9\xB2",
    interpretation: chunkMX4OLB7V_cjs.interpretEtaSq(omega2)
  };
}
function rankBiserial(U, n1, n2) {
  const r = 1 - 2 * U / (n1 * n2);
  return {
    value: r,
    name: "r (rank-biserial)",
    interpretation: chunkMX4OLB7V_cjs.interpretR(r)
  };
}
function rankBiserialWilcoxon(T, n) {
  const maxT = n * (n + 1) / 2;
  const r = maxT > 0 ? T / maxT * 2 - 1 : 0;
  return {
    value: r,
    name: "r (rank-biserial)",
    interpretation: chunkMX4OLB7V_cjs.interpretR(r)
  };
}
function etaSquaredKW(H, k, n) {
  const eta2 = Math.max(0, (H - k + 1) / (n - k));
  return {
    value: eta2,
    name: "\u03B7\xB2_H",
    interpretation: chunkMX4OLB7V_cjs.interpretEtaSq(eta2)
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
      const num = yatesCorrection ? Math.max(0, Math.abs(o - e) - 0.5) : o - e;
      chiSq += num * num / e;
    }
  }
  const df = (R - 1) * (C - 1);
  const pValue = chunkMX4OLB7V_cjs.chiSqPValue(chiSq, df);
  const minDim = Math.min(R, C) - 1;
  const cramersV = Math.sqrt(chiSq / (n * Math.max(1, minDim)));
  const effectSize = {
    value: cramersV,
    name: "Cram\xE9r's V",
    interpretation: chunkMX4OLB7V_cjs.interpretCramerV(cramersV, df)
  };
  const ci = [NaN, NaN];
  const formatted = chunkMX4OLB7V_cjs.formatChiSq(chiSq, df, pValue, cramersV, "V");
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
  const z = chunkMX4OLB7V_cjs.normalQuantile(0.975);
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
    formatted: `OR = ${oddsRatio.toFixed(2)}, ${chunkMX4OLB7V_cjs.formatP(pValue)}, 95% CI [${ci[0].toFixed(2)}, ${ci[1].toFixed(2)}]`
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
  const pValue = chunkMX4OLB7V_cjs.chiSqPValue(chiSq, df);
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
    formatted: chunkMX4OLB7V_cjs.formatChiSq(chiSq, df, pValue, w, "w")
  };
}

// src/stats/comparison.ts
function tTestIndependent(x1, x2, equalVariances = false, ciLevel = 0.95, alternative = "two.sided") {
  if (x1.length < 2 || x2.length < 2) throw new Error("tTestIndependent: need at least 2 per group");
  const n1 = x1.length, n2 = x2.length;
  const m1 = chunkMX4OLB7V_cjs.mean(x1), m2 = chunkMX4OLB7V_cjs.mean(x2);
  const v1 = chunkMX4OLB7V_cjs.variance(x1), v2 = chunkMX4OLB7V_cjs.variance(x2);
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
  const pFull = chunkMX4OLB7V_cjs.tDistPValue(t, df);
  const pValue = alternative === "two.sided" ? pFull : alternative === "less" ? t < 0 ? pFull / 2 : 1 - pFull / 2 : t > 0 ? pFull / 2 : 1 - pFull / 2;
  const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const diff = m1 - m2;
  const ci = [diff - tCrit * se2, diff + tCrit * se2];
  const effectSize = cohensD(x1, x2);
  const formatted = chunkMX4OLB7V_cjs.formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel);
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
  const mDiff = chunkMX4OLB7V_cjs.mean(diffs);
  const seDiff = chunkMX4OLB7V_cjs.se(diffs);
  const df = n - 1;
  const t = seDiff === 0 ? 0 : mDiff / seDiff;
  const pValue = chunkMX4OLB7V_cjs.tDistPValue(t, df);
  const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const ci = [mDiff - tCrit * seDiff, mDiff + tCrit * seDiff];
  const effectSize = cohensDPaired(diffs);
  const formatted = chunkMX4OLB7V_cjs.formatTTest(t, df, pValue, effectSize.value, "d", ci, ciLevel);
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
  const grandMean = chunkMX4OLB7V_cjs.mean(allValues);
  let ssBetween = 0, ssWithin = 0;
  const groupStats = groups.map((g) => {
    const gm = chunkMX4OLB7V_cjs.mean(g.values);
    const gn = g.values.length;
    ssBetween += gn * (gm - grandMean) ** 2;
    ssWithin += g.values.reduce((s, v) => s + (v - gm) ** 2, 0);
    return { label: g.label, n: gn, mean: gm, sd: chunkMX4OLB7V_cjs.sd(g.values) };
  });
  const ssTotal = ssBetween + ssWithin;
  const dfBetween = k - 1;
  const dfWithin = n - k;
  const msBetween = ssBetween / dfBetween;
  const msWithin = ssWithin / dfWithin;
  const F = msWithin === 0 ? Infinity : msBetween / msWithin;
  const pValue = chunkMX4OLB7V_cjs.fDistPValue(F, dfBetween, dfWithin);
  const omega = omegaSquared(ssBetween, ssTotal, dfBetween, msWithin);
  const formatted = chunkMX4OLB7V_cjs.formatANOVA(F, dfBetween, dfWithin, pValue, omega.value, "\u03C9\xB2");
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
  const allCombined = chunkMX4OLB7V_cjs.rank(combined.map((d) => d.v));
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
  const formatted = chunkMX4OLB7V_cjs.formatMannWhitney(U1, pValue, effect.value);
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
  const erf = 1 - poly * Math.exp(-x * x);
  return 0.5 * (1 + (z >= 0 ? erf : -erf));
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
  const ranks_ = chunkMX4OLB7V_cjs.rank(absDiffs);
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
    formatted: `V = ${W}, ${chunkMX4OLB7V_cjs.formatP(pValue)}, r = ${effect.value.toFixed(2)}`
  };
}
function kruskalWallis(groups) {
  if (groups.length < 2) throw new Error("kruskalWallis: need at least 2 groups");
  const k = groups.length;
  const allValues = groups.flatMap((g) => [...g.values]);
  const n = allValues.length;
  const allRanks_ = chunkMX4OLB7V_cjs.rank(allValues);
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
  const pValue = chunkMX4OLB7V_cjs.chiSqPValue(H, df);
  const effect = etaSquaredKW(H, k, n);
  const formatted = chunkMX4OLB7V_cjs.formatKruskalWallis(H, df, pValue, effect.value);
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
    const rowRanks = chunkMX4OLB7V_cjs.rank(row);
    for (const r of rowRanks) sumRankSq += r * r;
  }
  const colRankSums = Array.from(
    { length: k },
    (_, j) => data.reduce((s, row) => {
      const rowRanks = chunkMX4OLB7V_cjs.rank(row);
      return s + (rowRanks[j] ?? 0);
    }, 0)
  );
  const Rj2sum = colRankSums.reduce((s, r) => s + r * r, 0);
  const chi2 = 12 / (n * k * (k + 1)) * Rj2sum - 3 * n * (k + 1);
  const df = k - 1;
  const pValue = chunkMX4OLB7V_cjs.chiSqPValue(chi2, df);
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
    formatted: `\u03C7\xB2_F(${df}) = ${chi2.toFixed(2)}, ${chunkMX4OLB7V_cjs.formatP(pValue)}, W = ${w.toFixed(2)}`
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
      const m1 = chunkMX4OLB7V_cjs.mean(g1.values);
      const m2 = chunkMX4OLB7V_cjs.mean(g2.values);
      const diff = m1 - m2;
      const se2 = Math.sqrt(msWithin / 2 * (1 / n1 + 1 / n2));
      const q = se2 === 0 ? 0 : Math.abs(diff) / se2;
      const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - alpha / (k * (k - 1)), dfWithin);
      const pValue = pValueStudentizedRange(q, k);
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
function pValueStudentizedRange(q, k, _df) {
  const t = q / Math.SQRT2;
  const pOne = 2 * (1 - normCDF(t));
  return Math.min(1, k * (k - 1) / 2 * pOne);
}
function normCDF(z) {
  const x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const erf = 1 - poly * Math.exp(-x * x);
  return 0.5 * (1 + (z >= 0 ? erf : -erf));
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
      const m1 = chunkMX4OLB7V_cjs.mean(g1.values);
      const m2 = chunkMX4OLB7V_cjs.mean(g2.values);
      const v1 = chunkMX4OLB7V_cjs.variance(g1.values);
      const v2 = chunkMX4OLB7V_cjs.variance(g2.values);
      const diff = m1 - m2;
      const se2 = Math.sqrt(v1 / n1 + v2 / n2);
      const q = se2 === 0 ? 0 : Math.abs(diff) / se2;
      const dfNum = (v1 / n1 + v2 / n2) ** 2;
      const dfDen = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1);
      const df = dfDen > 0 ? dfNum / dfDen : n1 + n2 - 2;
      const pValue = pValueStudentizedRange(q, k);
      const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - alpha / (k * (k - 1)), df);
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
  const allRanks = chunkMX4OLB7V_cjs.rank(allValues);
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
  const rawPValues = [];
  const pairs = [];
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const diff = (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0);
      const se2 = Math.sqrt(
        (n * (n + 1) / 12 - tieAdj) * (1 / (groupNs[i] ?? 1) + 1 / (groupNs[j] ?? 1))
      );
      const z = se2 === 0 ? 0 : diff / se2;
      const p = 2 * (1 - normCDF(Math.abs(z)));
      rawPValues.push(p);
      pairs.push({ i, j });
    }
  }
  const adjPValues = chunkMX4OLB7V_cjs.adjustPValues(rawPValues, method);
  return pairs.map(({ i, j }, idx) => ({
    group1: groups[i].label,
    group2: groups[j].label,
    meanDiff: (groupRankMeans[i] ?? 0) - (groupRankMeans[j] ?? 0),
    se: 0,
    statistic: 0,
    pValue: rawPValues[idx],
    pValueAdj: adjPValues[idx],
    ci: [NaN, NaN],
    significant: (adjPValues[idx] ?? 1) < 0.05
  }));
}

// src/stats/correlation.ts
function pearsonCorrelation(x, y, ciLevel = 0.95) {
  if (x.length !== y.length) throw new Error("pearsonCorrelation: arrays must have equal length");
  const n = x.length;
  if (n < 3) throw new Error("pearsonCorrelation: need at least 3 observations");
  const sdX = chunkMX4OLB7V_cjs.sd(x), sdY = chunkMX4OLB7V_cjs.sd(y);
  if (sdX === 0 || sdY === 0) throw new Error("pearsonCorrelation: zero variance in input");
  const r = chunkMX4OLB7V_cjs.cov(x, y) / (sdX * sdY);
  const rClamped = Math.max(-1, Math.min(1, r));
  const df = n - 2;
  const t = Math.abs(rClamped) === 1 ? Infinity : rClamped * Math.sqrt(df / (1 - rClamped * rClamped));
  const pValue = chunkMX4OLB7V_cjs.tDistPValue(t, df);
  const ci = fisherZCI(rClamped, n, ciLevel);
  return {
    testName: "Pearson r",
    statistic: chunkMX4OLB7V_cjs.roundTo(rClamped, 4),
    df,
    pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4),
    effectSize: {
      value: rClamped,
      name: "Pearson r",
      interpretation: chunkMX4OLB7V_cjs.interpretR(rClamped)
    },
    ci,
    ciLevel,
    n,
    formatted: chunkMX4OLB7V_cjs.formatCorrelation(rClamped, df, pValue, ci, "r", ciLevel)
  };
}
function fisherZCI(r, n, ciLevel) {
  const z = Math.log((1 + r) / (1 - r)) / 2;
  const se2 = 1 / Math.sqrt(n - 3);
  const zCrit = chunkMX4OLB7V_cjs.normalQuantile(1 - (1 - ciLevel) / 2);
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
  const rx = chunkMX4OLB7V_cjs.rank(x), ry = chunkMX4OLB7V_cjs.rank(y);
  const rhoResult = pearsonCorrelation(rx, ry, ciLevel);
  return {
    ...rhoResult,
    testName: "Spearman's \u03C1",
    effectSize: {
      ...rhoResult.effectSize,
      name: "Spearman's \u03C1"
    },
    formatted: chunkMX4OLB7V_cjs.formatCorrelation(rhoResult.statistic, typeof rhoResult.df === "number" ? rhoResult.df : 0, rhoResult.pValue, rhoResult.ci, "\u03C1", ciLevel)
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
  const varTau = 2 * (2 * n + 5) / (9 * n * (n - 1));
  const z = tau / Math.sqrt(varTau);
  const pValue = 2 * (1 - chunkMX4OLB7V_cjs.normalCDF(Math.abs(z)));
  const ci = fisherZCI(tau, n, ciLevel);
  const df = n - 2;
  return {
    testName: "Kendall's \u03C4",
    statistic: chunkMX4OLB7V_cjs.roundTo(tau, 4),
    df,
    pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4),
    effectSize: {
      value: tau,
      name: "Kendall's \u03C4",
      interpretation: chunkMX4OLB7V_cjs.interpretR(tau)
    },
    ci,
    ciLevel,
    n,
    formatted: chunkMX4OLB7V_cjs.formatCorrelation(tau, df, pValue, ci, "\u03C4", ciLevel)
  };
}
function partialCorrelation(x, y, controls) {
  if (x.length !== y.length) throw new Error("partialCorrelation: arrays must have equal length");
  const xRes = residualize(x, controls);
  const yRes = residualize(y, controls);
  return pearsonCorrelation(xRes, yRes);
}
function residualize(y, predictors) {
  if (predictors.length === 0) return [...y];
  const n = y.length;
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((p) => p[i] ?? 0)])
  );
  const Xt = X.transpose();
  const XtX = Xt.multiply(X);
  const XtY = Xt.multiply(chunkFSSEZIKV_cjs.Matrix.colVec(y));
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
        return corrFn(data[i], data[j]).statistic;
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
  const XtY = Xt.multiply(chunkFSSEZIKV_cjs.Matrix.colVec(y));
  const XtXInv = XtX.inverse();
  const betaM = XtXInv.multiply(XtY);
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0));
  const fitted = Array.from({ length: n }, (_, i) => {
    let val = 0;
    for (let j = 0; j < p; j++) val += X.get(i, j) * (beta[j] ?? 0);
    return val;
  });
  const residuals = y.map((v, i) => v - (fitted[i] ?? 0));
  const yMean = chunkMX4OLB7V_cjs.mean(y);
  const ss_res = residuals.reduce((s, r) => s + r * r, 0);
  const ss_tot = y.reduce((s, v) => s + (v - yMean) ** 2, 0);
  const r2 = ss_tot > 0 ? Math.max(0, 1 - ss_res / ss_tot) : 0;
  const adjR2 = 1 - (1 - r2) * (n - 1) / (n - p);
  const dfRes = n - p;
  if (dfRes <= 0) throw new Error("fitOLS: not enough degrees of freedom");
  const sigma2 = ss_res / dfRes;
  const covBeta = XtXInv.scale(sigma2);
  const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, dfRes);
  const coefficients = beta.map((b, i) => {
    const se2 = Math.sqrt(Math.max(0, covBeta.get(i, i)));
    const t = se2 === 0 ? 0 : b / se2;
    const pVal = chunkMX4OLB7V_cjs.tDistPValue(t, dfRes);
    const ci = [b - tCrit * se2, b + tCrit * se2];
    return {
      name: coefNames[i] ?? `\u03B2${i}`,
      estimate: chunkMX4OLB7V_cjs.roundTo(b, 6),
      se: chunkMX4OLB7V_cjs.roundTo(se2, 6),
      tValue: chunkMX4OLB7V_cjs.roundTo(t, 4),
      pValue: chunkMX4OLB7V_cjs.roundTo(pVal, 4),
      ci
    };
  });
  const dfModel = p - 1;
  const ss_reg = ss_tot - ss_res;
  const F = sigma2 === 0 || dfModel === 0 ? 0 : ss_reg / dfModel / sigma2;
  const fPValue = chunkMX4OLB7V_cjs.fDistPValue(F, dfModel, dfRes);
  const rssSafe = Math.max(ss_res, 1e-15);
  const logLik = -n / 2 * (Math.log(2 * Math.PI) + Math.log(rssSafe / n) + 1);
  const aic = -2 * logLik + 2 * (p + 1);
  const bic = -2 * logLik + Math.log(n) * (p + 1);
  const formatted = chunkMX4OLB7V_cjs.formatRegression(r2, adjR2, F, dfModel, dfRes, fPValue);
  return {
    coefficients,
    r2: chunkMX4OLB7V_cjs.roundTo(r2, 6),
    adjR2: chunkMX4OLB7V_cjs.roundTo(adjR2, 6),
    fStatistic: chunkMX4OLB7V_cjs.roundTo(F, 4),
    fDf: [dfModel, dfRes],
    fPValue: chunkMX4OLB7V_cjs.roundTo(fPValue, 4),
    aic: chunkMX4OLB7V_cjs.roundTo(aic, 2),
    bic: chunkMX4OLB7V_cjs.roundTo(bic, 2),
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
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(Array.from({ length: n }, (_, i) => [1, x[i] ?? 0]));
  return fitOLS(X, y, ["(Intercept)", "x"], ciLevel);
}
function multipleRegression(y, predictors, ciLevel = 0.95) {
  if (predictors.length === 0) throw new Error("multipleRegression: need at least 1 predictor");
  const n = y.length;
  for (const p of predictors) {
    if (p.values.length !== n) throw new Error(`multipleRegression: predictor '${p.name}' length mismatch`);
  }
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((p) => p.values[i] ?? 0)])
  );
  const names = ["(Intercept)", ...predictors.map((p) => p.name)];
  return fitOLS(X, y, names, ciLevel);
}
function polynomialRegression(x, y, degree, ciLevel = 0.95) {
  if (degree < 1) throw new Error("polynomialRegression: degree must be \u2265 1");
  if (x.length !== y.length) throw new Error("polynomialRegression: arrays must match length");
  const n = x.length;
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
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
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((pr) => pr.values[i] ?? 0)])
  );
  const names = ["(Intercept)", ...predictors.map((pr) => pr.name)];
  let beta = new Array(p).fill(0);
  for (let iter = 0; iter < maxIter; iter++) {
    const eta2 = Array.from({ length: n }, (_, i) => {
      let v = 0;
      for (let j = 0; j < p; j++) v += X.get(i, j) * (beta[j] ?? 0);
      return v;
    });
    const mu2 = eta2.map((e) => 1 / (1 + Math.exp(-e)));
    const w2 = mu2.map((m) => Math.max(1e-10, m * (1 - m)));
    const Xw2 = chunkFSSEZIKV_cjs.Matrix.fromArray(
      Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_2, j) => X.get(i, j) * Math.sqrt(w2[i])))
    );
    const yAdj = Array.from({ length: n }, (_, i) => Math.sqrt(w2[i]) * ((y[i] ?? 0) - (mu2[i] ?? 0)));
    try {
      const Xwt = Xw2.transpose();
      const XwtXw = Xwt.multiply(Xw2);
      const XwtY = Xwt.multiply(chunkFSSEZIKV_cjs.Matrix.colVec(yAdj));
      const delta = XwtXw.inverse().multiply(XwtY);
      let maxChange = 0;
      for (let j = 0; j < p; j++) {
        const d = delta.get(j, 0);
        beta[j] = (beta[j] ?? 0) + d;
        maxChange = Math.max(maxChange, Math.abs(d));
      }
      if (maxChange < tol) break;
    } catch {
      break;
    }
  }
  const eta = Array.from({ length: n }, (_, i) => {
    let v = 0;
    for (let j = 0; j < p; j++) v += X.get(i, j) * (beta[j] ?? 0);
    return v;
  });
  const mu = eta.map((e) => 1 / (1 + Math.exp(-e)));
  const w = mu.map((m) => Math.max(1e-10, m * (1 - m)));
  const Xw = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => Array.from({ length: p }, (_2, j) => X.get(i, j) * Math.sqrt(w[i])))
  );
  let cov2;
  try {
    cov2 = Xw.transpose().multiply(Xw).inverse();
  } catch {
    cov2 = chunkFSSEZIKV_cjs.Matrix.identity(p);
  }
  const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, n - p);
  const coefficients = beta.map((b, i) => {
    const se2 = Math.sqrt(Math.max(0, cov2.get(i, i)));
    const z = se2 === 0 ? 0 : b / se2;
    const pVal = 2 * (1 - normCDFLocal(Math.abs(z)));
    const ci = [b - tCrit * se2, b + tCrit * se2];
    return {
      name: names[i] ?? `\u03B2${i}`,
      estimate: chunkMX4OLB7V_cjs.roundTo(b, 6),
      se: chunkMX4OLB7V_cjs.roundTo(se2, 6),
      tValue: chunkMX4OLB7V_cjs.roundTo(z, 4),
      pValue: chunkMX4OLB7V_cjs.roundTo(pVal, 4),
      ci
    };
  });
  const logLik = mu.reduce((s, m, i) => {
    const yi = y[i] ?? 0;
    return s + yi * Math.log(Math.max(1e-15, m)) + (1 - yi) * Math.log(Math.max(1e-15, 1 - m));
  }, 0);
  const pMeanRaw = chunkMX4OLB7V_cjs.mean([...y]);
  const pMean = Math.min(1 - 1e-12, Math.max(1e-12, pMeanRaw));
  const nullLogLik = n * (pMean * Math.log(Math.max(1e-15, pMean)) + (1 - pMean) * Math.log(Math.max(1e-15, 1 - pMean)));
  const r2 = Math.abs(nullLogLik) < 1e-12 ? NaN : 1 - logLik / nullLogLik;
  const aic = -2 * logLik + 2 * p;
  const bic = -2 * logLik + Math.log(n) * p;
  const residuals = y.map((v, i) => (v ?? 0) - (mu[i] ?? 0));
  return {
    coefficients,
    r2: chunkMX4OLB7V_cjs.roundTo(r2, 6),
    adjR2: chunkMX4OLB7V_cjs.roundTo(r2, 6),
    // McFadden's for logistic
    fStatistic: NaN,
    fDf: [p - 1, n - p],
    fPValue: NaN,
    aic: chunkMX4OLB7V_cjs.roundTo(aic, 2),
    bic: chunkMX4OLB7V_cjs.roundTo(bic, 2),
    residuals,
    fitted: mu,
    n,
    formatted: `McFadden R\xB2 = ${chunkMX4OLB7V_cjs.roundTo(r2, 3)}, AIC = ${chunkMX4OLB7V_cjs.roundTo(aic, 1)}`
  };
}
function normCDFLocal(z) {
  const x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const erf = 1 - poly * Math.exp(-x * x);
  return 0.5 * (1 + (z >= 0 ? erf : -erf));
}
function regressionDiagnostics(result, predictors) {
  const n = result.n;
  const p = result.coefficients.length;
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predictors.map((pr) => pr.values[i] ?? 0)])
  );
  const Xt = X.transpose();
  let XtXInv;
  try {
    XtXInv = Xt.multiply(X).inverse();
  } catch {
    XtXInv = chunkFSSEZIKV_cjs.Matrix.identity(p);
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
    let ss = 0;
    for (let i = 0; i < n; i++) {
      const diff = data[i][j] - colMeans[j];
      ss += diff * diff;
    }
    return Math.sqrt(ss / (n - 1));
  });
}

// src/stats/pca.ts
function runPCA(data, nComponents, scale = true) {
  if (data.length < 2) throw new Error("runPCA: need at least 2 observations");
  const n = data.length;
  const k = data[0].length;
  if (k < 2) throw new Error("runPCA: need at least 2 variables");
  const pp = preprocessData(data, { method: scale ? "standardize" : "center" });
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(pp.data);
  const Xs = X.scale(1 / Math.sqrt(n - 1));
  const { U: _U, S, V } = Xs.svd();
  const nc = nComponents ?? Math.min(n - 1, k);
  const eigenvalues = S.slice(0, nc).map((s) => s * s);
  const totalVar = S.reduce((sum, s) => sum + s * s, 0);
  const loadings = Array.from(
    { length: k },
    (_, varIdx) => Array.from({ length: nc }, (_2, compIdx) => V.get(varIdx, compIdx))
  );
  const Vk = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: k }, (_, i) => Array.from({ length: nc }, (_2, j) => V.get(i, j)))
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
    eigenvalues: eigenvalues.map((e) => chunkMX4OLB7V_cjs.roundTo(e, 6)),
    varianceExplained: varianceExplained.map((v) => chunkMX4OLB7V_cjs.roundTo(v, 6)),
    cumulativeVariance: cumulativeVariance.map((v) => chunkMX4OLB7V_cjs.roundTo(v, 6)),
    nComponents: nc
  };
}
function varimaxRotation(loadings, maxIter = 1e3, tol = 1e-6) {
  const k = loadings.length;
  const m = loadings[0].length;
  let L = loadings.map((row) => [...row]);
  let T = Array.from(
    { length: m },
    (_, i) => Array.from({ length: m }, (_2, j) => i === j ? 1 : 0)
  );
  for (let iter = 0; iter < maxIter; iter++) {
    let delta = 0;
    for (let p = 0; p < m - 1; p++) {
      for (let q = p + 1; q < m; q++) {
        const u = L.map((row) => (row[p] ?? 0) ** 2 - (row[q] ?? 0) ** 2);
        const v = L.map((row) => 2 * (row[p] ?? 0) * (row[q] ?? 0));
        const A = u.reduce((s, ui) => s + ui, 0);
        const B = v.reduce((s, vi) => s + vi, 0);
        const C = u.reduce((s, ui, i) => s + ui ** 2 - (v[i] ?? 0) ** 2, 0);
        const D = u.reduce((s, ui, i) => s + ui * (v[i] ?? 0), 0) * 2;
        const X_ = C - (A ** 2 - B ** 2) / k;
        const Y_ = D - 2 * A * B / k;
        const angle = Math.atan2(Y_, X_) / 4;
        if (Math.abs(angle) < 1e-12) continue;
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);
        delta += Math.abs(angle);
        const newLp = L.map((row) => (row[p] ?? 0) * cos + (row[q] ?? 0) * sin);
        const newLq = L.map((row) => -(row[p] ?? 0) * sin + (row[q] ?? 0) * cos);
        L.forEach((row, i) => {
          row[p] = newLp[i];
          row[q] = newLq[i];
        });
        for (let r = 0; r < m; r++) {
          const tp = (T[r]?.[p] ?? 0) * cos + (T[r]?.[q] ?? 0) * sin;
          const tq = -(T[r]?.[p] ?? 0) * sin + (T[r]?.[q] ?? 0) * cos;
          T[r][p] = tp;
          T[r][q] = tq;
        }
      }
    }
    if (delta < tol) break;
  }
  return { rotatedLoadings: L, rotationMatrix: T };
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
  const Dmat = ZtZ.add(chunkFSSEZIKV_cjs.Matrix.identity(q).scale(1 / psi));
  let DInv;
  let logDetD;
  try {
    DInv = Dmat.inverse();
    logDetD = Dmat.logDet();
  } catch {
    return { negLogLik: Infinity, sigmae2: 0, sigmab2: 0 };
  }
  const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose());
  const VpsiInv = chunkFSSEZIKV_cjs.Matrix.fromArray(
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
  const XtVinvY = XtVinv.multiply(chunkFSSEZIKV_cjs.Matrix.colVec(y));
  const beta = XtVinvXInv.multiply(XtVinvY);
  const Xbeta = X.multiply(beta);
  const e = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0));
  const eM = chunkFSSEZIKV_cjs.Matrix.colVec(e);
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
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [
      1,
      ...predNames.map((name) => (fixedPredictors[name] ?? [])[i] ?? 0)
    ])
  );
  const Z = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from(
      { length: n },
      (_, i) => groupLevels.map((g) => groupId[i] === g ? 1 : 0)
    )
  );
  const objFn = (theta) => remlProfileLogLik(theta[0] ?? 0, y, X, Z).negLogLik;
  const starts = [-4, -2, 0, 2, 4];
  let optResult = chunkMX4OLB7V_cjs.nelderMead(objFn, [starts[0]], { maxIter: 1e3, tol: 1e-8 });
  for (let si = 1; si < starts.length; si++) {
    const cand = chunkMX4OLB7V_cjs.nelderMead(objFn, [starts[si]], { maxIter: 1e3, tol: 1e-8 });
    if (cand.fval < optResult.fval) optResult = cand;
  }
  const finalModel = remlProfileLogLik(optResult.x[0] ?? 0, y, X, Z);
  const sigmab2 = finalModel.sigmab2;
  const sigmae2 = finalModel.sigmae2;
  const scale = sigmab2 / sigmae2;
  const ZtZ = Z.transpose().multiply(Z);
  let VinvScaled;
  if (scale < 1e-10) {
    VinvScaled = chunkFSSEZIKV_cjs.Matrix.identity(n);
  } else {
    const Dmat = ZtZ.add(chunkFSSEZIKV_cjs.Matrix.identity(nGroups).scale(1 / scale));
    let DInv;
    try {
      DInv = Dmat.inverse();
      const ZDinvZt = Z.multiply(DInv).multiply(Z.transpose());
      VinvScaled = chunkFSSEZIKV_cjs.Matrix.fromArray(
        Array.from(
          { length: n },
          (_, i) => Array.from(
            { length: n },
            (_2, j) => (i === j ? 1 : 0) - ZDinvZt.get(i, j)
          )
        )
      );
    } catch {
      VinvScaled = chunkFSSEZIKV_cjs.Matrix.identity(n);
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
    XtVinvXInv = chunkFSSEZIKV_cjs.Matrix.identity(p);
  }
  const XtVinvY = XtVinv.multiply(chunkFSSEZIKV_cjs.Matrix.colVec([...y]));
  const betaM = XtVinvXInv.multiply(XtVinvY);
  const beta = Array.from({ length: p }, (_, i) => betaM.get(i, 0));
  const df = Math.max(1, n - p - nGroups + 1);
  const tCrit = chunkMX4OLB7V_cjs.tDistQuantile(1 - (1 - ciLevel) / 2, df);
  const covBeta = XtVinvXInv.scale(sigmae2);
  const fixedEffectNames = ["(Intercept)", ...predNames];
  const fixedEffects = beta.map((b, i) => {
    const seVal = Math.sqrt(Math.max(0, covBeta.get(i, i)));
    const t = seVal === 0 ? 0 : b / seVal;
    const pVal = chunkMX4OLB7V_cjs.tDistPValue(t, df);
    return {
      name: fixedEffectNames[i] ?? `\u03B2${i}`,
      estimate: chunkMX4OLB7V_cjs.roundTo(b, 6),
      se: chunkMX4OLB7V_cjs.roundTo(seVal, 6),
      tValue: chunkMX4OLB7V_cjs.roundTo(t, 4),
      pValue: chunkMX4OLB7V_cjs.roundTo(pVal, 4),
      ci: [chunkMX4OLB7V_cjs.roundTo(b - tCrit * seVal, 6), chunkMX4OLB7V_cjs.roundTo(b + tCrit * seVal, 6)]
    };
  });
  const icc = sigmab2 / (sigmab2 + sigmae2);
  const remlConst = 0.5 * (n - p) * (1 + Math.log(2 * Math.PI));
  const logLik = -finalModel.negLogLik - remlConst;
  const aic = -2 * logLik + 2 * (p + 2);
  const bic = -2 * logLik + Math.log(n) * (p + 2);
  const formatted = chunkMX4OLB7V_cjs.formatLMM(icc, aic, bic, logLik);
  return {
    fixedEffects,
    varianceComponents: {
      intercept: chunkMX4OLB7V_cjs.roundTo(sigmab2, 6),
      residual: chunkMX4OLB7V_cjs.roundTo(sigmae2, 6)
    },
    icc: chunkMX4OLB7V_cjs.roundTo(icc, 6),
    logLik: chunkMX4OLB7V_cjs.roundTo(logLik, 4),
    aic: chunkMX4OLB7V_cjs.roundTo(aic, 2),
    bic: chunkMX4OLB7V_cjs.roundTo(bic, 2),
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
  const X = chunkFSSEZIKV_cjs.Matrix.fromArray(
    Array.from({ length: n }, (_, i) => [1, ...predNames.map((name) => (fixedPredictors[name] ?? [])[i] ?? 0)])
  );
  const beta = result.fixedEffects.map((fe) => fe.estimate);
  const Xbeta = X.multiply(chunkFSSEZIKV_cjs.Matrix.colVec(beta));
  const residuals = Array.from({ length: n }, (_, i) => (y[i] ?? 0) - Xbeta.get(i, 0));
  const psi = sigmab2 / sigmae2;
  return groupLevels.map((g) => {
    const indices = Array.from({ length: n }, (_, i) => i).filter((i) => groupId[i] === g);
    const sumResid = indices.reduce((s, i) => s + (residuals[i] ?? 0), 0);
    const nj = indices.length;
    const blup = psi / (1 + psi * nj) * sumResid;
    return { group: g, blup: chunkMX4OLB7V_cjs.roundTo(blup, 6) };
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
  const I = chunkFSSEZIKV_cjs.Matrix.identity(d);
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
        const { values, vectors } = chunkFSSEZIKV_cjs.Matrix.fromArray(cov2).eigen();
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
        const { values, vectors } = chunkFSSEZIKV_cjs.Matrix.fromArray(pool).eigen();
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
    const U = new chunkFSSEZIKV_cjs.Matrix(d, d, uFlats[j]);
    const D = chunkFSSEZIKV_cjs.Matrix.fromArray(Array.from(
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
      formatted: `GMM (K = ${k}, ${modelType}): BIC = ${chunkMX4OLB7V_cjs.roundTo(bic, 1)}, AIC = ${chunkMX4OLB7V_cjs.roundTo(aic, 1)}, LL = ${chunkMX4OLB7V_cjs.roundTo(logL, 1)}, AvePP = [${computeAvePP(resp, k).map((v) => chunkMX4OLB7V_cjs.roundTo(v, 2)).join(", ")}]`
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
      formatted: `LCA (K = ${k}): BIC = ${chunkMX4OLB7V_cjs.roundTo(bic, 1)}, AIC = ${chunkMX4OLB7V_cjs.roundTo(aic, 1)}, LL = ${chunkMX4OLB7V_cjs.roundTo(logL, 1)}`
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
      const logBeta = Array.from({ length: T }, () => new Array(k));
      for (let s = 0; s < k; s++) logBeta[T - 1][s] = 0;
      for (let t = T - 2; t >= 0; t--) {
        for (let s = 0; s < k; s++) {
          const combined = new Float64Array(k);
          for (let next = 0; next < k; next++) {
            combined[next] = Math.log(Math.max(tau[s][next], MIN_PROB)) + logB[t + 1][next] + logBeta[t + 1][next];
          }
          logBeta[t][s] = logSumExp(combined);
        }
      }
      const gamma = Array.from({ length: T }, () => new Array(k));
      for (let t = 0; t < T; t++) {
        const logRow = new Float64Array(k);
        for (let s = 0; s < k; s++) logRow[s] = logAlpha[t][s] + logBeta[t][s];
        const den = logSumExp(logRow);
        for (let s = 0; s < k; s++) gamma[t][s] = Math.exp(logRow[s] - den);
      }
      iterGamma[i] = gamma;
      for (let t = 0; t < T - 1; t++) {
        const logXi = new Float64Array(k * k);
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            logXi[j * k + l] = logAlpha[t][j] + Math.log(Math.max(tau[j][l], MIN_PROB)) + logB[t + 1][l] + logBeta[t + 1][l];
          }
        }
        const xiDen = logSumExp(logXi);
        for (let j = 0; j < k; j++) {
          for (let l = 0; l < k; l++) {
            tauNum[j][l] += Math.exp(logXi[j * k + l] - xiDen);
          }
        }
      }
      for (let s = 0; s < k; s++) piAcc[s] += gamma[0][s];
      for (let t = 0; t < T - 1; t++) {
        for (let s = 0; s < k; s++) tauDen[s] += gamma[t][s];
      }
      for (let t = 0; t < T; t++) {
        for (let s = 0; s < k; s++) {
          rhoDen[s] += gamma[t][s];
          for (let d = 0; d < m; d++) {
            if (data[i][t][d] === 1) rhoNum[s][d] += gamma[t][s];
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
      formatted: `LTA (K = ${k}, T = ${T}): BIC = ${chunkMX4OLB7V_cjs.roundTo(bic, 1)}, AIC = ${chunkMX4OLB7V_cjs.roundTo(aic, 1)}, LL = ${chunkMX4OLB7V_cjs.roundTo(logL, 1)}`
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
  const formatted = `DBSCAN (eps = ${chunkMX4OLB7V_cjs.roundTo(eps, 3)}, minPts = ${minPts}): ${nClusters} clusters, ${nNoise} noise points, silhouette = ${isNaN(sil.mean) ? "N/A" : chunkMX4OLB7V_cjs.roundTo(sil.mean, 3)}`;
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
  const formatted = `HAC (${linkage}): ${n} observations, cophenetic r = ${chunkMX4OLB7V_cjs.roundTo(cophCorr, 3)}`;
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
  return chunkFSSEZIKV_cjs.Matrix.fromArray(R);
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
  const pValue = df > 0 ? chunkMX4OLB7V_cjs.chiSqPValue(chiSq, df) : 1;
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
  const diagMat = chunkFSSEZIKV_cjs.Matrix.fromArray(diagArr);
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
    const adjM = chunkFSSEZIKV_cjs.Matrix.fromArray(adjR);
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
    const { values, vectors } = chunkFSSEZIKV_cjs.Matrix.fromArray(scaledR).eigen();
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
    const Sigma = chunkFSSEZIKV_cjs.Matrix.fromArray(sigmaArr);
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
    const { values } = chunkFSSEZIKV_cjs.Matrix.fromArray(scaledR).eigen();
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
  const nmResult = chunkMX4OLB7V_cjs.nelderMead(concentratedML, psi0, {
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
    let ss = 0;
    for (let j = 0; j < k; j++) ss += L[i][j] ** 2;
    sc[i] = Math.sqrt(ss || 1e-15);
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
    const { values: sigma2, vectors: Vmat } = chunkFSSEZIKV_cjs.Matrix.fromArray(BtB).eigen();
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
function criterionOblimin(L, gamma) {
  const p = L.length, k = L[0].length;
  const Gq = Array.from({ length: p }, () => new Array(k).fill(0));
  let f = 0;
  for (let i = 0; i < p; i++) {
    let rowSumSq = 0;
    for (let m = 0; m < k; m++) rowSumSq += L[i][m] ** 2;
    for (let j = 0; j < k; j++) {
      const lij = L[i][j];
      const otherSumSq = rowSumSq - lij * lij;
      Gq[i][j] = lij * (otherSumSq - gamma / p * rowSumSq);
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
  const Amat = chunkFSSEZIKV_cjs.Matrix.fromArray(A);
  let Tarr = Tinit ? Tinit.map((r) => [...r]) : Array.from(
    { length: k },
    (_, i) => Array.from({ length: k }, (_2, j) => i === j ? 1 : 0)
  );
  let Tmat = chunkFSSEZIKV_cjs.Matrix.fromArray(Tarr);
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
        let ss = 0;
        for (let i = 0; i < k; i++) ss += X[i][j] ** 2;
        const invNorm = 1 / Math.sqrt(ss || 1e-15);
        for (let i = 0; i < k; i++) X[i][j] = X[i][j] * invNorm;
      }
      const Ttmat = chunkFSSEZIKV_cjs.Matrix.fromArray(X);
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
  const Lmat = chunkFSSEZIKV_cjs.Matrix.fromArray(L);
  const GqMat = chunkFSSEZIKV_cjs.Matrix.fromArray(Gq);
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
  const Vmat = chunkFSSEZIKV_cjs.Matrix.fromArray(vari);
  const Qmat = chunkFSSEZIKV_cjs.Matrix.fromArray(Q);
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
  const Umat = chunkFSSEZIKV_cjs.Matrix.fromArray(Uarr);
  const rotatedNorm = Vmat.multiply(Umat).toArray();
  const rotated = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => rotatedNorm[i][j] * sqrtH2[i])
  );
  const TvarMat = chunkFSSEZIKV_cjs.Matrix.fromArray(T_varimax);
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
    const gamma = method === "quartimin" ? 0 : 0;
    const criterionFn = (L) => criterionOblimin(L, gamma);
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
    const Lmat = chunkFSSEZIKV_cjs.Matrix.fromArray(Larr);
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
  const pValue = chunkMX4OLB7V_cjs.chiSqPValue(Math.max(chiSq, 0), df);
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
  return chunkFSSEZIKV_cjs.Matrix.fromArray(sigma);
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
    const LMat = chunkFSSEZIKV_cjs.Matrix.fromArray(L);
    const PhiMat = chunkFSSEZIKV_cjs.Matrix.fromArray(Phi);
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
    covMat = chunkFSSEZIKV_cjs.Matrix.fromArray(infoArr).inverse();
  } catch {
    try {
      covMat = chunkFSSEZIKV_cjs.Matrix.fromArray(infoArr).pseudoInverse();
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
  const formatted = `EFA (${extraction}/${rotation}): ${chunkMX4OLB7V_cjs.formatCFAFit(fit)}`;
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
      const pValue = 2 * (1 - chunkMX4OLB7V_cjs.normalCDF(Math.abs(z)));
      const sigmaII = Math.max(impliedSigma.get(r, r), 1e-12);
      const stdAll = est * Math.sqrt(Phi[c][c]) / Math.sqrt(sigmaII);
      return { estimate: chunkMX4OLB7V_cjs.roundTo(est, 4), se: chunkMX4OLB7V_cjs.roundTo(se2, 4), z: chunkMX4OLB7V_cjs.roundTo(z, 3), pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4), stdAll: chunkMX4OLB7V_cjs.roundTo(stdAll, 4) };
    })
  );
  const paramUniqueness = Array.from(Theta).map((est, i) => {
    const se2 = Math.max(thetaSE[i], 1e-6);
    const z = est / se2;
    const pValue = 2 * (1 - chunkMX4OLB7V_cjs.normalCDF(Math.abs(z)));
    const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12);
    const stdAll = est / sigmaII;
    return { estimate: chunkMX4OLB7V_cjs.roundTo(est, 4), se: chunkMX4OLB7V_cjs.roundTo(se2, 4), z: chunkMX4OLB7V_cjs.roundTo(z, 3), pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4), stdAll: chunkMX4OLB7V_cjs.roundTo(stdAll, 4) };
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
      const pValue = 2 * (1 - chunkMX4OLB7V_cjs.normalCDF(Math.abs(z)));
      return { estimate: chunkMX4OLB7V_cjs.roundTo(est, 4), se: chunkMX4OLB7V_cjs.roundTo(se2, 4), z: chunkMX4OLB7V_cjs.roundTo(z, 3), pValue: chunkMX4OLB7V_cjs.roundTo(pValue, 4), stdAll: chunkMX4OLB7V_cjs.roundTo(est, 4) };
    })
  );
  const communalities = Array.from(Theta).map((t) => chunkMX4OLB7V_cjs.roundTo(1 - t, 4));
  const uniquenessArr = Array.from(Theta).map((t) => chunkMX4OLB7V_cjs.roundTo(t, 4));
  const stdLoadings = Array.from(
    { length: d },
    (_, i) => Array.from({ length: k }, (_2, j) => {
      const sigmaII = Math.max(impliedSigma.get(i, i), 1e-12);
      return chunkMX4OLB7V_cjs.roundTo(L[i][j] * Math.sqrt(Phi[j][j]) / Math.sqrt(sigmaII), 4);
    })
  );
  const eigenvalues = S.eigen().values.map((v) => chunkMX4OLB7V_cjs.roundTo(v, 4));
  const variableNames = options?.variableNames ?? Array.from({ length: d }, (_, i) => `V${i + 1}`);
  const factorNames = options?.factorNames ?? factors;
  const formatted = `CFA: ${chunkMX4OLB7V_cjs.formatCFAFit(fit)}`;
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

exports.analyze = analyze;
exports.chiSquareTest = chiSquareTest;
exports.ciMean = ciMean;
exports.cohensD = cohensD;
exports.cohensDCI = cohensDCI;
exports.cohensDPaired = cohensDPaired;
exports.computeBLUPs = computeBLUPs;
exports.contingencyTable = contingencyTable;
exports.correlationMatrix = correlationMatrix;
exports.cutTree = cutTree;
exports.cutTreeHeight = cutTreeHeight;
exports.describe = describe;
exports.detectFieldType = detectFieldType;
exports.dunnTest = dunnTest;
exports.etaSquared = etaSquared;
exports.etaSquaredKW = etaSquaredKW;
exports.euclideanDistMatrix = euclideanDistMatrix;
exports.findBestGMM = findBestGMM;
exports.fisherExactTest = fisherExactTest;
exports.fitGMM = fitGMM;
exports.fitGMMRange = fitGMMRange;
exports.fitKMeansRange = fitKMeansRange;
exports.fitLCA = fitLCA;
exports.fitLTA = fitLTA;
exports.frequencyTable = frequencyTable;
exports.friedmanTest = friedmanTest;
exports.gamesHowell = gamesHowell;
exports.goodnessOfFit = goodnessOfFit;
exports.hedgesG = hedgesG;
exports.inverseTransform = inverseTransform;
exports.kDistancePlot = kDistancePlot;
exports.kendallTau = kendallTau;
exports.kruskalWallis = kruskalWallis;
exports.kurtosis = kurtosis;
exports.linearRegression = linearRegression;
exports.logisticRegression = logisticRegression;
exports.mannWhitneyU = mannWhitneyU;
exports.multipleRegression = multipleRegression;
exports.omegaSquared = omegaSquared;
exports.oneWayANOVA = oneWayANOVA;
exports.partialCorrelation = partialCorrelation;
exports.pearsonCorrelation = pearsonCorrelation;
exports.phiCoefficient = phiCoefficient;
exports.polynomialRegression = polynomialRegression;
exports.predictGMM = predictGMM;
exports.predictKMeans = predictKMeans;
exports.preprocessData = preprocessData;
exports.rankBiserial = rankBiserial;
exports.rankBiserialWilcoxon = rankBiserialWilcoxon;
exports.regressionDiagnostics = regressionDiagnostics;
exports.runCFA = runCFA;
exports.runDBSCAN = runDBSCAN;
exports.runEFA = runEFA;
exports.runFADiagnostics = runFADiagnostics;
exports.runHierarchical = runHierarchical;
exports.runKMeans = runKMeans;
exports.runLMM = runLMM;
exports.runPCA = runPCA;
exports.screeData = screeData;
exports.shapiroWilk = shapiroWilk;
exports.silhouetteScores = silhouetteScores;
exports.skewness = skewness;
exports.spearmanCorrelation = spearmanCorrelation;
exports.tTestIndependent = tTestIndependent;
exports.tTestPaired = tTestPaired;
exports.trimmedMean = trimmedMean;
exports.tukeyHSD = tukeyHSD;
exports.varimaxRotation = varimaxRotation;
exports.wilcoxonSignedRank = wilcoxonSignedRank;
//# sourceMappingURL=chunk-RPPDUX5D.cjs.map
//# sourceMappingURL=chunk-RPPDUX5D.cjs.map