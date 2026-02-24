import { quantile, sortAsc, normalQuantile, mean, normalCDF } from './chunk-IRX4LIZX.js';

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
  });
  g.append("text").attr("transform", "rotate(-90)").attr("x", -height / 2).attr("y", -52).attr("text-anchor", "middle").attr("font-family", theme.fontFamily).attr("font-size", theme.fontSize).attr("fill", theme.text).text(label);
}
function renderGridLines(g, scale, width, theme = DEFAULT_THEME) {
  const ticks = scale.ticks(6);
  g.selectAll(".grid-line").data(ticks).join("line").attr("class", "grid-line").attr("x1", 0).attr("x2", width).attr("y1", scale).attr("y2", scale).attr("stroke", theme.gridLine).attr("stroke-width", 1);
}

// src/viz/plots/violin-box.ts
function renderViolinBox(container, data, config = {}) {
  import('d3').then((d3) => renderViolinBoxD3(d3, container, data, config));
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
  import('d3').then((d3) => renderScatterD3(d3, container, data, config));
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
  import('d3').then((d3) => renderHistogramD3(d3, container, data, config));
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
  import('d3').then((d3) => renderBarD3(d3, container, data, config));
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
  import('d3').then((d3) => renderCorrelogramD3(d3, container, data, config));
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
  const margin = { top: 60, left: 80 };
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
  import('d3').then((d3) => renderCoefD3(d3, container, coefficients, config));
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
  import('d3').then((d3) => renderQQD3(d3, container, values, config));
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
  import('d3').then((d3) => renderResidualD3(d3, container, result, leverage, config));
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
  import('d3').then((d3) => renderRaincloudD3(d3, container, data, config));
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
  import('d3').then((d3) => renderMixedD3(d3, container, result, blups, config));
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
  import('d3').then((d3) => {
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
  const margin = { top: 60, left: 100 };
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
  import('d3').then((d3) => renderDistD3(d3, container, params, config));
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
  import('d3').then((d3) => renderDensityD3(d3, container, data, config));
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
  import('d3').then((d3) => renderBoxplotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderLollipopD3(d3, container, data, config));
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
  import('d3').then((d3) => renderDotPlotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderGroupedBarD3(d3, container, data, config));
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
  import('d3').then((d3) => renderLineChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderBubbleChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderParetoD3(d3, container, data, config));
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
  import('d3').then((d3) => renderFunnelD3(d3, container, data, config));
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
  import('d3').then((d3) => renderPieChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderAreaChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderForestPlotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderROCCurveD3(d3, container, data, config));
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
  import('d3').then((d3) => renderStripPlotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderSwarmPlotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderMosaicPlotD3(d3, container, data, config));
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
  import('d3').then((d3) => renderPairPlotD3(d3, container, data, config));
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
        drawHistogramCell(d3, cellG, data.data[row], scales[col], cellW, cellH, cellPad, getColor(row, theme));
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
        const r = pearsonR(xVals.slice(0, n_obs), yVals.slice(0, n_obs));
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
function pearsonR(xs, ys) {
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
  import('d3').then((d3) => renderRadarChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderParallelCoordsD3(d3, container, data, config));
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
  import('d3').then((d3) => renderTreemapD3(d3, container, data, config));
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
  import('d3').then((d3) => renderWaffleChartD3(d3, container, data, config));
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
  import('d3').then((d3) => renderSparklineD3(d3, container, data, config));
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
  import('d3').then((d3) => renderSunburstD3(d3, container, data, config));
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
  import('d3').then((d3) => renderMarimekkoD3(d3, container, data, config));
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
  import('d3').then((d3) => renderChordDiagramD3(d3, container, data, config));
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
  import('d3').then((d3) => renderArcDiagramD3(d3, container, data, config));
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
  import('d3').then((d3) => renderAlluvialD3(d3, container, data, config));
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
  import('d3').then((d3) => renderEdgeBundlingD3(d3, container, data, config));
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

export { CARM_PALETTE, DARK_THEME, DEFAULT_THEME, OKABE_ITO, addCaption, addNLabel, addRegressionEquation, addStatBadge, addSubtitle, applyTheme, exportPNG, exportSVG, formatTooltipRow, getColor, hideTooltip, renderAlluvialPlot, renderArcDiagram, renderAreaChart, renderBarStats, renderBoxplot, renderBrackets, renderBubbleChart, renderChordDiagram, renderCoefPlot, renderCorrelogram, renderDensity, renderDistribution, renderDotPlot, renderEdgeBundling, renderForestPlot, renderFunnel, renderGridLines, renderGroupedBar, renderHistogram, renderLineChart, renderLollipop, renderMarimekko, renderMixedPlot, renderMosaicPlot, renderPCAPlot, renderPairPlot, renderParallelCoords, renderPareto, renderPieChart, renderQQPlot, renderROCCurve, renderRadarChart, renderRaincloud, renderResidualPanel, renderScatterStats, renderSparkline, renderStripPlot, renderSunburst, renderSwarmPlot, renderTreemap, renderViolinBox, renderWaffleChart, renderXAxis, renderYAxis, showTooltip, themeColorScale, totalBracketHeight };
//# sourceMappingURL=chunk-MPRPR3Z3.js.map
//# sourceMappingURL=chunk-MPRPR3Z3.js.map