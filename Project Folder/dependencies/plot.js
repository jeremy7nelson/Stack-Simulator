'use strict';

import { $, state, onDataUpdated, IMG_W, IMG_H, MAX16, STACK_LEAD_BASELINE_SEC, fmtConc } from './main.js';
import {
  getRawCurves, getNsaParams, getDriftParams, A_D,
  isNSAOn, isDriftOn, isPixelNoiseOn, getPixelNoiseParams, pixelNoiseValue
} from './noise.js';
import { regionBrightness16, pausePlayback } from './render.js';

function cssVar(name, fallback) {
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return v || fallback;
}

const COLORS = {
  rule: cssVar('--line', '#ddd6c8'),
  ruleSoft: cssVar('--line-soft', '#e8e2d6'),
  muted: cssVar('--muted', '#746f66'),
  ink: cssVar('--ink', '#1b1a17'),
  accent: cssVar('--accent', '#0f6b66'),
  accentSoft: cssVar('--accent-soft', '#e2efed'),
  warn: cssVar('--warn', '#b4541f'),
  panel: cssVar('--panel', '#fffdf8'),
  nsa: '#7a3e9d',
  drift: '#b0343f',
  jit: '#78838f'
};

function hexToRgb(hex) {
  const h = hex.replace('#', '');
  return [parseInt(h.substring(0, 2), 16), parseInt(h.substring(2, 4), 16), parseInt(h.substring(4, 6), 16)];
}

function hexToRgba(hex, alpha) {
  const [r, g, b] = hexToRgb(hex);
  return `rgba(${r},${g},${b},${alpha})`;
}

const RAMP_S = [[0, ...hexToRgb(COLORS.panel)], [0.55, 150, 196, 192], [1, ...hexToRgb(COLORS.accent)]];
const RAMP_W = [[0, ...hexToRgb(COLORS.panel)], [0.55, 214, 168, 110], [1, ...hexToRgb(COLORS.warn)]];

function parseConcs(str) {
  return (str || '').split(/[\s,;]+/).map(Number).filter(v => Number.isFinite(v) && v > 0);
}

function getTiming() {
  const tBase = +($('tBase')?.value ?? 0);
  const tAssoc = +($('tAssoc')?.value ?? 0);
  const tDissoc = +($('tDissoc')?.value ?? 0);
  const concs = parseConcs($('concSeries')?.value);
  const cyc = STACK_LEAD_BASELINE_SEC + tAssoc + tDissoc;
  return { tBase, tAssoc, tDissoc, concs, cyc };
}

function decodeContinuousT(t, timing) {
  const { tBase, cyc, concs } = timing;
  if (t < tBase || !concs.length) return { concIdx: -1, timeIdx: 0 };
  let concIdx = Math.floor((t - tBase) / cyc);
  if (concIdx >= concs.length) concIdx = concs.length - 1;
  const timeIdx = Math.round((t - tBase) - concIdx * cyc);
  return { concIdx, timeIdx };
}

function plotAxes(ctx, w, h, pad, yMin, yMax, yLabel, tMax) {
  ctx.clearRect(0, 0, w, h);
  ctx.font = '10px "IBM Plex Mono", monospace';
  ctx.strokeStyle = COLORS.rule;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, h - pad.b);
  ctx.lineTo(w - pad.r, h - pad.b);
  ctx.stroke();

  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let g = 0; g <= 4; g++) {
    const v = yMin + (yMax - yMin) * g / 4;
    const y = h - pad.b - (h - pad.t - pad.b) * g / 4;
    ctx.fillStyle = COLORS.muted;
    ctx.fillText(v.toFixed(Math.abs(yMax) < 1 ? 3 : 1), pad.l - 6, y);
    if (g > 0) {
      ctx.strokeStyle = COLORS.ruleSoft;
      ctx.beginPath();
      ctx.moveTo(pad.l, y);
      ctx.lineTo(w - pad.r, y);
      ctx.stroke();
    }
  }

  ctx.fillStyle = COLORS.ink;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';
  ctx.fillText(yLabel, pad.l, pad.t - 6);

  ctx.fillStyle = COLORS.muted;
  ctx.textAlign = 'left';
  ctx.fillText('0', pad.l, h - pad.b + 16);
  ctx.textAlign = 'right';
  ctx.fillText(Math.round(tMax) + ' s', w - pad.r, h - pad.b + 16);
}

function injectionBands(ctx, w, h, pad, timing) {
  const { tBase, tAssoc, cyc, concs } = timing;
  const tMax = tBase + concs.length * cyc;
  const px = t => pad.l + (w - pad.l - pad.r) * (t / (tMax || 1));

  ctx.fillStyle = COLORS.accentSoft;
  concs.forEach((c, q) => {
    const x0 = px(tBase + q * cyc + STACK_LEAD_BASELINE_SEC);
    const x1 = px(tBase + q * cyc + STACK_LEAD_BASELINE_SEC + tAssoc);
    ctx.fillRect(x0, pad.t, x1 - x0, h - pad.t - pad.b);
  });

  ctx.fillStyle = COLORS.muted;
  ctx.font = '9px "IBM Plex Mono", monospace';
  ctx.textAlign = 'left';
  concs.forEach((c, q) => {
    const x0 = px(tBase + q * cyc + STACK_LEAD_BASELINE_SEC);
    ctx.fillText(fmtConc(c), x0 + 2, pad.t + 9);
  });
}

function playheadLine(ctx, w, h, pad, tMax, tNow) {
  if (tNow == null) return;
  const x = pad.l + (w - pad.l - pad.r) * (tNow / (tMax || 1));
  ctx.strokeStyle = COLORS.warn;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(x, pad.t);
  ctx.lineTo(x, h - pad.b);
  ctx.stroke();
}

function drawBandChart(canvasId, grid, loArr, hiArr, colorHex, fillRgba, yMin, yMax, yLabel, timing, tNow) {
  const canvas = $(canvasId);
  if (!canvas || !grid.length) return;
  const w = canvas.width, h = canvas.height, pad = { l: 46, r: 14, t: 20, b: 26 };
  const ctx = canvas.getContext('2d');
  const tMax = grid[grid.length - 1];

  plotAxes(ctx, w, h, pad, yMin, yMax, yLabel, tMax);
  injectionBands(ctx, w, h, pad, timing);

  const X = t => pad.l + (w - pad.l - pad.r) * (t / (tMax || 1));
  const Y = v => h - pad.b - (h - pad.t - pad.b) * ((v - yMin) / ((yMax - yMin) || 1));

  ctx.beginPath();
  for (let i = 0; i < grid.length; i++) {
    const x = X(grid[i]), y = Y(hiArr[i]);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  for (let i = grid.length - 1; i >= 0; i--) ctx.lineTo(X(grid[i]), Y(loArr[i]));
  ctx.closePath();
  ctx.fillStyle = fillRgba;
  ctx.fill();

  ctx.strokeStyle = colorHex;
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  for (let i = 0; i < grid.length; i++) {
    const x = X(grid[i]), y = Y(hiArr[i]);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();

  playheadLine(ctx, w, h, pad, tMax, tNow);
}

function renderYnsChart(timing, tNow) {
  const { grid, Yns } = getRawCurves();
  if (!grid || !grid.length) return;
  const { a_n, b_n } = getNsaParams();
  const lower = Yns.map(y => a_n * y);
  const upper = Yns.map(y => (a_n + b_n) * y);
  const yMax = Math.max(...upper, 1e-6);
  drawBandChart('plot-yns', grid, lower, upper, COLORS.accent, hexToRgba(COLORS.accent, .14),
    0, yMax, 'Yns(t) — background a_n·Yns vs cell-margin (a_n+b_n)·Yns', timing, tNow);
}

function renderDriftChart(timing, tNow) {
  const { grid, driftCommon } = getRawCurves();
  if (!grid || !grid.length) return;
  const { b_d } = getDriftParams();
  const lower = Array.from(driftCommon, v => A_D * v);
  const upper = Array.from(driftCommon, v => (A_D + b_d) * v);
  const yMin = Math.min(0, ...lower, ...upper);
  const yMax = Math.max(...lower, ...upper, 1e-6);
  drawBandChart('plot-drift', grid, lower, upper, COLORS.drift, hexToRgba(COLORS.drift, .14),
    yMin, yMax, 'driftCommon(t) — background a_d·drift vs cell-margin (a_d+b_d)·drift', timing, tNow);
}

function rampColor(stops, u) {
  u = u < 0 ? 0 : (u > 1 ? 1 : u);
  for (let i = 1; i < stops.length; i++) {
    if (u <= stops[i][0]) {
      const a = stops[i - 1], b = stops[i];
      const f = (u - a[0]) / ((b[0] - a[0]) || 1);
      return [a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, a[3] + (b[3] - a[3]) * f];
    }
  }
  const L = stops[stops.length - 1];
  return [L[1], L[2], L[3]];
}

function rampCSS(stops) {
  const parts = [];
  for (let q = 0; q <= 10; q++) {
    const c = rampColor(stops, q / 10);
    parts.push(`rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0}) ${q * 10}%`);
  }
  return 'linear-gradient(to right,' + parts.join(',') + ')';
}

function paintField(canvasId, values, stops) {
  const canvas = $(canvasId);
  if (!canvas || !values) return null;
  canvas.width = IMG_W;
  canvas.height = IMG_H;
  let lo = Infinity, hi = -Infinity;
  for (let k = 0; k < values.length; k++) {
    if (values[k] < lo) lo = values[k];
    if (values[k] > hi) hi = values[k];
  }
  if (!isFinite(lo)) { lo = 0; hi = 1; }
  const span = (hi - lo) || 1;
  const ctx = canvas.getContext('2d');
  const img = ctx.createImageData(IMG_W, IMG_H);
  for (let k = 0, q = 0; k < values.length; k++, q += 4) {
    const c = rampColor(stops, (values[k] - lo) / span);
    img.data[q] = c[0]; img.data[q + 1] = c[1]; img.data[q + 2] = c[2]; img.data[q + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  return { lo, hi };
}

function setText(id, text) {
  const el = $(id);
  if (el) el.textContent = text;
}

function renderSpatialFields(field) {
  if (!field) return;
  const { edgeWeight, capacityOf, w } = field;
  const s = new Float32Array(edgeWeight.length);
  for (let k = 0; k < s.length; k++) s[k] = edgeWeight[k] * capacityOf[k];

  const sBar = $('s-bar'), wBar = $('w-bar');
  if (sBar) sBar.style.background = rampCSS(RAMP_S);
  if (wBar) wBar.style.background = rampCSS(RAMP_W);

  const sRange = paintField('plot-s', s, RAMP_S);
  const wRange = paintField('plot-w', w, RAMP_W);
  if (sRange) { setText('s-lo', sRange.lo.toFixed(2)); setText('s-hi', sRange.hi.toFixed(2)); }
  if (wRange) { setText('w-lo', wRange.lo.toFixed(2)); setText('w-hi', wRange.hi.toFixed(2)); }
}

function findProbePixel(field) {
  const { regionOf, edgeWeight, capacityOf } = field;
  let best = -1, bestVal = -1;
  for (let k = 0; k < regionOf.length; k++) {
    if (regionOf[k] < 0) continue;
    const v = edgeWeight[k] * capacityOf[k];
    if (v > bestVal) { bestVal = v; best = k; }
  }
  return best;
}

function computeDecompSeries(field, k, timing) {
  const { grid, Yns, driftCommon } = getRawCurves();
  const { nsaRate, gDrift, regionOf, edgeWeight, capacityOf } = field;
  const rIdx = regionOf[k];
  const region = rIdx >= 0 && state.parsed ? state.parsed.regions[rIdx] : null;
  const nsaOn = isNSAOn(), driftOn = isDriftOn(), pxOn = isPixelNoiseOn();
  let pxSigma = 0, pxSeed = 0;
  if (pxOn) ({ sigma: pxSigma, seed: pxSeed } = getPixelNoiseParams());
  const i = Math.floor(k / IMG_W), j = k % IMG_W;
  const weight = edgeWeight[k] * capacityOf[k];

  const spec = new Float64Array(grid.length);
  const nsa = new Float64Array(grid.length);
  const drift = new Float64Array(grid.length);
  const jit = new Float64Array(grid.length);

  for (let idx = 0; idx < grid.length; idx++) {
    const t = grid[idx];
    const { concIdx, timeIdx } = decodeContinuousT(t, timing);
    if (region && concIdx >= 0) spec[idx] = regionBrightness16(region, concIdx, timeIdx) * weight;
    if (nsaOn) nsa[idx] = nsaRate[k] * Yns[idx] * MAX16;
    if (driftOn) drift[idx] = gDrift[k] * driftCommon[idx] * MAX16;
    if (pxOn && concIdx >= 0) jit[idx] = pixelNoiseValue(i, j, timeIdx, pxSigma, pxSeed) * MAX16;
  }
  return { grid, spec, nsa, drift, jit, weight };
}

function clearReadout() {
  const tbody = document.querySelector('#decomp-readout tbody');
  if (tbody) tbody.innerHTML = '';
}

function updateReadout(field, k, series, tNow) {
  const tbody = document.querySelector('#decomp-readout tbody');
  if (!tbody) return;
  if (tNow == null) { tbody.innerHTML = ''; return; }

  const { grid } = getRawCurves();
  let idx = 0, bestDiff = Infinity;
  for (let i = 0; i < grid.length; i++) {
    const diff = Math.abs(grid[i] - tNow);
    if (diff < bestDiff) { bestDiff = diff; idx = i; }
  }

  const pct = v => (v / MAX16 * 100).toFixed(2) + '%';
  const { sigma } = getPixelNoiseParams();
  const rows = [
    ['drift', `gDrift = ${field.gDrift[k].toFixed(3)}`, pct(series.drift[idx])],
    ['NSA', `nsaRate = ${field.nsaRate[k].toFixed(3)}`, pct(series.nsa[idx])],
    ['specific', `edge×cap = ${series.weight.toFixed(3)}`, pct(series.spec[idx])],
    ['jitter', `σ = ${sigma.toFixed(3)}`, pct(series.jit[idx])]
  ];
  const total = series.spec[idx] + series.nsa[idx] + series.drift[idx] + series.jit[idx];

  tbody.innerHTML = rows.map(r =>
    `<tr><td>${r[0]}</td><td class="num">${r[1]}</td><td class="num">${r[2]}</td></tr>`
  ).join('') + `<tr><td><b>total</b></td><td></td><td class="num"><b>${pct(total)}</b></td></tr>`;
}

function renderDecomposition(field, timing, tNow) {
  const canvas = $('plot-decomp');
  if (!canvas) return;
  if (!field) { clearReadout(); return; }

  const kProbe = findProbePixel(field);
  if (kProbe < 0) { clearReadout(); return; }

  const series = computeDecompSeries(field, kProbe, timing);
  const { grid, spec, nsa, drift, jit } = series;

  const total = new Float64Array(grid.length);
  for (let idx = 0; idx < grid.length; idx++) total[idx] = spec[idx] + nsa[idx] + drift[idx] + jit[idx];
  let lo = 0, hi = -Infinity;
  for (let idx = 0; idx < grid.length; idx++) if (total[idx] > hi) hi = total[idx];
  if (!isFinite(hi) || hi <= 0) hi = MAX16 * 0.1;

  const w = canvas.width, h = canvas.height, pad = { l: 46, r: 14, t: 20, b: 26 };
  const ctx = canvas.getContext('2d');
  const tMax = grid[grid.length - 1];

  plotAxes(ctx, w, h, pad, lo / MAX16 * 100, hi / MAX16 * 100, 'stack composition — % of full range', tMax);
  injectionBands(ctx, w, h, pad, timing);

  const X = t => pad.l + (w - pad.l - pad.r) * (t / (tMax || 1));
  const Y = v => h - pad.b - (h - pad.t - pad.b) * ((v - lo) / ((hi - lo) || 1));

  const layers = [
    { arr: drift, color: COLORS.drift },
    { arr: nsa, color: COLORS.nsa },
    { arr: spec, color: COLORS.accent },
    { arr: jit, color: COLORS.jit }
  ];
  const lower = new Float64Array(grid.length);
  const upper = new Float64Array(grid.length);
  for (const layer of layers) {
    for (let idx = 0; idx < grid.length; idx++) upper[idx] = lower[idx] + layer.arr[idx];

    ctx.beginPath();
    for (let idx = 0; idx < grid.length; idx++) {
      const x = X(grid[idx]), y = Y(upper[idx]);
      if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    for (let idx = grid.length - 1; idx >= 0; idx--) ctx.lineTo(X(grid[idx]), Y(lower[idx]));
    ctx.closePath();
    ctx.fillStyle = hexToRgba(layer.color, .28);
    ctx.fill();

    ctx.strokeStyle = layer.color;
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let idx = 0; idx < grid.length; idx++) {
      const x = X(grid[idx]), y = Y(upper[idx]);
      if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    for (let idx = 0; idx < grid.length; idx++) lower[idx] = upper[idx];
  }

  playheadLine(ctx, w, h, pad, tMax, tNow);
  updateReadout(field, kProbe, series, tNow);
}

let lastFieldRef = null;
let plotsT = 0;

function renderAllPlots() {
  const timing = getTiming();
  const tNow = plotsT;

  renderYnsChart(timing, tNow);
  renderDriftChart(timing, tNow);

  const field = state.capacityField;
  if (field !== lastFieldRef) {
    renderSpatialFields(field);
    lastFieldRef = field;
  }
  renderDecomposition(field, timing, tNow);
}

const PLOTS_STEP = 2;
let plotsPlaying = false;
let plotsRaf = null;

function plotsTick() {
  if (!plotsPlaying) return;
  renderAllPlots();
  const { grid } = getRawCurves();
  const tMax = grid.length ? grid[grid.length - 1] : 0;
  plotsT += PLOTS_STEP;
  if (plotsT > tMax) plotsT = 0;
  plotsRaf = requestAnimationFrame(plotsTick);
}

function startPlotsTimer() {
  plotsT = 0;
  plotsPlaying = true;
  if (plotsRaf) cancelAnimationFrame(plotsRaf);
  plotsRaf = requestAnimationFrame(plotsTick);
}

function stopPlotsTimer() {
  plotsPlaying = false;
  if (plotsRaf) cancelAnimationFrame(plotsRaf);
  plotsRaf = null;
}

let activeTab = 'tab-simulator';

function showTab(tabId) {
  if (tabId === activeTab) return;

  if (activeTab === 'tab-simulator' && tabId === 'tab-plots') pausePlayback();
  if (activeTab === 'tab-plots' && tabId === 'tab-simulator') stopPlotsTimer();

  document.querySelectorAll('.tab-panel').forEach(el => {
    el.style.display = el.id === tabId ? '' : 'none';
  });
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.tab === tabId);
  });

  activeTab = tabId;

  if (tabId === 'tab-plots') startPlotsTimer();
}

document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => showTab(btn.dataset.tab));
});

onDataUpdated(renderAllPlots);

renderAllPlots();