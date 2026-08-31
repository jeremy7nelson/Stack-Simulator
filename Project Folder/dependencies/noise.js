'use strict';

import { $, on, notifyDataUpdated, STACK_LEAD_BASELINE_SEC } from './main.js';

const SIGMA_RU = 0.02;

const A_N_DEFAULT = 0.65, B_N_DEFAULT = 0.65;
const KA_NS_DEFAULT = 4e3, KD_NS_DEFAULT = 1e-4;
const B_D_DEFAULT = 0.15;

const DRIFT_D = 0.35;
const DRIFT_TAU = 500;
const SIGMA_OU = 0.02;
const THETA_OU = 0.005;
const DECAY_OU = true;
const A_D = 1;

function readNum(id, fallback) {
  const el = $(id);
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : fallback;
}

function getNsaParams() {
  return {
    a_n:  readNum("nsaBackground", A_N_DEFAULT),
    b_n:  readNum("nsaEdge", B_N_DEFAULT),
    kaNs: readNum("nsaBuildup", KA_NS_DEFAULT),
    kdNs: readNum("nsaFade", KD_NS_DEFAULT)
  };
}

function getDriftParams() {
  return { b_d: readNum("driftEdge", B_D_DEFAULT) };
}

export function getSpatialParamsKey() {
  const { a_n, b_n } = getNsaParams();
  const { b_d } = getDriftParams();
  return JSON.stringify([a_n, b_n, b_d]);
}

function getNoiseSeed() {
  const el = $("noiseSeed");
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : 0;
}

export const isPixelNoiseOn = () => { const el = $("noiseOn");  return el ? el.checked : false; };
export const isDriftOn      = () => { const el = $("driftOn");  return el ? el.checked : false; };
export const isNSAOn        = () => { const el = $("NSAOn");    return el ? el.checked : false; };

function refreshAndNotify() {
  cache.key = null;
  notifyDataUpdated();
}

on("genNoiseSeed", "click", () => {
  const el = $("noiseSeed");
  if (el) el.value = Math.floor(Math.random() * 2 ** 32);
  refreshAndNotify();
});
on("noiseSeed", "input", refreshAndNotify);

function toggleFieldGroup(fieldId, on_) {
  const el = $(fieldId);
  if (!el) return;
  el.style.opacity       = on_ ? "1" : ".45";
  el.style.pointerEvents = on_ ? "auto" : "none";
}

on("driftOn", "change", () => { toggleFieldGroup("driftField", isDriftOn()); notifyDataUpdated(); });
on("NSAOn",   "change", () => { toggleFieldGroup("nsaField", isNSAOn()); notifyDataUpdated(); });
on("noiseOn", "change", () => notifyDataUpdated());

["nsaBackground", "nsaEdge", "driftEdge"].forEach(id => on(id, "input", () => notifyDataUpdated()));

["nsaBuildup", "nsaFade"].forEach(id => on(id, "input", refreshAndNotify));

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function hashU01(i, j, t, seed) {
  let x = (Math.imul(i, 0x9E3779B1) ^ Math.imul(j, 0x85EBCA77)
         ^ Math.imul(t, 0xC2B2AE3D) ^ seed) >>> 0;
  x ^= x >>> 16; x = Math.imul(x, 0x7FEB352D) >>> 0;
  x ^= x >>> 15; x = Math.imul(x, 0x846CA68B) >>> 0;
  x ^= x >>> 16;
  return (x >>> 0) / 4294967296;
}

const HASH_GAUSS_OFFSETS = [0x5BF03635, 0x1B873593, 0x3AC5D673, 0x27D4EB2F];

function hashGauss(i, j, t, seed) {
  let sum = 0;
  for (let k = 0; k < HASH_GAUSS_OFFSETS.length; k++) {
    sum += hashU01(i, j, t, (seed ^ HASH_GAUSS_OFFSETS[k]) >>> 0);
  }
  return (sum - 2) * 1.7320508075688772;
}

export function pixelNoiseValue(i, j, timeIdx) {
  return SIGMA_RU * hashGauss(i, j, timeIdx, getNoiseSeed());
}

const parseConcs = str =>
  (str || "").split(/[\s,;]+/).map(Number).filter(v => Number.isFinite(v) && v > 0);

function buildCfun(tBase, tAssoc, cyc, concsM) {
  const nSpots = concsM.length;
  return t => {
    if (t < tBase) return 0;
    let k = Math.floor((t - tBase) / cyc);
    if (k >= nSpots) k = nSpots - 1;
    const inCyc = (t - tBase) - k * cyc;
    if (inCyc < STACK_LEAD_BASELINE_SEC) return 0;
    return (inCyc - STACK_LEAD_BASELINE_SEC) < tAssoc ? concsM[k] : 0;
  };
}

function simRK4Scalar(grid, deriv, y0, Cfun) {
  let y = y0;
  const out = [y];
  for (let i = 1; i < grid.length; i++) {
    const t0 = grid[i - 1], h = grid[i] - grid[i - 1];
    const C0 = Cfun(t0), Cm = Cfun(t0 + h / 2), C1 = Cfun(t0 + h);
    const k1 = deriv(y, C0);
    const k2 = deriv(y + k1 * h / 2, Cm);
    const k3 = deriv(y + k2 * h / 2, Cm);
    const k4 = deriv(y + k3 * h, C1);
    y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    out.push(y);
  }
  return out;
}

let cache = { key: null, YnsByConc: null, driftByConc: null };

function ensureCurves() {
  const tBase   = +$("tBase").value;
  const tAssoc  = +$("tAssoc").value;
  const tDissoc = +$("tDissoc").value;
  const concs   = parseConcs($("concSeries") ? $("concSeries").value : "");
  const seed    = getNoiseSeed();
  const { kaNs, kdNs } = getNsaParams();

  const key = JSON.stringify([tBase, tAssoc, tDissoc, concs, seed, kaNs, kdNs]);
  if (cache.key === key) return cache;

  const cyc     = STACK_LEAD_BASELINE_SEC + tAssoc + tDissoc;
  const nSpots  = concs.length;
  const total   = tBase + nSpots * cyc;
  const grid    = Array.from({ length: Math.round(total) + 1 }, (_, i) => i);
  const nFrames = Math.round(cyc) + 1;
  const concsM  = concs.map(c => c * 1e-9);
  const Cfun    = buildCfun(tBase, tAssoc, cyc, concsM);

  const Yns = simRK4Scalar(grid, (y, C) => kaNs * C * (1 - y) - kdNs * y, 0, Cfun);

  const rng = mulberry32(seed ^ 0x2545F491);
  function gaussStream() {
    let u = 0, v = 0;
    while (u === 0) u = rng();
    while (v === 0) v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  const driftCommon = new Float64Array(grid.length);
  {
    let w = 0;
    for (let t = 0; t < grid.length; t++) {
      const envelope = DECAY_OU ? Math.exp(-t / DRIFT_TAU) : 1;
      w += -THETA_OU * w + SIGMA_OU * envelope * gaussStream();
      driftCommon[t] = DRIFT_D * (1 - Math.exp(-t / DRIFT_TAU)) + w;
    }
  }

  const YnsByConc = concs.map((_, k) => {
    const startIdx = Math.round(tBase + k * cyc);
    return Yns.slice(startIdx, startIdx + nFrames);
  });
  const driftByConc = concs.map((_, k) => {
    const startIdx = Math.round(tBase + k * cyc);
    return Array.from(driftCommon.slice(startIdx, startIdx + nFrames));
  });

  cache = { key, YnsByConc, driftByConc };
  return cache;
}

export function getFrameArtifactScalars(concIdx, timeIdx) {
  const { YnsByConc, driftByConc } = ensureCurves();
  const yns      = YnsByConc[concIdx]   ? (YnsByConc[concIdx][timeIdx]   ?? 0) : 0;
  const driftVal = driftByConc[concIdx] ? (driftByConc[concIdx][timeIdx] ?? 0) : 0;
  return { yns, driftVal };
}

export function computeSpatialWeights(edgeWeight, m, n) {
  const P = m * n;
  const gradMag = new Float32Array(P);
  for (let i = 1; i < m - 1; i++) {
    for (let j = 1; j < n - 1; j++) {
      const k = i * n + j;
      const gi = ( edgeWeight[k + n - 1] + 2 * edgeWeight[k + n] + edgeWeight[k + n + 1]
                 - edgeWeight[k - n - 1] - 2 * edgeWeight[k - n] - edgeWeight[k - n + 1] ) / 8;
      const gj = ( edgeWeight[k - n + 1] + 2 * edgeWeight[k + 1] + edgeWeight[k + n + 1]
                 - edgeWeight[k - n - 1] - 2 * edgeWeight[k - 1] - edgeWeight[k + n - 1] ) / 8;
      gradMag[k] = Math.sqrt(gi * gi + gj * gj);
    }
  }

  let p95;
  {
    const nz = [];
    for (let k = 0; k < P; k++) if (gradMag[k] > 1e-12) nz.push(gradMag[k]);
    nz.sort((a, b) => a - b);
    p95 = nz.length ? nz[Math.min(nz.length - 1, Math.floor(0.95 * nz.length))] : 1;
    if (!(p95 > 0)) p95 = 1;
  }

  const { a_n, b_n } = getNsaParams();
  const { b_d } = getDriftParams();

  const nsaRate = new Float32Array(P);
  const gDrift  = new Float32Array(P);
  for (let k = 0; k < P; k++) {
    const ghat = gradMag[k] / p95;
    nsaRate[k] = a_n + b_n * ghat;
    gDrift[k]  = A_D + b_d * ghat;
  }
  return { nsaRate, gDrift };
}