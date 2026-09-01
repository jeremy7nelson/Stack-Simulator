'use strict';

import { $, on, notifyDataUpdated, STACK_LEAD_BASELINE_SEC } from './main.js';

const A_N_DEFAULT = 0.65, B_N_DEFAULT = 0.65;
const KA_NS_DEFAULT = 4e3, KD_NS_DEFAULT = 1e-4;
const B_D_DEFAULT = 0.15;
const D_DEFAULT = 0.35;
const TAU_DEFAULT = 500;
const SIGMA_OU_DEFAULT = 0.002;
const THETA_OU_DEFAULT = 0.005;
const SIGMA_DEFAULT = 0.01;
const LAMBDA_IN_DEFAULT = 8;
const LAMBDA_OUT_DEFAULT = 2;
const DECAY_OU = true;
export const A_D = 1;

function readNum(id, fallback) {
  const el = $(id);
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : fallback;
}

export function getNsaParams() {
  return {
    a_n:  readNum("nsaBackground", A_N_DEFAULT),
    b_n:  readNum("nsaEdge", B_N_DEFAULT),
    kaNs: readNum("nsaBuildup", KA_NS_DEFAULT),
    kdNs: readNum("nsaFade", KD_NS_DEFAULT)
  };
}

export function getDriftParams() {
  return {
    b_d: readNum("driftEdge", B_D_DEFAULT),
    D: readNum("Dmultiplier", D_DEFAULT),
    tau: readNum("tau", TAU_DEFAULT),
    sigmaOU: readNum("sigmaOU", SIGMA_OU_DEFAULT),
    thetaOU: readNum("thetaOU", THETA_OU_DEFAULT)
  };
}

function getPixelNoiseSigma() {
  return readNum("sigma", SIGMA_DEFAULT);
}

function getEdgeParams() {
  return {
    lambdaIn: readNum("lambdaIn", LAMBDA_IN_DEFAULT),
    lambdaOut: readNum("lambdaOut", LAMBDA_OUT_DEFAULT)
  };
}

export function getSpatialParamsKey() {
  const { a_n, b_n } = getNsaParams();
  const { b_d } = getDriftParams();
  const { lambdaIn, lambdaOut } = getEdgeParams();
  return JSON.stringify([a_n, b_n, b_d, lambdaIn, lambdaOut]);
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

function updateEdgeGeometryVisibility() {
  toggleFieldGroup("edgeGeometryField", isDriftOn() || isNSAOn());
}

on("driftOn", "change", () => { toggleFieldGroup("driftField", isDriftOn()); updateEdgeGeometryVisibility(); notifyDataUpdated(); });
on("NSAOn",   "change", () => { toggleFieldGroup("nsaField", isNSAOn()); updateEdgeGeometryVisibility(); notifyDataUpdated(); });
on("noiseOn", "change", () => { toggleFieldGroup("pixelNoiseField", isPixelNoiseOn()); notifyDataUpdated(); });

["nsaBackground", "nsaEdge", "driftEdge", "lambdaIn", "lambdaOut", "sigma"]
  .forEach(id => on(id, "input", () => notifyDataUpdated()));

["nsaBuildup", "nsaFade", "Dmultiplier", "tau", "sigmaOU", "thetaOU"]
  .forEach(id => on(id, "input", refreshAndNotify));

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

export function getPixelNoiseParams() {
  return { sigma: getPixelNoiseSigma(), seed: getNoiseSeed() };
}

export function pixelNoiseValue(i, j, timeIdx, sigma, seed) {
  return sigma * hashGauss(i, j, timeIdx, seed);
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
  const { D, tau, sigmaOU, thetaOU } = getDriftParams();

  const key = JSON.stringify([tBase, tAssoc, tDissoc, concs, seed, kaNs, kdNs, D, tau, sigmaOU, thetaOU]);
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
    let sum = 0;
    for (let i = 0; i < 4; i++) sum += rng();
    return (sum - 2) * 1.7320508075688772;
  }
  const driftCommon = new Float64Array(grid.length);
  {
    let w = 0;
    for (let t = 0; t < grid.length; t++) {
      const envelope = DECAY_OU ? Math.exp(-t / tau) : 1;
      w += -thetaOU * w + sigmaOU * envelope * gaussStream();
      driftCommon[t] = D * (1 - Math.exp(-t / tau)) + w;
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

  cache = { key, YnsByConc, driftByConc, grid, Yns, driftCommon };
  return cache;
}

export function getFrameArtifactScalars(concIdx, timeIdx) {
  const { YnsByConc, driftByConc } = ensureCurves();
  const yns      = YnsByConc[concIdx]   ? (YnsByConc[concIdx][timeIdx]   ?? 0) : 0;
  const driftVal = driftByConc[concIdx] ? (driftByConc[concIdx][timeIdx] ?? 0) : 0;
  return { yns, driftVal };
}

export function getRawCurves() {
  const { grid, Yns, driftCommon } = ensureCurves();
  return { grid, Yns, driftCommon };
}

const EDT_INF = 1e20;

function edt1d(f, N, d, v, z) {
  let k = 0;
  v[0] = 0; z[0] = -Infinity; z[1] = Infinity;
  for (let q = 1; q < N; q++) {
    let s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
    while (s <= z[k]) {
      k--;
      s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
    }
    k++; v[k] = q; z[k] = s; z[k + 1] = Infinity;
  }
  k = 0;
  for (let q = 0; q < N; q++) {
    while (z[k + 1] < q) k++;
    d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
  }
}

function edt2d(mask, m, n) {
  const L = Math.max(m, n);
  const f = new Float64Array(L), d = new Float64Array(L);
  const v = new Int32Array(L + 1), z = new Float64Array(L + 1);
  const out = new Float64Array(m * n);
  for (let k = 0; k < m * n; k++) out[k] = mask[k] ? 0 : EDT_INF;
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) f[i] = out[i * n + j];
    edt1d(f, m, d, v, z);
    for (let i = 0; i < m; i++) out[i * n + j] = d[i];
  }
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) f[j] = out[i * n + j];
    edt1d(f, n, d, v, z);
    for (let j = 0; j < n; j++) out[i * n + j] = d[j];
  }
  return out;
}

export function computeSpatialWeights(covered, m, n) {
  const P = m * n;
  const inverted = new Uint8Array(P);
  for (let k = 0; k < P; k++) inverted[k] = covered[k] ? 0 : 1;

  const dToCell = edt2d(covered, m, n);
  const dToVoid = edt2d(inverted, m, n);

  const signedDist = new Float32Array(P);
  for (let k = 0; k < P; k++) {
    signedDist[k] = covered[k] ? (Math.sqrt(dToVoid[k]) - 0.5)
                               : -(Math.sqrt(dToCell[k]) - 0.5);
  }

  const { lambdaIn, lambdaOut } = getEdgeParams();
  const w = new Float32Array(P);
  for (let k = 0; k < P; k++) {
    const d = signedDist[k];
    w[k] = (d >= 0) ? Math.exp(-d / lambdaIn) : Math.exp(d / lambdaOut);
  }

  const { a_n, b_n } = getNsaParams();
  const { b_d } = getDriftParams();

  const nsaRate = new Float32Array(P);
  const gDrift  = new Float32Array(P);
  for (let k = 0; k < P; k++) {
    nsaRate[k] = a_n + b_n * w[k];
    gDrift[k]  = A_D + b_d * w[k];
  }
  return { nsaRate, gDrift, w, signedDist };
}