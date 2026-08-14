'use strict';

/* ══════════════════════════════════════════════════════════
   NOISE / DRIFT / NSA
   Leaf module (like main.js): only imports from main.js, never
   from kinetics.js or render.js, so both of those can import
   from here with no circular dependency.

   Ported from the "SPRm 2-D DATA SIMULATOR" reference file:

     artifact[pixel][t] =   nsaRate[pixel]  * Yns[t]           (non-specific adsorption)
                           + gDrift[pixel]  * driftCommon[t]   (baseline drift)
                           + sigma          * hash(i,j,t)      (per-pixel jitter)

   Yns[t] and driftCommon[t] are single shared time-curves (temporal —
   computed here). nsaRate/gDrift are per-pixel weight fields derived
   from the capacity field's edge-weight map via a Sobel gradient
   (spatial — computeSpatialWeights() is called by render.js whenever
   it rebuilds the capacity field). Per-pixel noise needs no shared
   state at all: it's a pure hash of (i, j, t, seed).
   ══════════════════════════════════════════════════════════ */

import { $, on, notifyDataUpdated, STACK_LEAD_BASELINE_SEC } from './main.js';

/* ── HARD-CODED PARAMETERS ─────────────────────────────────
   Not yet exposed in the UI (only on/off + a shared seed are).
   Calibrated in RU at Rmax≈1 — verified against our real grid/Cfun/
   capacity-field via a standalone harness before wiring this in.
   If per-cell-type Rmax grows much past ~1, these may need to scale
   with it (the reference file scales everything off a single RmaxD;
   we don't have one global Rmax anymore, so this is a flat default). */
const SIGMA_RU = 0.10;      // per-pixel jitter sd, RU

const A_N = 0.667, B_N = 0.667;      // NSA: uniform floor / edge-gradient weight, RU
const KA_NS = 4e3, KD_NS = 1e-4;     // NSA pseudo-kinetics (slow, low-affinity, accumulates)

const DRIFT_D = 3.5;        // drift asymptote, RU
const DRIFT_TAU = 500;      // exponential time constant, s
const SIGMA_OU = 0.02;      // OU kick size, RU s^-1/2
const THETA_OU = 0.005;     // OU mean-reversion rate, s^-1
const DECAY_OU = true;      // OU kicks share the exp(-t/tau) envelope
const A_D = 1, B_D = 0.15;  // drift: uniform (cancels under reference subtraction) / edge-gradient weight

/* ── seed + toggles ───────────────────────────────────────── */
function getNoiseSeed() {
  const el = $("noiseSeed");
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : 0;
}

export const isPixelNoiseOn = () => { const el = $("noiseOn");  return el ? el.checked : false; };
export const isDriftOn      = () => { const el = $("driftOn");  return el ? el.checked : false; };
export const isNSAOn        = () => { const el = $("NSAOn");    return el ? el.checked : false; };

function refreshAndNotify() {
  cache.key = null; // force curve rebuild (seed may have changed)
  notifyDataUpdated();
}

on("genNoiseSeed", "click", () => {
  const el = $("noiseSeed");
  if (el) el.value = Math.floor(Math.random() * 2 ** 32);
  refreshAndNotify();
});
on("noiseSeed", "input", refreshAndNotify);
["noiseOn", "driftOn", "NSAOn"].forEach(id => on(id, "change", () => notifyDataUpdated()));

/* ── seeded RNG + coordinate hash (mirrors render.js's mulberry32) ── */
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Uniform in (0,1) determined entirely by (i, j, t, seed) — not a sequential
// stream, so any pixel/frame can be recomputed independently and reproducibly,
// regardless of render order.
function hashU01(i, j, t, seed) {
  let x = (Math.imul(i, 0x9E3779B1) ^ Math.imul(j, 0x85EBCA77)
         ^ Math.imul(t, 0xC2B2AE3D) ^ seed) >>> 0;
  x ^= x >>> 16; x = Math.imul(x, 0x7FEB352D) >>> 0;
  x ^= x >>> 15; x = Math.imul(x, 0x846CA68B) >>> 0;
  x ^= x >>> 16;
  return (x >>> 0) / 4294967296;
}

function hashGauss(i, j, t, seed) {
  const u = 1 - hashU01(i, j, t, seed);
  const v = hashU01(i, j, t, (seed ^ 0x5BF03635) >>> 0);
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// Per-pixel jitter for one (pixel, local-frame) pair. Called from render.js's
// pixel loop only when the per-pixel-noise toggle is on.
export function pixelNoiseValue(i, j, timeIdx) {
  return SIGMA_RU * hashGauss(i, j, timeIdx, getNoiseSeed());
}

/* ── temporal curves: NSA basis (Yns) + common-mode drift ──────────
   Same serial-injection Cfun as kinetics.js (including the 5s lead
   baseline), re-derived here from the DOM so this stays a leaf module
   with no import from kinetics.js. Cached per (timing, concs, seed),
   sliced into per-concentration local-frame segments exactly like
   kinetics.js slices region traces, so indexing lines up with
   decodeFrame()'s (concIdx, timeIdx). */

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

  const key = JSON.stringify([tBase, tAssoc, tDissoc, concs, seed]);
  if (cache.key === key) return cache;

  const cyc     = STACK_LEAD_BASELINE_SEC + tAssoc + tDissoc;
  const nSpots  = concs.length;
  const total   = tBase + nSpots * cyc;
  const grid    = Array.from({ length: Math.round(total) + 1 }, (_, i) => i);
  const nFrames = Math.round(cyc) + 1;
  const concsM  = concs.map(c => c * 1e-9);
  const Cfun    = buildCfun(tBase, tAssoc, cyc, concsM);

  // NSA basis curve: slow, low-affinity, nearly-irreversible — ratchets up
  // across the run rather than tracking individual injections.
  const Yns = simRK4Scalar(grid, (y, C) => KA_NS * C * (1 - y) - KD_NS * y, 0, Cfun);

  // Common-mode drift: deterministic exponential ramp to asymptote D, plus an
  // Ornstein-Uhlenbeck wobble whose kicks shrink under the same envelope.
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

// Called once per frame render (not per pixel) — returns the two scalar
// artifact values shared by every pixel that frame.
export function getFrameArtifactScalars(concIdx, timeIdx) {
  const { YnsByConc, driftByConc } = ensureCurves();
  const yns      = YnsByConc[concIdx]   ? (YnsByConc[concIdx][timeIdx]   ?? 0) : 0;
  const driftVal = driftByConc[concIdx] ? (driftByConc[concIdx][timeIdx] ?? 0) : 0;
  return { yns, driftVal };
}

/* ── spatial weights: Sobel gradient of the capacity field's edge-weight
   map -> ghat -> nsaRate/gDrift. Pure function; render.js calls this once
   whenever it rebuilds the capacity field and caches the result alongside
   regionOf/edgeWeight. Computed for every pixel (not just cell-owned ones):
   both NSA and drift affect the whole sensor surface, not just where cells
   are — a_n/a_d are flat floors, b_n/b_d add extra at edges. ── */
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

  // Robust normaliser: 95th percentile over nonzero gradients (the max is a
  // one-pixel rasterisation artifact and isn't a stable normaliser).
  let p95;
  {
    const nz = [];
    for (let k = 0; k < P; k++) if (gradMag[k] > 1e-12) nz.push(gradMag[k]);
    nz.sort((a, b) => a - b);
    p95 = nz.length ? nz[Math.min(nz.length - 1, Math.floor(0.95 * nz.length))] : 1;
    if (!(p95 > 0)) p95 = 1;
  }

  const nsaRate = new Float32Array(P);
  const gDrift  = new Float32Array(P);
  for (let k = 0; k < P; k++) {
    const ghat = gradMag[k] / p95; // deliberately unclamped (can exceed 1 at strong edges)
    nsaRate[k] = A_N + B_N * ghat;
    gDrift[k]  = A_D + B_D * ghat;
  }
  return { nsaRate, gDrift };
}