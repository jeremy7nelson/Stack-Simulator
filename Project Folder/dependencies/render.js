'use strict';

import {
  $, on, state, onDataUpdated, fmtConc, IMG_W, IMG_H, MAX16
} from './main.js';
import { simulate } from './kinetics.js';

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function buildCellTypeBag(frequencies) {
  const bag = [];
  frequencies.forEach((freq, idx) => {
    for (let i = 0; i < Math.round(freq); i++) bag.push(idx);
  });
  return bag.length ? bag : [0];
}

const EDGE_FLOOR = 0.15;

function generateCapacityField({
  m, n,
  cellTypeBag,
  targetConfluence,
  rMin, rMax,
  maxCircles = 5000,
  seed,
  edgeDominance
}) {
  const rng = mulberry32(seed);
  const noOverlap = !allowOverlap();

  const covered = new Uint8Array(m * n);
  let coveredCount = 0;

  const regionOf   = new Int16Array(m * n).fill(-1);
  const edgeWeight  = new Float32Array(m * n).fill(1);
  const bestPower   = new Float32Array(m * n).fill(Infinity);
  const circles = [];

  const MAX_OVERLAP_RETRIES = 30;

  function placeOneCircle() {
    let ci, cj, r, regionIdx, tries = 0;
    do {
      ci = rng() * m;
      cj = rng() * n;
      r  = rMin + rng() * (rMax - rMin);
      regionIdx = cellTypeBag[(rng() * cellTypeBag.length) | 0];
      tries++;
    } while (noOverlap && circleOverlapsAny(ci, cj, r, circles) && tries < MAX_OVERLAP_RETRIES);

    if (noOverlap && circleOverlapsAny(ci, cj, r, circles)) return;

    const r2 = r * r;
    const iLo = Math.max(0, Math.floor(ci - r)), iHi = Math.min(m - 1, Math.ceil(ci + r));
    const jLo = Math.max(0, Math.floor(cj - r)), jHi = Math.min(n - 1, Math.ceil(cj + r));
    for (let i = iLo; i <= iHi; i++) {
      for (let j = jLo; j <= jHi; j++) {
        const di = i - ci, dj = j - cj;
        const d2 = di * di + dj * dj;
        if (d2 <= r2) {
          const k = i * n + j;
          if (covered[k] === 0) { covered[k] = 1; coveredCount++; }

          // Power-diagram (Laguerre) ownership: whichever circle has the
          // smallest d^2 - r^2 at this pixel wins it, so a bigger disc can
          // rightfully claim pixels nearer a smaller disc's own centre.
          // Order-independent, so overlapping circles never simply erase
          // whichever was placed first (or last).
          const power = d2 - r2;
          if (power < bestPower[k]) {
            bestPower[k] = power;
            regionOf[k]  = regionIdx;
            // Edge dominance: brightest at the disc's own rim, dimmest at its
            // own centre — measured from THIS circle's outer boundary only,
            // so a cell fully enclosed by neighbours never shows a bright rim.
            // When disabled, every owned pixel stays at full weight (flat).
            edgeWeight[k] = edgeDominance
              ? EDGE_FLOOR + (1 - EDGE_FLOOR) * (d2 / r2)
              : 1;
          }
        }
      }
    }
    circles.push({ ci, cj, r, regionIdx });
  }

  let iter = 0;
  while (coveredCount / (m * n) < targetConfluence && iter < maxCircles) {
    placeOneCircle();
    iter++;
  }

  return { regionOf, edgeWeight, circles, achievedConfluence: coveredCount / (m * n) };
}

function circleOverlapsAny(ci, cj, r, circles) {
  for (const c of circles) {
    const dx = ci - c.ci, dy = cj - c.cj;
    const minDist = r + c.r;
    if (dx * dx + dy * dy < minDist * minDist) return true;
  }
  return false;
}

function allowOverlap() {
  const el = $("overlap");
  return el ? el.checked : true;
}

function allowEdgeDominance() {
  const el = $("dominance");
  return el ? el.checked : false;
}

on("overlap", "change", () => refreshCapacityFieldIfNeeded());
on("dominance", "change", () => refreshCapacityFieldIfNeeded());

function getSeed() {
  const el = $("inputSeed");
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : 0;
}

on("genSeed", "click", () => {
  const el = $("inputSeed");
  if (el) el.value = Math.floor(Math.random() * 2 ** 32);
  refreshCapacityFieldIfNeeded();
});

on("inputSeed", "input", refreshCapacityFieldIfNeeded);

const CIRCLE_R_MIN = 15, CIRCLE_R_MAX = 25;

let lastFieldSeed, lastFieldConfluence, lastFieldOverlap, lastFieldFreqKey, lastFieldDominance;

function getTargetConfluence() {
  const el = $("Confluency");
  const pct = el ? +el.value : NaN;
  return (Number.isFinite(pct) ? pct : 0) / 100;
}

function refreshCapacityFieldIfNeeded() {
  const seed = getSeed();
  const targetConfluence = getTargetConfluence();
  const overlap = allowOverlap();
  const dominance = allowEdgeDominance();
  const freqKey = (state.cellFrequencies || []).join(",");

  const unchanged = state.capacityField
    && seed === lastFieldSeed
    && targetConfluence === lastFieldConfluence
    && overlap === lastFieldOverlap
    && dominance === lastFieldDominance
    && freqKey === lastFieldFreqKey;
  if (unchanged) return;

  state.capacityField = generateCapacityField({
    m: IMG_H, n: IMG_W,
    cellTypeBag: buildCellTypeBag(state.cellFrequencies || [100]),
    targetConfluence,
    rMin: CIRCLE_R_MIN, rMax: CIRCLE_R_MAX,
    seed,
    edgeDominance: dominance
  });

  lastFieldSeed = seed;
  lastFieldConfluence = targetConfluence;
  lastFieldOverlap = overlap;
  lastFieldDominance = dominance;
  lastFieldFreqKey = freqKey;
}

on("Confluency", "input", refreshCapacityFieldIfNeeded);

function updateStackImage() {
  refreshCapacityFieldIfNeeded();

  const results = $("results");
  if (!state.parsed || !state.parsed.regions || state.parsed.nSpots === 0 || state.parsed.nFrames === 0) {
    if (results) results.style.display = "none";
    return;
  }

  const totalFrames = state.parsed.nFrames * state.parsed.nSpots;
  const slider = $("frame-slider");
  if (slider) {
    slider.max = totalFrames - 1;
    if (+slider.value > totalFrames - 1) slider.value = 0;
  }

  const totalEl = $("total-frames");
  if (totalEl) totalEl.textContent =
    `${totalFrames}  (${state.parsed.nFrames} time points × ${state.parsed.nSpots} concentrations)`;
  if (results) results.style.display = "block";
  renderPreview(slider ? +slider.value : 0);
}

export function decodeFrame(globalFrame) {
  return {
    concIdx: Math.floor(globalFrame / state.parsed.nFrames),
    timeIdx: globalFrame % state.parsed.nFrames
  };
}

function regionBrightness16(rg, concIdx, timeIdx) {
  if (concIdx >= rg.traces.length) return 0;
  const { globalMin, globalMax } = state.parsed;
  const denom = (globalMax - globalMin) || 1;
  const ru    = rg.traces[concIdx][timeIdx];
  const norm  = Math.max(0, (ru - globalMin) / denom);
  return Math.round(norm * MAX16);
}

function compositeCapacityField(mat, b16ByRegion) {
  const field = state.capacityField;
  if (!field) return;
  const { regionOf, edgeWeight } = field;
  for (let k = 0; k < mat.length; k++) {
    const rIdx = regionOf[k];
    let v = rIdx >= 0 ? Math.round((b16ByRegion[rIdx] ?? 0) * edgeWeight[k]) : 0;
    if (v < 0) v = 0; else if (v > MAX16) v = MAX16;
    mat[k] = v;
  }
}

export function getMatrix16(globalFrame) {
  const mat = new Uint16Array(IMG_H * IMG_W);
  if (!state.parsed || !state.parsed.regions) return mat;
  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  const b16ByRegion = state.parsed.regions.map(rg => regionBrightness16(rg, concIdx, timeIdx));
  compositeCapacityField(mat, b16ByRegion);
  return mat;
}

export function getMatrix16ForBrightness(b16) {
  const mat = new Uint16Array(IMG_H * IMG_W);
  const field = state.capacityField;
  if (!field) return mat;
  const { regionOf, edgeWeight } = field;
  for (let k = 0; k < mat.length; k++) {
    if (regionOf[k] < 0) continue;
    let v = Math.round(b16 * edgeWeight[k]);
    if (v < 0) v = 0; else if (v > MAX16) v = MAX16;
    mat[k] = v;
  }
  return mat;
}

function renderPreview(globalFrame) {
  if (!state.parsed || !state.parsed.regions) return;
  const canvas = $("img-canvas");
  if (!canvas) return;
  canvas.width = IMG_W; canvas.height = IMG_H;

  const mat = getMatrix16(globalFrame);
  const ctx = canvas.getContext("2d");
  const img = ctx.createImageData(IMG_W, IMG_H);
  for (let i = 0; i < mat.length; i++) {
    const g = mat[i] >> 8;
    const j = i * 4;
    img.data[j] = g; img.data[j + 1] = g; img.data[j + 2] = g; img.data[j + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);

  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  const frameVal = $("frame-val");
  if (frameVal) frameVal.textContent = globalFrame;

  const concTag = $("conc-tag");
  if (concTag) {
    const cLabel = concIdx < state.parsed.concs.length ? fmtConc(state.parsed.concs[concIdx])
                                                 : `spot ${concIdx + 1}`;
    concTag.textContent = `${cLabel}  —  time point ${timeIdx} / ${state.parsed.nFrames - 1}`;
  }
}

export function findPeakInjectionFrame() {
  const tA = 0;
  const tD = +$("tAssoc").value;
  const { regions, nFrames, nSpots, times } = state.parsed;

  let bestFrame = 0, bestRU = -Infinity;
  for (let f = 0; f < nFrames * nSpots; f++) {
    const { concIdx, timeIdx } = decodeFrame(f);
    const t = times[timeIdx];
    if (t < tA || t >= tD) continue;
    let ru = 0;
    for (const rg of regions)
      if (concIdx < rg.traces.length) ru = Math.max(ru, rg.traces[concIdx][timeIdx]);
    if (ru > bestRU) { bestRU = ru; bestFrame = f; }
  }
  return { frame: bestFrame, ru: bestRU };
}

on("frame-slider", "input", function () { renderPreview(+this.value); });

let playing = false, raf = null;

function stepFrame() {
  if (!playing || !state.parsed) { playing = false; return; }
  const slider = $("frame-slider");
  const totalFrames = state.parsed.nFrames * state.parsed.nSpots;
  if (!slider || totalFrames <= 0) { playing = false; return; }

  const t = (+slider.value + 4) % totalFrames;
  slider.value = t;
  renderPreview(t);
  raf = requestAnimationFrame(stepFrame);
}

on("play", "click", () => {
  const btn = $("play");
  playing = !playing;
  if (btn) btn.textContent = playing ? "Pause" : "Play";
  if (playing) stepFrame(); else cancelAnimationFrame(raf);
});

onDataUpdated(updateStackImage);

simulate();

refreshCapacityFieldIfNeeded();