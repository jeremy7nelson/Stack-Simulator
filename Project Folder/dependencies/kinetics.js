'use strict';

import {
  $, on, setStatus, state, notifyDataUpdated,
  STACK_LEAD_BASELINE_SEC
} from './main.js';

const vadd   = (a, b) => a.map((v, i) => v + b[i]);
const vscale = (a, s) => a.map(v => v * s);
const vsum   = a => a.reduce((x, y) => x + y, 0);

const simRK4 = (grid, deriv, y0, Cfun) => {
  const out = []; let y = [...y0]; out.push(vsum(y));
  for (let i = 1; i < grid.length; i++) {
    const [t0, h] = [grid[i - 1], grid[i] - grid[i - 1]];
    const [C0, Cm, C1] = [Cfun(t0), Cfun(t0 + h / 2), Cfun(t0 + h)];
    const k1 = deriv(y, C0);
    const k2 = deriv(vadd(y, vscale(k1, h / 2)), Cm);
    const k3 = deriv(vadd(y, vscale(k2, h / 2)), Cm);
    const k4 = deriv(vadd(y, vscale(k3, h)),     C1);
    y = vadd(y, vscale(vadd(vadd(k1, vscale(k2, 2)), vadd(vscale(k3, 2), k4)), h / 6));
    out.push(vsum(y));
  }
  return out;
};

function makeDeriv(model, Rmax, gv) {
  gv = gv || (id => +$(id).value);

  if (model === "langmuir") {
    const ka = gv("ka"), kd = gv("kd");
    return {
      size: 1,
      deriv: (y, C) => { const R = y[0]; return [ka * C * (Rmax - R) - kd * R]; },
      fluxCoef: (y) => ({ a: ka * (Rmax - y[0]), b: -kd * y[0] })
    };
  }
  if (model === "hetLigand") {
    const ka1 = gv("hetka1"), kd1 = gv("hetkd1"),
          ka2 = gv("hetka2"), kd2 = gv("hetkd2"), Rmax2 = gv("Rmax2");
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1];
        return [ka1 * C * (Rmax - R1) - kd1 * R1, ka2 * C * (Rmax2 - R2) - kd2 * R2]; },
      fluxCoef: (y) => ({ a: ka1 * (Rmax - y[0]) + ka2 * (Rmax2 - y[1]),
                          b: -(kd1 * y[0] + kd2 * y[1]) })
    };
  }
  if (model === "bivAnalyte") {
    const ka1 = gv("bivka1"), kd1 = gv("bivkd1"),
          ka2 = gv("bivka2"), kd2 = gv("bivkd2");
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1], free = Rmax - R1 - 2 * R2;
        return [2 * ka1 * C * free - kd1 * R1 - ka2 * R1 * free + 2 * kd2 * R2,
                2 * ka2 * R1 * free - 2 * kd2 * R2]; },
      fluxCoef: (y) => { const free = Rmax - y[0] - 2 * y[1];
        return { a: 2 * ka1 * free, b: -kd1 * y[0] }; }
    };
  }
  const ka1 = gv("ka1"), kd1 = gv("kd1"), ka2 = gv("ka2"), kd2 = gv("kd2");
  return {
    size: 2,
    deriv: (y, C) => { const AB = y[0], ABs = y[1], free = Rmax - AB - ABs;
      return [ka1 * C * free - kd1 * AB - ka2 * AB + kd2 * ABs, ka2 * AB - kd2 * ABs]; },
    fluxCoef: (y) => { const free = Rmax - y[0] - y[1];
      return { a: ka1 * free, b: -kd1 * y[0] }; }
  };
}

function withTransport(base, kt) {
  return {
    size: base.size,
    deriv: (y, C) => {
      const { a, b } = base.fluxCoef(y);
      let Cs = (kt * C - b) / (kt + a);
      if (!isFinite(Cs) || Cs < 0) Cs = 0;
      return base.deriv(y, Cs);
    }
  };
}

const parseConcs = str =>
  (str || "").split(/[\s,;]+/).map(Number).filter(v => Number.isFinite(v) && v > 0);

const MODEL_HINTS = {
  langmuir:   "Simplest case: one analyte binding one immobilised ligand.",
  hetLigand:  "Two available binding sites with two completely independent dynamics.",
  twostate:   "Binding followed by a conformational change that locks the complex — note the slow dissociation.",
  bivAnalyte: "The analyte may, at sufficient density, bind two membrane receptors simultaneously."
};

const MODEL_INPUT_IDS = [
  "model", "ka", "kd", "ka1", "kd1", "ka2", "kd2",
  "hetka1", "hetkd1", "hetka2", "hetkd2",
  "bivka1", "bivkd1", "bivka2", "bivkd2",
  "kt", "Rmax", "Rmax2"
];

const MODEL_FIELD_IDS = [...MODEL_INPUT_IDS, "mtlOn"];

function defaultCellModelConfig() {
  return {
    model: "langmuir",
    ka: 1e6, kd: 1e-3,
    ka1: 1e6, kd1: 1e-2, ka2: 5e-3, kd2: 1e-3,
    hetka1: 1e6, hetkd1: 1e-2, hetka2: 1e5, hetkd2: 1e-3,
    bivka1: 1e6, bivkd1: 1e-2, bivka2: 5e-3, bivkd2: 1e-3,
    kt: 1e8, Rmax: 1, Rmax2: 1,
    mtlOn: false
  };
}

let cellModelConfigs = [];
let regionCache = [];
let lastSharedKey = null;

function saveCellModelConfig(idx) {
  if (idx == null || idx < 0) return;
  const cfg = {};
  MODEL_FIELD_IDS.forEach(id => {
    const el = $(id);
    if (!el) return;
    cfg[id] = el.type === "checkbox" ? el.checked : el.value;
  });
  cellModelConfigs[idx] = cfg;
}

function loadCellModelConfig(idx) {
  const cfg = cellModelConfigs[idx] || defaultCellModelConfig();
  MODEL_FIELD_IDS.forEach(id => {
    const el = $(id);
    if (!el || !(id in cfg)) return;
    if (el.type === "checkbox") el.checked = cfg[id];
    else el.value = cfg[id];
  });
  setModelVisibility();
  const ktField = $("ktField");
  const mtlEl   = $("mtlOn");
  if (ktField) ktField.style.display = (mtlEl && mtlEl.checked) ? "" : "none";
}

function sharedParamsKey(tBase, tAssoc, tDissoc, concs) {
  return JSON.stringify([tBase, tAssoc, tDissoc, concs]);
}

export function simulate() {
  const tBase   = +$("tBase").value;
  const tAssoc  = +$("tAssoc").value;
  const tDissoc = +$("tDissoc").value;
  const cyc     = STACK_LEAD_BASELINE_SEC + tAssoc + tDissoc;

  const concs  = parseConcs($("concSeries").value);
  const nSpots = concs.length;
  const total  = tBase + nSpots * cyc;

  const grid   = Array.from({ length: Math.round(total) + 1 }, (_, i) => i);
  const concsM = concs.map(Cnm => Cnm * 1e-9);

  const nFrames   = Math.round(cyc) + 1;
  const localGrid = Array.from({ length: nFrames }, (_, i) => i);

  const sharedKey     = sharedParamsKey(tBase, tAssoc, tDissoc, concs);
  const sharedChanged = sharedKey !== lastSharedKey;
  lastSharedKey = sharedKey;

  let gMin = Infinity, gMax = -Infinity;

  const regions = cellModelConfigs.map((cfg, idx) => {
    const key    = JSON.stringify(cfg);
    const cached = regionCache[idx];

    let region;
    if (!sharedChanged && cached && cached.key === key) {
      region = cached.region;
    } else {
      const Rmax   = +cfg.Rmax || 0;
      const base   = makeDeriv(cfg.model, Rmax, id => +cfg[id]);
      const engine = cfg.mtlOn ? withTransport(base, +cfg.kt || 0) : base;

      const Cfun = t => {
        if (t < tBase) return 0;
        let k = Math.floor((t - tBase) / cyc);
        if (k >= nSpots) k = nSpots - 1;
        const inCyc = (t - tBase) - k * cyc;
        if (inCyc < STACK_LEAD_BASELINE_SEC) return 0;
        return (inCyc - STACK_LEAD_BASELINE_SEC) < tAssoc ? concsM[k] : 0;
      };

      const Y = simRK4(grid, engine.deriv, new Array(engine.size).fill(0), Cfun);

      const traces = concs.map((_, k) => {
        const startIdx = Math.round(tBase + k * cyc);
        return Y.slice(startIdx, startIdx + nFrames);
      });

      region = { idx: idx + 1, traces };
      regionCache[idx] = { key, region };
    }

    region.traces.forEach(tr => tr.forEach(v => {
      if (v < gMin) gMin = v;
      if (v > gMax) gMax = v;
    }));
    return region;
  });

  regionCache.length = cellModelConfigs.length;

  if (!isFinite(gMin)) { gMin = 0; gMax = 0; }

  state.parsed          = { times: new Float64Array(localGrid), grid: localGrid, nFrames, nSpots, concs,
                             globalMin: gMin, globalMax: gMax, regions };
  state.lastData         = { grid: localGrid };
  state.cellFrequencies  = cellFrequencies;

  setStatus(`${nFrames} pts × ${nSpots} conc · ${regions.length} cell type${regions.length > 1 ? "s" : ""} · serial, no regen`);

  notifyDataUpdated();
}

function setModelVisibility() {
  const m = $("model").value;
  const groupMap = {
    simple:     m === "langmuir",
    twostate:   m === "twostate",
    bivAnalyte: m === "bivAnalyte",
    hetLigand:  m === "hetLigand",
  };
  document.querySelectorAll("[data-group]").forEach(el => {
    el.style.display = groupMap[el.dataset.group] ? "" : "none";
  });
  const hint = $("modelHint");
  if (hint) hint.textContent = MODEL_HINTS[m] ?? "";
}

function genDilution() {
  const [top, f, n] = ["dilTop","dilFactor","dilN"].map(id => +$(id).value);
  const pts = Array.from({ length: Math.max(1, Math.round(n)) }, (_, i) =>
    +(top / f ** i).toPrecision(4)
  );
  $("concSeries").value = pts.reverse().join(", ");
  simulate();
}

let cellFrequencies = [];
let selectedCellIdx = 0;

function evenSplit(n) {
  const base = Math.floor(100 / n);
  const arr = new Array(n).fill(base);
  arr[n - 1] += 100 - base * n;
  return arr;
}

function rebuildCellTypes() {
  const n = Math.max(1, Math.round(+$("regionCount").value) || 1);

  saveCellModelConfig(selectedCellIdx);

  cellFrequencies = evenSplit(n);
  while (cellModelConfigs.length < n) cellModelConfigs.push(defaultCellModelConfig());
  cellModelConfigs.length = n;

  selectedCellIdx = Math.min(selectedCellIdx, n - 1);

  const select = $("regionSelect");
  if (select) {
    select.innerHTML = "";
    for (let i = 0; i < n; i++) {
      const opt = document.createElement("option");
      opt.value = i;
      opt.textContent = `Model #${i + 1}`;
      select.appendChild(opt);
    }
    select.selectedIndex = selectedCellIdx;
  }

  loadCellModelConfig(selectedCellIdx);
  refreshFrequencyField();
}

function refreshFrequencyField() {
  const field = $("frequncy");
  if (!field) return;
  const usedByOthers = cellFrequencies.reduce(
    (sum, f, i) => i === selectedCellIdx ? sum : sum + f, 0
  );
  field.min = 0;
  field.max = Math.max(0, 100 - usedByOthers);
  field.value = cellFrequencies[selectedCellIdx] ?? 0;
}

MODEL_INPUT_IDS.forEach(id => {
  const el = $(id);
  if (!el) { console.warn("Missing element:", id); return; }
  el.addEventListener("input", () => {
    if (id === "model") setModelVisibility();
    saveCellModelConfig(selectedCellIdx);
    simulate();
  });
});

on("mtlOn", "change", () => {
  const kt = $("ktField");
  if (kt) kt.style.display = $("mtlOn").checked ? "" : "none";
  saveCellModelConfig(selectedCellIdx);
  simulate();
});

["concSeries", "tBase", "tAssoc", "tDissoc"]
  .forEach(id => on(id, "input", simulate));

on("genDil", "click", genDilution);

on("regionCount", "input", () => { rebuildCellTypes(); simulate(); });

on("regionSelect", "change", () => {
  saveCellModelConfig(selectedCellIdx);
  const select = $("regionSelect");
  selectedCellIdx = select ? select.selectedIndex : 0;
  loadCellModelConfig(selectedCellIdx);
  refreshFrequencyField();
});

on("frequncy", "input", () => {
  const field = $("frequncy");
  if (!field) return;
  const usedByOthers = cellFrequencies.reduce(
    (sum, f, i) => i === selectedCellIdx ? sum : sum + f, 0
  );
  const max = Math.max(0, 100 - usedByOthers);
  let v = +field.value;
  if (!Number.isFinite(v)) v = 0;
  if (v < 0) v = 0;
  if (v > max) v = max;
  field.value = v;
  field.max = max;
  cellFrequencies[selectedCellIdx] = v;
  state.cellFrequencies = cellFrequencies;
  notifyDataUpdated();
});

setModelVisibility();
rebuildCellTypes();