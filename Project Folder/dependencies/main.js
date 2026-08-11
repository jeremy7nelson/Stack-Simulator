'use strict';

export const state = { parsed: null, lastData: null, capacityField: null };

export const $  = id => document.getElementById(id);
export const on = (id, evt, fn) => { const el = $(id); if (el) el.addEventListener(evt, fn); };

export const setStatus = msg => { const el = $("status"); if (el) el.textContent = msg; };

export const DATA_UPDATED_EVENT = 'data-updated';
export const notifyDataUpdated = () => document.dispatchEvent(new CustomEvent(DATA_UPDATED_EVENT));
export const onDataUpdated = fn => document.addEventListener(DATA_UPDATED_EVENT, fn);

export const IMG_W = 480, IMG_H = 640, MAX16 = 65535;

export const REGION_X = IMG_W / 2, REGION_Y = IMG_H / 2, REGION_R = 140;

export const fmtConc = Cnm => (Cnm >= 1 ? Cnm : Cnm.toPrecision(3)) + " nM";

export const STACK_LEAD_BASELINE_SEC = 5;