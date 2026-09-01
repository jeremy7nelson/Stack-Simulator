import { state, setStatus, $, fmtConc, IMG_W, IMG_H, MAX16, STACK_LEAD_BASELINE_SEC } from './main.js';
import { getMatrix16, decodeFrame, findPeakInjectionFrame, getMatrix16ForBrightness } from './render.js';

const CALIBRATION_RMAX = -0.575;
const BLANK_RMAX = 0;

function u8ToBase64(u8) {
  const CHUNK = 0x8000;
  let binary = "";
  for (let i = 0; i < u8.length; i += CHUNK) {
    binary += String.fromCharCode.apply(null, u8.subarray(i, i + CHUNK));
  }
  return btoa(binary);
}

async function deflateRawCompress(u8) {
  const cs = new CompressionStream('deflate-raw');
  const writer = cs.writable.getWriter();
  writer.write(u8);
  writer.close();
  const chunks = [];
  const reader = cs.readable.getReader();
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  const totalLen = chunks.reduce((n, c) => n + c.length, 0);
  const out = new Uint8Array(totalLen);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

async function deflateRawDecompress(u8) {
  const ds = new DecompressionStream('deflate-raw');
  const writer = ds.writable.getWriter();
  writer.write(u8);
  writer.close();
  const chunks = [];
  const reader = ds.readable.getReader();
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  const totalLen = chunks.reduce((n, c) => n + c.length, 0);
  const out = new Uint8Array(totalLen);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

async function encodeCompressedData(typedArrayOrBuffer) {
  const u8 = typedArrayOrBuffer instanceof ArrayBuffer
      ? new Uint8Array(typedArrayOrBuffer)
      : new Uint8Array(typedArrayOrBuffer.buffer, typedArrayOrBuffer.byteOffset, typedArrayOrBuffer.byteLength);
  const compressed = await deflateRawCompress(u8);
  return u8ToBase64(compressed);
}

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(u8) {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < u8.length; i++) crc = CRC_TABLE[(crc ^ u8[i]) & 0xFF] ^ (crc >>> 8);
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

function dosDateTime(date) {
  const time = ((date.getHours() & 0x1f) << 11) | ((date.getMinutes() & 0x3f) << 5) | ((date.getSeconds() >> 1) & 0x1f);
  const day  = (((date.getFullYear() - 1980) & 0x7f) << 9) | (((date.getMonth() + 1) & 0xf) << 5) | (date.getDate() & 0x1f);
  return { time, day };
}

class ZipWriter {
  constructor() {
    this.chunks = [];
    this.offset = 0;
    this.central = [];
  }

  _push(u8) {
    this.chunks.push(u8);
    this.offset += u8.length;
  }

  async addFile(name, uint8Data, when = new Date()) {
    const nameBytes = new TextEncoder().encode(name);
    const { time, day } = dosDateTime(when);
    const crc = crc32(uint8Data);
    const compressed = await deflateRawCompress(uint8Data);
    const usize = uint8Data.length, csize = compressed.length;
    const localOffset = this.offset;

    const header = new Uint8Array(30 + nameBytes.length);
    const dv = new DataView(header.buffer);
    dv.setUint32(0, 0x04034b50, true);
    dv.setUint16(4, 20, true);
    dv.setUint16(6, 0, true);
    dv.setUint16(8, 8, true);
    dv.setUint16(10, time, true);
    dv.setUint16(12, day, true);
    dv.setUint32(14, crc, true);
    dv.setUint32(18, csize, true);
    dv.setUint32(22, usize, true);
    dv.setUint16(26, nameBytes.length, true);
    dv.setUint16(28, 0, true);
    header.set(nameBytes, 30);

    this._push(header);
    this._push(compressed);
    this.central.push({ name: nameBytes, crc, csize, usize, time, day, offset: localOffset });
  }

  finalize() {
    const centralStart = this.offset;
    for (const e of this.central) {
      const header = new Uint8Array(46 + e.name.length);
      const dv = new DataView(header.buffer);
      dv.setUint32(0, 0x02014b50, true);
      dv.setUint16(4, 20, true);
      dv.setUint16(6, 20, true);
      dv.setUint16(8, 0, true);
      dv.setUint16(10, 8, true);
      dv.setUint16(12, e.time, true);
      dv.setUint16(14, e.day, true);
      dv.setUint32(16, e.crc, true);
      dv.setUint32(20, e.csize, true);
      dv.setUint32(24, e.usize, true);
      dv.setUint16(28, e.name.length, true);
      dv.setUint16(30, 0, true);
      dv.setUint16(32, 0, true);
      dv.setUint16(34, 0, true);
      dv.setUint16(36, 0, true);
      dv.setUint32(38, 0, true);
      dv.setUint32(42, e.offset, true);
      header.set(e.name, 46);
      this._push(header);
    }
    const centralSize = this.offset - centralStart;

    const eocd = new Uint8Array(22);
    const dv = new DataView(eocd.buffer);
    dv.setUint32(0, 0x06054b50, true);
    dv.setUint16(4, 0, true);
    dv.setUint16(6, 0, true);
    dv.setUint16(8, this.central.length, true);
    dv.setUint16(10, this.central.length, true);
    dv.setUint32(12, centralSize, true);
    dv.setUint32(16, centralStart, true);
    dv.setUint16(20, 0, true);
    this._push(eocd);

    const out = new Uint8Array(this.offset);
    let pos = 0;
    for (const c of this.chunks) { out.set(c, pos); pos += c.length; }
    return out;
  }
}

function formatTimestamp(when = new Date()) {
  const pad = n => String(n).padStart(2, '0');
  return `${pad(when.getMonth()+1)}/${pad(when.getDate())}/${when.getFullYear()} ` +
         `${pad(when.getHours())}:${pad(when.getMinutes())}:${pad(when.getSeconds())}`;
}


async function encodeTimeInput() {
  const { nFrames, nSpots } = state.parsed;
  const grid = state.lastData.grid;
  const nTotal = nSpots + 2;
  const f32  = new Float32Array(nFrames * nTotal);

  for (let seqIndex = 0; seqIndex < nTotal; seqIndex++) {
    const timeOffset = stkTimeOffset(seqIndex);
    for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
      f32[seqIndex * nFrames + timeIdx] = timeOffset + grid[timeIdx];
    }
  }

  const dataB64 = await encodeCompressedData(f32);
  return (
    `  <Input>\n` +
    `    <Name>Time</Name>\n` +
    `    <Data>${dataB64}</Data>\n` +
    `  </Input>`
  );
}

async function encodeResponseInput() {
  const { regions, nFrames, nSpots } = state.parsed;
  const tAssoc = +$("tAssoc").value;
  const nTotal = nSpots + 2;

  const calibTrace = stepResponseTrace(nFrames, STACK_LEAD_BASELINE_SEC, tAssoc, CALIBRATION_RMAX);
  const blankTrace  = stepResponseTrace(nFrames, STACK_LEAD_BASELINE_SEC, tAssoc, BLANK_RMAX);

  const parts = [];
  for (const rg of regions) {
    const f32 = new Float32Array(nFrames * nTotal);

    for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
      f32[0 * nFrames + timeIdx] = calibTrace[timeIdx] / 240;
    }
    for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
      f32[1 * nFrames + timeIdx] = blankTrace[timeIdx] / 240;
    }
    for (let concIdx = 0; concIdx < nSpots; concIdx++) {
      const seqIndex = concIdx + 2;
      for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
        const ru = concIdx < rg.traces.length ? rg.traces[concIdx][timeIdx] : 0;
        f32[seqIndex * nFrames + timeIdx] = ru / 240;
      }
    }

    const dataB64 = await encodeCompressedData(f32);
    parts.push(
      `  <Input>\n` +
      `    <Name>Roi${rg.idx}</Name>\n` +
      `    <Data>${dataB64}</Data>\n` +
      `  </Input>`
    );
  }
  return parts.join("\n");
}

function stkTimeOffset(concIdx) {
  const g = state.lastData.grid, n = state.parsed.nFrames;
  const span = g[n - 1] - g[0];
  return concIdx * span;
}

function stkFileName(concIdx) {
  const c = state.parsed.concs[concIdx];
  const tag = (c != null) ? fmtConc(c).replace(/[^0-9A-Za-z.]+/g, '') : `spot${concIdx + 1}`;
  return `spr_stack_${tag}.stk`;
}

function getInjectionWindow(seqIndex = 0) {
  const tAssoc = +$("tAssoc").value;
  const offset = stkTimeOffset(seqIndex);
  const tB = STACK_LEAD_BASELINE_SEC + offset;
  const tD = STACK_LEAD_BASELINE_SEC + tAssoc + offset;
  return { tB, tD };
}

function stepResponseTrace(nFrames, leadSec, tAssoc, Rmax) {
  const trace = new Float64Array(nFrames);
  for (let t = 0; t < nFrames; t++) {
    trace[t] = (t >= leadSec && t < leadSec + tAssoc) ? Rmax : 0;
  }
  return trace;
}

function specialTraceBrightness16(trace, timeIdx) {
  let lo = Infinity, hi = -Infinity;
  for (const v of trace) { if (v < lo) lo = v; if (v > hi) hi = v; }
  const denom = (hi - lo) || 1;
  const norm = (trace[timeIdx] - lo) / denom;
  return Math.round(norm * MAX16);
}

function buildSpecialStkBuffer(kind, baseDate, wallSeqIndex) {
  const FRAME_TYPE_SPR_GRAY16 = 101;
  const nFrames       = state.parsed.nFrames;
  const bytesPerFrame = IMG_W * IMG_H * 2;
  const tAssoc = +$("tAssoc").value;

  const isCalibration = kind === 'calibration';
  const trace = stepResponseTrace(
    nFrames, STACK_LEAD_BASELINE_SEC, tAssoc,
    isCalibration ? CALIBRATION_RMAX : BLANK_RMAX
  );

  const timeOffset   = 0;
  const wallOffset   = stkTimeOffset(wallSeqIndex);
  const startTimeStr = formatTimestamp(new Date(baseDate.getTime() + wallOffset * 1000));

  const enc = new TextEncoder();
  const concStr  = isCalibration ? 'C' : '0';
  const labelStr = 'SPR simulation';
  const descStr  = isCalibration ? 'calibration standard (step)' : 'blank / buffer injection';

  const strBytes = s => enc.encode(s);
  const strSize  = s => 1 + strBytes(s).length;

  const headerSize =
    4 + strSize(startTimeStr) + strSize(concStr) + strSize(labelStr) + strSize(descStr) + 4 * 12;
  const frameHeaderSize = 4 + 4 + 4 + 4;
  const markersSize = 6 * 4;
  const totalSize = headerSize + nFrames * (frameHeaderSize + bytesPerFrame) + markersSize;

  const buf  = new ArrayBuffer(totalSize);
  const view = new DataView(buf);
  const u8   = new Uint8Array(buf);
  let   pos  = 0;
  const LE   = true;

  const writeInt32   = v => { view.setInt32(pos, v, LE);   pos += 4; };
  const writeFloat32 = v => { view.setFloat32(pos, v, LE); pos += 4; };
  const writeString  = s => { const b = strBytes(s); u8[pos++] = b.length; u8.set(b, pos); pos += b.length; };

  writeInt32(2);
  writeString(startTimeStr);
  writeString(concStr);
  writeString(labelStr);
  writeString(descStr);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(25.0);
  writeFloat32(0.0);
  writeFloat32(1.0);
  writeFloat32(1.0);
  writeFloat32(0.0);
  writeFloat32(1.0);
  writeFloat32(1.0);
  writeFloat32(1.0);

  for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
    const timestamp = timeOffset + timeIdx;
    const b16 = specialTraceBrightness16(trace, timeIdx);
    const mat = getMatrix16ForBrightness(b16);

    writeInt32(FRAME_TYPE_SPR_GRAY16);
    writeFloat32(timestamp);
    writeInt32(IMG_W);
    writeInt32(IMG_H);
    for (let i = 0; i < mat.length; i++) { view.setUint16(pos, mat[i], LE); pos += 2; }
  }

  writeInt32(1000);
  writeFloat32(STACK_LEAD_BASELINE_SEC + timeOffset);
  writeInt32(101);
  writeInt32(1000);
  writeFloat32(STACK_LEAD_BASELINE_SEC + tAssoc + timeOffset);
  writeInt32(102);

  return buf;
}

function buildStkBuffer(baseDate = new Date(), concIdx = 0, seqIndex = concIdx, wallSeqIndex = seqIndex) {
  const FRAME_TYPE_SPR_GRAY16 = 101;
  const nFrames       = state.parsed.nFrames;
  const bytesPerFrame = IMG_W * IMG_H * 2;
  const timeOffset    = stkTimeOffset(seqIndex);
  const wallOffset    = stkTimeOffset(wallSeqIndex);

  const startTimeStr  = formatTimestamp(new Date(baseDate.getTime() + wallOffset * 1000));

  const enc = new TextEncoder();
  const c        = state.parsed.concs[concIdx];
  const concStr  = (c != null) ? fmtConc(c) : `spot ${concIdx + 1}`;
  const labelStr = 'SPR simulation';
  const descStr  = `model: ${$('model').value}`;

  const strBytes = s => enc.encode(s);
  const strSize  = s => 1 + strBytes(s).length;

  const headerSize =
    4 + strSize(startTimeStr) + strSize(concStr) + strSize(labelStr) + strSize(descStr) + 4 * 12;
  const frameHeaderSize = 4 + 4 + 4 + 4;
  const markersSize = 6 * 4;
  const totalSize = headerSize + nFrames * (frameHeaderSize + bytesPerFrame) + markersSize;

  const buf  = new ArrayBuffer(totalSize);
  const view = new DataView(buf);
  const u8   = new Uint8Array(buf);
  let   pos  = 0;
  const LE   = true;

  const writeInt32   = v => { view.setInt32(pos, v, LE);   pos += 4; };
  const writeFloat32 = v => { view.setFloat32(pos, v, LE); pos += 4; };
  const writeString  = s => { const b = strBytes(s); u8[pos++] = b.length; u8.set(b, pos); pos += b.length; };

  writeInt32(2);
  writeString(startTimeStr);
  writeString(concStr);
  writeString(labelStr);
  writeString(descStr);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(0.0);
  writeFloat32(25.0);
  writeFloat32(0.0);
  writeFloat32(1.0);
  writeFloat32(1.0);
  writeFloat32(0.0);
  writeFloat32(1.0);
  writeFloat32(1.0);
  writeFloat32(1.0);

  for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
    const f         = concIdx * nFrames + timeIdx;
    const timestamp = timeOffset + state.lastData.grid[timeIdx];
    const mat       = getMatrix16(f);

    writeInt32(FRAME_TYPE_SPR_GRAY16);
    writeFloat32(timestamp);
    writeInt32(IMG_W);
    writeInt32(IMG_H);
    for (let i = 0; i < mat.length; i++) { view.setUint16(pos, mat[i], LE); pos += 2; }
  }

  const { tB, tD } = getInjectionWindow(seqIndex);
  writeInt32(1000);
  writeFloat32(tB);
  writeInt32(101);
  writeInt32(1000);
  writeFloat32(tD);
  writeInt32(102);

  return buf;
}

async function addStkEntries(writer, baseDate = new Date()) {
  await writer.addFile('Stacks/spr_stack_calibration.stk', new Uint8Array(buildSpecialStkBuffer('calibration', baseDate, 0)), baseDate);
  await writer.addFile('Stacks/spr_stack_blank.stk', new Uint8Array(buildSpecialStkBuffer('blank', baseDate, 1)), baseDate);
  for (let c = 0; c < state.parsed.nSpots; c++) {
    await writer.addFile('Stacks/' + stkFileName(c), new Uint8Array(buildStkBuffer(baseDate, c, c, c + 2)), baseDate);
  }
}

async function buildRoiXml(timestamp = formatTimestamp()) {
  const { frame: peakFrame } = findPeakInjectionFrame();

  const mat16   = getMatrix16(peakFrame);
  const grayBuf = new ArrayBuffer(mat16.length * 2);
  const grayDV  = new DataView(grayBuf);
  for (let i = 0; i < mat16.length; i++) grayDV.setUint16(i * 2, mat16[i], true);
  const grayCompressed = await deflateRawCompress(new Uint8Array(grayBuf));
  const grayB64  = u8ToBase64(grayCompressed);
  const sprGrayW = IMG_W, sprGrayH = IMG_H;

  const BF_W = sprGrayW * 2, BF_H = sprGrayH * 2;
  const bfBuf = new Uint8Array(BF_W * BF_H * 3).fill(128);
  const bfCompressed = await deflateRawCompress(bfBuf);
  const bfB64 = u8ToBase64(bfCompressed);

  const winW = Math.round(sprGrayW / 2), winH = Math.round(sprGrayH / 2);
  const winX0 = Math.round((sprGrayW - winW) / 2), winY0 = Math.round((sprGrayH - winH) / 2);
  const sprWindow = `${winX0}, ${winY0}, ${winX0 + winW}, ${winY0 + winH}`;

  return `<?xml version="1.0" encoding="utf-8"?>
          <RoiGroup>
            <Timestamp>${timestamp}</Timestamp>
            <SprWindow>${sprWindow}</SprWindow>
            <Snapshot>
              <Width>${sprGrayW}</Width>
              <Height>${sprGrayH}</Height>
              <Type>SprGray16</Type>
              <Data>${grayB64}</Data>
            </Snapshot>
            <Snapshot>
              <Width>${BF_W}</Width>
              <Height>${BF_H}</Height>
              <Type>BrightFieldRgb24</Type>
              <Data>${bfB64}</Data>
            </Snapshot>
          </RoiGroup>`;
}

async function buildBiXml(timestamp = formatTimestamp()) {
  const timeBlock  = await encodeTimeInput();
  const roiEntries = await encodeResponseInput();

  return `<?xml version="1.0" encoding="utf-8"?>
          <SPRm-Realtime>
            <Version>2.8.2</Version>
            <StartTime>${timestamp}</StartTime>
          ${roiEntries}
          ${timeBlock}
          </SPRm-Realtime>
          `;
}

async function downloadAll() {
  if (!state.parsed || !state.lastData) return;
  setStatus("Building export…");

  const baseDate  = new Date();
  const startTime = formatTimestamp(baseDate);
  const enc = new TextEncoder();

  const roiXml = await buildRoiXml(startTime);
  const biXml  = await buildBiXml(startTime);

  const writer = new ZipWriter();
  await writer.addFile('ROI/spr.roi', enc.encode(roiXml), baseDate);
  await writer.addFile('Data/data.bi', enc.encode(biXml), baseDate);
  await addStkEntries(writer, baseDate);

  const zipBytes = writer.finalize();
  const blob = new Blob([zipBytes], { type: "application/zip" });
  const a = Object.assign(document.createElement("a"), {
    href:     URL.createObjectURL(blob),
    download: "spr_export.zip"
  });
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);

  setStatus("Export ready: spr_export.zip");
}

document.getElementById("downloadAll").addEventListener("click", downloadAll);