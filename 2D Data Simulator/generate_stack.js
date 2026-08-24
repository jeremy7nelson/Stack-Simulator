// =============================================================================
//  SPRm 2-D DATA SIMULATOR
// =============================================================================
//  Generates a synthetic SPR-microscopy image stack: an (m x n) spatial window
//  observed over a discretised time vector T, for a serial (no-regeneration)
//  injection series.
//
//  Assembly scheme:
//
//    Stack[k][t] =   s[k]        * Yhat_bin[k][t]     (specific binding)
//                  + nsaRate[k]  * Yns[t]             (non-specific adsorption)
//                  + gDrift[k]   * driftCommon[t]     (baseline drift)
//                  + sigma       * hash(i, j, t)      (per-pixel jitter)
//
//  The two spatial weights share one purely geometric edge-proximity field:
//
//    nsaRate[k] = a_n + b_n * w[k]
//    gDrift[k]  = a_d + b_d * w[k]
//    w[k]       = exp(-d/lambdaIn)  inside a cell
//                 exp( d/lambdaOut) outside it,  d = signed distance to the
//                                                cell/medium interface
//
//  w depends only on where the cell/medium interface is — NOT on the owning
//  cell's capacity bin. An earlier revision derived this weight from |grad s|,
//  which made edge drift scale ~2.5x with receptor expression: a high-capacity
//  cell drifted harder at its rim purely because it was brighter. Drift is an
//  optical and mechanical phenomenon and must not inherit that dependence.
//
//  Phase 0 : integrate.  All ODE work; produces a few short time vectors.
//  Phase 1 : build per-pixel weights.  No time axis.
//  Phase 2 : assemble the stack, frame-major.
//
//  ---------------------------------------------------------------------------
//  OBJECTS AVAILABLE FOR PLOTTING  (for whoever builds the interface)
//  ---------------------------------------------------------------------------
//  Four objects are exported specifically so the components can be inspected
//  separately from the assembled stack. All four are plain Float32Arrays and are
//  final once Phase 1 completes.
//
//  TWO TIME SERIES — plot against `grid`, which holds the frame times in
//  SECONDS (length T). Both are BEFORE any spatial weighting.
//
//    Yns          length T.  DIMENSIONLESS fractional occupancy, 0 .. ~0.15 at
//                 the default rate constants, because it is integrated at
//                 Rns = 1. To read it in RU, multiply by a pixel's nsaRate[k]:
//                 the background level is a_n * Yns[t] and the strongest cell
//                 margin is (a_n + b_n) * Yns[t]. Plotting those two bounds
//                 together brackets the whole NSA contribution.
//
//    driftCommon  length T.  Already in RU — no rescaling needed. This is the
//                 shared drift time course; pixel k receives gDrift[k] * this,
//                 with gDrift ranging from a_d (open background) to
//                 a_d + b_d (cell margin).
//
//  TWO SPATIAL FIELDS — both length m*n, row-major, index k = i*n + j, so
//  i = row = Math.floor(k/n) and j = column = k % n. Render directly as an
//  (m x n) image.
//
//    s            DIMENSIONLESS, in [0, 1]. The capacity field: 0 on background,
//                 rising to 1 only at the rim of a full-capacity (bin 1.0) cell.
//                 It is capacity x edge profile, so it carries BOTH how much
//                 receptor a cell has and where on that cell you are. The RU
//                 scale lives in Yhat, not here — a pixel's specific-binding
//                 contribution is s[k] * Yhat[binField[k]][t], and Yhat peaks
//                 near RmaxD. To display s in RU, multiply by RmaxD.
//
//    w            DIMENSIONLESS, in (0, 1]. The edge-proximity weight driving
//                 BOTH nsaRate and gDrift. Peaks just under 1 at the cell/medium
//                 interface and decays inward with lambdaIn, outward with
//                 lambdaOut. Note w is deliberately INDEPENDENT of capacity, so
//                 unlike s it looks identical around every cell — comparing the
//                 two images side by side shows exactly that.
//
//  Also exported and useful for overlays: signedDist (px, + inside a cell,
//  - in open medium, 0 at the interface), binField, covered, circles, gDrift,
//  nsaRate, Yhat (the per-bin kinetic curves, RU, length T each).
//  ---------------------------------------------------------------------------
//
//  NO USER INTERFACE.  Every parameter is hard-coded.  Parameters intended to
//  be exposed by an interface are collected in the USER-DEFINED PARAMETERS
//  block below and annotated with the control type they want.
//
//  MEMORY: the whole stack is held in memory as a Float32Array.
//    640 x 480 x 2551 frames x 4 bytes = 3.13 GB.
//  Run with an enlarged heap:
//      node --max-old-space-size=8192 generate_sprm_stack.js
//  A browser tab cannot allocate this; full-resolution runs need Node.
// =============================================================================


// #############################################################################
// ##  USER-DEFINED PARAMETERS                                                ##
// ##  Everything in this block is intended to be surfaced by an interface.   ##
// ##  Each entry is tagged with the control type it wants, and with whether  ##
// ##  it is ALWAYS shown or CONDITIONAL on the selected model.               ##
// #############################################################################

// ---- Model selection --------------------------------------------------------
const model = "langmuir";
// CONTROL: dropdown list.  ALWAYS SHOWN.
// Options: "langmuir" | "hetLigand" | "twoState" | "bivAnalyte"

const useMTL = false;
// CONTROL: checkbox.  ALWAYS SHOWN.
// Mass-transport limitation. Wraps whichever model is selected above.

// ---- Surface capacity -------------------------------------------------------
const RmaxD = 120;
// CONTROL: numeric entry box.  ALWAYS SHOWN.
// Total surface capacity in RU, at capacity = 1.
// NOTE: this is also the reference scale for sigma and D (see NOISE below), so
// that instrument noise does NOT track ligand density. Changing RmaxD alone
// therefore genuinely changes SNR, which is the intended behaviour.

// ---- Kinetic rate constants -------------------------------------------------
// CONTROL for all of these: numeric entry boxes.
// CONDITIONAL: show only the group matching the selected model.

// model === "langmuir"
const ka = 1e6;                // M^-1 s^-1
const kd = 1e-3;               // s^-1

// model === "hetLigand"   (two independent receptor populations)
const hetKa1 = 1e6, hetKd1 = 1e-2;
const hetKa2 = 1e5, hetKd2 = 1e-3;
const Rmax1Frac = 0.5;
// CONTROL: numeric entry box, 0 < value < 1.  CONDITIONAL on "hetLigand".
// Rmax_1 = Rmax * Rmax1Frac ; Rmax_2 = Rmax * (1 - Rmax1Frac).
// Default 0.5 gives Rmax_1 = Rmax_2 = 0.5 * RmaxD = 60 RU, preserving a TOTAL
// capacity of RmaxD. (draft8 exposed two independent Rmax boxes each defaulting
// to 120, i.e. a total of 240; the fractional form here keeps the total fixed.)
// NOTE: the model is invariant under simultaneously swapping
// (ka1,kd1) <-> (ka2,kd2) and Rmax1Frac -> 1 - Rmax1Frac, so site labelling is
// arbitrary unless a convention is imposed (e.g. site 1 = higher affinity).

// model === "twoState"    (conformational change: A+B <-> AB <-> AB*)
const tsKa1 = 1e6, tsKd1 = 1e-2;
const tsKa2 = 5e-3, tsKd2 = 1e-3;   // ka2, kd2 are s^-1 (surface isomerisation)

// model === "bivAnalyte"  (bivalent analyte)
const bivKa1 = 1e6, bivKd1 = 1e-2;
const bivKa2 = 5e-3, bivKd2 = 1e-3; // ka2 is RU^-1 s^-1 (crosslinking)

const ktr = 1e9;
// CONTROL: numeric entry box.  CONDITIONAL on useMTL === true.
// Transport coefficient, RU M^-1 s^-1. The regime is set by ktr/(ka*Rmax), so
// this must be re-considered whenever RmaxD changes. At RmaxD = 120 and
// ka = 1e6: ktr = 1e9 gives a mild ~3% curve distortion; 1e8 gives ~23%.

// ---- Injection schedule -----------------------------------------------------
const concSeries = "200, 100, 50, 25, 12.5, 6.25";
// CONTROL: multi-line text box (comma/space separated).  ALWAYS SHOWN.
// Analyte concentrations in nM, injected in the order given.
// NOTE: this series is DESCENDING and there is NO regeneration between cycles,
// so bound analyte carries over. Drift accrual is front-loaded, which means the
// high-concentration injections are the most corrupted.

const tBase = 30;      // CONTROL: numeric box. ALWAYS SHOWN. Baseline, s.
const tAssoc = 120;    // CONTROL: numeric box. ALWAYS SHOWN. Association, s.
const tDiss = 300;     // CONTROL: numeric box. ALWAYS SHOWN. Dissociation, s.

// ---- Canvas -----------------------------------------------------------------
const m = 640;         // CONTROL: numeric box. ALWAYS SHOWN. Number of ROWS.
const n = 480;         // CONTROL: numeric box. ALWAYS SHOWN. Number of COLUMNS.
// NOTE: 640 rows x 480 columns is deliberate (a portrait frame), per spec.

// ---- Surface morphology -----------------------------------------------------
const targetConfluence = 0.50;   // CONTROL: numeric box or slider, 0-1. ALWAYS.
const rMin = 15, rMax = 25;      // CONTROL: two numeric boxes. ALWAYS. Pixels.
const capacityBins = [0.4, 0.6, 0.8, 1.0];
// CONTROL: editable list of numbers in (0, 1].  ALWAYS SHOWN.
// Discrete capacities assignable to a cell. Discreteness is what makes the
// per-bin integration path finite for MTL and bivalent analyte.

// ---- Non-specific adsorption ------------------------------------------------
const a_n = 80;
// CONTROL: numeric entry box, > 0.  ALWAYS SHOWN.
// Space-independent NSA weight, RU. Yns peaks at ~0.149 fractional occupancy,
// so a_n = 80 puts background NSA at ~11.9 RU, about equal to the jitter sigma.

const b_n = 80;
// CONTROL: numeric entry box, >= 0.  ALWAYS SHOWN.
// Gradient-weighted NSA weight, RU. Set 0 to make NSA spatially uniform.
// b_n = a_n makes cell edges roughly twice background.

const kaNs = 4e3;      // CONTROL: numeric box. ALWAYS SHOWN. M^-1 s^-1.
const kdNs = 1e-4;     // CONTROL: numeric box. ALWAYS SHOWN. s^-1.
// Slow kdNs (half-life ~6900 s >> run length) makes NSA ACCUMULATE across the
// series rather than washing out. Low kaNs keeps peak occupancy ~0.15, i.e. in
// the quasi-linear regime where the "NSA is linear in C" approximation holds.
// The non-saturable ("unlimited sites") model is this same equation in the
// low-occupancy limit; it is not a separate code path.

// ---- Drift ------------------------------------------------------------------
const b_d = 0.15;
// CONTROL: numeric entry box, 0 - 0.3 suggested.  ALWAYS SHOWN.
// Edge-weighted drift. Set 0 for purely common-mode drift (which a cell-free
// reference region then cancels exactly). Non-zero b_d is what makes drift
// survive reference subtraction, concentrated at cluster edges.
// Because w is bounded in (0, 1] with w = 1 exactly ON the boundary, b_d is now
// literally "the peak fractional drift enhancement at the cell margin".

// ---- Edge-proximity weight w[k] ---------------------------------------------
// w[k] = exp(-d/lambdaIn)  inside a cell,   exp(-d/lambdaOut) outside it,
// where d is the distance to the nearest CELL/MEDIUM interface (the boundary of
// the UNION of all discs, not of any individual disc). This weight is purely
// GEOMETRIC: it is independent of a cell's capacity bin, which is the point.
// Instrument drift is an optical/mechanical phenomenon and has no reason to
// scale with how much receptor a given cell happens to express.
const lambdaIn = 8;
// CONTROL: numeric entry box, px, > 0.  ALWAYS SHOWN.
// Falloff INTO the cell. Larger values push the affected zone further under the
// cell body. Half-width of the ring on the cell side is lambdaIn * ln2 px.

const lambdaOut = 2;
// CONTROL: numeric entry box, px, > 0.  ALWAYS SHOWN.
// Falloff OUT into open medium. The asymmetry is deliberate and mechanism-
// dependent: margin-localised cell processes (lamellipodial dynamics, changes in
// the cell-substrate gap) act on the cell side and argue for lambdaIn >
// lambdaOut, whereas optical effects extending into the medium would argue for
// the reverse. Defaults favour the former.

// ---- Drift time course ------------------------------------------------------
// driftCommon(t) = D * (1 - exp(-t/tau)) + w_OU(t)
// This is the SHARED time course. Every pixel receives it, scaled by gDrift[k].

const D = 3.5 * RmaxD;
// CONTROL: numeric entry box, RU.  ALWAYS SHOWN.
// ASYMPTOTIC drift, i.e. the limit approached as t -> infinity. It is NOT the
// amount accrued over the run: at the default tau the series reaches ~99% of D
// by 2550 s. Expressed against RmaxD so instrument drift does not track ligand
// density; the multiplier, not RmaxD, is the thing to vary.

const tau = 500;
// CONTROL: numeric entry box, s, > 0.  ALWAYS SHOWN.
// EXPONENTIAL TIME CONSTANT — the pace of the approach, not its completion.
// One tau closes 63.2% of the remaining distance, so drift is ~95% complete at
// 3*tau and ~99% at 5*tau. Empirically anchored: tau = 500 s puts the drift RATE
// at ~18% of its initial value by injection 3 and ~8% by injection 4, matching
// the reported "largely settled by the third or fourth injection".

const sigmaOU = 0.02 * RmaxD;
// CONTROL: numeric entry box, RU s^-1/2.  ALWAYS SHOWN.
// Kick size of the Ornstein-Uhlenbeck wander superimposed on the smooth
// exponential. This is what makes the drift unfittable by a parametric baseline
// model; without it driftCommon is a clean analytic curve anyone could subtract.
// NOTE the RmaxD scaling: the source value (0.02 at Rmax = 1) would give a
// 0.2 RU wander at RmaxD = 120, i.e. 0.17% of full scale and invisible.

const thetaOU = 0.005;
// CONTROL: numeric entry box, s^-1, > 0.  ALWAYS SHOWN.
// OU mean-reversion rate. Together with sigmaOU it fixes the stationary spread
// of the wander at sigmaOU/sqrt(2*thetaOU) — 24 RU at the defaults, twice the
// per-pixel jitter. Raising thetaOU makes the wander faster and SMALLER; the two
// parameters are not independent in their effect on amplitude.

// ---- Per-pixel noise --------------------------------------------------------
const sigma = 0.10 * RmaxD;
// CONTROL: numeric entry box, RU.  ALWAYS SHOWN.
// Standard deviation of the independent Gaussian jitter added to every pixel at
// every frame. Expressed against RmaxD for the same reason as D: instrument
// noise is a property of the optics, not of how much ligand is on the surface.
// PROVISIONAL — not yet calibrated against instrument noise characterisation.

// ---- Seeds ------------------------------------------------------------------
const SURFACE_SEED = 12345;   // CONTROL: numeric box, integer 0 .. 2^32-1. ALWAYS.
const NOISE_SEED   = 67890;   // CONTROL: numeric box, integer 0 .. 2^32-1. ALWAYS.
// Two independent streams: changing NOISE_SEED leaves the surface bit-identical.


// #############################################################################
// ##  HARD-CODED PARAMETERS  (not intended for the interface)                ##
// #############################################################################

const dt = 1;                    // time resolution, s
const edgeFloor = 0.15;          // cell centre = 15% of rim weight
const MAX_CIRCLES = 5000;        // iteration cap, prevents a hang
const RnsUnit = 1;               // Yns is integrated at unit capacity

const decayOU = true;            // OU kicks share the exp(-t/tau) envelope
// decayOU = true is correct for this system: if the drift settles by injection
// 3-4, as reported, then its unpredictable component should settle too. Set
// false only to model a wander that persists after the baseline has settled.

const a_d = 1;                   // spatially uniform drift weight


// #############################################################################
// ##  UTILITIES                                                              ##
// #############################################################################

// --- seeded PRNG (surface placement only) ---
function mulberry32(seed){
    let a = seed >>> 0;
    return function(){
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// --- coordinate-keyed integer hash -------------------------------------------
// Returns a uniform in (0,1) determined ENTIRELY by (i, j, t, seed). Unlike a
// sequential stream, the value does not depend on traversal order, on m or n,
// or on how many draws preceded it. That is what allows frame-major assembly
// (and, later, streaming export) to reproduce the same field.
function hashU01(i, j, t, seed){
    let x = (Math.imul(i, 0x9E3779B1) ^ Math.imul(j, 0x85EBCA77)
    ^ Math.imul(t, 0xC2B2AE3D) ^ seed) >>> 0;
    x ^= x >>> 16; x = Math.imul(x, 0x7FEB352D) >>> 0;
    x ^= x >>> 15; x = Math.imul(x, 0x846CA68B) >>> 0;
    x ^= x >>> 16;
    return (x >>> 0) / 4294967296;
}

// --- Box-Muller on hashed uniforms -------------------------------------------
// Two decorrelated uniforms are obtained by hashing with two different seeds.
// The (1 - u) guard avoids log(0) without a rejection loop, so consumption is
// constant and the result stays a pure function of the coordinates.
function hashGauss(i, j, t, seed){
    const u = 1 - hashU01(i, j, t, seed);
    const v = hashU01(i, j, t, (seed ^ 0x5BF03635) >>> 0);
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

const vadd   = (a,b)=>a.map((v,i)=>v+b[i]);
const vscale = (a,s)=>a.map(v=>v*s);
const vsum   = a=>a.reduce((x,y)=>x+y,0);

function parseConcs(str){
    return str.split(/[\s,;]+/).map(s=>parseFloat(s))
    .filter(v=>Number.isFinite(v) && v>0);
}

// --- exact Euclidean distance transform (Felzenszwalb & Huttenlocher, 2012) --
// Returns, for every pixel, the SQUARED Euclidean distance to the nearest pixel
// where mask === 1. Exact, not an approximation, and O(P): two passes of a 1-D
// lower-envelope computation, once down the columns and once along the rows.
// Cost is ~6 ms at 320x240 and ~44 ms at 640x480 — a one-off Phase 1 expense.
const EDT_INF = 1e20;

function edt1d(f, N, d, v, z){
    let k = 0;
    v[0] = 0; z[0] = -Infinity; z[1] = Infinity;
    for (let q = 1; q < N; q++){
        let s = ((f[q] + q*q) - (f[v[k]] + v[k]*v[k])) / (2*q - 2*v[k]);
        while (s <= z[k]){
            k--;
            s = ((f[q] + q*q) - (f[v[k]] + v[k]*v[k])) / (2*q - 2*v[k]);
        }
        k++; v[k] = q; z[k] = s; z[k+1] = Infinity;
    }
    k = 0;
    for (let q = 0; q < N; q++){
        while (z[k+1] < q) k++;
        d[q] = (q - v[k])*(q - v[k]) + f[v[k]];
    }
}

function edt2d(mask, m, n){
    const L = Math.max(m, n);
    const f = new Float64Array(L), d = new Float64Array(L);
    const v = new Int32Array(L + 1), z = new Float64Array(L + 1);
    const out = new Float64Array(m * n);
    for (let k = 0; k < m*n; k++) out[k] = mask[k] ? 0 : EDT_INF;
    for (let j = 0; j < n; j++){                    // columns
        for (let i = 0; i < m; i++) f[i] = out[i*n + j];
        edt1d(f, m, d, v, z);
        for (let i = 0; i < m; i++) out[i*n + j] = d[i];
    }
    for (let i = 0; i < m; i++){                    // rows
        for (let j = 0; j < n; j++) f[j] = out[i*n + j];
        edt1d(f, n, d, v, z);
        for (let j = 0; j < n; j++) out[i*n + j] = d[j];
    }
    return out;
}


// #############################################################################
// ##  PHASE 0 : INTEGRATE                                                    ##
// ##  All ODE work. Output is a handful of short time vectors.               ##
// #############################################################################

const concsNm = parseConcs(concSeries);
const cyc     = tAssoc + tDiss;
const N       = concsNm.length;
const total   = tBase + N * cyc;
const npts    = Math.round(total / dt);
const grid    = new Array(npts + 1);
for (let i = 0; i <= npts; i++) grid[i] = +(i * dt).toFixed(4);
const T = grid.length;

const concsM = concsNm.map(c => c * 1e-9);

// Serial injection schedule. R is never reset between cycles (no regeneration).
function Cfun(t){
    if (t < tBase) return 0;
    let k = Math.floor((t - tBase) / cyc);
    if (k >= N) k = N - 1;
    return ((t - tBase) - k * cyc) < tAssoc ? concsM[k] : 0;
}

// --- generic fixed-step RK4 ---------------------------------------------------
// NOTE: the concentration schedule is discontinuous at every injection edge, and
// an RK4 step straddling such an edge evaluates a stage on the far side of the
// jump. This degrades the method to FIRST order at those instants: at dt = 1 s
// the response at the first association onset reads ~3% of Rmax where the exact
// value is 0. The transient decays within the association phase, but a fresh one
// appears at every onset. The exact fix is to restart the integrator at each
// event time; deferred, but it should not be mistaken for mechanism when reading
// residuals at injection fronts.
function simRK4(grid, deriv, y0, Cfun){
    const out = []; let y = y0.slice(); out.push(vsum(y));
    for (let i = 1; i < grid.length; i++){
        const t0 = grid[i-1], h = grid[i] - grid[i-1];
        const C0 = Cfun(t0), Cm = Cfun(t0 + h/2), C1 = Cfun(t0 + h);
        const k1 = deriv(y, C0);
        const k2 = deriv(vadd(y, vscale(k1, h/2)), Cm);
        const k3 = deriv(vadd(y, vscale(k2, h/2)), Cm);
        const k4 = deriv(vadd(y, vscale(k3, h  )), C1);
        y = vadd(y, vscale(vadd(vadd(k1, vscale(k2,2)), vadd(vscale(k3,2), k4)), h/6));
        out.push(vsum(y));
    }
    return out;
}

// --- model derivatives -------------------------------------------------------
// Each returns { size, deriv(y,C), fluxCoef(y) }.
// fluxCoef expresses the rate at which analyte is drawn from SOLUTION as
// J = a*C + b. Used only by the transport wrapper. Ported from draft8.
function makeDeriv(model, Rtot){
    if (model === "langmuir"){
        return { size: 1,
            deriv: (y,C) => { const R = y[0]; return [ka*C*(Rtot-R) - kd*R]; },
            fluxCoef: (y) => ({ a: ka*(Rtot-y[0]), b: -kd*y[0] }) };
    }
    if (model === "hetLigand"){
        const Rm1 = Rtot * Rmax1Frac, Rm2 = Rtot * (1 - Rmax1Frac);
        return { size: 2,
            deriv: (y,C) => { const R1 = y[0], R2 = y[1];
                return [ hetKa1*C*(Rm1-R1) - hetKd1*R1,
                hetKa2*C*(Rm2-R2) - hetKd2*R2 ]; },
                // Both sites draw from the SAME depleted pool, so under MTL they are
                // coupled and must be integrated jointly. Summing two independent
                // single-site integrations instead introduces a ~2-3% error.
                fluxCoef: (y) => ({ a: hetKa1*(Rm1-y[0]) + hetKa2*(Rm2-y[1]),
                    b: -(hetKd1*y[0] + hetKd2*y[1]) }) };
    }
    if (model === "bivAnalyte"){
        return { size: 2,
            deriv: (y,C) => { const R1 = y[0], R2 = y[1], free = Rtot - R1 - 2*R2;
                return [ 2*bivKa1*C*free - bivKd1*R1 - bivKa2*R1*free + 2*bivKd2*R2,
                2*bivKa2*R1*free - 2*bivKd2*R2 ]; },
                // Only the first step draws from solution; crosslinking does not.
                fluxCoef: (y) => { const free = Rtot - y[0] - 2*y[1];
                    return { a: 2*bivKa1*free, b: -bivKd1*y[0] }; } };
    }
    // two-state conformational change
    return { size: 2,
        deriv: (y,C) => { const AB = y[0], ABs = y[1], free = Rtot - AB - ABs;
            return [ tsKa1*C*free - tsKd1*AB - tsKa2*AB + tsKd2*ABs,
            tsKa2*AB - tsKd2*ABs ]; },
            fluxCoef: (y) => { const free = Rtot - y[0] - y[1];
                return { a: tsKa1*free, b: -tsKd1*y[0] }; } };
}

// --- mass-transport modifier (wraps ANY model) -------------------------------
// Quasi-steady two-compartment: kt*(C - Cs) = J(Cs) = a*Cs + b, so
//   Cs = (kt*C - b) / (kt + a).
function withTransport(base, kt){
    return { size: base.size,
        deriv: (y,C) => {
            const { a, b } = base.fluxCoef(y);
            let Cs = (kt*C - b) / (kt + a);
            if (!isFinite(Cs) || Cs < 0) Cs = 0;
            return base.deriv(y, Cs);
        } };
}

// --- does this configuration need one curve per capacity bin? ----------------
// A model factors (one shared curve suffices) when its solution is LINEAR in
// Rmax. Langmuir, heterogeneous ligand and conformational change all are:
// Rmax enters only the forcing term, never the operator.
//   - Bivalent analyte does NOT: the crosslinking term ka2*R1*free is bilinear
//     in surface species, so the avidity group ka2*Rmax changes with capacity.
//   - Mass transport does NOT, for ANY model: Rmax appears inside the Cs
//     denominator alongside R.
const perBin = useMTL || (model === "bivAnalyte");

// Integrate the bin bank. Each bin b is integrated at Rmax*cap_b and then
// NORMALISED by the dimensionless fraction cap_b (NOT by Rmax*cap_b).
//
//   Yhat_b = Y_b / cap_b       so that      Stack = s[k] * Yhat_bin[k]
//
// is correct, because s[k] already carries capacity (s = cap * profile). Without
// the division capacity would be applied twice and every bin below 1.0 would
// come out low by exactly its own capacity factor. Note that dividing by
// Rmax*cap_b instead would discard the RU scale entirely.
//
// When the model factors, all bins collapse to the same vector automatically --
// no branch is needed, the redundancy is simply not exploited.
function integrateBin(capFraction){
    const Rtot = RmaxD * capFraction;
    let engine = makeDeriv(model, Rtot);
    if (useMTL) engine = withTransport(engine, ktr);
    const Y = simRK4(grid, engine.deriv, new Array(engine.size).fill(0), Cfun);
    return Float32Array.from(Y, v => v / capFraction);   // <- normalise by cap_b
}

const Yhat = [];
if (perBin){
    for (let b = 0; b < capacityBins.length; b++) Yhat.push(integrateBin(capacityBins[b]));
} else {
    const shared = integrateBin(1.0);
    for (let b = 0; b < capacityBins.length; b++) Yhat.push(shared);  // same object
}

// --- NSA basis ---------------------------------------------------------------
// Integrated at Rns = 1, so Yns is DIMENSIONLESS fractional occupancy (peaking
// near 0.15 with the default rate constants). It is deliberately NOT normalised
// to unit peak -- a_n and b_n are calibrated against the raw curve.
const Yns = Float32Array.from(
    simRK4(grid, (y,C) => [ kaNs*C*(RnsUnit - y[0]) - kdNs*y[0] ], [0], Cfun)
);

// --- common-mode drift -------------------------------------------------------
// One value per frame, shared by every pixel. This is the only place the noise
// PRNG is used as a sequential stream; it is a length-T vector, independent of
// canvas size, so it stays reproducible.
const rngNoise = mulberry32(NOISE_SEED);
function gaussStream(){
    let u = 0, v = 0;
    while (u === 0) u = rngNoise();
    while (v === 0) v = rngNoise();
    return Math.sqrt(-2*Math.log(u)) * Math.cos(2*Math.PI*v);
}

const driftCommon = new Float32Array(T);
{
    let w = 0;
    for (let t = 0; t < T; t++){
        const envelope = decayOU ? Math.exp(-t*dt/tau) : 1;
        w += -thetaOU*w*dt + sigmaOU*envelope*Math.sqrt(dt)*gaussStream();
        driftCommon[t] = D*(1 - Math.exp(-t*dt/tau)) + w;
    }
}

console.log(`Phase 0 complete: model="${model}", MTL=${useMTL}, ${perBin ? capacityBins.length : 1} kinetic curve(s), T=${T} frames over ${total}s`);


// #############################################################################
// ##  PHASE 1 : BUILD PER-PIXEL WEIGHTS                                      ##
// ##  No time axis anywhere in this phase.                                   ##
// #############################################################################

const P = m * n;

const s        = new Float32Array(P);              // capacity * edge profile
const binField = new Int8Array(P).fill(-1);        // -1 = background sentinel
const covered  = new Uint8Array(P);
const bestPower = new Float32Array(P).fill(Infinity);
const circles  = [];
let coveredCount = 0;

const rngSurface = mulberry32(SURFACE_SEED);

// --- circle placement, overlap resolved by POWER DIAGRAM ---------------------
// Ownership goes to the least power (d^2 - r^2), which is a Laguerre
// tessellation, NOT a nearest-centre Voronoi: the comparison is radius-weighted,
// so a larger disc can win pixels closer to a smaller disc's centre. Internal
// boundaries are therefore straight radical axes rather than circular arcs.
//
// Edge weighting is radial in the OWNER's own frame, measured from the outer
// boundary only. A cell fully enclosed by neighbours never reaches d = r inside
// its own territory, so it shows no bright rim -- edge dominance fades as
// confluence rises.
function placeOneCircle(){
    const ci  = rngSurface() * m;
    const cj  = rngSurface() * n;
    const r   = rMin + rngSurface() * (rMax - rMin);
    const bin = (rngSurface() * capacityBins.length) | 0;   // keep the INDEX
    const cap = capacityBins[bin];                          // derive the value
    const r2  = r * r;

    const iLo = Math.max(0,   Math.floor(ci - r));
    const iHi = Math.min(m-1, Math.ceil (ci + r));
    const jLo = Math.max(0,   Math.floor(cj - r));
    const jHi = Math.min(n-1, Math.ceil (cj + r));

    for (let i = iLo; i <= iHi; i++){
        for (let j = jLo; j <= jHi; j++){
            const di = i - ci, dj = j - cj;
            const d2 = di*di + dj*dj;
            if (d2 <= r2){
                const k = i*n + j;
                if (covered[k] === 0){ covered[k] = 1; coveredCount++; }
                const power = d2 - r2;                      // <= 0 inside the disc
                if (power < bestPower[k]){
                    bestPower[k] = power;
                    const profile = edgeFloor + (1 - edgeFloor) * (d2 / r2);
                    s[k]        = cap * profile;
                    binField[k] = bin;
                }
            }
        }
    }
    circles.push({ ci, cj, r, capacity: cap, binIdx: bin });
}

let iter = 0;
while (coveredCount / P < targetConfluence && iter < MAX_CIRCLES){
    placeOneCircle();
    iter++;
}

// --- signed distance to the cell / medium interface --------------------------
// Two exact distance transforms give a SIGNED distance to the boundary of the
// UNION of all discs:
//   edt2d(covered)  -> for BACKGROUND pixels, distance to the nearest cell pixel
//   edt2d(inverted) -> for CELL pixels, distance to the nearest background pixel
//
// It must be the union boundary, not individual disc outlines. A disc buried
// inside a confluent cluster has an outline that touches no medium at all; the
// phenomenon being modelled is localised to where cells meet open medium, so a
// buried outline must carry no weight.
const inverted = new Uint8Array(P);
for (let k = 0; k < P; k++) inverted[k] = covered[k] ? 0 : 1;

const dToCell = edt2d(covered,  m, n);    // squared distance, background pixels
const dToVoid = edt2d(inverted, m, n);    // squared distance, cell pixels

// HALF-PIXEL CORRECTION. The transforms measure centre-to-centre distance to
// the nearest pixel of the opposite class, so the smallest value either side is
// 1 px, and the true interface — which lies midway between two adjacent pixel
// centres — would otherwise never be reached. Subtracting 0.5 px from the
// magnitude places d = 0 at the interface itself and makes the field symmetric
// about it. No pixel centre lands exactly on the interface, so w approaches 1
// without attaining it; that is correct, not an error.
const signedDist = new Float32Array(P);   // + inside a cell, - in open medium
for (let k = 0; k < P; k++){
    signedDist[k] = covered[k] ?  (Math.sqrt(dToVoid[k]) - 0.5)
    : -(Math.sqrt(dToCell[k]) - 0.5);
}

// --- edge-proximity weight ---------------------------------------------------
// w = exp(-d/lambdaIn) inside, exp(-d/lambdaOut) outside. Both branches give
// exp(0) = 1 at d = 0, so w is CONTINUOUS across the interface with a kink in
// slope rather than a step — there is no seam artifact to normalise away.
// w is bounded in (0, 1] and needs no percentile normaliser: unlike the old
// gradient field it is already dimensionless, already scaled, and by
// construction identical for every capacity bin.
const w = new Float32Array(P);
for (let k = 0; k < P; k++){
    const d = signedDist[k];
    w[k] = (d >= 0) ? Math.exp(-d / lambdaIn)
    : Math.exp( d / lambdaOut);
}

// --- per-pixel weight fields -------------------------------------------------
const nsaRate = new Float32Array(P);
const gDrift  = new Float32Array(P);
for (let k = 0; k < P; k++){
    nsaRate[k] = a_n + b_n * w[k];   // a_n > 0 required: a pure edge weight would
    // leave the cell-free reference with NO
    // background creep, contradicting the
    // phenomenon being modelled.
    gDrift[k]  = a_d + b_d * w[k];   // b_d > 0 is what survives reference
    // subtraction; the a_d part cancels exactly.
}

{
    let wMax = 0, ring = 0;
    for (let k = 0; k < P; k++){ if (w[k] > wMax) wMax = w[k]; if (w[k] > 0.5) ring++; }
    console.log(`Phase 1 complete: ${circles.length} circles, confluence=${(coveredCount/P).toFixed(4)}, `
    + `lambdaIn=${lambdaIn} lambdaOut=${lambdaOut}, peak w=${wMax.toFixed(4)}, `
    + `${(100*ring/P).toFixed(1)}% of frame with w>0.5`);
}


// #############################################################################
// ##  PHASE 2 : ASSEMBLE THE STACK, FRAME-MAJOR                              ##
// #############################################################################
//
// Storage layout is flat: index = k*T + t, with k = i*n + j. All time points of
// one pixel are contiguous, which is the fast layout for extracting a single
// pixel's sensorgram and the slow one for extracting a frame.
//
// The outer loop runs over t and the inner over pixels. Because the jitter is
// hashed on (i, j, t) rather than drawn from a stream, this ordering produces
// exactly the same field a pixel-major loop would -- but it is also the ordering
// a streaming exporter needs, so the loop can later emit and discard each frame
// without changing any arithmetic.

const stack = new Float32Array(P * T);

for (let t = 0; t < T; t++){
    const dc  = driftCommon[t];
    const yns = Yns[t];
    let k = 0;
    for (let i = 0; i < m; i++){
        for (let j = 0; j < n; j++, k++){
            const b = binField[k];
            // Background: binField = -1 and s[k] = 0, so the specific term
            // vanishes arithmetically. The guard exists only to keep the array
            // index in range, not to change the result.
            const specific = (b < 0) ? 0 : s[k] * Yhat[b][t];
            stack[k*T + t] = specific
            + nsaRate[k] * yns
            + gDrift[k]  * dc
            + sigma * hashGauss(i, j, t, NOISE_SEED);
        }
    }
}

console.log(`Phase 2 complete: stack ${m}x${n}x${T} = ${(stack.length*4/1e9).toFixed(2)} GB`);


// #############################################################################
// ##  EXPORTS                                                                ##
// #############################################################################

if (typeof module !== "undefined" && module.exports){
    module.exports = {
        stack, m, n, T, dt, grid,
        s, binField, capacityBins, circles,
        signedDist, w, lambdaIn, lambdaOut, nsaRate, gDrift, covered,
        Yhat, Yns, driftCommon,
        perBin, model, useMTL, RmaxD
    };
}
