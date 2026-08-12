// Very basic version of a different way of generating two-dimentionsal data.
// This prototype does not incorporate user-interface - all parameters hard-coded
// Advantages of this method:
// Integrates one time per model/set of parameters.
// Adding any number of randomly generated circles / overlapping circles is fast and computationally inexpensive
// Adding edge dominance is a simple extension
// Incorporating drift in kinetic parameters is also reasonable, and can be done inexpensively

// Outline of what this does, assuming one dynamic present:
// Step 1: Integrate one time (for whole time series, incorporating all concentrations) at RMax = 1. Set aside result Y
// Step 2: Uses seeded PRNG to distribute random circular regions, with random radii and capacity,
// which then forms the scalar capacity field
// Step 3: The final result is the outer product of the scalar field and Y, result is a flat vector of length m*n*t which corresponds to the
//  mxnxT array. Overlap is resolved by a POWER DIAGRAM (Laguerre tessellation) rather than MAX,
//  with edge dominance measured only from the outer boundary (d = r).
// Step 4: Generation and addition of noise and drift, using the same PRNG procedure with a different seed, according to specifications

// Next steps
// 1. Implement edge dominance - Complete (outer-boundary-only form, under the power diagram)
// 2. Implement all models EXCEPT Bivalent analyte
// 3. Implement Bivalent Analyte
// 4. Kinetic Heterogeneity
// 5. Model Heterogeneity

// QOL Steps
// It would be good to make dt declared in one place alone

// ═══════════════════════════════════════════════════════════════
// Step 1: Integrate once at Rmax = 1
// ═══════════════════════════════════════════════════════════════

const model = "langmuir";
const ka = 1e6, kd = 1e-3, Rmax = 1;
const str = "200, 100, 50, 25, 12.5, 6.25";
const tBase = 30, tAssoc = 120, tDiss = 300;

const vadd = (a,b)=>a.map((v,i)=>v+b[i]);
const vscale = (a,s)=>a.map(v=>v*s);
const vsum  = a=>a.reduce((x,y)=>x+y,0);

function simRK4(grid, deriv, y0, Cfun){
    const out=[]; let y=y0.slice(); out.push(vsum(y));
    for (let i=1;i<grid.length;i++){
        const t0=grid[i-1], h=grid[i]-grid[i-1];
        const C0=Cfun(t0), Cm=Cfun(t0+h/2), C1=Cfun(t0+h);
        const k1=deriv(y, C0);
        const k2=deriv(vadd(y,vscale(k1,h/2)), Cm);
        const k3=deriv(vadd(y,vscale(k2,h/2)), Cm);
        const k4=deriv(vadd(y,vscale(k3,h )), C1);
        y = vadd(y, vscale(vadd(vadd(k1,vscale(k2,2)), vadd(vscale(k3,2),k4)), h/6));
        out.push(vsum(y));
    }
    return out;
}

function makeDeriv(model, RmaxArg){
    return { size:1, deriv:(y,C)=>{ const R=y[0]; return [ka*C*(RmaxArg-R) - kd*R]; } };
}

function parseConcs(str){
    return str.split(/[\s,;]+/).map(s=>parseFloat(s)).filter(v=>Number.isFinite(v)&&v>0);
}

function simulate(){
    const concsNm = parseConcs(str);
    const cyc = tAssoc + tDiss;
    const N = concsNm.length;
    const total = tBase + N*cyc;
    const dt = 2; // SET TIME RESOLUTION HERE
    const npts = Math.round(total/dt);
    const grid=[]; for(let i=0;i<=npts;i++) grid.push(+(i*dt).toFixed(4));

    const concsM = concsNm.map(c=>c*1e-9);
    const Cfun = t=>{
        if (t < tBase) return 0;
        let k = Math.floor((t - tBase)/cyc);
        if (k >= N) k = N-1;
        return ((t - tBase) - k*cyc) < tAssoc ? concsM[k] : 0;
    };
    const base = makeDeriv(model, Rmax);
    return simRK4(grid, base.deriv, new Array(base.size).fill(0), Cfun);
}
const Y = simulate();
const T = Y.length;
console.log(Y.length, "points; final RU =", Y[Y.length-1]);

// ═══════════════════════════════════════════════════════════════
// Step 2: Implement Seeded Pseudo-randomness
// ═══════════════════════════════════════════════════════════════

// This is an alternative to math.random, but exactly reproducible given the same seed.
// seed can be any non-negative integer between 0 and 2^{32}-1
function mulberry32(seed){
    let a = seed >>> 0;
    return function(){
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
const SURFACE_SEED = 12345;          // change this to regenerate a different surface
const rngSurface = mulberry32(SURFACE_SEED);
// This is a function. call this to get a pseudo-random number /in [0,1).  Every call gets a number, and then the sequence advances to a new pseudo-random number.

// ═══════════════════════════════════════════════════════════════
// Step 3: build the scalar capacity field according to pseudo-random center coordinates, radii, and capacities, until a confluence is achieved.
// ═══════════════════════════════════════════════════════════════

const m = 320, n = 240;              // 320 rows (tall) x 240 cols (wide) demonstration window

// Quantize capacities - not necessary for Langmuir, but will be important for bivalent analyte.

const capacityBins = [0.4, 0.6, 0.8, 1.0];

// Hard-coded parameters governing placement
const targetConfluence = 0.50;       // representing a union of pixels
const rMin = 15, rMax = 25;          // radius range, pixels
const MAX_CIRCLES = 5000;            // a cap in iterations to prevent a hangup

// ─── EDGE DOMINANCE: OUTER BOUNDARY ONLY ───────────────────────
// Ownership is geometric (power diagram); weighting is radial in the OWNER's
// own frame. Because the free arc of a cell sits exactly where d = r for its
// own circle, "distance to the outer boundary" is still r - d and needs no
// distance transform. The two concerns decouple cleanly:
//     ownership  <- min power   (which cell holds this pixel)
//     weighting  <- (d/r)^2     (how bright, in the owner's own frame)
// CONSEQUENCE: a cell fully enclosed by neighbours never reaches d = r inside
// its own territory, so its profile tops out below 1 and it shows no bright
// rim. Only cells facing free space get the bright arc. That is arguably the
// right physics (an enclosed cell has no free edge), but it means edge
// dominance fades as confluence rises.
// Weighting on ALL boundaries (incl. flat internal ones) would instead require
// a Euclidean distance transform - deferred.
const edgeFloor = 0.15;              // center is 15% of rim weight, per cell
//   profile(d) = edgeFloor + (1 - edgeFloor) * (d/r)^2
// ───────────────────────────────────────────────────────────────

// covered[k] flags 1/0 whether a pixel is covered, so that we may know when confluence is attained.
// NOTE: the tessellation only redistributes OWNERSHIP inside the union of discs;
// it does not change the union itself, so confluence is unaffected by this change.
const covered = new Uint8Array(m * n);
let coveredCount = 0;

// s[k] = owning cell's capacity * its edge profile at this pixel.
// Continuous again (not four-valued), because the profile is continuous.
const s = new Float32Array(m * n);

// binField[k] = index of the bin of the owning cell at pixel k.
// This will be important for implementing the bivalent analyte model.
const binField = new Int8Array(m * n).fill(-1);

// bestPower[k] = smallest power seen so far at pixel k.
// POWER DIAGRAM (Laguerre tessellation): the power of a point p w.r.t. a circle
// of center c, radius r is |p - c|^2 - r^2. Ownership goes to the circle of
// LEAST power. The locus of equal power between two circles is their radical
// axis — a straight line through their intersection points — so cells meet along
// flat internal boundaries and keep circular arcs only where they face free space.
// min() is commutative, so this stays order-independent and seed-reproducible.
const bestPower = new Float32Array(m * n).fill(Infinity);

const circles = [];                  // record of what was placed

// this function takes the pseudo-random stream and sets the centers, radius, and capacity, then stamps that onto the field.
function placeOneCircle(){
    // draw geometry + capacity from the seeded stream
    const ci  = rngSurface() * m;                          // center row (real-valued)
    const cj  = rngSurface() * n;                          // center col
    const r   = rMin + rngSurface() * (rMax - rMin);       // radius
    const cap = capacityBins[(rngSurface() * capacityBins.length) | 0];
    const binIdx = capacityBins.indexOf(cap);

    const r2 = r * r;
    // bounding box, clamped to grid — only visit pixels near this circle
    const iLo = Math.max(0, Math.floor(ci - r)), iHi = Math.min(m-1, Math.ceil(ci + r));
    const jLo = Math.max(0, Math.floor(cj - r)), jHi = Math.min(n-1, Math.ceil(cj + r));

    for (let i = iLo; i <= iHi; i++){
        for (let j = jLo; j <= jHi; j++){
            const di = i - ci, dj = j - cj;
            const d2 = di*di + dj*dj;
            if (d2 <= r2){                         // clip to this disc: union stays the union
                const k = i*n + j;
                if (covered[k] === 0){ covered[k] = 1; coveredCount++; }

                // MIN POWER replaces the former MAX-on-weighted-value composition.
                // Because we only consider pixels with d2 <= r2, every candidate has
                // power <= 0, so the winner always genuinely contains the pixel.
                const power = d2 - r2;
                if (power < bestPower[k]){
                    bestPower[k]  = power;
                    // edge profile, radial in THIS (the owning) circle's own frame:
                    // floor at its center (d=0) rising to 1.0 at its rim (d=r).
                    const profile = edgeFloor + (1 - edgeFloor) * (d2 / r2);
                    s[k]          = cap * profile;     // capacity sets brightness; power set territory
                    binField[k]   = binIdx;            // owning cell sets the bin too
                }
            }
        }
    }
    circles.push({ ci, cj, r, capacity: cap, binIdx });
}

// continue to place until either union coverage hits the target or iter > MAX_CIRCLES
let iter = 0;
while (coveredCount / (m*n) < targetConfluence && iter < MAX_CIRCLES){
    placeOneCircle();
    iter++;
}
const achievedConfluence = coveredCount / (m*n);
console.log(`placed ${circles.length} circles; confluence target ${targetConfluence}, achieved ${achievedConfluence.toFixed(4)}${iter>=MAX_CIRCLES ? "  (HIT SAFETY CAP)" : ""}`);


// Step 3: outer product  s ⊗ Y  → flat m*n*T stack
// In this step we take the single integration and broadcast it across the scalar array
// Note that this is not valid for the bivalent analyte model.

const stack = new Float32Array(s.length * T);
for (let k = 0; k < s.length; k++){
    const sk = s[k];
    for (let t = 0; t < T; t++) stack[k*T + t] = sk * Y[t];
}

// ═══════════════════════════════════════════════════════════════
// Step 4: noise + drift
// ═══════════════════════════════════════════════════════════════
// This now incorporates its own seed for reproducible noise

// Hard-Coded Noise Parameters:
const D          = 3.5 * Rmax;   // TOTAL accumulated drift
const tau        = 500;          // settling time, s. Larger = more gradual.
const sigmaOU    = 0.02;         // OU step size — the hard-to-subtract wander
const thetaOU    = 0.005;        // OU mean-reversion rate (1/s). Small = long correlation
const decayOU    = true;         // scale OU jitter by the same exp(-t/tau) envelope
const sigmaPixel = 0.10 * Rmax;  // per-pixel noise: 5-30% of Rmax range

const NOISE_SEED = 67890;        // As above, may use any integer in the given range.
const rngNoise = mulberry32(NOISE_SEED);

// Box-Muller transform
function gauss(){
    let u = 0, v = 0;
    while (u === 0) u = rngNoise();
    while (v === 0) v = rngNoise();
    return Math.sqrt(-2*Math.log(u)) * Math.cos(2*Math.PI*v);
}


// drift vector: one value per frame, shared by all pixels
const dtGrid = 2; // THIS SHOULD MATCH TIME RESOLUTION
const driftCommon = new Float32Array(T);
{ let w = 0;    // OU state — persists across t for temporal correlation
    for (let t = 0; t < T; t++){
        const envelope = decayOU ? Math.exp(-t*dtGrid/tau) : 1;
        // OU update: restoring pull toward 0, plus a random kick
        w += -thetaOU*w*dtGrid + sigmaOU*envelope*Math.sqrt(dtGrid)*gauss();
        driftCommon[t] = D*(1 - Math.exp(-t*dtGrid/tau)) + w;
    }
}

// --- apply drift + noise to the stack ---
for (let k = 0; k < s.length; k++) {
    const sk = s[k];
    for (let t = 0; t < T; t++) {
        stack[k*T + t] = sk * Y[t]              // signal
        + driftCommon[t]          // SAME for every k — common-mode
        + sigmaPixel * gauss();   // fresh draw for every (k,t) — i.i.d.
    }
}
