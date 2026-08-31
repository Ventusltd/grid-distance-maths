/**
 * Adversarial verification of the nearest-feature search.
 *
 * WHY THIS FILE EXISTS
 * --------------------
 * verify.mjs already asserted "spatial index finds the same nearest as a
 * brute-force sweep". It passed. The search was still wrong in about one query
 * in nine, because that assertion used a single hand-made fixture whose layout
 * happened not to trip the defect.
 *
 * A fixture proves the code runs. Only a randomised sweep against exhaustive
 * truth proves the code is right, so this generates thousands of layouts across
 * the whole UK latitude band and fails on ANY disagreement.
 *
 * The two defects it was written to catch, both since fixed:
 *
 *   1. The stop test converted cells to km with ky, the latitude scale, on both
 *      axes. A cell is `cell` DEGREES on both axes, but a degree of longitude is
 *      shorter than a degree of latitude and narrows towards the pole -- 0.588
 *      of it at 54N, 0.500 at 60N. Using ky permitted stopping up to twice as
 *      early as the swept box justified. Measured before the fix: 657 of 6,000
 *      queries wrong, worst case a circuit 35.79 km away reported as 65.42 km.
 *
 *   2. addSegment registered only the two endpoint cells, so a segment crossing
 *      a cell without ending in it was invisible from that cell.
 *
 * Both failures overstate distance, never understate it, which is what makes
 * them dangerous: a project reads as further from the network than it is, and
 * nothing about the output looks wrong.
 *
 *   node test/verify_nearest.mjs
 */

import {
  SpatialIndex, distanceToSegmentKm, localScaleKm,
} from "../src/geodesy.mjs";

let passed = 0;
const failures = [];
function check(label, condition, detail = "") {
  if (condition) { passed += 1; console.log(`  [PASS] ${label}${detail ? `  ${detail}` : ""}`); }
  else { failures.push(`${label}${detail ? ` -- ${detail}` : ""}`); console.log(`  [FAIL] ${label}  ${detail}`); }
}

// A deterministic generator, so a failure is reproducible rather than a story
// about a run nobody can repeat.
let rng = 20260831;
const rand = () => { rng = (rng * 1103515245 + 12345) & 0x7fffffff; return rng / 0x7fffffff; };

const brute = (lon, lat, features) => {
  let best = null;
  features.forEach((f, id) => {
    const km = distanceToSegmentKm(lon, lat, f[0], f[1], f[2], f[3]).km;
    if (best === null || km < best.km) best = { id, km };
  });
  return best;
};

const viaIndex = (lon, lat, features, cell) => {
  const index = new SpatialIndex(cell);
  features.forEach((f, id) => index.addSegment(id, f[0], f[1], f[2], f[3]));
  return index.nearest(lon, lat, (id) => {
    const f = features[id];
    return distanceToSegmentKm(lon, lat, f[0], f[1], f[2], f[3]).km;
  });
};

function sweep(label, { trials, cell, spread, segment, latLo, latHi }) {
  let misses = 0;
  let worst = { over: 0, truth: 0, got: 0, lat: 0 };
  for (let t = 0; t < trials; t += 1) {
    const lat = latLo + rand() * (latHi - latLo);
    const lon = -7 + rand() * 9;
    const n = 2 + Math.floor(rand() * 5);
    const features = [];
    for (let i = 0; i < n; i += 1) {
      const flon = lon + (rand() - 0.5) * spread;
      const flat = lat + (rand() - 0.5) * spread;
      features.push([flon, flat, flon + (rand() - 0.5) * segment, flat + (rand() - 0.5) * segment]);
    }
    const truth = brute(lon, lat, features);
    const got = viaIndex(lon, lat, features, cell);
    if (!truth || !got) continue;
    if (Math.abs(truth.km - got.km) > 1e-9) {
      misses += 1;
      if (got.km - truth.km > worst.over) worst = { over: got.km - truth.km, truth: truth.km, got: got.km, lat };
    }
  }
  check(label, misses === 0,
    misses === 0 ? `${trials} layouts, exhaustive agreement`
      : `${misses}/${trials} wrong, worst ${worst.truth.toFixed(2)} km reported as ${worst.got.toFixed(2)} km at ${worst.lat.toFixed(1)}N`);
}

console.log("\nrandomised sweeps against exhaustive truth\n");
sweep("short segments, default 0.1 deg cells", { trials: 6000, cell: 0.1, spread: 2.0, segment: 0.02, latLo: 50, latHi: 61 });
sweep("far north, where a degree of longitude is narrowest", { trials: 3000, cell: 0.1, spread: 2.0, segment: 0.02, latLo: 58, latHi: 61 });
sweep("coarse 0.5 deg cells", { trials: 3000, cell: 0.5, spread: 4.0, segment: 0.05, latLo: 50, latHi: 61 });
sweep("fine 0.02 deg cells", { trials: 3000, cell: 0.02, spread: 0.5, segment: 0.01, latLo: 50, latHi: 61 });
sweep("long segments that cross whole cells", { trials: 3000, cell: 0.1, spread: 2.0, segment: 0.9, latLo: 50, latHi: 61 });
sweep("widely scattered features, sparse buckets", { trials: 2000, cell: 0.1, spread: 8.0, segment: 0.02, latLo: 50, latHi: 61 });

console.log("\nthe specific defects, as direct assertions\n");

// 1. A segment that crosses a cell without ending in it must still be found.
{
  const cell = 0.1;
  const index = new SpatialIndex(cell);
  // Spans four cells of longitude; neither end is in the middle two.
  index.addSegment(0, -1.35, 54.0, -0.95, 54.0);
  const hit = index.nearest(-1.15, 54.005, (id) =>
    distanceToSegmentKm(-1.15, 54.005, -1.35, 54.0, -0.95, 54.0).km);
  const truth = distanceToSegmentKm(-1.15, 54.005, -1.35, 54.0, -0.95, 54.0).km;
  check("a segment crossing a cell is findable from inside that cell",
    hit !== null && Math.abs(hit.km - truth) < 1e-9,
    hit === null ? "not found at all" : `${hit.km.toFixed(4)} km`);
}

// 2. The stop test must never claim more clearance than the box provides.
{
  const cell = 0.1;
  const index = new SpatialIndex(cell);
  for (const lat of [50, 54, 58, 60, 61]) {
    for (const lonOffset of [0.001, 0.037, 0.099]) {
      const lon = -1 + lonOffset;
      const { ky } = localScaleKm(lat);
      for (const ring of [1, 2, 5]) {
        const claimed = index.sweptClearanceKm(lon, lat, ring);
        // Re-derive the box edges independently and take the true distance to
        // the nearest one, on the narrower scale for the east-west pair.
        const ci = Math.floor(lat / cell);
        const cj = Math.floor(lon / cell);
        const latLo = (ci - ring) * cell; const latHi = (ci + ring + 1) * cell;
        const lonLo = (cj - ring) * cell; const lonHi = (cj + ring + 1) * cell;
        const kx = localScaleKm(Math.min(Math.max(Math.abs(latLo), Math.abs(latHi)), 89.9)).kx;
        const actual = Math.min((lat - latLo) * ky, (latHi - lat) * ky,
                                (lon - lonLo) * kx, (lonHi - lon) * kx);
        check(`clearance at ${lat}N ring ${ring} lon+${lonOffset} does not exceed the swept box`,
          claimed <= actual + 1e-9, `claims ${claimed.toFixed(3)} km, box gives ${actual.toFixed(3)} km`);
      }
    }
  }
  check("clearance is zero before any ring is complete", index.sweptClearanceKm(-1, 54, 0) === 0);
}

// 3. Every miss this class of defect produces overstates. Prove the direction.
{
  const cell = 0.1;
  let understated = 0;
  for (let t = 0; t < 2000; t += 1) {
    const lat = 50 + rand() * 11;
    const lon = -7 + rand() * 9;
    const features = [];
    for (let i = 0; i < 5; i += 1) {
      const flon = lon + (rand() - 0.5) * 3;
      const flat = lat + (rand() - 0.5) * 3;
      features.push([flon, flat, flon + 0.01, flat + 0.01]);
    }
    const truth = brute(lon, lat, features);
    const got = viaIndex(lon, lat, features, cell);
    if (got && truth && got.km < truth.km - 1e-9) understated += 1;
  }
  check("the index never returns a distance shorter than the true nearest",
    understated === 0, `${understated} understatements`);
}

console.log(`\n${passed}/${passed + failures.length} checks passed`);
if (failures.length) {
  console.error("\nFAILURES");
  for (const f of failures) console.error("  " + f);
  process.exit(1);
}
