/**
 * Grid Distance Maths — verification.
 *
 * Checks the canonical implementation against ellipsoidal truth and against
 * the specific defects that made this repository necessary. Run:
 *
 *   node test/verify.mjs
 *
 * Exits non-zero on any failure. No dependencies.
 */

import {
  R_ATLAS, R_MEAN, R_UK,
  distanceKm, distanceEllipsoidalKm, pathLengthKm,
  destinationPoint, initialBearingDeg,
  distanceToSegmentKm, distanceToLineKm,
  geodesicCircle, circleApproximationError, polygonAreaKm2, areaUnits,
  SpatialIndex, curvatureKm,
} from "../src/geodesy.mjs";

const checks = [];
const ok = (name, pass, detail = "") => checks.push({ name, pass: Boolean(pass), detail });
const near = (a, b, tol) => Math.abs(a - b) <= tol;

/* ---- constants ---------------------------------------------------------- */
ok("R_ATLAS is the WGS84 semi-major axis", R_ATLAS === 6378.137);
ok("R_MEAN is the IUGG mean radius", R_MEAN === 6371.0088);
ok("R_UK is the Gaussian radius at 54N", near(R_UK, 6384.7272, 1e-4));
ok("R_ATLAS and R_MEAN differ by 0.112 percent",
  near(((R_ATLAS - R_MEAN) / R_MEAN) * 100, 0.1119, 1e-4));

/* ---- agreement with the deployed tools ---------------------------------- */
// ventus-corev8engine.js haversine(), transcribed independently.
function engineHaversine(lon1, lat1, lon2, lat2) {
  const R = 6378.137; const r = Math.PI / 180;
  const dLat = (lat2 - lat1) * r; const dLon = (lon2 - lon1) * r;
  const a = Math.sin(dLat / 2) ** 2
    + Math.cos(lat1 * r) * Math.cos(lat2 * r) * Math.sin(dLon / 2) ** 2;
  return R * 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
}
let worstDrift = 0;
for (let i = 0; i < 500; i += 1) {
  const lo1 = -7 + Math.random() * 9; const la1 = 50 + Math.random() * 8.5;
  const lo2 = lo1 + (Math.random() - 0.5); const la2 = la1 + (Math.random() - 0.5);
  worstDrift = Math.max(worstDrift,
    Math.abs(distanceKm(lo1, la1, lo2, la2) - engineHaversine(lo1, la1, lo2, la2)));
}
ok("distanceKm agrees with the GridAtlas engine to the last bit",
  worstDrift < 1e-12, `worst drift ${worstDrift.toExponential(2)} km`);

/* ---- ellipsoidal reference ---------------------------------------------- */
// Land's End to John o' Groats, a published long UK baseline.
const lgj = distanceEllipsoidalKm(-5.7124, 50.0657, -3.0170, 58.6373);
ok("Vincenty reproduces the Land's End to John o' Groats baseline",
  near(lgj, 969.87, 0.05), `${lgj.toFixed(3)} km`);
ok("Vincenty is symmetric",
  near(lgj, distanceEllipsoidalKm(-3.0170, 58.6373, -5.7124, 50.0657), 1e-9));
ok("zero distance is zero", distanceEllipsoidalKm(-1.5, 54, -1.5, 54) === 0);

/* R_UK must beat R_ATLAS, which must beat R_MEAN, over UK pairs. */
let sumUk = 0; let sumAtlas = 0; let sumMean = 0; let n = 0;
for (let i = 0; i < 2000; i += 1) {
  const la1 = 50 + Math.random() * 8.5; const lo1 = -7 + Math.random() * 9;
  const brg = Math.random() * 360; const d = 0.5 + Math.random() * 59.5;
  const [lo2, la2] = destinationPoint(lo1, la1, d, brg);
  if (la2 < 50 || la2 > 58.6 || lo2 < -7.2 || lo2 > 1.9) continue;
  const truth = distanceEllipsoidalKm(lo1, la1, lo2, la2);
  if (!Number.isFinite(truth) || truth <= 0) continue;
  sumUk += Math.abs(distanceKm(lo1, la1, lo2, la2, R_UK) - truth) / truth;
  sumAtlas += Math.abs(distanceKm(lo1, la1, lo2, la2, R_ATLAS) - truth) / truth;
  sumMean += Math.abs(distanceKm(lo1, la1, lo2, la2, R_MEAN) - truth) / truth;
  n += 1;
}
ok("R_UK is more accurate than R_ATLAS over UK pairs",
  sumUk / n < sumAtlas / n, `${(1e6 * sumUk / n).toFixed(0)} vs ${(1e6 * sumAtlas / n).toFixed(0)} ppm`);
ok("R_ATLAS is more accurate than R_MEAN over UK pairs",
  sumAtlas / n < sumMean / n, `${(1e6 * sumAtlas / n).toFixed(0)} vs ${(1e6 * sumMean / n).toFixed(0)} ppm`);

/* ---- bearing and destination round-trip --------------------------------- */
const [dl, dp] = destinationPoint(-1.5, 54, 25, 73.4);
ok("destination then measure returns the distance travelled",
  near(distanceKm(-1.5, 54, dl, dp), 25, 1e-9), `${distanceKm(-1.5, 54, dl, dp).toFixed(9)} km`);
ok("destination then bearing returns the bearing travelled",
  near(initialBearingDeg(-1.5, 54, dl, dp), 73.4, 1e-6));

/* ---- point to segment: the defect this repo exists to prevent ------------ */
// A due-east span at 54N. A point opposite its midpoint is 1 km away from the
// LINE, but much further from either END. Measuring to endpoints overstates.
const A = [-1.5, 54.0];
const B = destinationPoint(A[0], A[1], 20, 90);
const mid = destinationPoint(A[0], A[1], 10, 90);
const off = destinationPoint(mid[0], mid[1], 1, 0);
const seg = distanceToSegmentKm(off[0], off[1], A[0], A[1], B[0], B[1]);
const toNearestEnd = Math.min(
  distanceKm(off[0], off[1], A[0], A[1]),
  distanceKm(off[0], off[1], B[0], B[1]),
);
// A GeoJSON segment is a CHORD, so the measured distance is to the chord, not
// to the great-circle arc through the same endpoints. The arc bulges away from
// the chord by L^2/8R -- 7.8 m over this deliberately long 20 km span. Real
// network segments are ~200 m, where the effect is nanometres, so the tolerance
// here is set to admit the known sagitta rather than to hide it.
ok("distance to a segment is the perpendicular, not the nearest end",
  near(seg.km, 1, 0.02), `${seg.km.toFixed(4)} km, chord sagitta included`);
{
  // The same geometry at a realistic segment length, where chord and arc agree.
  const a2 = [-1.5, 54.0];
  const b2 = destinationPoint(a2[0], a2[1], 0.5, 90);
  const m2 = destinationPoint(a2[0], a2[1], 0.25, 90);
  const o2 = destinationPoint(m2[0], m2[1], 1, 0);
  const s2 = distanceToSegmentKm(o2[0], o2[1], a2[0], a2[1], b2[0], b2[1]);
  ok("on a realistic 500 m segment the perpendicular is exact",
    near(s2.km, 1, 5e-5), `${s2.km.toFixed(7)} km`);
  ok("chord sagitta scales as L squared over 8R",
    near(seg.km - 1, (20 ** 2) / (8 * R_ATLAS), 0.004),
    `measured ${((seg.km - 1) * 1000).toFixed(1)} m, predicted ${((20 ** 2) / (8 * R_ATLAS) * 1000).toFixed(1)} m`);
}
ok("measuring to the nearest end would have overstated by an order of magnitude",
  toNearestEnd > 9, `${toNearestEnd.toFixed(2)} km vs ${seg.km.toFixed(2)} km`);
ok("the foot of the perpendicular lands inside the segment",
  seg.t > 0.45 && seg.t < 0.55, `t=${seg.t.toFixed(4)}`);

// Beyond an endpoint, the projection must clamp rather than run off the line.
const past = destinationPoint(B[0], B[1], 5, 90);
const clamped = distanceToSegmentKm(past[0], past[1], A[0], A[1], B[0], B[1]);
ok("projection clamps at the segment end", clamped.t === 1 && near(clamped.km, 5, 0.01));

const line = [A, mid, B];
ok("polyline search finds the same perpendicular",
  near(distanceToLineKm(off[0], off[1], line).km, 1, 0.01));

/* ---- circles ------------------------------------------------------------ */
const ring = geodesicCircle(-1.5, 54, 10, 128);
ok("geodesic circle closes", ring[0][0] === ring[ring.length - 1][0]);
const radii = ring.slice(0, -1).map((p) => distanceKm(-1.5, 54, p[0], p[1]));
ok("every circle vertex is exactly the radius away",
  radii.every((r) => near(r, 10, 1e-9)), `spread ${(Math.max(...radii) - Math.min(...radii)).toExponential(2)}`);
const e24 = circleApproximationError(24);
const e128 = circleApproximationError(128);
ok("24-vertex ring understates area by 1.14 percent",
  near(e24.areaShortFraction * 100, 1.138, 0.01), `${(e24.areaShortFraction * 100).toFixed(3)}%`);
ok("128-vertex ring understates area by 0.04 percent",
  near(e128.areaShortFraction * 100, 0.040, 0.005));
ok("24-vertex ring is 85.6 m inside a 10 km circle",
  near(e24.radialFraction * 10000, 85.6, 0.5), `${(e24.radialFraction * 10000).toFixed(1)} m`);

const area = polygonAreaKm2(ring.slice(0, -1));
ok("circle area is within the inscribed-polygon bound of pi r squared",
  area.areaKm2 < Math.PI * 100 && area.areaKm2 > Math.PI * 100 * 0.999,
  `${area.areaKm2.toFixed(4)} vs ${(Math.PI * 100).toFixed(4)} km2`);
ok("circle perimeter is close to 2 pi r",
  near(area.perimeterKm, 2 * Math.PI * 10, 0.05), `${area.perimeterKm.toFixed(4)} km`);
const u = areaUnits(1);
ok("a square kilometre is 100 hectares", near(u.hectares, 100, 1e-9));
ok("a square kilometre is 247.105 acres", near(u.acres, 247.105, 1e-3));

/* ---- path length -------------------------------------------------------- */
ok("path length sums its legs", near(pathLengthKm([A, mid, B]), 20, 1e-6));
ok("a one-point path has no length", pathLengthKm([A]) === 0);

/* ---- spatial index: the bounded-search defect ---------------------------- */
// 5,000 scattered points plus one deliberately far outside the first rings.
const pts = [];
for (let i = 0; i < 5000; i += 1) {
  pts.push([-7 + Math.random() * 9, 50 + Math.random() * 8.5]);
}
const target = [-3.0, 56.0];
const idx = new SpatialIndex(0.1);
pts.forEach((p, i) => idx.add(i, p[0], p[1]));
const measure = (i) => distanceKm(target[0], target[1], pts[i][0], pts[i][1]);
const viaIndex = idx.nearest(target[0], target[1], measure);
let brute = { id: -1, km: Infinity };
pts.forEach((p, i) => {
  const d = measure(i);
  if (d < brute.km) brute = { id: i, km: d };
});
ok("spatial index finds the same nearest as a brute-force sweep",
  viaIndex && viaIndex.id === brute.id,
  `index=${viaIndex && viaIndex.km.toFixed(4)} brute=${brute.km.toFixed(4)}`);

// A lone point far from everything must still be found, not silently missed.
const lonely = new SpatialIndex(0.1);
lonely.add(0, 2.0, 51.0);
const far = lonely.nearest(-6.5, 57.5, () => distanceKm(-6.5, 57.5, 2.0, 51.0));
ok("an isolated feature is still found after many rings",
  far !== null && far.id === 0, far ? `${far.km.toFixed(1)} km` : "not found");

/* ---- curvature ---------------------------------------------------------- */
const c54 = curvatureKm(54);
ok("prime vertical exceeds meridional at 54N", c54.primeVertical > c54.meridional);
ok("Gaussian radius at 54N matches R_UK",
  near(Math.sqrt(c54.meridional * c54.primeVertical), R_UK, 1e-3));

/* ---- report ------------------------------------------------------------- */
const failed = checks.filter((c) => !c.pass);
for (const c of checks) {
  console.log(`  [${c.pass ? "PASS" : "FAIL"}] ${c.name}${c.detail ? "  " + c.detail : ""}`);
}
console.log(`\n${checks.length - failed.length}/${checks.length} checks passed`);
process.exit(failed.length ? 1 : 0);
