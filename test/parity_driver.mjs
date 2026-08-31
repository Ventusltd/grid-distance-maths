/**
 * Parity driver: read a JSON case list on stdin, answer with the JavaScript
 * result for each case on stdout.
 *
 * It exists so verify_parity.py can compare the two implementations on
 * identical inputs without either side knowing the other's expected answers.
 * Nothing here computes anything itself — every case dispatches straight into
 * src/geodesy.mjs, so a divergence is a real divergence.
 */

import {
  R_ATLAS, R_MEAN, R_UK,
  curvatureKm, localScaleKm,
  distanceKm, distanceEllipsoidalKm, pathLengthKm,
  initialBearingDeg, destinationPoint,
  distanceToSegmentKm, distanceToLineKm,
  geodesicCircle, circleApproximationError,
  polygonAreaKm2, areaUnits,
  SpatialIndex, STRAIGHT_LINE_CAVEAT,
} from "../src/geodesy.mjs";

const OPS = {
  constants: () => ({ R_ATLAS, R_MEAN, R_UK }),
  curvature: (a) => curvatureKm(a.lat),
  localScale: (a) => localScaleKm(a.lat),
  distance: (a) => distanceKm(a.lon1, a.lat1, a.lon2, a.lat2, a.radius ?? R_ATLAS),
  ellipsoidal: (a) => {
    const v = distanceEllipsoidalKm(a.lon1, a.lat1, a.lon2, a.lat2);
    return Number.isNaN(v) ? null : v;
  },
  pathLength: (a) => pathLengthKm(a.coords, a.radius ?? R_ATLAS),
  bearing: (a) => initialBearingDeg(a.lon1, a.lat1, a.lon2, a.lat2),
  destination: (a) => destinationPoint(a.lon, a.lat, a.km, a.bearing, a.radius ?? R_ATLAS),
  segment: (a) => {
    const r = distanceToSegmentKm(a.lon, a.lat, a.aLon, a.aLat, a.bLon, a.bLat, a.radius ?? R_ATLAS);
    return { km: r.km, foot: r.foot, t: r.t };
  },
  line: (a) => {
    const r = distanceToLineKm(a.lon, a.lat, a.coords, a.radius ?? R_ATLAS);
    return { km: r.km, foot: r.foot, t: r.t, segmentIndex: r.segmentIndex };
  },
  circle: (a) => geodesicCircle(a.lon, a.lat, a.radiusKm, a.vertices, a.radius ?? R_ATLAS),
  circleError: (a) => circleApproximationError(a.vertices),
  polygon: (a) => polygonAreaKm2(a.coords, a.radius ?? R_ATLAS),
  areaUnits: (a) => areaUnits(a.areaKm2),
  caveat: () => STRAIGHT_LINE_CAVEAT,
  nearest: (a) => {
    const index = new SpatialIndex(a.cell ?? 0.1);
    a.features.forEach((f, i) => index.addSegment(i, f[0], f[1], f[2], f[3]));
    const hit = index.nearest(a.lon, a.lat, (i) => {
      const f = a.features[i];
      return distanceToSegmentKm(a.lon, a.lat, f[0], f[1], f[2], f[3]).km;
    });
    return hit === null ? null : { id: hit.id, km: hit.km };
  },
};

let raw = "";
process.stdin.setEncoding("utf8");
process.stdin.on("data", (chunk) => { raw += chunk; });
process.stdin.on("end", () => {
  const cases = JSON.parse(raw);
  const out = cases.map((c) => {
    const op = OPS[c.op];
    if (!op) throw new Error(`unknown op ${c.op}`);
    return op(c.args ?? {});
  });
  process.stdout.write(JSON.stringify(out));
});
