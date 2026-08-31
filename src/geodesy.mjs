/**
 * Grid Distance Maths — canonical geodesy.
 *
 * Every Ventus tool that measures a distance on the grid should call this, so
 * that the same two points always produce the same number. See
 * docs/EARTH-MODEL.md for why R_ATLAS is the default and when to reach past it.
 *
 * Nothing here returns a cable route. Every function returns straight-line or
 * great-circle geometry, with no wayleave, crossing, terrain or consent
 * content. Consumers must carry that caveat to the user.
 */

/* ---- constants ---------------------------------------------------------- */

/** WGS84 semi-major axis. The constant every deployed Ventus tool already uses. */
export const R_ATLAS = 6378.137;

/** IUGG mean radius. Turf.js default. Reads 0.112% shorter than R_ATLAS. */
export const R_MEAN = 6371.0088;

/**
 * Gaussian mean radius of curvature at 54°N, the UK centroid. The most accurate
 * single sphere for these latitudes: mean error 46 ppm against the ellipsoid,
 * against 1,078 ppm for R_ATLAS. Opt in where accuracy matters more than
 * agreeing with already-published numbers.
 */
export const R_UK = 6384.7272;

export const WGS84 = Object.freeze({
  a: 6378.137,
  f: 1 / 298.257223563,
  get b() { return this.a * (1 - this.f); },
  get e2() { return this.f * (2 - this.f); },
});

const DEG = Math.PI / 180;

/* ---- curvature ---------------------------------------------------------- */

/**
 * Meridional (north–south) and prime-vertical (east–west) radii of curvature
 * at a latitude. These are what make a local tangent plane correct rather than
 * merely convenient.
 */
export function curvatureKm(latDeg) {
  const s = Math.sin(latDeg * DEG);
  const t = 1 - WGS84.e2 * s * s;
  return {
    meridional: (WGS84.a * (1 - WGS84.e2)) / t ** 1.5,
    primeVertical: WGS84.a / Math.sqrt(t),
  };
}

/** Scale factors converting a degree of lon/lat to km at this latitude. */
export function localScaleKm(latDeg) {
  const { meridional, primeVertical } = curvatureKm(latDeg);
  return {
    kx: primeVertical * Math.cos(latDeg * DEG) * DEG,
    ky: meridional * DEG,
  };
}

/* ---- distance ----------------------------------------------------------- */

/**
 * Great-circle distance in km. The atan2 form, which stays stable at both very
 * small and near-antipodal separations where the asin form loses precision.
 *
 * Identical in form and default constant to ventus-corev8engine.js haversine()
 * and gis-sld-v5-drawing.js atlasHaversineKm().
 */
export function distanceKm(lon1, lat1, lon2, lat2, radius = R_ATLAS) {
  const dLat = (lat2 - lat1) * DEG;
  const dLon = (lon2 - lon1) * DEG;
  const x =
    Math.sin(dLat / 2) ** 2 +
    Math.cos(lat1 * DEG) * Math.cos(lat2 * DEG) * Math.sin(dLon / 2) ** 2;
  return radius * 2 * Math.atan2(Math.sqrt(x), Math.sqrt(1 - x));
}

/**
 * Vincenty inverse solution on the WGS84 ellipsoid — the reference this
 * repository measures its spheres against. Millimetre-grade, ~30x slower than
 * the sphere, and non-convergent for near-antipodal pairs (returns NaN rather
 * than a wrong answer). Use it to validate, or where a number will be quoted.
 */
export function distanceEllipsoidalKm(lon1, lat1, lon2, lat2) {
  const { a, f, b } = WGS84;
  const L = (lon2 - lon1) * DEG;
  const U1 = Math.atan((1 - f) * Math.tan(lat1 * DEG));
  const U2 = Math.atan((1 - f) * Math.tan(lat2 * DEG));
  const sU1 = Math.sin(U1); const cU1 = Math.cos(U1);
  const sU2 = Math.sin(U2); const cU2 = Math.cos(U2);
  let lambda = L; let prev; let iter = 0;
  let sinSigma; let cosSigma; let sigma; let cos2Alpha; let cos2SigmaM;
  do {
    const sl = Math.sin(lambda); const cl = Math.cos(lambda);
    sinSigma = Math.sqrt((cU2 * sl) ** 2 + (cU1 * sU2 - sU1 * cU2 * cl) ** 2);
    if (sinSigma === 0) return 0;
    cosSigma = sU1 * sU2 + cU1 * cU2 * cl;
    sigma = Math.atan2(sinSigma, cosSigma);
    const sinAlpha = (cU1 * cU2 * sl) / sinSigma;
    cos2Alpha = 1 - sinAlpha * sinAlpha;
    cos2SigmaM = cos2Alpha === 0 ? 0 : cosSigma - (2 * sU1 * sU2) / cos2Alpha;
    const C = (f / 16) * cos2Alpha * (4 + f * (4 - 3 * cos2Alpha));
    prev = lambda;
    lambda = L + (1 - C) * f * sinAlpha
      * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM ** 2)));
  } while (Math.abs(lambda - prev) > 1e-12 && ++iter < 200);
  if (iter >= 200) return NaN;
  const u2 = (cos2Alpha * (a * a - b * b)) / (b * b);
  const A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
  const B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
  const dSigma = B * sinSigma * (cos2SigmaM + (B / 4)
    * (cosSigma * (-1 + 2 * cos2SigmaM ** 2)
      - (B / 6) * cos2SigmaM * (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)));
  return b * A * (sigma - dSigma);
}

/** Cumulative length of a polyline, in km. */
export function pathLengthKm(coords, radius = R_ATLAS) {
  if (!Array.isArray(coords) || coords.length < 2) return 0;
  let total = 0;
  for (let i = 1; i < coords.length; i += 1) {
    total += distanceKm(coords[i - 1][0], coords[i - 1][1], coords[i][0], coords[i][1], radius);
  }
  return total;
}

/* ---- bearing and destination -------------------------------------------- */

export function initialBearingDeg(lon1, lat1, lon2, lat2) {
  const p1 = lat1 * DEG; const p2 = lat2 * DEG;
  const dl = (lon2 - lon1) * DEG;
  const y = Math.sin(dl) * Math.cos(p2);
  const x = Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl);
  return (Math.atan2(y, x) / DEG + 360) % 360;
}

/** Direct spherical destination: travel `km` from a point on `bearingDeg`. */
export function destinationPoint(lon, lat, km, bearingDeg, radius = R_ATLAS) {
  const ad = km / radius;
  const brg = bearingDeg * DEG;
  const p1 = lat * DEG;
  const p2 = Math.asin(Math.sin(p1) * Math.cos(ad) + Math.cos(p1) * Math.sin(ad) * Math.cos(brg));
  const l2 = lon * DEG + Math.atan2(
    Math.sin(brg) * Math.sin(ad) * Math.cos(p1),
    Math.cos(ad) - Math.sin(p1) * Math.sin(p2),
  );
  return [l2 / DEG, p2 / DEG];
}

/* ---- point to line ------------------------------------------------------ */

/**
 * Perpendicular distance from a point to a single segment, and the foot of that
 * perpendicular.
 *
 * This is the function whose absence caused the circuit_km defect: measuring to
 * a segment's endpoints instead of to the segment can only overstate, and by up
 * to half the vertex spacing.
 *
 * The projection runs on a local tangent plane built from the WGS84 radii of
 * curvature at the query point's own latitude, so the foot is geometrically
 * right; the returned distance is then measured with `distanceKm`, so it is
 * directly comparable with every other number this module produces.
 */
export function distanceToSegmentKm(lon, lat, aLon, aLat, bLon, bLat, radius = R_ATLAS) {
  const { kx, ky } = localScaleKm(lat);
  const ax = (aLon - lon) * kx; const ay = (aLat - lat) * ky;
  const bx = (bLon - lon) * kx; const by = (bLat - lat) * ky;
  const dx = bx - ax; const dy = by - ay;
  const len2 = dx * dx + dy * dy;
  let t = 0;
  if (len2 > 0) {
    t = -(ax * dx + ay * dy) / len2;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
  }
  const foot = [aLon + (bLon - aLon) * t, aLat + (bLat - aLat) * t];
  return { km: distanceKm(lon, lat, foot[0], foot[1], radius), foot, t };
}

/** The nearest point on a polyline, and how far away it is. */
export function distanceToLineKm(lon, lat, coords, radius = R_ATLAS) {
  let best = null;
  for (let i = 0; i < coords.length - 1; i += 1) {
    const hit = distanceToSegmentKm(
      lon, lat, coords[i][0], coords[i][1], coords[i + 1][0], coords[i + 1][1], radius,
    );
    if (!best || hit.km < best.km) best = { ...hit, segmentIndex: i };
  }
  return best;
}

/* ---- circles and area --------------------------------------------------- */

/**
 * A geodesic circle as a closed ring. The ring is INSCRIBED, so it understates:
 * at n vertices the worst radial error is R(1 − cos(π/n)) and the area is short
 * by 1 − (n/2π)·sin(2π/n).
 *
 *   n = 24   85.6 m inward at 10 km, area short 1.138%   (the zonedraw default)
 *   n = 64   12.0 m inward at 10 km, area short 0.161%
 *   n = 128   3.0 m inward at 10 km, area short 0.040%
 *
 * Default 128: the cost is trivial and it keeps the drawn ring honest against
 * the number printed beside it.
 */
export function geodesicCircle(lon, lat, radiusKm, vertices = 128, radius = R_ATLAS) {
  const ring = [];
  for (let i = 0; i < vertices; i += 1) {
    ring.push(destinationPoint(lon, lat, radiusKm, (i / vertices) * 360, radius));
  }
  ring.push(ring[0]);
  return ring;
}

/** Inscribed-polygon error for a given vertex count, so callers can report it. */
export function circleApproximationError(vertices) {
  return {
    radialFraction: 1 - Math.cos(Math.PI / vertices),
    areaShortFraction: 1 - (vertices / (2 * Math.PI)) * Math.sin((2 * Math.PI) / vertices),
  };
}

/** Spherical polygon area (km²) and perimeter (km) for a closed or open ring. */
export function polygonAreaKm2(coords, radius = R_ATLAS) {
  if (!Array.isArray(coords) || coords.length < 3) return { areaKm2: 0, perimeterKm: 0 };
  let sum = 0;
  for (let i = 0; i < coords.length; i += 1) {
    const j = (i + 1) % coords.length;
    const xi = coords[i][0] * DEG; const yi = coords[i][1] * DEG;
    const xj = coords[j][0] * DEG; const yj = coords[j][1] * DEG;
    sum += (xj - xi) * (2 + Math.sin(yi) + Math.sin(yj));
  }
  const areaKm2 = (Math.abs(sum) * radius * radius) / 2;
  return { areaKm2, perimeterKm: pathLengthKm([...coords, coords[0]], radius) };
}

/** Familiar units for an area, for panels that report land take. */
export function areaUnits(areaKm2) {
  const m2 = areaKm2 * 1e6;
  return {
    km2: areaKm2,
    hectares: m2 / 1e4,
    acres: m2 / 4046.85642,
    squareMiles: areaKm2 * 0.386102,
    squareMetres: m2,
  };
}

/* ---- nearest-feature search --------------------------------------------- */

/**
 * A uniform-cell spatial index. Built once over the network, it turns a
 * whole-layer sweep into a ring search over a handful of cells.
 *
 * The ring expands until the best hit found is provably inside the swept area —
 * stopping earlier is the bug that left ten projects reporting a circuit tens of
 * kilometres further away than the one actually nearest them.
 */
export class SpatialIndex {
  constructor(cellDegrees = 0.1) {
    this.cell = cellDegrees;
    this.buckets = new Map();
  }

  static key(i, j) { return `${i}|${j}`; }

  cellOf(lon, lat) {
    return [Math.floor(lat / this.cell), Math.floor(lon / this.cell)];
  }

  add(id, lon, lat) {
    const [i, j] = this.cellOf(lon, lat);
    const k = SpatialIndex.key(i, j);
    if (!this.buckets.has(k)) this.buckets.set(k, []);
    this.buckets.get(k).push(id);
  }

  /** Add every cell a segment's endpoints fall in, so long spans stay findable. */
  addSegment(id, aLon, aLat, bLon, bLat) {
    this.add(id, aLon, aLat);
    const [i1, j1] = this.cellOf(aLon, aLat);
    const [i2, j2] = this.cellOf(bLon, bLat);
    if (i1 !== i2 || j1 !== j2) this.add(id, bLon, bLat);
  }

  /**
   * @param {(id:number)=>number} measure distance in km for a candidate id
   * @returns {{id:number, km:number}|null}
   */
  nearest(lon, lat, measure, maxRings = 90) {
    const [ci, cj] = this.cellOf(lon, lat);
    const { ky } = localScaleKm(lat);
    let best = null;
    for (let ring = 0; ring < maxRings; ring += 1) {
      for (let i = ci - ring; i <= ci + ring; i += 1) {
        for (let j = cj - ring; j <= cj + ring; j += 1) {
          if (ring && Math.abs(i - ci) !== ring && Math.abs(j - cj) !== ring) continue;
          const bucket = this.buckets.get(SpatialIndex.key(i, j));
          if (!bucket) continue;
          for (const id of bucket) {
            const km = measure(id);
            if (best === null || km < best.km) best = { id, km };
          }
        }
      }
      // Only safe to stop when the best hit lies inside the area already swept.
      if (best && best.km <= ring * this.cell * ky * 0.999) return best;
    }
    return best;
  }
}

/* ---- guardrail ---------------------------------------------------------- */

/**
 * The caveat every consumer must show. Exported so no tool has to reword it,
 * and so a tool cannot quietly ship without it.
 *
 * `headroom` is the one a distance cannot answer at all. Fault level and
 * thermal headroom are a property of the network, not of the geometry, and the
 * only honest thing a proximity number can say about them is that they require
 * DNO data and a study.
 */
export const STRAIGHT_LINE_CAVEAT = Object.freeze({
  distance: "Straight-line distance to mapped geometry. Not a cable route, not a "
    + "connection length, and no wayleave, crossing, terrain or consent content.",
  substation: "A mapped substation point does not confirm capacity, voltage "
    + "suitability, connection rights, queue position or acceptance by any network party.",
  coverage: "Absence from a mapped layer is not absence on the ground.",
  headroom: "Fault level and thermal headroom cannot be inferred from distance. "
    + "They depend on DNO network data such as source impedance, fault infeed and "
    + "existing committed connections, and are established by a connection study, "
    + "not by geometry.",
});
