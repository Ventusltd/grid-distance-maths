"""Grid Distance Maths - canonical geodesy, Python twin of src/geodesy.mjs.

Every Ventus build script that measures a distance on the grid should import
this, so that a payload built in Python and a panel drawn in JavaScript answer
the same question with the same number by construction rather than by
inspection.

The two files are kept in lockstep by test/verify_parity.py, which runs the
same inputs through both and fails on any disagreement beyond floating-point
noise. If you change one, change the other, and let the parity test prove it.

Nothing here returns a cable route. Every function returns straight-line or
great-circle geometry, with no wayleave, crossing, terrain or consent content.
Consumers must carry that caveat to the user.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

# ---- constants -----------------------------------------------------------

#: WGS84 semi-major axis. The constant every deployed Ventus tool already uses.
R_ATLAS = 6378.137

#: IUGG mean radius. Turf.js default. Reads 0.112% shorter than R_ATLAS.
R_MEAN = 6371.0088

#: Gaussian mean radius of curvature at 54 degrees N, the UK centroid. The most
#: accurate single sphere for these latitudes: mean error 46 ppm against the
#: ellipsoid, against 1,078 ppm for R_ATLAS. Opt in where accuracy matters more
#: than agreeing with already-published numbers.
R_UK = 6384.7272

WGS84_A = 6378.137
WGS84_F = 1 / 298.257223563
WGS84_B = WGS84_A * (1 - WGS84_F)
WGS84_E2 = WGS84_F * (2 - WGS84_F)

DEG = math.pi / 180


# ---- curvature -----------------------------------------------------------

def curvature_km(lat_deg: float) -> dict:
    """Meridional (north-south) and prime-vertical (east-west) radii of
    curvature at a latitude. These are what make a local tangent plane correct
    rather than merely convenient."""
    s = math.sin(lat_deg * DEG)
    t = 1 - WGS84_E2 * s * s
    return {
        "meridional": (WGS84_A * (1 - WGS84_E2)) / t ** 1.5,
        "prime_vertical": WGS84_A / math.sqrt(t),
    }


def local_scale_km(lat_deg: float) -> dict:
    """Scale factors converting a degree of lon/lat to km at this latitude."""
    c = curvature_km(lat_deg)
    return {
        "kx": c["prime_vertical"] * math.cos(lat_deg * DEG) * DEG,
        "ky": c["meridional"] * DEG,
    }


# ---- distance ------------------------------------------------------------

def distance_km(lon1: float, lat1: float, lon2: float, lat2: float,
                radius: float = R_ATLAS) -> float:
    """Great-circle distance in km. The atan2 form, which stays stable at both
    very small and near-antipodal separations where the asin form loses
    precision.

    Identical in form and default constant to ventus-corev8engine.js
    haversine() and gis-sld-v5-drawing.js atlasHaversineKm()."""
    d_lat = (lat2 - lat1) * DEG
    d_lon = (lon2 - lon1) * DEG
    x = (math.sin(d_lat / 2) ** 2
         + math.cos(lat1 * DEG) * math.cos(lat2 * DEG) * math.sin(d_lon / 2) ** 2)
    return radius * 2 * math.atan2(math.sqrt(x), math.sqrt(1 - x))


def distance_ellipsoidal_km(lon1: float, lat1: float,
                            lon2: float, lat2: float) -> float:
    """Vincenty inverse solution on the WGS84 ellipsoid - the reference this
    repository measures its spheres against. Millimetre-grade, ~30x slower than
    the sphere, and non-convergent for near-antipodal pairs (returns nan rather
    than a wrong answer). Use it to validate, or where a number will be quoted.
    """
    a, f, b = WGS84_A, WGS84_F, WGS84_B
    L = (lon2 - lon1) * DEG
    U1 = math.atan((1 - f) * math.tan(lat1 * DEG))
    U2 = math.atan((1 - f) * math.tan(lat2 * DEG))
    sU1, cU1 = math.sin(U1), math.cos(U1)
    sU2, cU2 = math.sin(U2), math.cos(U2)
    lam = L
    sin_sigma = cos_sigma = sigma = cos2_alpha = cos2_sigma_m = 0.0
    converged = False
    for _ in range(200):
        sl, cl = math.sin(lam), math.cos(lam)
        sin_sigma = math.sqrt((cU2 * sl) ** 2 + (cU1 * sU2 - sU1 * cU2 * cl) ** 2)
        if sin_sigma == 0:
            return 0.0
        cos_sigma = sU1 * sU2 + cU1 * cU2 * cl
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = (cU1 * cU2 * sl) / sin_sigma
        cos2_alpha = 1 - sin_alpha * sin_alpha
        cos2_sigma_m = 0.0 if cos2_alpha == 0 else cos_sigma - (2 * sU1 * sU2) / cos2_alpha
        C = (f / 16) * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))
        prev = lam
        lam = L + (1 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (cos2_sigma_m + C * cos_sigma
                                     * (-1 + 2 * cos2_sigma_m ** 2)))
        if abs(lam - prev) <= 1e-12:
            converged = True
            break
    if not converged:
        return float("nan")

    u2 = (cos2_alpha * (a * a - b * b)) / (b * b)
    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    d_sigma = B * sin_sigma * (
        cos2_sigma_m + (B / 4) * (
            cos_sigma * (-1 + 2 * cos2_sigma_m ** 2)
            - (B / 6) * cos2_sigma_m * (-3 + 4 * sin_sigma ** 2)
            * (-3 + 4 * cos2_sigma_m ** 2)))
    return b * A * (sigma - d_sigma)


def path_length_km(coords: Sequence[Sequence[float]],
                   radius: float = R_ATLAS) -> float:
    """Cumulative length of a polyline, in km."""
    if not coords or len(coords) < 2:
        return 0.0
    total = 0.0
    for i in range(1, len(coords)):
        total += distance_km(coords[i - 1][0], coords[i - 1][1],
                             coords[i][0], coords[i][1], radius)
    return total


# ---- bearing and destination ---------------------------------------------

def initial_bearing_deg(lon1: float, lat1: float,
                        lon2: float, lat2: float) -> float:
    p1, p2 = lat1 * DEG, lat2 * DEG
    dl = (lon2 - lon1) * DEG
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.atan2(y, x) / DEG + 360) % 360


def destination_point(lon: float, lat: float, km: float, bearing_deg: float,
                      radius: float = R_ATLAS) -> list:
    """Direct spherical destination: travel `km` from a point on `bearing_deg`."""
    ad = km / radius
    brg = bearing_deg * DEG
    p1 = lat * DEG
    p2 = math.asin(math.sin(p1) * math.cos(ad)
                   + math.cos(p1) * math.sin(ad) * math.cos(brg))
    l2 = lon * DEG + math.atan2(
        math.sin(brg) * math.sin(ad) * math.cos(p1),
        math.cos(ad) - math.sin(p1) * math.sin(p2))
    return [l2 / DEG, p2 / DEG]


# ---- point to line -------------------------------------------------------

def distance_to_segment_km(lon: float, lat: float,
                           a_lon: float, a_lat: float,
                           b_lon: float, b_lat: float,
                           radius: float = R_ATLAS) -> dict:
    """Perpendicular distance from a point to a single segment, and the foot of
    that perpendicular.

    This is the function whose absence caused the circuit_km defect: measuring
    to a segment endpoint instead of to the segment can only overstate, and by
    up to half the vertex spacing.

    The projection runs on a local tangent plane built from the WGS84 radii of
    curvature at the query point own latitude, so the foot is geometrically
    right; the returned distance is then measured with distance_km, so it is
    directly comparable with every other number this module produces."""
    scale = local_scale_km(lat)
    kx, ky = scale["kx"], scale["ky"]
    ax, ay = (a_lon - lon) * kx, (a_lat - lat) * ky
    bx, by = (b_lon - lon) * kx, (b_lat - lat) * ky
    dx, dy = bx - ax, by - ay
    len2 = dx * dx + dy * dy
    t = 0.0
    if len2 > 0:
        t = -(ax * dx + ay * dy) / len2
        t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    foot = [a_lon + (b_lon - a_lon) * t, a_lat + (b_lat - a_lat) * t]
    return {"km": distance_km(lon, lat, foot[0], foot[1], radius),
            "foot": foot, "t": t}


def distance_to_line_km(lon: float, lat: float,
                        coords: Sequence[Sequence[float]],
                        radius: float = R_ATLAS):
    """The nearest point on a polyline, and how far away it is."""
    best = None
    for i in range(len(coords) - 1):
        hit = distance_to_segment_km(lon, lat, coords[i][0], coords[i][1],
                                     coords[i + 1][0], coords[i + 1][1], radius)
        if best is None or hit["km"] < best["km"]:
            best = dict(hit, segment_index=i)
    return best


# ---- circles and area ----------------------------------------------------

def geodesic_circle(lon: float, lat: float, radius_km: float,
                    vertices: int = 128, radius: float = R_ATLAS) -> list:
    """A geodesic circle as a closed ring. The ring is INSCRIBED, so it
    understates: at n vertices the worst radial error is R(1 - cos(pi/n)) and
    the area is short by 1 - (n/2pi)*sin(2pi/n).

        n = 24   85.6 m inward at 10 km, area short 1.138%   (zonedraw default)
        n = 64   12.0 m inward at 10 km, area short 0.161%
        n = 128   3.0 m inward at 10 km, area short 0.040%

    Default 128: the cost is trivial and it keeps the drawn ring honest against
    the number printed beside it."""
    ring = [destination_point(lon, lat, radius_km, (i / vertices) * 360, radius)
            for i in range(vertices)]
    ring.append(ring[0])
    return ring


def circle_approximation_error(vertices: int) -> dict:
    """Inscribed-polygon error for a vertex count, so callers can report it."""
    return {
        "radial_fraction": 1 - math.cos(math.pi / vertices),
        "area_short_fraction": 1 - (vertices / (2 * math.pi))
        * math.sin((2 * math.pi) / vertices),
    }


def polygon_area_km2(coords: Sequence[Sequence[float]],
                     radius: float = R_ATLAS) -> dict:
    """Spherical polygon area (km2) and perimeter (km) for a closed or open ring."""
    if not coords or len(coords) < 3:
        return {"area_km2": 0.0, "perimeter_km": 0.0}
    total = 0.0
    n = len(coords)
    for i in range(n):
        j = (i + 1) % n
        xi, yi = coords[i][0] * DEG, coords[i][1] * DEG
        xj, yj = coords[j][0] * DEG, coords[j][1] * DEG
        total += (xj - xi) * (2 + math.sin(yi) + math.sin(yj))
    area_km2 = (abs(total) * radius * radius) / 2
    return {"area_km2": area_km2,
            "perimeter_km": path_length_km(list(coords) + [coords[0]], radius)}


def area_units(area_km2: float) -> dict:
    """Familiar units for an area, for panels that report land take."""
    m2 = area_km2 * 1e6
    return {
        "km2": area_km2,
        "hectares": m2 / 1e4,
        "acres": m2 / 4046.85642,
        "square_miles": area_km2 * 0.386102,
        "square_metres": m2,
    }


# ---- nearest-feature search ----------------------------------------------

class SpatialIndex:
    """A uniform-cell spatial index. Built once over the network, it turns a
    whole-layer sweep into a ring search over a handful of cells.

    The ring expands until the best hit found is provably inside the swept area
    - stopping earlier is the bug that left ten projects reporting a circuit
    tens of kilometres further away than the one actually nearest them."""

    def __init__(self, cell_degrees: float = 0.1):
        self.cell = cell_degrees
        self.buckets: dict = {}

    def cell_of(self, lon: float, lat: float):
        return (math.floor(lat / self.cell), math.floor(lon / self.cell))

    def add(self, ident, lon: float, lat: float):
        self.buckets.setdefault(self.cell_of(lon, lat), []).append(ident)

    def add_segment(self, ident, a_lon, a_lat, b_lon, b_lat):
        """Register a segment in EVERY cell its bounding box covers.

        Indexing only the two endpoint cells is not enough: a segment that
        enters a cell and leaves it again without either end landing inside is
        invisible from that cell, so a query sitting right beside the conductor
        misses it until the ring search happens to reach an endpoint. That is a
        silent overstatement, which is the failure mode this module exists to
        prevent.

        The bounding box is deliberately conservative rather than an exact
        supercover walk: it can add a diagonal segment to a few cells it does
        not actually cross, which costs one extra distance measurement each and
        can never cost correctness, because nearest() measures true distance to
        every candidate it pulls out.
        """
        i1, j1 = self.cell_of(a_lon, a_lat)
        i2, j2 = self.cell_of(b_lon, b_lat)
        for i in range(min(i1, i2), max(i1, i2) + 1):
            for j in range(min(j1, j2), max(j1, j2) + 1):
                self.buckets.setdefault((i, j), []).append(ident)

    def swept_clearance_km(self, lon: float, lat: float, ring: int) -> float:
        """The radius around the query point that a Chebyshev ring sweep has
        PROVABLY covered. Stopping before the best hit is inside this is the
        bounded-search defect; stopping later is merely slow.

        Two things make this subtler than ring * cell:

        1. A cell is `cell` DEGREES on both axes, but a degree of longitude is
           shorter than a degree of latitude and narrows towards the pole --
           0.588 of it at 54N, 0.500 at 60N. Converting both axes with ky, as
           this did until it was measured, permits stopping up to twice as early
           as the box justifies. A randomised sweep over 6,000 layouts then
           returned the wrong nearest feature in 10.95% of them, the worst
           reporting 65.4 km for a circuit 35.8 km away. Every such error
           overstates distance.

        2. Sweeping to Chebyshev distance `ring` covers an axis-aligned box of
           cells, not a disc. The query point sits somewhere inside its own
           cell, not at its centre, so the guaranteed radius is the distance to
           the NEAREST edge of that box -- which can be almost a whole cell less
           than ring * cell on the side the point is closest to.

        Taking the true distance to each of the four edges is both correct and
        tighter than assuming the worst corner. This mirrors swept_radius_km in
        the pipelinenews grid-proximity builder, which reached the better
        formulation first; the canonical module should not be behind its own
        consumer.

        kx is evaluated at the highest latitude the box reaches, where a degree
        of longitude is narrowest, so the east-west guarantee holds across the
        whole box rather than only at the query latitude.
        """
        if ring <= 0:
            return 0.0
        ci, cj = self.cell_of(lon, lat)
        lat_lo = (ci - ring) * self.cell
        lat_hi = (ci + ring + 1) * self.cell
        lon_lo = (cj - ring) * self.cell
        lon_hi = (cj + ring + 1) * self.cell
        ky = local_scale_km(lat)["ky"]
        worst_lat = min(max(abs(lat_lo), abs(lat_hi)), 89.9)
        kx = local_scale_km(worst_lat)["kx"]
        return min((lat - lat_lo) * ky, (lat_hi - lat) * ky,
                   (lon - lon_lo) * kx, (lon_hi - lon) * kx)

    def nearest(self, lon: float, lat: float,
                measure: Callable[[object], float], max_rings: int = 90):
        ci, cj = self.cell_of(lon, lat)
        best = None
        seen = set()
        for ring in range(max_rings):
            for i in range(ci - ring, ci + ring + 1):
                for j in range(cj - ring, cj + ring + 1):
                    if ring and abs(i - ci) != ring and abs(j - cj) != ring:
                        continue
                    bucket = self.buckets.get((i, j))
                    if not bucket:
                        continue
                    for ident in bucket:
                        # A segment spans several cells, so the same id surfaces
                        # more than once. Measuring it once is the same answer
                        # for less work.
                        if ident in seen:
                            continue
                        seen.add(ident)
                        km = measure(ident)
                        if best is None or km < best["km"]:
                            best = {"id": ident, "km": km}
            if best and best["km"] <= self.swept_clearance_km(lon, lat, ring):
                return best
        return best


# ---- guardrail -----------------------------------------------------------

#: The caveat every consumer must show. Exported so no tool has to reword it,
#: and so a tool cannot quietly ship without it.
#:
#: `headroom` is the one a distance cannot answer at all. Fault level and
#: thermal headroom are a property of the network, not of the geometry, and the
#: only honest thing a proximity number can say about them is that they require
#: DNO data and a study.
STRAIGHT_LINE_CAVEAT = {
    "distance": ("Straight-line distance to mapped geometry. Not a cable route, "
                 "not a connection length, and no wayleave, crossing, terrain or "
                 "consent content."),
    "substation": ("A mapped substation point does not confirm capacity, voltage "
                   "suitability, connection rights, queue position or acceptance "
                   "by any network party."),
    "coverage": "Absence from a mapped layer is not absence on the ground.",
    "headroom": ("Fault level and thermal headroom cannot be inferred from "
                 "distance. They depend on DNO network data such as source "
                 "impedance, fault infeed and existing committed connections, "
                 "and are established by a connection study, not by geometry."),
}
