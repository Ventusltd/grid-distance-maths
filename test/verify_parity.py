"""Prove src/geodesy.py and src/geodesy.mjs are the same mathematics.

The Python twin exists so that a payload built by a script and a panel drawn in
a browser agree by construction. That claim is only worth anything if something
checks it, so this does: it builds a case list, runs it through both
implementations, and fails on any disagreement larger than floating-point noise.

    python test/verify_parity.py

Pure stdlib on the Python side; the JavaScript side is `node test/parity_driver.mjs`,
which does nothing but dispatch into src/geodesy.mjs.

Tolerance is 1e-9 relative for distances (nanometre-scale on a kilometre) and
exact for the caveat strings, because a consumer that shows a differently worded
caveat in one runtime than the other is a real defect, not a rounding one.
"""

import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import geodesy as G  # noqa: E402

REL_TOL = 1e-9
ABS_TOL = 1e-12

# A spread of real UK grid geometry plus the awkward cases: the equator, a
# meridian crossing, a pole-adjacent pair, and a zero-length segment.
SITES = [
    ("Thorpe Marsh BESS", -1.085062, 53.580258),
    ("Crimscote Solar", -1.663000, 52.140000),
    ("Beinn an Tuirc", -5.585000, 55.560000),
    ("Cleve Hill Solar", 0.900000, 51.340000),
    ("Dogger Bank A", 1.900000, 54.750000),
    ("Lands End", -5.716000, 50.066000),
    ("Unst", -0.850000, 60.760000),
]

CASES = []
EXPECT = []


def add(op, args, py):
    CASES.append({"op": op, "args": args})
    EXPECT.append(py)


def close(a, b, label, path=""):
    """Structural comparison. Numbers within tolerance, everything else exact."""
    if isinstance(a, dict) and isinstance(b, dict):
        # Key names differ by language convention (primeVertical / prime_vertical),
        # so compare by position-independent normalised key.
        na = {k.replace("_", "").lower(): v for k, v in a.items()}
        nb = {k.replace("_", "").lower(): v for k, v in b.items()}
        if set(na) != set(nb):
            return f"{label}{path}: key mismatch {sorted(na)} vs {sorted(nb)}"
        for k in na:
            bad = close(na[k], nb[k], label, f"{path}.{k}")
            if bad:
                return bad
        return None
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return f"{label}{path}: length {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            bad = close(x, y, label, f"{path}[{i}]")
            if bad:
                return bad
        return None
    if isinstance(a, bool) or isinstance(b, bool):
        return None if a == b else f"{label}{path}: {a!r} vs {b!r}"
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if math.isnan(a) and math.isnan(b):
            return None
        if math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            return None
        return f"{label}{path}: {a!r} vs {b!r}  (delta {abs(a - b):.3e})"
    if a is None and b is None:
        return None
    return None if a == b else f"{label}{path}: {a!r} vs {b!r}"


# ---- build the case list -------------------------------------------------

add("constants", {},
    {"R_ATLAS": G.R_ATLAS, "R_MEAN": G.R_MEAN, "R_UK": G.R_UK})

add("caveat", {}, G.STRAIGHT_LINE_CAVEAT)

for name, lon, lat in SITES:
    add("curvature", {"lat": lat}, G.curvature_km(lat))
    add("localScale", {"lat": lat}, G.local_scale_km(lat))

# Every ordered site pair, on all three radii.
for i, (n1, lo1, la1) in enumerate(SITES):
    for j, (n2, lo2, la2) in enumerate(SITES):
        if i == j:
            continue
        for r in (G.R_ATLAS, G.R_MEAN, G.R_UK):
            add("distance",
                {"lon1": lo1, "lat1": la1, "lon2": lo2, "lat2": la2, "radius": r},
                G.distance_km(lo1, la1, lo2, la2, r))
        ell = G.distance_ellipsoidal_km(lo1, la1, lo2, la2)
        add("ellipsoidal", {"lon1": lo1, "lat1": la1, "lon2": lo2, "lat2": la2},
            None if math.isnan(ell) else ell)
        add("bearing", {"lon1": lo1, "lat1": la1, "lon2": lo2, "lat2": la2},
            G.initial_bearing_deg(lo1, la1, lo2, la2))

# Degenerate and boundary distances.
for args in [
    {"lon1": 0.0, "lat1": 0.0, "lon2": 0.0, "lat2": 0.0},        # identical
    {"lon1": -0.001, "lat1": 51.5, "lon2": 0.001, "lat2": 51.5},  # metres apart
    {"lon1": -179.9, "lat1": 0.0, "lon2": 179.9, "lat2": 0.0},    # antimeridian
    {"lon1": 0.0, "lat1": 89.9, "lon2": 180.0, "lat2": 89.9},     # over the pole
]:
    add("distance", args,
        G.distance_km(args["lon1"], args["lat1"], args["lon2"], args["lat2"]))

# Point-to-segment, including the zero-length segment and both clamped ends.
SEGMENTS = [
    (-1.09, 53.57, -1.07, 53.59),
    (-1.09, 53.58, -1.09, 53.58),   # zero length
    (-2.00, 53.00, -1.00, 54.00),
    (0.90, 51.34, 0.91, 51.35),
]
for name, lon, lat in SITES:
    for (a_lon, a_lat, b_lon, b_lat) in SEGMENTS:
        py = G.distance_to_segment_km(lon, lat, a_lon, a_lat, b_lon, b_lat)
        add("segment", {"lon": lon, "lat": lat, "aLon": a_lon, "aLat": a_lat,
                        "bLon": b_lon, "bLat": b_lat},
            {"km": py["km"], "foot": py["foot"], "t": py["t"]})

# A 132 kV-shaped polyline, measured from every site.
LINE = [[-1.20, 53.50], [-1.15, 53.55], [-1.10, 53.58], [-1.02, 53.61], [-0.95, 53.66]]
for name, lon, lat in SITES:
    py = G.distance_to_line_km(lon, lat, LINE)
    add("line", {"lon": lon, "lat": lat, "coords": LINE},
        {"km": py["km"], "foot": py["foot"], "t": py["t"],
         "segmentIndex": py["segment_index"]})

add("pathLength", {"coords": LINE}, G.path_length_km(LINE))
add("pathLength", {"coords": [[0.0, 0.0]]}, G.path_length_km([[0.0, 0.0]]))

for v in (8, 24, 64, 128):
    add("circleError", {"vertices": v}, G.circle_approximation_error(v))

for name, lon, lat in SITES[:3]:
    for km in (0.5, 10.0, 50.0):
        add("circle", {"lon": lon, "lat": lat, "radiusKm": km, "vertices": 24},
            G.geodesic_circle(lon, lat, km, 24))
        for brg in (0.0, 47.5, 180.0, 359.9):
            add("destination", {"lon": lon, "lat": lat, "km": km, "bearing": brg},
                G.destination_point(lon, lat, km, brg))

RING = G.geodesic_circle(-1.085062, 53.580258, 3.0, 64)
add("polygon", {"coords": RING}, G.polygon_area_km2(RING))
add("polygon", {"coords": [[-1.1, 53.5], [-1.0, 53.5], [-1.0, 53.6], [-1.1, 53.6]]},
    G.polygon_area_km2([[-1.1, 53.5], [-1.0, 53.5], [-1.0, 53.6], [-1.1, 53.6]]))
add("areaUnits", {"areaKm2": 12.3456}, G.area_units(12.3456))

# Nearest-feature search, including the isolated feature many rings out that
# the bounded-search defect used to miss.
FEATURES = [
    [-1.20, 53.50, -1.15, 53.55],
    [-1.10, 53.58, -1.02, 53.61],
    [0.90, 51.34, 0.91, 51.35],
    [-5.60, 55.55, -5.58, 55.57],
    [-0.85, 60.75, -0.84, 60.77],
]
for name, lon, lat in SITES:
    index = G.SpatialIndex(0.1)
    for i, f in enumerate(FEATURES):
        index.add_segment(i, f[0], f[1], f[2], f[3])

    def measure(i, lon=lon, lat=lat):
        f = FEATURES[i]
        return G.distance_to_segment_km(lon, lat, f[0], f[1], f[2], f[3])["km"]

    hit = index.nearest(lon, lat, measure)
    add("nearest", {"lon": lon, "lat": lat, "features": FEATURES, "cell": 0.1},
        None if hit is None else {"id": hit["id"], "km": hit["km"]})


# ---- run the JavaScript side and compare ---------------------------------

def main():
    driver = os.path.join(HERE, "parity_driver.mjs")
    try:
        proc = subprocess.run(
            ["node", driver], input=json.dumps(CASES).encode("utf-8"),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError:
        print("SKIP  node is not on PATH; parity cannot be checked")
        return 2
    except subprocess.CalledProcessError as exc:
        print("FAIL  parity driver did not run")
        print(exc.stderr.decode("utf-8", "replace"))
        return 1

    js = json.loads(proc.stdout.decode("utf-8"))
    if len(js) != len(EXPECT):
        print("FAIL  case count mismatch: python %d, javascript %d"
              % (len(EXPECT), len(js)))
        return 1

    failures = []
    by_op = {}
    for case, py, jsv in zip(CASES, EXPECT, js):
        bad = close(py, jsv, case["op"])
        by_op.setdefault(case["op"], [0, 0])
        by_op[case["op"]][0] += 1
        if bad:
            failures.append(bad)
        else:
            by_op[case["op"]][1] += 1

    width = max(len(k) for k in by_op)
    for op in sorted(by_op):
        total, passed = by_op[op]
        mark = "PASS" if passed == total else "FAIL"
        print("  [%s] %-*s  %d/%d" % (mark, width, op, passed, total))

    print()
    if failures:
        print("%d of %d parity checks FAILED" % (len(failures), len(CASES)))
        for f in failures[:25]:
            print("   " + f)
        if len(failures) > 25:
            print("   ... and %d more" % (len(failures) - 25))
        return 1

    print("%d/%d parity checks passed - geodesy.py and geodesy.mjs agree"
          % (len(CASES), len(CASES)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
