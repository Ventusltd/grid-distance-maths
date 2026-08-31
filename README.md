# Grid Distance Maths

One implementation of the geometry the grid runs on, so every Ventus tool
answers the same question with the same number.

Distance to a circuit. Distance to a substation. Radius search. Geodesic
circles. Polygon area. Point-to-line projection. From 11 kV to 400 kV and
beyond, onshore and offshore.

The purpose is to support the energy transition through the maths of the grid:
if a developer, an engineer or an investor asks how far a site is from the
network, the answer should not depend on which of our applications they happen
to have open.

## Why this repo exists

The same calculation was implemented four times across the estate, on three
different Earth radii, and the answers disagreed.

| Where | Constant | Method |
|---|---|---|
| `gridatlas` &middot; `ventus-corev8engine.js` | 6378.137 | haversine, `atan2` form |
| `gis-sld-financial-sandbox` &middot; `atlasHaversineKm` | 6378.137 | haversine, `atan2` form |
| `gis-sld-financial-sandbox` &middot; `turf.distance` | 6371.0088 | Turf default, used for nearest-feature |
| `pipelinenews` &middot; `circuit_km` | 6371.0088 | haversine to a decimated point sample |

Two consequences, both measured, both real:

- The sandbox mixes radii **inside one tool** — route lengths on 6378.137,
  nearest-feature distances on 6371.0088.
- `circuit_km` measured to the nearest *sampled vertex* rather than to the line,
  so it could only ever equal or overstate the distance to the conductor, and it
  omitted 220 kV and 66 kV. 71 projects were pushed outside the 2 km band by the
  omission. One of them, Crimscote Solar Farm, read 18.21 km and is 0.89 km from
  a 66 kV line.

See [`docs/EARTH-MODEL.md`](docs/EARTH-MODEL.md) for the measurements and the
decision on which radius to standardise.

## What is here

```
src/geodesy.mjs      canonical implementation, browser and Node
docs/EARTH-MODEL.md  the audit, the error budget, and the decision
test/verify.mjs      34 checks against ellipsoidal truth and against the
                     specific defects this repo exists to prevent
```

Run the tests with `node test/verify.mjs`. No dependencies.

### Not here yet

- `src/geodesy.py`, a Python twin, so the build-time scripts and the browser
  agree by construction rather than by inspection. Every payload builder in the
  estate is Python; today they each carry their own copy of the formula.
- The 33 kV and 11 kV layers wired into consumers. See the inventory in
  `docs/EARTH-MODEL.md`: 33 kV is 104,557 segments and the voltage most
  utility-scale solar actually connects at, and it is currently unused.
- A published pin of the layer set, so a consumer can state which network
  vintage a distance was measured against.

## The rules this repo holds

1. **One constant, named, never inlined.** `R_ATLAS = 6378.137` is what every
   deployed Ventus tool already uses. It is not the most accurate sphere for UK
   latitudes — see the decision — but a silent change would move every published
   number, so it stays the default and the better answers are opt-in.
2. **Distance to a line means distance to the line**, not to whichever vertex
   somebody happened to sample. `distanceToLineKm` projects onto the segment.
3. **A sphere is a stated approximation.** `distanceEllipsoidalKm` is the
   reference; the sphere is the fast path, and its error is documented, not
   assumed negligible.
4. **Straight-line is not a route.** Nothing here returns a cable route, a
   wayleave, or a connection length. Every consumer must carry that caveat.
5. **Absence from a layer is not absence on the ground.** Coverage gaps are
   reported as null, never as a large number.

## Scope

In: geodesic distance, bearings, destination points, point-to-segment and
point-to-polyline projection, geodesic circles and buffers, spherical polygon
area and perimeter, nearest-feature search with a spatial index, voltage layer
handling from 11 kV upward, substation proximity.

Out: routing, wayleaves, terrain, thermal ratings, load flow, capacity. Those
are engineering, not geometry, and belong with the tools that own them.

## Licence and provenance

Network geometry used in validation is OpenStreetMap-derived, ODbL-1.0,
© OpenStreetMap contributors. Project data is DESNZ REPD under Open Government
Licence v3.0. This repository holds the mathematics, not the data.
