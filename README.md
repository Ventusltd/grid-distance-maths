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
src/geodesy.mjs        canonical implementation, browser and Node
src/geodesy.py         the Python twin, for build-time payload scripts
docs/EARTH-MODEL.md    the audit, the error budget, and the decision
test/verify.mjs        34 checks against ellipsoidal truth and against the
                       specific defects this repo exists to prevent
test/verify_parity.py  446 checks that the two implementations return the
                       same number for the same input
test/parity_driver.mjs the JavaScript side of that comparison
test/verify_nearest.mjs 20,000 randomised layouts checking the nearest-feature
                       search against exhaustive truth
```

Run all three suites. No dependencies, and no network.

```
node   test/verify.mjs          34/34
node   test/verify_nearest.mjs  54/54
python test/verify_parity.py    446/446
```

Every payload builder in the estate is Python and every panel is JavaScript.
The twin exists so those two agree by construction; the parity suite exists so
that claim is checked rather than asserted. Change one file, change the other,
and let `verify_parity.py` prove it — it compares distances, bearings,
destinations, segment projections, circles, areas and nearest-feature search
across a spread of real UK sites plus the awkward cases (identical points, the
antimeridian, over the pole, a zero-length segment), and it compares the caveat
strings exactly, because a consumer that shows different wording in one runtime
than the other is a real defect, not a rounding one.

### Not here yet

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
6. **A nearest search may only stop when it can prove it.** The ring sweep
   terminates when the best hit found is inside the radius the swept box
   provably covers -- distance to the nearest box edge, with the east-west pair
   converted on the longitude scale at the highest latitude the box reaches.
   Anything looser silently overstates. See below.
7. **Distance says nothing about headroom.** Fault level and thermal headroom
   are properties of the network, not of the geometry. They depend on DNO data
   — source impedance, fault infeed, existing committed connections — and are
   established by a connection study. `STRAIGHT_LINE_CAVEAT.headroom` carries
   that sentence so no consumer has to invent it.

## The nearest-search defect, and why the tests missed it

`verify.mjs` asserted "spatial index finds the same nearest as a brute-force
sweep" from the first commit. It passed. The search was wrong in about one query
in nine anyway, because that assertion used one hand-made fixture whose layout
happened not to trip it.

Two faults, both found on 31 Aug 2026 after an external review flagged the
class of problem, both since fixed:

1. **The stop test converted both axes with `ky`.** A cell is `cell` *degrees*
   on both axes, but a degree of longitude is shorter than a degree of latitude
   and narrows towards the pole -- 0.588 of it at 54N, 0.500 at 60N. Using the
   latitude scale for the east-west span permitted stopping up to twice as early
   as the swept box justified.
2. **`addSegment` registered only the two endpoint cells.** A segment crossing a
   cell without ending in it was invisible from inside that cell.

Measured on 6,000 randomised layouts before the fix: **657 wrong, 10.95%**, the
worst reporting **65.42 km for a circuit 35.79 km away**. Both faults overstate
and never understate, which is what makes them dangerous -- a site reads as
further from the network than it is and nothing about the output looks wrong.
It is the same failure this repository was created to end, reappearing inside
the module meant to end it.

**No published number was affected.** The distances on globalgrid2050 were
produced by `pipelinenews` `grid-proximity/build_payload.py`, which has its own
nearest search and had already got both of these right -- it indexes by
bounding box, and its `swept_radius_km` takes the true distance to each of the
four box edges on both scales. The canonical module was behind its own consumer.
That formulation has now been ported here, so the library is at least as good as
the code it is meant to replace, and `verify_nearest.mjs` holds the line with
20,000 randomised layouts across the UK latitude band rather than one fixture.

The lesson is in the test, not the arithmetic: a fixture proves the code runs;
only randomised comparison against exhaustive truth proves it is right. Parity
between the two implementations proves nothing here either -- both were wrong in
the same way, so they agreed. `verify_parity.py` now checks its nearest cases
against brute force on the Python side as well.

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
