# Earth model — the audit, the error budget, the decision

Measured 31 August 2026 against 4,000 random UK point pairs between 0.5 and
60 km, using the Vincenty inverse solution on WGS84 as truth.

## What the estate was doing

| Implementation | Radius | Error vs ellipsoid |
|---|---|---|
| `gridatlas` `ventus-corev8engine.js` | 6378.137 | −1,078 ppm mean, 141 m worst |
| `gis-sld-financial-sandbox` `atlasHaversineKm` | 6378.137 | as above |
| `gis-sld-financial-sandbox` `turf.distance` | 6371.0088 | −2,194 ppm mean, 206 m worst |
| `pipelinenews` `circuit_km` | 6371.0088 | as above, plus sampling error below |

Two implementations inside one tool, on two radii. 6378.137 is the WGS84
**equatorial** radius; 6371.0088 is the IUGG **mean**. They differ by 0.112%,
which is 11.2 m per 10 km and 112 m per 100 km.

## What is actually correct for the UK

The best single sphere for a latitude is the Gaussian mean radius of curvature,
√(M·N), where M is meridional and N prime-vertical:

| Latitude | M | N | Gaussian |
|---|---|---|---|
| 50°N | 6372.956 | 6390.702 | 6381.823 |
| 52°N | 6375.150 | 6391.435 | 6383.287 |
| 54°N | 6377.307 | 6392.156 | **6384.727** |
| 56°N | 6379.417 | 6392.861 | 6386.135 |
| 58°N | 6381.469 | 6393.546 | 6387.505 |

Both deployed constants are **too small** for these latitudes — 6378.137 by
0.103%, 6371.0088 by 0.215%.

| Sphere | Mean error | Worst error |
|---|---|---|
| 6371.0088 IUGG mean | −2,194 ppm | 206 m |
| 6378.137 WGS84 equatorial | −1,078 ppm | 141 m |
| 6383.635 UK minimax | −217 ppm | **91 m** |
| 6384.727 Gaussian at 54°N | **−46 ppm** | 102 m |

The error is also **anisotropic**: on 6378.137 a 10 km east–west leg at 54°N is
22.0 m short, while a 10 km north–south leg is 1.2 m long. A sphere cannot fix
that; only the ellipsoid can.

## Decision

**`R_ATLAS = 6378.137` stays the default.** Not because it is the most accurate
— it is not — but because it is what every deployed tool already uses, and every
published number depends on it. Changing it silently would move figures already
quoted from the Atlas, the sandbox and Pipeline News.

Accordingly:

- `distanceKm` defaults to `R_ATLAS`. Anything measuring alongside existing
  Ventus output must use it, so the numbers agree.
- `R_UK = 6384.7272` is exported for new work where accuracy matters more than
  backward agreement. It cuts mean error by a factor of 23.
- `distanceEllipsoidalKm` is the reference. Use it to validate, and use it for
  any figure that will be quoted to a third party.
- **`R_MEAN = 6371.0088` should not be used in new code.** It is exported only
  so existing Turf-based results can be reproduced and compared.

The gap between R_ATLAS and the truth — about 1 part in 1,000 — is well inside
the uncertainty of a REPD site centroid, which can sit hundreds of metres from
the point of connection. The constant is not the dominant error. What follows
is.

## The errors that actually dominate

Ranked by size, for the `circuit_km` case:

**1. Measuring to a vertex instead of to the line.** The published figure
measured to the nearest point of a sampled set, so it could only equal or
overstate. Vertex spacing on that sample: median 264 m, 90th percentile 872 m,
worst 1,115 m — a systematic overstatement of typically 130 m, up to 560 m.
That is five times the worst-case error from the radius choice.
Fixed by `distanceToSegmentKm`.

**2. Missing voltages.** The sample carried 400/275/132 only. 171 projects
(5.6%) are nearest to a 220 or 66 kV circuit, and 71 were pushed outside the
2 km band by the omission — Crimscote Solar Farm read 18.21 km against a true
0.89 km. That is a 17 km error from a data gap, not a maths one.

**3. A bounded nearest-neighbour search.** Ten projects reported a circuit that
was not the nearest, worst by 33.79 km, because the search stopped before it was
provably complete. `SpatialIndex.nearest` only stops once the best hit lies
inside the area already swept.

**4. Decimation.** The sample held 47,897 points where the source geometry holds
163,905 vertices over 149,340 segments — about a third.

**5. Coverage holes reported as numbers.** The Isle of Man has no vertices in
the layer, so Billown Solar Farm reported 76.55 km, its distance to the mainland
across open sea. Return null, never a large number.

The order matters: fixing the constant while leaving 1 to 5 in place would have
changed the fourth decimal place of a figure that was tens of kilometres wrong.

## Inscribed circles

Both the Atlas radius tool and the zonedraw elastic band draw a circle as an
inscribed polygon, which understates:

| Vertices | Radial error at 10 km | Area understated | Used by |
|---|---|---|---|
| 24 | 85.6 m | 1.138% | `ZONE_DRAW_VERTICES`, the elastic band |
| 64 | 12.0 m | 0.161% | `createGeoJSONCircle`, radius ≤ 500 km |
| 96 | 5.4 m | 0.071% | radius 500–5,000 km |
| 128 | 3.0 m | 0.040% | radius > 5,000 km |

The elastic band's 1.14% area understatement is larger than every distance error
discussed above. `geodesicCircle` defaults to 128 vertices, and
`circleApproximationError` returns the residual so a panel can report it.

## Layer inventory

The pinned voltage set, from the Overpass fetch scripts in
`repd_grid_atlasv8/scripts/`:

| Layer | Files | Features | Vertices | Segments |
|---|---|---|---|---|
| 11 kV (UKPN) | 1 | 15,126 | 15,126 | points |
| 33 kV (11 regions) | 11 | 5,748 | 110,305 | 104,557 |
| 66 kV | 1 | 1,171 | 18,341 | 17,170 |
| 132 kV | 1 | 6,227 | 81,327 | 75,100 |
| 220 kV | 1 | 126 | 1,875 | 1,749 |
| 275 kV | 1 | 2,935 | 24,856 | 21,921 |
| 400 kV | 1 | 4,106 | 37,506 | 33,400 |
| Substations | 1 | 5,800 | 5,800 | points |
| **Total** | | | **295,136** | **253,897** |

Pipeline News `202608311530` currently uses 66–400 kV, 149,340 segments.
**33 kV is the largest unused layer and the voltage most utility-scale solar
actually connects at.** Adding it is the next material improvement, worth more
than any further refinement of the Earth model.
