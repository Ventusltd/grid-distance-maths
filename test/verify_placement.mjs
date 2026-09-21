/* Holds placement.mjs to placement.py. Neither is authoritative.
 *
 *   node test/verify_placement.mjs
 *
 * It writes the JavaScript answers, runs the Python half over the same inputs,
 * and compares. A check that reached nothing FAILS: if the Python side does
 * not run, or if either side produced no rows, this exits non-zero rather than
 * reporting agreement over an empty set.
 */
import { execFileSync } from "node:child_process";
import { GOLDEN, golden, hilbert, hilbertOrder, lattice, morph } from "../src/placement.mjs";

const TOL = 1e-12;
const ORDER = 6;                       // 64 x 64 = 4096 points
const KEYS = 4096;
const RADICES = [3, 6, 15, 17, 101];

function jsRows() {
  const rows = [];
  for (let k = 0; k < KEYS; k++) {
    const [gx, gy] = golden(k, 1.7);
    const [hx, hy] = hilbert(k, ORDER);
    const lat = lattice(k, RADICES);
    const [mx, my] = morph([gx, gy], [hx, hy], (k % 11) / 10);
    rows.push([gx, gy, hx, hy, ...lat, mx, my]);
  }
  return rows;
}

const PY = `
import json, sys
sys.path.insert(0, "src")
from placement import golden, hilbert, lattice, morph
ORDER, KEYS = ${ORDER}, ${KEYS}
RAD = ${JSON.stringify(RADICES)}
out = []
for k in range(KEYS):
    g = golden(k, 1.7)
    h = hilbert(k, ORDER)
    lat = lattice(k, RAD)
    m = morph(g, h, (k % 11) / 10)
    out.append([g[0], g[1], h[0], h[1], *lat, m[0], m[1]])
print(json.dumps(out))
`;

const js = jsRows();
let py;
try {
  py = JSON.parse(execFileSync("python", ["-c", PY], {
    encoding: "utf8", maxBuffer: 1 << 28,
  }));
} catch (err) {
  console.error("FAIL: the Python half did not run, so nothing was compared.");
  console.error("A check that reached nothing is not a pass.");
  console.error(String(err.message).slice(0, 400));
  process.exit(1);
}

if (js.length === 0 || py.length === 0) {
  console.error("FAIL: one side produced no rows (js %d, py %d).", js.length, py.length);
  process.exit(1);
}
if (js.length !== py.length) {
  console.error("FAIL: row counts differ (js %d, py %d).", js.length, py.length);
  process.exit(1);
}

let worst = 0, worstAt = -1, differ = 0;
for (let i = 0; i < js.length; i++) {
  if (js[i].length !== py[i].length) {
    console.error("FAIL: row %d has different widths.", i);
    process.exit(1);
  }
  for (let c = 0; c < js[i].length; c++) {
    const d = Math.abs(js[i][c] - py[i][c]);
    if (d > worst) { worst = d; worstAt = i; }
    if (d > TOL) differ += 1;
  }
}

console.log("placement parity");
console.log("  rows compared      %d, of %d columns each", js.length, js[0].length);
console.log("  golden angle       %s", GOLDEN.toPrecision(12));
console.log("  hilbert order      %d  (side %d)", ORDER, 1 << ORDER);
console.log("  values over %s     %d", TOL, differ);
console.log("  worst difference   %s  (row %d)", worst.toExponential(3), worstAt);

if (differ > 0) {
  console.error("FAIL: the two implementations disagree.");
  process.exit(1);
}
console.log("PASS: placement.mjs and placement.py agree on every value.");
