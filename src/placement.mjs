/* Placement: where a thing goes on a plane, and what "near" means there.
 *
 * The JavaScript half. placement.py is the other half and neither is
 * authoritative: they must agree, and test/verify_placement.mjs is what holds
 * them to it. Same rule as geodesy.
 *
 *   golden(key)      an archive. Adjacency is ORDER. The Kuiper arrangement.
 *   hilbert(d, o)    a swept space. Adjacency is one parameter, one step.
 *   lattice(i, rad)  mixed-radix decode: the space being arranged.
 *   morph(a, b, t)   the same points on their way from one to the other.
 */

/* the golden angle, pi*(3-sqrt(5)) = 2.39996... rad. The most irrational
 * angle there is, which is the point: a rational one eventually makes spokes. */
export const GOLDEN = Math.PI * (3 - Math.sqrt(5));

export function golden(key, spacing = 1) {
  if (key < 0) throw new RangeError("key must be >= 0");
  const r = spacing * Math.sqrt(key), t = key * GOLDEN;
  return [r * Math.cos(t), r * Math.sin(t)];
}

export function hilbertOrder(n) {
  if (n <= 0) throw new RangeError("n must be > 0");
  let order = 0;
  while ((1 << order) * (1 << order) < n) order += 1;
  return order;
}

export function hilbert(d, order) {
  const side = 1 << order;
  if (d < 0 || d >= side * side) throw new RangeError("d out of range");
  let x = 0, y = 0, t = d, s = 1;
  while (s < side) {
    const rx = Math.floor(t / 2) % 2 ? 1 : 0;
    const ry = (t ^ rx) % 2 ? 1 : 0;
    if (ry === 0) {
      if (rx === 1) { x = s - 1 - x; y = s - 1 - y; }
      const tmp = x; x = y; y = tmp;
    }
    x += s * rx;
    y += s * ry;
    t = Math.floor(t / 4);
    s <<= 1;
  }
  return [x, y];
}

export function lattice(index, radices) {
  if (index < 0) throw new RangeError("index must be >= 0");
  const out = [];
  let r = index;
  for (const radix of radices) {
    if (radix < 1) throw new RangeError("every radix must be >= 1");
    out.push(r % radix);
    r = Math.floor(r / radix);
  }
  return out;
}

export function morph(a, b, t, ease = true) {
  if (!(t >= 0 && t <= 1)) throw new RangeError("t must be in [0, 1]");
  const k = ease ? t * t * (3 - 2 * t) : t;
  return [a[0] + (b[0] - a[0]) * k, a[1] + (b[1] - a[1]) * k];
}
