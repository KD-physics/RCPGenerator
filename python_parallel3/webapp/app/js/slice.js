// slice.js — 3D cross-section. A plane (normal = axis, position = z) intersects each
// sphere in a circle of radius sqrt(r^2 - d^2); the slice is a 2D packing the 2D
// renderer draws directly. Carries through any scalar overlay layer for the sliced
// particles. The slice area-fraction recovers phi (mean-chord theorem for a plane).
window.VIZ = window.VIZ || {};

VIZ.computeSlice = function (data, axis, z) {
  const N = data.N, p = data.positions, dia = data.diameters;
  const inPlane = [0, 1, 2].filter(a => a !== axis);
  const Laxis = data.box[axis];
  const xs = [], ys = [], ds = [], idx = [], d3 = [];   // d3 = TRUE 3D sphere diameter (for color-by-diameter)
  for (let i = 0; i < N; i++) {
    const r = dia[i] * 0.5;
    let dz = Math.abs(p[i * 3 + axis] - z); dz = Math.min(dz, Laxis - dz);  // periodic in z
    if (dz < r) { xs.push(p[i * 3 + inPlane[0]]); ys.push(p[i * 3 + inPlane[1]]);
                  ds.push(2 * Math.sqrt(r * r - dz * dz)); idx.push(i); d3.push(dia[i]); }
  }
  const n = xs.length, pos = new Float32Array(n * 2), diam = new Float32Array(n), diam3d = new Float32Array(n);
  for (let i = 0; i < n; i++) { pos[i * 2] = xs[i]; pos[i * 2 + 1] = ys[i]; diam[i] = ds[i]; diam3d[i] = d3[i]; }
  const out = { N: n, ndim: 2, box: [data.box[inPlane[0]], data.box[inPlane[1]]],
                phi: data.phi, positions: pos, diameters: diam, layers: {}, _slice: { axis, z },
                _srcIdx: idx, _diam3d: diam3d };
  // carry per-particle scalar layers (e.g. local_phi) onto the sliced subset
  for (const name in data.layers) {
    const L = data.layers[name];
    if (L.type === 'scalar_per_particle') {
      const v = new Float32Array(n); for (let i = 0; i < n; i++) v[i] = L.values[idx[i]];
      out.layers[name] = { type: L.type, values: v };
    }
  }
  // cross-section boundary: the container's intersection with the slice plane. Solid
  // = hard wall (curved container -> a circle that shrinks/grows with z for a sphere);
  // dashed = periodic. Absent walls -> leave undefined (renderer uses the dashed box).
  if (data.walls) {
    const w = data.walls, a = inPlane[0], b = inPlane[1], Lx = out.box[0], Ly = out.box[1];
    const t = w[0] < 0 ? -w[0] : 0;
    out._periodic = [a, b].map(dim => dim < t ? false : (w[dim] === 0));   // curved/hard -> no periodic ghost
    if (t === 0) {
      out.boundary = { segments: VIZ.boxEdges(Lx, Ly, w[a] === 0, w[b] === 0) };
    } else {
      const R = data.box[0] / 2, c = data.box[0] / 2;
      const aC = a < t, bC = b < t, axisC = axis < t;
      const hz = Math.sqrt(Math.max(0, R * R - (z - c) * (z - c)));  // sphere cross-section radius at z
      if (aC && bC) {                       // both in-plane dims curved -> a circle
        const Reff = axisC ? hz : R;        // sphere: shrinks with z; cylinder (perp slice): full R
        out.boundary = { segments: Reff > 1e-9 ? [{ path: VIZ.circlePath(c, c, Reff), dashed: false }] : [] };
      } else if (aC || bC) {                // cylinder sliced along its axis: one curved dim + one flat
        const h = axisC ? hz : R, segs = [];
        if (aC) { segs.push({ path: [[c - h, 0], [c - h, Ly]], dashed: false }, { path: [[c + h, 0], [c + h, Ly]], dashed: false });
                  const dY = w[b] === 0; segs.push({ path: [[0, 0], [Lx, 0]], dashed: dY }, { path: [[0, Ly], [Lx, Ly]], dashed: dY }); }
        else    { segs.push({ path: [[0, c - h], [Lx, c - h]], dashed: false }, { path: [[0, c + h], [Lx, c + h]], dashed: false });
                  const dX = w[a] === 0; segs.push({ path: [[0, 0], [0, Ly]], dashed: dX }, { path: [[Lx, 0], [Lx, Ly]], dashed: dX }); }
        out.boundary = { segments: segs };
      } else {
        out.boundary = { segments: VIZ.boxEdges(Lx, Ly, w[a] === 0, w[b] === 0) };
      }
    }
  }
  const res = VIZ.addGhosts(out);   // periodic in-plane edge wrap
  res.boundary = out.boundary;
  return res;
};
