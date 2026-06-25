// slice.js — 3D cross-section. A plane (normal = axis, position = z) intersects each
// sphere in a circle of radius sqrt(r^2 - d^2); the slice is a 2D packing the 2D
// renderer draws directly. Carries through any scalar overlay layer for the sliced
// particles. The slice area-fraction recovers phi (mean-chord theorem for a plane).
window.VIZ = window.VIZ || {};

VIZ.computeSlice = function (data, axis, z) {
  const N = data.N, p = data.positions, dia = data.diameters;
  const inPlane = [0, 1, 2].filter(a => a !== axis);
  const Laxis = data.box[axis];
  const xs = [], ys = [], ds = [], idx = [];
  for (let i = 0; i < N; i++) {
    const r = dia[i] * 0.5;
    let dz = Math.abs(p[i * 3 + axis] - z); dz = Math.min(dz, Laxis - dz);  // periodic in z
    if (dz < r) { xs.push(p[i * 3 + inPlane[0]]); ys.push(p[i * 3 + inPlane[1]]);
                  ds.push(2 * Math.sqrt(r * r - dz * dz)); idx.push(i); }
  }
  const n = xs.length, pos = new Float32Array(n * 2), diam = new Float32Array(n);
  for (let i = 0; i < n; i++) { pos[i * 2] = xs[i]; pos[i * 2 + 1] = ys[i]; diam[i] = ds[i]; }
  const out = { N: n, ndim: 2, box: [data.box[inPlane[0]], data.box[inPlane[1]]],
                phi: data.phi, positions: pos, diameters: diam, layers: {}, _slice: { axis, z },
                _srcIdx: idx };
  // carry per-particle scalar layers (e.g. local_phi) onto the sliced subset
  for (const name in data.layers) {
    const L = data.layers[name];
    if (L.type === 'scalar_per_particle') {
      const v = new Float32Array(n); for (let i = 0; i < n; i++) v[i] = L.values[idx[i]];
      out.layers[name] = { type: L.type, values: v };
    }
  }
  return VIZ.addGhosts(out);   // periodic in-plane edge wrap
};
