// slice.js — 3D cross-section. A plane (normal = axis, position = z) intersects each
// sphere in a circle of radius sqrt(r^2 - d^2); the slice is a 2D packing the 2D
// renderer draws directly. Carries through any scalar overlay layer for the sliced
// particles. The slice area-fraction recovers phi (mean-chord theorem for a plane).
window.VIZ = window.VIZ || {};

// --- ellipsoids (additive) ---------------------------------------------------------
// The plane x_axis = z cuts the ellipsoid (x-c)^T A^-1 (x-c) <= 1 in an ellipse, exactly. With n the normal
// axis, (a, b) the in-plane axes and h = z - c_n the signed distance of the plane from the centre:
//   it is cut iff h^2 < A_nn  (sqrt(A_nn) is the half-width of the ellipsoid along the normal);
//   centre of the cut  = (c_a, c_b) + h (A_an, A_bn) / A_nn      (shifted: the cut of a tilted body is off-centre);
//   shape of the cut   = (1 - h^2/A_nn) (A_uu - a a^T / A_nn),  a = (A_an, A_bn), A_uu the in-plane 2x2 block,
// so that the cut is { u : (u - centre)^T S^-1 (u - centre) <= 1 }. For a sphere this is the circle of radius
// sqrt(r^2 - h^2). Returned as [centre_a, centre_b, sxx, sxy, syy], or null when the plane misses the body.
VIZ.sliceEllipsoid = function (c, A6, axis, z, Laxis) {
  const inPlane = [0, 1, 2].filter(a => a !== axis), a = inPlane[0], b = inPlane[1];
  const ix = (i, j) => { if (i > j) { const t = i; i = j; j = t; } return i === 0 ? j : (i === 1 ? 2 + j : 5); };  // packed symmetric index
  const Ann = A6[ix(axis, axis)];
  let h = z - c[axis];
  if (Laxis) h -= Laxis * Math.round(h / Laxis);         // nearest periodic image along the normal
  if (!(h * h < Ann)) return null;
  const Aan = A6[ix(a, axis)], Abn = A6[ix(b, axis)], k = 1 - h * h / Ann;
  return [c[a] + h * Aan / Ann, c[b] + h * Abn / Ann,
          k * (A6[ix(a, a)] - Aan * Aan / Ann), k * (A6[ix(a, b)] - Aan * Abn / Ann), k * (A6[ix(b, b)] - Abn * Abn / Ann)];
};

VIZ.computeSliceEllipsoid = function (data, axis, z) {
  const N = data.N, p = data.positions, dia = data.diameters, A3 = data.A3;
  const inPlane = [0, 1, 2].filter(a => a !== axis), w = data.walls;
  const perAxis = w ? (w[axis] === 0) : true;            // no walls recorded -> periodic (as for spheres)
  const perIn = inPlane.map(dim => w ? (w[dim] === 0) : true);
  const Lin = [data.box[inPlane[0]], data.box[inPlane[1]]];
  const xs = [], ys = [], sh = [], idx = [];
  for (let i = 0; i < N; i++) {
    const cut = VIZ.sliceEllipsoid([p[i * 3], p[i * 3 + 1], p[i * 3 + 2]], A3.subarray(i * 6, i * 6 + 6), axis, z,
                                   perAxis ? data.box[axis] : 0);
    if (!cut) continue;
    for (let k = 0; k < 2; k++) if (perIn[k]) cut[k] -= Lin[k] * Math.floor(cut[k] / Lin[k]);   // centre back into the box
    xs.push(cut[0]); ys.push(cut[1]); sh.push(cut[2], cut[3], cut[4]); idx.push(i);
  }
  const n = xs.length, pos = new Float32Array(n * 2), diam = new Float32Array(n), diam3d = new Float32Array(n);
  const shape = new Float32Array(sh);
  for (let i = 0; i < n; i++) {
    pos[i * 2] = xs[i]; pos[i * 2 + 1] = ys[i]; diam3d[i] = dia[idx[i]];
    diam[i] = 2 * Math.pow(Math.max(0, shape[i * 3] * shape[i * 3 + 2] - shape[i * 3 + 1] * shape[i * 3 + 1]), 0.25);  // area-equivalent
  }
  const out = { N: n, ndim: 2, box: Lin, phi: data.phi, positions: pos, diameters: diam, shape, layers: {},
                _slice: { axis, z }, _srcIdx: idx, _diam3d: diam3d, _periodic: perIn };
  for (const name in data.layers) {
    const L = data.layers[name];
    if (L.type === 'scalar_per_particle') {
      const v = new Float32Array(n); for (let i = 0; i < n; i++) v[i] = L.values[idx[i]];
      out.layers[name] = { type: L.type, values: v };
    }
  }
  if (w) { const sb = VIZ.sliceBoundary(data, axis, z, inPlane, Lin[0], Lin[1]); out.boundary = sb.boundary; out._periodic = sb.periodic; }
  const res = VIZ.addGhosts(out);
  res.boundary = out.boundary;
  return res;
};

// The container's intersection with the slice plane, for spheres, ellipsoids and rods alike. Solid = hard wall
// (curved container -> a circle that shrinks/grows with z for a sphere); dashed = periodic. Returns the boundary
// segments and, for the two in-plane axes, whether they are periodic (curved and hard axes get no periodic ghost).
VIZ.sliceBoundary = function (data, axis, z, inPlane, Lx, Ly) {
  const w = data.walls, a = inPlane[0], b = inPlane[1];
  const t = w[0] < 0 ? (w[0] === -1 ? 2 : -w[0]) : 0;
  const periodic = [a, b].map(dim => dim < t ? false : (w[dim] === 0));
  let boundary;
  if (t === 0) {
    boundary = { segments: VIZ.boxEdges(Lx, Ly, w[a] === 0, w[b] === 0) };
  } else {
    const R = data.box[0] / 2, c = data.box[0] / 2;
    const aC = a < t, bC = b < t, axisC = axis < t;
    const hz = Math.sqrt(Math.max(0, R * R - (z - c) * (z - c)));  // sphere cross-section radius at z
    if (aC && bC) {                       // both in-plane dims curved -> a circle
      const Reff = axisC ? hz : R;        // sphere: shrinks with z; cylinder (perp slice): full R
      boundary = { segments: Reff > 1e-9 ? [{ path: VIZ.circlePath(c, c, Reff), dashed: false }] : [] };
    } else if (aC || bC) {                // cylinder sliced along its axis: one curved dim + one flat
      const h = axisC ? hz : R, segs = [];
      if (aC) { segs.push({ path: [[c - h, 0], [c - h, Ly]], dashed: false }, { path: [[c + h, 0], [c + h, Ly]], dashed: false });
                const dY = w[b] === 0; segs.push({ path: [[0, 0], [Lx, 0]], dashed: dY }, { path: [[0, Ly], [Lx, Ly]], dashed: dY }); }
      else    { segs.push({ path: [[0, c - h], [Lx, c - h]], dashed: false }, { path: [[0, c + h], [Lx, c + h]], dashed: false });
                const dX = w[a] === 0; segs.push({ path: [[0, 0], [0, Ly]], dashed: dX }, { path: [[Lx, 0], [Lx, Ly]], dashed: dX }); }
      boundary = { segments: segs };
    } else {
      boundary = { segments: VIZ.boxEdges(Lx, Ly, w[a] === 0, w[b] === 0) };
    }
  }
  return { boundary, periodic };
};

VIZ.computeSlice = function (data, axis, z) {
  if (data.A3) return VIZ.computeSliceEllipsoid(data, axis, z);
  if (data.rod) return VIZ.computeSliceRod(data, axis, z);          // rods: rod.js
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
  // cross-section boundary: the container's intersection with the slice plane (VIZ.sliceBoundary below)
  if (data.walls) {
    const sb = VIZ.sliceBoundary(data, axis, z, inPlane, out.box[0], out.box[1]);
    out._periodic = sb.periodic; out.boundary = sb.boundary;
  }
  const res = VIZ.addGhosts(out);   // periodic in-plane edge wrap
  res.boundary = out.boundary;
  return res;
};
