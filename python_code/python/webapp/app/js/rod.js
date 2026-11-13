// rod.js — rods (additive). A rod is a straight section of half-length hl about its centre c along the unit axis
// u, thickened by the radius R: CAPPED = every point within R of the axis segment (half-sphere ends; in 2D a
// stadium), UNCAPPED = the flat-ended cylinder |t| <= hl, distance from the axis <= R (in 2D a rectangle).
// A bundle carries manifest.rods = { capped, length, radius } (length and radius in units of each particle's D,
// each a number or a file of N float32) and the same `orientations` entry as ellipsoids; the axis is column 0
// of the particle's rotation matrix. The loader turns that into data.rod = { capped, halfLength, radius, axis }.
//
// THE CUT OF A ROD BY A PLANE, exactly. Let n be the plane's normal, h the signed distance of the plane from
// the rod's centre, un = u.n, s = |u - un n| (the length of the axis' projection into the plane), and in the
// plane take coordinates (x, y) from the projected centre, x along the projected axis, y across it. A point
// (x, y) of the plane has axial coordinate t = s x + h un, and squared distance from the axis LINE
// y^2 + (un x - s h)^2. Therefore the cut is { |y| <= w(x) }, with w(x)^2 = R^2 - f(x) and
//   straight section, |t| <= hl, i.e. x in [xlo, xhi] = [(-hl - h un)/s, (hl - h un)/s]:  f = (un x - s h)^2
//       -> an arc of the ellipse centred at x = s h/un with semi-axes R/|un| along x and R along y
//          (un = 0: the two straight lines y = +-sqrt(R^2 - h^2); s = 0: the circle of radius R);
//   beyond the + end (x > xhi), capped only:  f = (x - hl s)^2 + (h - hl un)^2
//       -> an arc of the circle the plane cuts from the end sphere: centre x = hl s, radius^2 = R^2 - (h - hl un)^2;
//   beyond the - end (x < xlo), capped only:  the same with hl -> -hl.
// f is continuous across xlo and xhi, so the three arcs join; an uncapped rod has only the middle piece, closed
// by the straight cuts of its two flat ends at x = xlo and x = xhi. The outline is drawn as a polygon with
// vertices ON these arcs and a vertex exactly at every corner. Circle arcs have vertices 2 pi / ROD_SEG apart
// (the resolution used for ellipses); on the ellipse arc the step in its parameter is smaller by
// sqrt(R/semi-major axis) = sqrt(|un|), so that every chord stays within (1 - cos(pi/ROD_SEG)) R = 5.4e-4 R of
// its arc however oblique the cut. A 2D rod is the case h = 0, un = 0, s = 1 (the plane contains the axis).
window.VIZ = window.VIZ || {};
VIZ.ROD_SEG = 96;

// Upper half (y >= 0) of the cut outline in the local frame above, as [[x, y], ...] from the largest x to the
// smallest; null when the plane misses the rod.
VIZ.rodSectionUpper = function (s, un, h, hl, R, capped) {
  const step = 2 * Math.PI / VIZ.ROD_SEG, up = [];
  const push = (x, y) => {
    const m = up.length;
    if (m && Math.abs(up[m - 1][0] - x) + Math.abs(up[m - 1][1] - y) < 1e-13 * (R + hl)) return;   // joins: no double vertex
    up.push([x, y]);
  };
  const arc = (xc, ax, ay, th0, th1, fine) => {   // (xc + ax cos th, ay sin th), th from th0 to th1 (0 <= th0 < th1 <= pi)
    const n = Math.max(1, Math.ceil((th1 - th0) / step * (fine || 1) - 1e-9));
    for (let k = 0; k <= n; k++) { const th = th0 + (th1 - th0) * k / n; push(xc + ax * Math.cos(th), ay * Math.sin(th)); }
  };
  const tc = h * un;                           // axial coordinate of the projected centre
  let xlo, xhi;                                // the straight section's strip in x
  if (s > 1e-12) { xlo = (-hl - tc) / s; xhi = (hl - tc) / s; }
  else if (tc > hl) { xlo = -Infinity; xhi = -Infinity; }      // axis along the normal: the whole plane is beyond the + end,
  else if (tc < -hl) { xlo = Infinity; xhi = Infinity; }       // ... beyond the - end,
  else { xlo = -Infinity; xhi = Infinity; }                    // ... or within the straight section
  if (capped) {                                // + end sphere, where x >= xhi
    const r2 = R * R - (h - hl * un) * (h - hl * un);
    if (r2 > 0) { const r = Math.sqrt(r2), q = (xhi - hl * s) / r; if (q < 1) arc(hl * s, r, r, 0, Math.acos(Math.max(-1, q))); }
  }
  if (xhi > xlo) {                             // straight section
    if (Math.abs(un) < 1e-8) {                 // axis in the plane: two straight sides
      const w2 = R * R - s * s * h * h;
      if (w2 > 0) { const w = Math.sqrt(w2); push(xhi, w); push(xlo, w); }
    } else {                                   // ellipse arc: un x - s h = sg R cos th, y = R sin th
      const sg = un > 0 ? 1 : -1, au = Math.abs(un), x0 = s * h / un;
      const ca = Math.min(1, (au * xhi - sg * s * h) / R), cb = Math.max(-1, (au * xlo - sg * s * h) / R);
      if (ca > cb) arc(x0, R / au, R, Math.acos(ca), Math.acos(cb), 1 / Math.sqrt(au));
    }
  }
  if (capped) {                                // - end sphere, where x <= xlo
    const r2 = R * R - (h + hl * un) * (h + hl * un);
    if (r2 > 0) { const r = Math.sqrt(r2), q = (xlo + hl * s) / r; if (q > -1) arc(-hl * s, r, r, Math.acos(Math.min(1, q)), Math.PI); }
  }
  return up.length >= 2 ? up : null;
};

// Closed outline (counter-clockwise, first point not repeated) in plane coordinates: c2 + x e1 + y e2.
VIZ.rodSectionPolygon = function (c2, e1, s, un, h, hl, R, capped) {
  const up = VIZ.rodSectionUpper(s, un, h, hl, R, capped);
  if (!up) return null;
  const e2 = [-e1[1], e1[0]], poly = [];
  for (let k = 0; k < up.length; k++) poly.push([c2[0] + up[k][0] * e1[0] + up[k][1] * e2[0], c2[1] + up[k][0] * e1[1] + up[k][1] * e2[1]]);
  for (let k = up.length - 1; k >= 0; k--) {
    if (up[k][1] < 1e-13 * (R + hl)) continue;                  // a point on the x axis is its own mirror image
    poly.push([c2[0] + up[k][0] * e1[0] - up[k][1] * e2[0], c2[1] + up[k][0] * e1[1] - up[k][1] * e2[1]]);
  }
  return poly.length >= 3 ? poly : null;
};

// The cut of one 3D rod (centre c, axis u, half-length hl, radius R) by the plane x_axis = z: the polygon in the
// plane's two coordinates (the remaining axes in increasing order), or null when the plane misses the rod.
// Laxis > 0: the rod's nearest periodic image along the normal is the one cut.
VIZ.sliceRod = function (c, u, hl, R, capped, axis, z, Laxis) {
  const inPlane = [0, 1, 2].filter(a => a !== axis), a = inPlane[0], b = inPlane[1];
  const nu = Math.hypot(u[0], u[1], u[2]) || 1, ua = u[a] / nu, ub = u[b] / nu, un = u[axis] / nu;
  let h = z - c[axis];
  if (Laxis) h -= Laxis * Math.round(h / Laxis);
  const s = Math.hypot(ua, ub);
  // quick reject: half-width of the rod along the normal
  const hn = capped ? hl * Math.abs(un) + R : hl * Math.abs(un) + R * s;
  if (!(Math.abs(h) < hn)) return null;
  const e1 = s > 1e-12 ? [ua / s, ub / s] : [1, 0];
  return VIZ.rodSectionPolygon([c[a], c[b]], e1, s, un, h, hl, R, capped);
};

// A polygon as the renderer keeps it: centre of its bounding box, vertices relative to that centre (flat
// [dx0, dy0, dx1, dy1, ...]), the half-widths of the box (for the periodic images) and the enclosed area.
VIZ.polyRecord = function (poly) {
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity, area = 0;
  const n = poly.length;
  for (let k = 0; k < n; k++) {
    const p = poly[k], q = poly[(k + 1) % n];
    if (p[0] < x0) x0 = p[0]; if (p[0] > x1) x1 = p[0]; if (p[1] < y0) y0 = p[1]; if (p[1] > y1) y1 = p[1];
    area += p[0] * q[1] - q[0] * p[1];
  }
  const cx = 0.5 * (x0 + x1), cy = 0.5 * (y0 + y1), rel = new Float64Array(2 * n);
  for (let k = 0; k < n; k++) { rel[2 * k] = poly[k][0] - cx; rel[2 * k + 1] = poly[k][1] - cy; }
  return { cx, cy, rel, hx: 0.5 * (x1 - x0), hy: 0.5 * (y1 - y0), area: 0.5 * Math.abs(area) };
};

// 2D rods: the outline of every particle (a stadium or a rectangle). Adds to a 2D packing object
//   outline[i] = vertices relative to positions[i];  hext = [hx, hy] per particle (bounding half-widths).
VIZ.rodOutlines2D = function (data) {
  const N = data.N, p = data.positions, rod = data.rod, outline = new Array(N), hext = new Float32Array(2 * N);
  for (let i = 0; i < N; i++) {
    const ux = rod.axis[2 * i], uy = rod.axis[2 * i + 1], nu = Math.hypot(ux, uy) || 1;
    const poly = VIZ.rodSectionPolygon([0, 0], [ux / nu, uy / nu], 1, 0, 0, rod.halfLength[i], rod.radius[i], rod.capped);
    const rec = VIZ.polyRecord(poly || [[0, 0], [0, 0], [0, 0]]);       // a rod is symmetric about its centre: cx = cy = 0
    outline[i] = rec.rel; hext[2 * i] = rec.hx; hext[2 * i + 1] = rec.hy;
  }
  return { outline, hext };
};

// 3D rods: the cross-section by the plane x_axis = z as a 2D packing of polygons the 2D renderer draws.
// positions = centre of each cut's bounding box; diameters = diameter of the disc of the same area.
VIZ.computeSliceRod = function (data, axis, z) {
  const N = data.N, p = data.positions, dia = data.diameters, rod = data.rod;
  const inPlane = [0, 1, 2].filter(a => a !== axis), w = data.walls;
  const perAxis = w ? (w[axis] === 0) : true;            // no walls recorded -> periodic (as for spheres)
  const perIn = inPlane.map(dim => w ? (w[dim] === 0) : true);
  const Lin = [data.box[inPlane[0]], data.box[inPlane[1]]], Ln = data.box[axis];
  const recs = [], idx = [];
  for (let i = 0; i < N; i++) {
    const u = [rod.axis[i * 3], rod.axis[i * 3 + 1], rod.axis[i * 3 + 2]];
    for (let k = (perAxis ? -1 : 0); k <= (perAxis ? 1 : 0); k++) {   // the rod and its two images along the normal
      const c = [p[i * 3], p[i * 3 + 1], p[i * 3 + 2]]; c[axis] += k * Ln;
      const poly = VIZ.sliceRod(c, u, rod.halfLength[i], rod.radius[i], rod.capped, axis, z, 0);
      if (!poly) continue;
      const rec = VIZ.polyRecord(poly);
      if (!(rec.area > 0)) continue;
      if (perIn[0]) rec.cx -= Lin[0] * Math.floor(rec.cx / Lin[0]);   // centre back into the box
      if (perIn[1]) rec.cy -= Lin[1] * Math.floor(rec.cy / Lin[1]);
      recs.push(rec); idx.push(i);
    }
  }
  const n = recs.length, pos = new Float32Array(n * 2), diam = new Float32Array(n), diam3d = new Float32Array(n);
  const outline = new Array(n), hext = new Float32Array(n * 2);
  for (let i = 0; i < n; i++) {
    const r = recs[i];
    pos[i * 2] = r.cx; pos[i * 2 + 1] = r.cy; diam3d[i] = dia[idx[i]]; diam[i] = 2 * Math.sqrt(r.area / Math.PI);
    outline[i] = r.rel; hext[i * 2] = r.hx; hext[i * 2 + 1] = r.hy;
  }
  const out = { N: n, ndim: 2, box: Lin, phi: data.phi, positions: pos, diameters: diam, outline, hext, layers: {},
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
