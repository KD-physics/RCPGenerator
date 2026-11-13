// wrap.js — periodic boundary handling. A particle whose circle crosses a box edge
// also appears on the opposite side (periodic image). We append ghost copies for
// edge-crossing particles so the render is correct (no edge gaps) and area-fraction
// recovers phi. Ghosts replicate scalar layer values; `_nreal` marks the real count.
window.VIZ = window.VIZ || {};

VIZ.addGhosts = function (data) {
  const [Lx, Ly] = data.box, N = data.N, p = data.positions, dia = data.diameters;
  // Periodic images (ghosts) are created ONLY for PERIODIC dimensions. Hard walls and
  // curved (disk/cylinder/sphere) boundaries have no periodic image, so wrapping them
  // drew spurious particles OUTSIDE the container. Periodicity: data._periodic (set per
  // in-plane dim by the slicer), else derived from data.walls, else all-periodic (legacy).
  const per = data._periodic ? data._periodic
    : (data.walls ? [0, 1].map(d => { const w = data.walls, t = (w[0] < 0 ? -w[0] : 0); return d < t ? false : (w[d] === 0); })
                  : [true, true]);
  const px = per[0], py = per[1];
  const sh = data.shape || null;                      // ellipses: [sxx, sxy, syy] per particle
  const ol = data.outline || null, he = data.hext || null;   // polygons (rods): outline + bounding half-widths per particle
  const gx = [], gy = [], gd = [], gi = [];           // ghost x,y,dia, source index
  for (let i = 0; i < N; i++) {
    const x = p[i * 2], y = p[i * 2 + 1], r = dia[i] * 0.5;
    // half-extent along each axis: the radius for a disc, sqrt of the shape matrix's diagonal for an ellipse
    const rx = he ? he[i * 2] : (sh ? Math.sqrt(sh[i * 3]) : r), ry = he ? he[i * 2 + 1] : (sh ? Math.sqrt(sh[i * 3 + 2]) : r);
    const sx = px ? (x < rx ? Lx : (x > Lx - rx ? -Lx : 0)) : 0;
    const sy = py ? (y < ry ? Ly : (y > Ly - ry ? -Ly : 0)) : 0;
    if (sx) { gx.push(x + sx); gy.push(y); gd.push(dia[i]); gi.push(i); }
    if (sy) { gx.push(x); gy.push(y + sy); gd.push(dia[i]); gi.push(i); }
    if (sx && sy) { gx.push(x + sx); gy.push(y + sy); gd.push(dia[i]); gi.push(i); }
  }
  const M = gx.length; if (!M) { data._nreal = N; data._srcIndex = null; return data; }
  const pos = new Float32Array((N + M) * 2), diam = new Float32Array(N + M);
  const src = new Int32Array(N + M);                  // real index of each particle (ghosts -> source)
  pos.set(p.subarray(0, N * 2)); diam.set(dia.subarray(0, N));
  for (let i = 0; i < N; i++) src[i] = i;
  for (let j = 0; j < M; j++) { pos[(N + j) * 2] = gx[j]; pos[(N + j) * 2 + 1] = gy[j]; diam[N + j] = gd[j]; src[N + j] = gi[j]; }
  const layers = {};
  for (const name in (data.layers || {})) {
    const L = data.layers[name];
    if (L.type === 'scalar_per_particle') {
      const v = new Float32Array(N + M); v.set(L.values.subarray(0, N));
      for (let j = 0; j < M; j++) v[N + j] = L.values[gi[j]];
      layers[name] = { type: L.type, values: v };
    } else layers[name] = L;
  }
  const res = Object.assign({}, data, { N: N + M, positions: pos, diameters: diam, layers, _nreal: N, _srcIndex: src });
  if (sh) {                                           // ghosts carry their source's shape
    const shape = new Float32Array((N + M) * 3); shape.set(sh.subarray(0, N * 3));
    for (let j = 0; j < M; j++) for (let k = 0; k < 3; k++) shape[(N + j) * 3 + k] = sh[gi[j] * 3 + k];
    res.shape = shape;
  }
  if (ol) {                                           // ghosts carry their source's outline
    const outline = ol.slice(), hext = new Float32Array((N + M) * 2); hext.set(he.subarray(0, N * 2));
    for (let j = 0; j < M; j++) { outline.push(ol[gi[j]]); hext[(N + j) * 2] = he[gi[j] * 2]; hext[(N + j) * 2 + 1] = he[gi[j] * 2 + 1]; }
    res.outline = outline; res.hext = hext;
  }
  return res;
};
