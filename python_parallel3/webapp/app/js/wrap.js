// wrap.js — periodic boundary handling. A particle whose circle crosses a box edge
// also appears on the opposite side (periodic image). We append ghost copies for
// edge-crossing particles so the render is correct (no edge gaps) and area-fraction
// recovers phi. Ghosts replicate scalar layer values; `_nreal` marks the real count.
window.VIZ = window.VIZ || {};

VIZ.addGhosts = function (data) {
  const [Lx, Ly] = data.box, N = data.N, p = data.positions, dia = data.diameters;
  const gx = [], gy = [], gd = [], gi = [];           // ghost x,y,dia, source index
  for (let i = 0; i < N; i++) {
    const x = p[i * 2], y = p[i * 2 + 1], r = dia[i] * 0.5;
    const sx = x < r ? Lx : (x > Lx - r ? -Lx : 0);
    const sy = y < r ? Ly : (y > Ly - r ? -Ly : 0);
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
  return Object.assign({}, data, { N: N + M, positions: pos, diameters: diam, layers, _nreal: N, _srcIndex: src });
};
