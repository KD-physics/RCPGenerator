// color.js — color-by-scalar (e.g. Voronoi local-phi or diameter). Computes a
// per-particle RGB array from a layer's values via a colormap, sets it on the
// renderer's current data, and refreshes. Reverting (layer=null) restores solid fill.
window.VIZ = window.VIZ || {};

VIZ.scalarFor = function (cur, name) {
  if (name === 'diameter') return cur.diameters;
  const L = cur.layers && cur.layers[name];
  return (L && L.type === 'scalar_per_particle') ? L.values : null;
};

VIZ.applyColorBy = function (renderer, name, cmapName, domain) {
  const cur = renderer.currentData();
  if (!name) { cur.colors = null; renderer.state.colorBy = null; renderer.refresh();
               return { domain: null }; }
  // render-class palette: 'p<k>' -> random assignment of palette[k]'s colors per
  // particle (ghost-consistent via _srcIndex). Matches rcpgenerator/render.py.
  if (name[0] === 'p' && VIZ.renderPalettes && VIZ.renderPalettes[name.slice(1)]) {
    const pal = VIZ.renderPalettes[name.slice(1)], src = cur._srcIndex;
    const cols = new Uint8Array(cur.N * 3);
    for (let i = 0; i < cur.N; i++) {
      const key = (src ? src[i] : i) >>> 0;
      const c = pal[((key * 2654435761) >>> 0) % pal.length];
      cols[i * 3] = c[0]; cols[i * 3 + 1] = c[1]; cols[i * 3 + 2] = c[2];
    }
    cur.colors = cols; renderer.state.colorBy = name; renderer.state.cmap = 'palette';
    renderer.refresh(); return { domain: null };
  }
  const vals = VIZ.scalarFor(cur, name);
  if (!vals) throw new Error('no scalar layer: ' + name);
  if (cmapName === 'random') {                    // random color per distinct value (ghost-consistent)
    const h = (v, k) => { const s = Math.sin((v + k) * 99.713) * 10000.0; return s - Math.floor(s); };
    const cols = new Uint8Array(cur.N * 3);
    for (let i = 0; i < cur.N; i++) { const v = vals[i];
      cols[i*3] = 55 + h(v,0.1)*200; cols[i*3+1] = 55 + h(v,1.7)*200; cols[i*3+2] = 55 + h(v,3.3)*200; }
    cur.colors = cols; renderer.state.colorBy = name; renderer.state.cmap = 'random';
    renderer.refresh(); return { domain: null };
  }
  const cmap = VIZ.colormaps[cmapName] || VIZ.colormaps.viridis;
  let lo, hi;
  if (domain) { [lo, hi] = domain; }
  else { lo = Infinity; hi = -Infinity; const n = cur._nreal != null ? cur._nreal : cur.N;
         for (let i = 0; i < n; i++) { const v = vals[i]; if (v < lo) lo = v; if (v > hi) hi = v; } }
  const span = (hi - lo) || 1;
  const cols = new Uint8Array(cur.N * 3);
  for (let i = 0; i < cur.N; i++) {
    const c = cmap((vals[i] - lo) / span);
    cols[i * 3] = c[0]; cols[i * 3 + 1] = c[1]; cols[i * 3 + 2] = c[2];
  }
  cur.colors = cols; renderer.state.colorBy = name; renderer.state.cmap = cmapName;
  renderer.state.domain = [lo, hi]; renderer.refresh();
  return { domain: [lo, hi] };
};
