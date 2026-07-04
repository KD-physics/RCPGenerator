// overlay_pores.js — segmented-pore map overlay (2D). Draws the porenet pore cells
// (Delaunay triangles of the scaffold) tinted by per-pore fill eta. Registered in the
// overlay registry; toggled via state.activeOverlays. Self-contained (modular add-on).
window.VIZ = window.VIZ || {};
VIZ.overlays = VIZ.overlays || {};
VIZ.overlays.pores = function (cur, state) {
  const L = cur.layers && cur.layers.pores;
  if (!L || L.type !== 'polygons' || !window.deck.SolidPolygonLayer) return null;
  const J = L.json, cmap = VIZ.colormaps[state.cmap] || VIZ.colormaps.viridis;
  let lo = Infinity, hi = -Infinity;
  for (const e of J.eta) { if (e < lo) lo = e; if (e > hi) hi = e; }
  const sp = (hi - lo) || 1;
  const data = J.tris.map((t, i) => ({ polygon: [[t[0], t[1]], [t[2], t[3]], [t[4], t[5]]], eta: J.eta[i] }));
  return new window.deck.SolidPolygonLayer({
    id: 'pores-overlay', data, getPolygon: d => d.polygon,
    getFillColor: d => { const c = cmap((d.eta - lo) / sp); return [c[0], c[1], c[2], 140]; },
    stroked: true, getLineColor: [15, 15, 22, 220], getLineWidth: 0.3, lineWidthUnits: 'common',
    parameters: { depthTest: false }
  });
};
VIZ.poreStats = function (cur) {
  const L = cur.layers && cur.layers.pores; if (!L) return null;
  return { n_bodies: L.json.n_bodies, n_cells: L.json.tris.length };
};
