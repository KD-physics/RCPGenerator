// render2d.js — renderer layer. deck.gl OrthographicView + ScatterplotLayer (binary
// attributes, crisp SDF circles in world units). Reused for 3D cross-sections: the
// caller swaps in a 2D slice via setData(). Overlays add layers via VIZ.overlays.
window.VIZ = window.VIZ || {};
VIZ.overlays = VIZ.overlays || {};   // name -> function(data,state,viewState)->deck.Layer

VIZ.Renderer2D = function (container, data0, opts) {
  opts = opts || {};
  const D = window.deck;
  const w = container.clientWidth || 800, h = container.clientHeight || 800;
  let cur = data0;
  let radii = halfRadii(cur);
  const [Lx, Ly] = cur.box;
  const fullFit = Math.log2(Math.min(w / Lx, h / Ly));
  const fitZoom = fullFit - 0.07;            // ~5% margin so boundary box + wrapped ghosts are visible
  let viewState = { target: [Lx / 2, Ly / 2, 0], zoom: fitZoom,
                    minZoom: fullFit - 4, maxZoom: fullFit + 34 };
  const state = { fillColor: opts.fillColor || [90, 150, 230], colorBy: null, activeOverlays: [],
                  outline: 1.0, outlineColor: [205, 205, 212, 255],
                  showBoundary: true, boundaryColor: [240, 240, 250, 255] };

  function halfRadii(d) { const r = new Float32Array(d.N); for (let i = 0; i < d.N; i++) r[i] = d.diameters[i] * 0.5; return r; }

  // Ellipses (cur.shape = [sxx, sxy, syy] per particle): each is drawn as a filled polygon of ELL_SEG points
  // on its exact outline, centre + L (cos t, sin t) with L the Cholesky factor of the shape matrix. Discs
  // (no cur.shape) keep the ScatterplotLayer below, unchanged.
  const ELL_SEG = 96;
  let ellCache = null, ellFor = null;
  function ellipsePolygons() {
    if (ellFor === cur && ellCache) return ellCache;
    const n = cur.N, sh = cur.shape, p = cur.positions, cs = new Float64Array(ELL_SEG), sn = new Float64Array(ELL_SEG);
    for (let k = 0; k < ELL_SEG; k++) { cs[k] = Math.cos(2 * Math.PI * k / ELL_SEG); sn[k] = Math.sin(2 * Math.PI * k / ELL_SEG); }
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const l11 = Math.sqrt(sh[i * 3]), l21 = l11 > 0 ? sh[i * 3 + 1] / l11 : 0;
      const l22 = Math.sqrt(Math.max(0, sh[i * 3 + 2] - l21 * l21)), x = p[i * 2], y = p[i * 2 + 1];
      const poly = new Array(ELL_SEG + 1);
      for (let k = 0; k < ELL_SEG; k++) poly[k] = [x + l11 * cs[k], y + l21 * cs[k] + l22 * sn[k]];
      poly[ELL_SEG] = poly[0];
      out[i] = { polygon: poly, i };
    }
    ellCache = out; ellFor = cur;
    return out;
  }
  // Polygons (cur.outline[i] = vertices relative to the particle's position; rods and their cuts, rod.js):
  // drawn by the same PolygonLayer as the ellipses.
  let outCache = null, outFor = null;
  function outlinePolygons() {
    if (outFor === cur && outCache) return outCache;
    const n = cur.N, p = cur.positions, out = new Array(n);
    for (let i = 0; i < n; i++) {
      const rel = cur.outline[i], m = rel.length / 2, x = p[i * 2], y = p[i * 2 + 1], poly = new Array(m + 1);
      for (let k = 0; k < m; k++) poly[k] = [x + rel[2 * k], y + rel[2 * k + 1]];
      poly[m] = poly[0];
      out[i] = { polygon: poly, i };
    }
    outCache = out; outFor = cur;
    return out;
  }
  function ellipseLayer(polygons) {
    const cols = (state.colorBy && cur.colors) ? cur.colors : null;
    return new D.PolygonLayer({
      id: 'particles', data: polygons || ellipsePolygons(), getPolygon: d => d.polygon, parameters: { depthTest: false },
      filled: true, getFillColor: cols ? (d => [cols[d.i * 3], cols[d.i * 3 + 1], cols[d.i * 3 + 2]]) : state.fillColor,
      stroked: state.outline > 0, getLineColor: state.outlineColor, lineWidthUnits: 'pixels',
      getLineWidth: state.outline, lineWidthMinPixels: 0, lineWidthMaxPixels: 3,
      updateTriggers: { getFillColor: [cols, state.fillColor], getLineColor: [state.outlineColor], getLineWidth: [state.outline] } });
  }

  function particleLayer() {
    if (cur.outline) return ellipseLayer(outlinePolygons());
    if (cur.shape) return ellipseLayer();
    const attrs = { getPosition: { value: cur.positions, size: cur.ndim },
                    getRadius: { value: radii, size: 1 } };
    const props = { id: 'particles', data: { length: cur.N, attributes: attrs },
                    radiusUnits: 'common', antialiasing: true, parameters: { depthTest: false },
                    getFillColor: state.fillColor,
                    // thin light-grey outline (visible on the dark bg), constant in SCREEN pixels
                    stroked: state.outline > 0, getLineColor: state.outlineColor, lineWidthUnits: 'pixels',
                    getLineWidth: state.outline, lineWidthMinPixels: 0, lineWidthMaxPixels: 3 };
    if (state.colorBy && cur.colors) { attrs.getFillColor = { value: cur.colors, size: 3 }; delete props.getFillColor; }
    return new D.ScatterplotLayer(props);
  }
  // Boundary: solid = hard wall, dashed = periodic; a solid circle is the cross-section
  // of a disk/cylinder/sphere container. cur.boundary (from the bundle's walls, built in
  // loader/slice) is a list of {path, dashed} segments. Absent -> legacy dashed box.
  function pathLayer(id, paths, dashed) {
    if (!paths.length) return null;
    const props = { id, data: paths.map(p => ({ path: p })), getPath: d => d.path,
                    getColor: state.boundaryColor, widthUnits: 'pixels', getWidth: 2,
                    widthMinPixels: 1.5, parameters: { depthTest: false } };
    if (dashed && D.PathStyleExtension) { props.extensions = [new D.PathStyleExtension({ dash: true })];
      props.getDashArray = [7, 6]; props.dashJustified = true; }
    return new D.PathLayer(props);
  }
  function boundaryLayers() {
    if (cur.boundary && cur.boundary.segments) {
      const dash = cur.boundary.segments.filter(s => s.dashed).map(s => s.path);
      const solid = cur.boundary.segments.filter(s => !s.dashed).map(s => s.path);
      return [pathLayer('boundary-dash', dash, true), pathLayer('boundary-solid', solid, false)].filter(Boolean);
    }
    const [Lx, Ly] = cur.box;   // legacy default: dashed periodic box
    return [pathLayer('boundary', [[[0, 0], [Lx, 0], [Lx, Ly], [0, Ly], [0, 0]]], true)].filter(Boolean);
  }
  function buildLayers() {
    const layers = [particleLayer()];
    for (const name of state.activeOverlays)
      if (VIZ.overlays[name]) { const l = VIZ.overlays[name](cur, state, viewState); if (l) layers.push(l); }
    if (state.showBoundary !== false) layers.push(...boundaryLayers());
    return layers;
  }

  const dk = new D.Deck({
    parent: container,
    views: [new D.OrthographicView({ flipY: false })],
    controller: true,
    initialViewState: viewState,
    layers: buildLayers(),
    glOptions: { preserveDrawingBuffer: true },
    onViewStateChange: ({ viewState: vs }) => {
      viewState = Object.assign({}, viewState, vs); dk.setProps({ viewState }); if (VIZ.onView) VIZ.onView(viewState);
    },
    onAfterRender: () => { window.__rendered = true; if (VIZ.diag) VIZ.diag.update();
      if (VIZ.stats && VIZ.stats.renderer) VIZ.stats.request();
      if (VIZ.__exportHook) VIZ.__exportHook(); }   // one-shot PNG capture (export.js)
  });

  const api = {
    deck: dk, box: cur.box, state,
    currentData: () => cur,
    getViewState: () => viewState,
    setViewState: (vs) => { viewState = Object.assign({}, viewState, vs); dk.setProps({ viewState }); },
    refresh: () => dk.setProps({ layers: buildLayers() }),
    setData: (nd) => { cur = nd; radii = halfRadii(cur); api.refresh(); },
    fit: () => api.setViewState({ target: [cur.box[0] / 2, cur.box[1] / 2, 0],
                                  zoom: (cur.box[0] === Lx && cur.box[1] === Ly) ? fitZoom : Math.log2(Math.min(w / cur.box[0], h / cur.box[1])) - 0.07 }),
    fitTight: () => api.setViewState({ target: [cur.box[0] / 2, cur.box[1] / 2, 0], zoom: fullFit })
  };
  return api;
};
