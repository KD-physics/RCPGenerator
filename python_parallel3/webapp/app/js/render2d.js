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

  function particleLayer() {
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
  function boundaryLayer() {                       // periodic box outline (dashed if available)
    const [Lx, Ly] = cur.box;
    const props = { id: 'boundary', data: [{ path: [[0, 0], [Lx, 0], [Lx, Ly], [0, Ly], [0, 0]] }],
                    getPath: d => d.path, getColor: state.boundaryColor, widthUnits: 'pixels',
                    getWidth: 2, widthMinPixels: 1.5, parameters: { depthTest: false } };
    if (D.PathStyleExtension) { props.extensions = [new D.PathStyleExtension({ dash: true })];
      props.getDashArray = [7, 6]; props.dashJustified = true; }
    return new D.PathLayer(props);
  }
  function buildLayers() {
    const layers = [particleLayer()];
    for (const name of state.activeOverlays)
      if (VIZ.overlays[name]) { const l = VIZ.overlays[name](cur, state, viewState); if (l) layers.push(l); }
    if (state.showBoundary !== false) layers.push(boundaryLayer());
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
      if (VIZ.stats && VIZ.stats.renderer) VIZ.stats.request(); }
  });

  const api = {
    deck: dk, box: cur.box, state,
    currentData: () => cur,
    getViewState: () => viewState,
    setViewState: (vs) => { viewState = Object.assign({}, viewState, vs); dk.setProps({ viewState }); },
    refresh: () => dk.setProps({ layers: buildLayers() }),
    setData: (nd) => { cur = nd; radii = halfRadii(cur); api.refresh(); },
    fit: () => api.setViewState({ target: [cur.box[0] / 2, cur.box[1] / 2, 0], zoom: fitZoom }),
    fitTight: () => api.setViewState({ target: [cur.box[0] / 2, cur.box[1] / 2, 0], zoom: fullFit })
  };
  return api;
};
