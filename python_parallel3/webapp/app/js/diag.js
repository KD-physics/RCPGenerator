// diag.js — diagnostic readout to window for the test harness. Reads the renderer's
// CURRENT data (the 2D packing, or the active 3D slice).
window.VIZ = window.VIZ || {};
VIZ.diag = {
  renderer: null, source: null, frames: 0,
  init(renderer, source) { this.renderer = renderer; this.source = source; this.update(); },
  countVisible() {
    const r = this.renderer; if (!r) return 0;
    const cur = r.currentData(); const vs = r.getViewState(); const pxPerWorld = Math.pow(2, vs.zoom);
    const cw = r.deck.width || 800, ch = r.deck.height || 800;
    const halfW = (cw / 2) / pxPerWorld, halfH = (ch / 2) / pxPerWorld;
    const cx = vs.target[0], cy = vs.target[1]; const p = cur.positions, d = cur.ndim; let n = 0;
    const nreal = cur._nreal != null ? cur._nreal : cur.N;   // exclude periodic ghosts
    for (let i = 0; i < nreal; i++) {
      if (Math.abs(p[i * d] - cx) <= halfW && Math.abs(p[i * d + 1] - cy) <= halfH) n++;
    }
    return n;
  },
  update() {
    const r = this.renderer; if (!r) return;
    const cur = r.currentData(); const vs = r.getViewState(); this.frames++;
    window.__diag = {
      ready: true, rendered: !!window.__rendered,
      N: (cur._nreal != null ? cur._nreal : cur.N), ndim: cur.ndim, box: cur.box, phi: cur.phi,
      sourceNdim: this.source ? this.source.ndim : cur.ndim,
      slice: cur._slice || null,
      zoom: vs.zoom, target: vs.target, visible: this.countVisible(),
      activeOverlays: r.state.activeOverlays.slice(), colorBy: r.state.colorBy, frames: this.frames
    };
  }
};
