// stats.js — live field-of-view statistics over the VISIBLE (real) particles:
// P(local-phi), size distribution, running phi of the view. Recomputed on pan/zoom.
// Exposes window.__stats for the test harness; draws lightweight canvas plots.
window.VIZ = window.VIZ || {};
VIZ.stats = {
  renderer: null, canvas: null, _raf: 0,
  init(renderer, canvas) { this.renderer = renderer; this.canvas = canvas; },
  request() { if (this._raf) return; this._raf = requestAnimationFrame(() => { this._raf = 0; this.draw(); }); },
  visibleIndices() {
    const r = this.renderer, cur = r.currentData(), vs = r.getViewState();
    const ppw = Math.pow(2, vs.zoom), cw = r.deck.width || 800, ch = r.deck.height || 800;
    const hw = cw / 2 / ppw, hh = ch / 2 / ppw, cx = vs.target[0], cy = vs.target[1];
    const n = cur._nreal != null ? cur._nreal : cur.N, p = cur.positions, d = cur.ndim, idx = [];
    for (let i = 0; i < n; i++)
      if (Math.abs(p[i * d] - cx) <= hw && Math.abs(p[i * d + 1] - cy) <= hh) idx.push(i);
    return idx;
  },
  hist(vals, lo, hi, nb) {
    const c = new Array(nb).fill(0), sp = (hi - lo) || 1;
    for (const v of vals) { let b = Math.floor((v - lo) / sp * nb); b = Math.max(0, Math.min(nb - 1, b)); c[b]++; }
    return { lo, hi, nb, counts: c };
  },
  compute() {
    const cur = this.renderer.currentData(), idx = this.visibleIndices();
    const lp = cur.layers && cur.layers.local_phi ? cur.layers.local_phi.values : null;
    const dia = cur.diameters; let sphi = 0; const phivals = [], dvals = [];
    let dmin = Infinity, dmax = -Infinity;
    for (const i of idx) {
      if (lp) { phivals.push(lp[i]); sphi += lp[i]; }
      const dd = dia[i]; dvals.push(dd); if (dd < dmin) dmin = dd; if (dd > dmax) dmax = dd;
    }
    const out = { n: idx.length,
      mean_local_phi: lp && phivals.length ? sphi / phivals.length : null,
      phi_hist: lp ? this.hist(phivals, 0, 1, 30) : null,
      diam_hist: idx.length ? this.hist(dvals, dmin, dmax, 30) : null };
    window.__stats = out; return out;
  },
  draw() {
    const s = this.compute(), cv = this.canvas; if (!cv) return;
    const ctx = cv.getContext('2d'), W = cv.width, H = cv.height;
    ctx.clearRect(0, 0, W, H); ctx.fillStyle = '#1a1a20'; ctx.fillRect(0, 0, W, H);
    ctx.fillStyle = '#bbb'; ctx.font = '11px system-ui';
    ctx.fillText(`view: ${s.n} particles` + (s.mean_local_phi != null ? `   <local phi>=${s.mean_local_phi.toFixed(3)}` : ''), 8, 16);
    const bar = (hist, y0, h, label, col) => {
      if (!hist) return; ctx.fillStyle = '#999'; ctx.fillText(label, 8, y0 - 4);
      const mx = Math.max(...hist.counts) || 1, w = (W - 16) / hist.nb;
      ctx.fillStyle = col;
      hist.counts.forEach((c, i) => { const bh = (c / mx) * h; ctx.fillRect(8 + i * w, y0 + h - bh, Math.max(1, w - 1), bh); });
    };
    bar(s.phi_hist, 36, 90, 'P(local phi)  [0..1]', '#5a96e6');
    bar(s.diam_hist, 150, 90, 'P(diameter)  [min..max]', '#7ad07a');
  }
};
