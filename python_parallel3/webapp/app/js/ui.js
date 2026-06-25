// ui.js — controls panel (right side) + live stats. Decoupled from the renderer:
// reads available layers, wires color-by / palette / overlay toggles / 3D slice /
// PNG export, and redraws the FOV stats on view change. VIZ.buildUI(renderer,data,source).
window.VIZ = window.VIZ || {};

// set the render-area background; flip outline/boundary to a contrasting ink so they
// stay visible (light bg -> dark ink, dark bg -> light ink).
VIZ.setBackground = function (renderer, color) {
  const deckEl = document.getElementById('deck'); if (deckEl) deckEl.style.background = color;
  let r = 13, g = 13, b = 16; const m = /^#?([0-9a-fA-F]{6})$/.exec(color);
  if (m) { const h = parseInt(m[1], 16); r = (h >> 16) & 255; g = (h >> 8) & 255; b = h & 255; }
  const lum = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
  renderer.state.outlineColor = lum > 0.55 ? [40, 40, 48, 255] : [205, 205, 212, 255];
  renderer.state.boundaryColor = lum > 0.55 ? [70, 70, 80, 255] : [240, 240, 250, 255];
  renderer.refresh(); window.__bg = color;
};

VIZ.BACKGROUNDS = { black: '#0d0d10', charcoal: '#22232a', slate: '#2b3a4a',
                    navy: '#10182e', 'light grey': '#e6e6ea', white: '#ffffff', beige: '#efe9dd' };

VIZ.buildUI = function (renderer, data, source) {
  const el = (t, p, kids) => { const e = document.createElement(t); Object.assign(e, p || {}); (kids || []).forEach(k => e.appendChild(k)); return e; };
  const panel = el('div', { id: 'panel' });
  panel.style.cssText = 'position:absolute;top:0;right:0;width:300px;height:100%;background:#15151b;' +
    'color:#ddd;font:12px system-ui;padding:10px;box-sizing:border-box;overflow:auto;border-left:1px solid #333';
  const row = (label, ctrl) => el('div', { style: 'margin:6px 0' }, [el('label', { textContent: label, style: 'display:block;color:#9aa;margin-bottom:2px' }), ctrl]);

  // color-by + palette
  const scalarLayers = ['(solid)', 'diameter', ...Object.keys(data.layers || {}).filter(k => (data.layers[k] || {}).type === 'scalar_per_particle')];
  const colorSel = el('select', { id: 'colorBy' }); scalarLayers.forEach(n => colorSel.appendChild(el('option', { value: n, textContent: n })));
  // render-class palettes (random per-particle assignment) appended to the same menu
  Object.keys(VIZ.renderPalettes || {}).forEach(k => colorSel.appendChild(
    el('option', { value: 'p' + k, textContent: 'palette ' + k + ' (' + VIZ.renderPalettes[k].length + ' colors)' })));
  const palSel = el('select', { id: 'palette' });
  [...Object.keys(VIZ.colormaps), 'random'].forEach(n => palSel.appendChild(
    el('option', { value: n, textContent: n + (VIZ.COLORBLIND_SAFE.includes(n) ? ' (cb-safe)' : '') })));
  function applyColor() {
    const n = colorSel.value === '(solid)' ? null : colorSel.value;
    VIZ.applyColorBy(renderer, n, palSel.value); VIZ.stats.draw();
  }
  colorSel.onchange = applyColor; palSel.onchange = applyColor;

  // 3D slice slider
  let sliceRow = null;
  if (source) {
    const Lz = source.box[2];
    const sl = el('input', { id: 'sliceSlider', type: 'range', min: 0, max: 1000, value: 500, style: 'width:100%' });
    const lbl = el('span', { textContent: 'z = 0.50' });
    sl.oninput = () => { const z = (sl.value / 1000) * Lz; lbl.textContent = 'z = ' + (sl.value / 1000).toFixed(2);
      renderer.setSlice(z); if (renderer.state.colorBy) applyColor(); VIZ.stats.draw(); };
    sliceRow = row('3D cross-section', el('div', {}, [sl, lbl]));
  }

  // PNG export
  const dpi = el('input', { type: 'number', value: 150, min: 72, max: 600, style: 'width:70px' });
  const expBtn = el('button', { id: 'exportBtn', textContent: 'Export PNG', style: 'margin-left:8px' });
  expBtn.onclick = () => VIZ.exportPNG(renderer, parseInt(dpi.value), true);

  // stats canvas
  const statsCanvas = el('canvas', { id: 'statsCanvas', width: 280, height: 260, style: 'background:#1a1a20;width:280px;height:260px;margin-top:8px' });

  panel.appendChild(el('div', { textContent: 'RCP Packing Visualizer', style: 'font-weight:600;margin-bottom:8px' }));
  panel.appendChild(row('color by', colorSel));
  panel.appendChild(row('palette', palSel));
  // background colour of the render area (presets + custom picker)
  const bgSel = el('select', { id: 'bgSelect' });
  Object.keys(VIZ.BACKGROUNDS).forEach(n => bgSel.appendChild(el('option', { value: VIZ.BACKGROUNDS[n], textContent: n })));
  const bgPick = el('input', { id: 'bgPick', type: 'color', value: VIZ.BACKGROUNDS.black, style: 'margin-left:6px;vertical-align:middle' });
  bgSel.onchange = () => { bgPick.value = bgSel.value; VIZ.setBackground(renderer, bgSel.value); };
  bgPick.oninput = () => VIZ.setBackground(renderer, bgPick.value);
  panel.appendChild(row('background', el('div', {}, [bgSel, bgPick])));
  if (sliceRow) panel.appendChild(sliceRow);
  panel.appendChild(row('export', el('div', {}, [dpi, expBtn])));

  // overlay toggles (polygon layers, e.g. segmented-pore map)
  Object.keys(data.layers || {}).forEach(name => {
    if ((data.layers[name] || {}).type !== 'polygons') return;
    const cb = el('input', { id: 'ov_' + name, type: 'checkbox' });
    cb.onchange = () => { const a = renderer.state.activeOverlays, i = a.indexOf(name);
      if (cb.checked && i < 0) a.push(name); else if (!cb.checked && i >= 0) a.splice(i, 1);
      renderer.refresh(); VIZ.stats.draw(); };
    panel.appendChild(row('overlay: ' + name + (VIZ.poreStats ? ` (${(VIZ.poreStats(renderer.currentData()) || {}).n_bodies || '?'} pores)` : ''), cb));
  });
  panel.appendChild(el('div', { textContent: 'field-of-view stats', style: 'color:#9aa;margin-top:10px' }));
  panel.appendChild(statsCanvas);
  document.body.appendChild(panel);

  VIZ.stats.init(renderer, statsCanvas);
  VIZ.onView = () => { VIZ.diag.update(); VIZ.stats.request(); };
  // default to color-by diameter so palettes recolor immediately (not "(solid)")
  colorSel.value = 'diameter'; applyColor();
  VIZ.stats.request();
  window.__ui = { colorSel, palSel, dpi, expBtn, statsCanvas };
};
