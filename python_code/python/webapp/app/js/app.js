// app.js — orchestration. Parse ?bundle= (fetch) or accept drag-drop / test hook,
// build the renderer, wire diagnostics. UI/controls layer is added separately.
window.VIZ = window.VIZ || {};
(function () {
  function start(data) {
    const container = document.getElementById('deck');
    let renderData = data, source = null, axis = 2, z = 0;
    if (data.ndim === 3) {                       // 3D -> render a cross-section slice
      source = data; axis = 2; z = data.box[2] / 2;
      renderData = VIZ.computeSlice(data, axis, z);
    } else {
      renderData = VIZ.addGhosts(data);          // 2D -> periodic edge wrap
    }
    const r = VIZ.Renderer2D(container, renderData);
    window.__viz = r; r.source = source;
    VIZ.diag.init(r, source);
    VIZ.onView = () => VIZ.diag.update();
    if (source) {
      r.setSlice = (zz, ax) => { if (ax != null) axis = ax; z = zz;
        r.setData(VIZ.computeSlice(source, axis, z)); VIZ.diag.update(); };
      r.sliceInfo = () => ({ axis, z, box: source.box });
    }
    // UI panel on by default (users); correctness tests pass ?ui=0 (full-canvas).
    if (VIZ.buildUI && new URLSearchParams(location.search).get('ui') !== '0') {
      document.getElementById('deck').style.right = '300px';   // leave room for panel
      VIZ.buildUI(r, data, source);
    }
    ['hint', 'loaders'].forEach(id => { const e = document.getElementById(id); if (e) e.style.display = 'none'; });
    return r;
  }
  VIZ.start = start;

  async function boot() {
    if (window.__EMBEDDED_BUNDLE) { try { start(VIZ.loadEmbedded(window.__EMBEDDED_BUNDLE)); } catch (e) { console.error(e); window.__loadError = String(e); } return; }
    const params = new URLSearchParams(location.search);
    const bundle = params.get('bundle');
    if (bundle) {
      try { start(await VIZ.loadFromUrl(bundle)); }
      catch (e) { console.error('load failed', e); window.__loadError = String(e); }
    }
  }

  // drag-drop a folder/bundle: folders arrive as directory ENTRIES, not files —
  // traverse them via webkitGetAsEntry. Falls back to a flat file list.
  async function collectFiles(entries, flat) {
    if (entries && entries.length) {
      const out = [];
      const walk = async (entry) => {
        if (entry.isFile) { await new Promise(res => entry.file(f => { out.push(f); res(); }, res)); }
        else if (entry.isDirectory) {
          const rd = entry.createReader();
          let batch; do { batch = await new Promise(res => rd.readEntries(res, () => res([])));
            for (const e of batch) await walk(e); } while (batch.length);
        }
      };
      for (const e of entries) await walk(e);
      if (out.length) return out;
    }
    return flat || [];
  }
  // wrapper kept for the test harness (captures entries+flat from a DataTransfer)
  async function filesFromDrop(dt) {
    if (!dt) return [];
    const entries = dt.items ? Array.from(dt.items).map(it => it.webkitGetAsEntry && it.webkitGetAsEntry()).filter(Boolean) : [];
    return collectFiles(entries, Array.from(dt.files || []));
  }
  // Listen on window in the CAPTURE phase so the deck canvas can't swallow the drop.
  // Capture the directory ENTRIES synchronously (they are neutered after any await).
  const onDragOver = e => { e.preventDefault(); e.stopPropagation(); if (e.dataTransfer) e.dataTransfer.dropEffect = 'copy'; };
  window.addEventListener('dragover', onDragOver, true);
  window.addEventListener('dragenter', onDragOver, true);
  window.addEventListener('drop', async e => {
    e.preventDefault(); e.stopPropagation();
    const dt = e.dataTransfer;
    const entries = dt && dt.items ? Array.from(dt.items).map(it => it.webkitGetAsEntry && it.webkitGetAsEntry()).filter(Boolean) : [];
    const flat = dt ? Array.from(dt.files || []) : [];
    try {
      const files = await collectFiles(entries, flat);
      if (files.length) start(await VIZ.loadFromFiles(files));
      else { window.__loadError = 'no files in drop'; alert('No bundle files found in the drop. Use the “Choose bundle folder” button.'); }
    } catch (err) { console.error(err); window.__loadError = String(err); alert('Could not load bundle: ' + err.message); }
  }, true);
  window.__filesFromDrop = filesFromDrop;   // exposed for the real-extraction audit test
  // file / directory picker buttons (robust fallback for fragile OS drag-drop)
  function wirePicker(id) {
    const inp = document.getElementById(id); if (!inp) return;
    inp.addEventListener('change', async e => {
      if (!e.target.files || !e.target.files.length) return;
      try { start(await VIZ.loadFromFiles(e.target.files)); }
      catch (err) { console.error(err); window.__loadError = String(err); alert('Load failed: ' + err.message); }
    });
  }
  wirePicker('dirpick'); wirePicker('filepick');

  // test hook: load from a FileList-like array
  window.__loadFiles = async (files) => start(await VIZ.loadFromFiles(files));

  if (document.readyState !== 'loading') boot(); else document.addEventListener('DOMContentLoaded', boot);
})();
