// export.js — export the current view to a PNG at a chosen DPI (re-renders the deck
// canvas at higher backing resolution). Returns a data URL; downloads if asked.
//
// The capture MUST happen in the same frame a deck render completes: deck.gl draws
// asynchronously, so reading the canvas on a fixed timer (the old approach) grabbed
// an empty (black) drawing buffer and left the live view un-repainted until a pan.
// We instead arm a one-shot hook that render2d's onAfterRender fires once the
// high-res frame is on screen, capture there (preserveDrawingBuffer is on), then
// restore resolution and force a redraw so the live view repaints immediately.
window.VIZ = window.VIZ || {};
VIZ.exportPNG = function (renderer, dpi, download) {
  return new Promise((resolve, reject) => {
    dpi = dpi || 150;
    const scale = dpi / 96;
    const dk = renderer.deck;
    const canvas = dk.canvas || document.querySelector('#deck canvas');
    let done = false;

    const finish = (url) => {
      if (done) return;
      done = true;
      VIZ.__exportHook = null;                 // disarm the one-shot
      dk.setProps({ useDevicePixels: true });  // back to screen resolution
      if (dk.redraw) dk.redraw('export-restore'); // repaint live view (fixes black-until-pan)
      window.__lastExport = { url, dpi, w: canvas.width, h: canvas.height };
      if (url) {
        if (download) {
          const a = document.createElement('a');
          a.href = url; a.download = `packing_${dpi}dpi.png`; a.click();
        }
        resolve(url);
      } else {
        reject(new Error('exportPNG: capture returned an empty image'));
      }
    };

    // Fires after the next render completes (at the export resolution).
    VIZ.__exportHook = () => {
      let url = null;
      try { url = canvas.toDataURL('image/png'); } catch (e) { /* tainted/oversize */ }
      finish(url);
    };

    // Bump backing resolution and request the high-res frame.
    dk.setProps({ useDevicePixels: scale });
    if (dk.redraw) dk.redraw('export');

    // Safety net: if no render fires (e.g. nothing changed), capture anyway.
    setTimeout(() => {
      if (done) return;
      let url = null;
      try { url = canvas.toDataURL('image/png'); } catch (e) {}
      finish(url);
    }, 2000);
  });
};
