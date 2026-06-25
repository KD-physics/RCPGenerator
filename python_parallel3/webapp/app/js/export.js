// export.js — export the current view to a PNG at a chosen DPI (re-renders the deck
// canvas at higher backing resolution). Returns a data URL; downloads if asked.
window.VIZ = window.VIZ || {};
VIZ.exportPNG = async function (renderer, dpi, download) {
  dpi = dpi || 150; const scale = dpi / 96;
  const dk = renderer.deck;
  dk.setProps({ useDevicePixels: scale }); dk.redraw && dk.redraw('export');
  await new Promise(r => setTimeout(r, 350));
  const canvas = dk.canvas || document.querySelector('#deck canvas');
  const url = canvas.toDataURL('image/png');
  dk.setProps({ useDevicePixels: true });
  if (download) {
    const a = document.createElement('a'); a.href = url;
    a.download = `packing_${dpi}dpi.png`; a.click();
  }
  window.__lastExport = { url, dpi, w: canvas.width, h: canvas.height };
  return url;
};
