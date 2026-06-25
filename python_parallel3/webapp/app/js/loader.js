// loader.js — data layer. Reads a bundle (manifest + raw float32 arrays) from a URL
// (fetch) or from dropped Files. Returns a normalized packing object. Modular: new
// layer types are read generically from the manifest.
window.VIZ = window.VIZ || {};

VIZ.normalize = function (man, buffers) {
  const N = man.N, d = man.ndim;
  const out = { N, ndim: d, box: man.box, phi: man.phi,
                positions: new Float32Array(buffers[man.positions.file]),
                diameters: new Float32Array(buffers[man.diameters.file]),
                layers: {} };
  for (const L of (man.layers || [])) {
    out.layers[L.name] = (L.type === 'scalar_per_particle')
      ? { type: L.type, values: new Float32Array(buffers[L.file]) }
      : { type: L.type, json: buffers[L.file] };
  }
  return out;
};

// Embedded bundle (CLI / single-file standalone): { manifest, files:{name:base64} }.
VIZ.loadEmbedded = function (bundle) {
  const buffers = {};
  for (const name in bundle.files) {
    const bin = atob(bundle.files[name]);
    if (name.endsWith('.json')) { buffers[name] = JSON.parse(bin); }
    else { const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); buffers[name] = u.buffer; }
  }
  return VIZ.normalize(bundle.manifest, buffers);
};

VIZ.loadFromUrl = async function (baseUrl) {
  const man = await (await fetch(baseUrl + '/manifest.json')).json();
  const need = [man.positions.file, man.diameters.file,
                ...(man.layers || []).map(L => L.file)];
  const buffers = {};
  for (const f of need) {
    const r = await fetch(baseUrl + '/' + f);
    buffers[f] = f.endsWith('.json') ? await r.json() : await r.arrayBuffer();
  }
  return VIZ.normalize(man, buffers);
};

// Drag-drop: a FileList (single manifest+files folder, or a lone .f32/.npz is not
// supported in-browser — folder/zip with manifest is the contract).
VIZ.loadFromFiles = async function (fileList) {
  const files = {}; let manFile = null;
  for (const f of Array.from(fileList)) {
    const name = (f.webkitRelativePath || f.name).split('/').pop();
    files[name] = f; if (name === 'manifest.json') manFile = f;
  }
  if (!manFile) throw new Error('bundle must contain manifest.json');
  const man = JSON.parse(await manFile.text());
  const buffers = {};
  const need = [man.positions.file, man.diameters.file, ...(man.layers || []).map(L => L.file)];
  for (const f of need) {
    if (!files[f]) throw new Error('missing ' + f);
    buffers[f] = f.endsWith('.json') ? JSON.parse(await files[f].text())
                                     : await files[f].arrayBuffer();
  }
  return VIZ.normalize(man, buffers);
};
