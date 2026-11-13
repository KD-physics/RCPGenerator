// loader.js — data layer. Reads a bundle (manifest + raw float32 arrays) from a URL
// (fetch) or from dropped Files. Returns a normalized packing object. Modular: new
// layer types are read generically from the manifest.
window.VIZ = window.VIZ || {};

// --- boundary helpers (shared by native-2D load and 3D cross-sections) -------------
VIZ.circlePath = function (cx, cy, R, n) {
  n = n || 96; const path = [];
  for (let i = 0; i <= n; i++) { const t = i / n * 2 * Math.PI; path.push([cx + R * Math.cos(t), cy + R * Math.sin(t)]); }
  return path;
};
VIZ.boxEdges = function (Lx, Ly, xDashed, yDashed) {   // per-edge; dashed = periodic, solid = hard
  return [{ path: [[0, 0], [0, Ly]], dashed: !!xDashed }, { path: [[Lx, 0], [Lx, Ly]], dashed: !!xDashed },
          { path: [[0, 0], [Lx, 0]], dashed: !!yDashed }, { path: [[0, Ly], [Lx, Ly]], dashed: !!yDashed }];
};
VIZ.boundaryForBox = function (box, walls) {   // native 2D container shape from walls
  const t = walls[0] < 0 ? -walls[0] : 0;      // walls[0] <= -2 -> disk of diameter box[0]
  if (t >= 2) return { segments: [{ path: VIZ.circlePath(box[0] / 2, box[0] / 2, box[0] / 2), dashed: false }] };
  return { segments: VIZ.boxEdges(box[0], box[1], walls[0] === 0, walls[1] === 0) };
};

// --- ellipsoids (additive). A bundle may carry `aspect_ratios` (ndim-1 numbers: axis k has length D x ratio)
// and `orientations` (one row-major ndim x ndim rotation per particle, columns = body axes). Without them, or
// with every ratio equal to 1, nothing below runs and the packing is drawn as discs/spheres exactly as before.
VIZ.bundleFiles = function (man) {     // every data file a bundle needs
  const rods = man.rods ? [man.rods.length, man.rods.radius].filter(v => v && typeof v === 'object').map(v => v.file) : [];
  return [man.positions.file, man.diameters.file, ...(man.orientations ? [man.orientations.file] : []), ...rods,
          ...(man.layers || []).map(L => L.file)];
};
VIZ.isEllipsoid = function (man) {
  return !!(man.aspect_ratios && man.orientations && man.aspect_ratios.some(a => a !== 1));
};
// Shape matrices A = R diag(semi^2) R^T, semi-axis k = D/2 x (1, ratios)[k]. 2D -> [axx, axy, ayy] per
// particle (what the renderer draws); 3D -> [a00, a01, a02, a11, a12, a22] (what the slicer cuts).
VIZ.shapeMatrices = function (N, d, dia, ori, ratios) {
  const rel = [1, ...ratios], nu = d === 2 ? 3 : 6, out = new Float32Array(N * nu), s2 = new Array(d);
  for (let i = 0; i < N; i++) {
    for (let k = 0; k < d; k++) { const a = 0.5 * dia[i] * rel[k]; s2[k] = a * a; }
    let c = 0;
    for (let a = 0; a < d; a++) for (let b = a; b < d; b++) {
      let sum = 0; for (let k = 0; k < d; k++) sum += ori[i * d * d + a * d + k] * s2[k] * ori[i * d * d + b * d + k];
      out[i * nu + c++] = sum;
    }
  }
  return out;
};

VIZ.normalize = function (man, buffers) {
  const N = man.N, d = man.ndim;
  const out = { N, ndim: d, box: man.box, phi: man.phi,
                walls: (man.walls != null ? man.walls : null),   // per-dim BC; null -> legacy dashed box
                positions: new Float32Array(buffers[man.positions.file]),
                diameters: new Float32Array(buffers[man.diameters.file]),
                layers: {} };
  if (d === 2 && out.walls) out.boundary = VIZ.boundaryForBox(man.box, out.walls);
  if (VIZ.isEllipsoid(man) && (d === 2 || d === 3)) {
    out.ratios = man.aspect_ratios.slice();
    out.orientations = new Float32Array(buffers[man.orientations.file]);
    const A = VIZ.shapeMatrices(N, d, out.diameters, out.orientations, out.ratios);
    if (d === 2) out.shape = A; else out.A3 = A;
  }
  // rods (rod.js): manifest.rods = { capped, length, radius }, length and radius in units of each particle's D,
  // each one number or a file of N float32; the axis is column 0 of the particle's rotation matrix.
  if (man.rods && man.orientations && (d === 2 || d === 3)) {
    const perParticle = v => (v && typeof v === 'object') ? new Float32Array(buffers[v.file]) : null;
    const Lp = perParticle(man.rods.length), Rp = perParticle(man.rods.radius);
    out.orientations = new Float32Array(buffers[man.orientations.file]);
    const hl = new Float32Array(N), rad = new Float32Array(N), ax = new Float32Array(N * d);
    for (let i = 0; i < N; i++) {
      hl[i] = 0.5 * out.diameters[i] * (Lp ? Lp[i] : man.rods.length);
      rad[i] = out.diameters[i] * (Rp ? Rp[i] : man.rods.radius);
      for (let a = 0; a < d; a++) ax[i * d + a] = out.orientations[i * d * d + a * d];
    }
    out.rod = { capped: !!man.rods.capped, halfLength: hl, radius: rad, axis: ax };
    if (d === 2) Object.assign(out, VIZ.rodOutlines2D(out));
  }
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
  const need = VIZ.bundleFiles(man);
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
  const need = VIZ.bundleFiles(man);
  for (const f of need) {
    if (!files[f]) throw new Error('missing ' + f);
    buffers[f] = f.endsWith('.json') ? JSON.parse(await files[f].text())
                                     : await files[f].arrayBuffer();
  }
  return VIZ.normalize(man, buffers);
};
