# RCP Packing Visualizer

A standalone, cross-platform webapp to inspect RCP packings — "Google Maps for
packings": smooth pan/zoom at constant crisp resolution from the whole box down to a
single fine particle, 2D directly or 3D via a slidable cross-section, with pore /
Voronoi overlays and live field-of-view statistics. Built to validate the solver and
inspect the power-law structural findings; useful for packing inspection generally.

## Use it (no install for the end user)
- **Drag-drop:** open `dist/packing_viz.html` in any browser; drag a bundle folder
  (or `.zip`) onto it. Standalone, offline, no server.
- **CLI:** `python launch.py <packing.npz | bundle_dir>` -> builds a self-contained
  HTML with the data embedded and opens it (no server). `--no-open`, `--out FILE`.
- **Dev/large data:** serve the folder (`python -m http.server` in `viz/`) and open
  `app/index.html?bundle=/bundles/<name>`.

## Make a bundle from a packing
`python bridge.py <packing.npz> <out_dir> [--layers]` — writes `manifest.json` + raw
float32 `pos.f32`/`dia.f32` (+ optional layers: Voronoi `local_phi` (3D), segmented
`pores` (2D), computed via the analysis toolbox). `--ground-truth <dir>` makes the
tiny test fixture.

## Features
pan/zoom (crisp SDF circles, periodic wrap) · 3D cross-section slider · color-by
diameter / Voronoi local-phi / any scalar layer · colormaps (viridis = colorblind-
safe, turbo, gray, coolwarm) · segmented-pore-map overlay (per-pore eta) · live FOV
stats (P(local-phi), P(diameter), running phi) · PNG export at chosen DPI.

## Architecture (modular — features/scale/tweaks are add-ons, not refactors)
`app/` classic-script modules: `loader` (bundle/embedded/drop) · `render2d` (deck.gl
OrthographicView + ScatterplotLayer, binary attrs) · `wrap` (periodic ghosts) ·
`slice` (3D->2D) · `colormap`/`color` · `overlay_pores` (registry: `VIZ.overlays`) ·
`stats` · `export` · `diag` · `ui` · `app`. Vendored `vendor/deck.gl.min.js`
(offline). `bridge.py` (npz->bundle), `build.py` (-> single `dist/packing_viz.html`),
`launch.py` (CLI). Add an overlay/stat/colormap = a module + a manifest entry.

## Tests (Playwright headless Chromium, WebGL via swiftshader)
`python tests/test_{core,nav,3d,overlay,ui,standalone,pores}.py` — 44 checks:
ground-truth pixels, area-fraction=phi (2D + 3D slice), zoom scaling, pan, periodic
wrap, large-N load, color-by, palettes, FOV stats, PNG-DPI, standalone file://, pore
overlay vs toolbox. `python gallery.py` -> `gallery/gallery.png` (visual review).
Note: software-GL renders the same pixels (correctness valid); fps is a floor —
on-GPU (e.g. RTX) is faster.
