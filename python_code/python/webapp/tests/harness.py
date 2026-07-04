"""Test harness: serve viz/ over http + drive headless Chromium (WebGL via the
verified swiftshader recipe) + read screenshots. Shared by all viz tests."""
import threading, functools, http.server, socketserver, pathlib, contextlib
import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright

VIZ = pathlib.Path(__file__).resolve().parent.parent
GL_ARGS = ["--use-gl=angle", "--use-angle=swiftshader"]
BG = (13, 13, 16)

class _Reuse(socketserver.TCPServer):
    allow_reuse_address = True
    def log_message(self, *a): pass

@contextlib.contextmanager
def serve(port=0):
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(VIZ))
    httpd = _Reuse(("127.0.0.1", port), handler)
    realport = httpd.server_address[1]
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    try: yield f"http://127.0.0.1:{realport}"
    finally: httpd.shutdown()

@contextlib.contextmanager
def page(viewport=(800, 800)):
    with sync_playwright() as p:
        b = p.chromium.launch(args=GL_ARGS)
        pg = b.new_page(viewport={"width": viewport[0], "height": viewport[1]})
        pg.on("console", lambda m: pg.__dict__.setdefault("_console", []).append(m))
        try: yield pg
        finally: b.close()

def load(pg, base, bundle, timeout=20000, ui=False):
    errs = []
    pg.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
    pg.goto(f"{base}/app/index.html?bundle=/bundles/{bundle}&ui={1 if ui else 0}")
    pg.wait_for_function("window.__diag && window.__diag.ready && window.__rendered === true", timeout=timeout)
    pg.wait_for_timeout(300)
    pg._errs = errs
    return pg.evaluate("window.__diag")

def shot(pg, name):
    out = VIZ / "tests" / "shots"; out.mkdir(exist_ok=True)
    path = out / f"{name}.png"; pg.screenshot(path=str(path))
    return np.asarray(Image.open(path).convert("RGB")), path

FILL = (90, 150, 230)   # default particle color in render2d.js

def colored_fraction(arr, bg=BG, tol=40):
    diff = np.abs(arr.astype(int) - np.array(bg)).sum(axis=2)
    return float((diff > tol).mean())

def coverage_fraction(arr, fill=FILL, bg=BG):
    """Anti-alias-aware area estimate: fractional coverage per pixel (blended edge
    pixels count proportionally). Unbiased vs binary thresholding -> area = phi."""
    full = abs(np.array(fill) - np.array(bg)).sum()
    diff = np.abs(arr.astype(int) - np.array(bg)).sum(axis=2)
    return float(np.clip(diff / full, 0, 1).mean())

def is_colored(arr, x, y, bg=BG, tol=40):
    px = arr[y, x].astype(int)
    return abs(px - np.array(bg)).sum() > tol

def no_outline(pg):
    """Disable the particle outline + boundary so area-fraction is measured on the
    fill color alone (outlines/boundary are real features, just not part of phi)."""
    pg.evaluate("window.__viz && (window.__viz.state.outline=0, window.__viz.state.showBoundary=false, window.__viz.refresh(), window.__viz.fitTight())")
    pg.wait_for_timeout(250)

def set_view(pg, vs):
    pg.evaluate("vs => window.__viz.setViewState(vs)", vs)
    pg.wait_for_timeout(250)

def colored_radius(arr, cx, cy, bg=BG, tol=40):
    """colored extent along +x from (cx,cy) until background."""
    h, w = arr.shape[:2]
    for r in range(1, w - cx):
        if not is_colored(arr, cx + r, cy, bg, tol):
            return r
    return -1
