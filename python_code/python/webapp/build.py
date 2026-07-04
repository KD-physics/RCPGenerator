#!/usr/bin/env python3
"""Inline vendor + js modules + index.html into a single self-contained HTML
(viz/dist/packing_viz.html) — the standalone deliverable (drag-drop a bundle onto it).
Modular dev (separate files) -> single distributable build."""
import re, pathlib
APP = pathlib.Path(__file__).resolve().parent / "app"
DIST = pathlib.Path(__file__).resolve().parent / "dist"

def build():
    html = (APP / "index.html").read_text()
    def inline(m):
        src = m.group(1); code = (APP / src).read_text()
        code = code.replace("</script>", "<\\/script>")   # avoid premature tag close
        return "<script>\n" + code + "\n</script>"
    html = re.sub(r'<script src="([^"]+)"></script>', inline, html)
    DIST.mkdir(exist_ok=True)
    out = DIST / "packing_viz.html"; out.write_text(html)
    print(f"built -> {out}  ({len(html)//1024} KB, self-contained)")
    return out

if __name__ == "__main__":
    build()
