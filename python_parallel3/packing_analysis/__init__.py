"""packing_analysis — structure-analysis library for RCP packings.

Canonical toolbox modules (importable):
  voidfill_method  — Farr-Groot fg_rcp, greedy void-fill kernel, radical Voronoi
  porenet          — periodic Delaunay pore network (throats, gap-merge)
  porefill         — per-cavity fill
  chord_cut / chord_study / chord_boxbin — chord-gap / FG decomposition

Plus an expandable ANALYSIS REGISTRY (registry.py) that turns named analyses
("local_phi", "pores", ...) into webapp bundle layers + summary metrics. Add a
new analysis by registering one function — no API change anywhere downstream.

Importing this package is cheap and dependency-free: the heavy deps (pyvoro for
Voronoi, scipy for the pore network) load LAZILY, only when an analysis actually
runs. A missing dep degrades gracefully (skip-with-warning), never an import error.
This is a general-purpose library — no environment- or path-specific assumptions.
"""
import os as _os, sys as _sys

# Put this package dir on sys.path so the toolbox modules' absolute intra-imports
# (`import voidfill_method`, `import porenet as pn`, ...) resolve, and the modules
# still work when run as standalone CLI scripts. Idempotent.
_pkg_dir = _os.path.dirname(_os.path.abspath(__file__))
if _pkg_dir not in _sys.path:
    _sys.path.insert(0, _pkg_dir)

from . import registry  # lightweight: numpy only at import time
from .registry import run as run_analyses, available as available_analyses, ANALYSES

__all__ = ["registry", "run_analyses", "available_analyses", "ANALYSES"]
