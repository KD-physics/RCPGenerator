"""Thin Python bindings and a small stateful wrapper for rcpgenerator."""

from ._rcpgenerator import initialize_particles, run_packing
from .render import animate_packing_2d, render_packing
from ._wrapper import Packing

__all__ = ["Packing", "animate_packing_2d", "initialize_particles", "render_packing", "run_packing"]

