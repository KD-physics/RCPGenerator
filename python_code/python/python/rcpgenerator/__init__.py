"""Thin Python bindings and a small stateful wrapper for rcpgenerator."""

from ._rcpgenerator import (
    initialize_particles,
    run_packing,
    set_num_threads,
    get_num_threads,
)
from .render import animate_packing_2d, render_packing
from ._wrapper import Packing

__all__ = [
    "Packing",
    "animate_packing_2d",
    "get_num_threads",
    "initialize_particles",
    "render_packing",
    "run_packing",
    "set_num_threads",
]

