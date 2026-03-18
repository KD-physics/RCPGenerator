"""Pure Python wrapper around the thin rcpgenerator bindings."""

from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path
from typing import Any

from ._rcpgenerator import initialize_particles as _initialize_particles
from ._rcpgenerator import run_packing_observed as _run_packing_observed
from ._rcpgenerator import run_packing as _run_packing
from .render import animate_packing_2d as _animate_packing_2d
from .render import render_packing as _render_packing


class Packing:
    """Stateful v1 public API around the validated dict-based rcpgenerator core."""

    _INIT_KEYS = {"phi", "N", "Ndim", "box", "walls", "fix_height", "dist"}
    _PACK_KEYS = {"box", "walls", "neighbor_max", "seed", "fix_height"}
    _REALIZED_STATE_KEYS = {"positions", "diameters"}
    _RESULT_KEYS = {
        "steps",
        "phi_final",
        "max_min_dist",
        "force_magnitude",
        "phi_history",
        "force_history",
        "energy_history",
    }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Packing":
        init_config = deepcopy(data.get("init_config", {}))
        pack_config = deepcopy(data.get("pack_config", {}))
        realized_state = deepcopy(data.get("realized_state", {}))
        final_result = deepcopy(data.get("final_result", {}))
        trajectory = deepcopy(data.get("trajectory", {}))
        status = deepcopy(data.get("status", {}))

        obj = cls.__new__(cls)
        object.__setattr__(obj, "_mutating_internal_state", True)
        object.__setattr__(obj, "_initialized", False)
        object.__setattr__(obj, "_packed", False)
        object.__setattr__(obj, "_needs_initialize", True)
        object.__setattr__(obj, "_needs_pack", True)

        obj.verbose = init_config.pop("verbose", pack_config.pop("verbose", False))
        obj.phi = init_config.get("phi", 0.05)
        obj.N = init_config.get("N", 0)
        obj.Ndim = init_config.get("Ndim", 0)
        obj.box = deepcopy(init_config.get("box", pack_config.get("box", [])))
        obj.walls = deepcopy(init_config.get("walls", pack_config.get("walls", [])))
        obj.fix_height = init_config.get("fix_height", pack_config.get("fix_height", False))
        obj.dist = deepcopy(init_config.get("dist", {}))
        obj.neighbor_max = pack_config.get("neighbor_max", 0)
        obj.seed = pack_config.get("seed", 0)

        obj.positions = deepcopy(realized_state.get("positions", []))
        obj.diameters = deepcopy(realized_state.get("diameters", []))
        obj.steps = final_result.get("steps", 0)
        obj.phi_final = final_result.get("phi_final")
        obj.max_min_dist = final_result.get("max_min_dist")
        obj.force_magnitude = final_result.get("force_magnitude")
        obj.phi_history = deepcopy(final_result.get("phi_history", []))
        obj.force_history = deepcopy(final_result.get("force_history", []))
        obj.energy_history = deepcopy(final_result.get("energy_history", []))
        obj.trajectory_positions = deepcopy(trajectory.get("positions", []))
        obj.trajectory_diameters = deepcopy(trajectory.get("diameters", []))
        obj.trajectory_steps = deepcopy(trajectory.get("steps", []))
        obj.trajectory_phi = deepcopy(trajectory.get("phi", []))
        obj.trajectory_force = deepcopy(trajectory.get("force", []))
        obj.trajectory_energy = deepcopy(trajectory.get("energy", []))
        obj.trajectory_max_min_dist = deepcopy(trajectory.get("max_min_dist", []))

        object.__setattr__(obj, "_initialized", bool(status.get("initialized", bool(obj.positions) and bool(obj.diameters))))
        object.__setattr__(obj, "_packed", bool(status.get("packed", bool(obj.steps))))
        object.__setattr__(obj, "_needs_initialize", bool(status.get("needs_initialize", not obj._initialized)))
        object.__setattr__(obj, "_needs_pack", bool(status.get("needs_pack", not obj._packed)))
        object.__setattr__(obj, "_mutating_internal_state", False)
        return obj

    def __init__(self, **kwargs: Any) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        object.__setattr__(self, "_initialized", False)
        object.__setattr__(self, "_packed", False)
        object.__setattr__(self, "_needs_initialize", True)
        object.__setattr__(self, "_needs_pack", True)

        self.verbose = kwargs.pop("verbose", False)
        self.phi = kwargs.pop("phi", 0.05)
        self.N = kwargs.pop("N", 0)
        self.Ndim = kwargs.pop("Ndim", 0)
        self.box = deepcopy(kwargs.pop("box", []))
        default_walls = [0] * len(self.box) if self.box else []
        self.walls = deepcopy(kwargs.pop("walls", default_walls))    
        self.fix_height = kwargs.pop("fix_height", False)    
        self.dist = deepcopy(kwargs.pop("dist", {"type": "mono", "d": 1.0}))
        self.neighbor_max = kwargs.pop("neighbor_max", 0)
        self.seed = kwargs.pop("seed", 0)
        
        # self.box = deepcopy(kwargs.pop("box", []))
        # self.walls = deepcopy(kwargs.pop("walls", []))
        # self.fix_height = kwargs.pop("fix_height", False)
        # self.dist = deepcopy(kwargs.pop("dist", {}))
        # self.neighbor_max = kwargs.pop("neighbor_max", 0)
        # self.seed = kwargs.pop("seed", 0)

        self.positions = []
        self.diameters = []
        self.steps = 0
        self.phi_final = None
        self.max_min_dist = None
        self.force_magnitude = None
        self.phi_history = []
        self.force_history = []
        self.energy_history = []
        self.trajectory_positions = []
        self.trajectory_diameters = []
        self.trajectory_steps = []
        self.trajectory_phi = []
        self.trajectory_force = []
        self.trajectory_energy = []
        self.trajectory_max_min_dist = []

        if kwargs:
            unknown = ", ".join(sorted(kwargs))
            raise TypeError(f"Unexpected keyword argument(s): {unknown}")

        object.__setattr__(self, "_mutating_internal_state", False)
        self.initialize()

    def __repr__(self) -> str:
        return (
            "Packing("
            f"N={self.N}, "
            f"Ndim={self.Ndim}, "
            f"initialized={self._initialized}, "
            f"packed={self._packed}, "
            f"needs_initialize={self._needs_initialize}, "
            f"needs_pack={self._needs_pack}"
            ")"
        )

    @property
    def needs_initialize(self) -> bool:
        return self._needs_initialize

    @property
    def needs_pack(self) -> bool:
        return self._needs_pack

    @property
    def initialized(self) -> bool:
        return self._initialized

    @property
    def packed(self) -> bool:
        return self._packed

    def __setattr__(self, name: str, value: Any) -> None:
        object.__setattr__(self, name, value)
        if name.startswith("_"):
            return
        if getattr(self, "_mutating_internal_state", False):
            return
        if name in self._INIT_KEYS:
            self._mark_needs_initialize()
        elif name in self._PACK_KEYS or name in self._REALIZED_STATE_KEYS:
            self._mark_needs_pack()

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message)

    def _mark_needs_initialize(self) -> None:
        object.__setattr__(self, "_needs_initialize", True)
        object.__setattr__(self, "_needs_pack", True)
        object.__setattr__(self, "_packed", False)

    def _mark_needs_pack(self) -> None:
        object.__setattr__(self, "_needs_pack", True)
        object.__setattr__(self, "_packed", False)

    def _collect_init_config(self) -> dict[str, Any]:
        return {
            "phi": self.phi,
            "N": self.N,
            "Ndim": self.Ndim,
            "box": deepcopy(self.box),
            "walls": deepcopy(self.walls),
            "fix_height": self.fix_height,
            "dist": deepcopy(self.dist),
        }

    def _collect_pack_input(self) -> dict[str, Any]:
        return {
            "positions": deepcopy(self.positions),
            "diameters": deepcopy(self.diameters),
        }

    def _collect_pack_config(self) -> dict[str, Any]:
        return {
            "box": deepcopy(self.box),
            "walls": deepcopy(self.walls),
            "neighbor_max": self.neighbor_max,
            "seed": self.seed,
            "fix_height": self.fix_height,
        }

    def _collect_relax_options(
        self,
        n_steps: int,
        mu: float | None,
        fix_diameter: bool,
        target_phi: float | None,
    ) -> dict[str, Any]:
        options = {
            "max_steps": int(n_steps),
            "fix_diameter": bool(fix_diameter),
        }
        if mu is not None:
            options["mu"] = float(mu)
        if target_phi is not None:
            options["target_phi"] = float(target_phi)
        return options

    def _clear_results(self) -> None:
        self.steps = 0
        self.phi_final = None
        self.max_min_dist = None
        self.force_magnitude = None
        self.phi_history = []
        self.force_history = []
        self.energy_history = []

    def _clear_trajectory(self) -> None:
        self.trajectory_positions = []
        self.trajectory_diameters = []
        self.trajectory_steps = []
        self.trajectory_phi = []
        self.trajectory_force = []
        self.trajectory_energy = []
        self.trajectory_max_min_dist = []

    def _current_phi(self) -> float:
        if not self.positions or not self.diameters or not self.box:
            raise ValueError("Current realized state is unavailable")
        box_volume = 1.0
        for length in self.box:
            box_volume *= length
        particle_volume = 0.0
        coeff = math.pi ** (self.Ndim / 2.0) / math.gamma(self.Ndim / 2.0 + 1.0)
        for diameter in self.diameters:
            particle_volume += coeff * (diameter / 2.0) ** self.Ndim
        phi_modifier = 1.0
        if self.walls and self.walls[0] < 0:
            phi_modifier = math.pi ** ((-self.walls[0]) / 2.0) / math.gamma(((-self.walls[0]) / 2.0) + 1.0) * (0.5 ** (-self.walls[0]))
        return particle_volume / box_volume / phi_modifier

    def _mark_realized_state_updated(self) -> None:
        self._clear_results()
        self._clear_trajectory()
        self._mark_needs_pack()

    def _clear_realized_state(self) -> None:
        self.positions = []
        self.diameters = []
        self._clear_results()
        self._clear_trajectory()

    def _apply_initializer_result(self, result: dict[str, Any]) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.positions = result["positions"]
            self.diameters = result["diameters"]
            self.box = result["box"]
            self.walls = result["walls"]
            self._clear_results()
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

        object.__setattr__(self, "_initialized", True)
        object.__setattr__(self, "_needs_initialize", False)
        object.__setattr__(self, "_needs_pack", True)
        object.__setattr__(self, "_packed", False)

    def _apply_packing_result(self, result: dict[str, Any]) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.positions = result["positions"]
            self.diameters = result["diameters"]
            self.box = result["box"]
            self.walls = result["walls"]
            self.steps = result["steps"]
            self.phi_final = result["phi"]
            self.max_min_dist = result["max_min_dist"]
            self.force_magnitude = result["force_magnitude"]
            self.phi_history = result["phi_history"]
            self.force_history = result["force_history"]
            self.energy_history = result["energy_history"]
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

        object.__setattr__(self, "_packed", True)
        object.__setattr__(self, "_needs_pack", False)

    def _apply_trajectory(self, trace: dict[str, Any]) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.trajectory_positions = trace.get("positions", [])
            self.trajectory_diameters = trace.get("diameters", [])
            self.trajectory_steps = trace.get("steps", [])
            self.trajectory_phi = trace.get("phi", [])
            self.trajectory_force = trace.get("force", [])
            self.trajectory_energy = trace.get("energy", [])
            self.trajectory_max_min_dist = trace.get("max_min_dist", [])
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

    def _execute_packing(
        self,
        *,
        verbose: bool,
        progress_interval: int,
        capture_trajectory: bool,
        trajectory_interval: int,
        options: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        self._clear_trajectory()
        if verbose:
            print("rcpgenerator: running packing")

        if verbose or capture_trajectory or (options is not None and options):
            def on_progress(update: dict[str, Any]) -> None:
                if verbose:
                    print(
                        f"step={update['step']} "
                        f"phi={update['phi']:.6f} "
                        f"force={update['force_magnitude']:.6g} "
                        f"energy={update['energy']:.6g}"
                    )

            observed = _run_packing_observed(
                self._collect_pack_input(),
                self._collect_pack_config(),
                progress_interval if verbose else 0,
                capture_trajectory,
                trajectory_interval if capture_trajectory else 0,
                on_progress if verbose else None,
                options or {},
            )
            result = observed["result"]
            self._apply_packing_result(result)
            self._apply_trajectory(observed["trace"])
        else:
            result = _run_packing(self._collect_pack_input(), self._collect_pack_config())
            self._apply_packing_result(result)

        if verbose:
            print("rcpgenerator: packing complete")
        return result

    def initialize(self) -> dict[str, Any]:
        self._log("rcpgenerator: initializing particles")
        result = _initialize_particles(self._collect_init_config())
        self._apply_initializer_result(result)
        self._log("rcpgenerator: initialization complete")
        return result

    def pack(
        self,
        verbose: bool | None = None,
        progress_interval: int = 1000,
        capture_trajectory: bool = False,
        trajectory_interval: int = 1000,
    ) -> dict[str, Any]:
        effective_verbose = self.verbose if verbose is None else verbose
        if self._needs_initialize or not self._initialized:
            if effective_verbose:
                print("rcpgenerator: initialization required before packing")
            self.initialize()
        elif self._needs_pack and effective_verbose:
            print("rcpgenerator: rerunning packing on current realized state")

        return self._execute_packing(
            verbose=effective_verbose,
            progress_interval=progress_interval,
            capture_trajectory=capture_trajectory,
            trajectory_interval=trajectory_interval,
        )

    def relax(
        self,
        n_steps: int,
        mu: float | None = None,
        fix_diameter: bool = False,
        target_phi: float | None = None,
        verbose: bool = False,
        progress_interval: int = 1000,
        capture_trajectory: bool = False,
        trajectory_interval: int = 1000,
    ) -> dict[str, Any]:
        if n_steps <= 0:
            raise ValueError("n_steps must be positive")
        if target_phi is not None and fix_diameter:
            raise ValueError("target_phi cannot be combined with fix_diameter=True")
        if self._needs_initialize or not self._initialized:
            if verbose:
                print("rcpgenerator: initialization required before relaxation")
            self.initialize()
        elif self._needs_pack and verbose:
            print("rcpgenerator: relaxing current realized state")

        return self._execute_packing(
            verbose=verbose,
            progress_interval=progress_interval,
            capture_trajectory=capture_trajectory,
            trajectory_interval=trajectory_interval,
            options=self._collect_relax_options(
                n_steps=n_steps,
                mu=mu,
                fix_diameter=fix_diameter,
                target_phi=target_phi,
            ),
        )

    def render(
        self,
        path: str | Path | None = None,
        show: bool = False,
        palette_choice: int = 1,
        figsize=(4, 4),
        dpi=110,
    ) -> str | None:
        return _render_packing(
            self,
            path=path,
            show=show,
            palette_choice=palette_choice,
            figsize=figsize,
            dpi=dpi,
        )

    def savefig(
        self,
        path: str | Path,
        show: bool = False,
        palette_choice: int = 1,
    ) -> str:
        return self.render(path=path, show=show, palette_choice=palette_choice)
    
    def show_packing(
        self,
        palette_choice: int = 1,
        figsize=(4, 4),
        dpi=110,
    ) -> None:
        self.render(
            path=None,
            show=True,
            palette_choice=palette_choice,
            figsize=figsize,
            dpi=dpi,
        )
        
    def animate_2d(
        self,
        path: str | Path | None = None,
        show: bool = True,
        palette_choice: int = 1,
        interval_ms: int = 80,
        repeat: bool = False,
    ):
        return _animate_packing_2d(
            self,
            path=path,
            show=show,
            palette_choice=palette_choice,
            interval_ms=interval_ms,
            repeat=repeat,
        )

    def summary(self) -> str:
        return (
            f"Packing summary: N={self.N}, Ndim={self.Ndim}, phi={self.phi}, "
            f"initialized={self.initialized}, packed={self.packed}, "
            f"needs_initialize={self.needs_initialize}, needs_pack={self.needs_pack}, "
            f"phi_final={self.phi_final}, steps={self.steps}, "
            f"force_magnitude={self.force_magnitude}, "
            f"trajectory_samples={len(self.trajectory_steps)}"
        )

    def reset(self) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self._clear_realized_state()
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)
        object.__setattr__(self, "_initialized", False)
        object.__setattr__(self, "_packed", False)
        object.__setattr__(self, "_needs_initialize", True)
        object.__setattr__(self, "_needs_pack", True)

    def update_phi(self, target_phi: float) -> float:
        if target_phi <= 0:
            raise ValueError("target_phi must be positive")
        current_phi = self._current_phi()
        exponent = self.Ndim - 1 if self.fix_height else self.Ndim
        scale = (target_phi / current_phi) ** (1.0 / exponent)

        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.diameters = [diameter * scale for diameter in self.diameters]
            if self.fix_height:
                updated_box = list(self.box)
                updated_box[self.Ndim - 1] *= scale
                self.box = updated_box
                self.positions = [
                    [
                        coord * scale if dim == self.Ndim - 1 else coord
                        for dim, coord in enumerate(position)
                    ]
                    for position in self.positions
                ]
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

        self._mark_realized_state_updated()
        return self._current_phi()

    def update_box(self, box: list[float]) -> list[float]:
        if len(box) != self.Ndim:
            raise ValueError("box length must match Ndim")
        if not self.positions or not self.diameters or not self.box:
            raise ValueError("Current realized state is unavailable")
        old_box = list(self.box)
        old_volume = math.prod(old_box)
        new_box = [float(length) for length in box]
        new_volume = math.prod(new_box)
        if old_volume <= 0 or new_volume <= 0:
            raise ValueError("box dimensions must be positive")
        diameter_scale = (new_volume / old_volume) ** (1.0 / self.Ndim)
        position_scales = [new_box[dim] / old_box[dim] for dim in range(self.Ndim)]

        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.box = new_box
            self.positions = [
                [coord * position_scales[dim] for dim, coord in enumerate(position)]
                for position in self.positions
            ]
            self.diameters = [diameter * diameter_scale for diameter in self.diameters]
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

        self._mark_realized_state_updated()
        return list(self.box)

    def update_height(self, height: float) -> list[float]:
        new_box = list(self.box)
        if not new_box:
            raise ValueError("Current box is unavailable")
        new_box[self.Ndim - 1] = float(height)
        return self.update_box(new_box)

    def copy(self) -> "Packing":
        return type(self).from_dict(self.to_dict())

    def clone(self) -> "Packing":
        return self.copy()

    def to_dict(self) -> dict[str, Any]:
        return {
            "init_config": self._collect_init_config(),
            "pack_config": self._collect_pack_config(),
            "realized_state": self._collect_pack_input(),
            "final_result": {
                "steps": self.steps,
                "phi_final": self.phi_final,
                "max_min_dist": self.max_min_dist,
                "force_magnitude": self.force_magnitude,
                "phi_history": deepcopy(self.phi_history),
                "force_history": deepcopy(self.force_history),
                "energy_history": deepcopy(self.energy_history),
            },
            "trajectory": {
                "positions": deepcopy(self.trajectory_positions),
                "diameters": deepcopy(self.trajectory_diameters),
                "steps": deepcopy(self.trajectory_steps),
                "phi": deepcopy(self.trajectory_phi),
                "force": deepcopy(self.trajectory_force),
                "energy": deepcopy(self.trajectory_energy),
                "max_min_dist": deepcopy(self.trajectory_max_min_dist),
            },
            "status": {
                "initialized": self._initialized,
                "packed": self._packed,
                "needs_initialize": self._needs_initialize,
                "needs_pack": self._needs_pack,
            },
        }
