"""Pure Python wrapper around the thin rcpgenerator bindings."""

from __future__ import annotations

from copy import deepcopy
import math
import numbers
import warnings
from pathlib import Path
from typing import Any

from ._rcpgenerator import initialize_particles as _initialize_particles
from ._rcpgenerator import run_packing_observed as _run_packing_observed
from ._rcpgenerator import run_packing as _run_packing
from ._rcpgenerator import contact_scale as _contact_scale
from ._rcpgenerator import random_orientations as _random_orientations


def _container_fraction(walls):
    """Box fraction filled by a curved (hyperspherical) container: walls[0] = -t marks the
    first t dims as one t-ball of diameter box[0]; returns 1.0 for a flat / periodic box.
    Reported phi is divided by this so a confined packing is referenced to the actual
    container volume (disk / cylinder / sphere), not the bounding box."""
    t = -int(walls[0]) if (walls and int(walls[0]) < 0) else 0
    if t < 2:
        return 1.0
    return math.pi ** (t / 2.0) / math.gamma(t / 2.0 + 1.0) * (0.5 ** t)
from .render import animate_packing_2d as _animate_packing_2d
from .render import render_packing as _render_packing


# Packing.quench() outcome, by the core's quench_status code
_QUENCH_STATUS = {1: "converged", 2: "max_iter", 3: "stalled", 4: "unjammed", 5: "step_collapse", 6: "non_finite"}


class Packing:
    """Stateful v1 public API around the validated dict-based rcpgenerator core.

    Diameter permutation: pass ``diameter_permute=True`` to enable the optional move
    that permutes which particle carries which diameter during the run (positions are
    never moved by it, the diameter multiset is preserved exactly). Two knobs set the
    rate: ``permute_budget`` (cap on exchanges per pass as a fraction of N, default
    0.25) and ``permute_interval`` (steps between passes, default 5); a third,
    ``permute_start`` (step at which permutation begins, default 4000), is rarely
    needed. The defaults are the validated recipe (capacity-matched exchanges with
    3-cycles); a smaller budget, a longer interval or a later start gives less uplift
    on a smooth scale down to the pristine result. ``diameter_permute=False`` (default)
    is the original algorithm bit for bit. See README "Diameter permutation". After
    ``pack()``, ``permute_counts``, ``permute_passes`` and ``permute_exchanges`` hold
    the permutation diagnostics. (The recipe is assembled
    from internal fields under a legacy ``swap_`` prefix; they are not part of the
    public API.)

    Ellipsoids: pass ``aspect_ratios`` (Ndim - 1 numbers) to pack ellipsoids instead of spheres. The ratios
    are written in DECREASING order, each at most 1 (``1 >= ratio_2 >= ratio_3 >= ... > 0``): ``D`` is the full
    length of a particle along its first body axis, which is its LONGEST axis, and axis k has length
    ``D * aspect_ratios[k-2]`` (2D 2:1 ellipse ``[0.5]``; 3D prolate ``[0.5, 0.5]``, oblate ``[1, 0.5]``).
    Input that breaks the order raises ``ValueError`` and is never reordered. One list applies to every
    particle; a list of N lists gives a shape per particle. Row i goes to the i-th largest particle when the
    diameters are generated (the order of ``diameters`` after construction), and with entry i of the user's own
    diameter list when ``dist`` is of type ``"custom"`` (rod lists follow the same rule). After ``pack()``, ``aspect_ratios`` is the table in the
    packing's own particle order.

    Rods: pass ``rod_length`` and ``rod_radius`` (both, each one number or a list of N, in units of the
    particle's ``D``) and ``rod_capped``. Capped (``True``, the default): a straight section of length
    ``D * rod_length`` with a half-sphere of radius ``D * rod_radius`` on each end, so ``rod_length + 2 *
    rod_radius`` from tip to tip, and length 0 is a sphere. Uncapped (``False``): flat ends, ``rod_length`` is
    the full length; 2D and 3D only. With ``rod_radius=0.5``, ``D`` is the rod's thickness. Rods and aspect
    ratios cannot be combined.

    For both shapes the initializer generates random orientations (``orientations``: one Ndim x Ndim rotation
    matrix per particle, columns = body axes in the lab frame; column 0 is the longest axis, or the rod's
    axis); pass ``orientations=`` to supply your own. After ``pack()``: ``orientations`` holds the final
    rotations, ``phi_final`` the raw (slightly overlapping) packing fraction as for spheres, and
    ``phi_contact`` the packing fraction after the exact uniform shrink to zero overlap (``shrink_factor``):
    the value to report. Without ``aspect_ratios`` or the rod keywords nothing changes. Available with both
    shapes: periodic boxes, flat hard walls, curved containers and ``fix_height``; ``diameter_permute=True`` (an exchange trades diameter and
    shape; orientations stay in place by default, ``ell_permute_orient=1`` makes them travel); ``quench()``
    (force and torque balance; defaults ``kappa_mult=1.01`` and ``dt_max=0.2``; rods are balanced as slightly
    rounded rods, ``rod_round=0.01``, in 2D and 3D only); ``render`` / ``savefig`` / ``show_packing`` /
    ``animate_2d``, which draw the true shapes and orientations. Not available with either shape (each
    raises): the research inputs (distribution morph, fixed particles, internal ``swap_*`` settings). Flat-ended
    rods and the quench of rods are 2D and 3D only. See README "Ellipsoids and rods".
    """

    _INIT_KEYS = {"phi", "N", "Ndim", "box", "walls", "fix_height", "dist", "aspect_ratios",
                  "rod_length", "rod_radius", "rod_capped"}
    # Ellipsoid controls that travel in the pack config (sent to the core only when aspect_ratios is set).
    _ELL_DEFAULTS = {"ell_iters": 5, "ell_rot_scale": 1.0, "ell_track_tol": 1e-6,
                     "ell_permute_orient": 0, "ell_permute_partners": 64, "ell_permute_guard": 1,
                     "ell_capture": 0}
    # Public diameter-permutation controls (translated onto the internal swap_* fields in
    # _collect_pack_config; inert - and bitwise identical to the pristine core - when off).
    _PERMUTE_DEFAULTS = {"diameter_permute": False, "permute_budget": 0.25, "permute_interval": 5, "permute_start": 4000}
    # [219] Distribution morph + mu pin: an undocumented research feature (travels in the pack config; inert when
    # morph_target is empty and mu_pin <= 0). morph_target: N relative diameters (any order).
    _MORPH_DEFAULTS = {"morph_target": [], "morph_start": 0, "morph_interval": 10, "morph_rate0": 1e-3,
                       "morph_band": 2.0, "morph_grow": 1.25, "morph_shrink": 0.7, "morph_rate_min": 1e-5, "morph_grow_after": 5,
                       "morph_rate_max": 3e-5, "morph_settle_flat": 0.0, "morph_settle_window": 1000,
                       "morph_kappa_damp": 1.0, "morph_band_floor": 5e-4, "morph_margin": 0.5, "morph_max_hold_streak": 2000, "morph_kappa_drop": 1e-3, "morph_kappa_tau": 100,
                       "mu_pin": 0.0}
    # [221] quench to force balance: core fields behind Packing.quench() (inert when quench_mode == 0). The public
    # interface is the quench() method; these fields are internal, kept for checkpoint-path runs and the research record.
    _QUENCH_DEFAULTS = {"quench_stall_window": 0, "quench_stall_factor": 0.5, "quench_dt_min": 0.0, "quench_mode": 0, "quench_tol": 1e-10, "quench_max_iter": 200000, "quench_dt0": 0.05, "quench_dt_max": 0.5,
                        "quench_alpha0": 0.1, "quench_mass_scaled": 0, "quench_report_every": 500, "quench_kappa_mult": 1.0, "quench_mix_mode": 0, "quench_two_level": 0, "quench_light_frac": 0.3, "quench_light_omega": 1.0, "quench_warm_steps": 0, "quench_warm_amp": 0.1, "quench_freeze_damp": 0.0, "quench_freeze_frac": 0.1, "quench_class_every": 100, "quench_held": [], "quench_rod_round": 0.0}
    # Internal diameter-permutation fields (legacy ``swap_`` prefix, not public API); all
    # inert at swap == 0. They travel in the pack config so they reach the core through
    # pack() (PackingRunOptions does not).
    _SWAP_DEFAULTS = {
        "swap": 0.0, "swap_interval": 1, "swap_mode": 5, "swap_rank_window": 8,
        "swap_size_tol": 0.2, "swap_eps": 0.0, "swap_start": 0, "swap_stop": -1,
        "swap_phase": 0, "swap_adam": 0, "swap_inject_adam": -1, "swap_inject_freeze": 0, "swap_inject_accept": 0, "swap_inject_tol": 0.0, "swap_inject_freeze_shell": 1, "swap_max_per_particle": 4, "swap_refresh": 0, "swap_backoff": 64, "swap_exact": 0, "swap_window_rel": 0.016, "swap_relayout_ratio": 1.3, "swap_rank_space": 1, "swap_idle_frac": 1e-2,
        "swap_ctrl": 0, "swap_a": 0.3, "swap_R": 0.20, "swap_eta": 0.2,
        "swap_wmin": 5e-4, "swap_wmax": 0.10, "swap_smax": 32.0,
        "swap_floor_budget": 1.0,
        "swap_rescue": 0.0, "swap_rescue_w": 0.05, "swap_rescue_interval": 4,
        "swap_cycle3": 0.0,
        # session survivability (not swap knobs, but they travel in the pack config)
        "checkpoint_path": "", "checkpoint_interval_s": 1800.0, "checkpoint_every_steps": 0, "mu_hold_start": -1, "mu_hold_steps": 0, "mu_hold_mult": 1.0, "mu_hold_shift_deadline": 1, "record_visibility": 0, "swap_mu_ctrl": 0, "mu_ctrl_v_target": 0.30, "mu_ctrl_deadband": 0.05, "mu_ctrl_eta": 0.5, "mu_ctrl_interval": 50, "mu_ctrl_min_mult": 0.1, "mu_ctrl_max_mult": 1.0, "swap_geo": 0, "swap_geo_T": 0.0, "swap_geo_T_end": -1.0, "swap_geo_T_decay": 1.0, "swap_geo_w": 0.05, "swap_geo_B": 2000, "swap_geo_k": 200, "swap_geo_start": 0, "swap_geo_stop": -1, "swap_geo_reopen": 0, "swap_geo_cycle": 0, "swap_geo_cycle_cap": 12, "swap_geo_cycle4": 0, "swap_geo_cycle4_cap": 6, "swap_geo_capk": 0, "swap_geo_capmode": 3, "swap_geo_capk_start": 0, "swap_geo_cage": 0, "swap_geo_cage_move": 0, "swap_geo_bottleneck": 0, "swap_geo_cage_start": 0, "swap_geo_nonsep": 0, "swap_geo_ctarget": 0.5, "swap_geo_nsk": 4, "swap_geo_nonsep_start": 0, "swap_geo_nonsep_stop": -1, "max_steps": 0,
        "swap_event_log_stride": 0, "swap_event_log_cap": 1000000, "swap_event_log_geom": 0,
        "swap_sweeps_off": 0, "swap_inject": [],
    }
    _PACK_KEYS = {"box", "walls", "neighbor_max", "seed", "fix_height"} | set(_ELL_DEFAULTS) | set(_SWAP_DEFAULTS) | set(_PERMUTE_DEFAULTS) | set(_MORPH_DEFAULTS) | set(_QUENCH_DEFAULTS)
    _REALIZED_STATE_KEYS = {"positions", "diameters", "orientations"}
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
        obj._init_extension_state(pack_config, consume=False)
        obj._init_ellipsoid_state(init_config.get("aspect_ratios", pack_config.get("aspect_ratios", [])),
                                  pack_config)
        obj._init_rod_state(init_config.get("rod_length"), init_config.get("rod_radius"), init_config.get("rod_capped", True))
        for key in ("rod_length", "rod_radius"):   # per-particle lists in the order of the realized state
            if obj.is_rod and isinstance(getattr(obj, key), list) and realized_state.get(key):
                setattr(obj, key, [float(v) for v in realized_state[key]])

        obj.positions = deepcopy(realized_state.get("positions", []))
        obj.diameters = deepcopy(realized_state.get("diameters", []))
        obj.orientations = obj._nested_orientations(realized_state.get("orientations", []))
        if realized_state.get("aspect_ratios"):   # per-particle shapes in the order of the realized state
            obj.aspect_ratios = [[float(a) for a in row] for row in realized_state["aspect_ratios"]]
        obj.phi_contact = final_result.get("phi_contact")
        obj.shrink_factor = final_result.get("shrink_factor")
        obj.torque_magnitude = final_result.get("torque_magnitude")
        obj.trajectory_orientations = deepcopy(trajectory.get("orientations", []))
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
        obj.trajectory_origin = deepcopy(trajectory.get("origins", []))

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
        self._init_extension_state(kwargs, consume=True)
        self._init_ellipsoid_state(kwargs.pop("aspect_ratios", None),
                                   {k: kwargs.pop(k) for k in list(self._ELL_DEFAULTS) if k in kwargs})
        self._init_rod_state(kwargs.pop("rod_length", None), kwargs.pop("rod_radius", None), kwargs.pop("rod_capped", True))
        user_orientations = kwargs.pop("orientations", None)
        
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
        self.max_overlap_history = []
        self.mu_flag_history = []
        self.mu_history = []
        self.alpha_history = []
        self.trajectory_positions = []
        self.trajectory_diameters = []
        self.trajectory_steps = []
        self.trajectory_phi = []
        self.trajectory_force = []
        self.trajectory_energy = []
        self.trajectory_max_min_dist = []
        self.trajectory_origin = []

        if kwargs:
            unknown = ", ".join(sorted(kwargs))
            raise TypeError(f"Unexpected keyword argument(s): {unknown}")

        object.__setattr__(self, "_mutating_internal_state", False)
        self.initialize()
        if user_orientations is not None:
            if not self._nonsphere:
                raise ValueError("orientations were given without aspect_ratios (or rod_length and rod_radius): "
                                 "supply both (ellipsoids or rods) or neither (spheres)")
            object.__setattr__(self, "_mutating_internal_state", True)
            try:
                self.orientations = self._nested_orientations(user_orientations, check=True)
            finally:
                object.__setattr__(self, "_mutating_internal_state", False)

    # ------------------------------------------------------------------ ellipsoids ----------------------------
    def _init_ellipsoid_state(self, aspect_ratios, config: dict[str, Any]) -> None:
        """Shape and orientation state (all empty / None for spheres)."""
        self.aspect_ratios = self._normalised_ratios(aspect_ratios)
        # the shapes as the user gave them: for a per-particle table, row i belongs to particle i of the
        # INITIALISED packing (the order of `diameters` after initialize(): largest first). pack() returns the
        # particles in another order and the table with them, so the table for a re-initialisation is kept apart
        object.__setattr__(self, "_aspect_init", deepcopy(self.aspect_ratios))
        for k, v in self._ELL_DEFAULTS.items():
            setattr(self, k, type(v)(config.get(k, v)))
        self.orientations = []
        self.phi_contact = None       # packing fraction after the exact shrink to zero overlap (ellipsoid packs)
        self.shrink_factor = None     # the uniform factor of that shrink (multiply every diameter by it)
        self.torque_magnitude = None
        self.trajectory_orientations = []

    # ------------------------------------------------------------------ rods ----------------------------------
    def _init_rod_state(self, rod_length, rod_radius, rod_capped) -> None:
        """Rod shape (None for spheres and ellipsoids): rod_length and rod_radius in units of each particle's D,
        each one number or a list of N (entry i belongs to particle i of the INITIALISED packing, largest
        first); rod_capped one flag for the packing."""
        self.rod_length = self._normalised_rod(rod_length)
        self.rod_radius = self._normalised_rod(rod_radius)
        self.rod_capped = bool(rod_capped)
        object.__setattr__(self, "_rod_init", (deepcopy(self.rod_length), deepcopy(self.rod_radius)))

    @staticmethod
    def _normalised_rod(value):
        if value is None:
            return None
        if hasattr(value, "tolist"):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return [float(v) for v in value]
        return float(value)

    @property
    def is_rod(self) -> bool:
        return self.rod_length is not None or self.rod_radius is not None

    @property
    def _nonsphere(self) -> bool:
        """Particles with an orientation: ellipsoids (aspect_ratios) or rods."""
        return bool(self.aspect_ratios) or self.is_rod

    def _rod_rows(self, n: int | None = None):
        """(lengths, radii) per particle, in the order of ``diameters``."""
        n = len(self.diameters) if n is None else n
        L = self.rod_length if isinstance(self.rod_length, list) else [self.rod_length] * n
        R = self.rod_radius if isinstance(self.rod_radius, list) else [self.rod_radius] * n
        return L, R

    def _rod_bounding_factor(self, L: float, R: float) -> float:
        """Largest extent of a rod in units of its D."""
        return L + 2.0 * R if self.rod_capped else 2.0 * math.sqrt(0.25 * L * L + R * R)

    def _check_rod_inputs(self) -> None:
        if self.aspect_ratios:
            raise ValueError("rod_length / rod_radius were given together with aspect_ratios: a packing is of rods "
                             "or of ellipsoids, not both")
        if self.rod_length is None or self.rod_radius is None:
            raise ValueError("rods need both rod_length and rod_radius (each one number or one per particle); "
                             "rod_radius=0.5 makes D the rod's thickness")
        n_part = len(self.diameters) if self.diameters else int(self.N)
        for name, val in (("rod_length", self.rod_length), ("rod_radius", self.rod_radius)):
            vals = val if isinstance(val, list) else [val]
            if isinstance(val, list) and len(val) != n_part:
                raise ValueError(f"{name} as a list needs one entry per particle ({n_part}), got {len(val)}")
            for i, v in enumerate(vals):
                who = f" (particle {i})" if isinstance(val, list) else ""
                if not math.isfinite(v):
                    raise ValueError(f"{name}{who} must be finite")
                if name == "rod_radius" and not v > 0.0:
                    raise ValueError(f"rod_radius{who} must be positive")
                if name == "rod_length" and (v < 0.0 or (not self.rod_capped and not v > 0.0)):
                    raise ValueError(f"rod_length{who} must be " + ("zero or positive (0 is a sphere)" if self.rod_capped
                                     else "positive: an uncapped rod of length 0 has no volume"))
        if not self.rod_capped and int(self.Ndim) > 3:
            raise NotImplementedError("uncapped (flat-ended) rods are implemented in 2D and 3D only; capped rods work "
                                      "in any dimension")

    @staticmethod
    def _normalised_ratios(aspect_ratios):
        """One list of floats (the same shape for every particle) or a list of N such lists (a shape per particle)."""
        if aspect_ratios is None:
            return []
        if hasattr(aspect_ratios, "tolist"):
            aspect_ratios = aspect_ratios.tolist()
        ratios = list(aspect_ratios)
        if ratios and (isinstance(ratios[0], (list, tuple)) or hasattr(ratios[0], "__len__")):
            return [[float(a) for a in row] for row in ratios]
        return [float(a) for a in ratios]

    @property
    def is_ellipsoid(self) -> bool:
        return bool(self.aspect_ratios)

    @property
    def has_shape_table(self) -> bool:
        """True when ``aspect_ratios`` is a table with one row per particle."""
        return bool(self.aspect_ratios) and isinstance(self.aspect_ratios[0], list)

    def _shape_rows(self, n: int | None = None) -> list:
        """The shape of every particle: N rows of Ndim-1 ratios, in the order of ``diameters``."""
        if self.has_shape_table:
            return self.aspect_ratios
        n = len(self.diameters) if n is None else n
        return [self.aspect_ratios] * n

    def _all_ratios_one(self) -> bool:
        rows = self.aspect_ratios if self.has_shape_table else [self.aspect_ratios]
        return all(float(a) == 1.0 for row in rows for a in row)

    def _reach_factor(self) -> float:
        """Largest extent of a particle in units of D: 1 for ellipsoids (D is the longest axis); a rod's tip-to-tip
        (capped) or corner-to-corner (uncapped) length."""
        if self.is_rod and self.rod_length is not None and self.rod_radius is not None:
            L, R = self._rod_rows(max(1, len(self.diameters)))
            return max(self._rod_bounding_factor(a, b) for a, b in zip(L, R))
        return 1.0

    _ORDER_RULE = ("aspect ratios are listed in decreasing order, 1 >= ratio_2 >= ratio_3 >= ... > 0: D is the "
                   "LONGEST axis and axis k has length D x ratio_k (3D prolate [0.5, 0.5], oblate [1, 0.5]). The "
                   "input is not reordered for you, because that would rescale D and change what the size "
                   "distribution means: divide D by the largest ratio and list the ratios longest first")

    def _check_shape_rows(self) -> None:
        nd = int(self.Ndim)
        table = self.has_shape_table
        rows = self.aspect_ratios if table else [self.aspect_ratios]
        n_part = len(self.diameters) if self.diameters else int(self.N)
        if table and len(rows) != n_part:
            raise ValueError(f"a per-particle aspect_ratios table needs one row per particle ({n_part}), got {len(rows)}")
        for i, row in enumerate(rows):
            who = f" (particle {i})" if table else ""
            if len(row) != nd - 1:
                raise ValueError(f"aspect_ratios{who} needs Ndim-1 = {nd - 1} entries (axis k has length D x ratio_k; "
                                 f"the first axis is D), got {len(row)}")
            if any(not (a > 0.0) or not math.isfinite(a) for a in row):
                raise ValueError(f"every aspect ratio{who} must be a positive finite number")
            prev = 1.0
            for a in row:
                if a > prev:
                    raise ValueError(f"aspect_ratios{who} = {list(row)} breaks the order rule: {self._ORDER_RULE}")
                prev = a

    def _nested_orientations(self, value, check: bool = False):
        """Orientations as N nested Ndim x Ndim matrices; accepts nested matrices or flat row-major lists."""
        nd = int(self.Ndim)
        out = []
        for item in (value if value is not None else []):
            rows = [list(map(float, row)) for row in item] if (len(item) and isinstance(item[0], (list, tuple))
                    or (hasattr(item, "ndim") and getattr(item, "ndim") == 2)) else None
            if rows is None:
                flat = [float(v) for v in item]
                if len(flat) != nd * nd:
                    raise ValueError(f"each orientation needs Ndim*Ndim = {nd * nd} numbers (a rotation matrix), "
                                     f"got {len(flat)}")
                rows = [flat[a * nd:(a + 1) * nd] for a in range(nd)]
            if len(rows) != nd or any(len(r) != nd for r in rows):
                raise ValueError(f"each orientation must be an Ndim x Ndim = {nd} x {nd} rotation matrix")
            out.append(rows)
        if check and len(out) != len(self.diameters):
            raise ValueError(f"orientations needs one rotation per particle ({len(self.diameters)}), got {len(out)}")
        return out

    def _flat_orientations(self) -> list:
        return [[v for row in R for v in row] for R in self.orientations]

    def _check_ellipsoid_supported(self) -> None:
        """Features not yet extended to ellipsoids refuse aspect ratios; they never run on the sphere assumption."""
        if not self._nonsphere:
            return
        if self.is_rod:
            self._check_rod_inputs()
        else:
            self._check_shape_rows()
        what = "rods (rod_length / rod_radius)" if self.is_rod else "aspect_ratios (ellipsoids)"

        def refuse(feature: str) -> None:
            raise NotImplementedError(f"{feature} is not supported with {what} yet; "
                                      f"remove them or do not use it")
        if float(self.swap) > 0.0 and not self.diameter_permute:
            refuse("the internal swap_* permutation settings (use diameter_permute=True)")
        if self.morph_target:
            refuse("the distribution morph (morph_target)")

    def shrink_to_contact(self) -> dict[str, Any]:
        """Exact shrink-to-contact of the current configuration (spheres or ellipsoids; flat or periodic box).

        Returns ``scale`` (the largest uniform factor <= 1 by which every particle can be scaled about its
        centre so that no pair and no wall overlaps, computed from the converged contact function over all
        pairs, less a 1e-12 safety margin), ``phi`` (the packing fraction at that scale), ``n_overlapping``
        (contacts that overlapped before the shrink) and ``pair`` (the tightest pair). The packing itself is
        not modified: multiply ``diameters`` by ``scale`` to get the non-overlapping configuration."""
        cs = _contact_scale(self._collect_pack_input(), self._collect_pack_config(translate=False))
        scale = float(cs["scale"])
        if scale < 1.0:
            scale *= (1.0 - 1e-12)
        return {"scale": scale, "phi": self._current_phi() * scale ** int(self.Ndim),
                "n_overlapping": int(cs["n_overlapping"]), "pair": (int(cs["pair_i"]), int(cs["pair_j"]))}

    def _init_extension_state(self, config: dict[str, Any], *, consume: bool) -> None:
        """Set the fields added after v1.0 (diameter permutation, morph, quench, internal permutation fields)
        and reset their diagnostics. Shared by __init__ (config = the keyword arguments, consumed) and from_dict
        (config = a saved pack config; absent keys - e.g. a dict saved by v1.0 - take their defaults)."""
        take = config.pop if consume else config.get
        for k, v in self._SWAP_DEFAULTS.items():
            setattr(self, k, type(v)(take(k, v)))
        for k, v in self._PERMUTE_DEFAULTS.items():
            setattr(self, k, type(v)(take(k, v)))
        for k, v in self._MORPH_DEFAULTS.items():
            val = take(k, v)
            setattr(self, k, [float(t) for t in val] if k == "morph_target" else type(v)(val))
        for k, v in self._QUENCH_DEFAULTS.items():
            val = take(k, v)
            setattr(self, k, [float(t) for t in val] if k == "quench_held" else type(v)(val))
        self.quench_residual = 0.0; self.quench_global_residual = 0.0; self.quench_iters = 0; self.quench_converged = 0
        self.quench_log = []; self.quench_held_force = []; self.quench_status = ""; self.quench_message = ""
        self.morph_log = []         # [219] per-update morph diagnostics: 6 values per update [step, lambda, increment, max_overlap, held, sum D0^3]
        self.morph_lambda = 0.0
        self.morph_updates = 0
        self.morph_holds = 0
        self.morph_stalled = 0
        self.morph_rate_ceiling = 0.0
        self.swap_log = []          # per-sweep event log from the core (8 values per sweep)
        self.particle_origin = []   # slot k -> original particle index (through spatial reorders)
        self.swap_counts = []
        self.swap_relayouts = 0
        self.swap_rescue_candidates = 0
        self.swap_rescue_accepted = 0
        self.swap_rescue_false_rej = 0
        self.swap_rescue_false_acc = 0
        self.swap_cycle3_candidates = 0
        self.swap_cycle3_accepted = 0
        self.swap_cycle3_blocked = 0
        self.capped = False
        self.resumed_from_step = 0
        self.checkpoint_last_step = 0
        self.swap_events = []
        self.swap_event_width = 8

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
        if name in ("rod_length", "rod_radius"):   # new rod numbers from the user: what initialize() will use
            object.__setattr__(self, name, self._normalised_rod(value))
            object.__setattr__(self, "_rod_init", (deepcopy(self.rod_length), deepcopy(self.rod_radius)))
        if name == "aspect_ratios":   # a new shape (or table) from the user: it is what initialize() will use
            object.__setattr__(self, "aspect_ratios", self._normalised_ratios(value))
            object.__setattr__(self, "_aspect_init", deepcopy(self.aspect_ratios))
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
            "seed": self.seed,   # deterministic, seed-varying initial configuration
            **({"aspect_ratios": deepcopy(self._aspect_init)} if self.aspect_ratios else {}),
            **({"rod_length": deepcopy(self._rod_init[0]), "rod_radius": deepcopy(self._rod_init[1]),
                "rod_capped": self.rod_capped} if self.is_rod else {}),
        }

    def _collect_pack_input(self) -> dict[str, Any]:
        state = {
            "positions": deepcopy(self.positions),
            "diameters": deepcopy(self.diameters),
        }
        if self.is_rod:   # rods: length and radius are part of the state (one number or one per particle)
            state["orientations"] = self._flat_orientations()
            state["rod_length"] = list(self.rod_length) if isinstance(self.rod_length, list) else [self.rod_length]
            state["rod_radius"] = list(self.rod_radius) if isinstance(self.rod_radius, list) else [self.rod_radius]
        if self.aspect_ratios:
            state["orientations"] = self._flat_orientations()
            if self.has_shape_table:   # a shape per particle is part of the state, like the diameters
                state["aspect_ratios"] = deepcopy(self.aspect_ratios)
        return state

    def _collect_pack_config(self, translate: bool = True) -> dict[str, Any]:
        cfg = {
            "box": deepcopy(self.box),
            "walls": deepcopy(self.walls),
            "neighbor_max": self.neighbor_max,
            "seed": self.seed,
            "fix_height": self.fix_height,
            **{k: getattr(self, k) for k in self._SWAP_DEFAULTS},
            **{k: getattr(self, k) for k in self._PERMUTE_DEFAULTS},
            **{k: deepcopy(getattr(self, k)) for k in self._MORPH_DEFAULTS},
            **{k: deepcopy(getattr(self, k)) for k in self._QUENCH_DEFAULTS},
        }
        if self.aspect_ratios:   # ellipsoids only: a sphere pack config is exactly what it always was
            if not self.has_shape_table:   # one shape for all travels in the config, a table in the input
                cfg["aspect_ratios"] = list(self.aspect_ratios)
            cfg.update({k: getattr(self, k) for k in self._ELL_DEFAULTS})
        if self.is_rod:
            cfg["rod_capped"] = int(self.rod_capped)
            cfg.update({k: getattr(self, k) for k in self._ELL_DEFAULTS})
        if translate and self.diameter_permute:   # to_dict() keeps the user's settings (translate=False)
            cfg.update(self._permute_recipe())
        return cfg

    def _permute_recipe(self) -> dict[str, Any]:
        """Internal swap_* settings for ``diameter_permute=True``: capacity-matched
        exchanges within a 5 % rank window plus 3-cycles, quota-controlled, from
        ``permute_start`` on, with ``permute_budget * N`` exchanges per pass. This is the
        recipe validated at S = 80 and S = 100 (research log [149]/[150]); the two public
        knobs are the only things a user needs to vary."""
        n = len(self.diameters) if self.diameters else int(self.N)
        budget = float(self.permute_budget)
        if not (0.0 < budget <= 0.5):
            raise ValueError("permute_budget must lie in (0, 0.5]")
        if int(self.permute_start) < 0:
            raise ValueError("permute_start must be >= 0")
        if int(self.permute_interval) < 1:
            raise ValueError("permute_interval must be >= 1 (steps between permutation passes)")
        return {
            "swap": 1.0, "swap_ctrl": 1, "swap_a": 0.3, "swap_R": 0.05, "swap_sweeps_off": 1,
            "swap_geo": 1, "swap_geo_T": 0.0, "swap_geo_w": 0.05, "swap_geo_k": int(self.permute_interval),
            "swap_geo_cycle": 1, "swap_geo_B": max(1, int(budget * n)), "swap_geo_start": int(self.permute_start),
        }

    # ---- diameter-permutation diagnostics (public names; the underlying fields keep their legacy prefix) ----
    @property
    def permute_counts(self) -> list:
        """Accepted diameter exchanges per particle (empty when permutation is off)."""
        return self.swap_counts

    @property
    def permute_passes(self) -> int:
        """Number of permutation passes run."""
        return int(getattr(self, "swap_geo_passes", 0) or 0)

    @property
    def permute_exchanges(self) -> int:
        """Total accepted diameter exchanges (pairwise + 3-cycle rotations)."""
        return int(getattr(self, "swap_geo_trades", 0) or 0) + int(getattr(self, "swap_geo_cycle_trades", 0) or 0)

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
        self.swap_log = []
        self.swap_counts = []
        self.swap_relayouts = 0
        self.swap_rescue_candidates = 0
        self.swap_rescue_accepted = 0
        self.swap_rescue_false_rej = 0
        self.swap_rescue_false_acc = 0
        self.swap_cycle3_candidates = 0
        self.swap_cycle3_accepted = 0
        self.swap_cycle3_blocked = 0
        self.capped = False
        self.resumed_from_step = 0
        self.checkpoint_last_step = 0
        self.particle_origin = []
        self.steps = 0
        self.phi_final = None
        self.max_min_dist = None
        self.force_magnitude = None
        self.phi_history = []
        self.force_history = []
        self.energy_history = []
        self.max_overlap_history = []
        self.mu_flag_history = []
        self.mu_history = []
        self.alpha_history = []
        self.phi_contact = None
        self.shrink_factor = None
        self.torque_magnitude = None

    def _clear_trajectory(self) -> None:
        self.trajectory_positions = []
        self.trajectory_diameters = []
        self.trajectory_steps = []
        self.trajectory_phi = []
        self.trajectory_force = []
        self.trajectory_energy = []
        self.trajectory_max_min_dist = []
        self.trajectory_orientations = []
        self.trajectory_origin = []

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
        if self.aspect_ratios:   # ellipsoid volume = sphere volume for D x product of the ratios
            particle_volume = 0.0
            for diameter, row in zip(self.diameters, self._shape_rows()):
                particle_volume += coeff * (diameter / 2.0) ** self.Ndim * math.prod(row)
        if self.is_rod:   # straight section x cross-section, plus one ball of the rod's radius when capped
            nd = int(self.Ndim)
            c1 = math.pi ** ((nd - 1) / 2.0) / math.gamma((nd - 1) / 2.0 + 1.0)
            particle_volume = 0.0
            L, R = self._rod_rows()
            for diameter, a, b in zip(self.diameters, L, R):
                rad = diameter * b
                particle_volume += c1 * rad ** (nd - 1) * diameter * a + (coeff * rad ** nd if self.rod_capped else 0.0)
        return particle_volume / box_volume / phi_modifier

    def _mark_realized_state_updated(self) -> None:
        self._clear_results()
        self._clear_trajectory()
        self._mark_needs_pack()

    def _clear_realized_state(self) -> None:
        self.positions = []
        self.diameters = []
        self.orientations = []
        self._clear_results()
        self._clear_trajectory()

    def _user_list_order(self, n: int):
        """Which entry of the user's own diameter list each particle carries, or None.

        The initializer sorts particles by decreasing diameter. When the code GENERATES the diameters (any
        ``dist`` but "custom") there is no user order: row i of a per-particle shape table or rod list goes to the
        i-th largest particle. When the user SUPPLIES the diameters (``dist = {"type": "custom", "custom": [...]}``)
        entry i of a shape table or rod list belongs to entry i of that diameter list and stays with it through
        the sort: particle k (the k-th largest) carries the user's entry ``order[k]``. Equal diameters keep the
        user's order among themselves."""
        dist = self.dist if isinstance(self.dist, dict) else {}
        custom = dist.get("custom") if dist.get("type") == "custom" else None
        if custom is None or len(custom) != n:
            return None
        return sorted(range(n), key=lambda i: -float(custom[i]))

    def _apply_initializer_result(self, result: dict[str, Any]) -> None:
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            self.positions = result["positions"]
            self.diameters = result["diameters"]
            self.box = result["box"]
            self.walls = result["walls"]
            order = self._user_list_order(len(self.diameters))
            if self.is_rod:
                self.rod_length, self.rod_radius = deepcopy(self._rod_init[0]), deepcopy(self._rod_init[1])
                if order is not None:   # the user's own diameter list: entry i of a rod list goes with entry i of it
                    if isinstance(self.rod_length, list):
                        self.rod_length = [self.rod_length[k] for k in order]
                    if isinstance(self.rod_radius, list):
                        self.rod_radius = [self.rod_radius[k] for k in order]
                self.orientations = self._nested_orientations(
                    _random_orientations(len(self.diameters), int(self.Ndim), int(self.seed)))
                if self.fix_height:
                    # fix_height: the height is a multiple of the largest particle EXTENT. The initializer set it
                    # from the largest D; a rod's extent is its bounding diameter D x (tip to tip, or the diagonal),
                    # so the last axis is stretched by (largest extent) / (largest D), positions with it. The
                    # packing loop then keeps height / size fixed as every particle grows by the same factor.
                    L, R = self._rod_rows()
                    extent = max(d * self._rod_bounding_factor(a, b) for d, a, b in zip(self.diameters, L, R))
                    stretch = extent / max(self.diameters)
                    self.box = list(self.box[:-1]) + [self.box[-1] * stretch]
                    self.positions = [list(x[:-1]) + [x[-1] * stretch] for x in self.positions]
            elif self.aspect_ratios:   # the initializer generates orientations, as it generates positions
                self.aspect_ratios = deepcopy(self._aspect_init)   # table: row i <-> particle i as initialised
                if order is not None and self.has_shape_table:   # the user's own diameter list: row i goes with entry i of it
                    self.aspect_ratios = [self.aspect_ratios[k] for k in order]
                self.orientations = self._nested_orientations(
                    _random_orientations(len(self.diameters), int(self.Ndim), int(self.seed)))
            else:
                self.orientations = []
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
            _cf = _container_fraction(result["walls"])   # curved container -> phi vs disk/sphere volume, not the box
            self.phi_final = result["phi"] / _cf
            self.max_min_dist = result["max_min_dist"]
            self.force_magnitude = result["force_magnitude"]
            self.phi_history = ([v / _cf for v in result["phi_history"]] if _cf != 1.0 else result["phi_history"])
            self.force_history = result["force_history"]
            self.energy_history = result["energy_history"]
            self.max_overlap_history = result.get("max_overlap_history", [])
            self.mu_flag_history = result.get("mu_flag_history", [])
            self.zpos_history = result.get("zpos_history", [])
            self.mu_eff_history = result.get("mu_eff_history", [])
            self.mu_hold_applied = int(result.get("mu_hold_applied", 0))
            self.mu_ctrl_final_mult = float(result.get("mu_ctrl_final_mult", 0.0))
            self.swap_geo_trades = int(result.get("swap_geo_trades", 0))
            self.swap_geo_passes = int(result.get("swap_geo_passes", 0))
            self.swap_geo_reopens = int(result.get("swap_geo_reopens", 0))
            self.swap_geo_cycle_trades = int(result.get("swap_geo_cycle_trades", 0))
            self.swap_geo_cycle4_trades = int(result.get("swap_geo_cycle4_trades", 0))
            self.swap_geo_capk_trades = int(result.get("swap_geo_capk_trades", 0))
            self.swap_geo_cage_trades = int(result.get("swap_geo_cage_trades", 0))
            self.swap_geo_bottleneck_trades = int(result.get("swap_geo_bottleneck_trades", 0))
            self.mu_history = result.get("mu_history", [])
            self.alpha_history = result.get("alpha_history", [])
            self.swap_log = result.get("reorder_log", [])
            self.particle_origin = result.get("particle_origin", [])
            self.swap_counts = result.get("spike_log", [])   # per-particle accepted swaps when swap > 0
            self.swap_relayouts = int(result.get("swap_relayouts", 0))
            self.swap_rescue_candidates = int(result.get("swap_rescue_candidates", 0))
            self.swap_rescue_accepted = int(result.get("swap_rescue_accepted", 0))
            self.swap_rescue_false_rej = int(result.get("swap_rescue_false_rej", 0))
            self.swap_rescue_false_acc = int(result.get("swap_rescue_false_acc", 0))
            self.swap_cycle3_candidates = int(result.get("swap_cycle3_candidates", 0))
            self.swap_cycle3_accepted = int(result.get("swap_cycle3_accepted", 0))
            self.swap_cycle3_blocked = int(result.get("swap_cycle3_blocked", 0))
            self.capped = bool(result.get("capped", False))
            self.resumed_from_step = int(result.get("resumed_from_step", 0))
            self.checkpoint_last_step = int(result.get("checkpoint_last_step", 0))
            self.morph_log = list(result.get("morph_log", []))          # [219]
            self.morph_lambda = float(result.get("morph_lambda", 0.0))
            self.morph_updates = int(result.get("morph_updates", 0))
            self.morph_holds = int(result.get("morph_holds", 0))
            self.morph_stalled = int(result.get("morph_stalled", 0))
            self.quench_residual = float(result.get("quench_residual", 0.0)); self.quench_global_residual = float(result.get("quench_global_residual", 0.0))
            self.quench_iters = int(result.get("quench_iters", 0)); self.quench_converged = int(result.get("quench_converged", 0))
            self.quench_log = list(result.get("quench_log", [])); self.quench_held_force = list(result.get("quench_held_force", []))
            self.quench_status = _QUENCH_STATUS.get(int(result.get("quench_status", 0)), "")
            self.quench_message = ""
            self.morph_rate_ceiling = float(result.get("morph_rate_ceiling", 0.0))
            self.swap_events = list(result.get("swap_events", []))
            self.swap_event_width = int(result.get("swap_event_width", 8))
            if self.is_rod:
                for key in ("rod_length", "rod_radius"):   # per-particle lists come back with their particles
                    if isinstance(getattr(self, key), list) and result.get(key):
                        setattr(self, key, [float(v) for v in result[key]])
            if self._nonsphere:   # ellipsoids and rods: orientations out, and the reported phi after the exact shrink
                self.orientations = self._nested_orientations(result.get("orientations", []))
                if self.has_shape_table and result.get("aspect_ratios"):   # shapes come back with their particles
                    self.aspect_ratios = [[float(a) for a in row] for row in result["aspect_ratios"]]
                self.torque_magnitude = float(result.get("torque_magnitude", 0.0))
                self.ell_capture_record = result.get("ell_capture")
                self.ell_permute_candidates = int(result.get("ell_permute_candidates", 0))
                shrunk = self.shrink_to_contact()
                self.shrink_factor = shrunk["scale"]
                self.phi_contact = shrunk["phi"]
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
            self.trajectory_orientations = [self._nested_orientations(o) for o in trace.get("orientations", [])]
            # which particle is in which slot of each frame (the core renumbers particles while it runs):
            # trajectory_origin[f][k] = index at the start of this run of the particle in slot k of frame f
            self.trajectory_origin = [list(o) for o in trace.get("origins", [])]
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
        trajectory_per_decade: int = 0,
    ) -> dict[str, Any]:
        self._clear_trajectory()
        self._check_ellipsoid_supported()
        if self._nonsphere and len(self.orientations) != len(self.diameters):
            raise ValueError("aspect_ratios (or rod inputs) are set but orientations do not match the particles: supply both "
                             "(ellipsoids) or neither (spheres); initialize() generates orientations")
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

            pack_config = self._collect_pack_config()
            if trajectory_per_decade > 0:      # frames evenly in log(step); the interval is then not used
                pack_config = dict(pack_config, trajectory_per_decade=int(trajectory_per_decade))
            observed = _run_packing_observed(
                self._collect_pack_input(),
                pack_config,
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
        self._warn_small_periodic_box()
        return result

    def _min_periodic_box(self, dmax: float) -> float:
        """Smallest periodic box length, in largest diameters, below which a pair can be counted wrongly.
        2D/3D: two particles touch through more than one periodic copy when L < D_i + D_j <= 2 Dmax (forces use the
        nearest copy only); 2.1 leaves 5 % for overlaps. 4D and up: the neighbour search can also store a pair
        twice (double-counted force) once its reach, about 1.01 Dmax plus the neighbour skin, passes half the box."""
        if int(self.Ndim) <= 3:
            return 2.1
        d = sorted(self.diameters)
        skin = 0.5 * min(max(1.02 * d[len(d) // 2], 2.25 * d[0]), d[-1])   # the core's neighbour-list shell
        return 2.0 * (1.01 + skin / dmax)

    def _warn_small_periodic_box(self) -> None:
        if not self.diameters or not self.box or not self.walls:
            return
        dmax = max(self.diameters)
        if self._nonsphere:   # an elongated particle reaches further than its D
            dmax *= self._reach_factor()
        if not dmax > 0:
            return
        w0 = int(self.walls[0])
        t = (2 if w0 == -1 else -w0) if w0 < 0 else 0      # axes inside a curved container are not periodic
        need = self._min_periodic_box(dmax)
        short = [(k, self.box[k] / dmax) for k in range(t, len(self.box))
                 if int(self.walls[k]) == 0 and self.box[k] < need * dmax]
        if not short:
            return
        axes = ", ".join(f"axis {k}: {r:.2f}" for k, r in short)
        warnings.warn(
            f"periodic box is small for the largest particle (box length / largest diameter = {axes}; needs at "
            f"least {need:.2f} here). Particles can then interact through more than one periodic copy, which the "
            f"forces do not account for. Use more particles, or hard walls (walls entry 1) on that axis.",
            UserWarning, stacklevel=4)

    _WALLS_HELP = ("one entry per axis: 0 = periodic, 1 = hard walls at both ends (e.g. [0, 0, 1]); a curved "
                   "container is a negative FIRST entry over the first t axes: [-2, 0] disk, [-2, 0, 0] cylinder "
                   "along z, [-2, 0, 1] capped cylinder, [-3, 0, 0] sphere. See rcpgenerator.example(), "
                   "'Boundary Conditions Examples'.")

    def _validated_walls(self) -> list[int]:
        """Check `walls` against Ndim and return it as integers. Everything the core handles correctly is
        accepted unchanged (including -1 as an alias of -2 and any subordinate entry of a curved container);
        forms the core would silently turn into an inconsistent boundary raise ValueError."""
        nd = int(self.Ndim)
        if len(self.box) != nd:
            raise ValueError(f"box={list(self.box)} needs exactly Ndim={nd} lengths")
        w = list(self.walls) if self.walls is not None else []
        if len(w) != nd:
            raise ValueError(f"walls={w} has {len(w)} entries but Ndim={nd}: {self._WALLS_HELP}")
        out = []
        for v in w:
            if isinstance(v, bool) or not isinstance(v, numbers.Real) or float(v) != int(v):
                raise ValueError(f"walls={w}: entries must be integers; {self._WALLS_HELP}")
            out.append(int(v))
        t = 0
        if out[0] < 0:
            t = 2 if out[0] == -1 else -out[0]
            if t > nd:
                raise ValueError(f"walls={w}: walls[0]={out[0]} asks for a {t}-dimensional container in {nd} "
                                 f"dimensions; {self._WALLS_HELP}")
        for k in range(t, nd):
            if out[k] not in (0, 1):
                where = "a curved container must be the FIRST entry" if out[k] < 0 else f"walls[{k}]={out[k]} is not 0 or 1"
                raise ValueError(f"walls={w}: {where}; {self._WALLS_HELP}")
        for k in range(1, t):
            if out[k] not in (-1, 0, 1):
                raise ValueError(f"walls={w}: walls[{k}] lies inside the curved container (ignored there) and "
                                 f"must be 0; {self._WALLS_HELP}")
        return out

    def initialize(self) -> dict[str, Any]:
        self._log("rcpgenerator: initializing particles")
        config = self._collect_init_config()
        config["walls"] = self._validated_walls()
        self._check_ellipsoid_supported()
        result = _initialize_particles(config)
        self._apply_initializer_result(result)
        self._log("rcpgenerator: initialization complete")
        return result

    def pack(
        self,
        verbose: bool | None = None,
        progress_interval: int = 1000,
        capture_trajectory: bool = False,
        trajectory_interval: int = 1000,
        trajectory_spacing: str = "linear",
        trajectory_per_decade: int = 20,
    ) -> dict[str, Any]:
        """Run the packing. ``capture_trajectory=True`` records frames for ``animate_2d``: every
        ``trajectory_interval`` steps (``trajectory_spacing="linear"``), or evenly in the LOGARITHM of the step,
        ``trajectory_per_decade`` frames per factor of 10 in steps (``trajectory_spacing="log"``). Most of the
        motion happens in the first few hundred steps of some 16,000, so log spacing gives a movie in which
        things happen throughout; linear spacing gives one that is almost still after the first frames."""
        if trajectory_spacing not in ("linear", "log"):
            raise ValueError("trajectory_spacing must be 'linear' or 'log'")
        if trajectory_spacing == "log" and int(trajectory_per_decade) < 1:
            raise ValueError("trajectory_per_decade must be at least 1")
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
            trajectory_per_decade=int(trajectory_per_decade) if (capture_trajectory and trajectory_spacing == "log") else 0,
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

    def quench(
        self,
        *,
        tol: float = 1e-10,
        max_iter: int = 200000,
        kappa_mult: float | None = None,
        resume: bool = False,
        held=None,
        stall_window: int | None = None,
        verbose: bool | None = None,
        mass_scaled: bool = True,
        dt0: float = 0.05,
        dt_max: float | None = None,
        alpha0: float = 0.03,
        mix_mode: int = 1,
        report_every: int = 500,
        rod_round: float | None = None,
    ) -> dict[str, Any]:
        """Drive the current packing to mechanical force balance at fixed diameters.

        With ``aspect_ratios`` (ellipsoids) the same minimisation moves positions and orientations together and
        the residual measures force and torque balance: |(Σf, Στ/a)| / Σ|f| per load-bearing particle, with ``a``
        the particle's largest semi-axis. ``orientations`` come back in the caller's order like the positions.

        A FIRE minimisation of the overlap energy replaces the packing update: every
        diameter is first multiplied once by ``kappa_mult`` and then frozen, the μ schedule
        is suspended and no diameter permutation runs. The packing is modified in place,
        whether or not the quench converges, so the state it stopped in is always
        available afterwards; particle ``i`` keeps index ``i`` (positions, diameters,
        ``particle_origin`` and ``permute_counts`` stay aligned).

        Parameters
        ----------
        tol : float
            Convergence threshold on the largest per-particle residual |Σf| / Σ|f| over
            the load-bearing particles (Ndim + 1 or more contacts, wall contacts
            included, and at least 1 % of the typical contact load).
        max_iter : int
            Iteration cap for this call.
        kappa_mult : float, optional
            One-time rescaling of all diameters before the quench. A packing run ends at
            the edge of jamming; at exactly its own diameters it carries no load and cannot
            be balanced, so a small compression is needed. Default for a fresh quench: 1.003
            (0.3 % larger diameters) in 2D and 3D, which converges for size ratios up to
            about 10; 1.01 in 4D and above, where broad size distributions need the larger
            compression to converge reliably, and 1.01 for ellipsoids in any dimension (their
            packings have fewer contacts than degrees of freedom and are held stable by the
            compression itself). Not applied when ``resume=True``.
        resume : bool
            Continue from the state a previous ``quench()`` left (for example with a
            tighter ``tol`` or after a stop): the diameters are not rescaled again.
        held : sequence of int, optional
            Indices (into ``positions`` as they are now) of particles that must not move;
            the net force on them is reported in ``quench_held_force``.
        rod_round : float, optional
            Rods only. The quench balances forces and torques on slightly ROUNDED rods (edges rounded by this
            relative amount, default 0.01): a sharp rod's contact energy has kinks where two flat faces rest on
            each other, and no balance point exists for a descent method to find. 0 uses the sharp rods (the
            quench then usually stalls). The independent check of the result must use the same value
            (``evaluate_once`` config key ``rod_round``).
        stall_window : int, optional
            Give up when the best residual has not halved within this many iterations.
            Default ``max(20000, 20 * N)``; 0 disables the check. With ellipsoids a fall of the
            energy by more than 0.1 % within the window also counts as progress (slow rearrangements
            lower the energy for tens of thousands of iterations before the residual drops).
        verbose : bool, optional
            Print the outcome (defaults to the packing's ``verbose``).
        mass_scaled, dt0, dt_max, alpha0, mix_mode, report_every
            Advanced FIRE settings; the defaults are the measured best. ``mix_mode`` 0 uses
            global velocity mixing (Bitzek et al.), 1 per-particle mixing.

        Returns
        -------
        dict
            The run result, including ``quench_status`` and ``quench_message``.

        After the call, ``quench_status`` is one of ``"converged"``, ``"max_iter"``
        (cap reached while still improving), ``"stalled"``, ``"unjammed"`` (no
        load-bearing particles left), ``"step_collapse"`` or ``"non_finite"``, and
        ``quench_message`` explains it. Anything other than ``"converged"`` also raises a
        ``RuntimeWarning``. ``quench_residual``, ``quench_global_residual``,
        ``quench_iters`` and ``quench_log`` (per report: iteration, max residual, global
        residual, energy, load-bearing count, time step) hold the details.
        """
        if self._needs_initialize or not self._initialized or not self.positions:
            raise RuntimeError("quench() needs a realized state: call pack() first")
        if resume and kappa_mult is not None and float(kappa_mult) != 1.0:
            raise ValueError("quench(resume=True) continues at the current diameters; kappa_mult cannot be applied "
                             "(call quench(kappa_mult=...) without resume to rescale)")
        default_kmult = 1.003 if int(self.Ndim) <= 3 else 1.01   # [228] 4D+: x1.003 gives up in ~40 % of broad-distribution runs
        if self.is_rod or (self.aspect_ratios and not self._all_ratios_one()):
            # ellipsoids: x1.003 leaves 2D ratio 3, 3D ratio 1.05 and the 1.25:1:0.8 shape wandering at a residual
            # of 1e-3 to 1e-2 for 200,000 iterations; x1.01 converges them (EllipsoidalPacking LOG [020], [021])
            default_kmult = 1.01
            ell_shaped = True
        else:
            ell_shaped = False
        if dt_max is None:
            # spheres 0.5 as before. Ellipsoids 0.2: at 0.5 the step keeps overshooting and being cut back on the
            # more elongated or flattened shapes (2D ratio 4, 3D ratio 3 and 0.33 did not converge); at 0.2 it holds
            # (EllipsoidalPacking LOG [022])
            dt_max = 0.2 if ell_shaped else 0.5
        kmult = 1.0 if resume else (default_kmult if kappa_mult is None else float(kappa_mult))
        if not (tol > 0.0) or int(max_iter) < 1 or not (kmult > 0.0):
            raise ValueError("quench(): tol and kappa_mult must be > 0 and max_iter >= 1")
        if int(mix_mode) not in (0, 1):
            raise ValueError("quench(): mix_mode must be 0 (global mixing) or 1 (per-particle mixing)")
        n = len(self.diameters)
        window = max(20000, 20 * n) if stall_window is None else int(stall_window)
        if window < 0:
            raise ValueError("quench(): stall_window must be >= 0 (0 disables the stall check)")
        held_idx = [] if held is None else [int(i) for i in held]
        if any(i < 0 or i >= n for i in held_idx):
            raise ValueError("quench(): held indices must lie in [0, N)")
        saved = {k: deepcopy(getattr(self, k)) for k in self._QUENCH_DEFAULTS}
        saved["checkpoint_path"] = self.checkpoint_path
        saved["diameter_permute"] = self.diameter_permute
        saved["swap"] = self.swap
        # the pack run's permutation diagnostics and particle identity map survive the quench
        kept = {k: deepcopy(getattr(self, k, 0)) for k in ("swap_counts", "particle_origin", "swap_geo_trades",
                                                            "swap_geo_passes", "swap_geo_cycle_trades")}
        try:
            # the quench acts on the realized state held here, never on a checkpoint: with checkpoint_path set the
            # run would resume the blob (dropping kappa_mult) and overwrite it at exit
            self.checkpoint_path = ""
            diameters_before = list(self.diameters)
            run_finished = False
            if self.is_rod and int(self.Ndim) > 3 and (0.01 if rod_round is None else float(rod_round)) > 0.0:
                raise NotImplementedError("quench() of rods uses the rounded-rod contact, which is implemented in 2D and "
                                          "3D only")
            if kmult != 1.0:
                self.diameters = [float(d) * kmult for d in self.diameters]
            for k, v in self._QUENCH_DEFAULTS.items():   # research-only fields stay at their (off) defaults
                setattr(self, k, deepcopy(v))
            self.quench_mode = 1
            self.quench_tol = float(tol)
            self.quench_max_iter = int(max_iter)
            self.quench_dt0 = float(dt0)
            self.quench_dt_max = float(dt_max)
            self.quench_alpha0 = float(alpha0)
            self.quench_mass_scaled = int(bool(mass_scaled))
            self.quench_report_every = int(report_every)
            self.quench_kappa_mult = 1.0          # applied above, on the diameters themselves
            self.quench_mix_mode = int(mix_mode)
            self.quench_stall_window = window
            self.quench_dt_min = 1e-8
            # rods: the quench works on slightly rounded rods (smooth contact; see rod.hpp), default rounding 1 %
            self.quench_rod_round = (0.01 if rod_round is None else float(rod_round)) if self.is_rod else 0.0
            self.quench_held = [float(i) for i in held_idx]
            self.diameter_permute = False          # frozen diameters: no permutation during the quench
            self.swap = 0.0
            effective_verbose = self.verbose if verbose is None else verbose
            result = self._execute_packing(
                verbose=effective_verbose,
                progress_interval=max(1, int(report_every)),
                capture_trajectory=False,
                trajectory_interval=0,
                options={"max_steps": int(max_iter) + 1, "fix_diameter": True},
            )
            run_finished = True
        finally:
            object.__setattr__(self, "_mutating_internal_state", True)
            try:
                if not run_finished and "diameters_before" in locals():   # refused or failed: nothing was done
                    self.diameters = diameters_before
                for k, v in saved.items():
                    setattr(self, k, v)
            finally:
                object.__setattr__(self, "_mutating_internal_state", False)
        self._restore_order_and_diagnostics(kept, result)
        message = self._quench_message(tol=float(tol), max_iter=int(max_iter), window=window)
        object.__setattr__(self, "quench_message", message)
        result["quench_status"] = self.quench_status
        result["quench_message"] = message
        if effective_verbose:
            print(f"rcpgenerator: quench {message}")
        if self.quench_status != "converged":
            warnings.warn(f"quench() did not converge: {message}", RuntimeWarning, stacklevel=2)
        return result

    def _restore_order_and_diagnostics(self, kept: dict[str, Any], result: dict[str, Any]) -> None:
        """After a quench: particle i keeps index i. The core reorders particles spatially during any run; the
        quench's own particle_origin (final slot -> slot at quench start) undoes that, so positions and diameters
        come back in the caller's order and the pack run's particle_origin and permutation diagnostics stay valid."""
        q_origin = list(self.particle_origin)
        n = len(self.diameters)
        object.__setattr__(self, "_mutating_internal_state", True)
        try:
            if len(q_origin) == n and sorted(q_origin) == list(range(n)):
                pos, dia = [None] * n, [0.0] * n
                for k, i in enumerate(q_origin):
                    pos[i] = self.positions[k]
                    dia[i] = self.diameters[k]
                self.positions, self.diameters = pos, dia
                result["positions"], result["diameters"] = pos, dia
                if self.has_shape_table and len(self.aspect_ratios) == n:   # so does a per-particle shape
                    tab = [None] * n
                    for k in range(n):
                        tab[q_origin[k]] = self.aspect_ratios[k]
                    self.aspect_ratios = tab
                    result["aspect_ratios"] = tab
                for key in ("rod_length", "rod_radius"):   # rods: per-particle lists follow their particle too
                    val = getattr(self, key, None)
                    if isinstance(val, list) and len(val) == n:
                        lst = [0.0] * n
                        for k in range(n):
                            lst[q_origin[k]] = val[k]
                        setattr(self, key, lst)
                        result[key] = lst
                if self._nonsphere and len(self.orientations) == n:   # orientations follow their particle
                    ori = [None] * n
                    for k, i in enumerate(q_origin):
                        ori[i] = self.orientations[k]
                    self.orientations = ori
                    result["orientations"] = [[v for row in R for v in row] for R in ori]
            for k, v in kept.items():
                setattr(self, k, v)
            result["particle_origin"] = self.particle_origin
        finally:
            object.__setattr__(self, "_mutating_internal_state", False)

    def _quench_message(self, *, tol: float, max_iter: int, window: int) -> str:
        """One-paragraph explanation of how the last quench ended, with the likely cause."""
        st, res, it = self.quench_status, self.quench_residual, self.quench_iters
        d = [float(v) for v in self.diameters if v > 0]
        S = (max(d) / min(d)) if d else 1.0
        log = self.quench_log
        e0, e1 = (log[3], log[-3]) if len(log) >= 12 else (0.0, 0.0)
        size_note = (f" The size ratio here is {S:.3g}; above about 20 this method often cannot find a balanced "
                     f"state of the packing (see README, 'Quench to force balance').") if S > 15.0 else ""
        if st == "converged":
            return f"converged: max residual {res:.2e} <= tol {tol:.0e} after {it} iterations."
        if st == "max_iter":
            return (f"stopped at max_iter={max_iter} with max residual {res:.2e} (tol {tol:.0e}) while still improving; "
                    f"call quench(resume=True) to continue from here.{size_note}")
        if st == "stalled":
            if res < 1e-8:
                return (f"stopped at the resolution floor: the max residual has not halved in {window} iterations "
                        f"(now {res:.2e}, tol {tol:.0e}). Residuals this small are limited by floating-point rounding "
                        f"of the positions (roughly 6e-14 divided by the mean relative overlap), so a tighter tol is "
                        f"not reachable at this compression; the state is balanced to the precision available.")
            if e0 > 0.0 and e1 < 1e-2 * e0:
                cause = (" The overlap energy fell by a factor of %.0e, so the packing is relaxing apart rather than "
                         "balancing: it is under-compressed. Quench the original state with a larger kappa_mult." % (e0 / max(e1, 1e-300)))
            else:
                cause = (" Quench the ORIGINAL (packed) state again with a lower dt_max (for example half the present one; much lower "
                         "than needed costs many more iterations) or with slightly more compression (a larger kappa_mult); when a "
                         "lower dt_max does not help, more compression is the next thing to try. "
                         "Resuming with unchanged settings does not help unless the residual was still falling.")
                if self.is_rod and not self.rod_capped and int(self.Ndim) == 3 and self.walls and int(self.walls[0]) < 0:
                    cause += " (Flat-ended rods in a 3D curved container often need this retry.)"
            return (f"gave up: the max residual has not halved in {window} iterations (now {res:.2e}, tol {tol:.0e})."
                    f"{cause}{size_note}")
        if st == "unjammed":
            return ("unjammed: no load-bearing particles remain, so the diameters are too small for this structure to "
                    "carry load. Quench the original state with a larger kappa_mult.")
        if st == "step_collapse":
            return (f"gave up: the minimiser's time step collapsed after {it} iterations (max residual {res:.2e}), so it "
                    f"found no step it could follow downhill. Quench the original (packed) state again with a lower dt_max "
                    f"(for example half the present one) or with slightly more compression (a larger kappa_mult); "
                    f"resuming with unchanged settings does not help.{size_note}")
        if st == "non_finite":
            return (f"failed: a non-finite force or energy appeared after {it} iterations; the stored state may contain "
                    "invalid values and should not be used.")
        return "no quench has been run."

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
        data = self._to_dict_base()
        if self._nonsphere:   # ellipsoids and rods only: a sphere packing gets none of these keys (its only new key is pack_config.quench_rod_round)
            data["final_result"].update(phi_contact=self.phi_contact, shrink_factor=self.shrink_factor,
                                        torque_magnitude=self.torque_magnitude)
            data["trajectory"]["orientations"] = deepcopy(self.trajectory_orientations)
        if self.trajectory_origin:   # only with a captured trajectory: a packing without one serialises exactly as before
            data["trajectory"]["origins"] = deepcopy(self.trajectory_origin)
        return data

    def _to_dict_base(self) -> dict[str, Any]:
        return {
            "init_config": self._collect_init_config(),
            "pack_config": self._collect_pack_config(translate=False),
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
