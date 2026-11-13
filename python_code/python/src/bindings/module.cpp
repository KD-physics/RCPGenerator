#include "rcpgenerator/initialize_particles.hpp"
#include "rcpgenerator/rcp_generator.hpp"
#include "rcpgenerator/ellipsoid.hpp"
#include "rcpgenerator/rod.hpp"
#include "rcpgenerator/curved_wall.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

template <typename T>
void assign_if_present(const py::dict& source, const char* key, T& target) {
    if (source.contains(py::str(key))) {
        target = source[py::str(key)].cast<T>();
    }
}

std::vector<std::int8_t> cast_walls(const py::handle& value) {
    std::vector<std::int8_t> walls;
    for (const auto item : value) {
        walls.push_back(py::cast<std::int8_t>(item));
    }
    return walls;
}

py::list walls_to_python(const std::vector<std::int8_t>& walls) {
    py::list result;
    for (std::int8_t wall : walls) {
        result.append(py::int_(wall));
    }
    return result;
}

rcpgenerator::Distribution distribution_from_dict(const py::dict& source) {
    rcpgenerator::Distribution dist;
    assign_if_present(source, "type", dist.type);
    assign_if_present(source, "d", dist.d);
    assign_if_present(source, "mu", dist.mu);
    assign_if_present(source, "sigma", dist.sigma);
    assign_if_present(source, "mu1", dist.mu1);
    assign_if_present(source, "sigma1", dist.sigma1);
    assign_if_present(source, "mu2", dist.mu2);
    assign_if_present(source, "sigma2", dist.sigma2);
    assign_if_present(source, "p", dist.p);
    assign_if_present(source, "d1", dist.d1);
    assign_if_present(source, "d2", dist.d2);
    assign_if_present(source, "d_min", dist.d_min);
    assign_if_present(source, "d_max", dist.d_max);
    assign_if_present(source, "exponent", dist.exponent);
    assign_if_present(source, "scale", dist.scale);
    assign_if_present(source, "shape", dist.shape);
    assign_if_present(source, "custom", dist.custom);
    return dist;
}

rcpgenerator::InitializerConfig initializer_config_from_dict(const py::dict& source) {
    rcpgenerator::InitializerConfig config;
    assign_if_present(source, "phi", config.phi);
    assign_if_present(source, "N", config.N);
    assign_if_present(source, "Ndim", config.Ndim);
    assign_if_present(source, "box", config.box);
    assign_if_present(source, "fix_height", config.fix_height);
    assign_if_present(source, "seed", config.seed);
    if (source.contains(py::str("walls"))) {
        config.walls = cast_walls(source[py::str("walls")]);
    }
    if (source.contains(py::str("dist"))) {
        config.dist = distribution_from_dict(source[py::str("dist")].cast<py::dict>());
    }
    return config;
}

rcpgenerator::PackingInput packing_input_from_dict(const py::dict& source) {
    rcpgenerator::PackingInput input;
    assign_if_present(source, "positions", input.positions);
    assign_if_present(source, "diameters", input.diameters);
    assign_if_present(source, "fixed", input.fixed);   // cascade: frozen-diameter mask
    assign_if_present(source, "orientations", input.orientations);
    assign_if_present(source, "aspect_ratios", input.aspect_ratios);
    assign_if_present(source, "rod_length", input.rod_length);   // rods: 1 or N numbers each
    assign_if_present(source, "rod_radius", input.rod_radius);   // ellipsoids: per-particle table, N x (Ndim-1)   // ellipsoids: N x (Ndim*Ndim), row-major
    return input;
}

rcpgenerator::PackingConfig packing_config_from_dict(const py::dict& source) {
    rcpgenerator::PackingConfig config;
    assign_if_present(source, "box", config.box);
    assign_if_present(source, "neighbor_max", config.neighbor_max);
    assign_if_present(source, "seed", config.seed);
    assign_if_present(source, "fix_height", config.fix_height);
    if (source.contains(py::str("walls"))) {
        config.walls = cast_walls(source[py::str("walls")]);
    }
    // Diameter-SWAP knobs (all inert when swap == 0).
    assign_if_present(source, "swap", config.swap);
    assign_if_present(source, "swap_interval", config.swap_interval);
    assign_if_present(source, "swap_mode", config.swap_mode);
    assign_if_present(source, "swap_rank_window", config.swap_rank_window);
    assign_if_present(source, "swap_size_tol", config.swap_size_tol);
    assign_if_present(source, "swap_eps", config.swap_eps);
    assign_if_present(source, "swap_start", config.swap_start);
    assign_if_present(source, "swap_stop", config.swap_stop);
    assign_if_present(source, "swap_phase", config.swap_phase);
    assign_if_present(source, "swap_adam", config.swap_adam);
    assign_if_present(source, "swap_inject_adam", config.swap_inject_adam);
    assign_if_present(source, "swap_inject_freeze", config.swap_inject_freeze);
    assign_if_present(source, "swap_inject_accept", config.swap_inject_accept);
    assign_if_present(source, "swap_inject_tol", config.swap_inject_tol);
    assign_if_present(source, "swap_inject_freeze_shell", config.swap_inject_freeze_shell);
    assign_if_present(source, "swap_max_per_particle", config.swap_max_per_particle);
    assign_if_present(source, "swap_refresh", config.swap_refresh);
    assign_if_present(source, "swap_backoff", config.swap_backoff);
    assign_if_present(source, "swap_exact", config.swap_exact);
    assign_if_present(source, "swap_window_rel", config.swap_window_rel);
    assign_if_present(source, "swap_relayout_ratio", config.swap_relayout_ratio);
    assign_if_present(source, "swap_rank_space", config.swap_rank_space);
    assign_if_present(source, "swap_idle_frac", config.swap_idle_frac);
    assign_if_present(source, "swap_ctrl", config.swap_ctrl);
    assign_if_present(source, "swap_a", config.swap_a);
    assign_if_present(source, "swap_R", config.swap_R);
    assign_if_present(source, "swap_eta", config.swap_eta);
    assign_if_present(source, "swap_wmin", config.swap_wmin);
    assign_if_present(source, "swap_wmax", config.swap_wmax);
    assign_if_present(source, "swap_smax", config.swap_smax);
    assign_if_present(source, "swap_floor_budget", config.swap_floor_budget);
    assign_if_present(source, "swap_rescue", config.swap_rescue);
    assign_if_present(source, "swap_rescue_w", config.swap_rescue_w);
    assign_if_present(source, "swap_rescue_interval", config.swap_rescue_interval);
    assign_if_present(source, "swap_cycle3", config.swap_cycle3);
    assign_if_present(source, "checkpoint_path", config.checkpoint_path);
    assign_if_present(source, "checkpoint_interval_s", config.checkpoint_interval_s);
    assign_if_present(source, "max_steps", config.max_steps);
    assign_if_present(source, "swap_event_log_stride", config.swap_event_log_stride);
    assign_if_present(source, "swap_event_log_cap", config.swap_event_log_cap);
    assign_if_present(source, "swap_event_log_geom", config.swap_event_log_geom);
    assign_if_present(source, "checkpoint_every_steps", config.checkpoint_every_steps);
    assign_if_present(source, "mu_hold_start", config.mu_hold_start);
    assign_if_present(source, "mu_hold_steps", config.mu_hold_steps);
    assign_if_present(source, "mu_hold_mult", config.mu_hold_mult);
    assign_if_present(source, "mu_hold_shift_deadline", config.mu_hold_shift_deadline);
    assign_if_present(source, "record_visibility", config.record_visibility);
    assign_if_present(source, "swap_mu_ctrl", config.swap_mu_ctrl);
    assign_if_present(source, "mu_ctrl_v_target", config.mu_ctrl_v_target);
    assign_if_present(source, "mu_ctrl_deadband", config.mu_ctrl_deadband);
    assign_if_present(source, "mu_ctrl_eta", config.mu_ctrl_eta);
    assign_if_present(source, "mu_ctrl_interval", config.mu_ctrl_interval);
    assign_if_present(source, "mu_ctrl_min_mult", config.mu_ctrl_min_mult);
    assign_if_present(source, "mu_ctrl_max_mult", config.mu_ctrl_max_mult);
    assign_if_present(source, "swap_geo", config.swap_geo);
    assign_if_present(source, "swap_geo_T", config.swap_geo_T);
    assign_if_present(source, "swap_geo_T_end", config.swap_geo_T_end);
    assign_if_present(source, "swap_geo_T_decay", config.swap_geo_T_decay);
    assign_if_present(source, "swap_geo_w", config.swap_geo_w);
    assign_if_present(source, "swap_geo_B", config.swap_geo_B);
    assign_if_present(source, "swap_geo_k", config.swap_geo_k);
    assign_if_present(source, "swap_geo_start", config.swap_geo_start);
    assign_if_present(source, "swap_geo_stop", config.swap_geo_stop);
    assign_if_present(source, "swap_geo_reopen", config.swap_geo_reopen);
    assign_if_present(source, "swap_geo_cycle", config.swap_geo_cycle);
    assign_if_present(source, "swap_geo_cycle_cap", config.swap_geo_cycle_cap);
    assign_if_present(source, "swap_geo_cycle4", config.swap_geo_cycle4);
    assign_if_present(source, "swap_geo_cycle4_cap", config.swap_geo_cycle4_cap);
    assign_if_present(source, "swap_geo_capk", config.swap_geo_capk);
    assign_if_present(source, "swap_geo_capmode", config.swap_geo_capmode);
    assign_if_present(source, "swap_geo_capk_start", config.swap_geo_capk_start);
    assign_if_present(source, "swap_geo_bottleneck", config.swap_geo_bottleneck);
    assign_if_present(source, "swap_geo_cage", config.swap_geo_cage);
    assign_if_present(source, "swap_geo_cage_move", config.swap_geo_cage_move);
    assign_if_present(source, "swap_geo_nonsep", config.swap_geo_nonsep);
    assign_if_present(source, "swap_geo_ctarget", config.swap_geo_ctarget);
    assign_if_present(source, "swap_geo_nsk", config.swap_geo_nsk);
    assign_if_present(source, "swap_geo_nonsep_start", config.swap_geo_nonsep_start);
    assign_if_present(source, "swap_geo_nonsep_stop", config.swap_geo_nonsep_stop);
    assign_if_present(source, "swap_geo_cage_start", config.swap_geo_cage_start);
    assign_if_present(source, "swap_sweeps_off", config.swap_sweeps_off);
    assign_if_present(source, "swap_inject", config.swap_inject);
    // [219] distribution morph + mu pin (default off)
    assign_if_present(source, "morph_target", config.morph_target);
    assign_if_present(source, "morph_start", config.morph_start);
    assign_if_present(source, "morph_interval", config.morph_interval);
    assign_if_present(source, "morph_rate0", config.morph_rate0);
    assign_if_present(source, "morph_band", config.morph_band);
    assign_if_present(source, "morph_grow", config.morph_grow);
    assign_if_present(source, "morph_shrink", config.morph_shrink);
    assign_if_present(source, "morph_grow_after", config.morph_grow_after);
    assign_if_present(source, "morph_rate_max", config.morph_rate_max);
    assign_if_present(source, "morph_settle_flat", config.morph_settle_flat);
    assign_if_present(source, "morph_settle_window", config.morph_settle_window);
    assign_if_present(source, "morph_kappa_damp", config.morph_kappa_damp);
    assign_if_present(source, "morph_band_floor", config.morph_band_floor);
    assign_if_present(source, "morph_margin", config.morph_margin);
    assign_if_present(source, "morph_max_hold_streak", config.morph_max_hold_streak);
    assign_if_present(source, "morph_kappa_drop", config.morph_kappa_drop);
    assign_if_present(source, "morph_kappa_tau", config.morph_kappa_tau);
    assign_if_present(source, "morph_rate_min", config.morph_rate_min);
    assign_if_present(source, "mu_pin", config.mu_pin);
    assign_if_present(source, "quench_mode", config.quench_mode);
    assign_if_present(source, "quench_tol", config.quench_tol);
    assign_if_present(source, "quench_max_iter", config.quench_max_iter);
    assign_if_present(source, "quench_dt0", config.quench_dt0);
    assign_if_present(source, "quench_dt_max", config.quench_dt_max);
    assign_if_present(source, "quench_alpha0", config.quench_alpha0);
    assign_if_present(source, "quench_mass_scaled", config.quench_mass_scaled);
    assign_if_present(source, "quench_report_every", config.quench_report_every);
    assign_if_present(source, "quench_kappa_mult", config.quench_kappa_mult);
    assign_if_present(source, "quench_mix_mode", config.quench_mix_mode);
    assign_if_present(source, "quench_two_level", config.quench_two_level);
    assign_if_present(source, "quench_light_frac", config.quench_light_frac);
    assign_if_present(source, "quench_light_omega", config.quench_light_omega);
    assign_if_present(source, "quench_warm_steps", config.quench_warm_steps);
    assign_if_present(source, "quench_freeze_damp", config.quench_freeze_damp);
    assign_if_present(source, "quench_freeze_frac", config.quench_freeze_frac);
    assign_if_present(source, "quench_warm_amp", config.quench_warm_amp);
    assign_if_present(source, "quench_class_every", config.quench_class_every);
    assign_if_present(source, "quench_held", config.quench_held);
    assign_if_present(source, "quench_stall_window", config.quench_stall_window);
    assign_if_present(source, "quench_stall_factor", config.quench_stall_factor);
    assign_if_present(source, "quench_dt_min", config.quench_dt_min);
    // Ellipsoids (absent -> spheres).
    assign_if_present(source, "aspect_ratios", config.aspect_ratios);
    assign_if_present(source, "ell_iters", config.ell_iters);
    assign_if_present(source, "ell_rot_scale", config.ell_rot_scale);
    assign_if_present(source, "ell_track_tol", config.ell_track_tol);
    assign_if_present(source, "ell_permute_orient", config.ell_permute_orient);
    assign_if_present(source, "rod_capped", config.rod_capped);
    assign_if_present(source, "rod_round", config.rod_round);
    assign_if_present(source, "quench_rod_round", config.quench_rod_round);
    assign_if_present(source, "trajectory_per_decade", config.trajectory_per_decade);
    assign_if_present(source, "ell_permute_partners", config.ell_permute_partners);
    assign_if_present(source, "ell_permute_guard", config.ell_permute_guard);
    assign_if_present(source, "ell_capture", config.ell_capture);
    return config;
}

rcpgenerator::PackingRunOptions packing_run_options_from_dict(const py::dict& source) {
    rcpgenerator::PackingRunOptions options;
    assign_if_present(source, "max_steps", options.max_steps);
    assign_if_present(source, "fix_diameter", options.fix_diameter);
    assign_if_present(source, "reorder_interval", options.reorder_interval);
    if (source.contains(py::str("mu"))) {
        options.has_mu_override = true;
        options.mu_override = source[py::str("mu")].cast<double>();
    }
    if (source.contains(py::str("target_phi"))) {
        options.has_target_phi = true;
        options.target_phi = source[py::str("target_phi")].cast<double>();
    }
    if (source.contains(py::str("initial_phi_target"))) {
        options.has_initial_phi_target = true;
        options.initial_phi_target =
            source[py::str("initial_phi_target")].cast<double>();
    }
    return options;
}

py::dict distribution_to_dict(const rcpgenerator::Distribution& dist) {
    py::dict result;
    result["type"] = dist.type;
    result["d"] = dist.d;
    result["mu"] = dist.mu;
    result["sigma"] = dist.sigma;
    result["mu1"] = dist.mu1;
    result["sigma1"] = dist.sigma1;
    result["mu2"] = dist.mu2;
    result["sigma2"] = dist.sigma2;
    result["p"] = dist.p;
    result["d1"] = dist.d1;
    result["d2"] = dist.d2;
    result["d_min"] = dist.d_min;
    result["d_max"] = dist.d_max;
    result["exponent"] = dist.exponent;
    result["scale"] = dist.scale;
    result["shape"] = dist.shape;
    result["custom"] = dist.custom;
    return result;
}

py::dict initializer_result_to_dict(const rcpgenerator::InitializerResult& result) {
    py::dict output;
    output["positions"] = result.positions;
    output["diameters"] = result.diameters;
    output["box"] = result.box;
    output["walls"] = walls_to_python(result.walls);
    output["diameter_scale_factor"] = result.diameter_scale_factor;
    output["phi_modifier"] = result.phi_modifier;
    return output;
}

py::dict packing_result_to_dict(const rcpgenerator::PackingResult& result) {
    py::dict output;
    output["positions"] = result.positions;
    output["diameters"] = result.diameters;
    output["box"] = result.box;
    output["walls"] = walls_to_python(result.walls);
    output["phi_history"] = result.phi_history;
    output["force_history"] = result.force_history;
    output["energy_history"] = result.energy_history;
    // Cycle 16: per-step diagnostics for ADAM convergence study.
    output["max_overlap_history"] = result.max_overlap_history;
    output["mu_history"] = result.mu_history;
    output["alpha_history"] = result.alpha_history;
    output["mu_flag_history"] = result.mu_flag_history;
    output["zpos_history"] = result.zpos_history;
    output["mu_eff_history"] = result.mu_eff_history;
    output["mu_hold_applied"] = result.mu_hold_applied;
    output["mu_ctrl_final_mult"] = result.mu_ctrl_final_mult;
    output["swap_geo_trades"] = result.swap_geo_trades;
    output["swap_geo_passes"] = result.swap_geo_passes;
    output["swap_geo_reopens"] = result.swap_geo_reopens;
    output["swap_geo_cycle_trades"] = result.swap_geo_cycle_trades;
    output["swap_geo_cycle4_trades"] = result.swap_geo_cycle4_trades;
    output["swap_geo_capk_trades"] = result.swap_geo_capk_trades;
    output["swap_geo_cage_trades"] = result.swap_geo_cage_trades;
    output["swap_geo_bottleneck_trades"] = result.swap_geo_bottleneck_trades;
    // Cycle 18: per-particle troublemaker stats from last window.
    output["pp_F_sq_sum"] = result.pp_F_sq_sum;
    output["pp_F_sign_flips"] = result.pp_F_sign_flips;
    output["pp_pos_range_max"] = result.pp_pos_range_max;
    output["pp_contact_sum"] = result.pp_contact_sum;
    output["pp_signflip_resets"] = result.pp_signflip_resets;
    output["pp_window_steps"] = result.pp_window_steps;
    // Cycle 18 Phase A: spike event log (flat, 8 doubles per event).
    output["spike_log"] = result.spike_log;
    output["reorder_log"] = result.reorder_log;
    output["swap_relayouts"] = result.swap_relayouts;
    output["swap_rescue_candidates"] = result.swap_rescue_candidates;
    output["swap_rescue_accepted"] = result.swap_rescue_accepted;
    output["swap_rescue_false_rej"] = result.swap_rescue_false_rej;
    output["swap_rescue_false_acc"] = result.swap_rescue_false_acc;
    output["swap_cycle3_candidates"] = result.swap_cycle3_candidates;
    output["swap_cycle3_accepted"] = result.swap_cycle3_accepted;
    output["swap_cycle3_blocked"] = result.swap_cycle3_blocked;
    output["capped"] = result.capped;
    output["resumed_from_step"] = result.resumed_from_step;
    output["checkpoint_last_step"] = result.checkpoint_last_step;
    output["swap_events"] = result.swap_events;
    output["swap_event_width"] = result.swap_event_width;
    output["steps"] = result.steps;
    output["phi"] = result.phi;
    output["max_min_dist"] = result.max_min_dist;
    output["force_magnitude"] = result.force_magnitude;
    output["particle_origin"] = result.particle_origin;
    // [219] distribution morph diagnostics
    output["morph_log"] = result.morph_log;
    output["morph_lambda"] = result.morph_lambda;
    output["morph_updates"] = result.morph_updates;
    output["morph_holds"] = result.morph_holds;
    output["morph_stalled"] = result.morph_stalled;
    output["morph_rate_ceiling"] = result.morph_rate_ceiling;
    output["quench_residual"] = result.quench_residual;
    output["quench_global_residual"] = result.quench_global_residual;
    output["quench_iters"] = result.quench_iters;
    output["quench_converged"] = result.quench_converged;
    output["quench_status"] = result.quench_status;
    output["quench_log"] = result.quench_log;
    output["quench_held_force"] = result.quench_held_force;
    // Ellipsoids only: a sphere run's result dict has exactly the keys it always had.
    if (!result.aspect_ratios.empty()) output["aspect_ratios"] = result.aspect_ratios;
    if (!result.rod_length.empty()) { output["rod_length"] = result.rod_length; output["rod_radius"] = result.rod_radius; }
    if (!result.orientations.empty()) output["ell_permute_candidates"] = result.ell_permute_candidates;
    if (result.ell_capture_step >= 0) {
        py::dict cap;
        cap["step"] = result.ell_capture_step;
        cap["energy"] = result.ell_capture_energy;
        cap["positions"] = result.ell_capture_positions;       // flat N x Ndim
        cap["diameters"] = result.ell_capture_diameters;
        cap["orientations"] = result.ell_capture_orientations; // flat N x Ndim^2
        cap["aspect_ratios"] = result.ell_capture_ratios;      // flat N x (Ndim-1)
        cap["forces"] = result.ell_capture_forces;             // flat N x Ndim (the gradient, as in evaluate_once)
        cap["torques"] = result.ell_capture_torques;           // flat N x nrot
        cap["exchanges"] = result.ell_capture_exchanges;       // (i, j, gain) of the pass just before
        cap["z"] = result.ell_capture_z;                       // pair contacts per particle as the engine counted them
        output["ell_capture"] = cap;
    }
    if (!result.orientations.empty()) {
        output["orientations"] = result.orientations;
        output["torque_magnitude"] = result.torque_magnitude;
    }
    return output;
}

py::dict packing_trace_to_dict(const rcpgenerator::PackingTrace& trace) {
    py::dict output;
    output["steps"] = trace.steps;
    output["positions"] = trace.positions;
    output["diameters"] = trace.diameters;
    output["phi"] = trace.phi;
    output["force"] = trace.force;
    output["energy"] = trace.energy;
    output["max_min_dist"] = trace.max_min_dist;
    output["origins"] = trace.origins;
    if (!trace.orientations.empty()) {
        output["orientations"] = trace.orientations;
    }
    return output;
}

}  // namespace

PYBIND11_MODULE(_rcpgenerator, m) {
    m.doc() = "Thin pybind11 bindings for rcpgenerator";

    // Default to single-threaded so behavior matches the legacy serial build
    // unless the user explicitly opts into threading via set_num_threads().
    // OMP_NUM_THREADS environment variable, if set, overrides this default.
    if (std::getenv("OMP_NUM_THREADS") == nullptr) {
        omp_set_num_threads(1);
    }

    m.def(
        "set_num_threads",
        [](int n) {
            if (n < 1) {
                throw std::invalid_argument("num_threads must be >= 1");
            }
            omp_set_num_threads(n);
        },
        py::arg("n"),
        "Set the number of OpenMP threads used by the parallel kernels. "
        "Default is 1 (single-threaded — bit-equivalent to the legacy build). "
        "Pass n > 1 to enable parallel force loop, neighbor build, etc.");

    m.def(
        "get_num_threads",
        []() { return omp_get_max_threads(); },
        "Return the current OpenMP thread count.");

    m.def(
        "initialize_particles",
        [](const py::dict& config_dict) {
            return initializer_result_to_dict(
                rcpgenerator::initialize_particles(initializer_config_from_dict(config_dict)));
        },
        py::arg("config"));

    m.def(
        "run_packing",
        [](const py::dict& input_dict, const py::dict& config_dict) {
            return packing_result_to_dict(
                rcpgenerator::run_packing(
                    packing_input_from_dict(input_dict),
                    packing_config_from_dict(config_dict)));
        },
        py::arg("input"),
        py::arg("config"));

    m.def(
        "run_packing_observed",
        [](const py::dict& input_dict,
           const py::dict& config_dict,
           std::size_t progress_interval,
           bool capture_positions,
           std::size_t trajectory_interval,
           py::object callback,
           py::dict options_dict) {
            py::function py_callback;
            if (!callback.is_none()) {
                py_callback = callback.cast<py::function>();
            }

            auto observed = rcpgenerator::run_packing_observed(
                packing_input_from_dict(input_dict),
                packing_config_from_dict(config_dict),
                progress_interval,
                capture_positions,
                trajectory_interval,
                py_callback
                    ? rcpgenerator::PackingObserver(
                          [&](const rcpgenerator::PackingObservation& observation,
                              const std::vector<std::vector<double>>* positions) {
                              py::dict update;
                              update["step"] = observation.step;
                              update["phi"] = observation.phi;
                              update["force_magnitude"] = observation.force_magnitude;
                              update["energy"] = observation.energy;
                              update["max_min_dist"] = observation.max_min_dist;
                              if (positions != nullptr) {
                                  update["positions"] = *positions;
                              }
                              py_callback(update);
                          })
                    : rcpgenerator::PackingObserver(),
                packing_run_options_from_dict(options_dict));

            py::dict output;
            output["result"] = packing_result_to_dict(observed.first);
            output["trace"] = packing_trace_to_dict(observed.second);
            return output;
        },
        py::arg("input"),
        py::arg("config"),
        py::arg("progress_interval") = 0,
        py::arg("capture_positions") = false,
        py::arg("trajectory_interval") = 0,
        py::arg("callback") = py::none(),
        py::arg("options") = py::dict());

    // ---- ellipsoids and single evaluations (additive) -------------------------------------------------------
    m.def(
        "evaluate_once",
        [](const py::dict& input_dict, const py::dict& config_dict, double mu) {
            const auto ev = rcpgenerator::evaluate_once(
                packing_input_from_dict(input_dict), packing_config_from_dict(config_dict), mu);
            py::dict out;
            out["forces"] = ev.forces;
            out["torques"] = ev.torques;
            out["U"] = ev.U;
            out["dkappa"] = ev.dkappa;
            out["max_overlap"] = ev.max_overlap;
            out["energy"] = ev.energy;
            out["Fmean"] = ev.Fmean;
            out["Lc"] = ev.Lc;
            out["z"] = ev.z;
            return out;
        },
        py::arg("input"), py::arg("config"), py::arg("mu") = 0.0,
        "One evaluation of forces, torques, energy and the kappa gradient on a configuration as given. "
        "No aspect_ratios and no orientations -> the sphere kernel; both -> the ellipsoid kernel.");

    m.def(
        "contact_scale",
        [](const py::dict& input_dict, const py::dict& config_dict) {
            const auto cs = rcpgenerator::contact_scale(
                packing_input_from_dict(input_dict), packing_config_from_dict(config_dict));
            py::dict out;
            out["scale"] = cs.scale;
            out["pair_i"] = cs.pair_i;
            out["pair_j"] = cs.pair_j;
            out["n_overlapping"] = cs.n_overlapping;
            out["particle_volume"] = cs.particle_volume;
            return out;
        },
        py::arg("input"), py::arg("config"),
        "Exact shrink-to-contact factor of a configuration (converged contact function, all pairs).");

    auto permute_out = [](const rcpgenerator::PermuteOnce& r) {
        py::dict out;
        out["diameters"] = r.diameters;
        out["aspect_ratios"] = r.aspect_ratios;
        out["rod_length"] = r.rod_length;
        out["rod_radius"] = r.rod_radius;
        out["orientations"] = r.orientations;
        out["exchanges"] = r.exchanges;
        out["candidates"] = r.candidates;
        return out;
    };
    m.def(
        "permute_once",
        [permute_out](const py::dict& input_dict, const py::dict& config_dict) {
            return permute_out(rcpgenerator::permute_once(packing_input_from_dict(input_dict), packing_config_from_dict(config_dict)));
        },
        py::arg("input"), py::arg("config"),
        "One ellipsoid permutation pass on a configuration as given (all pairs listed). Returns the new diameters, "
        "per-particle aspect ratios and orientations, and the accepted exchanges as flat (i, j, gain) triples.");
    m.def(
        "exchange_once",
        [permute_out](const py::dict& input_dict, const py::dict& config_dict, std::size_t i, std::size_t j) {
            return permute_out(rcpgenerator::exchange_once(packing_input_from_dict(input_dict), packing_config_from_dict(config_dict), i, j));
        },
        py::arg("input"), py::arg("config"), py::arg("i"), py::arg("j"),
        "The exchange move itself, forced on particles i and j: diameters, shapes and (ell_permute_orient) orientations trade places.");
    m.def(
        "exchange_admissible",
        [](double Di, const std::vector<double>& ri, double Dj, const std::vector<double>& rj, double w) {
            const std::size_t n = ri.size() + 1;
            std::vector<double> qi(n, 1.0), qj(n, 1.0);
            for (std::size_t k = 1; k < n; ++k) { qi[k] = ri[k - 1]; qj[k] = rj.at(k - 1); }
            return rcpgenerator::ell_exchange_admissible(n, Di, qi.data(), Dj, qj.data(), w);
        },
        py::arg("Di"), py::arg("ratios_i"), py::arg("Dj"), py::arg("ratios_j"), py::arg("window"),
        "May two particles be exchanged: every axis length within the relative window.");

    m.def("rod_kernel_stats", [](bool reset) {
        const auto st = rcpgenerator::rod_kernel_stats(reset);
        py::dict out;
        out["pairs"] = st.pairs; out["quick_apart"] = st.quick_apart; out["closed_form"] = st.closed_form;
        out["iterated"] = st.iterated; out["iterations"] = st.iterations; out["contacts"] = st.contacts; out["warm"] = st.warm;
        out["warm_failed"] = st.warm_failed; out["backtracks"] = st.backtracks; out["descents"] = st.descents;
        out["tracked_ok"] = st.tracked_ok; out["tracked_rejected"] = st.tracked_rejected; out["sep_cached_ok"] = st.sep_cached_ok;
        out["sep_cached_failed"] = st.sep_cached_failed; out["lean_ok"] = st.lean_ok; out["lean_failed"] = st.lean_failed;
        return out;
    }, py::arg("reset") = true, "Counters of the 3D uncapped rod kernel since the last reset (diagnostic).");
    m.def(
        "curved_wall_probe",
        [](std::size_t n, std::size_t t, const std::string& shape, const std::vector<double>& xc, const std::vector<double>& A,
           const std::vector<double>& u, double hl, double rad, double eps, double Rc, const std::vector<double>& Rot) {
            if (n < 2 || n > rcpgenerator::ELL_MAX_DIM || t < 2 || t > n || xc.size() != n) throw std::runtime_error("curved_wall_probe: bad n, t or xc");
            if (!Rot.empty() && Rot.size() != n * n) throw std::runtime_error("curved_wall_probe: Rot must have n*n entries");
            rcpgenerator::CurvedWallHit hit;
            double scale = 0.0;
            if (shape == "ellipsoid") {
                if (A.size() != n * n) throw std::runtime_error("curved_wall_probe: A must have n*n entries");
                rcpgenerator::ell_farthest_from_axis(n, t, A.data(), xc.data(), hit);
                scale = rcpgenerator::ell_wall_scale(n, t, A.data(), xc.data(), Rc, 10.0);
            } else {
                if (u.size() != n) throw std::runtime_error("curved_wall_probe: u must have n entries");
                const int kind = shape == "capped" ? rcpgenerator::SHAPE_ROD_CAPPED : rcpgenerator::SHAPE_ROD_UNCAPPED;
                if (kind == rcpgenerator::SHAPE_ROD_UNCAPPED && n > 3) throw std::runtime_error("curved_wall_probe: uncapped rods are 2D and 3D only");
                rcpgenerator::rod_farthest_from_axis(n, t, kind, u.data(), hl, rad, eps, xc.data(), hit);
                scale = rcpgenerator::rod_wall_scale(n, t, kind, u.data(), hl, rad, xc.data(), Rc, 10.0);
            }
            py::dict out;
            out["g"] = hit.g;
            out["nvec"] = std::vector<double>(hit.nvec, hit.nvec + n);
            out["p"] = std::vector<double>(hit.p, hit.p + n);
            out["scale"] = scale;
            rcpgenerator::CurvedWallContact c;
            const bool touching = rcpgenerator::curved_wall_contact(n, Rc, 1.0, hit, c);
            out["touching"] = touching;
            if (shape != "ellipsoid") {      // rods: the signed sum of point contacts the passes use
                std::vector<double> F(n, 0.0), T(n * (n - 1) / 2, 0.0);
                rcpgenerator::CurvedWallSum s;
                const int kind = shape == "capped" ? rcpgenerator::SHAPE_ROD_CAPPED : rcpgenerator::SHAPE_ROD_UNCAPPED;
                const bool any = rcpgenerator::rod_curved_wall(n, t, kind, u.data(), hl, rad, xc.data(), Rc, 1.0, F.data(), T.data(), s, Rot.empty() ? nullptr : Rot.data());
                out["contact"] = any;
                out["energy"] = s.energy; out["force"] = F; out["torque"] = T; out["dkappa"] = s.dkappa; out["U_add"] = s.U;
                out["F_mag"] = s.Fmean; out["count"] = s.count; out["f_sum"] = s.f;
                if (touching) { out["h"] = c.h; out["delta"] = c.delta; }
                return out;
            }
            if (touching) {
                std::vector<double> F(n, 0.0), T(n * (n - 1) / 2, 0.0);
                rcpgenerator::curved_wall_apply(n, hit, c, F.data(), T.data());
                out["h"] = c.h; out["delta"] = c.delta; out["energy"] = c.energy; out["force"] = F; out["torque"] = T;
                out["dkappa"] = c.dkappa; out["U_add"] = c.U_add; out["F_mag"] = c.F_mag;
            }
            return out;
        }, py::arg("n"), py::arg("t"), py::arg("shape"), py::arg("xc"), py::arg("A") = std::vector<double>(),
        py::arg("u") = std::vector<double>(), py::arg("hl") = 0.0, py::arg("rad") = 0.0, py::arg("eps") = 0.0, py::arg("Rc") = 0.5,
        py::arg("Rot") = std::vector<double>(),
        "One particle against a curved container wall (test hook): farthest point from the axis, contact quantities, scale to touch.");
    m.def(
        "rod_pair",
        [](const std::vector<double>& r, const std::vector<double>& ui, double hli, double radi,
           const std::vector<double>& uj, double hlj, double radj, bool capped, int max_iter) {
            const std::size_t n = r.size();
            rcpgenerator::RodPairOut o;
            std::vector<double> Ti(n * (n - 1) / 2, 0.0), Tj(n * (n - 1) / 2, 0.0);
            const int kind = capped ? rcpgenerator::SHAPE_ROD_CAPPED : rcpgenerator::SHAPE_ROD_UNCAPPED;
            const bool contact = rcpgenerator::rod_pair(n, kind, r.data(), ui.data(), hli, radi, uj.data(), hlj, radj, 1.0, max_iter, o, Ti.data(), Tj.data());
            py::dict out;
            out["contact"] = contact;
            out["scale_exact"] = rcpgenerator::rod_pair_scale_exact(n, kind, r.data(), ui.data(), hli, radi, uj.data(), hlj, radj);
            if (contact) {
                out["s"] = o.s; out["overlap"] = o.overlap; out["energy"] = o.energy; out["sigma"] = o.sigma;
                out["dkappa"] = o.dkappa; out["converged"] = o.converged;
                out["dEdr"] = std::vector<double>(o.dEdr, o.dEdr + n);
                out["normal"] = std::vector<double>(o.nvec, o.nvec + n);
                out["point_i"] = std::vector<double>(o.pi, o.pi + n);
                out["point_j"] = std::vector<double>(o.pj, o.pj + n);
                out["torque_i"] = Ti; out["torque_j"] = Tj;
            }
            return out;
        },
        py::arg("r"), py::arg("ui"), py::arg("hli"), py::arg("radi"), py::arg("uj"), py::arg("hlj"), py::arg("radj"),
        py::arg("capped"), py::arg("max_iter") = 400,
        "One pair of rods (actual sizes: half-lengths and radii; r = x_j - x_i): contact, scale, energy, dE/dr, "
        "torques, contact normal and the points where the pair force acts.");

    m.def(
        "random_orientations",
        [](std::size_t N, std::size_t Ndim, std::uint64_t seed) {
            std::vector<double> flat;
            rcpgenerator::ell_random_orientations(N, Ndim, seed, flat);
            std::vector<std::vector<double>> out(N, std::vector<double>(Ndim * Ndim));
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t k = 0; k < Ndim * Ndim; ++k) out[i][k] = flat[i * Ndim * Ndim + k];
            return out;
        },
        py::arg("N"), py::arg("Ndim"), py::arg("seed") = 0,
        "N uniformly random rotation matrices (row-major Ndim x Ndim each), deterministic in seed.");

    m.def(
        "ell_pair_contact",
        [](const std::vector<double>& Ai, const std::vector<double>& Aj, const std::vector<double>& r, int n_iter) {
            const std::size_t n = r.size();
            if (n < 2 || n > rcpgenerator::ELL_MAX_DIM || Ai.size() != n * n || Aj.size() != n * n)
                throw std::invalid_argument("ell_pair_contact: Ai, Aj must be flat n x n and r length n, 2 <= n <= 8");
            double s2 = 0.0, lam = 0.0, xv[rcpgenerator::ELL_MAX_DIM] = {0};
            const bool contact = rcpgenerator::ell_pair_contact(Ai.data(), Aj.data(), r.data(), n, n_iter, s2, lam, xv);
            py::dict out;
            out["contact"] = contact;          // false: certainly apart (a lower bound on s^2 reached 1)
            out["s2"] = s2;                    // meaningful only when the iteration ran to the end
            out["lam"] = lam;
            out["s2_exact"] = rcpgenerator::ell_pair_s2_exact(Ai.data(), Aj.data(), r.data(), n);
            return out;
        },
        py::arg("Ai"), py::arg("Aj"), py::arg("r"), py::arg("n_iter") = 5,
        "The pair contact function as the packing loop evaluates it (fixed n_iter), and its converged value.");

    m.def(
        "ell_pair_cost_probe",
        [](std::size_t Ndim, const std::vector<double>& aspect_ratios, int n_iter, std::size_t n_pairs,
           std::size_t repeats, std::uint64_t seed) {
            double te = 0.0, ts = 0.0, tt = 0.0, chk = 0.0;
            rcpgenerator::ell_pair_cost_probe(Ndim, aspect_ratios, n_iter, n_pairs, repeats, seed, te, ts, tt, chk);
            py::dict out;
            out["ellipsoid_s_per_pair"] = te;
            out["sphere_s_per_pair"] = ts;
            out["tracked_s_per_pair"] = tt;
            out["checksum"] = chk;
            return out;
        },
        py::arg("Ndim"), py::arg("aspect_ratios"), py::arg("n_iter") = 5, py::arg("n_pairs") = 4096,
        py::arg("repeats") = 200, py::arg("seed") = 1,
        "Wall time per overlapping pair of the ellipsoid pair kernel and of the sphere pair body.");
}
