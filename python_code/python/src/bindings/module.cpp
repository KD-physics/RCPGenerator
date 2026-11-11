#include "rcpgenerator/initialize_particles.hpp"
#include "rcpgenerator/rcp_generator.hpp"

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
}
