#include "rcpgenerator/initialize_particles.hpp"
#include "rcpgenerator/rcp_generator.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

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
    return config;
}

rcpgenerator::PackingRunOptions packing_run_options_from_dict(const py::dict& source) {
    rcpgenerator::PackingRunOptions options;
    assign_if_present(source, "max_steps", options.max_steps);
    assign_if_present(source, "fix_diameter", options.fix_diameter);
    if (source.contains(py::str("mu"))) {
        options.has_mu_override = true;
        options.mu_override = source[py::str("mu")].cast<double>();
    }
    if (source.contains(py::str("target_phi"))) {
        options.has_target_phi = true;
        options.target_phi = source[py::str("target_phi")].cast<double>();
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
    output["steps"] = result.steps;
    output["phi"] = result.phi;
    output["max_min_dist"] = result.max_min_dist;
    output["force_magnitude"] = result.force_magnitude;
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
