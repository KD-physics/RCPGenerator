#include "rcpgenerator/initialize_particles.hpp"
#include "rcpgenerator/rcp_generator.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require_finite(double value, const char* label) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string(label) + " is not finite");
    }
}

template <typename Matrix>
void require_all_finite(const Matrix& values, const char* label) {
    for (const auto& row : values) {
        for (double value : row) {
            require_finite(value, label);
        }
    }
}

void require_all_finite(const std::vector<double>& values, const char* label) {
    for (double value : values) {
        require_finite(value, label);
    }
}

}  // namespace

int main() {
    try {
        rcpgenerator::InitializerConfig init_config;
        init_config.phi = 0.11;
        init_config.N = 64;
        init_config.Ndim = 2;
        init_config.box = {1.0, 1.0};
        init_config.walls = {1, 1};
        init_config.fix_height = false;
        init_config.dist.type = "mono";
        init_config.dist.d = 1.0;

        auto initialized = rcpgenerator::initialize_particles(init_config);
        require_all_finite(initialized.positions, "initialized position");
        require_all_finite(initialized.diameters, "initialized diameter");
        require_all_finite(initialized.box, "initialized box");
        require_finite(initialized.diameter_scale_factor, "diameter scale factor");
        require_finite(initialized.phi_modifier, "phi modifier");

        rcpgenerator::PackingInput input;
        input.positions = initialized.positions;
        input.diameters = initialized.diameters;

        rcpgenerator::PackingConfig pack_config;
        pack_config.box = initialized.box;
        pack_config.walls = initialized.walls;
        pack_config.fix_height = false;
        pack_config.neighbor_max = 0;
        pack_config.seed = 123u;

        rcpgenerator::PackingRunOptions options;
        options.max_steps = 500;

        auto observed = rcpgenerator::run_packing_observed(
            input,
            pack_config,
            0,
            false,
            0,
            rcpgenerator::PackingObserver(),
            options);

        const auto& result = observed.first;
        require_all_finite(result.positions, "packed position");
        require_all_finite(result.diameters, "packed diameter");
        require_all_finite(result.box, "packed box");
        require_finite(result.phi, "packed phi");
        require_finite(result.max_min_dist, "packed max_min_dist");
        require_finite(result.force_magnitude, "packed force_magnitude");

        std::cout << "Hard-wall regression checks passed.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Hard-wall regression failure: " << ex.what() << '\n';
        return 1;
    }
}
