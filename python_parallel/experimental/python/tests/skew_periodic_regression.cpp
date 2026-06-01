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

void require_close(double value, double expected, double tolerance, const char* label) {
    if (std::abs(value - expected) > tolerance) {
        throw std::runtime_error(std::string(label) + " is outside tolerance");
    }
}

}  // namespace

int main() {
    try {
        rcpgenerator::InitializerConfig finite_init;
        finite_init.phi = 0.11;
        finite_init.N = 32;
        finite_init.Ndim = 2;
        finite_init.box = {1.0, 1.0};
        finite_init.box_skew = {0.0, -20.0};
        finite_init.walls = {0, 0};
        finite_init.fix_height = false;
        finite_init.dist.type = "mono";
        finite_init.dist.d = 1.0;

        auto initialized = rcpgenerator::initialize_particles(finite_init);
        require_all_finite(initialized.positions, "finite-skew initialized position");
        require_all_finite(initialized.diameters, "finite-skew initialized diameter");
        require_close(initialized.box_skew[0], 0.0, 1e-12, "finite-skew corrected x skew");
        require_close(initialized.box_skew[1], -20.0, 1e-12, "finite-skew corrected y skew");

        rcpgenerator::PackingInput input;
        input.positions = initialized.positions;
        input.diameters = initialized.diameters;

        rcpgenerator::PackingConfig finite_pack;
        finite_pack.box = initialized.box;
        finite_pack.box_skew = initialized.box_skew;
        finite_pack.walls = initialized.walls;
        finite_pack.fix_height = false;
        finite_pack.neighbor_max = 0;
        finite_pack.seed = 123u;

        rcpgenerator::PackingRunOptions finite_options;
        finite_options.max_steps = 2000;

        auto finite_observed = rcpgenerator::run_packing_observed(
            input,
            finite_pack,
            0,
            false,
            0,
            rcpgenerator::PackingObserver(),
            finite_options);

        const auto& finite_result = finite_observed.first;
        require_all_finite(finite_result.positions, "finite-skew packed position");
        require_all_finite(finite_result.diameters, "finite-skew packed diameter");
        require_finite(finite_result.phi, "finite-skew packed phi");
        require_finite(finite_result.max_min_dist, "finite-skew packed max_min_dist");
        require_finite(finite_result.force_magnitude, "finite-skew packed force_magnitude");
        require_close(finite_result.box_skew[0], 0.0, 1e-12, "finite-skew packed x skew");
        require_close(finite_result.box_skew[1], -20.0, 1e-12, "finite-skew packed y skew");

        rcpgenerator::InitializerConfig guardrail_init = finite_init;
        guardrail_init.box_skew = {0.0, -75.0};

        auto guarded = rcpgenerator::initialize_particles(guardrail_init);
        require_all_finite(guarded.positions, "guardrail initialized position");
        require_all_finite(guarded.diameters, "guardrail initialized diameter");
        require_close(guarded.box_skew[0], 0.0, 1e-12, "guardrail corrected x skew");
        require_close(guarded.box_skew[1], -60.0, 0.2, "guardrail corrected y skew");

        std::cout << "Skew periodic regression checks passed.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Skew periodic regression failure: " << ex.what() << '\n';
        return 1;
    }
}
