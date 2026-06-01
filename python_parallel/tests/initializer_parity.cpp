#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "rcpgenerator/initialize_particles.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define main legacy_initialize_particles_main
#include "../../InitializeParticles.cpp"
#undef main

namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

bool doubles_equal(double expected, double actual) {
    return (std::isnan(expected) && std::isnan(actual)) || expected == actual;
}

template <typename T>
void compare_scalar(const std::string& label, const T& expected, const T& actual) {
    if (!(expected == actual)) {
        std::ostringstream oss;
        oss << label << " mismatch: expected=" << expected << " actual=" << actual;
        fail(oss.str());
    }
}

void compare_double(const std::string& label, double expected, double actual) {
    if (!doubles_equal(expected, actual)) {
        std::ostringstream oss;
        oss << std::setprecision(17)
            << label << " mismatch: expected=" << expected << " actual=" << actual;
        fail(oss.str());
    }
}

void compare_vector(const std::string& label,
                    const std::vector<double>& expected,
                    const std::vector<double>& actual) {
    if (expected.size() != actual.size()) {
        std::ostringstream oss;
        oss << label << " size mismatch: expected=" << expected.size()
            << " actual=" << actual.size();
        fail(oss.str());
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!doubles_equal(expected[i], actual[i])) {
            std::ostringstream oss;
            oss << std::setprecision(17)
                << label << "[" << i << "] mismatch: expected=" << expected[i]
                << " actual=" << actual[i];
            fail(oss.str());
        }
    }
}

void compare_vector_i8(const std::string& label,
                       const std::vector<std::int8_t>& expected,
                       const std::vector<std::int8_t>& actual) {
    if (expected.size() != actual.size()) {
        std::ostringstream oss;
        oss << label << " size mismatch: expected=" << expected.size()
            << " actual=" << actual.size();
        fail(oss.str());
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (expected[i] != actual[i]) {
            std::ostringstream oss;
            oss << label << "[" << i << "] mismatch: expected="
                << static_cast<int>(expected[i]) << " actual="
                << static_cast<int>(actual[i]);
            fail(oss.str());
        }
    }
}

Distribution to_legacy_distribution(const rcpgenerator::Distribution& dist) {
    Distribution legacy;
    legacy.type = dist.type;
    legacy.d = dist.d;
    legacy.mu = dist.mu;
    legacy.sigma = dist.sigma;
    legacy.mu1 = dist.mu1;
    legacy.sigma1 = dist.sigma1;
    legacy.mu2 = dist.mu2;
    legacy.sigma2 = dist.sigma2;
    legacy.p = dist.p;
    legacy.d1 = dist.d1;
    legacy.d2 = dist.d2;
    legacy.d_min = dist.d_min;
    legacy.d_max = dist.d_max;
    legacy.exponent = dist.exponent;
    legacy.scale = dist.scale;
    legacy.shape = dist.shape;
    legacy.custom = dist.custom;
    return legacy;
}

struct LegacyInitializerDeterministicResult {
    std::vector<double> diameters;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    double diameter_scale_factor = 0.0;
    double phi_modifier = 0.0;
};

LegacyInitializerDeterministicResult run_legacy_initializer_deterministic(
    const rcpgenerator::InitializerConfig& config) {
    double phi = config.phi;
    std::size_t N = config.N;
    std::size_t Ndim = config.Ndim;
    std::vector<double> box = config.box;
    std::vector<std::int8_t> walls = config.walls;
    bool fix_height = config.fix_height;
    Distribution dist = to_legacy_distribution(config.dist);

    if (Ndim < 2) {
        fail("legacy deterministic initializer setup requires Ndim >= 2");
    }
    if (box.empty()) box.assign(Ndim, 1.0);
    if (phi <= 0 || phi >= 1 || N == 0) {
        fail("legacy deterministic initializer setup requires valid phi/N");
    }

    double Llast = box[0];
    double h = 1;
    if (fix_height)
    {
        h = box[Ndim - 1];
        box[Ndim - 1] = box[0];
    }

    double phi_modifier = 1.0;
    if (walls[0] < 0); { phi_modifier = sphereVolume(1 / 2., -walls[0]); }

    std::mt19937 rng(123456u);
    auto rawD = generate_diameter_distribution(N, dist, rng);
    auto scaled = scale_diametersND(rawD, phi * phi_modifier, box, Ndim);
    std::vector<double> D_scaled = std::move(scaled.first);
    double factor = scaled.second;

    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
        [&](std::size_t a, std::size_t b) { return D_scaled[a] > D_scaled[b]; });

    std::vector<double> D_sorted(N);
    for (std::size_t i = 0; i < N; ++i) {
        D_sorted[i] = D_scaled[idx[i]];
    }
    D_scaled = std::move(D_sorted);

    if (fix_height)
    {
        double s = std::pow(h * D_scaled[0] / Llast, 1.0 / (Ndim - 1));
        double eps = std::pow(s, Ndim);

        for (std::size_t i = 0; i < D_scaled.size(); ++i) {
            D_scaled[i] *= s;
        }

        box[Ndim - 1] *= eps;
    }

    LegacyInitializerDeterministicResult result;
    result.diameters = std::move(D_scaled);
    result.box = std::move(box);
    result.walls = std::move(walls);
    result.diameter_scale_factor = factor;
    result.phi_modifier = phi_modifier;
    return result;
}

void test_distribution_parity() {
    struct DistributionCase {
        std::string name;
        std::size_t N;
        rcpgenerator::Distribution dist;
    };

    std::vector<DistributionCase> cases;

    {
        rcpgenerator::Distribution dist;
        dist.type = "mono";
        dist.d = 1.25;
        cases.push_back({"mono", 6, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "gaussian";
        dist.mu = 1.0;
        dist.sigma = 0.2;
        cases.push_back({"gaussian", 8, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "bigaussian";
        dist.mu1 = 1.0;
        dist.sigma1 = 0.2;
        dist.mu2 = 2.0;
        dist.sigma2 = 0.3;
        dist.p = 0.375;
        cases.push_back({"bigaussian", 8, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "bidisperse";
        dist.d1 = 0.8;
        dist.d2 = 1.7;
        dist.p = 0.4;
        cases.push_back({"bidisperse", 10, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "lognormal";
        dist.mu = 0.0;
        dist.sigma = 0.5;
        cases.push_back({"lognormal", 7, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "flat";
        dist.d_min = 0.2;
        dist.d_max = 1.6;
        cases.push_back({"flat", 9, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "powerlaw";
        dist.d_min = 0.3;
        dist.d_max = 2.0;
        dist.exponent = -2.5;
        cases.push_back({"powerlaw", 9, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "exponential";
        dist.d_min = 0.25;
        dist.d_max = 1.75;
        cases.push_back({"exponential", 9, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "weibull";
        dist.scale = 1.5;
        dist.shape = 2.25;
        cases.push_back({"weibull", 9, dist});
    }
    {
        rcpgenerator::Distribution dist;
        dist.type = "custom";
        dist.custom = {0.4, 0.6, 0.9, 1.1, 1.3};
        cases.push_back({"custom", 5, dist});
    }

    for (const auto& distribution_case : cases) {
        std::mt19937 legacy_rng(123456u);
        std::mt19937 new_rng(123456u);
        Distribution legacy_dist = to_legacy_distribution(distribution_case.dist);
        auto legacy_values = generate_diameter_distribution(distribution_case.N, legacy_dist, legacy_rng);
        auto new_values = rcpgenerator::generate_diameter_distribution(distribution_case.N, distribution_case.dist, new_rng);
        compare_vector("distribution " + distribution_case.name, legacy_values, new_values);
    }
}

void test_scale_parity() {
    std::vector<double> diameters = {0.4, 0.8, 1.2, 1.6};
    std::vector<double> box = {1.0, 1.5, 2.0};
    auto legacy_scaled = scale_diametersND(diameters, 0.17, box, 3);
    auto new_scaled = rcpgenerator::scale_diameters_nd(diameters, 0.17, box, 3);
    compare_vector("scale_diameters_nd values", legacy_scaled.first, new_scaled.first);
    compare_double("scale_diameters_nd factor", legacy_scaled.second, new_scaled.second);
}

void test_parse_walls_parity() {
    std::vector<double> legacy_box = {2.0, 3.0, 4.0, 5.0};
    std::vector<double> new_box = legacy_box;
    auto legacy_walls = parse_walls("-1,true,false,0", 4, legacy_box);
    auto new_walls = rcpgenerator::parse_walls("-1,true,false,0", 4, new_box);
    compare_vector_i8("parse_walls walls", legacy_walls, new_walls);
    compare_vector("parse_walls box", legacy_box, new_box);
}

void test_initializer_deterministic_outputs() {
    {
        rcpgenerator::InitializerConfig config;
        config.phi = 0.11;
        config.N = 4;
        config.Ndim = 2;
        config.box = {1.0, 1.0};
        config.walls = {0, 0};
        config.fix_height = false;
        config.dist.type = "mono";
        config.dist.d = 1.0;

        auto legacy = run_legacy_initializer_deterministic(config);
        auto current = rcpgenerator::initialize_particles(config);
        compare_vector("initializer mono diameters", legacy.diameters, current.diameters);
        compare_vector("initializer mono box", legacy.box, current.box);
        compare_vector_i8("initializer mono walls", legacy.walls, current.walls);
        compare_double("initializer mono scale factor", legacy.diameter_scale_factor, current.diameter_scale_factor);
        compare_double("initializer mono phi modifier", legacy.phi_modifier, current.phi_modifier);
    }

    {
        rcpgenerator::InitializerConfig config;
        config.phi = 0.08;
        config.N = 3;
        config.Ndim = 3;
        config.box = {1.0, 1.0, 0.75};
        config.walls = {0, 0, 0};
        config.fix_height = true;
        config.dist.type = "custom";
        config.dist.custom = {0.4, 0.7, 1.1};

        auto legacy = run_legacy_initializer_deterministic(config);
        auto current = rcpgenerator::initialize_particles(config);
        compare_vector("initializer custom diameters", legacy.diameters, current.diameters);
        compare_vector("initializer custom box", legacy.box, current.box);
        compare_vector_i8("initializer custom walls", legacy.walls, current.walls);
        compare_double("initializer custom scale factor", legacy.diameter_scale_factor, current.diameter_scale_factor);
        compare_double("initializer custom phi modifier", legacy.phi_modifier, current.phi_modifier);
    }
}

}  // namespace

int main() {
    try {
        test_distribution_parity();
        test_scale_parity();
        test_parse_walls_parity();
        test_initializer_deterministic_outputs();
        std::cout << "Initializer parity checks passed.\n";
        std::cout << "Initializer position parity is not checked exactly because both legacy placement paths seed internal std::mt19937 instances from std::random_device.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Initializer parity failure: " << ex.what() << '\n';
        return 1;
    }
}
