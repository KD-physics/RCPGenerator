#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace rcpgenerator {

struct Distribution {
    std::string type = "mono";
    double d = 1.0;
    double mu = 1.0, sigma = 0.2;
    double mu1 = 1.0, sigma1 = 0.2;
    double mu2 = 2.0, sigma2 = 0.3;
    double p = 0.5;
    double d1 = 1.0, d2 = 2.0;
    double d_min = 0.5, d_max = 1.5;
    double exponent = -2.5;
    double scale = 1.0, shape = 2.0;
    std::vector<double> custom;
};

struct InitializerConfig {
    double phi = 0.05;
    std::size_t N = 0;
    std::size_t Ndim = 0;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    bool fix_height = false;
    Distribution dist;
};

struct InitializerResult {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    double diameter_scale_factor = 1.0;
    double phi_modifier = 1.0;
};

std::vector<double> generate_diameter_distribution(
    std::size_t N,
    const Distribution& dist,
    std::mt19937& rng);

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim);

void initialize_positions_binned(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    std::vector<std::int8_t>& walls,
    bool fix_height);

void initialize_positions_naive(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    std::vector<std::int8_t>& walls,
    bool fix_height);

std::vector<std::int8_t> parse_walls(
    const std::string& s,
    std::size_t Ndim,
    std::vector<double>& box);

#ifndef RCPGENERATOR_SPHERE_VOLUME_INLINE
#define RCPGENERATOR_SPHERE_VOLUME_INLINE
inline double sphere_volume(double r, std::size_t Ndim) {
    return std::pow(M_PI, Ndim / 2.0) * std::pow(r, Ndim) /
           std::tgamma(Ndim / 2.0 + 1.0);
}
#endif

InitializerResult initialize_particles(const InitializerConfig& config);

}  // namespace rcpgenerator
