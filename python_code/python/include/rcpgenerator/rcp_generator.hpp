#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace rcpgenerator {

struct PackingConfig {
    std::vector<double> box;
    // Legacy packer behavior: empty walls defaults to Ndim zeros before first use.
    std::vector<std::int8_t> walls;
    std::uint32_t neighbor_max = 0;
    std::uint32_t seed = 0;
    bool fix_height = false;
};

struct PackingInput {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
};

struct PackingResult {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    std::vector<double> phi_history;
    std::vector<double> force_history;
    std::vector<double> energy_history;
    std::size_t steps = 0;
    double phi = 0.0;
    double max_min_dist = 0.0;
    double force_magnitude = 0.0;
};

struct PackingObservation {
    std::size_t step = 0;
    double phi = 0.0;
    double force_magnitude = 0.0;
    double energy = 0.0;
    double max_min_dist = 0.0;
};

struct PackingTrace {
    std::vector<std::size_t> steps;
    std::vector<std::vector<std::vector<double>>> positions;
    std::vector<std::vector<double>> diameters;
    std::vector<double> phi;
    std::vector<double> force;
    std::vector<double> energy;
    std::vector<double> max_min_dist;
};

struct PackingRunOptions {
    std::size_t max_steps = 0;
    bool has_mu_override = false;
    double mu_override = 0.0;
    bool fix_diameter = false;
    bool has_target_phi = false;
    double target_phi = 0.0;
};

using PackingObserver = std::function<void(
    const PackingObservation& observation,
    const std::vector<std::vector<double>>* positions)>;

double delta_x(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& x_old,
    const std::vector<double>& D);

double norm(
    const std::vector<double>& a,
    const std::vector<double>& b);

double compute_mean_force(
    const std::vector<std::vector<double>>& F,
    double Lc,
    const std::vector<std::size_t>& z,
    std::size_t Ndim);

double mean(
    const std::vector<double>& v,
    std::size_t start,
    std::size_t end);

double max_elem(const std::vector<std::vector<double>>& M);
double compute_max_min_dist(const std::vector<std::vector<double>>& min_dist);
std::size_t compute_max_neighbors(const std::vector<std::vector<std::uint32_t>>& pairs);
std::size_t compute_num_changes(const std::vector<std::uint32_t>& refresh);

std::vector<std::size_t> sort_indices_by_column(
    const std::vector<std::vector<double>>& x,
    std::size_t col);

#ifndef RCPGENERATOR_SPHERE_VOLUME_INLINE
#define RCPGENERATOR_SPHERE_VOLUME_INLINE
inline double sphere_volume(double r, std::size_t Ndim) {
    return std::pow(M_PI, Ndim / 2.0) * std::pow(r, Ndim) /
           std::tgamma(Ndim / 2.0 + 1.0);
}
#endif

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim,
    bool fix_height);

void get_pairs_nd_3(
    std::size_t N,
    std::size_t Ndim,
    const std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    const std::vector<std::uint32_t>& refresh,
    std::size_t neighborMax,
    std::vector<std::vector<std::uint32_t>>& pairs);

void get_forces_nd_3(
    const std::vector<std::vector<std::uint32_t>>& pairs,
    const std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    std::vector<std::vector<double>>& F,
    double& U,
    std::vector<std::vector<double>>& min_dist,
    double& max_min_dist,
    double& Lc,
    double& Fmean,
    double mu,
    double& dkappa,
    std::vector<std::size_t>& z);

void adam_update(
    const std::string& method,
    std::size_t N,
    std::size_t Ndim,
    const std::vector<std::vector<double>>& F,
    double dkappa,
    double beta1,
    double beta2,
    std::size_t t,
    double dt,
    double verlet_drag,
    std::vector<std::vector<double>>& m,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& v_max,
    std::vector<std::vector<double>>& v_update,
    std::vector<std::vector<double>>& m_hat,
    std::vector<std::vector<double>>& v_hat,
    std::vector<std::vector<double>>& a,
    std::vector<std::vector<double>>& v_verlet,
    std::vector<std::vector<double>>& a_old,
    double& m_kappa,
    double& v_kappa,
    double& v_update_kappa,
    double& m_hat_kappa,
    double& v_hat_kappa);

PackingResult run_packing(
    const PackingInput& input,
    const PackingConfig& config);

std::pair<PackingResult, PackingTrace> run_packing_observed(
    const PackingInput& input,
    const PackingConfig& config,
    std::size_t progress_interval,
    bool capture_positions,
    std::size_t trajectory_interval,
    const PackingObserver& observer,
    const PackingRunOptions& options = PackingRunOptions{});

}  // namespace rcpgenerator
