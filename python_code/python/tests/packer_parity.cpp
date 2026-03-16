#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "rcpgenerator/rcp_generator.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#define main legacy_rcp_generator_main
#include "../../RCPGenerator.cpp"
#undef main

namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

bool doubles_equal(double expected, double actual) {
    return (std::isnan(expected) && std::isnan(actual)) || expected == actual;
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

void compare_matrix(const std::string& label,
                    const std::vector<std::vector<double>>& expected,
                    const std::vector<std::vector<double>>& actual) {
    if (expected.size() != actual.size()) {
        std::ostringstream oss;
        oss << label << " row count mismatch: expected=" << expected.size()
            << " actual=" << actual.size();
        fail(oss.str());
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (expected[i].size() != actual[i].size()) {
            std::ostringstream oss;
            oss << label << "[" << i << "] column count mismatch: expected="
                << expected[i].size() << " actual=" << actual[i].size();
            fail(oss.str());
        }
        for (std::size_t j = 0; j < expected[i].size(); ++j) {
            if (!doubles_equal(expected[i][j], actual[i][j])) {
                std::ostringstream oss;
                oss << std::setprecision(17)
                    << label << "[" << i << "][" << j << "] mismatch: expected="
                    << expected[i][j] << " actual=" << actual[i][j];
                fail(oss.str());
            }
        }
    }
}

struct LegacyPackingResult {
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

LegacyPackingResult run_legacy_packing(
    const rcpgenerator::PackingInput& input,
    const rcpgenerator::PackingConfig& config) {
    std::size_t N = input.positions.size();
    std::vector<std::vector<double>> x = input.positions;
    std::vector<double> D = input.diameters;
    std::vector<double> D0;
    std::size_t Ndim = x.empty() ? 0 : x[0].size();
    std::vector<std::vector<double>> x_last = x;

    std::vector<double> box = config.box;
    std::vector<std::int8_t> walls = config.walls;
    std::uint32_t neighborMax = config.neighbor_max;
    std::uint32_t seed = config.seed;
    bool fix_height = config.fix_height;

    if (box.empty()) box.assign(Ndim, 1.0);
    else if (box.size() != Ndim) {
        fail("legacy packer parity setup: invalid box length");
    }
    if (walls.empty()) {
        walls.assign(Ndim, 0);
    }
    else if (walls.size() != Ndim) {
        fail("legacy packer parity setup: invalid walls length");
    }

    double phi0 = 0.025;
    double phi = phi0;
    double phi_modifier = 1;
    if (walls[0] < 0)
    {
        if (walls[0] == -1) { walls[0] = -2; }
        for (std::size_t k = 1; k <= static_cast<std::size_t>(-walls[0]); ++k) { walls[k] = -1; box[k] = box[0]; }

        phi_modifier = sphereVolume(1 / 2., -walls[0]);
    }

    if (fix_height)
    {
        box[Ndim - 1] = box[Ndim - 1] * D[0];
    }

    if (Ndim == 2) phi = 0.575;
    else if (Ndim == 3) phi = 0.33;
    else if (Ndim == 4) phi = 0.15;
    else if (Ndim == 5) phi = 0.10;
    else if (Ndim > 5) phi = 0.10 / std::pow(2, (Ndim - 4));
    double delta_phi0 = DELTA_PHI0;
    double delta_phi = DELTA_PHI0;
    double dphi = 0.0;
    std::size_t count = 0;
    bool update_flag = false;
    int direction_flag = 0;
    double F_magnitude = 0.0;
    auto scaled = scale_diametersND(D, phi * phi_modifier, box, Ndim, fix_height);
    std::vector<double> D_scaled = std::move(scaled.first);
    double factor = scaled.second;
    if (fix_height) {
        box[Ndim - 1] = box[Ndim - 1] * factor;
        for (std::size_t i = 0; i < N; ++i)
        {
            x[i][Ndim - 1] = x[i][Ndim - 1] * factor;
        }
    }
    D = D_scaled;
    double coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
    double current_volume = 0.0;
    for (std::size_t i = 0; i < D.size(); ++i) {
        current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
    }

    double box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());
    phi = current_volume / box_volume / phi_modifier;
    phi0 = phi;

    double Dmin = 0;
    auto it = std::min_element(D.begin(), D.end());
    Dmin = *std::min_element(D.begin(), D.end());
    double Dmin_last = Dmin;
    double max_min_dist = 0;

    if (neighborMax == 0) {
        std::size_t idx = (Ndim >= 2 ? std::min<std::size_t>(Ndim - 2, MAX_NEIGHBORS.size() - 1) : 0);
        neighborMax = MAX_NEIGHBORS[idx];
    }

    std::vector<std::uint32_t> refresh(N, 1);
    std::vector<std::vector<std::uint32_t>> pairs(
        N, std::vector<std::uint32_t>(neighborMax + 1, 0)
    );

    std::vector<std::vector<double>> F(N, std::vector<double>(Ndim, 0.0));
    std::vector<std::vector<double>> min_dist(
        N, std::vector<double>(neighborMax, 0.0)
    );
    std::vector<std::size_t> z(N, 0);
    double U = 0.0, Lc = 0.0;

    std::string method = METHOD;
    double beta1 = BETA1, beta2 = BETA2, t = 0, R = 1.00005, alpha = ALPHA_MAX, alpha_max = ALPHA_MAX, N_steps = N_STEPS;
    double dt = DT;
    double verlet_drag = 2;
    std::vector<std::vector<double>>
        m(N, std::vector<double>(Ndim, 0.0)),
        v(N, std::vector<double>(Ndim, 0.0)),
        v_max(N, std::vector<double>(Ndim, 0.0)),
        v_update(N, std::vector<double>(Ndim, 0.0)),
        m_hat(N, std::vector<double>(Ndim, 0.0)),
        v_hat(N, std::vector<double>(Ndim, 0.0)),
        a(N, std::vector<double>(Ndim, 0.0)),
        v_verlet(N, std::vector<double>(Ndim, 0.0)),
        a_old(N, std::vector<double>(Ndim, 0.0));

    double m_kappa = 0, v_kappa = 0, v_update_kappa = 0, m_hat_kappa = 0, v_hat_kappa = 0;

    std::size_t LastPhiUpdate = 0;
    std::size_t LastAlphaUpdate = 0;
    std::vector<double> U_history(N_steps, 0.0);

    std::vector<double> phi_history(N_steps, 0.0);
    std::vector<double> F_history(N_steps, 0.0);

    double phi_min = 0.8;
    if (Ndim == 2) { phi_min = 0.76; }
    if (Ndim == 3) { phi_min = 0.45; }
    if (Ndim == 4) { phi_min = 0.15; }
    if (Ndim == 5) { phi_min = 0.1; }
    if (Ndim > 5) { phi_min = 0.1 / std::pow(2, (Ndim - 4)); }

    double dkappa = 1.0;
    double kappa = 1.0;
    double mu = 5E-4;
    int mu_flag = 1;
    double Fmean = 0.0;
    double mu_change = 0.;

    std::pair<double, std::size_t> phi_max = {0.0, 0};

    D0 = D;

    alpha = 0.005;
    alpha_max = 0.005;
    if (Ndim == 2) { alpha = 0.005; alpha_max = 0.005; }
    if (Ndim == 3) { alpha = 0.0045; alpha_max = 0.0045; }
    if (Ndim == 4) { alpha = 0.0035; alpha_max = 0.0035; }
    if (Ndim == 5) { alpha = 0.0025; alpha_max = 0.0025; }
    if (Ndim > 5) { alpha = 0.0025; alpha_max = 0.0025; }

    GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
    GetForcesND_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);

    std::size_t steps = 0;
    for (std::size_t step = 1; step <= N_steps; ++step) {
        steps = step;

        t += 1;

        if ((step - LastAlphaUpdate > 500) &&
            (method == "ADAM") &&
            ((phi - phi_history[std::max<std::size_t>(step - 250, 1)]) < 5e-4) &&
            (step > 2500) && step > mu_change + 500) {

            double trend1 = mean(F_history, step - 450, step - 250) / mean(F_history, step - 75, step - 1);
            double trend2 = (phi - phi_history[std::max<std::size_t>(step - 250, 1)]);

            if ((mu_flag < 1 && trend1 < 0.85) || trend2 < -5E-4) {
                LastAlphaUpdate = step;
                alpha /= 1.5;
            }
            else if (mu_flag < 1 && step - LastAlphaUpdate > 1500) {
                trend1 = mean(F_history, step - 200, step - 1) / mean(F_history, step - 3000, step - 2500);
                if (trend1 > 0.5) {
                    LastAlphaUpdate = step;
                    alpha *= 1.15;
                }
            }

            if (mu_flag == 1) {
                alpha = std::max(alpha, alpha_max / 10.0);
            }
        }

        alpha = std::min(alpha * R, alpha_max);

        if (step % 500 == 0) {
            alpha = std::min(alpha_max, 1.1 * alpha);
        }

        for (std::size_t i = 0; i < D0.size(); ++i) {
            D[i] = D0[i] * kappa;
        }
        double Dmin_lastest = *std::min_element(D.begin(), D.end());

        coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
        current_volume = 0.0;
        for (std::size_t i = 0; i < D.size(); ++i) {
            current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
        }

        box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());

        if (fix_height) {
            box[Ndim - 1] = box[Ndim - 1] * kappa;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i][Ndim - 1] = x[i][Ndim - 1] * kappa;
            }
        }

        phi = current_volume / box_volume;
        phi_history[step - 1] = phi / phi_modifier;

        if (phi > phi_max.first) {
            phi_max = std::make_pair(phi, step);
        }

        if ((step > 500 && mu_flag == 1)) {
            std::vector<double> XX(phi_history.begin() + step - 250, phi_history.begin() + step - 1);
            double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

            if (delta_XX < 5e-6 || step - phi_max.second > 3500 || (step > std::max(15000, static_cast<int>(N / 3)))) {
                mu /= 10.0;
                alpha /= 2.0;
                mu_flag = 0;
                mu_change = step;
            }
        }

        if (mu_flag == 0 && step == 45000) { alpha = alpha / 10; }
        if (mu_flag == 0 && step == 55000) { alpha = alpha / 10; }
        if (mu_flag == 1 && step == 60000) { alpha = alpha / 10; }
        if (step > 500 && mu_flag == 0 && step > mu_change + 250) {
            std::vector<double> XX(phi_history.begin() + step - 250, phi_history.begin() + step - 1);
            double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

            if (delta_XX < 5e-6) {
                mu /= 10.0;
                alpha /= 2.0;
                mu_flag = -1;
                mu_change = step;
            }
        }

        if (Dmin_lastest / Dmin > 1.05) {
            std::fill(refresh.begin(), refresh.end(), 1u);
            GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
            x_last = x;
            Dmin = Dmin_lastest;
        } else {
            std::fill(refresh.begin(), refresh.end(), 0u);
            for (std::size_t k = 0; k < N; ++k) {
                if (norm(x[k], x_last[k]) > Dmin / 4) {
                    refresh[k] = 1;
                    x_last[k] = x[k];
                    update_flag = true;
                }
            }
            if (update_flag) {
                GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
            }
        }
        update_flag = false;

        GetForcesND_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);
        F_magnitude = computeMeanForce(F, Lc, z, Ndim) / Fmean;
        F_history[step - 1] = F_magnitude;

        AdamUpdate(method, N, Ndim,
                   F, dkappa, beta1, beta2, t, dt, verlet_drag,
                   m, v, v_max, v_update,
                   m_hat, v_hat,
                   a, v_verlet, a_old,
                   m_kappa, v_kappa, v_update_kappa,
                   m_hat_kappa, v_hat_kappa);

        if (step > 1100) {
            if (mu_flag == -1 &&
                F_history[step - 1] < 5.e-3 &&
                max_min_dist > 1e-16 &&
                step > mu_change + 250) {

                std::vector<double> XX(phi_history.begin() + step - 1000, phi_history.begin() + step - 1);
                double delta_XX = *std::max_element(XX.begin(), XX.end()) - *std::min_element(XX.begin(), XX.end());

                if (delta_XX < 2.5e-6) {
                    break;
                }
            }
        }

        kappa = kappa - alpha * m_hat_kappa / (std::sqrt(v_hat_kappa) + EPSILON);
        if (method != "Verlet") {
            for (std::size_t k = 0; k < N; ++k)
                for (std::size_t d = 0; d < Ndim; ++d)
                    x[k][d] -= alpha * m_hat[k][d] / (std::sqrt(v_hat[k][d]) + EPSILON);
        } else {
            for (std::size_t k = 0; k < N; ++k)
                for (std::size_t d = 0; d < Ndim; ++d)
                    x[k][d] += v_verlet[k][d] * dt + 0.5 * a_old[k][d] * dt * dt;
        }

        for (std::size_t k = 0; k < N; ++k)
            for (std::size_t d = 0; d < Ndim; ++d)
                x[k][d] = std::fmod(x[k][d] + box[d], box[d]);

        U_history[step - 1] = U;

        if (step % 5000 == 0) {
            std::this_thread::yield();
        }

        (void)delta_phi0;
        (void)delta_phi;
        (void)dphi;
        (void)count;
        (void)direction_flag;
        (void)phi0;
        (void)Dmin_last;
        (void)LastPhiUpdate;
        (void)phi_min;
        (void)seed;
        (void)it;
    }

    LegacyPackingResult result;
    result.positions = std::move(x);
    result.diameters = std::move(D);
    result.box = std::move(box);
    result.walls = std::move(walls);
    result.phi_history = std::move(phi_history);
    result.force_history = std::move(F_history);
    result.energy_history = std::move(U_history);
    result.steps = steps;
    result.phi = phi;
    result.max_min_dist = max_min_dist;
    result.force_magnitude = F_magnitude;
    return result;
}

void compare_packing_results(const std::string& case_name,
                             const LegacyPackingResult& legacy,
                             const rcpgenerator::PackingResult& current) {
    compare_vector(case_name + " phi_history", legacy.phi_history, current.phi_history);
    compare_vector(case_name + " force_history", legacy.force_history, current.force_history);
    compare_vector(case_name + " energy_history", legacy.energy_history, current.energy_history);
    compare_matrix(case_name + " final positions", legacy.positions, current.positions);
    compare_vector(case_name + " final diameters", legacy.diameters, current.diameters);
    compare_vector(case_name + " final box", legacy.box, current.box);
    compare_vector_i8(case_name + " final walls", legacy.walls, current.walls);
    if (legacy.steps != current.steps) {
        std::ostringstream oss;
        oss << case_name << " steps mismatch: expected=" << legacy.steps
            << " actual=" << current.steps;
        fail(oss.str());
    }
    compare_double(case_name + " final phi", legacy.phi, current.phi);
    compare_double(case_name + " final max_min_dist", legacy.max_min_dist, current.max_min_dist);
    compare_double(case_name + " final force_magnitude", legacy.force_magnitude, current.force_magnitude);
}

rcpgenerator::PackingInput make_case_2d_input() {
    rcpgenerator::PackingInput input;
    input.positions = {
        {0.10, 0.10},
        {0.35, 0.25},
        {0.70, 0.55},
        {0.82, 0.85}
    };
    input.diameters = {0.12, 0.10, 0.09, 0.08};
    return input;
}

rcpgenerator::PackingInput make_case_3d_input() {
    rcpgenerator::PackingInput input;
    input.positions = {
        {0.12, 0.15, 0.10},
        {0.46, 0.22, 0.31},
        {0.73, 0.68, 0.42},
        {0.28, 0.74, 0.81}
    };
    input.diameters = {0.10, 0.09, 0.08, 0.07};
    return input;
}

void run_case(const std::string& case_name,
              const rcpgenerator::PackingInput& input,
              const rcpgenerator::PackingConfig& config) {
    auto legacy = run_legacy_packing(input, config);
    auto current = rcpgenerator::run_packing(input, config);
    compare_packing_results(case_name, legacy, current);
}

}  // namespace

int main() {
    try {
        {
            rcpgenerator::PackingConfig config;
            config.box = {1.0, 1.0};
            config.walls = {};
            config.neighbor_max = 0;
            config.seed = 123u;
            config.fix_height = false;
            run_case("packer_2d_empty_walls_default", make_case_2d_input(), config);
        }

        {
            rcpgenerator::PackingConfig config;
            config.box = {1.0, 1.0, 0.8};
            config.walls = {0, 0, 0};
            config.neighbor_max = 64;
            config.seed = 999u;
            config.fix_height = true;
            run_case("packer_3d_fix_height", make_case_3d_input(), config);
        }

        std::cout << "Packer parity checks passed.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Packer parity failure: " << ex.what() << '\n';
        return 1;
    }
}
