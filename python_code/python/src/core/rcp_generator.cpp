#include "rcpgenerator/rcp_generator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace rcpgenerator {

namespace {

static constexpr double ALPHA_MAX = 0.0025;
static constexpr double BETA1 = 0.9;
static constexpr double BETA2 = 0.999;
static constexpr double EPSILON = 1e-8;
static constexpr std::size_t N_STEPS = 60000;
static constexpr double DT = 0.1;
static const std::string METHOD = "ADAM";
static const std::vector<std::uint32_t> MAX_NEIGHBORS = {300, 750, 5500, 5500};
static constexpr double DELTA_PHI0 = 1.5 * 1.5e-3;

}  // namespace

double delta_x(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& x_old,
    const std::vector<double>& D)
{
    std::size_t N = x.size();
    if (x_old.size() != N || D.size() != N) {
        throw std::invalid_argument("delta_x: mismatched sizes");
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double d2 = 0.0;
        for (std::size_t d = 0; d < x[i].size(); ++d) {
            double diff = x_old[i][d] - x[i][d];
            d2 += diff * diff;
        }
        double dist = std::sqrt(d2);
        sum += dist / D[i];
    }

    return sum / double(N);
}

double norm(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("norm: vector sizes differ");
    }
    double sum2 = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        sum2 += d * d;
    }
    return std::sqrt(sum2);
}

double mean(
    const std::vector<double>& v,
    std::size_t start,
    std::size_t end)
{
    double sum = 0.0;
    for (std::size_t k = start; k <= end; ++k) {
        sum += v[k];
    }
    return sum / double(end - start + 1);
}

double compute_mean_force(
    const std::vector<std::vector<double>>& F,
    double Lc,
    const std::vector<std::size_t>& z,
    std::size_t Ndim)
{
    std::size_t N = F.size();
    double sum_mag = 0.0;
    double sum_z = 0.0;

    std::uint32_t count = 0;

    for (std::size_t i = 0; i < N; ++i) {
        double d2 = 0.0;
        for (std::size_t d = 0; d < Ndim; ++d) {
            d2 += F[i][d] * F[i][d];
        }
        if (d2 > 1.0E-16) {
            sum_mag += std::sqrt(d2);
            ++count;
        }
    }

    double mean_mag = sum_mag / double(count);

    (void)Lc;
    (void)z;
    return mean_mag / std::sqrt(double(Ndim));
}

double max_elem(const std::vector<std::vector<double>>& M) {
    if (M.empty()) {
        throw std::invalid_argument("maxElem: empty matrix");
    }
    double m = M[0][0];
    for (const auto& row : M) {
        if (row.empty()) continue;
        double row_max = *std::max_element(row.begin(), row.end());
        if (row_max > m) m = row_max;
    }
    return m;
}

double compute_max_min_dist(const std::vector<std::vector<double>>& min_dist) {
    double m = 0.0;
    for (auto& row : min_dist)
        for (double v : row)
            if (v > m) m = v;
    return m;
}

std::size_t compute_max_neighbors(const std::vector<std::vector<std::uint32_t>>& pairs) {
    std::size_t m = 0;
    for (auto& row : pairs)
        if (row[0] > m) m = row[0];
    return m;
}

std::size_t compute_num_changes(const std::vector<std::uint32_t>& refresh) {
    std::size_t sum = 0;
    for (auto v : refresh)
        sum += v;
    return sum;
}

std::vector<std::size_t> sort_indices_by_column(
    const std::vector<std::vector<double>>& x,
    std::size_t col)
{
    std::size_t N = x.size();
    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b) { return x[a][col] < x[b][col]; });
    return idx;
}

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim,
    bool fix_height)
{
    double V_box = 1.0;
    for (double L : box) {
        V_box *= L;
    }

    double sumVol = 0.0;
    for (double Di : D) {
        double radius = Di * 0.5;
        sumVol += sphere_volume(radius, Ndim);
    }

    double factor = 1.0;
    if (fix_height)
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim - 1));
    }
    else
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim));
    }

    std::vector<double> D_scaled = D;
    for (double& Di : D_scaled) {
        Di *= factor;
    }

    return {D_scaled, factor};
}

void get_pairs_nd_3(
    std::size_t N,
    std::size_t Ndim,
    const std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    const std::vector<std::uint32_t>& refresh,
    std::size_t neighborMax,
    std::vector<std::vector<std::uint32_t>>& pairs)
{
    double Dmax = *std::max_element(D.begin(), D.end());
    double Dmin = *std::min_element(D.begin(), D.end());
    std::vector<double> Dsorted = D;
    std::sort(Dsorted.begin(), Dsorted.end());
    double median = Dsorted[N / 2];
    double t = std::min(std::max(median * 1.02, Dmin * 2.25), Dmax) * 0.5;

    auto sort_idx = sort_indices_by_column(x, 0);
    std::vector<std::size_t> sort_loc(N);
    for (std::size_t k = 0; k < N; ++k)
        sort_loc[sort_idx[k]] = k;

    bool warningIssued = false;

    for (std::size_t i = 0; i < N; ++i) {
        if (refresh[i] != 1) continue;
        pairs[i][0] = 0;

        double r_c_max = (D[i] + 1.01 * Dmax) / 2 + t;
        for (int dir : {-1, 1}) {
            std::size_t jdx = 0;
            bool go = true;
            while (go) {
                ++jdx;
                if (jdx > N / 2) { go = false; break; }

                std::size_t pos = (sort_loc[i] + dir * static_cast<int>(jdx) + N) % N;
                std::size_t j = sort_idx[pos];

                double dx = x[j][0] - x[i][0];
                dx -= std::round(dx / box[0]) * box[0];
                if (std::abs(dx) > r_c_max) { go = false; break; }

                double r_c = (D[i] + D[j]) / 2 + t;
                if (std::abs(dx) < r_c) {
                    bool proceed = true;
                    double d2 = dx * dx;
                    for (std::size_t d = 1; d < Ndim; ++d) {
                        double dz = x[j][d] - x[i][d];
                        dz -= std::round(dz / box[d]) * box[d];
                        if (std::abs(dz) > r_c) { proceed = false; }
                        d2 += dz * dz;
                    }
                    if (!proceed) continue;

                    std::uint32_t& cnt = pairs[i][0];
                    if (cnt < neighborMax) {
                        ++cnt;
                        pairs[i][cnt] = static_cast<std::uint32_t>(j);
                    } else if (!warningIssued) {
                        std::cerr << "Error: found more neighbors than allotted ("
                                  << neighborMax
                                  << "). Restart with a larger --NeighborMax.\n";
                        std::exit(EXIT_FAILURE);
                    }

                    if (cnt <= neighborMax) {
                        bool duplicate = false;
                        for (std::size_t k = 1; k <= pairs[j][0]; ++k)
                            if (pairs[j][k] == i) { duplicate = true; break; }
                        if (!duplicate) {
                            std::uint32_t& cntj = pairs[j][0];
                            if (cntj < neighborMax) {
                                ++cntj;
                                pairs[j][cntj] = static_cast<std::uint32_t>(i);
                            } else if (!warningIssued) {
                                std::cerr << "Error: found more neighbors than allotted ("
                                          << neighborMax
                                          << "). Restart with a larger --NeighborMax.\n";
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
    }

    (void)walls;
}

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
    std::vector<std::size_t>& z)
{
    std::size_t N = pairs.size();
    std::size_t Ndim = x[0].size();
    std::size_t Mplus1 = pairs[0].size();
    double K = 1.0;

    F.assign(N, std::vector<double>(Ndim, 0.0));
    min_dist.assign(N, std::vector<double>(Mplus1, 0.0));
    z.assign(N, 0);
    U = 0.0;
    Lc = 0.0;
    std::size_t count = 1;

    std::vector<double> dx(Ndim);

    max_min_dist = 0;
    bool circle_flag = (walls[0] < 0);

    bool no_walls = true;
    for (std::size_t d = 0; d < Ndim; ++d) {
        if (walls[d] != 0) {
            no_walls = false;
        }
    }

    dkappa = 0;
    Fmean = 0;

    for (std::size_t i = 0; i < N; ++i) {
        std::uint32_t numNbr = pairs[i][0];
        for (std::uint32_t jdx = 1; jdx <= numNbr; ++jdx) {
            std::size_t j = pairs[i][jdx];
            if (j <= i) continue;

            double r_ij = 0.5 * (D[i] + D[j]);
            bool flag = true;
            double d2 = 0.0;

            for (std::size_t d = 0; d < Ndim; ++d) {
                double delta = x[j][d] - x[i][d];
                delta -= std::floor(delta / box[d] + 0.5) * box[d];
                if (delta > r_ij | delta < -r_ij) {
                    flag = false;
                }
                dx[d] = delta;
                d2 += delta * delta;
            }

            if (flag) {
                double dist = std::sqrt(d2);
                if (dist < r_ij) {
                    z[i]++;
                    z[j]++;
                    double F_mag = -K * (r_ij / dist - 1.0);
                    U += (1.0 - dist / r_ij);
                    ++count;
                    min_dist[i][jdx] = (1.0 - dist / r_ij);
                    Lc += r_ij;
                    if (min_dist[i][jdx] > max_min_dist) { max_min_dist = min_dist[i][jdx]; }

                    for (std::size_t d = 0; d < Ndim; ++d) {
                        double fcomp = F_mag * dx[d];
                        F[i][d] += fcomp;
                        F[j][d] -= fcomp;
                    }

                    dkappa = dkappa + K * r_ij * (dist - r_ij);
                    Fmean = Fmean + F_mag;
                }
            }
        }
    }

    if (~no_walls) {

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t d = 0; d < Ndim; ++d) {
                if (walls[d] == 1)
                {
                    double r_ij = D[i] / 2;
                    for (std::size_t wall = 0; wall < 2; ++wall) {
                        dx[d] = box[d] * wall - x[i][d];
                        double dist = std::abs(dx[d]);
                        if (dist < r_ij) {
                            double F_mag = -2.0 * K * (r_ij / dist - 1.0);
                            U += 2.0 * (1.0 - dist / r_ij);
                            ++count;
                            double fcomp = F_mag * dx[d];
                            F[i][d] += fcomp;
                            dkappa = dkappa + 2 * K * r_ij * (dist - r_ij);
                            Fmean = Fmean + F_mag;
                        }
                    }
                }
            }

            if (circle_flag)
            {
                double r_ij = D[i] / 2;
                double d_ij = 0;
                double R = box[0] / 2;
                for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                {
                    d_ij = d_ij + (x[i][M] - R) * (x[i][M] - R);
                }
                d_ij = std::sqrt(d_ij);
                double delta = R - d_ij;
                if (delta < r_ij)
                {
                    delta = d_ij;
                    d_ij = 0;
                    for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                    {
                        dx[M] = std::abs(R / delta - 1) * (x[i][M] - R);
                        d_ij = d_ij + dx[M] * dx[M];
                    }
                    d_ij = sqrt(d_ij);

                    double F_mag = -2.0 * K * std::abs(r_ij / d_ij - 1.0);

                    U = U + (1 - d_ij / r_ij);
                    ++count;

                    for (std::size_t M = 0; M < static_cast<std::size_t>(-walls[0]); ++M)
                    {
                        F[i][M] = F[i][M] + F_mag * dx[M];
                    }
                    dkappa = dkappa + 2 * K * r_ij * (d_ij - r_ij);
                    Fmean = Fmean + F_mag;
                }
            }
        }
    }

    U = std::pow(U / double(count), 2);
    Lc = Lc / double(count);

    double sumD = 0.0;
    for (double Di : D) {
        sumD += Di;
    }

    dkappa = dkappa + mu * sumD;
    Fmean = -Fmean / count;
}

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
    double& v_hat_kappa)
{
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t kk = 0; kk < Ndim; ++kk) { std::size_t idx = k * Ndim + kk; double Fkd = F[k][kk];
            m[k][kk] = beta1 * m[k][kk] - (1.0 - beta1) * Fkd;
            v[k][kk] = beta2 * v[k][kk] + (1.0 - beta2) * Fkd * Fkd;
            v_update[k][kk] = v[k][kk];
            m_hat[k][kk] = m[k][kk] / (1.0 - std::pow(beta1, double(t)));
            v_hat[k][kk] = v_update[k][kk] / (1.0 - std::pow(beta2, double(t)));
            (void)idx;
        }
    }

    m_kappa = beta1 * m_kappa - (1 - beta1) * dkappa;
    v_kappa = beta2 * v_kappa + (1 - beta2) * (dkappa * dkappa);
    v_update_kappa = v_kappa;
    m_hat_kappa = m_kappa / (1.0 - std::pow(beta1, double(t)));
    v_hat_kappa = v_update_kappa / (1.0 - std::pow(beta2, double(t)));

    (void)method;
    (void)dt;
    (void)verlet_drag;
    (void)v_max;
    (void)a;
    (void)v_verlet;
    (void)a_old;
}

PackingResult run_packing(
    const PackingInput& input,
    const PackingConfig& config)
{
    return run_packing_observed(input, config, 0, false, 0, nullptr).first;
}

std::pair<PackingResult, PackingTrace> run_packing_observed(
    const PackingInput& input,
    const PackingConfig& config,
    std::size_t progress_interval,
    bool capture_positions,
    std::size_t trajectory_interval,
    const PackingObserver& observer,
    const PackingRunOptions& options)
{
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
        throw std::runtime_error("Error: --box length mismatch");
    }
    if (walls.empty()) {
        walls.assign(Ndim, 0);
    }
    else if (walls.size() != Ndim) {
        throw std::runtime_error("Error: --walls expects Ndim entries");
    }

    double phi0 = 0.025;
    double phi = phi0;
    double phi_modifier = 1;
    if (walls[0] < 0)
    {
        if (walls[0] == -1) { walls[0] = -2; }
        for (std::size_t k = 1; k < static_cast<std::size_t>(-walls[0]) && k < Ndim; ++k) { walls[k] = -1; box[k] = box[0]; }

        phi_modifier = sphere_volume(1 / 2., -walls[0]);
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
    if (!options.fix_diameter) {
        auto scaled_pair = scale_diameters_nd(D, phi * phi_modifier, box, Ndim, fix_height);
        auto D_scaled = std::move(scaled_pair.first);
        double factor = scaled_pair.second;
        if (fix_height) {
            box[Ndim - 1] = box[Ndim - 1] * factor;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i][Ndim - 1] = x[i][Ndim - 1] * factor;
            }
        }
        D = D_scaled;
    }
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
    double beta1 = BETA1, beta2 = BETA2, t = 0, R = 1.00005, alpha = ALPHA_MAX, alpha_max = ALPHA_MAX, N_steps = N_STEPS, epsilon = EPSILON;
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
    double mu = options.has_mu_override ? options.mu_override : 5E-4;
    bool runtime_fix_diameter = options.fix_diameter;
    bool target_phi_locked = false;
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

    get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);

    get_forces_nd_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);

    PackingTrace trace;
    auto record_observation = [&](std::size_t step_value, bool capture_sample_positions) {
        PackingObservation observation;
        observation.step = step_value;
        observation.phi = phi / phi_modifier;
        observation.force_magnitude = F_magnitude;
        observation.energy = U;
        observation.max_min_dist = max_min_dist;

        if (capture_sample_positions) {
            trace.steps.push_back(step_value);
            trace.diameters.push_back(D);
            trace.phi.push_back(observation.phi);
            trace.force.push_back(observation.force_magnitude);
            trace.energy.push_back(observation.energy);
            trace.max_min_dist.push_back(observation.max_min_dist);
            trace.positions.push_back(x);
        }

        if (observer) {
            observer(observation, capture_sample_positions ? &trace.positions.back() : nullptr);
        }
    };

    const bool capture_initial =
        capture_positions && trajectory_interval > 0;
    const bool report_initial =
        observer && progress_interval > 0;
    if (capture_initial || report_initial) {
        record_observation(0, capture_initial);
    }

    std::size_t steps = 0;
    const std::size_t max_steps = options.max_steps > 0 ? options.max_steps : static_cast<std::size_t>(N_steps);
    for (std::size_t step = 1; step <= max_steps; ++step) {
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

        if (!runtime_fix_diameter) {
            for (std::size_t i = 0; i < D0.size(); ++i) {
                D[i] = D0[i] * kappa;
            }
        }
        double Dmin_lastest = *std::min_element(D.begin(), D.end());

        coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
        current_volume = 0.0;
        for (std::size_t i = 0; i < D.size(); ++i) {
            current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
        }

        box_volume = std::accumulate(box.begin(), box.end(), 1.0, std::multiplies<double>());

        if (fix_height && !runtime_fix_diameter) {
            box[Ndim - 1] = box[Ndim - 1] * kappa;
            for (std::size_t i = 0; i < N; ++i)
            {
                x[i][Ndim - 1] = x[i][Ndim - 1] * kappa;
            }
        }

        phi = current_volume / box_volume;

        if (options.has_target_phi && !target_phi_locked && !runtime_fix_diameter) {
            const double user_phi = phi / phi_modifier;
            if (user_phi >= options.target_phi) {
                auto target_scaled = scale_diameters_nd(
                    D,
                    options.target_phi * phi_modifier,
                    box,
                    Ndim,
                    fix_height);
                D = std::move(target_scaled.first);
                D0 = D;
                kappa = 1.0;
                mu = 0.0;
                runtime_fix_diameter = true;
                target_phi_locked = true;

                Dmin_lastest = *std::min_element(D.begin(), D.end());
                coeff = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
                current_volume = 0.0;
                for (std::size_t i = 0; i < D.size(); ++i) {
                    current_volume += coeff * std::pow(D[i] / 2.0, Ndim);
                }
                phi = current_volume / box_volume;
                std::fill(refresh.begin(), refresh.end(), 1u);
                get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
                x_last = x;
                Dmin = Dmin_lastest;
            }
        }

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
            get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
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
                get_pairs_nd_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
            }
        }
        update_flag = false;

        get_forces_nd_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, Fmean, mu, dkappa, z);
        F_magnitude = compute_mean_force(F, Lc, z, Ndim) / Fmean;
        F_history[step - 1] = F_magnitude;

        adam_update(method, N, Ndim,
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

        if (!runtime_fix_diameter) {
            kappa = kappa - alpha * m_hat_kappa / (std::sqrt(v_hat_kappa) + epsilon);
        }
        auto x_old = x;
        if (method != "Verlet") {
            for (std::size_t k = 0; k < N; ++k)
                for (std::size_t d = 0; d < Ndim; ++d)
                    x[k][d] -= alpha * m_hat[k][d] / (std::sqrt(v_hat[k][d]) + epsilon);
        } else {
            for (std::size_t k = 0; k < N; ++k)
                for (std::size_t d = 0; d < Ndim; ++d)
                    x[k][d] += v_verlet[k][d] * dt + 0.5 * a_old[k][d] * dt * dt;
        }

        for (std::size_t k = 0; k < N; ++k)
            for (std::size_t d = 0; d < Ndim; ++d)
                x[k][d] = std::fmod(x[k][d] + box[d], box[d]);

        U_history[step - 1] = U;

        const bool should_report =
            observer && progress_interval > 0 && (step % progress_interval == 0);
        const bool should_capture =
            capture_positions && trajectory_interval > 0 && (step % trajectory_interval == 0);
        if (should_report || should_capture) {
            record_observation(step, should_capture);
        }

        if (step % 5000 == 0) {
            std::this_thread::yield();
        }

        (void)x_old;
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
    }

    PackingResult result;
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

    const bool final_already_captured =
        capture_positions && trajectory_interval > 0 && steps > 0 && (steps % trajectory_interval == 0);
    const bool final_already_reported =
        observer && progress_interval > 0 && steps > 0 && (steps % progress_interval == 0);
    if ((capture_positions && !final_already_captured) || (observer && !final_already_reported)) {
        trace.steps.push_back(steps);
        trace.diameters.push_back(result.diameters);
        trace.phi.push_back(phi);
        trace.force.push_back(F_magnitude);
        trace.energy.push_back(U);
        trace.max_min_dist.push_back(max_min_dist);
        if (capture_positions) {
            trace.positions.push_back(result.positions);
        }
        if (observer) {
            PackingObservation observation;
            observation.step = steps;
            observation.phi = phi;
            observation.force_magnitude = F_magnitude;
            observation.energy = U;
            observation.max_min_dist = max_min_dist;
            observer(observation, capture_positions ? &trace.positions.back() : nullptr);
        }
    }

    return {std::move(result), std::move(trace)};
}

}  // namespace rcpgenerator
