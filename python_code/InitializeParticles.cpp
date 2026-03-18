#include "rcpgenerator/initialize_particles.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rcpgenerator {

std::vector<double> generate_diameter_distribution(
    std::size_t N,
    const Distribution& dist,
    std::mt19937& rng)
{
    std::vector<double> D(N);
    std::uniform_real_distribution<> ur(0.0, 1.0);

    if (dist.type == "mono") {
        std::fill(D.begin(), D.end(), dist.d);
    }
    else if (dist.type == "gaussian") {
        std::normal_distribution<> nd(dist.mu, dist.sigma);
        for (auto& v : D) v = std::abs(nd(rng));
    }
    else if (dist.type == "bigaussian") {
        std::size_t N1 = static_cast<std::size_t>(std::round(dist.p * N));
        std::normal_distribution<> g1(dist.mu1, dist.sigma1);
        std::normal_distribution<> g2(dist.mu2, dist.sigma2);
        for (std::size_t i = 0; i < N1; ++i) D[i] = std::abs(g1(rng));
        for (std::size_t i = N1; i < N; ++i) D[i] = std::abs(g2(rng));
        std::shuffle(D.begin(), D.end(), rng);
    }
    else if (dist.type == "bidisperse") {
        std::size_t N1 = static_cast<std::size_t>(std::round(dist.p * N));
        for (std::size_t i = 0; i < N1; ++i) D[i] = dist.d1;
        for (std::size_t i = N1; i < N; ++i) D[i] = dist.d2;
        std::shuffle(D.begin(), D.end(), rng);
    }
    else if (dist.type == "lognormal") {
        std::lognormal_distribution<> ld(dist.mu, dist.sigma);
        for (auto& v : D) {
            double tmp = ld(rng);
            v = std::clamp(tmp, 0.0025, 30.0);
        }
    }
    else if (dist.type == "flat") {
        std::uniform_real_distribution<> ud(dist.d_min, dist.d_max);
        for (auto& v : D) v = ud(rng);
        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {
            for (std::size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        }
        else {
            std::fill(D.begin(), D.end(), dist.d_min);
        }
    }
    else if (dist.type == "powerlaw") {
        std::uniform_real_distribution<> ud(0.0, 1.0);
        double a = dist.exponent + 1.0;
        for (std::size_t i = 0; i < N; ++i) {
            double u = ud(rng);
            if (std::abs(a) < 1e-12) {
                D[i] = dist.d_min * std::exp(u * std::log(dist.d_max / dist.d_min));
            }
            else {
                double top = std::pow(dist.d_max, a) - std::pow(dist.d_min, a);
                D[i] = std::pow(top * u + std::pow(dist.d_min, a), 1.0 / a);
            }
        }

        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {
            for (std::size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        }
        else {
            std::fill(D.begin(), D.end(), dist.d_min);
        }
    }
    else if (dist.type == "exponential") {
        double lambda = -std::log(0.5) / (dist.d_max - dist.d_min);
        std::uniform_real_distribution<> ud(0.0, 1.0);
        for (std::size_t i = 0; i < N; ++i) {
            double u = ud(rng);
            D[i] = dist.d_min - (1.0 / lambda) * std::log(1.0 - u);
            if (D[i] > dist.d_max) D[i] = dist.d_max;
        }

        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {
            for (std::size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        }
        else {
            std::fill(D.begin(), D.end(), dist.d_min);
        }
    }
    else if (dist.type == "weibull") {
        std::weibull_distribution<> wd(dist.shape, dist.scale);
        for (auto& v : D) v = wd(rng);
    }
    else if (dist.type == "custom") {
        if (dist.custom.size() != N)
            throw std::runtime_error("custom distribution needs N entries");
        D = dist.custom;
    }
    else {
        std::cerr << "Unsupported distribution: " << dist.type << std::endl;
        std::exit(1);
    }
    (void)ur;
    return D;
}

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim)
{
    double V_box = 1.0;
    for (double L : box) V_box *= L;
    double sumVol = 0.0;
    for (double Di : D) {
        double r = Di * 0.5;
        double vol_unit = std::pow(M_PI, Ndim / 2.0) / std::tgamma(Ndim / 2.0 + 1.0);
        sumVol += vol_unit * std::pow(r, double(Ndim));
    }
    double factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim));
    std::vector<double> scaled = D;
    for (auto& v : scaled) v *= factor;
    return {scaled, factor};
}

namespace {

double sqDist(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        s += (a[i] - b[i]) * (a[i] - b[i]);
    return s;
}

std::string trim(const std::string& s) {
    auto wsfront = std::find_if_not(s.begin(), s.end(),
                                    [](int c) { return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(),
                                   [](int c) { return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

}  // namespace

void initialize_positions_binned(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    std::vector<std::int8_t>& walls,
    bool fix_height)
{
    std::size_t N = x.size();
    std::size_t Ndim = box.size();
    std::size_t NB = 1u << Ndim;
    std::vector<std::vector<std::size_t>> bins(NB);

    std::vector<double> bin_size(Ndim);
    for (std::size_t d = 0; d < Ndim; ++d)
        bin_size[d] = box[d] * 0.5;

    auto getBin = [&](const std::vector<double>& pt) {
        std::size_t id = 0;
        for (std::size_t d = 0; d < Ndim; ++d)
            if (pt[d] > box[d] * 0.5) id |= (1u << d);
        return id;
    };

    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<> ud(0.0, 1.0);

    std::vector<double> xt0(Ndim);
    std::size_t bin_id0 = rng() % NB;
    for (std::size_t d = 0; d < Ndim; ++d) {
        std::size_t in_bin = (bin_id0 >> d) & 1;
        double origin = in_bin * bin_size[d];
        xt0[d] = origin + D[0] / 2 + ud(rng) * (bin_size[d] - D[0]);
        x[0][d] = xt0[d];
    }
    bins[bin_id0].push_back(0);

    std::size_t count = 1;
    while (count < N) {
        std::vector<double> xt(Ndim);
        std::size_t bin_id = rng() % NB;
        for (std::size_t d = 0; d < Ndim; ++d) {
            std::size_t in_bin = (bin_id >> d) & 1;
            double origin = in_bin * bin_size[d];
            xt[d] = origin + D[count] / 2 + ud(rng) * (bin_size[d] - D[count]);
        }

        std::size_t b = getBin(xt);
        bool ok = true;

        if (walls[0] < 0) {
            std::size_t K = static_cast<std::size_t>(-walls[0]);
            double R = box[0] * 0.5;
            double sumsq = 0.0;
            for (std::size_t d = 0; d < K; ++d) {
                double diff = xt[d] - R;
                sumsq += diff * diff;
            }
            if (sumsq > R * R) {
                ok = false;
            }
        }

        if (ok)
        {
            for (std::size_t idx : bins[b]) {
                double r_ij = 0.5 * (D[idx] + D[count]);
                if (sqDist(x[idx], xt) <= r_ij * r_ij) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                x[count] = std::move(xt);
                bins[b].push_back(count);
                ++count;
            }
        }
    }

    (void)fix_height;
}

void initialize_positions_naive(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    std::vector<std::int8_t>& walls,
    bool fix_height)
{
    std::size_t N = x.size(), Ndim = box.size();
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<> ud(0.0, 1.0);
    for (std::size_t d = 0; d < Ndim; ++d) x[0][d] = ud(rng) * (box[d] - D[0]) + D[0] / 2;
    std::size_t count = 1;
    while (count < N) {
        std::vector<double> xt(Ndim);
        for (std::size_t d = 0; d < Ndim; ++d) xt[d] = ud(rng) * (box[d] - D[count]) + D[count] / 2;
        bool ok = true;

        if (walls[0] < 0) {
            std::size_t K = static_cast<std::size_t>(-walls[0]);
            double R = box[0] * 0.5;
            double sumsq = 0.0;
            for (std::size_t d = 0; d < K; ++d) {
                double diff = xt[d] - R;
                sumsq += diff * diff;
            }
            if (sumsq > (R - D[count] / 2) * (R - D[count] / 2)) {
                ok = false;
            }
        }

        if (ok)
        {
            for (std::size_t i = 0; i < count; ++i) {
                double r_ij = 0.5 * (D[i] + D[count]);
                if (sqDist(x[i], xt) <= r_ij * r_ij) { ok = false; break; }
            }
            if (ok) { x[count] = std::move(xt); ++count; }
        }
    }

    (void)fix_height;
}

std::vector<std::int8_t> parse_walls(
    const std::string& s,
    std::size_t Ndim,
    std::vector<double>& box)
{
    std::vector<std::int8_t> walls;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trim(token);
        std::string lower = token;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (lower == "true") {
            walls.push_back(1);
        }
        else if (lower == "false") {
            walls.push_back(0);
        }
        else {
            try {
                int v = std::stoi(lower);
                if (v < -128 || v > 127) throw std::out_of_range("int8_t");
                walls.push_back(static_cast<std::int8_t>(v));
            }
            catch (...) {
                std::cerr << "Invalid --walls entry: \"" << token << "\"\n";
                std::exit(EXIT_FAILURE);
            }
        }
    }

    if (walls.empty()) {
        walls.assign(Ndim, 0);
    }

    if (walls[0] < 0) {
        if (walls[0] == -1) {
            walls[0] = -2;
        }
        std::size_t limit = static_cast<std::size_t>(-walls[0]);
        for (std::size_t k = 1; k <= limit && k < Ndim; ++k) {
            walls[k] = -1;
            box[k] = box[0];
        }
    }

    if (walls.size() < Ndim) {
        walls.resize(Ndim, 0);
    }

    return walls;
}

InitializerResult initialize_particles(const InitializerConfig& config) {
    double phi = config.phi;
    std::size_t N = config.N;
    std::size_t Ndim = config.Ndim;
    std::vector<double> box = config.box;
    std::vector<std::int8_t> walls = config.walls;
    bool fix_height = config.fix_height;
    Distribution dist = config.dist;

    if (Ndim < 2) {
        throw std::runtime_error("--Ndim must be >=2");
    }
    if (box.empty()) box.assign(Ndim, 1.0);
    if (phi <= 0 || phi >= 1 || N == 0) {
        throw std::runtime_error("Missing or invalid phi/N");
    }

    double Llast = box[0];
    double h = 1;
    if (fix_height)
    {
        h = box[Ndim - 1];
        box[Ndim - 1] = box[0];
    }

    double phi_modifier = 1.0;
    if (walls[0] < 0) {
        phi_modifier = sphere_volume(1 / 2., static_cast<std::size_t>(-walls[0]));
    }

    std::mt19937 rng{std::random_device{}()};
    auto rawD = generate_diameter_distribution(N, dist, rng);
    auto scaled_pair = scale_diameters_nd(rawD, phi * phi_modifier, box, Ndim);
    auto D_scaled = std::move(scaled_pair.first);
    double factor = scaled_pair.second;

    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
        [&](std::size_t a, std::size_t b) { return D_scaled[a] > D_scaled[b]; });

    std::vector<double> D_sorted(N);
    for (std::size_t i = 0; i < N; ++i)
        D_sorted[i] = D_scaled[idx[i]];
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

    std::vector<std::vector<double>> x(N, std::vector<double>(Ndim));
    initialize_positions_binned(x, D_scaled, box, walls, fix_height);

    InitializerResult result;
    result.positions = std::move(x);
    result.diameters = std::move(D_scaled);
    result.box = std::move(box);
    result.walls = std::move(walls);
    result.diameter_scale_factor = factor;
    result.phi_modifier = phi_modifier;
    return result;
}

}  // namespace rcpgenerator
