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
    std::size_t Ndim,
    const std::vector<double>& box_skew)
{
    double V_box = 1.0;
    for (double L : box) V_box *= L;
    if (!box_skew.empty()) {
        std::vector<std::vector<double>> basis(box.size(), std::vector<double>(box.size(), 0.0));
        for (std::size_t axis = 0; axis < box.size(); ++axis) {
            const double theta = (axis < box_skew.size() ? box_skew[axis] : 0.0) * M_PI / 180.0;
            if (axis == 0) {
                basis[0][0] = box[0] * std::cos(theta);
                if (box.size() > 1) {
                    basis[1][0] = box[0] * std::sin(theta);
                }
            } else {
                basis[0][axis] = -box[axis] * std::sin(theta);
                basis[axis][axis] = box[axis] * std::cos(theta);
            }
        }
        double det = 1.0;
        int sign = 1;
        for (std::size_t col = 0; col < basis.size(); ++col) {
            std::size_t pivot = col;
            double pivot_abs = std::abs(basis[pivot][col]);
            for (std::size_t row = col + 1; row < basis.size(); ++row) {
                const double value_abs = std::abs(basis[row][col]);
                if (value_abs > pivot_abs) {
                    pivot = row;
                    pivot_abs = value_abs;
                }
            }
            if (pivot_abs <= 1e-14) {
                det = 0.0;
                break;
            }
            if (pivot != col) {
                std::swap(basis[pivot], basis[col]);
                sign *= -1;
            }
            det *= basis[col][col];
            for (std::size_t row = col + 1; row < basis.size(); ++row) {
                const double factor = basis[row][col] / basis[col][col];
                for (std::size_t k = col; k < basis.size(); ++k) {
                    basis[row][k] -= factor * basis[col][k];
                }
            }
        }
        V_box = std::abs(det * static_cast<double>(sign));
    }
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

struct CellGeometry {
    std::vector<std::vector<double>> basis;
    std::vector<std::vector<double>> inverse;
    std::vector<double> corrected_skew;
    double volume = 0.0;
    bool skewed = false;
};

double deg_to_rad(double value) {
    return value * M_PI / 180.0;
}

bool is_effectively_zero(double value) {
    return std::abs(value) < 1e-12;
}

bool is_skewed(const std::vector<double>& box_skew) {
    for (double value : box_skew) {
        if (!is_effectively_zero(value)) {
            return true;
        }
    }
    return false;
}

std::vector<double> mat_vec_mul(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector)
{
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        for (std::size_t j = 0; j < vector.size(); ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

double vec_dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

double vec_norm(const std::vector<double>& a) {
    return std::sqrt(vec_dot(a, a));
}

double angle_between_degrees(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    const double denom = vec_norm(a) * vec_norm(b);
    if (denom <= 1e-18) {
        return 0.0;
    }
    const double cosine = std::clamp(vec_dot(a, b) / denom, -1.0, 1.0);
    return std::acos(cosine) * 180.0 / M_PI;
}

std::vector<double> basis_vector_from_skew(
    const std::vector<double>& box,
    const std::vector<double>& box_skew,
    std::size_t axis)
{
    const std::size_t Ndim = box.size();
    std::vector<double> result(Ndim, 0.0);
    const double length = box[axis];
    const double theta = deg_to_rad(box_skew[axis]);

    if (axis == 0) {
        result[0] = length * std::cos(theta);
        if (Ndim > 1) {
            result[1] = length * std::sin(theta);
        }
        return result;
    }

    result[0] = -length * std::sin(theta);
    result[axis] = length * std::cos(theta);
    return result;
}

double determinant(std::vector<std::vector<double>> matrix) {
    const std::size_t N = matrix.size();
    double det = 1.0;
    int sign = 1;
    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivot = col;
        double pivot_abs = std::abs(matrix[pivot][col]);
        for (std::size_t row = col + 1; row < N; ++row) {
            const double value_abs = std::abs(matrix[row][col]);
            if (value_abs > pivot_abs) {
                pivot = row;
                pivot_abs = value_abs;
            }
        }
        if (pivot_abs <= 1e-14) {
            return 0.0;
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            sign *= -1;
        }
        det *= matrix[col][col];
        for (std::size_t row = col + 1; row < N; ++row) {
            const double factor = matrix[row][col] / matrix[col][col];
            for (std::size_t k = col; k < N; ++k) {
                matrix[row][k] -= factor * matrix[col][k];
            }
        }
    }
    return det * static_cast<double>(sign);
}

std::vector<std::vector<double>> invert_matrix(std::vector<std::vector<double>> matrix) {
    const std::size_t N = matrix.size();
    std::vector<std::vector<double>> inverse(N, std::vector<double>(N, 0.0));
    for (std::size_t i = 0; i < N; ++i) {
        inverse[i][i] = 1.0;
    }

    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivot = col;
        double pivot_abs = std::abs(matrix[pivot][col]);
        for (std::size_t row = col + 1; row < N; ++row) {
            const double value_abs = std::abs(matrix[row][col]);
            if (value_abs > pivot_abs) {
                pivot = row;
                pivot_abs = value_abs;
            }
        }
        if (pivot_abs <= 1e-14) {
            throw std::runtime_error("Skew cell basis is singular");
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(inverse[pivot], inverse[col]);
        }

        const double diag = matrix[col][col];
        for (std::size_t k = 0; k < N; ++k) {
            matrix[col][k] /= diag;
            inverse[col][k] /= diag;
        }

        for (std::size_t row = 0; row < N; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = matrix[row][col];
            for (std::size_t k = 0; k < N; ++k) {
                matrix[row][k] -= factor * matrix[col][k];
                inverse[row][k] -= factor * inverse[col][k];
            }
        }
    }

    return inverse;
}

CellGeometry build_cell_geometry(
    const std::vector<double>& box,
    const std::vector<double>& requested_box_skew,
    bool apply_guardrail,
    const char* warning_context)
{
    const std::size_t Ndim = box.size();
    CellGeometry geometry;
    geometry.corrected_skew = requested_box_skew;
    if (geometry.corrected_skew.empty()) {
        geometry.corrected_skew.assign(Ndim, 0.0);
    }
    if (geometry.corrected_skew.size() != Ndim) {
        throw std::runtime_error("box_skew length must match Ndim");
    }

    geometry.skewed = is_skewed(geometry.corrected_skew);
    geometry.basis.assign(Ndim, std::vector<double>(Ndim, 0.0));

    auto rebuild_basis = [&]() {
        for (std::size_t col = 0; col < Ndim; ++col) {
            const auto basis_vec = basis_vector_from_skew(box, geometry.corrected_skew, col);
            for (std::size_t row = 0; row < Ndim; ++row) {
                geometry.basis[row][col] = basis_vec[row];
            }
        }
    };

    rebuild_basis();

    if (apply_guardrail && Ndim >= 2) {
        const auto anchor = basis_vector_from_skew(box, geometry.corrected_skew, 0);
        for (std::size_t axis = 1; axis < Ndim; ++axis) {
            double angle = angle_between_degrees(anchor, basis_vector_from_skew(box, geometry.corrected_skew, axis));
            if (angle >= 30.0) {
                continue;
            }

            double theta = geometry.corrected_skew[axis];
            const auto trial_angle = [&](double new_theta) {
                auto trial_skew = geometry.corrected_skew;
                trial_skew[axis] = new_theta;
                return angle_between_degrees(anchor, basis_vector_from_skew(box, trial_skew, axis));
            };

            const double forward = trial_angle(theta + 0.1);
            const double backward = trial_angle(theta - 0.1);
            const double step = forward >= backward ? 0.1 : -0.1;

            std::size_t guard_steps = 0;
            while (angle < 30.0 && guard_steps < 4000) {
                theta += step;
                angle = trial_angle(theta);
                ++guard_steps;
            }

            if (angle < 30.0) {
                throw std::runtime_error("Unable to enforce the minimum skew basis angle");
            }

            std::cerr << warning_context
                      << ": requested box_skew[" << axis << "] produced a basis angle below 30 degrees; "
                      << "adjusted from " << geometry.corrected_skew[axis] << " to " << theta << ".\n";
            geometry.corrected_skew[axis] = theta;
            geometry.skewed = is_skewed(geometry.corrected_skew);
            rebuild_basis();
        }
    }

    geometry.volume = std::abs(determinant(geometry.basis));
    if (geometry.volume <= 1e-14) {
        throw std::runtime_error("Skew cell has near-zero volume");
    }
    geometry.inverse = invert_matrix(geometry.basis);
    return geometry;
}

std::vector<double> minimum_image_delta(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const CellGeometry& geometry)
{
    std::vector<double> delta(a.size(), 0.0);
    for (std::size_t d = 0; d < a.size(); ++d) {
        delta[d] = a[d] - b[d];
    }
    auto fractional = mat_vec_mul(geometry.inverse, delta);
    for (double& value : fractional) {
        value -= std::round(value);
    }
    return mat_vec_mul(geometry.basis, fractional);
}

std::vector<double> wrap_point(
    const std::vector<double>& point,
    const CellGeometry& geometry)
{
    auto fractional = mat_vec_mul(geometry.inverse, point);
    for (double& value : fractional) {
        value -= std::floor(value);
    }
    return mat_vec_mul(geometry.basis, fractional);
}

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

void initialize_positions_periodic_skew(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const CellGeometry& geometry)
{
    const std::size_t N = x.size();
    const std::size_t Ndim = geometry.basis.size();
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<> ud(0.0, 1.0);

    std::vector<double> fractional0(Ndim);
    for (std::size_t d = 0; d < Ndim; ++d) {
        fractional0[d] = ud(rng);
    }
    x[0] = mat_vec_mul(geometry.basis, fractional0);

    std::size_t count = 1;
    while (count < N) {
        std::vector<double> fractional(Ndim);
        for (std::size_t d = 0; d < Ndim; ++d) {
            fractional[d] = ud(rng);
        }
        std::vector<double> xt = mat_vec_mul(geometry.basis, fractional);
        bool ok = true;
        for (std::size_t i = 0; i < count; ++i) {
            const double r_ij = 0.5 * (D[i] + D[count]);
            const auto delta = minimum_image_delta(xt, x[i], geometry);
            if (vec_dot(delta, delta) <= r_ij * r_ij) {
                ok = false;
                break;
            }
        }
        if (ok) {
            x[count] = wrap_point(xt, geometry);
            ++count;
        }
    }
}

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
    std::vector<double> box_skew = config.box_skew;
    std::vector<std::int8_t> walls = config.walls;
    bool fix_height = config.fix_height;
    Distribution dist = config.dist;

    if (Ndim < 2) {
        throw std::runtime_error("--Ndim must be >=2");
    }
    if (box.empty()) box.assign(Ndim, 1.0);
    else if (box.size() != Ndim) {
        throw std::runtime_error("box length must match Ndim");
    }
    if (box_skew.empty()) box_skew.assign(Ndim, 0.0);
    if (box_skew.size() != Ndim) {
        throw std::runtime_error("box_skew length must match Ndim");
    }
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

    if (is_skewed(box_skew)) {
        const bool all_periodic = std::all_of(
            walls.begin(),
            walls.end(),
            [](std::int8_t wall) { return wall == 0; });
        if (!all_periodic) {
            throw std::runtime_error("Non-zero box_skew is currently supported only for fully periodic boundaries");
        }
        if (fix_height) {
            throw std::runtime_error("Non-zero box_skew is currently unsupported with fix_height");
        }
    }

    CellGeometry geometry = build_cell_geometry(box, box_skew, true, "rcpgenerator initialize_particles");
    box_skew = geometry.corrected_skew;

    std::mt19937 rng{std::random_device{}()};
    auto rawD = generate_diameter_distribution(N, dist, rng);
    auto scaled_pair = scale_diameters_nd(rawD, phi * phi_modifier, box, Ndim, box_skew);
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
    if (geometry.skewed) {
        geometry = build_cell_geometry(box, box_skew, false, "rcpgenerator initialize_particles");
    }

    std::vector<std::vector<double>> x(N, std::vector<double>(Ndim));
    if (geometry.skewed) {
        initialize_positions_periodic_skew(x, D_scaled, geometry);
    } else {
        initialize_positions_binned(x, D_scaled, box, walls, fix_height);
    }

    InitializerResult result;
    result.positions = std::move(x);
    result.diameters = std::move(D_scaled);
    result.box = std::move(box);
    result.box_skew = std::move(box_skew);
    result.walls = std::move(walls);
    result.diameter_scale_factor = factor;
    result.phi_modifier = phi_modifier;
    return result;
}

}  // namespace rcpgenerator
