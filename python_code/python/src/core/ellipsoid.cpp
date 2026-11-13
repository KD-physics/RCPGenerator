// Ellipsoid extension: contact kernel, force pass, orientation update, exact shrink-to-contact.
// See include/rcpgenerator/ellipsoid.hpp for the conventions. Nothing here is specific to 2D or 3D.

#include "rcpgenerator/ellipsoid.hpp"
#include "rcpgenerator/rod.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include <omp.h>

#include "rcpgenerator/curved_wall.hpp"

#if defined(__SSE2__)
#include <emmintrin.h>
#endif

namespace rcpgenerator {

namespace {

constexpr std::size_t MD = ELL_MAX_DIM;

// The helpers below are force-inlined so that, when they are reached through the dimension-specialised
// entry points (contact_fixed<2>, <3>, <4>, ...), the dimension is a compile-time constant and every small
// loop is unrolled. The arithmetic and its order are the same for every dimension and for the generic path.
#if defined(__GNUC__) || defined(__clang__)
#define ELL_INLINE inline __attribute__((always_inline))
#else
#define ELL_INLINE inline
#endif

// ---- small dense helpers on fixed [MD][MD] storage, loops run to n ------------------------------------------

// Cholesky factor of the symmetric positive-definite n x n matrix Y (row-major, stride n) into L (lower).
// Linv[i] = 1 / L[i][i]: the solves multiply by it instead of dividing (divisions dominate the pair cost).
ELL_INLINE void chol(const double* Y, std::size_t n, double L[MD][MD], double Linv[MD]) {
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = Y[i * n + j];
            for (std::size_t k = 0; k < j; ++k) sum -= L[i][k] * L[j][k];
            if (i == j) { L[i][i] = std::sqrt(sum); Linv[i] = 1.0 / L[i][i]; }
            else        L[i][j] = sum * Linv[j];
        }
    }
}

// Solve (L L^T) out = b.
ELL_INLINE void chol_solve(const double L[MD][MD], const double Linv[MD], std::size_t n, const double* b,
                           double* out) {
    double y[MD];
    for (std::size_t i = 0; i < n; ++i) {
        double sum = b[i];
        for (std::size_t k = 0; k < i; ++k) sum -= L[i][k] * y[k];
        y[i] = sum * Linv[i];
    }
    for (std::size_t ii = n; ii-- > 0;) {
        double sum = y[ii];
        for (std::size_t k = ii + 1; k < n; ++k) sum -= L[k][ii] * out[k];
        out[ii] = sum * Linv[ii];
    }
}

ELL_INLINE void matvec(const double* A, std::size_t n, const double* v, double* out) {
    for (std::size_t a = 0; a < n; ++a) {
        double sum = 0.0;
        for (std::size_t b = 0; b < n; ++b) sum += A[a * n + b] * v[b];
        out[a] = sum;
    }
}

ELL_INLINE double dot(const double* a, const double* b, std::size_t n) {
    double sum = 0.0;
    for (std::size_t k = 0; k < n; ++k) sum += a[k] * b[k];
    return sum;
}

// One evaluation of the contact function at lam. Fills x = Y^-1 r, Ax, Bx and returns r.x.
// One evaluation of the contact function at lam: x = Y^-1 r, Ax, Bx, r.x, x.A x, x.B x, and the factorisation
// of Y (explicit symmetric inverse in 2D and 3D, Cholesky otherwise) for a second solve.
struct LamEval {
    double x[MD], Ax[MD], Bx[MD];
    double rx, p, q;       // r.x, x.A x, x.B x
    double L[MD][MD];
    double Linv[MD];
    double Yi[6];          // 2D: i00 i01 i11;  3D: i00 i01 i02 i11 i12 i22
};

ELL_INLINE void factor(const double* Y, std::size_t n, LamEval& e) {
    if (n == 2) {
        const double inv = 1.0 / (Y[0] * Y[3] - Y[1] * Y[1]);
        e.Yi[0] = Y[3] * inv;
        e.Yi[1] = -Y[1] * inv;
        e.Yi[2] = Y[0] * inv;
    } else if (n == 3) {
        const double c00 = Y[4] * Y[8] - Y[5] * Y[5];
        const double c01 = Y[2] * Y[5] - Y[1] * Y[8];
        const double c02 = Y[1] * Y[5] - Y[2] * Y[4];
        const double inv = 1.0 / (Y[0] * c00 + Y[1] * c01 + Y[2] * c02);
        e.Yi[0] = c00 * inv;
        e.Yi[1] = c01 * inv;
        e.Yi[2] = c02 * inv;
        e.Yi[3] = (Y[0] * Y[8] - Y[2] * Y[2]) * inv;
        e.Yi[4] = (Y[1] * Y[2] - Y[0] * Y[5]) * inv;
        e.Yi[5] = (Y[0] * Y[4] - Y[1] * Y[1]) * inv;
    } else {
        chol(Y, n, e.L, e.Linv);
    }
}

ELL_INLINE void solve(const LamEval& e, std::size_t n, const double* b, double* out) {
    if (n == 2) {
        out[0] = e.Yi[0] * b[0] + e.Yi[1] * b[1];
        out[1] = e.Yi[1] * b[0] + e.Yi[2] * b[1];
    } else if (n == 3) {
        out[0] = e.Yi[0] * b[0] + e.Yi[1] * b[1] + e.Yi[2] * b[2];
        out[1] = e.Yi[1] * b[0] + e.Yi[3] * b[1] + e.Yi[4] * b[2];
        out[2] = e.Yi[2] * b[0] + e.Yi[4] * b[1] + e.Yi[5] * b[2];
    } else {
        chol_solve(e.L, e.Linv, n, b, out);
    }
}

// Stage 1: x = Y^-1 r and r.x, enough to know G(lam) (and so whether the pair is certainly apart).
ELL_INLINE void eval_lam_x(const double* A, const double* B, const double* r, std::size_t n, double lam, LamEval& e) {
    double Y[MD * MD];
    const double wa = 1.0 - lam, wb = lam;
    for (std::size_t k = 0; k < n * n; ++k) Y[k] = wa * A[k] + wb * B[k];
    factor(Y, n, e);
    solve(e, n, r, e.x);
    e.rx = dot(r, e.x, n);
}

// Stage 2: Ax, Bx, x.A x, x.B x, needed for G', the torques and the next lam.
ELL_INLINE void eval_lam_rest(const double* A, const double* B, std::size_t n, LamEval& e) {
    matvec(A, n, e.x, e.Ax);
    matvec(B, n, e.x, e.Bx);
    e.p = dot(e.x, e.Ax, n);
    e.q = dot(e.x, e.Bx, n);
}

// The lam iteration. G(lam) = lam(1-lam) r.x is concave on [0,1] with G' = (1-lam)^2 p - lam^2 q, so the root
// of G' is bracketed by the signs seen so far.
//   Cold start (lam_start outside (0,1)): lam = rho_A/(rho_A+rho_B), rho = sqrt(r.A r) (exact for spheres);
//     evaluation 0 is followed by the fixed-point step lam = sqrt(p)/(sqrt(p)+sqrt(q)), later ones by Newton
//     steps kept inside the bracket (fallback: fixed-point step, then bisection).
//   Warm start (lam_start in (0,1), the pair's lam from the previous step): Newton from the first evaluation.
// At most n_iter evaluations plus one final evaluation: a fixed bound. Any lam gives G(lam) <= s^2.
//   early_out: return false as soon as G >= 1 (a lower bound on s^2 reached 1 => certainly apart).
//   track_tol > 0: accept the CURRENT evaluation, without evaluating again, once the Newton step predicts a
//     relative gain in G below track_tol (G'^2 / (2|G''|) <= track_tol * G). The step itself is returned in
//     lam_next for the next call. track_tol == 0: exactly n_iter + 1 evaluations for every pair in contact.
//   conv_tol > 0: iterate until the step is below conv_tol (reporting only).
// On return true, e is the evaluation at lam; lam_next is the best estimate of the maximiser.
ELL_INLINE bool lam_iterate(const double* A, const double* B, const double* r, std::size_t n, int n_iter,
                            bool early_out, double conv_tol, double lam_start, double track_tol,
                            double& lam, double& lam_next, LamEval& e) {
    const bool warm = lam_start > 0.0 && lam_start < 1.0;
    if (warm) {
        lam = lam_start;
    } else {
        double Ar[MD], Br[MD];
        matvec(A, n, r, Ar);
        matvec(B, n, r, Br);
        const double ra = std::sqrt(dot(r, Ar, n)), rb = std::sqrt(dot(r, Br, n));
        lam = ra / (ra + rb);
    }
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < n_iter; ++it) {
        eval_lam_x(A, B, r, n, lam, e);
        const double G = lam * (1.0 - lam) * e.rx;
        if (early_out && G >= 1.0) {
            lam_next = lam;      // this lam certifies "apart"; it is the start for the next step
            return false;
        }
        eval_lam_rest(A, B, n, e);
        const double G1 = (1.0 - lam) * (1.0 - lam) * e.p - lam * lam * e.q;
        if (G1 > 0.0) lo = lam; else hi = lam;
        double next = -1.0;
        double gain = std::numeric_limits<double>::infinity();
        if (warm || it >= 1) {
            double w[MD], y[MD];
            for (std::size_t k = 0; k < n; ++k) w[k] = e.Bx[k] - e.Ax[k];
            solve(e, n, w, y);
            const double q1 = dot(e.x, w, n), q2 = dot(w, y, n);
            const double G2 = -2.0 * e.rx - 2.0 * (1.0 - 2.0 * lam) * q1 + 2.0 * lam * (1.0 - lam) * q2;
            if (G2 < 0.0) {
                next = lam - G1 / G2;
                gain = 0.5 * G1 * G1 / (-G2);
            }
        }
        if (!(next >= lo && next <= hi)) {
            const double sp = std::sqrt(e.p), sq = std::sqrt(e.q);
            next = sp / (sp + sq);
            if (!(next >= lo && next <= hi)) next = 0.5 * (lo + hi);
            gain = std::numeric_limits<double>::infinity();
        }
        if (track_tol > 0.0 && gain <= track_tol * G) {
            lam_next = next;
            return true;                       // e is the evaluation at lam
        }
        const double step = std::abs(next - lam);
        lam = next;
        if (conv_tol > 0.0 && step <= conv_tol) break;
    }
    eval_lam_x(A, B, r, n, lam, e);
    eval_lam_rest(A, B, n, e);
    lam_next = lam;
    return true;
}

// Physics of one overlapping pair, given the contact evaluation. Adds to Fi/Fj (forces), Ti/Tj (torques).
struct PairPhysics {
    double overlap;   // 1 - s
    double sigma;     // contact distance along r
    double F_mag;     // -K (sigma/d - 1)   (the sphere code's F_mag)
    double dkappa;    // K sigma (d - sigma)
    double energy;    // (K/2)(sigma - d)^2
    double f_norm;    // |force on i|
};

ELL_INLINE PairPhysics pair_physics(const double* r, double d2, std::size_t n, double s2, double lam,
                                    const LamEval& e, double K, double* Fi, double* Fj, double* Ti, double* Tj) {
    const double d = std::sqrt(d2);
    const double s = std::sqrt(s2);
    const double inv_d = 1.0 / d, inv_s = 1.0 / s;
    const double sigma = d * inv_s;
    const double pref = K * (sigma - d);
    PairPhysics out;
    out.overlap = 1.0 - s;
    out.sigma = sigma;
    out.F_mag = -K * (sigma * inv_d - 1.0);
    out.dkappa = K * sigma * (d - sigma);
    out.energy = 0.5 * K * (sigma - d) * (sigma - d);
    // F_i = pref * (dsigma/dr - r/d),  dsigma/dr = r/(d s) - d lam(1-lam) x / s^3.
    const double inv_s3 = inv_s * inv_s * inv_s;
    const double ll = lam * (1.0 - lam);
    const double c_r = pref * (inv_d * inv_s - inv_d);
    const double c_x = -pref * d * ll * inv_s3;
    double f2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double f = c_r * r[k] + c_x * e.x[k];
        Fi[k] += f;
        Fj[k] -= f;
        f2 += f * f;
    }
    out.f_norm = std::sqrt(f2);
    // dE/dA_i = ci x x^T, dE/dA_j = cj x x^T;  T_ab = -dE/dOmega_ab = -2 c (x_a (Ax)_b - (Ax)_a x_b).
    const double base = pref * d * ll * inv_s3;       // = 2 x (pref d ll / (2 s^3))
    const double ci2 = -base * (1.0 - lam), cj2 = -base * lam;
    std::size_t c = 0;
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = a + 1; b < n; ++b, ++c) {
            Ti[c] += ci2 * (e.x[a] * e.Ax[b] - e.Ax[a] * e.x[b]);
            Tj[c] += cj2 * (e.x[a] * e.Bx[b] - e.Bx[a] * e.x[b]);
        }
    }
    return out;
}

// Dimension-specialised entry points (NN == 0: the generic, runtime-dimension instance).
template <std::size_t NN>
bool contact_fixed(const double* Ai, const double* Aj, const double* r, std::size_t n_runtime, int n_iter,
                   double& s2, double& lam, double* xvec) {
    const std::size_t n = NN ? NN : n_runtime;
    LamEval e;
    double lam_next;
    if (!lam_iterate(Ai, Aj, r, n, n_iter, true, 0.0, -1.0, 0.0, lam, lam_next, e)) return false;
    s2 = lam * (1.0 - lam) * e.rx;
    for (std::size_t k = 0; k < n; ++k) xvec[k] = e.x[k];
    return s2 < 1.0;
}

using ContactFn = bool (*)(const double*, const double*, const double*, std::size_t, int, double&, double&, double*);

ContactFn pick_contact(std::size_t n) {
    switch (n) {
        case 2: return contact_fixed<2>;
        case 3: return contact_fixed<3>;
        case 4: return contact_fixed<4>;
        case 5: return contact_fixed<5>;
        case 6: return contact_fixed<6>;
        default: return contact_fixed<0>;
    }
}

// The whole pair: contact evaluation (cold or warm) and, if in contact, its physics. Returns 0 = certainly
// apart, 1 = evaluated but s >= 1, 2 = contact (ph filled, forces and torques accumulated).
template <std::size_t NN>
int pair_full(const double* Ai, const double* Aj, const double* r, double d2, std::size_t n_runtime, int n_iter,
              double lam_start, double track_tol, double K, double* Fi, double* Fj, double* Ti, double* Tj,
              double& lam_next, PairPhysics& ph) {
    const std::size_t n = NN ? NN : n_runtime;
    LamEval e;
    double lam;
    if (!lam_iterate(Ai, Aj, r, n, n_iter, true, 0.0, lam_start, track_tol, lam, lam_next, e)) return 0;
    const double s2 = lam * (1.0 - lam) * e.rx;
    if (!(s2 < 1.0)) return 1;
    ph = pair_physics(r, d2, n, s2, lam, e, K, Fi, Fj, Ti, Tj);
    return 2;
}

using PairFn = int (*)(const double*, const double*, const double*, double, std::size_t, int, double, double,
                       double, double*, double*, double*, double*, double&, PairPhysics&);

PairFn pick_pair(std::size_t n) {
    switch (n) {
        case 2: return pair_full<2>;
        case 3: return pair_full<3>;
        case 4: return pair_full<4>;
        case 5: return pair_full<5>;
        case 6: return pair_full<6>;
        default: return pair_full<0>;
    }
}

}  // namespace

// ------------------------------------------------------------------------------------------------------------

void ell_check_ratios(const double* ratios, std::size_t count, std::size_t Ndim, const std::string& who) {
    if (Ndim < 2 || Ndim > ELL_MAX_DIM) {
        throw std::runtime_error("ellipsoids: Ndim must be between 2 and " + std::to_string(ELL_MAX_DIM));
    }
    if (count != Ndim - 1) {
        throw std::runtime_error("ellipsoids: aspect_ratios" + who + " needs Ndim-1 = " + std::to_string(Ndim - 1) +
                                 " entries (axis k has length D x ratio_k; the first axis is D), got " +
                                 std::to_string(count));
    }
    double prev = 1.0;
    for (std::size_t k = 0; k < count; ++k) {
        const double a = ratios[k];
        if (!(a > 0.0) || !std::isfinite(a)) {
            throw std::runtime_error("ellipsoids: every aspect ratio" + who + " must be a positive finite number");
        }
        if (a > prev) {
            throw std::runtime_error(
                "ellipsoids: aspect_ratios" + who + " must be in decreasing order, 1 >= ratio_2 >= ratio_3 >= ... > 0 "
                "(D is the LONGEST axis; e.g. 3D prolate [0.5, 0.5], oblate [1, 0.5]). Got " + std::to_string(a) +
                " after " + std::to_string(prev) + ". The input is not reordered for you, because that would "
                "rescale D and change what the size distribution means: divide D by the largest ratio and list "
                "the ratios longest first");
        }
        prev = a;
    }
}

void ell_build_shape_table(const std::vector<double>& single, const std::vector<std::vector<double>>& rows,
                           std::size_t N, std::size_t Ndim, std::vector<double>& qtab, std::vector<double>& deq,
                           bool& all_one) {
    if (!single.empty() && !rows.empty()) {
        throw std::runtime_error("ellipsoids: aspect ratios were given both as one list and as a per-particle "
                                 "table; supply one of them");
    }
    if (!rows.empty() && rows.size() != N) {
        throw std::runtime_error("ellipsoids: the per-particle aspect_ratios table needs one row per particle (" +
                                 std::to_string(N) + "), got " + std::to_string(rows.size()));
    }
    if (rows.empty()) ell_check_ratios(single.data(), single.size(), Ndim, "");
    qtab.assign(N * Ndim, 1.0);
    deq.assign(N, 1.0);
    all_one = true;
    for (std::size_t i = 0; i < N; ++i) {
        const std::vector<double>& r = rows.empty() ? single : rows[i];
        if (!rows.empty()) ell_check_ratios(r.data(), r.size(), Ndim, " (particle " + std::to_string(i) + ")");
        double vf = 1.0; bool one = true;
        for (std::size_t k = 1; k < Ndim; ++k) {
            qtab[i * Ndim + k] = r[k - 1];
            vf *= r[k - 1];
            if (r[k - 1] != 1.0) one = false;
        }
        deq[i] = one ? 1.0 : std::pow(vf, 1.0 / static_cast<double>(Ndim));
        if (!one) all_one = false;
    }
}

double ell_volume_factor(std::size_t N, std::size_t Ndim, const double* D, const double* qtab) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double v = 1.0, vf = 1.0;
        for (std::size_t k = 0; k < Ndim; ++k) { v *= D[i]; vf *= qtab[i * Ndim + k]; }
        num += v * vf; den += v;
    }
    return den > 0.0 ? num / den : 1.0;
}

EllShape make_ell_shape(const std::vector<double>& aspect_ratios, std::size_t Ndim, int n_iter) {
    if (Ndim < 2 || Ndim > ELL_MAX_DIM) {
        throw std::runtime_error("ellipsoids: Ndim must be between 2 and " + std::to_string(ELL_MAX_DIM));
    }
    if (aspect_ratios.size() != Ndim - 1) {
        throw std::runtime_error("ellipsoids: aspect_ratios needs Ndim-1 = " + std::to_string(Ndim - 1) +
                                 " entries (axis k has length D x ratio_k; the first axis is D), got " +
                                 std::to_string(aspect_ratios.size()));
    }
    ell_check_ratios(aspect_ratios.data(), aspect_ratios.size(), Ndim, "");
    EllShape sh;
    sh.Ndim = Ndim;
    sh.n_iter = n_iter < 1 ? 1 : n_iter;
    sh.rmax = 1.0;             // with the order rule the first axis is the longest: Db == D
    sh.all_one = true;
    for (double a : aspect_ratios) if (a != 1.0) sh.all_one = false;
    sh.q[0] = 1.0 / sh.rmax;
    for (std::size_t k = 1; k < Ndim; ++k) sh.q[k] = aspect_ratios[k - 1] / sh.rmax;
    sh.vol_factor = 1.0;
    for (std::size_t k = 0; k < Ndim; ++k) sh.vol_factor *= sh.q[k];
    sh.deq_factor = sh.all_one ? 1.0 : std::pow(sh.vol_factor, 1.0 / static_cast<double>(Ndim));
    return sh;
}

void ell_check_orientation(const double* R, std::size_t Ndim, std::size_t index) {
    for (std::size_t a = 0; a < Ndim; ++a) {
        for (std::size_t b = a; b < Ndim; ++b) {
            double sum = 0.0;
            for (std::size_t k = 0; k < Ndim; ++k) sum += R[k * Ndim + a] * R[k * Ndim + b];
            const double want = (a == b) ? 1.0 : 0.0;
            if (!(std::abs(sum - want) < 1e-8)) {
                throw std::runtime_error("ellipsoids: orientation of particle " + std::to_string(index) +
                                         " is not a rotation matrix (columns must be orthonormal)");
            }
        }
    }
}

namespace {

// Modified Gram-Schmidt on the columns of the row-major n x n matrix R.
ELL_INLINE void orthonormalise_columns(double* R, std::size_t n) {
    for (std::size_t c = 0; c < n; ++c) {
        for (std::size_t p = 0; p < c; ++p) {
            double proj = 0.0;
            for (std::size_t k = 0; k < n; ++k) proj += R[k * n + c] * R[k * n + p];
            for (std::size_t k = 0; k < n; ++k) R[k * n + c] -= proj * R[k * n + p];
        }
        double norm = 0.0;
        for (std::size_t k = 0; k < n; ++k) norm += R[k * n + c] * R[k * n + c];
        norm = std::sqrt(norm);
        for (std::size_t k = 0; k < n; ++k) R[k * n + c] /= norm;
    }
}

inline double determinant(const double* R, std::size_t n) {
    double M[MD][MD];
    for (std::size_t a = 0; a < n; ++a) for (std::size_t b = 0; b < n; ++b) M[a][b] = R[a * n + b];
    double det = 1.0;
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t a = c + 1; a < n; ++a) if (std::abs(M[a][c]) > std::abs(M[piv][c])) piv = a;
        if (M[piv][c] == 0.0) return 0.0;
        if (piv != c) { for (std::size_t b = 0; b < n; ++b) std::swap(M[piv][b], M[c][b]); det = -det; }
        det *= M[c][c];
        for (std::size_t a = c + 1; a < n; ++a) {
            const double f = M[a][c] / M[c][c];
            for (std::size_t b = c; b < n; ++b) M[a][b] -= f * M[c][b];
        }
    }
    return det;
}

}  // namespace

void ell_random_orientations(std::size_t N, std::size_t Ndim, std::uint64_t seed, std::vector<double>& out) {
    out.assign(N * Ndim * Ndim, 0.0);
    std::mt19937_64 rng(seed ^ 0xE11150D5EEDULL);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (std::size_t i = 0; i < N; ++i) {
        double* R = out.data() + i * Ndim * Ndim;
        for (std::size_t k = 0; k < Ndim * Ndim; ++k) R[k] = gauss(rng);
        orthonormalise_columns(R, Ndim);                       // Gaussian matrix -> Haar-distributed frame
        if (determinant(R, Ndim) < 0.0) {
            for (std::size_t k = 0; k < Ndim; ++k) R[k * Ndim] = -R[k * Ndim];   // make it a proper rotation
        }
    }
}

bool ell_pair_contact(const double* Ai, const double* Aj, const double* r, std::size_t Ndim, int n_iter,
                      double& s2, double& lam, double* xvec) {
    return pick_contact(Ndim)(Ai, Aj, r, Ndim, n_iter, s2, lam, xvec);
}

double ell_pair_s2_exact(const double* Ai, const double* Aj, const double* r, std::size_t Ndim) {
    LamEval e;
    double lam = 0.5, lam_next = 0.5;
    lam_iterate(Ai, Aj, r, Ndim, 200, false, 1e-16, -1.0, 0.0, lam, lam_next, e);
    return lam * (1.0 - lam) * e.rx;
}

void ell_shape_matrices(std::size_t N, const EllShape& sh, const double* Db, const double* Rot, double* Amat) {
    const std::size_t n = sh.Ndim, nn = n * n;
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        const double* R = Rot + i * nn;
        double* A = Amat + i * nn;
        double a2[MD];
        const double half = 0.5 * Db[i];
        const double* q = sh.qtab ? sh.qtab + i * n : sh.q;
        for (std::size_t k = 0; k < n; ++k) { const double ak = half * q[k]; a2[k] = ak * ak; }
        for (std::size_t a = 0; a < n; ++a) {
            for (std::size_t b = a; b < n; ++b) {
                double sum = 0.0;
                for (std::size_t k = 0; k < n; ++k) sum += R[a * n + k] * a2[k] * R[b * n + k];
                A[a * n + b] = sum;
                A[b * n + a] = sum;
            }
        }
    }
}

namespace {

struct LoopCtx {
    const std::uint32_t* pairs_data;
    const std::uint64_t* pair_offsets;
    const double* x;
    const double* Dd;
    const double* Amat;
    const double* boxd;
    const double* inv_box;
    const double* half_box;
    std::size_t n_runtime;
    int n_iter;
    double track_tol;
    double conv_tol;
    double K;
    std::uint32_t* c_ids;
    double* c_lam;
    std::uint32_t* c_cnt;
    bool want_gsum;
};

struct LoopAcc {
    double U = 0.0, Lc = 0.0, dkappa = 0.0, Fmean = 0.0, E = 0.0, max_min = 0.0;
    std::size_t count = 0;
};

// Rows [i_begin, i_end) of the half neighbour list. Compiled once per dimension (NN == 0: generic), so the
// scan over list entries, the contact evaluation and the physics are all unrolled for that dimension.
template <std::size_t NN>
void process_rows(const LoopCtx& c, std::size_t i_begin, std::size_t i_end, double* F_t, double* T_t,
                  std::size_t* z_t, double* g_t, LoopAcc& acc) {
    const std::size_t n = NN ? NN : c.n_runtime, nn = n * n, nrot = n * (n - 1) / 2;
    constexpr std::size_t CAP = EllPairCache::CAP;
    double half_box[MD], boxl[MD];
    for (std::size_t d = 0; d < n; ++d) { half_box[d] = c.half_box[d]; boxl[d] = c.boxd[d]; }
    const double* x = c.x;
    const double* Dd = c.Dd;
    const double* Amat = c.Amat;
    std::vector<std::uint32_t> keep(256);
    std::uint32_t keep_cap = 255;
    for (std::size_t i = i_begin; i < i_end; ++i) {
        const std::uint32_t* pi = c.pairs_data + c.pair_offsets[i];
        const std::uint32_t numNbr = pi[0];
        double xi[MD];
        for (std::size_t d = 0; d < n; ++d) xi[d] = x[i * n + d];
        const double Di = Dd[i];
        const double* Ai = Amat + i * nn;
        // this row's lam cache from the previous pass (ascending j, like the row itself), and the new one
        const std::uint32_t* old_ids = c.c_ids ? c.c_ids + i * CAP : nullptr;
        const double* old_lam = c.c_lam ? c.c_lam + i * CAP : nullptr;
        const std::uint32_t old_n = c.c_cnt ? c.c_cnt[i] : 0;
        std::uint32_t po = 0, new_n = 0;
        std::uint32_t new_ids[CAP];
        double new_lam[CAP];
        // Pass 1 over the row, branch-free: which entries have overlapping bounding spheres. Positions are
        // kept inside the box, so the minimum image is at most one box length (same value as the rounding
        // form the sphere kernel uses); a wall axis has half_box = inf and is never wrapped.
        std::uint32_t n_keep = 0;
        if (numNbr > keep_cap) { keep_cap = numNbr + 64; keep.resize(keep_cap); }
        std::uint32_t* keep_p = keep.data();
#if defined(__SSE2__)
        if (n == 3) {
            // The same operations as the scalar form below, x and y in one SSE2 register (always available
            // on x86-64), in the same order, so d2 is bit-identical.
            const __m128d xi01 = _mm_loadu_pd(xi), h01 = _mm_loadu_pd(half_box), b01 = _mm_loadu_pd(boxl);
            const __m128d nh01 = _mm_sub_pd(_mm_setzero_pd(), h01);
            const __m128d xi2 = _mm_load_sd(xi + 2), h2 = _mm_load_sd(half_box + 2), b2 = _mm_load_sd(boxl + 2);
            const __m128d nh2 = _mm_sub_sd(_mm_setzero_pd(), h2);
            const __m128d halfDi = _mm_set_sd(0.5);
            for (std::uint32_t jdx = 1; jdx <= numNbr; ++jdx) {
                const std::uint32_t j = pi[jdx];
                const double* xj = x + static_cast<std::size_t>(j) * 3;
                __m128d d01 = _mm_sub_pd(_mm_loadu_pd(xj), xi01);
                d01 = _mm_sub_pd(d01, _mm_and_pd(_mm_cmpgt_pd(d01, h01), b01));
                d01 = _mm_add_pd(d01, _mm_and_pd(_mm_cmplt_pd(d01, nh01), b01));
                __m128d dz = _mm_sub_sd(_mm_load_sd(xj + 2), xi2);
                dz = _mm_sub_sd(dz, _mm_and_pd(_mm_cmpgt_sd(dz, h2), b2));
                dz = _mm_add_sd(dz, _mm_and_pd(_mm_cmplt_sd(dz, nh2), b2));
                const __m128d sq = _mm_mul_pd(d01, d01);
                __m128d d2v = _mm_add_sd(_mm_add_sd(sq, _mm_unpackhi_pd(sq, sq)), _mm_mul_sd(dz, dz));
                __m128d rbv = _mm_mul_sd(halfDi, _mm_add_sd(_mm_set_sd(Di), _mm_load_sd(Dd + j)));
                keep_p[n_keep] = j;
                n_keep += static_cast<std::uint32_t>(_mm_comilt_sd(d2v, _mm_mul_sd(rbv, rbv)));
            }
        } else
#endif
        for (std::uint32_t jdx = 1; jdx <= numNbr; ++jdx) {
            const std::size_t j = pi[jdx];
            const double rb = 0.5 * (Di + Dd[j]);          // bounding spheres
            const double* xj = x + j * n;
            double d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = xj[d] - xi[d];
                delta -= (delta > half_box[d] ? boxl[d] : 0.0);
                delta += (delta < -half_box[d] ? boxl[d] : 0.0);
                d2 += delta * delta;
            }
            keep_p[n_keep] = static_cast<std::uint32_t>(j);
            n_keep += (d2 < rb * rb) ? 1u : 0u;
        }
        // Pass 2: only the pairs whose bounding spheres overlap.
        for (std::uint32_t kk = 0; kk < n_keep; ++kk) {
            const std::size_t j = keep_p[kk];
            const double* xj = x + j * n;
            double r[MD];
            double d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = xj[d] - xi[d];
                delta -= (delta > half_box[d] ? boxl[d] : 0.0);
                delta += (delta < -half_box[d] ? boxl[d] : 0.0);
                r[d] = delta;
                d2 += delta * delta;
            }
            // The pair's cache entry: +lam if it was in contact at the last pass, -lam if it was apart.
            double cached = 0.0;
            while (po < old_n && old_ids[po] < j) ++po;
            if (po < old_n && old_ids[po] == j) cached = old_lam[po];
            const double* Aj = Amat + j * nn;
            if (!(cached > 0.0)) {
                // Not in contact at the last pass: try the separating axis along r first. The projections
                // of the two ellipsoids on r/|r| have half-widths sqrt(r.A r)/|r|; if they do not overlap
                // the bodies are apart (exact; the converse does not hold). No solve, no lam.
                double qa = 0.0, qb = 0.0;
                for (std::size_t a = 0; a < n; ++a) {
                    double sa = 0.0, sb = 0.0;
                    for (std::size_t b = 0; b < n; ++b) { sa += Ai[a * n + b] * r[b]; sb += Aj[a * n + b] * r[b]; }
                    qa += r[a] * sa;
                    qb += r[a] * sb;
                }
                if (d2 >= std::sqrt(qa) + std::sqrt(qb)) {
                    if (c.c_ids && new_n < CAP) { new_ids[new_n] = static_cast<std::uint32_t>(j); new_lam[new_n] = cached < 0.0 ? cached : -2.0; ++new_n; }
                    continue;
                }
            }
            LamEval e;
            double lam, lam_next = -1.0;
            const bool touching = lam_iterate(Ai, Aj, r, n, c.n_iter, true, c.conv_tol, std::abs(cached), c.track_tol, lam, lam_next, e);
            const double s2 = touching ? lam * (1.0 - lam) * e.rx : 2.0;
            const bool contact = s2 < 1.0;
            if (c.c_ids && new_n < CAP) { new_ids[new_n] = static_cast<std::uint32_t>(j); new_lam[new_n] = contact ? lam_next : -lam_next; ++new_n; }
            if (!contact) continue;
            const PairPhysics ph = pair_physics(r, d2, n, s2, lam, e, c.K, F_t + i * n, F_t + j * n, T_t + i * nrot, T_t + j * nrot);
            z_t[i]++;
            z_t[j]++;
            if (c.want_gsum) { g_t[i] += ph.f_norm; g_t[j] += ph.f_norm; }
            acc.U += ph.overlap;
            ++acc.count;
            acc.Lc += ph.sigma;
            if (ph.overlap > acc.max_min) acc.max_min = ph.overlap;
            acc.dkappa += ph.dkappa;
            acc.Fmean += ph.F_mag;
            acc.E += ph.energy;
        }
        if (c.c_ids) {
            for (std::uint32_t k = 0; k < new_n; ++k) { c.c_ids[i * CAP + k] = new_ids[k]; c.c_lam[i * CAP + k] = new_lam[k]; }
            c.c_cnt[i] = new_n;
        }
    }
}

}  // namespace

void ell_forces(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                std::size_t N, const EllShape& sh, const std::vector<double>& Db,
                const std::vector<double>& box, const std::vector<std::int8_t>& walls, const double* Rot,
                double* F, double* Tq, double& U, double& max_min_dist, double& Lc, double& Fmean, double mu,
                double& dkappa, std::vector<std::size_t>& z, double& energy, std::vector<double>& Amat_work,
                EllPairCache* cache, double* gsum) {
    const double K = 1.0;
    const std::size_t n = sh.Ndim, nn = n * n, nrot = ell_nrot(n);
    const std::size_t F_size = N * n, T_size = N * nrot;
    Amat_work.resize(N * nn);
    double* Amat = Amat_work.data();
    ell_shape_matrices(N, sh, Db.data(), Rot, Amat);

    for (std::size_t k = 0; k < F_size; ++k) F[k] = 0.0;
    for (std::size_t k = 0; k < T_size; ++k) Tq[k] = 0.0;
    if (z.size() == N) std::fill(z.begin(), z.end(), 0); else z.assign(N, 0);

    double inv_box[MD] = {0.0}, half_box[MD] = {0.0};
    for (std::size_t d = 0; d < n; ++d) {
        inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;
        half_box[d] = walls[d] == 0 ? 0.5 * box[d] : std::numeric_limits<double>::infinity();
    }

    const int nthreads = omp_get_max_threads();
    // Per-thread accumulation buffers (thread 0 writes straight into the outputs); summed in thread order.
    std::vector<double> Fbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * F_size : 0, 0.0);
    std::vector<double> Tbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * T_size : 0, 0.0);
    std::vector<std::size_t> zbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * N : 0, 0);
    std::vector<double> gbuf((gsum && nthreads > 1) ? static_cast<std::size_t>(nthreads - 1) * N : 0, 0.0);
    if (gsum) for (std::size_t k = 0; k < N; ++k) gsum[k] = 0.0;

    double U_red = 0.0, Lc_red = 0.0, dkappa_red = 0.0, Fmean_red = 0.0, E_red = 0.0, max_min_red = 0.0;
    std::size_t count_red = 0;
    const double* Dd = Db.data();
    const double* boxd = box.data();
    const int n_iter = sh.n_iter;
    LoopCtx ctx;
    ctx.pairs_data = pairs_data; ctx.pair_offsets = pair_offsets; ctx.x = x; ctx.Dd = Dd; ctx.Amat = Amat;
    ctx.boxd = boxd; ctx.inv_box = inv_box; ctx.half_box = half_box; ctx.n_runtime = n; ctx.n_iter = n_iter;
    ctx.track_tol = sh.track_tol; ctx.conv_tol = sh.conv_tol; ctx.K = K; ctx.want_gsum = gsum != nullptr;
    ctx.c_ids = cache ? cache->ids.data() : nullptr;
    ctx.c_lam = cache ? cache->lam.data() : nullptr;
    ctx.c_cnt = cache ? cache->count.data() : nullptr;

    #pragma omp parallel reduction(+:U_red, Lc_red, dkappa_red, Fmean_red, E_red, count_red) \
                         reduction(max:max_min_red)
    {
        const int tid = omp_get_thread_num();
        const int nth = omp_get_num_threads();
        double* F_t = tid == 0 ? F : Fbuf.data() + static_cast<std::size_t>(tid - 1) * F_size;
        double* T_t = tid == 0 ? Tq : Tbuf.data() + static_cast<std::size_t>(tid - 1) * T_size;
        std::size_t* z_t = tid == 0 ? z.data() : zbuf.data() + static_cast<std::size_t>(tid - 1) * N;
        double* g_t = !gsum ? nullptr : (tid == 0 ? gsum : gbuf.data() + static_cast<std::size_t>(tid - 1) * N);
        const std::size_t i_begin = N * static_cast<std::size_t>(tid) / static_cast<std::size_t>(nth);
        const std::size_t i_end = N * static_cast<std::size_t>(tid + 1) / static_cast<std::size_t>(nth);
        LoopAcc acc;
        switch (n) {
            case 2: process_rows<2>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
            case 3: process_rows<3>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
            case 4: process_rows<4>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
            case 5: process_rows<5>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
            case 6: process_rows<6>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
            default: process_rows<0>(ctx, i_begin, i_end, F_t, T_t, z_t, g_t, acc); break;
        }
        U_red += acc.U;
        Lc_red += acc.Lc;
        dkappa_red += acc.dkappa;
        Fmean_red += acc.Fmean;
        E_red += acc.E;
        count_red += acc.count;
        if (acc.max_min > max_min_red) max_min_red = acc.max_min;
    }
    for (int tt = 1; tt < nthreads; ++tt) {
        const double* F_t = Fbuf.data() + static_cast<std::size_t>(tt - 1) * F_size;
        const double* T_t = Tbuf.data() + static_cast<std::size_t>(tt - 1) * T_size;
        const std::size_t* z_t = zbuf.data() + static_cast<std::size_t>(tt - 1) * N;
        for (std::size_t k = 0; k < F_size; ++k) F[k] += F_t[k];
        for (std::size_t k = 0; k < T_size; ++k) Tq[k] += T_t[k];
        for (std::size_t k = 0; k < N; ++k) z[k] += z_t[k];
        if (gsum) { const double* g_o = gbuf.data() + static_cast<std::size_t>(tt - 1) * N; for (std::size_t k = 0; k < N; ++k) gsum[k] += g_o[k]; }
    }

    U = U_red;
    Lc = Lc_red;
    dkappa = dkappa_red;
    Fmean = Fmean_red;
    energy = E_red;
    std::size_t count = 1 + count_red;
    max_min_dist = max_min_red;

    // Flat hard walls at 0 and box[d]: the half-width of the particle along axis d is h = sqrt(A_dd) (its
    // support function), in the role of the sphere's D/2. Same spring as the sphere code (twice a pair's).
    for (std::size_t d = 0; d < n; ++d) {
        if (walls[d] != 1) continue;
        for (std::size_t i = 0; i < N; ++i) {
            const double* A = Amat + i * nn;
            const double h = std::sqrt(A[d * n + d]);
            for (std::size_t wall = 0; wall < 2; ++wall) {
                const double dxw = box[d] * static_cast<double>(wall) - x[i * n + d];
                const double dist = std::abs(dxw);
                if (dist < h) {
                    const double F_mag = -2.0 * K * (h / dist - 1.0);
                    U += 2.0 * (1.0 - dist / h);
                    ++count;
                    F[i * n + d] += F_mag * dxw;
                    dkappa += 2.0 * K * h * (dist - h);
                    Fmean += F_mag;
                    energy += K * (h - dist) * (h - dist);
                    // E = K (h - dist)^2, dE/dA = K (h - dist)/h * e e^T = c e e^T;
                    // T_ab = -2 c (e_a (Ae)_b - (Ae)_a e_b): non-zero only when a == d or b == d.
                    const double cw = K * (h - dist) / h;
                    std::size_t c = 0;
                    for (std::size_t a = 0; a < n; ++a) {
                        for (std::size_t b = a + 1; b < n; ++b, ++c) {
                            const double ea = (a == d) ? 1.0 : 0.0, eb = (b == d) ? 1.0 : 0.0;
                            Tq[i * nrot + c] += -2.0 * cw * (ea * A[b * n + d] - A[a * n + d] * eb);
                        }
                    }
                }
            }
        }
    }

    // Curved container (the first -walls[0] axes): the particle's farthest point from the axis in the role of the
    // sphere's d + r; the sphere code's formulas with r -> the extent along the radial direction there
    // (curved_wall.hpp). Reached only with walls[0] < 0.
    if (walls[0] < 0) {
        const std::size_t t = curved_axes(walls[0], n);
        const double Rc = 0.5 * box[0];
        for (std::size_t i = 0; i < N; ++i) {
            CurvedWallSum cw;
            if (!shape_curved_wall(sh, n, t, Rc, K, x + i * n, Dd[i], Amat + i * nn, nullptr, nullptr, F + i * n, Tq + i * nrot, cw)) continue;
            U += cw.U; count += cw.count; dkappa += cw.dkappa; Fmean += cw.Fmean; energy += cw.energy;
        }
    }

    U = std::pow(U / double(count), 2);
    Lc = Lc / double(count);
    // The pressure-like term acts on the volume-equivalent diameter (factor 1 for spheres).
    if (sh.deqtab) {
        double sumDeq = 0.0;
        for (std::size_t i = 0; i < N; ++i) sumDeq += Db[i] * sh.deqtab[i];
        dkappa = dkappa + mu * sumDeq;
    } else {
        double sumD = 0.0;
        for (double Di : Db) sumD += Di;
        dkappa = dkappa + mu * (sumD * sh.deq_factor);
    }
    Fmean = -Fmean / count;
}

namespace {

// One particle's orientation update. Force-inlined into the per-dimension loops below.
// R <- (I - W/2)^-1 (I + W/2) R for a skew matrix W: an exact rotation (angle ~ |W| when small).
ELL_INLINE void apply_cayley(std::size_t n, double W[MD][MD], double* R, bool clean) {
    // Cayley map: R <- (I - W/2)^-1 (I + W/2) R, an exact rotation for skew W (angle ~ |W| when small).
    if (n == 2) {
        const double h = 0.5 * W[0][1];                      // (I - W/2)^-1 (I + W/2) = [[c, s], [-s, c]]
        const double inv = 1.0 / (1.0 + h * h);
        const double cs = (1.0 - h * h) * inv, sn = 2.0 * h * inv;
        const double r00 = R[0], r01 = R[1], r10 = R[2], r11 = R[3];
        R[0] = cs * r00 + sn * r10;
        R[1] = cs * r01 + sn * r11;
        R[2] = -sn * r00 + cs * r10;
        R[3] = -sn * r01 + cs * r11;
        const double nrm = 1.0 / std::sqrt(R[0] * R[0] + R[2] * R[2]);   // keep it a rotation against rounding drift
        R[0] *= nrm; R[2] *= nrm;
        R[1] = -R[2]; R[3] = R[0];
        return;
    }
    if (n == 3) {
        // (I - W/2)^-1 (I + W/2) = [(1 - |a|^2) I + 2 a a^T - 2 [a]x] / (1 + |a|^2) with W/2 = -[a]x, i.e.
        // a = -(W12, -W02, W01)/2 ... written out through the skew matrix H = W/2 directly:
        // Q = (I - H)^-1 (I + H) = I + 2 (H + H^2) / (1 + |h|^2),  |h|^2 = h01^2 + h02^2 + h12^2.
        const double h01 = 0.5 * W[0][1], h02 = 0.5 * W[0][2], h12 = 0.5 * W[1][2];
        const double f = 2.0 / (1.0 + h01 * h01 + h02 * h02 + h12 * h12);
        // H = [[0, h01, h02], [-h01, 0, h12], [-h02, -h12, 0]];  H^2 is symmetric.
        const double s00 = -(h01 * h01 + h02 * h02), s11 = -(h01 * h01 + h12 * h12), s22 = -(h02 * h02 + h12 * h12);
        const double s01 = -h02 * h12, s02 = h01 * h12, s12 = -h01 * h02;
        const double q00 = 1.0 + f * s00, q01 = f * (h01 + s01), q02 = f * (h02 + s02);
        const double q10 = f * (-h01 + s01), q11 = 1.0 + f * s11, q12 = f * (h12 + s12);
        const double q20 = f * (-h02 + s02), q21 = f * (-h12 + s12), q22 = 1.0 + f * s22;
        for (std::size_t b = 0; b < 3; ++b) {
            const double r0 = R[b], r1 = R[3 + b], r2 = R[6 + b];
            R[b] = q00 * r0 + q01 * r1 + q02 * r2;
            R[3 + b] = q10 * r0 + q11 * r1 + q12 * r2;
            R[6 + b] = q20 * r0 + q21 * r1 + q22 * r2;
        }
        if (clean) orthonormalise_columns(R, 3);           // remove rounding drift now and then
        return;
    }
    double M[MD][MD], Bm[MD][MD];
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = 0; b < n; ++b) {
            M[a][b] = (a == b ? 1.0 : 0.0) - 0.5 * W[a][b];
            double sum = R[a * n + b];
            for (std::size_t k = 0; k < n; ++k) sum += 0.5 * W[a][k] * R[k * n + b];
            Bm[a][b] = sum;
        }
    }
    // Gaussian elimination without pivoting: I - W/2 has unit diagonal and a skew off-diagonal part, so every
    // pivot is >= 1.
    for (std::size_t col = 0; col < n; ++col) {
        const double inv = 1.0 / M[col][col];
        for (std::size_t a = col + 1; a < n; ++a) {
            const double f = M[a][col] * inv;
            for (std::size_t b = col; b < n; ++b) M[a][b] -= f * M[col][b];
            for (std::size_t b = 0; b < n; ++b) Bm[a][b] -= f * Bm[col][b];
        }
    }
    for (std::size_t ai = n; ai-- > 0;) {
        const double inv = 1.0 / M[ai][ai];
        for (std::size_t b = 0; b < n; ++b) {
            double sum = Bm[ai][b];
            for (std::size_t k = ai + 1; k < n; ++k) sum -= M[ai][k] * R[k * n + b];
            R[ai * n + b] = sum * inv;
        }
    }
    orthonormalise_columns(R, n);                          // remove rounding drift
}

ELL_INLINE void rotate_one(std::size_t n, std::size_t nrot, const double* Tq, double* m, double* v, double* R,
                           double alpha_rot, double beta1, double beta2, double inv_b1_corr, double inv_b2_corr,
                           double epsilon, bool clean) {
    double W[MD][MD];
    bool any = false;
    std::size_t c = 0;
    for (std::size_t a = 0; a < n; ++a) {
        W[a][a] = 0.0;
        for (std::size_t b = a + 1; b < n; ++b, ++c) {
            const double Tk = Tq[c];
            const double mk = beta1 * m[c] - (1.0 - beta1) * Tk;
            const double vk = beta2 * v[c] + (1.0 - beta2) * Tk * Tk;
            m[c] = mk;
            v[c] = vk;
            const double w = -alpha_rot * (mk * inv_b1_corr) / (std::sqrt(vk * inv_b2_corr) + epsilon);
            W[a][b] = w;
            W[b][a] = -w;
            if (w != 0.0) any = true;
        }
    }
    (void)nrot;
    if (!any) return;   // spheres: zero torque -> the orientation is never touched
    apply_cayley(n, W, R, clean);
}

template <std::size_t NN>
void rotation_loop(std::size_t N, std::size_t n_runtime, const double* Tq, double* m, double* v, double* Rot,
                   const double* Db, double alpha, double beta1, double beta2, double inv_b1_corr,
                   double inv_b2_corr, double epsilon, bool clean) {
    const std::size_t n = NN ? NN : n_runtime, nn = n * n, nrot = n * (n - 1) / 2;
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
        rotate_one(n, nrot, Tq + i * nrot, m + i * nrot, v + i * nrot, Rot + i * nn, alpha / (0.5 * Db[i]),
                   beta1, beta2, inv_b1_corr, inv_b2_corr, epsilon, clean);
    }
}

}  // namespace

void ell_rotate_by(std::size_t Ndim, double* R, const double* w) {
    double W[MD][MD];
    bool any = false;
    std::size_t c = 0;
    for (std::size_t a = 0; a < Ndim; ++a) {
        W[a][a] = 0.0;
        for (std::size_t b = a + 1; b < Ndim; ++b, ++c) { W[a][b] = w[c]; W[b][a] = -w[c]; if (w[c] != 0.0) any = true; }
    }
    if (any) apply_cayley(Ndim, W, R, true);
}

void ell_rotation_step(std::size_t N, std::size_t Ndim, const double* Tq, double* m, double* v, double* Rot,
                       const double* Db, double alpha, double beta1, double beta2, double t, double epsilon) {
    const double inv_b1_corr = 1.0 / (1.0 - std::pow(beta1, t));
    const double inv_b2_corr = 1.0 / (1.0 - std::pow(beta2, t));
    const bool clean = (static_cast<long long>(t) % 16) == 0;   // 3D fast path: Gram-Schmidt every 16th step
    switch (Ndim) {
        case 2: rotation_loop<2>(N, Ndim, Tq, m, v, Rot, Db, alpha, beta1, beta2, inv_b1_corr, inv_b2_corr, epsilon, clean); break;
        case 3: rotation_loop<3>(N, Ndim, Tq, m, v, Rot, Db, alpha, beta1, beta2, inv_b1_corr, inv_b2_corr, epsilon, clean); break;
        case 4: rotation_loop<4>(N, Ndim, Tq, m, v, Rot, Db, alpha, beta1, beta2, inv_b1_corr, inv_b2_corr, epsilon, clean); break;
        default: rotation_loop<0>(N, Ndim, Tq, m, v, Rot, Db, alpha, beta1, beta2, inv_b1_corr, inv_b2_corr, epsilon, clean); break;
    }
}

// ---- diameter permutation with ellipsoids ------------------------------------------------------------------------

void ell_exchange(std::size_t Ndim, std::size_t i, std::size_t j, double* D, double* D0, double* qtab, double* deq,
                  double* Rot, bool move_orientations) {
    std::swap(D[i], D[j]);
    if (D0) std::swap(D0[i], D0[j]);
    for (std::size_t k = 0; k < Ndim; ++k) std::swap(qtab[i * Ndim + k], qtab[j * Ndim + k]);
    if (deq) std::swap(deq[i], deq[j]);
    if (move_orientations && Rot) {
        const std::size_t nn = Ndim * Ndim;
        for (std::size_t k = 0; k < nn; ++k) std::swap(Rot[i * nn + k], Rot[j * nn + k]);
    }
}

bool ell_exchange_admissible(std::size_t Ndim, double Di, const double* qi, double Dj, const double* qj, double w) {
    for (std::size_t k = 0; k < Ndim; ++k) {
        const double a = Di * qi[k], b = Dj * qj[k];
        if (!(a > 0.0) || !(b > 0.0)) return false;
        if ((a > b ? a / b : b / a) - 1.0 > w) return false;
    }
    return true;
}

std::size_t ell_permute_pass(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                             std::size_t N, std::size_t Ndim, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, double* D, double* D0, double* qtab,
                             double* deq, double* Rot, const double* Amat, const EllPermuteParams& prm,
                             std::vector<EllExchangeRecord>& done) {
    const std::size_t n = Ndim, nn = n * n, fw = n + 1;   // a facet: unit normal (n numbers) and its distance b
    if (N < 2) return 0;
    double inv_box[MD] = {0.0};
    for (std::size_t d = 0; d < n; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;
    auto support2 = [n](const double* A, const double* v) {   // v.A v
        double s = 0.0;
        for (std::size_t a = 0; a < n; ++a) { double t = 0.0; for (std::size_t b = 0; b < n; ++b) t += A[a * n + b] * v[b]; s += v[a] * t; }
        return s;
    };
    // 1. the cage of every site, from the pairs of the neighbour lists inside the 1.5 shell
    const double shell = 1.5;
    std::vector<std::vector<double>> facets(N);
    for (std::size_t i = 0; i < N; ++i) {
        const std::uint32_t* row = pairs_data + pair_offsets[i];   // row[0] = count, then the neighbours j > i
        const std::size_t cnt = row[0];
        const double* Ai = Amat + i * nn;
        for (std::size_t q = 1; q <= cnt; ++q) {
            const std::size_t j = row[q];
            const double rb = shell * 0.5 * (D[i] + D[j]);
            double r[MD], rs[MD], d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = x[j * n + d] - x[i * n + d];
                delta -= std::nearbyint(delta * inv_box[d]) * box[d];
                r[d] = delta; rs[d] = delta / shell; d2 += delta * delta;
            }
            if (!(d2 < rb * rb) || !(d2 > 0.0)) continue;
            const double* Aj = Amat + j * nn;
            double nr = 0.0, nv[MD], hi, hj;
            if (prm.kind != 0) {   // rods: the normal of the pair's contact (or closest approach), same 1.5 shell
                double ua[MD], ub[MD];
                for (std::size_t d = 0; d < n; ++d) { ua[d] = Rot[i * nn + d * n]; ub[d] = Rot[j * nn + d * n]; }
                const double hla = D[i] * qtab[i * n], ra = D[i] * qtab[i * n + 1], hlb = D[j] * qtab[j * n], rbb = D[j] * qtab[j * n + 1];
                // uncapped: the maximising normal does not depend on the scale, so r / 1.5 gives the shell test and
                // the normal at once. Capped: the closest points of the two segments DO depend on the separation,
                // so the normal is taken at the true r and the shell is "segment distance < 1.5 (R_i + R_k)".
                if (prm.kind == SHAPE_ROD_CAPPED) { if (!(rod_pair_normal(n, prm.kind, r, ua, hla, ra, ub, hlb, rbb, 30, nv) < shell)) continue; }
                else if (!(rod_pair_normal(n, prm.kind, rs, ua, hla, ra, ub, hlb, rbb, 30, nv) < 1.0)) continue;
                for (std::size_t d = 0; d < n; ++d) nr += nv[d] * r[d];
                hi = rod_support(n, prm.kind, ua, hla, ra, nv);
                hj = rod_support(n, prm.kind, ub, hlb, rbb, nv);
            } else {
            double s2 = 0.0, lam = 0.5, xv[MD];
            if (!ell_pair_contact(Ai, Aj, rs, n, prm.n_iter, s2, lam, xv)) continue;   // apart even at 1.5 x: not a facet
            double xn = 0.0;
            for (std::size_t d = 0; d < n; ++d) xn += xv[d] * xv[d];
            if (!(xn > 0.0)) continue;
            xn = 1.0 / std::sqrt(xn);
            for (std::size_t d = 0; d < n; ++d) { nv[d] = xv[d] * xn; nr += nv[d] * r[d]; }
            hi = std::sqrt(support2(Ai, nv)); hj = std::sqrt(support2(Aj, nv));
            }
            std::vector<double>& fi = facets[i];
            std::vector<double>& fj = facets[j];
            for (std::size_t d = 0; d < n; ++d) { fi.push_back(nv[d]); fj.push_back(-nv[d]); }
            fi.push_back(nr - hj);   // room left at site i by neighbour j
            fj.push_back(nr - hi);
        }
    }
    // room of particle p (shape matrix A, or size and shape placed in the site's orientation) at site a
    const double INF = std::numeric_limits<double>::infinity();
    double room_overlap = 0.0;   // set by room(): sum over the planes of max(0, h - b)^2
    auto room = [&](std::size_t a, std::size_t p) {
        const std::vector<double>& f = facets[a];
        double best = INF, ov = 0.0;
        if (prm.kind != 0) {   // rods: particle p's length and radius, on its own axis or on the site's
            const std::size_t o = (!prm.decide_in_place || a == p) ? p : a;
            double ua[MD];
            for (std::size_t d = 0; d < n; ++d) ua[d] = Rot[o * nn + d * n];
            const double hlp = D[p] * qtab[p * n], rp = D[p] * qtab[p * n + 1];
            for (std::size_t e = 0; e + fw <= f.size(); e += fw) {
                const double h = rod_support(n, prm.kind, ua, hlp, rp, f.data() + e);
                const double t = f[e + n] / h;
                if (t < best) best = t;
                if (h > f[e + n]) ov += (h - f[e + n]) * (h - f[e + n]);
            }
            room_overlap = ov;
            return best;
        }
        if (!prm.decide_in_place || a == p) {
            const double* A = Amat + p * nn;
            for (std::size_t e = 0; e + fw <= f.size(); e += fw) {
                const double h = std::sqrt(support2(A, f.data() + e));
                const double t = f[e + n] / h;
                if (t < best) best = t;
                if (h > f[e + n]) ov += (h - f[e + n]) * (h - f[e + n]);
            }
        } else {
            const double* R = Rot + a * nn;   // the site keeps its orientation; p brings size and shape
            double a2[MD];
            for (std::size_t k = 0; k < n; ++k) { const double ak = 0.5 * D[p] * qtab[p * n + k]; a2[k] = ak * ak; }
            for (std::size_t e = 0; e + fw <= f.size(); e += fw) {
                double h2 = 0.0;
                for (std::size_t k = 0; k < n; ++k) {
                    double c = 0.0;
                    for (std::size_t d = 0; d < n; ++d) c += R[d * n + k] * f[e + d];
                    h2 += a2[k] * c * c;
                }
                const double h = std::sqrt(h2);
                const double t = f[e + n] / h;
                if (t < best) best = t;
                if (h > f[e + n]) ov += (h - f[e + n]) * (h - f[e + n]);
            }
        }
        room_overlap = ov;
        return best;
    };
    std::vector<double> self(N, INF), self_ov(N, 0.0);
    for (std::size_t p = 0; p < N; ++p) if (!facets[p].empty()) { self[p] = room(p, p); self_ov[p] = room_overlap; }
    // 2. candidates: partners in the window, in order of the longest axis D (ties by index: deterministic)
    std::vector<std::uint32_t> order(N);
    for (std::size_t k = 0; k < N; ++k) order[k] = static_cast<std::uint32_t>(k);
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return D[a] < D[b] || (D[a] == D[b] && a < b); });
    std::vector<EllExchangeRecord> cands;
    const long cap = std::max(1, prm.max_partners);
    std::size_t examined = 0, q1 = 0;
    for (std::size_t r = 0; r < N; ++r) {
        const std::uint32_t i = order[r];
        if (q1 < r) q1 = r;
        while (q1 + 1 < N && D[order[q1 + 1]] <= D[i] * (1.0 + prm.window)) ++q1;
        if (!(self[i] < INF) || q1 <= r) continue;
        const long span = static_cast<long>(q1 - r);
        const long stride = std::max<long>(1, (span + cap - 1) / cap);
        for (long q = static_cast<long>(r) + 1 + static_cast<long>(prm.pass_index % static_cast<std::uint64_t>(stride));
             q <= static_cast<long>(q1); q += stride) {
            const std::uint32_t j = order[q];
            if (!(self[j] < INF)) continue;
            if (prm.kind != 0) { if (!rod_exchange_admissible(D[i], qtab + i * n, D[j], qtab + j * n, prm.window)) continue; }
            else
            if (!ell_exchange_admissible(n, D[i], qtab + i * n, D[j], qtab + j * n, prm.window)) continue;
            if (!prm.move_orientations && D[i] == D[j] && std::equal(qtab + i * n, qtab + i * n + n, qtab + j * n)) continue;   // nothing would change
            ++examined;
            const double Ri2 = 0.25 * D[i] * D[i], Rj2 = 0.25 * D[j] * D[j];
            const double tji = room(j, i), oji = room_overlap, tij = room(i, j), oij = room_overlap;
            const double gain = 2.0 * (Ri2 * tji + Rj2 * tij - Ri2 * self[i] - Rj2 * self[j]);
            // the guard concerns orientation: for two spheres nothing can stick out differently after the exchange
            // than the room already measures, and the sphere engine's rule is kept as it is
            if (prm.overlap_guard && oji + oij > self_ov[i] + self_ov[j]) {
                bool spheres = true;
                for (std::size_t k = 1; k < n; ++k) if (qtab[i * n + k] != 1.0 || qtab[j * n + k] != 1.0) spheres = false;
                if (!spheres) continue;
            }
            if (gain > 1e-9 * (Ri2 + Rj2)) cands.push_back({i, j, gain});   // strictly better, clear of rounding
        }
    }
    std::sort(cands.begin(), cands.end(), [](const EllExchangeRecord& a, const EllExchangeRecord& b) {
        return a.gain > b.gain || (a.gain == b.gain && (a.i < b.i || (a.i == b.i && a.j < b.j)));
    });
    // 3. accept greedily: each particle once per pass, cages as they were at the start of the pass
    std::vector<unsigned char> used(N, 0);
    long acc = 0;
    for (const EllExchangeRecord& c : cands) {
        if (acc >= prm.budget) break;
        if (used[c.i] || used[c.j]) continue;
        used[c.i] = used[c.j] = 1;
        ell_exchange(n, c.i, c.j, D, D0, qtab, deq, Rot, prm.move_orientations);
        done.push_back(c);
        ++acc;
    }
    return examined;
}

double ell_min_contact_scale(std::size_t N, const EllShape& sh, const double* x, const double* Db,
                             const double* Rot, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, std::size_t& pair_i, std::size_t& pair_j,
                             std::size_t& n_overlapping) {
    const std::size_t n = sh.Ndim, nn = n * n;
    std::vector<double> Amat(N * nn);
    ell_shape_matrices(N, sh, Db, Rot, Amat.data());
    double inv_box[MD] = {0.0};
    for (std::size_t d = 0; d < n; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;

    double best = 1.0;
    std::size_t bi = N, bj = N, nover = 0;
    #pragma omp parallel
    {
        double tbest = 1.0;
        std::size_t tbi = N, tbj = N, tn = 0;
        #pragma omp for schedule(dynamic, 16)
        for (std::size_t i = 0; i < N; ++i) {
            const double* xi = x + i * n;
            const double* Ai = Amat.data() + i * nn;
            for (std::size_t j = i + 1; j < N; ++j) {
                const double rb = 0.5 * (Db[i] + Db[j]);
                double r[MD];
                double d2 = 0.0;
                for (std::size_t d = 0; d < n; ++d) {
                    double delta = x[j * n + d] - xi[d];
                    delta -= std::nearbyint(delta * inv_box[d]) * box[d];
                    r[d] = delta;
                    d2 += delta * delta;
                }
                if (!(d2 < rb * rb)) continue;
                const double s2 = ell_pair_s2_exact(Ai, Amat.data() + j * nn, r, n);
                if (s2 < 1.0) {
                    ++tn;
                    const double s = std::sqrt(s2);
                    if (s < tbest) { tbest = s; tbi = i; tbj = j; }
                }
            }
            for (std::size_t d = 0; d < n; ++d) {
                if (walls[d] != 1) continue;
                const double h = std::sqrt(Ai[d * n + d]);
                for (std::size_t wall = 0; wall < 2; ++wall) {
                    const double dist = std::abs(box[d] * static_cast<double>(wall) - xi[d]);
                    if (dist < h) {
                        ++tn;
                        const double s = dist / h;
                        if (s < tbest) { tbest = s; tbi = i; tbj = N; }
                    }
                }
            }
            if (walls[0] < 0) {   // curved container: the factor that brings the farthest point back onto the wall
                const std::size_t t = curved_axes(walls[0], n);
                const double Rc = 0.5 * box[0];
                double xc[MD];
                for (std::size_t k = 0; k < n; ++k) xc[k] = k < t ? xi[k] - Rc : 0.0;
                const double s = ell_wall_scale(n, t, Ai, xc, Rc, 1.0);
                if (s < 1.0) { ++tn; if (s < tbest) { tbest = s; tbi = i; tbj = N; } }
            }
        }
        #pragma omp critical
        {
            nover += tn;
            if (tbest < best || (tbest == best && tbi < bi)) { best = tbest; bi = tbi; bj = tbj; }
        }
    }
    pair_i = bi;
    pair_j = bj;
    n_overlapping = nover;
    return best;
}

void ell_pair_cost_probe(std::size_t Ndim, const std::vector<double>& aspect_ratios, int n_iter,
                         std::size_t n_pairs, std::size_t repeats, std::uint64_t seed,
                         double& t_ellipsoid, double& t_sphere, double& t_tracked, double& checksum) {
    const EllShape sh = make_ell_shape(aspect_ratios, Ndim, n_iter);
    const std::size_t n = Ndim, nn = n * n, nrot = ell_nrot(n);
    std::vector<double> Rot;
    ell_random_orientations(2 * n_pairs, n, seed, Rot);
    std::vector<double> Db(2 * n_pairs, 1.0), Amat(2 * n_pairs * nn), rvec(n_pairs * n), lam_conv(n_pairs, -1.0);
    ell_shape_matrices(2 * n_pairs, sh, Db.data(), Rot.data(), Amat.data());
    std::mt19937_64 rng(seed + 17);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::uniform_real_distribution<double> uni(0.80, 0.98);
    for (std::size_t p = 0; p < n_pairs; ++p) {            // overlapping pairs: s between 0.80 and 0.98
        double u[MD], norm = 0.0;
        for (std::size_t k = 0; k < n; ++k) { u[k] = gauss(rng); norm += u[k] * u[k]; }
        norm = std::sqrt(norm);
        for (std::size_t k = 0; k < n; ++k) u[k] /= norm;
        const double s_unit = std::sqrt(ell_pair_s2_exact(Amat.data() + (2 * p) * nn, Amat.data() + (2 * p + 1) * nn, u, n));
        const double dist = uni(rng) / s_unit;
        for (std::size_t k = 0; k < n; ++k) rvec[p * n + k] = u[k] * dist;
        LamEval e;
        double lam, lam_next;
        lam_iterate(Amat.data() + (2 * p) * nn, Amat.data() + (2 * p + 1) * nn, rvec.data() + p * n, n, 200, false,
                    1e-16, -1.0, 0.0, lam, lam_next, e);
        lam_conv[p] = lam;
    }
    std::vector<double> F(2 * n_pairs * n, 0.0), T(2 * n_pairs * nrot, 0.0);
    const PairFn pair = pick_pair(n);
    auto run = [&](bool tracked) {
        double acc = 0.0;
        for (std::size_t rep = 0; rep < repeats; ++rep) {
            for (std::size_t p = 0; p < n_pairs; ++p) {
                const double* r = rvec.data() + p * n;
                double d2 = 0.0;
                for (std::size_t k = 0; k < n; ++k) d2 += r[k] * r[k];
                double lam_next;
                PairPhysics ph;
                const int state = pair(Amat.data() + (2 * p) * nn, Amat.data() + (2 * p + 1) * nn, r, d2, n, n_iter,
                                       tracked ? lam_conv[p] : -1.0, tracked ? 1e-6 : 0.0, 1.0,
                                       F.data() + (2 * p) * n, F.data() + (2 * p + 1) * n,
                                       T.data() + (2 * p) * nrot, T.data() + (2 * p + 1) * nrot, lam_next, ph);
                if (state == 2) acc += ph.dkappa;
            }
        }
        return acc;
    };
    auto t0 = std::chrono::steady_clock::now();
    const double acc = run(false);
    auto t1 = std::chrono::steady_clock::now();
    // The sphere pair body of get_forces_nd_3's scalar kernel, on the same separations (contact distance set
    // so that every pair overlaps, as above).
    std::vector<double> Fs(2 * n_pairs * n, 0.0);
    double acc_s = 0.0;
    for (std::size_t rep = 0; rep < repeats; ++rep) {
        for (std::size_t p = 0; p < n_pairs; ++p) {
            const double* r = rvec.data() + p * n;
            double d2 = 0.0;
            for (std::size_t k = 0; k < n; ++k) d2 += r[k] * r[k];
            const double r_ij = 1.1 * std::sqrt(d2) + 1e-3 * static_cast<double>(rep % 2);
            if (d2 < r_ij * r_ij) {
                const double dist = std::sqrt(d2);
                const double F_mag = -(r_ij / dist - 1.0);
                for (std::size_t k = 0; k < n; ++k) {
                    const double fcomp = F_mag * r[k];
                    Fs[(2 * p) * n + k] += fcomp;
                    Fs[(2 * p + 1) * n + k] -= fcomp;
                }
                acc_s += r_ij * (dist - r_ij);
            }
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    const double acc_t = run(true);
    auto t3 = std::chrono::steady_clock::now();
    const double denom = static_cast<double>(n_pairs) * static_cast<double>(repeats);
    t_ellipsoid = std::chrono::duration<double>(t1 - t0).count() / denom;
    t_sphere = std::chrono::duration<double>(t2 - t1).count() / denom;
    t_tracked = std::chrono::duration<double>(t3 - t2).count() / denom;
    checksum = acc + acc_s + acc_t + F[0] + Fs[0] + T[0];
}

}  // namespace rcpgenerator
