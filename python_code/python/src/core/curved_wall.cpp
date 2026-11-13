// Ellipsoids and rods against the wall of a curved container. See rcpgenerator/curved_wall.hpp.
// Additive: nothing in this file is reached for spheres, or for shapes in periodic or flat-wall boxes.

#include "rcpgenerator/curved_wall.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "rcpgenerator/rod.hpp"

namespace rcpgenerator {

namespace {

constexpr std::size_t MD = ELL_MAX_DIM;
constexpr double PI = 3.14159265358979323846;

// Cyclic Jacobi for a symmetric t x t matrix S (row-major, destroyed): eigenvalues in w, eigenvectors in the
// COLUMNS of V. t <= ELL_MAX_DIM, so a fixed number of sweeps reaches rounding.
void jacobi(std::size_t t, double* S, double* w, double* V) {
    for (std::size_t a = 0; a < t; ++a)
        for (std::size_t b = 0; b < t; ++b) V[a * t + b] = a == b ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (std::size_t a = 0; a < t; ++a) {
            diag += S[a * t + a] * S[a * t + a];
            for (std::size_t b = a + 1; b < t; ++b) off += S[a * t + b] * S[a * t + b];
        }
        if (!(off > 1e-32 * diag)) break;
        for (std::size_t p = 0; p < t; ++p) {
            for (std::size_t q = p + 1; q < t; ++q) {
                const double apq = S[p * t + q];
                if (apq == 0.0) continue;
                const double theta = (S[q * t + q] - S[p * t + p]) / (2.0 * apq);
                const double tt = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(tt * tt + 1.0), s = tt * c;
                for (std::size_t k = 0; k < t; ++k) {      // columns p and q
                    const double skp = S[k * t + p], skq = S[k * t + q];
                    S[k * t + p] = c * skp - s * skq; S[k * t + q] = s * skp + c * skq;
                }
                for (std::size_t k = 0; k < t; ++k) {      // rows p and q
                    const double spk = S[p * t + k], sqk = S[q * t + k];
                    S[p * t + k] = c * spk - s * sqk; S[q * t + k] = s * spk + c * sqk;
                }
                for (std::size_t k = 0; k < t; ++k) {
                    const double vkp = V[k * t + p], vkq = V[k * t + q];
                    V[k * t + p] = c * vkp - s * vkq; V[k * t + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    for (std::size_t a = 0; a < t; ++a) w[a] = S[a * t + a];
}

// A unit radial direction when the natural one is undefined (a point exactly on the axis): the first curved axis.
void fallback_radial(std::size_t n, double* nv) {
    for (std::size_t k = 0; k < n; ++k) nv[k] = k == 0 ? 1.0 : 0.0;
}

// hit from a candidate point q (relative to the centre): distance of xc + q from the axis.
double radial(std::size_t n, std::size_t t, const double* xc, const double* q, double* nv) {
    double g2 = 0.0;
    for (std::size_t k = 0; k < t; ++k) { nv[k] = xc[k] + q[k]; g2 += nv[k] * nv[k]; }
    for (std::size_t k = t; k < n; ++k) nv[k] = 0.0;
    const double g = std::sqrt(g2);
    if (g > 0.0) for (std::size_t k = 0; k < t; ++k) nv[k] /= g; else fallback_radial(n, nv);
    return g;
}

void take_point(std::size_t n, std::size_t t, const double* xc, const double* q, CurvedWallHit& hit, bool& have) {
    double nv[MD];
    const double g = radial(n, t, xc, q, nv);
    if (!have || g > hit.g) {
        hit.g = g; have = true;
        for (std::size_t k = 0; k < n; ++k) { hit.nvec[k] = nv[k]; hit.p[k] = q[k]; }
    }
}

// Support point of the rounded rod (rod.hpp) along the unit vector m: the gradient of its support function.
void rounded_point(std::size_t n, int kind, double eps, const double* u, double hl, double rad, const double* m, double* p) {
    const double e2 = eps * eps;
    double c = 0.0;
    for (std::size_t k = 0; k < n; ++k) c += u[k] * m[k];
    const double S1 = std::sqrt(c * c + e2);
    for (std::size_t k = 0; k < n; ++k) p[k] = hl * (c * u[k] + e2 * m[k]) / S1;
    if (kind == SHAPE_ROD_CAPPED) {
        for (std::size_t k = 0; k < n; ++k) p[k] += rad * m[k];
        return;
    }
    double p2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) { const double v = m[k] - c * u[k]; p2 += v * v; }
    const double S2 = std::sqrt(p2 + e2);
    for (std::size_t k = 0; k < n; ++k) p[k] += rad * ((m[k] - c * u[k]) + e2 * m[k]) / S2;
}

}  // namespace

void ell_farthest_from_axis(std::size_t n, std::size_t t, const double* A, const double* xc, CurvedWallHit& hit) {
    // Maximise |b + q_t|^2 over q^T A^-1 q <= 1, b = the first t entries of xc. Stationarity gives q = A w / lam with
    // w = P(b + q) the radial vector, hence (lam I - A_tt) w = lam b and the constraint
    //     phi(lam) = sum_i L_i c_i^2 / (lam - L_i)^2 = 1,   lam > L_max,
    // in the eigenbasis of the t x t block A_tt = V diag(L) V^T, c = V^T b. phi falls monotonically on (L_max, inf).
    double S[MD * MD], L[MD], V[MD * MD], c[MD];
    for (std::size_t a = 0; a < t; ++a)
        for (std::size_t b = 0; b < t; ++b) S[a * t + b] = A[a * n + b];
    jacobi(t, S, L, V);
    std::size_t top = 0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < t; ++i) {
        if (L[i] > L[top]) top = i;
        c[i] = 0.0;
        for (std::size_t k = 0; k < t; ++k) c[i] += V[k * t + i] * xc[k];
        b2 += xc[i] * xc[i];
    }
    const double Lmax = L[top];
    // Terms whose eigenvalue equals L_max to rounding blow up together at lam -> L_max ("top" terms).
    auto is_top = [&](std::size_t i) { return L[i] >= Lmax * (1.0 - 1e-14); };
    double ctop2 = 0.0;
    for (std::size_t i = 0; i < t; ++i) if (is_top(i)) ctop2 += c[i] * c[i];
    auto phi_rest = [&](double mu) {       // the non-top terms at lam = L_max + mu
        double s = 0.0;
        for (std::size_t i = 0; i < t; ++i) if (!is_top(i)) { const double d = mu + Lmax - L[i]; s += L[i] * c[i] * c[i] / (d * d); }
        return s;
    };
    double w[MD] = {0.0};                 // the radial vector P(b + q), in the lab frame
    double lam;
    const double tiny = 1e-28 * (Lmax + b2);
    if (!(ctop2 > tiny) && !(phi_rest(0.0) > 1.0)) {
        // Hard case: the centre has no component along the longest projected axis (for example it is on the
        // container's axis). lam = L_max and the farthest point leans along that axis by tau; both signs are
        // equally far, the + one is taken.
        lam = Lmax;
        const double tau = std::sqrt(std::max(0.0, Lmax * (1.0 - phi_rest(0.0))));
        for (std::size_t i = 0; i < t; ++i) {
            const double coef = is_top(i) ? (i == top ? tau : 0.0) : lam * c[i] / (lam - L[i]);
            for (std::size_t k = 0; k < t; ++k) w[k] += coef * V[k * t + i];
        }
    } else {
        // Root of phi(L_max + mu) = 1 in mu > 0. Bracket: phi > 1 at mu_lo, phi <= 1 at mu_hi. Newton on
        // psi = phi^(-1/2) - 1 (nearly linear in mu), kept inside the bracket, else bisection.
        auto phi = [&](double mu, double* dphi) {
            double s = 0.0, ds = 0.0;
            for (std::size_t i = 0; i < t; ++i) {
                const double d = mu + Lmax - L[i];
                if (c[i] == 0.0) continue;
                const double term = L[i] * c[i] * c[i] / (d * d);
                s += term; ds += -2.0 * term / d;
            }
            if (dphi) *dphi = ds;
            return s;
        };
        double mu_hi = std::sqrt(Lmax * b2) + 1e-300, mu_lo = 0.0;      // phi(mu_hi) <= L_max b^2 / mu_hi^2 = 1
        double mu = mu_hi;
        for (int it = 0; it < 200; ++it) {
            double dphi;
            const double f = phi(mu, &dphi);
            if (f > 1.0) mu_lo = mu; else mu_hi = mu;
            if (std::abs(f - 1.0) <= 4e-16) break;
            double next = -1.0;
            if (f > 0.0 && dphi < 0.0) {
                const double psi = 1.0 / std::sqrt(f) - 1.0, dpsi = -0.5 * dphi / (f * std::sqrt(f));
                next = mu - psi / dpsi;
            }
            if (!(next > mu_lo && next < mu_hi)) next = 0.5 * (mu_lo + mu_hi);
            if (next == mu || !(mu_hi - mu_lo > 1e-16 * mu_hi)) { mu = next; break; }
            mu = next;
        }
        lam = Lmax + mu;
        // Near the hard case (centre close to the axis) mu is tiny and the top coefficients lam c_i / mu amplify
        // the rounding of the eigenvalues: an error of 1e-16 L_max in L_max is a RELATIVE error 1e-16 L_max / mu in
        // them. There the size of the top part is taken from the constraint instead (as in the hard case), with
        // the direction of c inside the top eigenspace; the other coefficients are well conditioned.
        const double share = 1.0 - phi_rest(mu);          // what the top terms must contribute to phi = 1
        const bool top_from_constraint = ctop2 > 0.0 && share > 0.25;
        const double top_scale = top_from_constraint ? lam * std::sqrt(share / (Lmax * ctop2)) : 0.0;
        for (std::size_t i = 0; i < t; ++i) {
            if (c[i] == 0.0) continue;
            // mu + (L_max - L_i), not lam - L_i: the sum L_max + mu loses mu
            const double coef = (top_from_constraint && is_top(i)) ? top_scale * c[i] : lam * c[i] / (mu + (Lmax - L[i]));
            for (std::size_t k = 0; k < t; ++k) w[k] += coef * V[k * t + i];
        }
    }
    double g2 = 0.0;
    for (std::size_t k = 0; k < t; ++k) g2 += w[k] * w[k];
    hit.g = std::sqrt(g2);
    for (std::size_t k = 0; k < n; ++k) {
        hit.nvec[k] = (k < t && hit.g > 0.0) ? w[k] / hit.g : 0.0;
        double q = 0.0;
        for (std::size_t m = 0; m < t; ++m) q += A[k * n + m] * w[m];
        hit.p[k] = q / lam;                // q = A w / lam (lam >= L_max > 0)
    }
    if (!(hit.g > 0.0)) fallback_radial(n, hit.nvec);
}

namespace {
// End rim of a 3D flat-ended rod: the circle of radius rad about the point e (relative to the container's centre),
// across the axis u. Its distance from the container axis, squared, is a degree-2 trigonometric polynomial of the
// angle, so it has at most two local maxima. They are found by sampling and Newton on f', and returned as rim points
// RELATIVE TO e. Each is a stationary point, so a contact placed there has an exact gradient (envelope theorem).
// A rim that is equally far all round (end-on in a sphere container) returns two opposite points.
int rim_local_maxima(std::size_t t, const double* u, double rad, const double* e, double pts[2][3]) {
    double a[3], b[3];
    {
        const double ax[3] = {std::abs(u[0]) < 0.6 ? 1.0 : 0.0, std::abs(u[0]) < 0.6 ? 0.0 : 1.0, 0.0};
        const double d = ax[0] * u[0] + ax[1] * u[1] + ax[2] * u[2];
        double nrm = 0.0;
        for (int k = 0; k < 3; ++k) { a[k] = ax[k] - d * u[k]; nrm += a[k] * a[k]; }
        nrm = std::sqrt(nrm);
        for (int k = 0; k < 3; ++k) a[k] /= nrm;
        b[0] = u[1] * a[2] - u[2] * a[1]; b[1] = u[2] * a[0] - u[0] * a[2]; b[2] = u[0] * a[1] - u[1] * a[0];
    }
    // f(th) = |P(e + rad (cos th a + sin th b))|^2 - |P e|^2 = 2 rad (ea cos + eb sin) + rad^2 (aa cos^2 + 2 ab sin cos + bb sin^2)
    double ea = 0.0, eb = 0.0, aa = 0.0, ab = 0.0, bb = 0.0, e2 = 0.0;
    for (std::size_t k = 0; k < t; ++k) { ea += e[k] * a[k]; eb += e[k] * b[k]; aa += a[k] * a[k]; ab += a[k] * b[k]; bb += b[k] * b[k]; e2 += e[k] * e[k]; }
    auto f = [&](double th) {
        const double cs = std::cos(th), sn = std::sin(th);
        return 2.0 * rad * (ea * cs + eb * sn) + rad * rad * (aa * cs * cs + 2.0 * ab * sn * cs + bb * sn * sn);
    };
    auto point = [&](double th, double* p) {
        const double cs = std::cos(th), sn = std::sin(th);
        for (int k = 0; k < 3; ++k) p[k] = rad * (cs * a[k] + sn * b[k]);
    };
    constexpr int NS = 24;
    double fs[NS], fmin = std::numeric_limits<double>::infinity(), fmax = -fmin;
    for (int s = 0; s < NS; ++s) { fs[s] = f(2.0 * PI * s / NS); fmin = std::min(fmin, fs[s]); fmax = std::max(fmax, fs[s]); }
    if (!(fmax - fmin > 1e-13 * (e2 + rad * rad))) {       // equally far all round: two opposite points
        point(0.0, pts[0]); point(PI, pts[1]);
        return 2;
    }
    double best[2] = {0.0, 0.0}, fbest[2] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (int s = 0; s < NS; ++s) {
        if (!(fs[s] >= fs[(s + 1) % NS] && fs[s] > fs[(s + NS - 1) % NS])) continue;     // local maximum of the samples
        const int slot = fs[s] > fbest[0] ? 0 : (fs[s] > fbest[1] ? 1 : -1);
        if (slot == 0) { best[1] = best[0]; fbest[1] = fbest[0]; }
        if (slot >= 0) { best[slot] = 2.0 * PI * s / NS; fbest[slot] = fs[s]; }
    }
    int m = 0;
    double found[2] = {0.0, 0.0};
    for (int cand = 0; cand < 2; ++cand) {
        if (!std::isfinite(fbest[cand])) continue;
        double th = best[cand], fth = fbest[cand];
        for (int it = 0; it < 40; ++it) {
            const double cs = std::cos(th), sn = std::sin(th);
            const double f1 = 2.0 * rad * (-ea * sn + eb * cs) + rad * rad * ((bb - aa) * 2.0 * sn * cs + 2.0 * ab * (cs * cs - sn * sn));
            const double f2 = 2.0 * rad * (-ea * cs - eb * sn) + rad * rad * ((bb - aa) * 2.0 * (cs * cs - sn * sn) - 8.0 * ab * sn * cs);
            if (!(f2 < 0.0)) break;
            double step = -f1 / f2;
            const double lim = PI / NS;                 // stay within one sample spacing
            step = std::max(-lim, std::min(lim, step));
            const double fn = f(th + step);
            if (!(fn >= fth)) break;
            th += step; fth = fn;
            if (std::abs(step) < 1e-15) break;
        }
        if (m == 1 && std::abs(std::remainder(th - found[0], 2.0 * PI)) < 1e-8) continue;   // the same maximum twice
        found[m] = th; point(th, pts[m]); ++m;
    }
    return m;
}
}  // namespace

void rod_farthest_from_axis(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, double eps,
                            const double* xc, CurvedWallHit& hit) {
    bool have = false;
    double q[MD];
    if (eps > 0.0) {
        // Rounded rod: ascent on the support function. From a radial direction m, the body's support point p(m)
        // is at least as far from the axis as the point that defined m, so g never decreases; the fixed point is
        // a farthest point. Both ends are tried as starts (the two lobes), the farther result is kept.
        for (int side = 0; side < 2; ++side) {
            double m[MD], p[MD];
            for (std::size_t k = 0; k < n; ++k) q[k] = (side == 0 ? hl : -hl) * u[k];
            double g = radial(n, t, xc, q, m), g_prev = -1.0;
            for (int it = 0; it < 500; ++it) {
                rounded_point(n, kind, eps, u, hl, rad, m, p);
                g_prev = g;
                g = radial(n, t, xc, p, m);
                if (!(g > g_prev * (1.0 + 1e-15))) break;
            }
            rounded_point(n, kind, eps, u, hl, rad, m, p);
            take_point(n, t, xc, p, hit, have);
        }
        return;
    }
    if (kind == SHAPE_ROD_CAPPED) {
        // The farthest point of a segment is one of its ends; the cap adds the radius along the radial direction.
        for (int side = 0; side < 2; ++side) {
            double nv[MD];
            for (std::size_t k = 0; k < n; ++k) q[k] = (side == 0 ? hl : -hl) * u[k];
            const double ge = radial(n, t, xc, q, nv);
            for (std::size_t k = 0; k < n; ++k) q[k] += rad * nv[k];
            if (!have || ge + rad > hit.g) {
                hit.g = ge + rad; have = true;
                for (std::size_t k = 0; k < n; ++k) { hit.nvec[k] = nv[k]; hit.p[k] = q[k]; }
            }
        }
        return;
    }
    // Uncapped (2D and 3D only): the farthest point is an extreme point: a corner (2D) or a point of an end rim (3D).
    if (n == 2) {
        const double v[2] = {-u[1], u[0]};
        for (int side = 0; side < 2; ++side)
            for (int cs = 0; cs < 2; ++cs) {
                for (std::size_t k = 0; k < 2; ++k) q[k] = (side == 0 ? hl : -hl) * u[k] + (cs == 0 ? rad : -rad) * v[k];
                take_point(n, t, xc, q, hit, have);
            }
        return;
    }
    // 3D: the farthest point is a local maximum of one of the two end rims.
    for (int side = 0; side < 2; ++side) {
        double e[3], pts[2][3];
        for (int k = 0; k < 3; ++k) e[k] = xc[k] + (side == 0 ? hl : -hl) * u[k];
        const int m = rim_local_maxima(t, u, rad, e, pts);
        for (int c = 0; c < m; ++c) {
            for (int k = 0; k < 3; ++k) q[k] = (side == 0 ? hl : -hl) * u[k] + pts[c][k];
            take_point(n, t, xc, q, hit, have);
        }
    }
}

bool curved_wall_contact(std::size_t n, double Rc, double K, const CurvedWallHit& hit, CurvedWallContact& c) {
    if (!(hit.g > Rc)) return false;
    double h = 0.0;
    for (std::size_t k = 0; k < n; ++k) h += hit.nvec[k] * hit.p[k];
    const double delta = h - (hit.g - Rc);
    c.h = h; c.delta = delta;
    c.F_mag = -2.0 * K * std::abs(h / delta - 1.0);
    c.U_add = 1.0 - delta / h;
    c.dkappa = 2.0 * K * h * (delta - h);
    c.energy = K * (h - delta) * (h - delta);
    c.f = 2.0 * K * (h - delta);
    return true;
}

void curved_wall_apply(std::size_t n, const CurvedWallHit& hit, const CurvedWallContact& c, double* F, double* Tq) {
    for (std::size_t k = 0; k < n; ++k) F[k] -= c.f * hit.nvec[k];
    std::size_t idx = 0;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b, ++idx) Tq[idx] -= c.f * (hit.nvec[a] * hit.p[b] - hit.nvec[b] * hit.p[a]);
}

namespace {
template <class G>
double wall_scale(G g_of, double Rc, double s_max) {
    if (g_of(s_max) <= Rc) return s_max;
    double lo = 0.0, hi = s_max;           // g(0) = distance of the centre (< Rc for a centre inside), g(hi) > Rc
    if (!(g_of(0.0) < Rc)) return 0.0;
    for (int it = 0; it < 200 && hi - lo > 1e-16 * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (g_of(mid) <= Rc) lo = mid; else hi = mid;
    }
    return lo;
}
}  // namespace

double ell_wall_scale(std::size_t n, std::size_t t, const double* A, const double* xc, double Rc, double s_max) {
    return wall_scale([&](double s) {
        double As[MD * MD];
        for (std::size_t k = 0; k < n * n; ++k) As[k] = s * s * A[k];
        CurvedWallHit hit;
        if (s == 0.0) { double q[MD] = {0.0}, nv[MD]; return radial(n, t, xc, q, nv); }
        ell_farthest_from_axis(n, t, As, xc, hit);
        return hit.g;
    }, Rc, s_max);
}

namespace {
// One point contact: the point q (relative to the centre) carries a ball of radius `ball` (0: a bare point).
void point_contact(std::size_t n, std::size_t t, const double* xc, const double* q, double ball, int sign, double Rc, double K,
                   double* F, double* Tq, CurvedWallSum& sum, bool& any) {
    CurvedWallHit hit;
    const double ge = radial(n, t, xc, q, hit.nvec);
    hit.g = ge + ball;
    if (!(hit.g > Rc)) return;
    for (std::size_t k = 0; k < n; ++k) hit.p[k] = q[k] + ball * hit.nvec[k];
    CurvedWallContact c;
    // h (the sphere code's r) for a point contact: the ball's radius for a cap; for a bare corner or rim point the
    // extent of the point along the radial direction, as for the ellipsoid.
    if (!curved_wall_contact(n, Rc, K, hit, c)) return;
    any = true;
    sum.U += sign * c.U_add; sum.dkappa += sign * c.dkappa; sum.Fmean += sign * c.F_mag; sum.energy += sign * c.energy;
    sum.f += sign * c.f; sum.count += sign;
    if (F) for (std::size_t k = 0; k < n; ++k) F[k] -= sign * c.f * hit.nvec[k];
    if (Tq) {
        std::size_t idx = 0;
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = a + 1; b < n; ++b, ++idx) Tq[idx] -= sign * c.f * (hit.nvec[a] * hit.p[b] - hit.nvec[b] * hit.p[a]);
    }
}
}  // namespace

namespace {
// Number of contact points on each end rim of a 3D flat-ended rod at a curved wall: 24. The environment variable
// RCP_ROD_RIM_POINTS (4 to 256, read once) changes it for MEASUREMENTS of its effect (LOG [038]); it is not an
// interface and is not set in normal use.
constexpr int RIM_MAX = 256;
struct RimTable {
    int n = 24;
    double c[RIM_MAX], s[RIM_MAX];
    RimTable() {
        if (const char* e = std::getenv("RCP_ROD_RIM_POINTS")) { const int v = std::atoi(e); if (v >= 4 && v <= RIM_MAX) n = v; }
        for (int m = 0; m < n; ++m) { c[m] = std::cos(2.0 * PI * m / n); s[m] = std::sin(2.0 * PI * m / n); }
    }
};
const RimTable RIM_TABLE;
const int RIM_POINTS = RIM_TABLE.n;
const double* const RIM_COS = RIM_TABLE.c;
const double* const RIM_SIN = RIM_TABLE.s;
}  // namespace

bool rod_curved_wall(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, const double* xc,
                     double Rc, double K, double* F, double* Tq, CurvedWallSum& sum, const double* Rot) {
    bool any = false;
    double q[MD];
    if (kind == SHAPE_ROD_CAPPED) {
        for (int side = 0; side < 2; ++side) {
            for (std::size_t k = 0; k < n; ++k) q[k] = (side == 0 ? hl : -hl) * u[k];
            point_contact(n, t, xc, q, rad, +1, Rc, K, F, Tq, sum, any);
        }
        for (std::size_t k = 0; k < n; ++k) q[k] = 0.0;
        point_contact(n, t, xc, q, rad, -1, Rc, K, F, Tq, sum, any);
        return any;
    }
    if (n == 2) {
        const double v[2] = {-u[1], u[0]};
        for (int side = 0; side < 2; ++side)
            for (int cs = 0; cs < 2; ++cs) {
                for (std::size_t k = 0; k < 2; ++k) q[k] = (side == 0 ? hl : -hl) * u[k] + (cs == 0 ? rad : -rad) * v[k];
                point_contact(n, t, xc, q, 0.0, +1, Rc, K, F, Tq, sum, any);
            }
        return any;
    }
    // 3D uncapped: RIM_POINTS points on each end rim, fixed in the rod's own frame (body axes 2 and 3), on the
    // circle that circumscribes the rim's regular polygon through them (radius rad / cos(pi / RIM_POINTS)), so the
    // rim lies inside the polygon and no part of it can be outside the wall unseen.
    double a[3], b[3];
    if (Rot) { for (int k = 0; k < 3; ++k) { a[k] = Rot[k * 3 + 1]; b[k] = Rot[k * 3 + 2]; } }
    else {      // no frame given (single-rod tests that pass only the axis): any pair across the axis
        const double ax[3] = {std::abs(u[0]) < 0.6 ? 1.0 : 0.0, std::abs(u[0]) < 0.6 ? 0.0 : 1.0, 0.0};
        const double d = ax[0] * u[0] + ax[1] * u[1] + ax[2] * u[2];
        double nrm = 0.0;
        for (int k = 0; k < 3; ++k) { a[k] = ax[k] - d * u[k]; nrm += a[k] * a[k]; }
        nrm = std::sqrt(nrm);
        for (int k = 0; k < 3; ++k) a[k] /= nrm;
        b[0] = u[1] * a[2] - u[2] * a[1]; b[1] = u[2] * a[0] - u[0] * a[2]; b[2] = u[0] * a[1] - u[1] * a[0];
    }
    const double rr = rad / std::cos(PI / RIM_POINTS);
    for (int side = 0; side < 2; ++side) {
        const double hs = side == 0 ? hl : -hl;
        for (int m = 0; m < RIM_POINTS; ++m) {
            const double cs = RIM_COS[m], sn = RIM_SIN[m];
            for (int k = 0; k < 3; ++k) q[k] = hs * u[k] + rr * (cs * a[k] + sn * b[k]);
            point_contact(n, t, xc, q, 0.0, +1, Rc, K, F, Tq, sum, any);
        }
    }
    return any;
}

bool shape_curved_wall(const EllShape& sh, std::size_t n, std::size_t t, double Rc, double K, const double* x_i, double Db_i,
                       const double* Amat_i, const double* Rot_i, const double* q_i, double* F_i, double* Tq_i,
                       CurvedWallSum& sum) {
    double xc[MD], d2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) { xc[k] = k < t ? x_i[k] - Rc : 0.0; d2 += xc[k] * xc[k]; }
    if (std::sqrt(d2) + 0.5 * Db_i <= Rc) return false;      // bounding sphere inside
    if (sh.kind == SHAPE_ELLIPSOID) {
        CurvedWallHit hit;
        ell_farthest_from_axis(n, t, Amat_i, xc, hit);
        CurvedWallContact c;
        if (!curved_wall_contact(n, Rc, K, hit, c)) return false;
        sum.U += c.U_add; sum.dkappa += c.dkappa; sum.Fmean += c.F_mag; sum.energy += c.energy; sum.f += c.f; sum.count += 1;
        if (F_i && Tq_i) curved_wall_apply(n, hit, c, F_i, Tq_i);
        return true;
    }
    double u[MD];
    for (std::size_t k = 0; k < n; ++k) u[k] = Rot_i[k * n];                              // the rod's axis: column 0
    return rod_curved_wall(n, t, sh.kind, u, Db_i * q_i[0], Db_i * q_i[1], xc, Rc, K, F_i, Tq_i, sum, Rot_i);
}

double rod_wall_scale(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, const double* xc,
                      double Rc, double s_max) {
    return wall_scale([&](double s) {
        CurvedWallHit hit;
        rod_farthest_from_axis(n, t, kind, u, s * hl, s * rad, 0.0, xc, hit);
        return hit.g;
    }, Rc, s_max);
}

}  // namespace rcpgenerator
