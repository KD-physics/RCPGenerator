// Rods: capped (spherocylinders) and uncapped (flat-ended cylinders). See rcpgenerator/rod.hpp for the conventions.
// Additive: nothing in this file is reached for spheres or for ellipsoids.

#include "rcpgenerator/rod.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

#include <omp.h>

#include "rcpgenerator/curved_wall.hpp"

// MSVC's <cmath> does not define M_PI; the sphere code's headers define it the same way (rcp_generator.hpp), but
// this file does not include them.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace rcpgenerator {

namespace {

constexpr std::size_t MD = ELL_MAX_DIM;
thread_local RodKernelStats g_stats;   // per thread: the force pass may run on several; the diagnostic reads the caller's

#if defined(__GNUC__) || defined(__clang__)
#define ROD_INLINE inline __attribute__((always_inline))
#else
#define ROD_INLINE inline
#endif

ROD_INLINE double dot(const double* a, const double* b, std::size_t n) {
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) s += a[k] * b[k];
    return s;
}
ROD_INLINE double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }
inline double sgn(double v) { return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0); }

double ball_coeff(std::size_t k) {   // volume of the unit k-ball
    return std::pow(M_PI, 0.5 * static_cast<double>(k)) / std::tgamma(0.5 * static_cast<double>(k) + 1.0);
}

// Generalised torque of a force f applied at the point p (relative to the particle's centre):
// T_ab = f_a p_b - f_b p_a for a < b (T_ab = -dE/dOmega_ab with dp = Omega p and f = -dE/dp).
inline void add_torque(std::size_t n, const double* f, const double* p, double scale, double* T) {
    std::size_t c = 0;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b, ++c) T[c] += scale * (f[a] * p[b] - f[b] * p[a]);
}

// ---- capped rods: closest points of two segments -----------------------------------------------------------------
// Segment i: t_i u_i, |t_i| <= hli (centre at the origin); segment j: r + t_j u_j, |t_j| <= hlj.
// Returns the squared distance and the parameters. Exactly parallel segments: the middle of the stretch over
// which they face each other (so the pair force acts where a uniform line contact would put it).
ROD_INLINE double segment_closest(std::size_t n, const double* r, const double* ui, double hli, const double* uj, double hlj,
                       double& ti, double& tj) {
    const double b = dot(ui, uj, n), c = dot(ui, r, n), f = dot(uj, r, n);
    const double denom = 1.0 - b * b;
    if (denom > 1e-12) {
        ti = clampd((c - b * f) / denom, -hli, hli);
        tj = b * ti - f;
        if (tj < -hlj) { tj = -hlj; ti = clampd(b * tj + c, -hli, hli); }
        else if (tj > hlj) { tj = hlj; ti = clampd(b * tj + c, -hli, hli); }
    } else {
        // (Nearly) parallel: the general formula divides by ~0. The closest approach is at one end of the stretch
        // over which the segments face each other (or anywhere on it when exactly parallel). Take the closer end
        // for the DISTANCE, so that it is the true minimum also for a tilt below the threshold; when both ends are
        // equally close (exactly parallel) the points are put at the middle of the stretch.
        const double lo = std::max(-hli, c - hlj), hi = std::min(hli, c + hlj);
        if (lo <= hi) {
            auto dist2 = [&](double a, double& bj) {
                bj = clampd(b * a - f, -hlj, hlj);
                double q = 0.0;
                for (std::size_t k = 0; k < n; ++k) { const double w = r[k] + bj * uj[k] - a * ui[k]; q += w * w; }
                return q;
            };
            double tjl, tjh;
            const double dl = dist2(lo, tjl), dh = dist2(hi, tjh);
            if (std::abs(dl - dh) <= 1e-14 * (dl + dh)) { ti = 0.5 * (lo + hi); tj = clampd(b * ti - f, -hlj, hlj); return 0.5 * (dl + dh) < std::min(dl, dh) ? 0.5 * (dl + dh) : std::min(dl, dh); }
            if (dl < dh) { ti = lo; tj = tjl; return dl; }
            ti = hi; tj = tjh; return dh;
        }
        ti = clampd(c, -hli, hli);
        tj = clampd(b * ti - f, -hlj, hlj);
    }
    double d2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) { const double w = r[k] + tj * uj[k] - ti * ui[k]; d2 += w * w; }
    return d2;
}

// 2D, crossing segments: among the four ends, the one closest to the other segment. Returns minus that distance,
// sets (ti, tj) to the two points (the end, and its foot on the other segment) and nh so that the signed distance
// obeys d = nh . (r + tj uj - ti ui).
double signed_crossing_2d(const double* r, const double* ui, double hli, const double* uj, double hlj, double& ti,
                          double& tj, double* nh) {
    double best = std::numeric_limits<double>::infinity();
    for (int which = 0; which < 4; ++which) {
        const bool end_of_i = which < 2;
        const double sg = (which & 1) ? 1.0 : -1.0;
        double e[2], t_other, w[2];
        if (end_of_i) {            // end of i against segment j: foot parameter on j
            e[0] = sg * hli * ui[0] - r[0]; e[1] = sg * hli * ui[1] - r[1];           // the end, relative to centre j
            t_other = clampd(e[0] * uj[0] + e[1] * uj[1], -hlj, hlj);
            w[0] = t_other * uj[0] - e[0]; w[1] = t_other * uj[1] - e[1];             // foot on j minus the end of i
        } else {                   // end of j against segment i
            e[0] = r[0] + sg * hlj * uj[0]; e[1] = r[1] + sg * hlj * uj[1];           // the end, relative to centre i
            t_other = clampd(e[0] * ui[0] + e[1] * ui[1], -hli, hli);
            w[0] = e[0] - t_other * ui[0]; w[1] = e[1] - t_other * ui[1];             // end of j minus its foot on i
        }
        const double dist = std::sqrt(w[0] * w[0] + w[1] * w[1]);
        if (dist < best) {
            best = dist;
            if (end_of_i) { ti = sg * hli; tj = t_other; } else { ti = t_other; tj = sg * hlj; }
            if (dist > 0.0) { nh[0] = -w[0] / dist; nh[1] = -w[1] / dist; } else { nh[0] = -ui[1]; nh[1] = ui[0]; }
        }
    }
    return -best;
}

// ---- uncapped rods ------------------------------------------------------------------------------------------------
// Support function of the cylinder {t u + rho: |t| <= hl, rho perpendicular to u, |rho| <= rad} along m (any length).
ROD_INLINE double cyl_support(std::size_t n, const double* u, double hl, double rad, const double* m) {
    const double c = dot(u, m, n);
    double p2 = 0.0;   // |m - (m.u) u|^2 from the components: m.m - c^2 would lose half the digits when m is along u
    for (std::size_t k = 0; k < n; ++k) { const double v = m[k] - c * u[k]; p2 += v * v; }
    return hl * std::abs(c) + rad * std::sqrt(p2);
}

// The support set of a cylinder along the unit normal nv is a point, a segment (nv perpendicular to the axis) or
// an end disc (nv along the axis). project() returns the point of that set closest to y.
struct SupportSet {
    int type = 0;                // 0 point, 1 segment (side), 2 disc (end face)
    double base[MD] = {0};       // point; segment: its middle (rad * v); disc: its centre (hl * sgn * u)
    const double* u = nullptr;
    double hl = 0.0, rad = 0.0;
    std::size_t n = 0;
    void project(const double* y, double* out) const {
        if (type == 0) { for (std::size_t k = 0; k < n; ++k) out[k] = base[k]; return; }
        const double uy = dot(u, y, n);
        if (type == 1) {
            const double t = clampd(uy, -hl, hl);
            for (std::size_t k = 0; k < n; ++k) out[k] = base[k] + t * u[k];
            return;
        }
        double rho[MD], r2 = 0.0;
        for (std::size_t k = 0; k < n; ++k) { rho[k] = y[k] - uy * u[k]; r2 += rho[k] * rho[k]; }
        const double sc = r2 > rad * rad ? rad / std::sqrt(r2) : 1.0;
        for (std::size_t k = 0; k < n; ++k) out[k] = base[k] + sc * rho[k];
    }
};

SupportSet support_set(std::size_t n, const double* u, double hl, double rad, const double* nv) {
    SupportSet s;
    s.u = u; s.hl = hl; s.rad = rad; s.n = n;
    const double c = dot(u, nv, n);
    double perp[MD], p2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) { perp[k] = nv[k] - c * u[k]; p2 += perp[k] * perp[k]; }
    const double sp = std::sqrt(p2);
    const double tol = 1e-9;
    if (sp <= tol) {                       // normal along the axis: the end face
        s.type = 2;
        for (std::size_t k = 0; k < n; ++k) s.base[k] = hl * sgn(c) * u[k];
    } else if (std::abs(c) <= tol) {       // normal perpendicular to the axis: a line of the side
        s.type = 1;
        for (std::size_t k = 0; k < n; ++k) s.base[k] = rad * perp[k] / sp;
    } else {
        for (std::size_t k = 0; k < n; ++k) s.base[k] = hl * sgn(c) * u[k] + rad * perp[k] / sp;
    }
    return s;
}

// Contact points: p_i in the support set of i and q_j in that of j (both along +nv; j is symmetric) with
// p_i + q_j = target (= r / s): where the two scaled rods touch. Unique points are taken as they are; where a
// set has room (a face, a line of the side) the point is the one consistent with the other rod.
void contact_points(std::size_t n, const double* ui, double hli, double radi, const double* uj, double hlj,
                    double radj, const double* nv, const double* target, double* pi, double* qj) {
    const SupportSet A = support_set(n, ui, hli, radi, nv), B = support_set(n, uj, hlj, radj, nv);
    double y[MD];
    for (std::size_t k = 0; k < n; ++k) { pi[k] = A.base[k]; qj[k] = B.base[k]; }
    if (A.type == 0 && B.type == 0) return;
    if (A.type == 1 && B.type == 1) {
        // side against side: t_i u_i + t_j u_j = target - (the two side offsets), solved in the plane of the axes;
        // parallel axes: the middle of the stretch over which the two sides face each other
        double t[MD];
        for (std::size_t k = 0; k < n; ++k) t[k] = target[k] - A.base[k] - B.base[k];
        const double b = dot(ui, uj, n), ta_ = dot(t, ui, n), tb_ = dot(t, uj, n), den = 1.0 - b * b;
        double ta, tb;
        if (den > 1e-12) {
            ta = clampd((ta_ - b * tb_) / den, -hli, hli);
            tb = clampd((tb_ - b * ta_) / den, -hlj, hlj);
        } else {
            const double lo = std::max(-hli, ta_ - hlj), hi = std::min(hli, ta_ + hlj);
            ta = lo <= hi ? 0.5 * (lo + hi) : clampd(ta_, -hli, hli);
            tb = clampd(b * (ta_ - ta), -hlj, hlj);
        }
        for (std::size_t k = 0; k < n; ++k) { pi[k] = A.base[k] + ta * ui[k]; qj[k] = B.base[k] + tb * uj[k]; }
        return;
    }
    const int sweeps = (A.type != 0 && B.type != 0) ? 12 : 1;
    for (int it = 0; it < sweeps; ++it) {
        if (A.type != 0) { for (std::size_t k = 0; k < n; ++k) y[k] = target[k] - qj[k]; A.project(y, pi); }
        if (B.type != 0) { for (std::size_t k = 0; k < n; ++k) y[k] = target[k] - pi[k]; B.project(y, qj); }
    }
}

// 2D: the Minkowski sum of two rectangles is a polygon whose edge normals are the four edge normals of the two
// rectangles, so the maximum over n of n.r / (h_i + h_j) is at one of them. Exact.
double scale_uncapped_2d(const double* r, const double* ui, double hli, double radi, const double* uj, double hlj,
                         double radj, double* nv) {
    const double cand[4][2] = {{ui[0], ui[1]}, {-ui[1], ui[0]}, {uj[0], uj[1]}, {-uj[1], uj[0]}};
    double best = -1.0;
    for (int k = 0; k < 4; ++k) {
        double m[2] = {cand[k][0], cand[k][1]};
        double nr = m[0] * r[0] + m[1] * r[1];
        if (nr < 0.0) { nr = -nr; m[0] = -m[0]; m[1] = -m[1]; }
        const double H = cyl_support(2, ui, hli, radi, m) + cyl_support(2, uj, hlj, radj, m);
        const double ratio = nr / H;
        if (ratio > best) { best = ratio; nv[0] = m[0]; nv[1] = m[1]; }
    }
    return best;
}

inline void cross3(const double* a, const double* b, double* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

// Lean re-check of a closed-form candidate that held at the last pass (same conditions as in scale_uncapped_3d).
// End face of a (axis ua) against b: true when the contact point still falls inside a's end disc.
ROD_INLINE bool face_still_holds(const double* r, const double* ua, double hla, double rada, const double* ub, double hlb,
                                 double radb, double* nv, double& s) {
    const double ar = ua[0] * r[0] + ua[1] * r[1] + ua[2] * r[2];
    const double sg = ar < 0.0 ? -1.0 : 1.0;
    const double nn[3] = {sg * ua[0], sg * ua[1], sg * ua[2]};
    const double c = ub[0] * nn[0] + ub[1] * nn[1] + ub[2] * nn[2];
    const double perp[3] = {nn[0] - c * ub[0], nn[1] - c * ub[1], nn[2] - c * ub[2]};
    const double sp = std::sqrt(perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]);
    if (sp <= 1e-9 || std::abs(c) <= 1e-9) return false;          // b's support is not a single point: the general route
    const double H = hla + hlb * std::abs(c) + radb * sp;         // h_a(n) + h_b(n), unit n
    const double aar = std::abs(ar);
    if (!(aar > 0.0)) return false;
    const double inv_s = H / aar, kq = radb / sp, sc = c < 0.0 ? -hlb : hlb;
    double rho2 = 0.0, along = 0.0, y[3];
    for (int k = 0; k < 3; ++k) { y[k] = r[k] * inv_s - hla * nn[k] - (sc * ub[k] + kq * perp[k]); along += ua[k] * y[k]; }
    for (int k = 0; k < 3; ++k) { const double v = y[k] - along * ua[k]; rho2 += v * v; }
    if (!(rho2 <= rada * rada * (1.0 + 1e-10))) return false;
    nv[0] = nn[0]; nv[1] = nn[1]; nv[2] = nn[2];
    s = aar / H;
    return true;
}
// Side against side: true when the two contact lines still cross inside both rods.
ROD_INLINE bool side_still_holds(const double* r, const double* ui, double hli, double radi, const double* uj, double hlj,
                                 double radj, double* nv, double& s) {
    const double w[3] = {ui[1] * uj[2] - ui[2] * uj[1], ui[2] * uj[0] - ui[0] * uj[2], ui[0] * uj[1] - ui[1] * uj[0]};
    const double w2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
    if (w2 < 1e-12) return false;
    const double wr = w[0] * r[0] + w[1] * r[1] + w[2] * r[2];
    const double iw = 1.0 / std::sqrt(w2), awr = std::abs(wr) * iw;
    if (!(awr > 0.0)) return false;
    const double sg = wr < 0.0 ? -iw : iw, rho = radi + radj, inv_s = rho / awr;
    double ta_ = 0.0, tb_ = 0.0, b = 0.0;
    for (int k = 0; k < 3; ++k) { const double t = r[k] * inv_s - rho * sg * w[k]; ta_ += t * ui[k]; tb_ += t * uj[k]; b += ui[k] * uj[k]; }
    const double ta = (ta_ - b * tb_) / w2, tb = (tb_ - b * ta_) / w2;
    if (!(std::abs(ta) <= hli * (1.0 + 1e-10) && std::abs(tb) <= hlj * (1.0 + 1e-10))) return false;
    nv[0] = sg * w[0]; nv[1] = sg * w[1]; nv[2] = sg * w[2];
    s = awr / rho;
    return true;
}

// 3D: minimise H(m) = h_i(m) + h_j(m) subject to m.r = 1 (convex; s = 1 / H at the minimum). Closed-form
// candidates first, each accepted as THE minimum when its contact point falls inside the feature it assumes (the
// optimality condition of the convex problem): an end face of i or of j (m along that axis) and side against side
// (m along u_i x u_j). Otherwise a Newton iteration on the two free coordinates, stepping no further than the next
// kink (a normal becoming perpendicular to an axis) and deciding there, from the one-sided slopes, whether to stay
// on the kink. Every m visited gives a valid lower bound on s; the best is returned. `converged` reports whether
// an optimality condition was met.
double scale_uncapped_3d(const double* r, const double* ui, double hli, double radi, const double* uj, double hlj,
                         double radj, int max_iter, double* nv, bool& converged, bool stop_apart = false,
                         int* tag = nullptr, double* nc = nullptr, double tol = 0.0) {
    const std::size_t n = 3;
    if (tag && *tag >= 1 && *tag <= 3) {   // the closed form that held last time, re-checked on its own
        double s_fast;
        const bool held = *tag == 1 ? face_still_holds(r, ui, hli, radi, uj, hlj, radj, nv, s_fast)
                        : *tag == 2 ? face_still_holds(r, uj, hlj, radj, ui, hli, radi, nv, s_fast)
                                    : side_still_holds(r, ui, hli, radi, uj, hlj, radj, nv, s_fast);
        if (held) { converged = true; ++g_stats.closed_form; ++g_stats.lean_ok; return s_fast; }
        ++g_stats.lean_failed;
    }
    const double d2 = dot(r, r, n), d = std::sqrt(d2);
    converged = false;
    auto Hof = [&](const double* m) { return cyl_support(n, ui, hli, radi, m) + cyl_support(n, uj, hlj, radj, m); };
    auto take_unit = [&](const double* m) {
        const double nm = std::sqrt(dot(m, m, n));
        for (std::size_t k = 0; k < n; ++k) nv[k] = m[k] / nm;
    };
    double w[3];
    cross3(ui, uj, w);
    const double w2 = dot(w, w, n);
    const double eps = 1e-10;
    // parallel axes: the Minkowski sum is one cylinder (half-length hli + hlj, radius radi + radj)
    if (w2 < 1e-20) {
        const double c = dot(ui, r, n);
        double pr[3], p2 = 0.0;
        for (std::size_t k = 0; k < n; ++k) { pr[k] = r[k] - c * ui[k]; p2 += pr[k] * pr[k]; }
        const double s_end = std::abs(c) / (hli + hlj), s_side = std::sqrt(p2) / (radi + radj);
        converged = true;
        if (s_end >= s_side) { for (std::size_t k = 0; k < n; ++k) nv[k] = sgn(c) * ui[k]; return s_end; }
        const double ip = 1.0 / std::sqrt(p2);
        for (std::size_t k = 0; k < n; ++k) nv[k] = pr[k] * ip;
        return s_side;
    }
    double best_m[3] = {r[0] / d2, r[1] / d2, r[2] / d2};
    double best_H = Hof(best_m);
    double tilt_m[3] = {0, 0, 0}, tilt_H = std::numeric_limits<double>::infinity();   // start near an end face that failed
    auto consider = [&](const double* m) { const double H = Hof(m); if (H < best_H) { best_H = H; for (std::size_t k = 0; k < n; ++k) best_m[k] = m[k]; } return H; };
    // stop_apart: the caller only needs to know the pair is apart; any m with H(m) <= 1 proves it (s >= 1)
    auto remember = [&](int t) {   // what the next pass can start from
        if (!tag) return;
        *tag = t;
        if (nc && t >= 4) { const double nm = std::sqrt(dot(best_m, best_m, n)); for (std::size_t k = 0; k < n; ++k) nc[k] = best_m[k] / nm; }
    };
    #define ROD_APART_EXIT() if (stop_apart && best_H <= 1.0) { converged = true; remember(5); take_unit(best_m); return 1.0 / best_H; }
    // warm start: the normal the iteration ended on at the last pass. The closed-form candidates are skipped; if
    // the iteration does not confirm an optimum from there, everything is redone from scratch.
    bool warm = false;
    double warm_m[3] = {0, 0, 0};
    bool warm_onA = false, warm_onB = false;   // the cached normal sat on a kink (perpendicular to an axis)
    if (tag && nc && *tag >= 4) {
        if (*tag != 5) { warm_onA = ((*tag) & 8) != 0; warm_onB = ((*tag) & 16) != 0; }
        const double nr = dot(nc, r, n);
        if (nr > 1e-12 * d) { warm = true; for (std::size_t k = 0; k < n; ++k) warm_m[k] = nc[k] / nr; ++g_stats.warm; }
    }
    ROD_APART_EXIT()
    // closed-form candidates, the one that held at the last pass first: end face of i (0), of j (1), side against side (2)
    int order[3] = {0, 1, 2};
    if (tag && *tag == 2) { order[0] = 1; order[1] = 0; }
    else if (tag && *tag == 3) { order[0] = 2; order[1] = 0; order[2] = 1; }
    for (int oi = 0; oi < 3 && !warm; ++oi) {
        if (order[oi] < 2) {
        const int which = order[oi];
        const double* ua = which == 0 ? ui : uj;
        const double* ub = which == 0 ? uj : ui;
        const double hla = which == 0 ? hli : hlj, rada = which == 0 ? radi : radj;
        const double hlb = which == 0 ? hlj : hli, radb = which == 0 ? radj : radi;
        const double ar = dot(ua, r, n);
        if (std::abs(ar) < eps * d) continue;
        double m[3], nn[3];
        for (std::size_t k = 0; k < n; ++k) { m[k] = ua[k] / ar; nn[k] = sgn(ar) * ua[k]; }
        const double H = consider(m);
        ROD_APART_EXIT()
        // the other rod's support point along nn; the face's point is r H - that point; inside the end disc?
        const SupportSet B = support_set(n, ub, hlb, radb, nn);
        double y[3], qb[3], rho2 = 0.0;
        for (std::size_t k = 0; k < n; ++k) y[k] = r[k] * H - hla * nn[k];     // target minus the face centre
        B.project(y, qb);                                                   // the point of B's set closest to it
        const double along = dot(ua, y, n) - dot(ua, qb, n);
        for (std::size_t k = 0; k < n; ++k) { const double v = y[k] - qb[k] - along * ua[k]; rho2 += v * v; }
        if (rho2 <= rada * rada * (1.0 + 1e-10)) { converged = true; ++g_stats.closed_form; remember(1 + which); take_unit(m); return 1.0 / H; }
        // not the face but possibly its rim: a start with the normal tilted towards where the contact point fell
        const double irho = 1.0 / std::sqrt(rho2);
        double tm[3], tr = 0.0;
        for (std::size_t k = 0; k < n; ++k) { tm[k] = nn[k] + 0.05 * (y[k] - qb[k] - along * ua[k]) * irho; tr += tm[k] * r[k]; }
        if (tr > 0.0) {
            for (std::size_t k = 0; k < n; ++k) tm[k] /= tr;
            const double Ht = Hof(tm);
            if (Ht < tilt_H) { tilt_H = Ht; for (std::size_t k = 0; k < n; ++k) tilt_m[k] = tm[k]; }
        }
        continue;
        }
    // side against side
    {
        const double wr = dot(w, r, n);
        if (!warm && std::abs(wr) > eps * d * std::sqrt(w2)) {
            double m[3], nn[3];
            const double iw = 1.0 / std::sqrt(w2);
            for (std::size_t k = 0; k < n; ++k) { m[k] = w[k] / wr; nn[k] = sgn(wr) * w[k] * iw; }
            const double H = consider(m);
            ROD_APART_EXIT()
            double t[3];
            for (std::size_t k = 0; k < n; ++k) t[k] = r[k] * H - (radi + radj) * nn[k];
            const double b = dot(ui, uj, n), ta_ = dot(t, ui, n), tb_ = dot(t, uj, n);
            const double ta = (ta_ - b * tb_) / w2, tb = (tb_ - b * ta_) / w2;       // 1 - b^2 = |u_i x u_j|^2
            if (std::abs(ta) <= hli * (1.0 + 1e-10) && std::abs(tb) <= hlj * (1.0 + 1e-10)) { converged = true; ++g_stats.closed_form; remember(3); take_unit(m); return 1.0 / H; }
        }
    }
    }
    // Newton on the chart m = m0 + a e1 + b e2 (e1, e2 perpendicular to r)
    double e1[3], e2[3];
    {
        const double rh[3] = {r[0] / d, r[1] / d, r[2] / d};
        const double ax[3] = {std::abs(rh[0]) < 0.6 ? 1.0 : 0.0, std::abs(rh[0]) < 0.6 ? 0.0 : 1.0, 0.0};
        cross3(rh, ax, e1);
        const double i1 = 1.0 / std::sqrt(dot(e1, e1, n));
        for (std::size_t k = 0; k < n; ++k) e1[k] *= i1;
        cross3(rh, e1, e2);
    }
    double m[3] = {r[0] / d2, r[1] / d2, r[2] / d2};
    if (best_H < Hof(m) && std::abs(dot(best_m, ui, n)) < 1e-14 / d && std::abs(dot(best_m, uj, n)) < 1e-14 / d)
        for (std::size_t k = 0; k < n; ++k) m[k] = best_m[k];      // start from the side-side point when it is lower
    if (tilt_H < Hof(m)) for (std::size_t k = 0; k < n; ++k) m[k] = tilt_m[k];
    const double aA_[2] = {dot(ui, e1, n), dot(ui, e2, n)}, aB_[2] = {dot(uj, e1, n), dot(uj, e2, n)};
    if (warm) {
        for (std::size_t k = 0; k < n; ++k) m[k] = warm_m[k];
        // back onto the kink it was on (the axis has turned a little since): the iteration then starts along it
        for (int rep = 0; rep < 2; ++rep) {
            if (warm_onA) { const double c = dot(ui, m, n), a2 = aA_[0] * aA_[0] + aA_[1] * aA_[1]; if (a2 > 0.0) for (std::size_t k = 0; k < n; ++k) m[k] -= c / a2 * (aA_[0] * e1[k] + aA_[1] * e2[k]); }
            if (warm_onB) { const double c = dot(uj, m, n), a2 = aB_[0] * aB_[0] + aB_[1] * aB_[1]; if (a2 > 0.0) for (std::size_t k = 0; k < n; ++k) m[k] -= c / a2 * (aB_[0] * e1[k] + aB_[1] * e2[k]); }
            if (!(warm_onA && warm_onB)) break;
        }
    }
    double Hcur = Hof(m);
    if (Hcur < best_H) { best_H = Hcur; for (std::size_t k = 0; k < n; ++k) best_m[k] = m[k]; }
    const double aA[2] = {dot(ui, e1, n), dot(ui, e2, n)}, aB[2] = {dot(uj, e1, n), dot(uj, e2, n)};
    const double aA2 = aA[0] * aA[0] + aA[1] * aA[1], aB2 = aB[0] * aB[0] + aB[1] * aB[1];
    // put m back exactly on a kink it is meant to stay on (c = 0), moving inside the chart
    auto snap = [&](double* mm, bool toA, bool toB) {
        for (int rep = 0; rep < 2; ++rep) {
            if (toA && aA2 > 0.0) { const double c = dot(ui, mm, n); for (std::size_t k = 0; k < n; ++k) mm[k] -= c / aA2 * (aA[0] * e1[k] + aA[1] * e2[k]); }
            if (toB && aB2 > 0.0) { const double c = dot(uj, mm, n); for (std::size_t k = 0; k < n; ++k) mm[k] -= c / aB2 * (aB[0] * e1[k] + aB[1] * e2[k]); }
            if (!(toA && toB)) break;
        }
    };
    ++g_stats.iterated;
    // A warm start that is still good needs one or two steps. If it needs more, the contact has changed kind (for
    // instance to an end face, which the iteration can only creep towards): stop and redo from scratch, where the
    // closed-form candidates are tried first.
    const int iter_cap = warm ? std::min(max_iter, 4) : max_iter;
    for (int it = 0; it < iter_cap; ++it) {
        ++g_stats.iterations;
        const double mn = std::sqrt(dot(m, m, n));
        const double cA = dot(ui, m, n), cB = dot(uj, m, n);
        const bool onA = std::abs(cA) <= 1e-13 * mn, onB = std::abs(cB) <= 1e-13 * mn;
        double vA[3], vB[3], sA2 = 0.0, sB2 = 0.0;
        for (std::size_t k = 0; k < n; ++k) { vA[k] = m[k] - cA * ui[k]; vB[k] = m[k] - cB * uj[k]; sA2 += vA[k] * vA[k]; sB2 += vB[k] * vB[k]; }
        const double sA = std::sqrt(sA2), sB = std::sqrt(sB2);
        if (sA <= 1e-11 * mn || sB <= 1e-11 * mn) break;          // at an end-face normal: the face candidates cover it
        for (std::size_t k = 0; k < n; ++k) { vA[k] /= sA; vB[k] /= sB; }
        // gradient of the smooth terms, and of the |c| terms away from their kinks
        const double vA1 = dot(vA, e1, n), vA2 = dot(vA, e2, n), vB1 = dot(vB, e1, n), vB2 = dot(vB, e2, n);
        const double kA = onA ? 0.0 : hli * sgn(cA), kB = onB ? 0.0 : hlj * sgn(cB);
        const double gs[2] = {radi * vA1 + radj * vB1 + kA * aA[0] + kB * aB[0], radi * vA2 + radj * vB2 + kA * aA[1] + kB * aB[1]};
        // On a kink the slope of hl |c| may be anything in [-hl, hl] (its subgradient). The steepest descent
        // direction is minus the SMALLEST gradient over that set; a multiplier strictly inside its range means
        // "stay on this kink", one at an end means "leave on that side". Zero smallest gradient: optimal.
        // Hessian of the smooth terms in the chart
        const double wA = radi / sA, wB = radj / sB;   // e1, e2 orthonormal: e_a.e_b = delta_ab
        double h11 = wA * (1.0 - aA[0] * aA[0] - vA1 * vA1) + wB * (1.0 - aB[0] * aB[0] - vB1 * vB1);
        double h12 = wA * (-aA[0] * aA[1] - vA1 * vA2) + wB * (-aB[0] * aB[1] - vB1 * vB2);
        double h22 = wA * (1.0 - aA[1] * aA[1] - vA2 * vA2) + wB * (1.0 - aB[1] * aB[1] - vB2 * vB2);
        const double ridge = 1e-12 * (h11 + h22) + 1e-300;
        h11 += ridge; h22 += ridge;
        double g[2] = {gs[0], gs[1]}, lamA = 0.0, lamB = 0.0;
        bool fixA = false, fixB = false;
        // On ONE kink: first go to the best point along it (a smooth one-dimensional problem); only there is the
        // multiplier asked whether to leave. Leaving earlier, from a first-order slope, fails near an end-face
        // normal, where the other rod's term bends sharply.
        auto best1 = [&](const double* base, const double* a, double a2, double hl) {
            return a2 > 0.0 ? clampd(-(base[0] * a[0] + base[1] * a[1]) / a2, -hl, hl) : 0.0;
        };
        auto along_kink = [&](const double* a, double a2, double hl, double& lam) {
            const double tn = std::sqrt(a2);
            const double t[2] = {-a[1] / tn, a[0] / tn};
            const double gt = gs[0] * t[0] + gs[1] * t[1];
            const double ht = h11 * t[0] * t[0] + 2.0 * h12 * t[0] * t[1] + h22 * t[1] * t[1];
            if (gt * gt / ht > std::max(tol, 1e-13) * Hcur) {       // not yet at rest along the kink: stay on it
                g[0] = gt * t[0]; g[1] = gt * t[1];
                lam = 0.0;
                return true;
            }
            lam = best1(gs, a, a2, hl);
            g[0] = gs[0] + lam * a[0]; g[1] = gs[1] + lam * a[1];
            return std::abs(lam) < hl;
        };
        if (onA && onB) {
            double bestn = std::numeric_limits<double>::infinity();
            double cand[5][2];
            int nc = 0;
            const double det = aA[0] * aB[1] - aA[1] * aB[0];
            if (std::abs(det) > 1e-12 * std::sqrt(aA2 * aB2)) {
                const double la = (-gs[0] * aB[1] + gs[1] * aB[0]) / det, lb = (-aA[0] * gs[1] + aA[1] * gs[0]) / det;
                if (std::abs(la) <= hli && std::abs(lb) <= hlj) { cand[nc][0] = la; cand[nc][1] = lb; ++nc; }
            }
            for (int e = 0; e < 4; ++e) {      // one multiplier at an end, the other the best for it
                double base[2];
                if (e < 2) { const double la = e == 0 ? hli : -hli; base[0] = gs[0] + la * aA[0]; base[1] = gs[1] + la * aA[1]; cand[nc][0] = la; cand[nc][1] = best1(base, aB, aB2, hlj); }
                else { const double lb = e == 2 ? hlj : -hlj; base[0] = gs[0] + lb * aB[0]; base[1] = gs[1] + lb * aB[1]; cand[nc][1] = lb; cand[nc][0] = best1(base, aA, aA2, hli); }
                ++nc;
            }
            for (int c = 0; c < nc; ++c) {
                const double gx = gs[0] + cand[c][0] * aA[0] + cand[c][1] * aB[0], gy = gs[1] + cand[c][0] * aA[1] + cand[c][1] * aB[1];
                const double nrm = gx * gx + gy * gy;
                if (nrm < bestn) { bestn = nrm; lamA = cand[c][0]; lamB = cand[c][1]; g[0] = gx; g[1] = gy; }
            }
            fixA = std::abs(lamA) < hli; fixB = std::abs(lamB) < hlj;
        } else if (onA && aA2 > 0.0) {
            fixA = along_kink(aA, aA2, hli, lamA);
        } else if (onB && aB2 > 0.0) {
            fixB = along_kink(aB, aB2, hlj, lamB);
        }
        const double gnorm = std::sqrt(g[0] * g[0] + g[1] * g[1]);
        if (gnorm <= 1e-13 * Hcur * mn || (fixA && fixB)) { converged = true; break; }   // H ~ Hcur, lengths ~ 1/mn... slope scale Hcur * mn
        double dl[2] = {0.0, 0.0};
        if (fixA || fixB) {
            const double* a = fixA ? aA : aB;                        // stay on the kink: move along it
            double t[2] = {-a[1], a[0]};
            const double tn = std::sqrt(t[0] * t[0] + t[1] * t[1]);
            t[0] /= tn; t[1] /= tn;
            const double gt = g[0] * t[0] + g[1] * t[1];
            const double ht = h11 * t[0] * t[0] + 2.0 * h12 * t[0] * t[1] + h22 * t[1] * t[1];
            const double st = -gt / ht;
            dl[0] = st * t[0]; dl[1] = st * t[1];
        } else {
            const double det = h11 * h22 - h12 * h12;
            if (det > 1e-14 * h11 * h22) {
                dl[0] = -(h22 * g[0] - h12 * g[1]) / det;
                dl[1] = -(-h12 * g[0] + h11 * g[1]) / det;
            }
            // leaving a kink: the step must go to the side its multiplier chose; otherwise plain steepest descent
            const bool wrongA = onA && (aA[0] * dl[0] + aA[1] * dl[1]) * lamA <= 0.0;
            const bool wrongB = onB && (aB[0] * dl[0] + aB[1] * dl[1]) * lamB <= 0.0;
            if (!(det > 1e-14 * h11 * h22) || wrongA || wrongB) {
                const double sc = 1.0 / (h11 + h22);
                dl[0] = -sc * g[0]; dl[1] = -sc * g[1];
            }
        }
        double step[3];
        for (std::size_t k = 0; k < n; ++k) step[k] = dl[0] * e1[k] + dl[1] * e2[k];
        const double pred = -(g[0] * dl[0] + g[1] * dl[1]);          // twice the predicted decrease
        if (!(pred > 1e-16 * Hcur)) { converged = true; break; }
        // do not step across a kink: stop on it and decide there
        double tau = 1.0;
        int hit = 0;
        const double dA = dot(ui, step, n), dB = dot(uj, step, n);
        if (!onA && cA * (cA + dA) < 0.0) { const double tt = -cA / dA; if (tt < tau) { tau = tt; hit = 1; } }
        if (!onB && cB * (cB + dB) < 0.0) { const double tt = -cB / dB; if (tt < tau) { tau = tt; hit = 2; } }
        double mt[3], Hn = 0.0;
        bool ok = false;
        if (hit == 0 && pred <= std::max(tol, 1e-13) * Hcur) {
            // The step is within the tolerance (inside a run: the tracking tolerance; otherwise rounding): take
            // it as it is and stop. Checking that H went down would compare numbers that differ by rounding.
            for (std::size_t k = 0; k < n; ++k) mt[k] = m[k] + step[k];
            snap(mt, fixA, fixB);
            Hn = Hof(mt);
            if (Hn < best_H) { best_H = Hn; for (std::size_t k = 0; k < n; ++k) best_m[k] = mt[k]; }
            converged = true;
            ROD_APART_EXIT()
            break;
        }
        for (int bt = 0; bt < 8; ++bt) {
            for (std::size_t k = 0; k < n; ++k) mt[k] = m[k] + tau * step[k];
            snap(mt, fixA || (hit == 1 && bt == 0), fixB || (hit == 2 && bt == 0));
            Hn = Hof(mt);
            if (Hn <= Hcur * (1.0 + 1e-13)) { ok = true; break; }
            tau *= 0.5; ++g_stats.backtracks;
        }
        if (!ok) {
            ++g_stats.descents;
            // the quadratic model failed (close to an end-face normal the function is a cone): plain descent along
            // minus the gradient, shortened until it goes down
            double len = 1.0 / (h11 + h22);
            if (warm) break;
            for (int bt = 0; bt < 12 && !ok; ++bt, len *= 0.125) {
                for (std::size_t k = 0; k < n; ++k) mt[k] = m[k] - len * (g[0] * e1[k] + g[1] * e2[k]);
                snap(mt, fixA, fixB);
                const double cA2 = dot(ui, mt, n), cB2 = dot(uj, mt, n);
                if ((!onA && cA * cA2 < 0.0) || (!onB && cB * cB2 < 0.0)) continue;   // would cross a kink: shorter
                Hn = Hof(mt);
                if (Hn < Hcur) ok = true;
            }
            hit = 1;   // not a Newton step: do not read convergence from it
            if (!ok) break;
        }
        const double gain = Hcur - Hn;
        for (std::size_t k = 0; k < n; ++k) m[k] = mt[k];
        Hcur = Hn;
        if (Hcur < best_H) { best_H = Hcur; for (std::size_t k = 0; k < n; ++k) best_m[k] = m[k]; }
        ROD_APART_EXIT()
        if (hit == 0 && tau == 1.0 && gain <= 1e-15 * Hcur && pred <= 1e-13 * Hcur) { converged = true; break; }
        // tracking (tol > 0, inside a run): a full Newton step whose predicted gain was already below the tolerance
        // is accepted without another look; the next pass starts from it
        if (tol > 0.0 && hit == 0 && tau == 1.0 && (pred <= tol * Hcur || (warm && pred * pred <= tol * Hcur * Hcur))) { converged = true; break; }
    }
    if (warm && !converged) {   // the warm start did not confirm an optimum: from scratch
        ++g_stats.warm_failed;
        int t0 = 0;
        const double s0 = scale_uncapped_3d(r, ui, hli, radi, uj, hlj, radj, max_iter, nv, converged, stop_apart, &t0, nc, tol);
        *tag = t0;
        return s0;
    }
    {
        const double mn = std::sqrt(dot(best_m, best_m, n));
        remember(4 + (std::abs(dot(ui, best_m, n)) <= 1e-13 * mn ? 8 : 0) + (std::abs(dot(uj, best_m, n)) <= 1e-13 * mn ? 16 : 0));
    }
    take_unit(best_m);
    return 1.0 / best_H;
}

// 2D pair test for the force pass: false as soon as one of the four edge normals separates the rectangles.
ROD_INLINE bool uncapped_2d_contact(const double* r, const double* ui, double hli, double radi, const double* uj,
                                    double hlj, double radj, double& s, double* nv) {
    const double cs = std::abs(ui[0] * uj[0] + ui[1] * uj[1]), sn = std::abs(ui[0] * uj[1] - ui[1] * uj[0]);
    const double a = ui[0] * r[0] + ui[1] * r[1], b = -ui[1] * r[0] + ui[0] * r[1];
    const double c = uj[0] * r[0] + uj[1] * r[1], e = -uj[1] * r[0] + uj[0] * r[1];
    const double Ha = hli + hlj * cs + radj * sn, Hb = radi + hlj * sn + radj * cs;
    const double Hc = hlj + hli * cs + radi * sn, He = radj + hli * sn + radi * cs;
    const double aa = std::abs(a), ab = std::abs(b), ac = std::abs(c), ae = std::abs(e);
    if (aa >= Ha || ab >= Hb || ac >= Hc || ae >= He) return false;
    const double ra = aa / Ha, rb = ab / Hb, rc = ac / Hc, re = ae / He;
    s = ra; nv[0] = a < 0.0 ? -ui[0] : ui[0]; nv[1] = a < 0.0 ? -ui[1] : ui[1];
    if (rb > s) { s = rb; const double sg = b < 0.0 ? -1.0 : 1.0; nv[0] = -sg * ui[1]; nv[1] = sg * ui[0]; }
    if (rc > s) { s = rc; nv[0] = c < 0.0 ? -uj[0] : uj[0]; nv[1] = c < 0.0 ? -uj[1] : uj[1]; }
    if (re > s) { s = re; const double sg = e < 0.0 ? -1.0 : 1.0; nv[0] = -sg * uj[1]; nv[1] = sg * uj[0]; }
    return true;
}

// 3D tracking (inside a run, tol > 0): is the contact normal of the last pass still good enough, without solving?
// The normal nc is first put back on the kink it sat on (perpendicular to axis i or j: the axis has turned a
// little). With the support points p_i, p_j along it, s = n.r / (h_i + h_j) is a valid lower bound on the scale,
// and the normal is the maximiser exactly when p_i + p_j = r / s. It is accepted when the mismatch is below
// sqrt(tol) of |r / s| (on a kink: the mismatch across the side line, and the point within the rod's length).
// Fills nv, s and the two contact points. apart: the bound already shows the pair is apart.
ROD_INLINE bool tracked_normal_holds(const double* r, const double* ui, double hli, double radi, const double* uj,
                                     double hlj, double radj, const double* nc, bool onA, bool onB, double tol,
                                     double* nv, double& s, double* pi, double* qj, bool& apart) {
    apart = false;
    if (onA && onB) return false;
    double nn[3] = {nc[0], nc[1], nc[2]};
    const double* uk = onA ? ui : uj;             // the rod whose side line carries the contact (if any)
    if (onA || onB) {
        const double c = uk[0] * nn[0] + uk[1] * nn[1] + uk[2] * nn[2];
        for (int k = 0; k < 3; ++k) nn[k] -= c * uk[k];
        const double l2 = nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2];
        if (!(l2 > 1e-18)) return false;
        const double il = 1.0 / std::sqrt(l2);
        for (int k = 0; k < 3; ++k) nn[k] *= il;
    }
    // unique support point of a cylinder along nn (false when nn is along or perpendicular to its axis)
    auto point = [&](const double* u, double hl, double rad, double* p, double& h) {
        const double c = u[0] * nn[0] + u[1] * nn[1] + u[2] * nn[2];
        const double q[3] = {nn[0] - c * u[0], nn[1] - c * u[1], nn[2] - c * u[2]};
        const double sp = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
        if (sp <= 1e-9 || std::abs(c) <= 1e-9) return false;
        const double a = c < 0.0 ? -hl : hl, bq = rad / sp;
        for (int k = 0; k < 3; ++k) p[k] = a * u[k] + bq * q[k];
        h = hl * std::abs(c) + rad * sp;
        return true;
    };
    const double nr = nn[0] * r[0] + nn[1] * r[1] + nn[2] * r[2];
    if (!(nr > 0.0)) return false;
    double hA = radi, hB = radj;
    if (!onA && !point(ui, hli, radi, pi, hA)) return false;
    if (!onB && !point(uj, hlj, radj, qj, hB)) return false;
    const double H = hA + hB;
    if (nr >= H) { apart = true; nv[0] = nn[0]; nv[1] = nn[1]; nv[2] = nn[2]; s = nr / H; return true; }
    const double inv_s = H / nr;
    double y[3], y2 = 0.0, t2 = 0.0;
    if (onA || onB) {
        const double* pk = onA ? qj : pi;         // the other rod's (unique) point
        const double radk = onA ? radi : radj, hlk = onA ? hli : hlj;
        double along = 0.0;
        for (int k = 0; k < 3; ++k) { y[k] = r[k] * inv_s - pk[k] - radk * nn[k]; along += uk[k] * y[k]; }
        if (std::abs(along) > hlk) return false;
        for (int k = 0; k < 3; ++k) { const double v = y[k] - along * uk[k]; y2 += v * v; t2 += r[k] * r[k]; }
        if (!(y2 <= tol * t2 * inv_s * inv_s)) return false;
        double* po = onA ? pi : qj;
        for (int k = 0; k < 3; ++k) po[k] = along * uk[k] + radk * nn[k];
    } else {
        for (int k = 0; k < 3; ++k) { const double v = r[k] * inv_s - pi[k] - qj[k]; y2 += v * v; t2 += r[k] * r[k]; }
        if (!(y2 <= tol * t2 * inv_s * inv_s)) return false;
    }
    nv[0] = nn[0]; nv[1] = nn[1]; nv[2] = nn[2];
    s = nr / H;
    return true;
}

// 3D: exact separation by the centre line, the two axes and their cross product (cheap; most listed pairs).
ROD_INLINE bool uncapped_3d_separated(const double* r, double d2, const double* ui, double hli, double radi,
                                      const double* uj, double hlj, double radj) {
    double w[3];
    cross3(ui, uj, w);
    const double w2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];     // = 1 - (u_i.u_j)^2
    if (w2 > 1e-20) { const double wr = w[0] * r[0] + w[1] * r[1] + w[2] * r[2]; if (wr * wr >= w2 * (radi + radj) * (radi + radj)) return true; }   // normal along u_i x u_j
    const double b = std::abs(ui[0] * uj[0] + ui[1] * uj[1] + ui[2] * uj[2]);
    const double a = std::abs(ui[0] * r[0] + ui[1] * r[1] + ui[2] * r[2]), c = std::abs(uj[0] * r[0] + uj[1] * r[1] + uj[2] * r[2]);
    const double ea = a - hli - hlj * b, ec = c - hlj - hli * b;   // normal along an axis: a >= hl_i + hl_j b + R_j sn, sn = sqrt(w2)
    if (ea >= 0.0 && ea * ea >= radj * radj * w2) return true;
    if (ec >= 0.0 && ec * ec >= radi * radi * w2) return true;
    const double pa = std::sqrt(std::max(0.0, d2 - a * a)), pc = std::sqrt(std::max(0.0, d2 - c * c));
    return d2 >= hli * a + radi * pa + hlj * c + radj * pc;        // normal along r
}

double scale_uncapped(std::size_t n, const double* r, const double* ui, double hli, double radi, const double* uj,
                      double hlj, double radj, int max_iter, double* nv) {
    if (n == 2) return scale_uncapped_2d(r, ui, hli, radi, uj, hlj, radj, nv);
    bool conv = false;
    return scale_uncapped_3d(r, ui, hli, radi, uj, hlj, radj, max_iter, nv, conv);
}

}  // namespace

// --------------------------------------------------------------------------------------------------------------------

// ---- rounded rods ------------------------------------------------------------------------------------------------
namespace {
// Support function of the rounded rod along m (any length), its gradient (the support point, homogeneous of degree 0)
// and, when Hv is given, its Hessian applied to the vectors in V (nv of them, each of length n), ADDED to Hv.
ROD_INLINE double rounded_support(std::size_t n, int kind, double eps, const double* u, double hl, double rad,
                                  const double* m, double* grad, const double* V, std::size_t nvec, double* Hv) {
    const double e2 = eps * eps;
    const double c = dot(u, m, n), m2 = dot(m, m, n);
    double p2 = 0.0, pm[MD];
    for (std::size_t k = 0; k < n; ++k) { pm[k] = m[k] - c * u[k]; p2 += pm[k] * pm[k]; }
    const double S1 = std::sqrt(c * c + e2 * m2);
    // first term: hl sqrt(m.Q1 m), Q1 = u u^T + eps^2 I
    double q1[MD];
    for (std::size_t k = 0; k < n; ++k) q1[k] = c * u[k] + e2 * m[k];
    double h = hl * S1;
    if (grad) for (std::size_t k = 0; k < n; ++k) grad[k] = hl * q1[k] / S1;
    if (Hv) for (std::size_t a = 0; a < nvec; ++a) {
        const double* v = V + a * n;
        const double uv = dot(u, v, n), qv = dot(q1, v, n);
        for (std::size_t k = 0; k < n; ++k) Hv[a * n + k] += hl * (uv * u[k] + e2 * v[k] - q1[k] * qv / (S1 * S1)) / S1;
    }
    if (kind == SHAPE_ROD_CAPPED) {     // + R |m|
        const double mn = std::sqrt(m2);
        h += rad * mn;
        if (grad) for (std::size_t k = 0; k < n; ++k) grad[k] += rad * m[k] / mn;
        if (Hv) for (std::size_t a = 0; a < nvec; ++a) {
            const double* v = V + a * n;
            const double mv = dot(m, v, n);
            for (std::size_t k = 0; k < n; ++k) Hv[a * n + k] += rad * (v[k] - m[k] * mv / m2) / mn;
        }
        return h;
    }
    // + R sqrt(m.Q2 m), Q2 = (1 + eps^2) I - u u^T
    const double S2 = std::sqrt(p2 + e2 * m2);
    double q2[MD];
    for (std::size_t k = 0; k < n; ++k) q2[k] = pm[k] + e2 * m[k];
    h += rad * S2;
    if (grad) for (std::size_t k = 0; k < n; ++k) grad[k] += rad * q2[k] / S2;
    if (Hv) for (std::size_t a = 0; a < nvec; ++a) {
        const double* v = V + a * n;
        const double uv = dot(u, v, n), qv = dot(q2, v, n);
        for (std::size_t k = 0; k < n; ++k) Hv[a * n + k] += rad * ((1.0 + e2) * v[k] - uv * u[k] - q2[k] * qv / (S2 * S2)) / S2;
    }
    return h;
}
}  // namespace

bool rod_pair_rounded(std::size_t n, int kind, double eps, const double* r, const double* ui, double hli,
                      double radi, const double* uj, double hlj, double radj, double K, RodPairOut& out, double* Ti,
                      double* Tj) {
    if (n != 2 && n != 3) throw std::runtime_error("rods: the rounded contact (quench) is implemented in 2D and 3D");
    const double d2 = dot(r, r, n), d = std::sqrt(d2);
    const std::size_t nc = n - 1;                 // chart: m = m0 + sum_k a_k e_k, e_k perpendicular to r
    double E[2 * MD];
    if (n == 2) { E[0] = -r[1] / d; E[1] = r[0] / d; }
    else {
        const double rh[3] = {r[0] / d, r[1] / d, r[2] / d};
        const double ax[3] = {std::abs(rh[0]) < 0.6 ? 1.0 : 0.0, std::abs(rh[0]) < 0.6 ? 0.0 : 1.0, 0.0};
        cross3(rh, ax, E);
        const double i1 = 1.0 / std::sqrt(dot(E, E, 3));
        for (int k = 0; k < 3; ++k) E[k] *= i1;
        cross3(rh, E, E + 3);
    }
    double m[MD];
    for (std::size_t k = 0; k < n; ++k) m[k] = r[k] / d2;
    auto Hof = [&](const double* mm) { return rounded_support(n, kind, eps, ui, hli, radi, mm, nullptr, nullptr, 0, nullptr) +
                                              rounded_support(n, kind, eps, uj, hlj, radj, mm, nullptr, nullptr, 0, nullptr); };
    double Hcur = Hof(m);
    for (int it = 0; it < 100; ++it) {
        double g3[MD], Hv[2 * MD] = {0.0}, gb[MD];
        rounded_support(n, kind, eps, ui, hli, radi, m, g3, E, nc, Hv);
        rounded_support(n, kind, eps, uj, hlj, radj, m, gb, E, nc, Hv);
        for (std::size_t k = 0; k < n; ++k) g3[k] += gb[k];
        double step[MD] = {0.0}, pred;
        if (nc == 1) {
            const double g = dot(g3, E, n), hh = dot(Hv, E, n);
            if (!(hh > 0.0)) break;
            const double a = -g / hh;
            for (std::size_t k = 0; k < n; ++k) step[k] = a * E[k];
            pred = g * g / hh;
        } else {
            const double g0 = dot(g3, E, 3), g1 = dot(g3, E + 3, 3);
            const double h00 = dot(Hv, E, 3), h01 = dot(Hv, E + 3, 3), h11 = dot(Hv + 3, E + 3, 3);
            const double det = h00 * h11 - h01 * h01;
            if (!(det > 0.0)) break;
            const double a0 = -(h11 * g0 - h01 * g1) / det, a1 = -(-h01 * g0 + h00 * g1) / det;
            for (int k = 0; k < 3; ++k) step[k] = a0 * E[k] + a1 * E[3 + k];
            pred = -(g0 * a0 + g1 * a1);
        }
        if (!(pred > 1e-30 * Hcur)) break;
        double tau = 1.0, mt[MD], Hn = Hcur;
        bool ok = false;
        // Sufficient decrease (Armijo, 0.1 of the predicted fall), not just any decrease: the support function has a
        // crease of width eps at n.u = 0, and a full Newton step taken from outside it lands on the far side. With
        // "any decrease" the iteration then zig-zagged across the crease, a quarter step at a time, and left after
        // 100 iterations with H 4 % too large, i.e. a pair energy 4 times too large that came and went with a
        // displacement of 1e-7 (two crossed capped rods; EllipsoidalPacking LOG [033]). The rounding slack 2e-16
        // keeps the last step at the minimum acceptable.
        for (int bt = 0; bt < 40; ++bt, tau *= 0.5) {
            for (std::size_t k = 0; k < n; ++k) mt[k] = m[k] + tau * step[k];
            Hn = Hof(mt);
            if (Hn <= Hcur * (1.0 + 2e-16) - 0.1 * tau * pred) { ok = true; break; }
        }
        if (!ok) break;
        const bool done = tau == 1.0 && pred <= 1e-15 * Hcur;
        for (std::size_t k = 0; k < n; ++k) m[k] = mt[k];
        Hcur = Hn;
        if (done) break;
    }
    const double s = 1.0 / Hcur;
    if (!(s < 1.0)) return false;
    const double mn = std::sqrt(dot(m, m, n));
    double nv[MD], pi[MD], qj[MD];
    for (std::size_t k = 0; k < n; ++k) nv[k] = m[k] / mn;
    rounded_support(n, kind, eps, ui, hli, radi, nv, pi, nullptr, 0, nullptr);     // unique support points
    rounded_support(n, kind, eps, uj, hlj, radj, nv, qj, nullptr, 0, nullptr);
    const double nr = dot(nv, r, n), H = nr / s, sigma = d / s, pen = sigma - d;
    out.s = s; out.overlap = 1.0 - s; out.sigma = sigma; out.energy = 0.5 * K * pen * pen;
    out.F_mag = -K * (sigma / d - 1.0); out.dkappa = K * sigma * (d - sigma); out.converged = true;
    double f2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        out.dEdr[k] = K * pen * (r[k] / (d * s) - d * nv[k] / (s * nr) - r[k] / d);
        out.nvec[k] = nv[k]; out.pi[k] = pi[k]; out.pj[k] = qj[k];
        f2 += out.dEdr[k] * out.dEdr[k];
    }
    out.f_norm = std::sqrt(f2);
    const double cT = -K * pen * d / (s * H);
    if (Ti) add_torque(n, nv, pi, cT, Ti);
    if (Tj) add_torque(n, nv, qj, cT, Tj);
    return true;
}

RodKernelStats rod_kernel_stats(bool reset) {
    const RodKernelStats out = g_stats;
    if (reset) g_stats = RodKernelStats();
    return out;
}

double rod_row_volume_factor(std::size_t Ndim, double q0, double q1, bool capped) {
    const double cn = ball_coeff(Ndim), cn1 = ball_coeff(Ndim - 1);
    double v = 2.0 * q0 * cn1 * std::pow(q1, static_cast<double>(Ndim - 1));
    if (capped) v += cn * std::pow(q1, static_cast<double>(Ndim));
    return v / (cn * std::pow(0.5, static_cast<double>(Ndim)));
}

void rod_build_shape_table(const std::vector<double>& lengths, const std::vector<double>& radii, bool capped,
                           std::size_t N, std::size_t Ndim, std::vector<double>& qtab, std::vector<double>& deq,
                           std::vector<double>& bfac) {
    if (Ndim < 2 || Ndim > ELL_MAX_DIM) throw std::runtime_error("rods: Ndim must be between 2 and " + std::to_string(ELL_MAX_DIM));
    if (!capped && Ndim > 3) {
        throw std::runtime_error("rods: uncapped (flat-ended) rods are implemented in 2D and 3D only; capped rods "
                                 "work in any dimension");
    }
    if (lengths.empty() || radii.empty()) {
        throw std::runtime_error("rods: supply both rod_length and rod_radius (each one number or one per particle)");
    }
    if ((lengths.size() != 1 && lengths.size() != N) || (radii.size() != 1 && radii.size() != N)) {
        throw std::runtime_error("rods: rod_length and rod_radius each need one number or one per particle (" +
                                 std::to_string(N) + "); got " + std::to_string(lengths.size()) + " and " +
                                 std::to_string(radii.size()));
    }
    qtab.assign(N * Ndim, 1.0);
    deq.assign(N, 1.0);
    bfac.assign(N, 1.0);
    for (std::size_t i = 0; i < N; ++i) {
        const double L = lengths.size() == 1 ? lengths[0] : lengths[i];
        const double R = radii.size() == 1 ? radii[0] : radii[i];
        if (!(R > 0.0) || !std::isfinite(R)) throw std::runtime_error("rods: rod_radius must be a positive finite number (particle " + std::to_string(i) + ")");
        if (!std::isfinite(L) || L < 0.0 || (!capped && !(L > 0.0))) {
            throw std::runtime_error(std::string("rods: rod_length must be ") + (capped ? "zero or positive" : "positive (an uncapped rod of length 0 has no volume)") +
                                     " and finite (particle " + std::to_string(i) + ")");
        }
        const double b = capped ? L + 2.0 * R : 2.0 * std::sqrt(0.25 * L * L + R * R);
        bfac[i] = b;
        qtab[i * Ndim + 0] = 0.5 * L / b;
        qtab[i * Ndim + 1] = R / b;
        deq[i] = std::pow(rod_row_volume_factor(Ndim, qtab[i * Ndim], qtab[i * Ndim + 1], capped), 1.0 / static_cast<double>(Ndim));
    }
}

double rod_volume_factor(std::size_t N, std::size_t Ndim, const double* Db, const double* qtab, bool capped) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double v = 1.0;
        for (std::size_t k = 0; k < Ndim; ++k) v *= Db[i];
        num += v * rod_row_volume_factor(Ndim, qtab[i * Ndim], qtab[i * Ndim + 1], capped);
        den += v;
    }
    return den > 0.0 ? num / den : 1.0;
}

double rod_support(std::size_t Ndim, int kind, const double* u, double hl, double rad, const double* nv) {
    if (kind == SHAPE_ROD_CAPPED) return hl * std::abs(dot(u, nv, Ndim)) + rad * std::sqrt(dot(nv, nv, Ndim));
    return cyl_support(Ndim, u, hl, rad, nv);
}

double rod_half_width(std::size_t Ndim, int kind, const double* R, double hl, double rad, std::size_t d) {
    const double ud = R[d * Ndim + 0];
    if (kind == SHAPE_ROD_CAPPED) return hl * std::abs(ud) + rad;
    double p2 = 0.0;   // |e_d - u_d u|^2 from the components (see cyl_support)
    for (std::size_t k = 0; k < Ndim; ++k) { const double v = (k == d ? 1.0 : 0.0) - ud * R[k * Ndim + 0]; p2 += v * v; }
    return hl * std::abs(ud) + rad * std::sqrt(p2);
}

namespace {
ROD_INLINE bool rod_pair_impl(std::size_t n, int kind, const double* r, double d2, const double* ui, double hli,
                              double radi, const double* uj, double hlj, double radj, double K, int max_iter,
                              RodPairOut& out, double* Ti, double* Tj, int* tag = nullptr, double* nc = nullptr,
                              double tol = 0.0);
}
bool rod_pair(std::size_t n, int kind, const double* r, const double* ui, double hli, double radi,
              const double* uj, double hlj, double radj, double K, int max_iter, RodPairOut& out, double* Ti,
              double* Tj) {
    return rod_pair_impl(n, kind, r, dot(r, r, n), ui, hli, radi, uj, hlj, radj, K, max_iter, out, Ti, Tj);
}
namespace {
ROD_INLINE bool rod_pair_impl(std::size_t n, int kind, const double* r, double d2, const double* ui, double hli,
                              double radi, const double* uj, double hlj, double radj, double K, int max_iter,
                              RodPairOut& out, double* Ti, double* Tj, int* tag, double* nc, double tol) {
    const double d = std::sqrt(d2);
    if (kind == SHAPE_ROD_CAPPED) {
        double ti, tj;
        const double ds2 = segment_closest(n, r, ui, hli, uj, hlj, ti, tj);
        const double rho = radi + radj;
        if (!(ds2 < rho * rho)) return false;
        double ds = std::sqrt(ds2);
        double nh[MD];
        bool crossing = false;
        if (n == 2 && hli > 0.0 && hlj > 0.0) {   // do the two axis segments cross? each has its ends on both sides of the other
            const double ci = ui[0] * r[1] - ui[1] * r[0], cij = ui[0] * uj[1] - ui[1] * uj[0];
            const double cj = uj[0] * (-r[1]) - uj[1] * (-r[0]);
            crossing = (ci + hlj * cij) * (ci - hlj * cij) < 0.0 && (cj - hli * cij) * (cj + hli * cij) < 0.0;
        }
        if (crossing) {
            // 2D only: the two axis segments CROSS. Their distance is zero over a whole range of positions, which
            // would leave the energy flat there. The distance is continued to negative values: minus the smallest
            // displacement that uncrosses them, i.e. minus the smallest distance from an end of one segment to
            // the other segment. Same formulas with the signed distance and the normal pointing the other way.
            ds = signed_crossing_2d(r, ui, hli, uj, hlj, ti, tj, nh);
        } else
        if (ds > 0.0) { for (std::size_t k = 0; k < n; ++k) nh[k] = (r[k] + tj * uj[k] - ti * ui[k]) / ds; }
        else if (d > 0.0) { for (std::size_t k = 0; k < n; ++k) nh[k] = r[k] / d; }
        else { for (std::size_t k = 0; k < n; ++k) nh[k] = k == 0 ? 1.0 : 0.0; }
        const double pen = rho - ds;
        out.s = ds > 0.0 ? ds / rho : 0.0;
        out.overlap = 1.0 - ds / rho;
        out.energy = 0.5 * K * pen * pen;
        out.f_norm = K * pen;
        out.sigma = d + pen;
        out.F_mag = d > 0.0 ? -K * (out.sigma / d - 1.0) : -K;
        out.dkappa = -K * pen * (pen + dot(nh, r, n));
        double pi[MD], pj[MD];
        for (std::size_t k = 0; k < n; ++k) { out.dEdr[k] = -K * pen * nh[k]; out.nvec[k] = nh[k]; pi[k] = ti * ui[k]; pj[k] = tj * uj[k]; out.pi[k] = pi[k]; out.pj[k] = pj[k]; }
        if (Ti) add_torque(n, out.dEdr, pi, 1.0, Ti);       // force on i = +dEdr at p_i
        if (Tj) add_torque(n, out.dEdr, pj, -1.0, Tj);      // force on j = -dEdr at p_j
        return true;
    }
    double nv[MD];
    bool conv = true, have_points = false;
    double s, track_pi[3] = {0, 0, 0}, track_qj[3] = {0, 0, 0};
    if (n == 2) {
        if (!uncapped_2d_contact(r, ui, hli, radi, uj, hlj, radj, s, nv)) return false;
    } else {
        ++g_stats.pairs;
        if (tag && *tag == 5) {   // the normal that separated the pair last time: does it still?
            const double nr = dot(nc, r, n);
            if (nr >= cyl_support(n, ui, hli, radi, nc) + cyl_support(n, uj, hlj, radj, nc)) { ++g_stats.quick_apart; ++g_stats.sep_cached_ok; return false; }
            ++g_stats.sep_cached_failed;
        }
        if (uncapped_3d_separated(r, d2, ui, hli, radi, uj, hlj, radj)) { ++g_stats.quick_apart; if (tag) *tag = 0; return false; }
        if (tag && tol > 0.0 && *tag >= 4 && *tag != 5) {
            bool apart = false;
            if (tracked_normal_holds(r, ui, hli, radi, uj, hlj, radj, nc, ((*tag) & 8) != 0, ((*tag) & 16) != 0, tol, nv, s,
                                     track_pi, track_qj, apart)) {
                for (int k = 0; k < 3; ++k) nc[k] = nv[k];
                if (apart) { *tag = 5; ++g_stats.quick_apart; return false; }
                have_points = true;
                ++g_stats.tracked_ok;
            } else ++g_stats.tracked_rejected;
        }
        if (!have_points) s = scale_uncapped_3d(r, ui, hli, radi, uj, hlj, radj, max_iter, nv, conv, true, tag, nc, tol);
    }
    out.converged = conv;
    if (!(s < 1.0)) return false;
    if (n == 3) ++g_stats.contacts;
    const double nr = dot(nv, r, n);
    const double H = nr / s;                                // h_i(n) + h_j(n) for the unit normal
    const double sigma = d / s, pen = sigma - d;
    out.s = s;
    out.overlap = 1.0 - s;
    out.sigma = sigma;
    out.energy = 0.5 * K * pen * pen;
    out.F_mag = -K * (sigma / d - 1.0);
    out.dkappa = K * sigma * (d - sigma);
    // dsigma/dr = r/(d s) - d n / (s n.r);  dE/dr = K pen (dsigma/dr - r/d)
    double f2 = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double ds_dr = r[k] / (d * s) - d * nv[k] / (s * nr);
        out.dEdr[k] = K * pen * (ds_dr - r[k] / d);
        out.nvec[k] = nv[k];
        f2 += out.dEdr[k] * out.dEdr[k];
    }
    out.f_norm = std::sqrt(f2);
    // rotations: dH = n.(Omega p) for the support point p, dsigma = d/(s H) dH, T_ab = -K pen d/(s H) (n_a p_b - n_b p_a)
    double target[MD], pi[MD], qj[MD];
    for (std::size_t k = 0; k < n; ++k) target[k] = r[k] / s;
    if (have_points) { for (int k = 0; k < 3; ++k) { pi[k] = track_pi[k]; qj[k] = track_qj[k]; } }
    else contact_points(n, ui, hli, radi, uj, hlj, radj, nv, target, pi, qj);
    for (std::size_t k = 0; k < n; ++k) { out.pi[k] = pi[k]; out.pj[k] = qj[k]; }
    const double cT = -K * pen * d / (s * H);
    if (Ti) add_torque(n, nv, pi, cT, Ti);
    if (Tj) add_torque(n, nv, qj, cT, Tj);
    return true;
}

// The pair loop of the force pass, compiled once per dimension (NN = 0: any dimension).
struct RowAcc { double U = 0.0, L = 0.0, dk = 0.0, Fm = 0.0, E = 0.0, mx = 0.0; std::size_t count = 0; };
template <std::size_t NN>
void rod_rows(std::size_t n_runtime, int kind, int max_iter, double K, std::size_t N, const std::uint32_t* pairs_data,
              const std::uint64_t* pair_offsets, const double* x, const double* Db, const double* box,
              const double* half_box, const double* axis, const double* hl, const double* rad, double* F, double* Tq,
              std::size_t* z, double* gsum, RowAcc& acc, RodPairCache* cache, double tol, std::size_t i_begin,
              std::size_t i_end) {
    (void)N;
    const std::size_t n = NN ? NN : n_runtime, nrot = n * (n - 1) / 2;
    constexpr std::size_t CAP = RodPairCache::CAP;
    const bool use_cache = cache && n == 3 && kind == SHAPE_ROD_UNCAPPED;
    RodPairOut po;
    std::vector<std::uint32_t> keep;
    for (std::size_t i = i_begin; i < i_end; ++i) {
        const std::uint32_t* row = pairs_data + pair_offsets[i];
        const std::size_t cnt = row[0];
        std::uint32_t new_ids[CAP], new_n = 0;
        std::int8_t new_tag[CAP];
        double new_nv[CAP * 3];
        const std::uint32_t old_n = use_cache ? cache->count[i] : 0;
        const std::uint32_t* old_ids = use_cache ? cache->ids.data() + i * CAP : nullptr;
        std::uint32_t po_ = 0;                      // walks the old row: both rows are in ascending j
        const double* xi = x + i * n;
        const double Di = Db[i];
        // Pass 1, branch-free: keep the entries whose bounding spheres overlap (most do not).
        keep.resize(cnt + 1);
        std::uint32_t* keep_p = keep.data();
        std::uint32_t n_keep = 0;
        for (std::size_t q = 1; q <= cnt; ++q) {
            const std::uint32_t j = row[q];
            const double rb = 0.5 * (Di + Db[j]);
            double d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = x[static_cast<std::size_t>(j) * n + d] - xi[d];
                delta -= (delta > half_box[d] ? box[d] : 0.0);
                delta += (delta < -half_box[d] ? box[d] : 0.0);
                d2 += delta * delta;
            }
            keep_p[n_keep] = j;
            n_keep += (d2 < rb * rb) ? 1u : 0u;
        }
        // Pass 2: the pair kernel on those.
        for (std::uint32_t kk = 0; kk < n_keep; ++kk) {
            const std::size_t j = keep_p[kk];
            double r[MD], d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = x[j * n + d] - xi[d];
                delta -= (delta > half_box[d] ? box[d] : 0.0);
                delta += (delta < -half_box[d] ? box[d] : 0.0);
                r[d] = delta; d2 += delta * delta;
            }
            int tag = 0;
            double nc[3] = {0, 0, 0};
            if (use_cache) {
                while (po_ < old_n && old_ids[po_] < j) ++po_;
                if (po_ < old_n && old_ids[po_] == j) { tag = cache->tag[i * CAP + po_]; for (int k = 0; k < 3; ++k) nc[k] = cache->nvec[(i * CAP + po_) * 3 + k]; }
            }
            const bool contact = rod_pair_impl(n, kind, r, d2, axis + i * n, hl[i], rad[i], axis + j * n, hl[j], rad[j], K,
                                               max_iter, po, Tq + i * nrot, Tq + j * nrot, use_cache ? &tag : nullptr, nc,
                                               use_cache ? tol : 0.0);
            if (use_cache && new_n < CAP && tag >= 1) {
                new_ids[new_n] = static_cast<std::uint32_t>(j); new_tag[new_n] = static_cast<std::int8_t>(tag);
                for (int k = 0; k < 3; ++k) new_nv[new_n * 3 + k] = nc[k];
                ++new_n;
            }
            if (!contact) continue;
            for (std::size_t d = 0; d < n; ++d) { F[i * n + d] += po.dEdr[d]; F[j * n + d] -= po.dEdr[d]; }
            z[i]++; z[j]++;
            if (gsum) { gsum[i] += po.f_norm; gsum[j] += po.f_norm; }
            acc.U += po.overlap; ++acc.count; acc.L += po.sigma;
            if (po.overlap > acc.mx) acc.mx = po.overlap;
            acc.dk += po.dkappa; acc.Fm += po.F_mag; acc.E += po.energy;
        }
        if (use_cache) {
            for (std::uint32_t c = 0; c < new_n; ++c) {
                cache->ids[i * CAP + c] = new_ids[c]; cache->tag[i * CAP + c] = new_tag[c];
                for (int k = 0; k < 3; ++k) cache->nvec[(i * CAP + c) * 3 + k] = new_nv[c * 3 + k];
            }
            cache->count[i] = new_n;
        }
    }
}
}  // namespace

double rod_pair_normal(std::size_t n, int kind, const double* r, const double* ui, double hli, double radi,
                       const double* uj, double hlj, double radj, int max_iter, double* nv) {
    if (kind == SHAPE_ROD_CAPPED) {
        double ti, tj;
        const double ds = std::sqrt(segment_closest(n, r, ui, hli, uj, hlj, ti, tj));
        const double d = std::sqrt(dot(r, r, n));
        for (std::size_t k = 0; k < n; ++k) nv[k] = ds > 0.0 ? (r[k] + tj * uj[k] - ti * ui[k]) / ds : (d > 0.0 ? r[k] / d : (k == 0 ? 1.0 : 0.0));
        return ds / (radi + radj);
    }
    return scale_uncapped(n, r, ui, hli, radi, uj, hlj, radj, max_iter, nv);
}

double rod_pair_scale_exact(std::size_t n, int kind, const double* r, const double* ui, double hli, double radi,
                            const double* uj, double hlj, double radj) {
    if (kind == SHAPE_ROD_CAPPED) {
        const double rho = radi + radj;
        double ti, tj;
        auto f = [&](double lam) { return std::sqrt(segment_closest(n, r, ui, lam * hli, uj, lam * hlj, ti, tj)) - lam * rho; };
        double lo = 0.0, hi = 1.0;
        while (f(hi) > 0.0 && hi < 1e6) { lo = hi; hi *= 2.0; }     // apart at scale 1: the scale lies above 1
        for (int it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (f(mid) > 0.0) lo = mid; else hi = mid;
            if (hi - lo <= 4e-16 * hi) break;
        }
        return 0.5 * (lo + hi);
    }
    double nv[MD];
    if (n == 2) return scale_uncapped_2d(r, ui, hli, radi, uj, hlj, radj, nv);
    bool conv = false;
    double best = scale_uncapped_3d(r, ui, hli, radi, uj, hlj, radj, 400, nv, conv);
    // Always finished by hand, converged flag or not: near a cone tip the iteration can stop on a tiny predicted
    // gain without being at the optimum.
    // The iteration may not have confirmed an optimum (two nearly parallel flat ends, a flat end resting on a side: the
    // function is close to a cone tip there). Reporting only: finish with a pattern search over directions around
    // the normal reached, which needs no derivative; every direction tried is a valid lower bound on the scale.
    auto ratio = [&](const double* nn) { return dot(nn, r, n) / (cyl_support(n, ui, hli, radi, nn) + cyl_support(n, uj, hlj, radj, nn)); };
    {   // the natural normals first: along either axis, across either axis towards the other rod, along r
        const double d = std::sqrt(dot(r, r, 3));
        double cand[5][3];
        int nc = 0;
        for (int which = 0; which < 2; ++which) {
            const double* u = which == 0 ? ui : uj;
            const double c = dot(u, r, 3);
            if (c != 0.0) { for (int k = 0; k < 3; ++k) cand[nc][k] = (c < 0.0 ? -u[k] : u[k]); ++nc; }
            double q[3], q2 = 0.0;
            for (int k = 0; k < 3; ++k) { q[k] = r[k] - c * u[k]; q2 += q[k] * q[k]; }
            if (q2 > 1e-24 * d * d) { const double iq = 1.0 / std::sqrt(q2); for (int k = 0; k < 3; ++k) cand[nc][k] = q[k] * iq; ++nc; }
        }
        for (int k = 0; k < 3; ++k) cand[nc][k] = r[k] / d;
        ++nc;
        for (int c = 0; c < nc; ++c) { const double v = ratio(cand[c]); if (v > best) { best = v; for (int k = 0; k < 3; ++k) nv[k] = cand[c][k]; } }
    }
    double step = 0.05;
    for (int sweep = 0; sweep < 20000 && step > 1e-16; ++sweep) {
        double t1[3], t2[3];
        const double ax[3] = {std::abs(nv[0]) < 0.6 ? 1.0 : 0.0, std::abs(nv[0]) < 0.6 ? 0.0 : 1.0, 0.0};
        cross3(nv, ax, t1);
        const double i1 = 1.0 / std::sqrt(dot(t1, t1, 3));
        for (int k = 0; k < 3; ++k) t1[k] *= i1;
        cross3(nv, t1, t2);
        bool moved = false;
        for (int dir = 0; dir < 12 && !moved; ++dir) {      // the set of directions turns from sweep to sweep (ridges)
            const double ang = M_PI / 6.0 * dir + 2.399963229728653 * sweep;
            const double a = std::cos(ang), b = std::sin(ang);
            double tn[3], l2 = 0.0;
            for (int k = 0; k < 3; ++k) { tn[k] = nv[k] + step * (a * t1[k] + b * t2[k]); l2 += tn[k] * tn[k]; }
            const double il = 1.0 / std::sqrt(l2);
            for (int k = 0; k < 3; ++k) tn[k] *= il;
            const double v = ratio(tn);
            if (v > best) { best = v; for (int k = 0; k < 3; ++k) nv[k] = tn[k]; moved = true; }
        }
        if (!moved) step *= 0.5;
    }
    return best;
}

void rod_forces(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                std::size_t N, const EllShape& sh, const std::vector<double>& Db,
                const std::vector<double>& box, const std::vector<std::int8_t>& walls, const double* Rot,
                double* F, double* Tq, double& U, double& max_min_dist, double& Lc, double& Fmean, double mu,
                double& dkappa, std::vector<std::size_t>& z, double& energy, std::vector<double>& Amat_work,
                double* gsum, RodPairCache* cache) {
    const double K = 1.0;
    const std::size_t n = sh.Ndim, nn = n * n, nrot = ell_nrot(n);
    const int kind = sh.kind;
    for (std::size_t k = 0; k < N * n; ++k) F[k] = 0.0;
    for (std::size_t k = 0; k < N * nrot; ++k) Tq[k] = 0.0;
    if (z.size() == N) std::fill(z.begin(), z.end(), 0); else z.assign(N, 0);
    if (gsum) for (std::size_t k = 0; k < N; ++k) gsum[k] = 0.0;
    // half-widths along the lab axes: read by the wall contacts below and by the quench (which asks for gsum)
    bool need_h = gsum != nullptr;
    for (std::size_t d = 0; d < n; ++d) if (walls[d] == 1) need_h = true;
    if (need_h || Amat_work.size() != N * nn) Amat_work.assign(N * nn, 0.0);
    // axes (column 0 of each rotation), half-lengths and radii
    std::vector<double> axis(N * n), hl(N), rad(N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t k = 0; k < n; ++k) axis[i * n + k] = Rot[i * nn + k * n];
        hl[i] = Db[i] * sh.qtab[i * n];
        rad[i] = Db[i] * sh.qtab[i * n + 1];
        if (need_h) for (std::size_t d = 0; d < n; ++d) {
            double h;
            if (sh.rod_round > 0.0) { double e[MD]; for (std::size_t k = 0; k < n; ++k) e[k] = k == d ? 1.0 : 0.0; h = rounded_support(n, kind, sh.rod_round, axis.data() + i * n, hl[i], rad[i], e, nullptr, nullptr, 0, nullptr); }
            else h = rod_half_width(n, kind, Rot + i * nn, hl[i], rad[i], d);
            Amat_work[i * nn + d * n + d] = h * h;
        }
    }
    double half_box[MD];
    for (std::size_t d = 0; d < n; ++d) half_box[d] = walls[d] == 0 ? 0.5 * box[d] : std::numeric_limits<double>::infinity();
    const double eps = sh.rod_round;
    // refused HERE, before the parallel region: an exception must not leave a thread of the region (it would end the process)
    if (eps > 0.0 && n != 2 && n != 3) throw std::runtime_error("rods: the rounded contact (quench) is implemented in 2D and 3D");
    // Rows [i_begin, i_end) of the half neighbour list on each thread, as ell_forces does it: thread 0 writes straight
    // into the outputs, the others into their own buffers, which are added in thread order afterwards. With one
    // thread nothing is buffered and the sums run in the order they always did. The cache rows of the 3D uncapped
    // kernel are per particle i, so each is written by the one thread that owns row i.
    const int nthreads = omp_get_max_threads();
    const std::size_t F_size = N * n, T_size = N * nrot;
    std::vector<double> Fbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * F_size : 0, 0.0);
    std::vector<double> Tbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * T_size : 0, 0.0);
    std::vector<std::size_t> zbuf(nthreads > 1 ? static_cast<std::size_t>(nthreads - 1) * N : 0, 0);
    std::vector<double> gbuf((gsum && nthreads > 1) ? static_cast<std::size_t>(nthreads - 1) * N : 0, 0.0);
    std::vector<RowAcc> accs(static_cast<std::size_t>(nthreads));
    const double* Dbd = Db.data();
    const double* boxd = box.data();
    #pragma omp parallel num_threads(nthreads)
    {
        const int tid = omp_get_thread_num();
        const int nth = omp_get_num_threads();
        double* F_t = tid == 0 ? F : Fbuf.data() + static_cast<std::size_t>(tid - 1) * F_size;
        double* T_t = tid == 0 ? Tq : Tbuf.data() + static_cast<std::size_t>(tid - 1) * T_size;
        std::size_t* z_t = tid == 0 ? z.data() : zbuf.data() + static_cast<std::size_t>(tid - 1) * N;
        double* g_t = !gsum ? nullptr : (tid == 0 ? gsum : gbuf.data() + static_cast<std::size_t>(tid - 1) * N);
        const std::size_t i_begin = N * static_cast<std::size_t>(tid) / static_cast<std::size_t>(nth);
        const std::size_t i_end = N * static_cast<std::size_t>(tid + 1) / static_cast<std::size_t>(nth);
        RowAcc acc;
        if (eps > 0.0) {
            // rounded rods (the quench): every listed pair whose rounded bounding spheres overlap, cold Newton each time
            RodPairOut po;
            const double grow = (1.0 + eps) * (1.0 + eps);
            for (std::size_t i = i_begin; i < i_end; ++i) {
                const std::uint32_t* row = pairs_data + pair_offsets[i];
                const std::size_t cnt = row[0];
                for (std::size_t q = 1; q <= cnt; ++q) {
                    const std::size_t j = row[q];
                    const double rb = 0.5 * (Dbd[i] + Dbd[j]);
                    double r[MD], d2 = 0.0;
                    for (std::size_t d = 0; d < n; ++d) {
                        double delta = x[j * n + d] - x[i * n + d];
                        delta -= (delta > half_box[d] ? boxd[d] : 0.0);
                        delta += (delta < -half_box[d] ? boxd[d] : 0.0);
                        r[d] = delta; d2 += delta * delta;
                    }
                    if (!(d2 < rb * rb * grow)) continue;
                    if (!rod_pair_rounded(n, kind, eps, r, axis.data() + i * n, hl[i], rad[i], axis.data() + j * n, hl[j], rad[j], K, po,
                                          T_t + i * nrot, T_t + j * nrot)) continue;
                    for (std::size_t d = 0; d < n; ++d) { F_t[i * n + d] += po.dEdr[d]; F_t[j * n + d] -= po.dEdr[d]; }
                    z_t[i]++; z_t[j]++;
                    if (g_t) { g_t[i] += po.f_norm; g_t[j] += po.f_norm; }
                    acc.U += po.overlap; ++acc.count; acc.L += po.sigma;
                    if (po.overlap > acc.mx) acc.mx = po.overlap;
                    acc.dk += po.dkappa; acc.Fm += po.F_mag; acc.E += po.energy;
                }
            }
        } else
        switch (n) {
            case 2: rod_rows<2>(n, kind, sh.n_iter, K, N, pairs_data, pair_offsets, x, Dbd, boxd, half_box, axis.data(), hl.data(), rad.data(), F_t, T_t, z_t, g_t, acc, cache, sh.track_tol, i_begin, i_end); break;
            case 3: rod_rows<3>(n, kind, sh.n_iter, K, N, pairs_data, pair_offsets, x, Dbd, boxd, half_box, axis.data(), hl.data(), rad.data(), F_t, T_t, z_t, g_t, acc, cache, sh.track_tol, i_begin, i_end); break;
            default: rod_rows<0>(n, kind, sh.n_iter, K, N, pairs_data, pair_offsets, x, Dbd, boxd, half_box, axis.data(), hl.data(), rad.data(), F_t, T_t, z_t, g_t, acc, cache, sh.track_tol, i_begin, i_end); break;
        }
        accs[static_cast<std::size_t>(tid)] = acc;
    }
    RowAcc acc = accs[0];
    for (int tt = 1; tt < nthreads; ++tt) {
        const double* F_o = Fbuf.data() + static_cast<std::size_t>(tt - 1) * F_size;
        const double* T_o = Tbuf.data() + static_cast<std::size_t>(tt - 1) * T_size;
        const std::size_t* z_o = zbuf.data() + static_cast<std::size_t>(tt - 1) * N;
        for (std::size_t k = 0; k < F_size; ++k) F[k] += F_o[k];
        for (std::size_t k = 0; k < T_size; ++k) Tq[k] += T_o[k];
        for (std::size_t k = 0; k < N; ++k) z[k] += z_o[k];
        if (gsum) { const double* g_o = gbuf.data() + static_cast<std::size_t>(tt - 1) * N; for (std::size_t k = 0; k < N; ++k) gsum[k] += g_o[k]; }
        const RowAcc& a = accs[static_cast<std::size_t>(tt)];
        acc.U += a.U; acc.L += a.L; acc.dk += a.dk; acc.Fm += a.Fm; acc.E += a.E; acc.count += a.count;
        if (a.mx > acc.mx) acc.mx = a.mx;
    }
    std::size_t count = 1 + acc.count;
    U = acc.U; Lc = acc.L; dkappa = acc.dk; Fmean = acc.Fm; energy = acc.E; max_min_dist = acc.mx;
    // flat hard walls: half-width h along the wall normal in the role of the sphere's radius, the sphere code's spring
    for (std::size_t d = 0; d < n; ++d) {
        if (walls[d] != 1) continue;
        for (std::size_t i = 0; i < N; ++i) {
            const double h = std::sqrt(Amat_work[i * nn + d * n + d]);
            for (std::size_t wall = 0; wall < 2; ++wall) {
                const double dxw = box[d] * static_cast<double>(wall) - x[i * n + d];
                const double dist = std::abs(dxw);
                if (!(dist < h)) continue;
                const double F_mag = -2.0 * K * (h / dist - 1.0);
                U += 2.0 * (1.0 - dist / h);
                ++count;
                F[i * n + d] += F_mag * dxw;
                dkappa += 2.0 * K * h * (dist - h);
                Fmean += F_mag;
                energy += K * (h - dist) * (h - dist);
                // E = K (h - dist)^2, dh = e.(Omega p) with p the support point along e_d: T_ab = -2K(h - dist)(e_a p_b - e_b p_a)
                const double* u = axis.data() + i * n;
                const double ud = u[d];
                double p[MD], e[MD];
                for (std::size_t k = 0; k < n; ++k) { e[k] = k == d ? 1.0 : 0.0; p[k] = hl[i] * sgn(ud) * u[k]; }
                if (sh.rod_round > 0.0) rounded_support(n, kind, sh.rod_round, u, hl[i], rad[i], e, p, nullptr, 0, nullptr);
                else if (kind == SHAPE_ROD_CAPPED) p[d] += rad[i];
                else {
                    double sp = 0.0;
                    for (std::size_t k = 0; k < n; ++k) { const double v = e[k] - ud * u[k]; sp += v * v; }
                    sp = std::sqrt(sp);
                    if (sp > 1e-12) for (std::size_t k = 0; k < n; ++k) p[k] += rad[i] * (e[k] - ud * u[k]) / sp;
                }
                add_torque(n, e, p, -2.0 * K * (h - dist), Tq + i * nrot);
            }
        }
    }
    // Curved container: point contacts at the rod's ends (curved_wall.hpp); the same in a quench (no rounding).
    if (walls[0] < 0) {
        const std::size_t t = curved_axes(walls[0], n);
        const double Rc = 0.5 * box[0];
        for (std::size_t i = 0; i < N; ++i) {
            CurvedWallSum cw;
            if (!shape_curved_wall(sh, n, t, Rc, K, x + i * n, Db[i], nullptr, Rot + i * nn, sh.qtab + i * n, F + i * n, Tq + i * nrot, cw)) continue;
            U += cw.U; count += cw.count; dkappa += cw.dkappa; Fmean += cw.Fmean; energy += cw.energy;
        }
    }
    U = std::pow(U / double(count), 2);
    Lc = Lc / double(count);
    double sumDeq = 0.0;
    for (std::size_t i = 0; i < N; ++i) sumDeq += Db[i] * sh.deqtab[i];
    dkappa = dkappa + mu * sumDeq;
    Fmean = -Fmean / count;
}

double rod_min_contact_scale(std::size_t N, const EllShape& sh, const double* x, const double* Db,
                             const double* Rot, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, std::size_t& pair_i, std::size_t& pair_j,
                             std::size_t& n_overlapping) {
    const std::size_t n = sh.Ndim, nn = n * n;
    const int kind = sh.kind;
    double inv_box[MD] = {0.0};
    for (std::size_t d = 0; d < n; ++d) inv_box[d] = walls[d] == 0 ? 1.0 / box[d] : 0.0;
    std::vector<double> axis(N * n), hl(N), rad(N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t k = 0; k < n; ++k) axis[i * n + k] = Rot[i * nn + k * n];
        hl[i] = Db[i] * sh.qtab[i * n];
        rad[i] = Db[i] * sh.qtab[i * n + 1];
    }
    double best = 1.0;
    std::size_t bi = N, bj = N, nover = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            const double rb = 0.5 * (Db[i] + Db[j]);
            double r[MD], d2 = 0.0;
            for (std::size_t d = 0; d < n; ++d) {
                double delta = x[j * n + d] - x[i * n + d];
                delta -= std::nearbyint(delta * inv_box[d]) * box[d];
                r[d] = delta; d2 += delta * delta;
            }
            if (!(d2 < rb * rb)) continue;
            const double s = rod_pair_scale_exact(n, kind, r, axis.data() + i * n, hl[i], rad[i], axis.data() + j * n, hl[j], rad[j]);
            if (s < 1.0) { ++nover; if (s < best) { best = s; bi = i; bj = j; } }
        }
        for (std::size_t d = 0; d < n; ++d) {
            if (walls[d] != 1) continue;
            const double h = rod_half_width(n, kind, Rot + i * nn, hl[i], rad[i], d);
            for (std::size_t wall = 0; wall < 2; ++wall) {
                const double dist = std::abs(box[d] * static_cast<double>(wall) - x[i * n + d]);
                if (dist < h) { ++nover; const double s = dist / h; if (s < best) { best = s; bi = i; bj = N; } }
            }
        }
        if (walls[0] < 0) {   // curved container: the factor that brings the farthest point back onto the wall
            const std::size_t t = curved_axes(walls[0], n);
            const double Rc = 0.5 * box[0];
            double xc[MD];
            for (std::size_t k = 0; k < n; ++k) xc[k] = k < t ? x[i * n + k] - Rc : 0.0;
            const double s = rod_wall_scale(n, t, kind, axis.data() + i * n, hl[i], rad[i], xc, Rc, 1.0);
            if (s < 1.0) { ++nover; if (s < best) { best = s; bi = i; bj = N; } }
        }
    }
    pair_i = bi; pair_j = bj; n_overlapping = nover;
    return best;
}

bool rod_exchange_admissible(double Dbi, const double* qi, double Dbj, const double* qj, double w) {
    const double Ri = Dbi * qi[1], Rj = Dbj * qj[1], Li = Dbi * qi[0], Lj = Dbj * qj[0];
    if (!(Ri > 0.0) || !(Rj > 0.0)) return false;
    if ((Ri > Rj ? Ri / Rj : Rj / Ri) - 1.0 > w) return false;
    if (Li == 0.0 && Lj == 0.0) return true;
    if (!(Li > 0.0) || !(Lj > 0.0)) return false;
    return (Li > Lj ? Li / Lj : Lj / Li) - 1.0 <= w;
}

}  // namespace rcpgenerator
