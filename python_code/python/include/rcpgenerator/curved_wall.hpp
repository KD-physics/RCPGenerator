#pragma once
// Ellipsoids and rods against the wall of a curved container (BUILD.md section 11).
//
// The container is the sphere code's: the first t = -walls[0] axes form a t-ball of radius Rc = box[0]/2 centred
// at (Rc, ..., Rc); the remaining axes are periodic or flat-walled (t = 2 in 2D: a disk; t = 2 in 3D: a cylinder;
// t = 3 in 3D: a sphere; any t <= Ndim above). A sphere of radius r at distance d from the axis touches the wall
// when Rc - d < r. For a convex particle the same statement is about its FARTHEST point from the axis:
//
//     g = max over points p of the particle of |P (p - c)|      (P: projection on the curved axes, c: the axis)
//
// and the particle sticks out iff g > Rc. The support function along the radial direction through the particle's
// CENTRE (what a flat wall uses) is only a lower bound on g and would report "inside" for a particle that sticks
// out. With p* the farthest point, nvec = P(p* - c)/g the outward radial direction there, h = nvec.(p* - x) the
// particle's extent along nvec (r for a sphere) and delta = h - (g - Rc) (Rc - d for a sphere), the sphere code's
// wall formulas are kept with r -> h:
//     energy K (h - delta)^2 = K (g - Rc)^2,   force on the centre -2K (h - delta) nvec,   torque of that force at p*,
//     U += 1 - delta / h,   dkappa += 2K h (delta - h),   Fmean += -2K (h / delta - 1).
// By the envelope theorem the force and torque are the exact gradient of K (g - Rc)^2 wherever the farthest point
// is unique. Near the axis two lobes of a particle can be equally far: the energy is continuous there and its
// gradient has a crease (the max of two smooth functions).
//
// Nothing here is reached for spheres, or for ellipsoids and rods in periodic or flat-wall boxes.

#include <cstddef>

#include "rcpgenerator/ellipsoid.hpp"

namespace rcpgenerator {

struct CurvedWallHit {
    double g = 0.0;                    // largest distance from the container axis over the particle
    double nvec[ELL_MAX_DIM] = {0};    // outward radial unit vector at the farthest point (0 on the other axes)
    double p[ELL_MAX_DIM] = {0};       // the farthest point, relative to the particle's centre
};

// Ellipsoid with shape matrix A = R diag(semi-axes^2) R^T (n x n, row-major). xc: the centre minus the container's
// centre (n entries; only the first t matter). Exact: eigen-decomposition of the t x t block of A and the
// one-dimensional secular equation, solved to rounding with a bracketed Newton iteration (bounded).
void ell_farthest_from_axis(std::size_t n, std::size_t t, const double* A, const double* xc, CurvedWallHit& hit);

// Rod with axis u (unit), half-length hl and radius rad. kind: SHAPE_ROD_CAPPED or SHAPE_ROD_UNCAPPED (rod.hpp).
// eps = 0: the sharp rod (capped: an end point plus the radius; uncapped: a corner in 2D, a point of an end rim in
// 3D). eps > 0: the rounded rod of the quench (rod.hpp), by ascent on its support function.
void rod_farthest_from_axis(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, double eps,
                            const double* xc, CurvedWallHit& hit);

// The sphere code's wall quantities from a hit. Returns false (and touches nothing) when the particle is inside.
struct CurvedWallContact {
    double h = 0.0, delta = 0.0;       // extent along nvec, and h minus the depth (the sphere's r and Rc - d)
    double F_mag = 0.0;                // -2K (h/delta - 1)
    double U_add = 0.0;                // 1 - delta/h
    double dkappa = 0.0;               // 2K h (delta - h)
    double energy = 0.0;               // K (h - delta)^2
    double f = 0.0;                    // |force| = 2K (h - delta): what the quench adds to the contact-force sum
};
bool curved_wall_contact(std::size_t n, double Rc, double K, const CurvedWallHit& hit, CurvedWallContact& c);

// Adds the wall force to F (n entries of this particle) and the torque to Tq (n(n-1)/2 entries, the generator
// convention of ellipsoid.hpp: T_ab = -dE/dOmega_ab, a < b).
void curved_wall_apply(std::size_t n, const CurvedWallHit& hit, const CurvedWallContact& c, double* F, double* Tq);

// Largest factor s <= s_max by which the particle can be scaled about its centre and still lie inside the container
// (g(s) = Rc; bisection on the routines above; for the exact shrink). Returns s_max if it fits at s_max.
double ell_wall_scale(std::size_t n, std::size_t t, const double* A, const double* xc, double Rc, double s_max);
// Number of curved axes from walls[0] (< 0): the sphere code reads -1 as 2.
inline std::size_t curved_axes(int w0, std::size_t n) {
    const std::size_t t = w0 == -1 ? 2 : static_cast<std::size_t>(-w0);
    return t < n ? t : n;
}

// ---- the wall contact of one particle, as the force passes, the quench and the tests use it ---------------------------
// Ellipsoid: ONE contact, at its farthest point (a strictly convex smooth body has one, away from the axis).
// Rod: a rod lying along a concave wall touches it with BOTH ends, and a flat end standing on it touches with its
// whole rim. One contact at "the" farthest point would give an energy with a crease exactly at those resting
// states (the farther end keeps changing), which a minimiser cannot balance. A rod's wall energy is therefore a
// signed sum of point contacts, each K (depth)^2 with the sphere code's formulas:
//     capped:        + end 1 (a ball of the rod's radius at the end of the axis)  + end 2  - the middle ball.
//                    Minus the middle so that a rod of length 0 is exactly one sphere contact (1 + 1 - 1), and so
//                    that the energy stays continuous when a short rod is out along its whole length. It is never
//                    negative: the distance from the axis is convex along the rod, so depth_mid^2 <= (d1^2 + d2^2)/2.
//     uncapped, 2D:  + each of the four corners.
//     uncapped, 3D:  + 24 points on each end rim, fixed in the rod's own frame. A flat end standing on a spherical
//                    wall touches along its whole rim, so "the farthest rim point" is not a function of the
//                    configuration there; a fixed ring of points is, and each point contact is smooth. The points
//                    lie on the circle of radius rad / cos(pi/24) (0.9 % outside the rim), so the rim is inside
//                    their polygon and no part of it can be outside the wall unseen. The price: the wall sees a
//                    24-sided prism, and turning the rod about its own axis changes the wall energy slightly.
// The farthest point of a capped rod or of a rectangle is among its points, and a 3D rim is inside its polygon, so a
// rod that sticks out is never reported inside. Point contacts are smooth in the rod's position and orientation, so
// the same energy serves the packing loop and the quench (no rounding is involved in the wall contact of a rod).
struct CurvedWallSum {
    double U = 0.0, dkappa = 0.0, Fmean = 0.0, energy = 0.0;   // to add to the pass's sums
    double f = 0.0;        // sum of the contact-force magnitudes (the quench's per-particle sum)
    int count = 0;         // contacts (signed sum: 1 for a sphere-like contact)
};
// x_i: the particle's position; Db_i: its bounding diameter; ellipsoid: Amat_i (shape matrix); rod: Rot_i (rotation
// matrix; the axis is column 0) and q_i (its row of the shape table). Adds force to F_i and torque to Tq_i when they
// are not null. Returns false, having computed nothing, when the bounding sphere is inside, or nothing touches.
bool shape_curved_wall(const EllShape& sh, std::size_t n, std::size_t t, double Rc, double K, const double* x_i, double Db_i,
                       const double* Amat_i, const double* Rot_i, const double* q_i, double* F_i, double* Tq_i,
                       CurvedWallSum& sum);
// The rod part alone (tests): xc = centre minus container centre.
// Rot: the rod's rotation matrix (row-major, columns = body axes), used by 3D flat-ended rods for the frame of
// their rim points; may be null (a frame is then built from u, which is enough for one evaluation but does not turn
// with the rod about its own axis).
bool rod_curved_wall(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, const double* xc,
                     double Rc, double K, double* F, double* Tq, CurvedWallSum& sum, const double* Rot = nullptr);
double rod_wall_scale(std::size_t n, std::size_t t, int kind, const double* u, double hl, double rad, const double* xc,
                      double Rc, double s_max);

}  // namespace rcpgenerator
