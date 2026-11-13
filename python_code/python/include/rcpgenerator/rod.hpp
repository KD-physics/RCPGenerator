#pragma once
// Rods (BUILD.md section 9): capped (a straight section with a half-sphere on each end; 2D: a stadium) and
// uncapped (a flat-ended cylinder; 2D: a rectangle). Additive to the ellipsoid extension: the packing loop treats
// a rod packing as a non-spherical packing (the `ell_on` path: orientations, torques, per-particle shape rows) and
// only the pair kernel, the wall contact, the volume and the permutation's support function differ, selected by
// EllShape::kind. Nothing here is called for spheres or for ellipsoids.
//
// Conventions
//   * User input: rod_length L and rod_radius R in units of the particle's D: straight section D*L, radius D*R.
//     Capped: tip to tip D*(L + 2R); L = 0 is a sphere of radius D*R. Uncapped: full length D*L.
//   * The rod's axis is the FIRST body axis: column 0 of the particle's rotation matrix.
//   * Inside the loop the size array holds the BOUNDING diameter Db = D * bfac, bfac = L + 2R (capped) or
//     2 sqrt((L/2)^2 + R^2) (uncapped), so the sphere neighbour machinery works unchanged. The shape row of a
//     particle (the same N x Ndim table ellipsoids use) holds q[0] = half-length / Db and q[1] = radius / Db
//     (other entries 1). `bfac` is kept per particle to give the user's D back.
//   * Capped pair energy: (K/2)(R_i + R_j - d)^2 with d the distance between the two axis segments: the sphere
//     spring on the closest points. Closed form, any dimension.
//   * Uncapped pair energy: the ellipsoid form (K/2)(sigma - |r|)^2, sigma = |r| / s, with s the factor by which
//     both rods would have to be scaled about their centres to just touch:
//         s = max over unit n of  n.r / (h_i(n) + h_j(n)),   h(n) = hl |n.u| + R |n - (n.u) u|   (support function).
//     Any n gives a lower bound on s (a phantom contact at worst, never a missed overlap). 2D: the maximum is at
//     one of the four edge normals (exact). 3D: closed-form candidates (end faces, side against side) verified
//     by where their contact point falls, else a bounded Newton iteration on the convex problem
//     min h_i(n) + h_j(n) subject to n.r = 1. Dimensions above 3 are refused for uncapped rods.
//   * Force, torque and the kappa gradient are analytic: with n the maximising normal and p_i, p_j the support
//     points of the two rods along n, d s = ... (see rod.cpp). Torques use the generator convention of
//     ellipsoid.hpp: T_ab = -dE/dOmega_ab for the rotation dR = Omega R, a < b.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "rcpgenerator/ellipsoid.hpp"

namespace rcpgenerator {

constexpr int SHAPE_ELLIPSOID = 0, SHAPE_ROD_CAPPED = 1, SHAPE_ROD_UNCAPPED = 2;

// Validate the rod inputs and build the per-particle rows. lengths and radii have 1 or N entries each.
// qtab: N * Ndim (q[0] = half-length / Db, q[1] = radius / Db, rest 1); deq: volume-equivalent diameter / Db;
// bfac: Db / D. Throws std::runtime_error with the rule that was broken.
void rod_build_shape_table(const std::vector<double>& lengths, const std::vector<double>& radii, bool capped,
                           std::size_t N, std::size_t Ndim, std::vector<double>& qtab, std::vector<double>& deq,
                           std::vector<double>& bfac);

// Volume of one rod / volume of its bounding sphere, from its row (q0 = half-length / Db, q1 = radius / Db).
double rod_row_volume_factor(std::size_t Ndim, double q0, double q1, bool capped);

// Total rod volume / total bounding-sphere volume (the loop's one volume factor).
double rod_volume_factor(std::size_t N, std::size_t Ndim, const double* Db, const double* qtab, bool capped);

// Half-width of a rod along lab axis d (its support function along e_d): the wall contact's "radius".
double rod_half_width(std::size_t Ndim, int kind, const double* R, double hl, double rad, std::size_t d);

// Support function of a rod along the unit vector nv (axis u, half-length hl, radius rad).
double rod_support(std::size_t Ndim, int kind, const double* u, double hl, double rad, const double* nv);

// One pair. r = x_j - x_i (minimum image applied); ui, uj unit axes; hl half-lengths, rad radii (actual sizes).
// max_iter bounds the 3D uncapped iteration (>= 1). Returns false when the pair does not overlap. Otherwise
// fills: s (contact scale < 1 for uncapped; d / (R_i + R_j) for capped), overlap = 1 - s, sigma, F_mag, dkappa,
// energy, f_norm, dEdr (force on i = +dEdr, on j = -dEdr), Ti, Tj (generalised torques, ADDED to the arrays),
// and nvec (the contact normal from i to j).
struct RodPairOut {
    double s = 1.0, overlap = 0.0, sigma = 0.0, F_mag = 0.0, dkappa = 0.0, energy = 0.0, f_norm = 0.0;
    double dEdr[ELL_MAX_DIM] = {0}, nvec[ELL_MAX_DIM] = {0};
    double pi[ELL_MAX_DIM] = {0}, pj[ELL_MAX_DIM] = {0};   // where the pair force acts, relative to each centre
                                                           // (uncapped: the support points along +nvec of both)
    bool converged = true;                                 // 3D uncapped: an optimality condition was met
};
bool rod_pair(std::size_t Ndim, int kind, const double* r, const double* ui, double hli, double radi,
              const double* uj, double hlj, double radj, double K, int max_iter, RodPairOut& out, double* Ti,
              double* Tj);

// Normal of a pair for the permutation's cage planes, at any separation: capped, the direction between the closest
// points of the two axis segments; uncapped, the maximising normal. Returns the scale (capped: segment distance /
// (R_i + R_j)); a return value >= 1 means the pair is apart at this r.
double rod_pair_normal(std::size_t Ndim, int kind, const double* r, const double* ui, double hli, double radi,
                       const double* uj, double hlj, double radj, int max_iter, double* nv);

// Contact scale of one pair, converged: the factor by which both rods must be scaled about their centres to just
// touch (< 1: overlap; may exceed 1). Used for the exact shrink and for tests.
double rod_pair_scale_exact(std::size_t Ndim, int kind, const double* r, const double* ui, double hli, double radi,
                            const double* uj, double hlj, double radj);

// What the 3D uncapped kernel found for a pair at the last pass, so that the next pass starts from it: for each
// particle i up to CAP entries (j, tag, normal) for pairs (i, j > i). tag 4: the iteration's normal (the next
// iteration starts there and usually needs one step); tag 5: a normal that proved the pair apart (checked first);
// other tags: nothing to reuse. Cleared when particle indices change; entries of exchanged particles are dropped.
struct RodPairCache {
    static constexpr std::size_t CAP = 32;
    std::vector<std::uint32_t> ids, count;
    std::vector<std::int8_t> tag;
    std::vector<double> nvec;   // 3 per entry
    void reset(std::size_t N) { ids.assign(N * CAP, 0); count.assign(N, 0); tag.assign(N * CAP, 0); nvec.assign(N * CAP * 3, 0.0); }
    void clear() { std::fill(count.begin(), count.end(), 0u); }
    void purge(const std::vector<char>& touched) {
        const std::size_t N = count.size();
        for (std::size_t i = 0; i < N; ++i) {
            if (touched[i]) { count[i] = 0; continue; }
            std::uint32_t keep = 0;
            for (std::uint32_t c = 0; c < count[i]; ++c) {
                const std::size_t a = i * CAP + c, b = i * CAP + keep;
                if (touched[ids[a]]) continue;
                ids[b] = ids[a]; tag[b] = tag[a];
                for (int k = 0; k < 3; ++k) nvec[b * 3 + k] = nvec[a * 3 + k];
                ++keep;
            }
            count[i] = keep;
        }
    }
};

// Full force pass, same outputs and conventions as ell_forces (ellipsoid.hpp). sh.kind selects the rod kind,
// sh.qtab the rows, sh.n_iter the iteration bound of the 3D uncapped kernel. Amat_work receives, for every
// particle, a matrix whose DIAGONAL holds the squared half-widths along the lab axes (what the quench's wall
// sums read); the off-diagonal entries are zero. Threaded over the rows of the neighbour list like ell_forces
// (per-thread buffers added in thread order); with one thread the sums run in one fixed order.
void rod_forces(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                std::size_t N, const EllShape& sh, const std::vector<double>& Db,
                const std::vector<double>& box, const std::vector<std::int8_t>& walls, const double* Rot,
                double* F, double* Tq, double& U, double& max_min_dist, double& Lc, double& Fmean, double mu,
                double& dkappa, std::vector<std::size_t>& z, double& energy, std::vector<double>& Amat_work,
                double* gsum = nullptr, RodPairCache* cache = nullptr);

// Exact shrink-to-contact for rods: as ell_min_contact_scale.
double rod_min_contact_scale(std::size_t N, const EllShape& sh, const double* x, const double* Db,
                             const double* Rot, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, std::size_t& pair_i, std::size_t& pair_j,
                             std::size_t& n_overlapping);

// ---- rounded rods (the quench) ---------------------------------------------------------------------------------------
// A sharp rod's contact energy has kinks exactly where rod packings like to sit: two flat sides or ends resting on
// each other. There the pair force may act anywhere on the shared flat region, and a descent method never finds the
// balance. For the quench the rod is therefore replaced by a slightly ROUNDED rod, a smooth strictly convex body
// that contains it, with support function (eps = EllShape::rod_round, |n| = 1)
//     uncapped: h(n) = hl sqrt(c^2 + eps^2) + R sqrt(1 - c^2 + eps^2),     capped: h(n) = hl sqrt(c^2 + eps^2) + R,
// c = n.u: the axis segment becomes a needle-shaped ellipsoid of half-width eps hl and the end disc a lens of
// half-thickness eps R. The energy is the scale form (K/2)(sigma - |r|)^2, sigma = |r| / s, for BOTH kinds,
// s = max_n n.r / (h_i + h_j), found by Newton on a smooth strictly convex problem; the support points are unique,
// so force and torque are well defined everywhere and the contact point of two flat faces sits inside them. 2D and
// 3D. With eps = 0 this is not used (the sharp kernels above are).
bool rod_pair_rounded(std::size_t Ndim, int kind, double eps, const double* r, const double* ui, double hli,
                      double radi, const double* uj, double hlj, double radj, double K, RodPairOut& out, double* Ti,
                      double* Tj);

// Counters of the 3D uncapped kernel since the last reset (single-threaded diagnostic): pairs entering it after the
// bounding test, separated by the four quick axes, settled by a closed-form candidate, sent to the iteration, and
// the iterations made. rod_kernel_stats(true) resets.
struct RodKernelStats { std::uint64_t pairs = 0, quick_apart = 0, closed_form = 0, iterated = 0, iterations = 0, contacts = 0, warm = 0,
                        warm_failed = 0, backtracks = 0, descents = 0,
                        tracked_ok = 0, tracked_rejected = 0, sep_cached_ok = 0, sep_cached_failed = 0, lean_ok = 0, lean_failed = 0; };
RodKernelStats rod_kernel_stats(bool reset);

// Which pairs may be exchanged by the permutation: straight length AND radius within the relative window w
// (two zero lengths agree; a zero and a non-zero length do not).
bool rod_exchange_admissible(double Dbi, const double* qi, double Dbj, const double* qj, double w);

}  // namespace rcpgenerator
