#pragma once
// Ellipsoid extension of the sphere packer (EllipsoidalPacking project).
//
// Everything specific to ellipsoids lives in this header and in src/core/ellipsoid.cpp. The packing loop in
// rcp_generator.cpp calls into it only when aspect ratios AND orientations were supplied; with neither, no
// function declared here is called and the sphere path is untouched.
//
// Conventions
//   * User input: D = full length of a particle along its FIRST body axis; axis k has length D * ratio_k,
//     with Ndim-1 ratios supplied (ratio_1 = 1 implied).
//   * Inside the loop the per-particle size array holds the BOUNDING diameter  Db = D * rmax,
//     rmax = max(1, ratios). All neighbour-list machinery of the sphere code (shell, refresh budgets,
//     capacities, tree pruning) then works unchanged and conservatively: two ellipsoids can touch only if
//     their bounding spheres overlap, and a rotation never changes a bounding sphere.
//   * Semi-axis k of particle i = (Db_i / 2) * q[k],  q[k] = ratio_k / rmax  (so max q = 1).
//   * Orientation: one rotation matrix per particle, row-major Ndim x Ndim, COLUMNS = body axes in the lab
//     frame. Shape matrix A = R diag(semi-axis^2) R^T.
//   * Contact function (Perram-Wertheim):  s^2 = max_{lam in [0,1]} lam(1-lam) r^T [(1-lam)A_i + lam A_j]^-1 r.
//     s < 1 <=> overlap. sigma = |r| / s is the centre distance at which the pair would touch in the same
//     direction with the same orientations (sigma = (D_i+D_j)/2 for spheres).
//   * Energy: (K/2)(sigma - |r|)^2 for s < 1 — the sphere code's spring on absolute overlap with sigma in
//     the role of r_ij. Force, generalised torque and the kappa gradient are analytic at the lam used.
//   * Generalised torque: for the rotation dR = Omega R (Omega skew), T_ab = -dE/dOmega_ab, stored for
//     a < b in lexicographic order: nrot = Ndim (Ndim-1)/2 numbers per particle.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rcpgenerator {

constexpr std::size_t ELL_MAX_DIM = 8;   // same hard limit as the sphere kernel's fixed arrays

struct EllShape {
    std::size_t Ndim = 0;
    double q[ELL_MAX_DIM] = {1, 1, 1, 1, 1, 1, 1, 1};   // semi-axis k = (Db/2) * q[k]; max q == 1
    double rmax = 1.0;         // Db = D_user * rmax
    double vol_factor = 1.0;   // particle volume / bounding-sphere volume = prod_k q[k]
    double deq_factor = 1.0;   // volume-equivalent diameter / Db = vol_factor^(1/Ndim)
    int n_iter = 5;            // bound on the lam iterations per pair per step (n_iter + 1 evaluations)
    double track_tol = 0.0;    // > 0: accept a pair's first evaluation once the Newton step predicts a relative
                               //      gain in s^2 below this (pairs start from their previous lam);
                               // 0: every pair in contact gets exactly n_iter + 1 evaluations from a cold start
    double conv_tol = 0.0;     // > 0: iterate lam until its step is below this (still at most n_iter evaluations);
                               //      used by the quench, where forces must be exact to rounding
    bool all_one = true;       // every ratio == 1 (the sphere limit run through this code)
    double rod_round = 0.0;    // rods, > 0: the ROUNDED rod of rod.hpp (the quench's smooth contact), relative rounding
    int kind = 0;              // 0 ellipsoid; 1 capped rod; 2 uncapped rod (rod.hpp: the rows of qtab then hold
                               // half-length / Db and radius / Db, and the pair kernel is rod_forces)
    // Per-particle shapes. When qtab is set (the packing loop, evaluate_once and contact_scale always set it)
    // particle i has semi-axis k = (Db_i/2) * qtab[i*Ndim + k] and volume-equivalent diameter Db_i * deqtab[i];
    // q, vol_factor and deq_factor above are then not used. The arrays belong to the caller and follow the
    // particles (spatial reorder, exchanges); the pointers must be refreshed if the arrays are reallocated.
    const double* qtab = nullptr;     // N * Ndim
    const double* deqtab = nullptr;   // N
};

// The order rule on one shape: Ndim-1 finite numbers with 1 >= ratio_2 >= ratio_3 >= ... > 0 (D, the first
// axis, is the longest). Throws std::runtime_error naming the rule; `who` says which particle (or "").
void ell_check_ratios(const double* ratios, std::size_t count, std::size_t Ndim, const std::string& who);

// Per-particle shape table from either one list (applied to every particle) or N rows. Fills qtab (N * Ndim,
// q[0] = 1) and deq (N; exactly 1.0 for a row of ones); all_one = every entry is 1.
void ell_build_shape_table(const std::vector<double>& single, const std::vector<std::vector<double>>& rows,
                           std::size_t N, std::size_t Ndim, std::vector<double>& qtab, std::vector<double>& deq,
                           bool& all_one);

// Total particle volume / total bounding-sphere volume = sum_i prod_k q_ik D_i^Ndim / sum_i D_i^Ndim. It does
// not change under a uniform rescale, a reorder, or an exchange of whole particles.
double ell_volume_factor(std::size_t N, std::size_t Ndim, const double* D, const double* qtab);

// Throws std::runtime_error on a wrong count or a non-positive ratio.
EllShape make_ell_shape(const std::vector<double>& aspect_ratios, std::size_t Ndim, int n_iter);

// Per-pair lam carried from one step to the next: for each particle i, up to CAP entries (j, lam) for the
// pairs (i, j > i) evaluated in the last force pass, in ascending j (the order of the neighbour row). A pair
// not found here starts cold. Must be cleared whenever particle indices change (the spatial reorder).
struct EllPairCache {
    static constexpr std::size_t CAP = 32;
    std::vector<std::uint32_t> ids;     // N * CAP
    std::vector<double> lam;            // N * CAP
    std::vector<std::uint32_t> count;   // N
    void reset(std::size_t N) { ids.assign(N * CAP, 0); lam.assign(N * CAP, 0.0); count.assign(N, 0); }
    void clear() { std::fill(count.begin(), count.end(), 0u); }
};

inline std::size_t ell_nrot(std::size_t Ndim) { return Ndim * (Ndim - 1) / 2; }

// Validate / orthonormalise a caller-supplied rotation (row-major Ndim x Ndim). Throws if it is not close to
// orthonormal.
void ell_check_orientation(const double* R, std::size_t Ndim, std::size_t index);

// Uniform random rotations (Haar), deterministic in seed. out is N * Ndim * Ndim.
void ell_random_orientations(std::size_t N, std::size_t Ndim, std::uint64_t seed, std::vector<double>& out);

// One pair, fixed cost. Ai, Aj: shape matrices (row-major Ndim x Ndim). r = x_j - x_i (minimum image applied).
// Returns false when the pair is certainly apart (a lower bound on s^2 reached 1). Otherwise fills
// s2 (a LOWER bound on the true s^2 — it can create a phantom contact, never hide an overlap), lam and
// xvec = Y^-1 r at that lam.
bool ell_pair_contact(const double* Ai, const double* Aj, const double* r, std::size_t Ndim, int n_iter,
                      double& s2, double& lam, double* xvec);

// Converged s^2 for one pair (bracketed iteration to machine precision; used only for reporting).
double ell_pair_s2_exact(const double* Ai, const double* Aj, const double* r, std::size_t Ndim);

// Shape matrices for all particles: Amat[i] = R_i diag(((Db_i/2) q_k)^2) R_i^T.
void ell_shape_matrices(std::size_t N, const EllShape& sh, const double* Db, const double* Rot, double* Amat);

// Full force pass over the half neighbour lists; same outputs and conventions as get_forces_nd_3, plus the
// generalised torques Tq (N * nrot) and the true energy. Flat hard walls (walls[d] == 1) are included.
void ell_forces(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                std::size_t N, const EllShape& sh, const std::vector<double>& Db,
                const std::vector<double>& box, const std::vector<std::int8_t>& walls, const double* Rot,
                double* F, double* Tq, double& U, double& max_min_dist, double& Lc, double& Fmean, double mu,
                double& dkappa, std::vector<std::size_t>& z, double& energy, std::vector<double>& Amat_work,
                EllPairCache* cache = nullptr, double* gsum = nullptr);
// gsum (optional, length N): per particle, the sum over its pair contacts of the contact force magnitude |f|
// (the quench's load measure; the sphere kernel's is the sum of overlaps, the same thing at K = 1).

// Rotate one particle: R <- (I - W/2)^-1 (I + W/2) R with the skew matrix W built from w (Ndim(Ndim-1)/2
// components, a < b in lexicographic order, W_ab = w, W_ba = -w): the same map the Adam orientation step uses.
void ell_rotate_by(std::size_t Ndim, double* R, const double* w);

// Adam step on the orientations, mirroring the fused position step (same beta1, beta2, t, epsilon). The angle
// step of particle i is alpha / (Db_i/2): the tip of the bounding sphere moves by at most alpha per component.
void ell_rotation_step(std::size_t N, std::size_t Ndim, const double* Tq, double* m, double* v, double* Rot,
                       const double* Db, double alpha, double beta1, double beta2, double t, double epsilon);

// Exact shrink-to-contact: the smallest converged s over all pairs (minimum image on periodic axes) and all
// flat-wall contacts, or 1 if nothing overlaps. Scaling every particle by this factor leaves the tightest
// pair exactly touching. O(N^2) with a bounding-sphere prefilter. pair_i/pair_j receive the tightest pair
// (pair_j == N for a wall).
double ell_min_contact_scale(std::size_t N, const EllShape& sh, const double* x, const double* Db,
                             const double* Rot, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, std::size_t& pair_i, std::size_t& pair_j,
                             std::size_t& n_overlapping);

// ---- diameter permutation with ellipsoids (BUILD.md section 8) -------------------------------------------------
// The move: two particles trade places. At fixed positions, slot i and slot j exchange diameter (D and D0), shape
// (row of qtab, deq) and, when move_orientations, the orientation. With move_orientations == false the
// orientations stay with the sites (the comparison case of target D6).
void ell_exchange(std::size_t Ndim, std::size_t i, std::size_t j, double* D, double* D0, double* qtab, double* deq,
                  double* Rot, bool move_orientations);

// Which pairs may be exchanged: every axis length must agree within the relative window w,
//     max(L_ik / L_jk, L_jk / L_ik) - 1 <= w   for k = 1..Ndim,   L_pk = D_p * q_pk
// (axis lengths line up entry by entry because ratios are listed in decreasing order). For spheres this is the
// sphere engine's window |D_j / D_i - 1| <= w on the diameter.
bool ell_exchange_admissible(std::size_t Ndim, double Di, const double* qi, double Dj, const double* qj, double w);

struct EllPermuteParams {
    double window = 0.05;        // relative window on every axis length (the sphere recipe's 5 %)
    long budget = 1;             // at most this many exchanges per pass
    int max_partners = 64;       // partners examined per particle per pass (strided through its window; the
                                 // stride start advances with pass_index so that every partner gets its turn)
    bool move_orientations = true;
    bool decide_in_place = false;   // judge the room of an arriving particle in the SITE's orientation (always the
                                    // case when orientations stay; with move_orientations it is a diagnostic:
                                    // decide as if the particle could take the site's orientation, then move it
                                    // with its own)
    int kind = 0;                   // EllShape::kind: for rods the cage normals, the support function and the
                                    // window (straight length and radius) come from rod.hpp
    bool overlap_guard = false;     // also require that the exchange does not raise the cage overlap of the two
                                    // particles, sum over cage planes of max(0, h_p(n) - b)^2 (how far each
                                    // particle would stick out through the planes of its site)
    std::uint64_t pass_index = 0;
    int n_iter = 5;              // contact-function iterations used to find a facet normal
};
struct EllExchangeRecord { std::uint32_t i, j; double gain; };

// One permutation pass at fixed positions (no random numbers). What plays the sphere engine's roles:
//   * "the radius that fits at a site" (c_i = min_k (r_ik - R_k) for spheres) becomes the site's CAGE: for every
//     neighbour k within the shell s_ik < 1.5 (the pair would touch if both were scaled by 1.5; for spheres
//     r < 1.5 (R_i + R_k), the sphere engine's shell) the tangent plane of k facing i, with unit normal n_ik
//     taken from the contact function, at distance b_ik = n_ik.r_ik - h_k(n_ik) from the centre of i
//     (h_k(n) = sqrt(n.A_k n) is k's support function). The planes bound a convex region free of every neighbour.
//   * the ROOM of particle p at site a is the largest factor by which p could be scaled about the centre of a
//     and still fit inside that cage: tau(a, p) = min_k b_ak / h_p(n_ak), with p in the orientation it would
//     have there (its own when orientations travel, the site's when they stay). Spheres: tau = c_a / R_p.
//   * the gain of exchanging i and j is the change of J = 2 sum_p R_p^2 tau(site(p), p), R_p = D_p / 2:
//         gain = 2 [R_i^2 tau(j,i) + R_j^2 tau(i,j) - R_i^2 tau(i,i) - R_j^2 tau(j,j)].
//     For spheres this is the sphere engine's 2 (R_j - R_i)(c_i - c_j), exactly.
// Candidates: partners within the window on every axis; accepted greedily by decreasing gain (gain > 0), each
// particle at most once per pass, at most `budget` exchanges. Walls are not facets (as in the sphere engine).
// Amat must be the shape matrices of the current state (ell_forces leaves them in its work array). Applies the
// exchanges through ell_exchange and appends them to `done`. Returns the number of candidates examined.
std::size_t ell_permute_pass(const std::uint32_t* pairs_data, const std::uint64_t* pair_offsets, const double* x,
                             std::size_t N, std::size_t Ndim, const std::vector<double>& box,
                             const std::vector<std::int8_t>& walls, double* D, double* D0, double* qtab,
                             double* deq, double* Rot, const double* Amat, const EllPermuteParams& prm,
                             std::vector<EllExchangeRecord>& done);

// Timing probe for requirement 3: average wall time per pair evaluation (seconds) of the ellipsoid pair kernel
// (energy + force + torque + kappa term, overlapping pairs) and of the sphere pair body, on the same random
// overlapping configurations. track_tol == 0: cold start, exactly n_iter + 1 evaluations (t_ellipsoid).
// t_tracked: the same pairs started from their converged lam with tracking tolerance 1e-6, the steady state
// of a pair inside a run (one evaluation).
void ell_pair_cost_probe(std::size_t Ndim, const std::vector<double>& aspect_ratios, int n_iter,
                         std::size_t n_pairs, std::size_t repeats, std::uint64_t seed,
                         double& t_ellipsoid, double& t_sphere, double& t_tracked, double& checksum);

}  // namespace rcpgenerator
