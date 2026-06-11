#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace rcpgenerator {

struct PackingConfig {
    std::vector<double> box;
    // Legacy packer behavior: empty walls defaults to Ndim zeros before first use.
    std::vector<std::int8_t> walls;
    std::uint32_t neighbor_max = 0;
    std::uint32_t seed = 0;
    bool fix_height = false;
};

struct PackingInput {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
};

struct PackingResult {
    std::vector<std::vector<double>> positions;
    std::vector<double> diameters;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    std::vector<double> phi_history;
    std::vector<double> force_history;
    std::vector<double> energy_history;
    // Cycle 16: per-step diagnostics for ADAM convergence study. All four
    // have length = hist_size (allocated alongside phi_history/force_history).
    std::vector<double> max_overlap_history;  // max pair overlap each step
    std::vector<double> mu_history;           // chemical-potential rung value
    std::vector<double> alpha_history;        // current learning rate
    std::vector<std::int8_t> mu_flag_history; // -1, 0, 1 (rung index)
    std::size_t steps = 0;
    double phi = 0.0;
    double max_min_dist = 0.0;
    double force_magnitude = 0.0;
    // Cycle 10: particle_origin[k] = the index this particle had at the
    // start of the run, before any spatial reorder. Lets callers compare
    // a candidate's positions against an oracle's positions on an
    // apples-to-apples basis (oracle.positions[particle_origin[k]] should
    // ~equal candidate.positions[k] within FP noise).
    std::vector<std::size_t> particle_origin;
    // Cycle 18: per-particle aggregate stats over the LAST window of the
    // run (window length controlled by RCP_PP_WINDOW env var, default 5000
    // steps). Used to identify "troublemakers" — particles whose force
    // oscillates, whose position swings, whose contacts churn. Empty if
    // window param is 0 or run ended too quickly.
    std::vector<double> pp_F_sq_sum;       // sum of |F_i|^2 over window
    std::vector<std::uint64_t> pp_F_sign_flips;   // total sign flips across dims
    std::vector<double> pp_pos_range_max;  // max range across dims (max-min)
    std::vector<std::uint64_t> pp_contact_sum;    // sum of pair counts (Cycle 13 own-writes)
    std::size_t pp_window_steps = 0;       // actual # steps accumulated
    // Cycle 18 Method B: count of sign-flip momentum resets per particle
    // over the window. High count = oscillating particle; low = stable.
    std::vector<std::uint64_t> pp_signflip_resets;
    // Cycle 18 Phase A: F-spike events (one row per event). Each row:
    //   [step, F_ratio_vs_smoothed, F_value, max_overlap, mu_flag,
    //    alpha, step_mod_reorder, step_mod_500]
    // Stored flat as 8 doubles per event in spike_log.
    std::vector<double> spike_log;
    // Cycle 18 Phase A2: per-reorder-event pair-count diagnostics.
    // Tracks pair list size before/after each reorder to test Hypothesis A
    // (reorder finds previously-missed pairs). Format per event (5 doubles):
    //   [step, mu_flag, total_pairs_before, total_pairs_after, F_at_step]
    std::vector<double> reorder_log;
};

struct PackingObservation {
    std::size_t step = 0;
    double phi = 0.0;
    double force_magnitude = 0.0;
    double energy = 0.0;
    double max_min_dist = 0.0;
};

struct PackingTrace {
    std::vector<std::size_t> steps;
    std::vector<std::vector<std::vector<double>>> positions;
    std::vector<std::vector<double>> diameters;
    std::vector<double> phi;
    std::vector<double> force;
    std::vector<double> energy;
    std::vector<double> max_min_dist;
};

struct PackingRunOptions {
    std::size_t max_steps = 0;
    bool has_mu_override = false;
    double mu_override = 0.0;
    bool fix_diameter = false;
    bool has_target_phi = false;
    double target_phi = 0.0;
    // Cycle 10/12: 0 disables periodic spatial reorder entirely (baseline
    // behavior, same particle indices for the whole run). Any positive
    // value sets the reorder period in steps. Default 2000 chosen by sweep
    // at N=100K — balances the L1-residence win on forces.loop against
    // the full-pair-rebuild cost incurred by each reorder.
    std::size_t reorder_interval = 2000;
    // Cycle 15o: override the dim-specific hardcoded initial packing
    // fraction (0.575 for 2D, 0.33 for 3D). Production workflow scales
    // diameters to phi=0.55 at ADAM start (ADAM eats the initial overlap
    // cleanly). The hardcoded 0.33 makes the bench start far from dense
    // regime. When has_initial_phi_target is true, use initial_phi_target
    // in scale_diameters_nd instead.
    bool has_initial_phi_target = false;
    double initial_phi_target = 0.55;
};

using PackingObserver = std::function<void(
    const PackingObservation& observation,
    const std::vector<std::vector<double>>* positions)>;

double delta_x(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& x_old,
    const std::vector<double>& D);

double norm(
    const std::vector<double>& a,
    const std::vector<double>& b);

// Cycle 7: F is now flat (vector<double> of length N*Ndim, addressed as
// F[i*Ndim + d]). Removes a pointer chase per particle in adam_update and
// in the forces reduce step, and makes the data contiguous for the
// compute_mean_force scan and for any SIMD inner loops downstream.
double compute_mean_force(
    const double* F,
    std::size_t N,
    std::size_t Ndim,
    double Lc,
    const std::vector<std::size_t>& z);

double mean(
    const std::vector<double>& v,
    std::size_t start,
    std::size_t end);

double max_elem(const std::vector<std::vector<double>>& M);
double compute_max_min_dist(const std::vector<std::vector<double>>& min_dist);
std::size_t compute_max_neighbors(const std::vector<std::vector<std::uint32_t>>& pairs);
std::size_t compute_num_changes(const std::vector<std::uint32_t>& refresh);

std::vector<std::size_t> sort_indices_by_column(
    const std::vector<std::vector<double>>& x,
    std::size_t col);

#ifndef RCPGENERATOR_SPHERE_VOLUME_INLINE
#define RCPGENERATOR_SPHERE_VOLUME_INLINE
inline double sphere_volume(double r, std::size_t Ndim) {
    return std::pow(M_PI, Ndim / 2.0) * std::pow(r, Ndim) /
           std::tgamma(Ndim / 2.0 + 1.0);
}
#endif

std::pair<std::vector<double>, double> scale_diameters_nd(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    std::size_t Ndim,
    bool fix_height);

// Cycle 13: flat pairs. Cycle 14: flat x (positions).
// x is a length-N*Ndim double*; x[i*Ndim + d] is particle i's d-th coord.
// pairs_data is one contiguous vector<uint32_t>; pair_offsets[i] is the start
// of particle i's slot. pairs_data[pair_offsets[i]] is its current neighbor
// count, pairs_data[pair_offsets[i] + 1 .. + count] are the neighbor indices.
// max_neighbors_per_particle[i] is the cap for that slot (derived from D[i]).
// Total storage = pair_offsets[N].
void get_pairs_nd_3(
    std::size_t N,
    std::size_t Ndim,
    const double* x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    const std::vector<std::uint32_t>& refresh,
    const std::uint32_t* max_neighbors_per_particle,
    const std::uint64_t* pair_offsets,
    std::uint32_t* pairs_data,
    std::uint32_t* out_max_observed,   // optional, may be nullptr
    std::uint32_t* out_min_observed,   // optional, may be nullptr
    // Cycle 21k packed-D0: parallel per-slot copy of D0[j], refilled in the
    // canonicalize pass so the force loop can replace its D[j] gather with
    // a contiguous load (Dj = packed*kappa is the same single rounding as
    // D[j] = D0[j]*kappa, so forces are byte-identical). Defaults nullptr:
    // sentinel/shadow scratch builds MUST leave these null so they never
    // overwrite the live packed array.
    double* packed_d0 = nullptr,
    const double* D0_data = nullptr);

void get_forces_nd_3(
    const std::uint32_t* pairs_data,
    const std::uint64_t* pair_offsets,
    const double* x,                    // flat [N*Ndim] — Cycle 14
    std::size_t N,
    std::size_t Ndim,
    const std::vector<double>& D,
    const std::vector<double>& box,
    const std::vector<std::int8_t>& walls,
    double* F,                                  // flat [N*Ndim] — Cycle 7
    double& U,
    std::vector<std::vector<double>>& min_dist,
    double& max_min_dist,
    double& Lc,
    double& Fmean,
    double mu,
    double& dkappa,
    std::vector<std::size_t>& z,
    // Cycle 21k packed-D0 (see get_pairs_nd_3): when non-null, the 2D/3D
    // SIMD reject paths load D0[j] contiguously from this array and
    // multiply by kappa instead of gathering D[j]. kappa MUST be the same
    // value the D-update pass used (D[j] == D0[j]*kappa bit-for-bit).
    const double* packed_d0 = nullptr,
    double kappa = 1.0);

// Cycle 5a: m, v, v_update, m_hat, v_hat moved to flat double* layout.
// They are addressed as [k * Ndim + d] for particle k, dimension d. The
// previously-nested ones used vector<vector<double>>, which forces a pointer
// chase per particle and was capping adam_update T-scaling at ~1.08×.
// v_max, a, v_verlet, a_old remain nested — they are touched only by the
// (unused-in-ADAM) Verlet code path and aren't on the hot path.
void adam_update(
    const std::string& method,
    std::size_t N,
    std::size_t Ndim,
    const double* F,                            // flat [N*Ndim] — Cycle 7
    double dkappa,
    double beta1,
    double beta2,
    std::size_t t,
    double dt,
    double verlet_drag,
    double* m,
    double* v,
    std::vector<std::vector<double>>& v_max,
    double* v_update,
    double* m_hat,
    double* v_hat,
    std::vector<std::vector<double>>& a,
    std::vector<std::vector<double>>& v_verlet,
    std::vector<std::vector<double>>& a_old,
    double& m_kappa,
    double& v_kappa,
    double& v_update_kappa,
    double& m_hat_kappa,
    double& v_hat_kappa);

PackingResult run_packing(
    const PackingInput& input,
    const PackingConfig& config);

std::pair<PackingResult, PackingTrace> run_packing_observed(
    const PackingInput& input,
    const PackingConfig& config,
    std::size_t progress_interval,
    bool capture_positions,
    std::size_t trajectory_interval,
    const PackingObserver& observer,
    const PackingRunOptions& options = PackingRunOptions{});

}  // namespace rcpgenerator
